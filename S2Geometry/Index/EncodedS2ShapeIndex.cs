// EncodedS2ShapeIndex is an S2ShapeIndex implementation that works directly
// with encoded data.  Rather than decoding everything in advance, geometry is
// decoded incrementally (down to individual edges) as needed.  It can be
// initialized from a single block of data in nearly constant time (about 1.3
// microseconds per million edges).  This saves large amounts of memory and is
// also much faster in the common situation where geometric data is loaded
// from somewhere, decoded, and then only a single operation is performed on
// it.  It supports all S2ShapeIndex operations including boolean operations,
// measuring distances, etc.
//
// The speedups can be over 1000x for large geometric objects.  For example
// vertices and 50,000 loops.  If this geometry is represented as an
// S2Polygon, then simply decoding it takes ~250ms and building its internal
// S2ShapeIndex takes a further ~1500ms.  These times are much longer than the
// time needed for many operations, e.g. e.g. measuring the distance from the
// polygon to one of its vertices takes only about 0.001ms.
//
// If the same geometry is represented using EncodedLaxPolygonShape and
// EncodedS2ShapeIndex, initializing the index takes only 0.005ms.  The
// distance measuring operation itself takes slightly longer than before
// (0.0013ms vs. the original 0.001ms) but the overall time is now much lower
// (~0.007ms vs. 1750ms).  This is possible because the new classes decode
// data lazily (only when it is actually needed) and with fine granularity
// (down to the level of individual edges).  The overhead associated with this
// incremental decoding is small; operations are typically 25% slower compared
// to fully decoding the MutableS2ShapeIndex and its underlying shapes.
//
// EncodedS2ShapeIndex also uses less memory than MutableS2ShapeIndex.  The
// encoded data is contiguous block of memory that is typically between 4-20%
// of the original index size (see MutableS2ShapeIndex.Encode for examples).
// Constructing the EncodedS2ShapeIndex uses additional memory, but even so
// the total memory usage immediately after construction is typically 25-35%
// of the corresponding MutableS2ShapeIndex size.
//
// Note that MutableS2ShapeIndex will still be faster and use less memory if
// you need to decode the entire index.  Similarly MutableS2ShapeIndex will be
// faster if you plan to execute a large number of operations on it.  The main
// advantage of EncodedS2ShapeIndex is that it is much faster and uses less
// memory when only a small portion of the data needs to be decoded.
//
// Example code showing how to create an encoded index:
//
//   Encoder encoder;
//   s2shapeutil.CompactEncodeTaggedShapes(index, encoder);
//   index.Encode(encoder);
//   string encoded(encoder.base(), encoder.length());  // Encoded data.
//
// Example code showing how to use an encoded index:
//
//   Decoder decoder(encoded.data(), encoded.size());
//   EncodedS2ShapeIndex index;
//   index.Init(&decoder, s2shapeutil.LazyDecodeShapeFactory(&decoder));
//   S2ClosestEdgeQuery query(&index);
//   S2ClosestEdgeQuery.PointTarget target(test_point);
//   if (query.IsDistanceLessOrEqual(&target, limit)) {
//     ...
//   }
//
// Note that EncodedS2ShapeIndex does not make a copy of the encoded data, and
// therefore the client must ensure that this data outlives the
// EncodedS2ShapeIndex object.
//
// There are a number of built-in classes that work with S2ShapeIndex objects.
// Generally these classes accept any collection of geometry that can be
// represented by an S2ShapeIndex, i.e. any combination of points, polylines,
// and polygons.  Such classes include:
//
// - S2ContainsPointQuery: returns the shape(s) that contain a given point.
//
// - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
//                       S2CellId, or S2ShapeIndex.
//
// - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
//
// - S2BooleanOperation: computes boolean operations such as union,
//                       and boolean predicates such as containment.
//
// - S2ShapeIndexRegion: can be used together with S2RegionCoverer to
//                       approximate geometry as a set of S2CellIds.
//
// - S2ShapeIndexBufferedRegion: computes approximations that have been
//                               expanded by a given radius.
//
// EncodedS2ShapeIndex is thread-compatible, meaning that const methods are
// thread safe, and non-const methods are not thread safe.  The only non-const
// method is Minimize(), so if you plan to call Minimize() while other threads
// are actively using the index that you must use an external reader-writer
// lock such as absl.Mutex to guard access to it.  (There is no global state
// and therefore each index can be guarded independently.)

namespace S2Geometry;

using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Runtime.InteropServices;
using Options = MutableS2ShapeIndex.Options;

public sealed class EncodedS2ShapeIndex : S2ShapeIndex, IDisposable
{
    public readonly Options Options_;

    private readonly ShapeFactory ShapeFactory_;

    // A vector containing all shapes in the index.  Initially all shapes are
    // set to kUndecodedShape(); as shapes are decoded, they are added to the
    // vector using atomic.compare_exchange_strong.
    private readonly ConcurrentDictionary<int, S2Shape?> Shapes;

    // A vector containing the S2CellIds of each cell in the index.
    private readonly EncodedS2CellIdVector CellIds;

    // A vector containing the encoded contents of each cell in the index.
    private readonly EncodedStringVector EncodedCells;

    // A raw array containing the decoded contents of each cell in the index.
    // Initially all values are *uninitialized memory*.  The cells_decoded_
    // field below keeps track of which elements are present.
    private readonly Lazy<S2ShapeIndexCell>[] Cells;

    // In order to minimize destructor time when very few cells of a large
    // S2ShapeIndex are needed, we keep track of the indices of the first few
    // cells to be decoded.  This lets us avoid scanning the cells_decoded_
    // vector when the number of cells decoded is very small.
    private readonly InputEdgeLoop CellCache;

    // Creates an index ~that must be initialized by calling Init()~.
    //
    // Initializes the EncodedS2ShapeIndex, returning true on success.
    //
    // This method does not decode the S2Shape objects in the index; this is
    // the responsibility of the client-provided function "shape_factory"
    // (see s2shapeutil_coding.h).  Example usage:
    //
    //   index.Init(decoder, S2ShapeUtil.LazyDecodeShapeFactory(decoder));
    //
    // Note that the encoded shape vector must *precede* the encoded S2ShapeIndex
    // in the Decoder's data buffer in this example.
    private EncodedS2ShapeIndex(int maxEdgesPerCell, ShapeFactory shapeFactory, EncodedS2CellIdVector cellIds,
                    Lazy<S2ShapeIndexCell>[] cells, EncodedStringVector encodedCells)
    {
        Minimize();
        Options_ = new()
        {
            MaxEdgesPerCell = maxEdgesPerCell
        };
        ShapeFactory_ = shapeFactory;
        CellIds = cellIds;
        Cells = cells;
        EncodedCells = encodedCells;
        Shapes = new();
        CellCache = [];
    }

    public static (bool, EncodedS2ShapeIndex?) Factory(Decoder decoder, ShapeFactory shape_factory)
    {
        if (!decoder.TryGetVarUInt64(out var max_edges_version)) return (false, null);
        int version = (int)(max_edges_version & 3);
        if (version != MutableS2ShapeIndex.kCurrentEncodingVersionNumber)
        {
            return (false, null);
        }
        var maxEdgesPerCell = (int)(max_edges_version >> 2);

        // AtomicShape is a subtype of atomic<S2Shape> that changes the
        // default constructor value to kUndecodedShape().  This saves the effort of
        // initializing all the elements twice.
        var shapeFactory = (ShapeFactory)shape_factory.CustomClone();
        var (success, cellIds) = EncodedS2CellIdVector.Init(decoder);
        if (!success) return (false, null);

        // The cells_ elements are *uninitialized memory*.  Instead we have bit
        // vector (cells_decoded_) to indicate which elements of cells_ are valid.
        // This reduces constructor times by more than a factor of 50, since rather
        // than needing to initialize one 64-bit pointer per cell to zero, we only
        // need to initialize one bit per cell to zero.
        //
        // For very large S2ShapeIndexes the internal memset() call to initialize
        // cells_decoded_ still takes about 1.3 microseconds per million cells
        // (assuming an optimized implementation that writes 32 bytes per cycle),
        // but this seems reasonable relative to other likely costs (I/O, etc).
        //
        // NOTE(ericv): DO NOT use make_unique<> here! make_unique<> allocates memory
        // using "new T[n]()", which initializes all elements of the array.  This
        // slows down some benchmarks by over 100x.
        //
        // cells_ = make_unique<S2ShapeIndexCell*>[](cell_ids_.size());
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        //                                NO NO NO
        var cells = new Lazy<S2ShapeIndexCell>[cellIds!.Count];

        var (success2, shape) = EncodedStringVector.Init(decoder);
        if (success2)
        {
            return (true, new(maxEdgesPerCell, shapeFactory, cellIds, cells, shape!));
        }
        return (false, null);
    }

    public void Dispose()
    {
        // Although Minimize() does slightly more than required for destruction
        // (i.e., it resets vector elements to their default values), this does not
        // affect benchmark times.
        Minimize();
    }

    // The number of distinct shape ids in the index.  This equals the number of
    // shapes in the index provided that no shapes have ever been removed.
    // (Shape ids are not reused.)
    public override int NumShapeIds() { return Shapes.Count; }

    // Return a pointer to the shape with the given id, or null if the shape
    // has been removed from the index.
    public override S2Shape? Shape(int id)
    {
        var shape = Shapes[id];
        if (shape is not null) return shape;
        return GetShape(id);
    }

    // Minimizes memory usage by requesting that any data structures that can be
    // rebuilt should be discarded.  This method invalidates all iterators.
    //
    // Like all non-methods, this method is not thread-safe.
    public override void Minimize()
    {
        if (Cells is null) return;  // Not initialized yet.

        foreach (var atomic_shape in Shapes)
        {
            if (atomic_shape.Value is not null)
                Shapes[atomic_shape.Key] = null;
        }

        if (CellCache.Count < MaxCellCacheSize())
        {
            // When only a tiny fraction of the cells are decoded, we keep track of
            // those cells in cell_cache_ to avoid the cost of scanning the
            // cells_decoded_ vector.  (The cost is only about 1 cycle per 64 cells,
            // but for a huge polygon with 1 million cells that's still 16000 cycles.)
            foreach (int pos in CellCache)
                Cells[pos] = new Lazy<S2ShapeIndexCell>();
        }
        else
        {
            for (var i = 0; i < Cells.Length; i++)
            {
                if (Cells[i].IsValueCreated)
                {
                    Cells[i] = new Lazy<S2ShapeIndexCell>();
                }
            }
        }
        CellCache.Clear();
    }

    public override (int pos, bool found) SeekCell(S2CellId target)
    {
        var pos = CellIds.LowerBound(target);

        if (pos < 0) return (~pos, false);
        return (pos, true);
    }

    // Returns the number of bytes currently occupied by the index (including any
    // unused space at the end of vectors, etc). It has the same thread safety
    // as the other "const" methods (see introduction).
    public override int SpaceUsed()
    {
        // TODO(ericv): Add SpaceUsed() method to S2Shape base class,and include
        // memory owned by the allocated S2Shapes (here and in S2ShapeIndex).
        int size = SizeHelper.SizeOf(this);
        size += Shapes.Count * SizeHelper.SizeOf<S2Shape>();
        size += CellIds.Count * Marshal.SizeOf<S2ShapeIndexCell>();  // cells_
        size += CellCache.Count * sizeof(int);
        return size;
    }

    private S2Shape? GetShape(int id)
    {
        // This method is called when a shape has not been decoded yet.
        var shape = ShapeFactory_[id];
        shape?.SetId(id);
        return Shapes[id];
    }

    public override S2CellId? GetCellId(int index) => CellIds[index];

    public override S2ShapeIndexCell? GetCell(int index)
    {
        // memory_order_release ensures that no reads or writes in the current
        // thread can be reordered after this store, and all writes in the current
        // thread are visible to other threads that acquire the same atomic
        // variable.
        //
        // memory_order_acquire ensures that no reads or writes in the current
        // thread can be reordered before this load, and all writes in other threads
        // that release the same atomic variable are visible in this thread.
        //
        // We use this to implement lock-free synchronization on the read path as
        // follows:
        //
        //  1. cells_decoded(i) is updated using acquire/release semantics
        //  2. cells_[i] is written before cells_decoded(i)
        //  3. cells_[i] is read after cells_decoded(i)
        //
        // Note that we do still use a lock for the write path to ensure that
        // cells_[i] and cell_decoded(i) are updated together atomically.
        if (Cells[index].IsValueCreated)
        {
            var cell2 = Cells[index].Value;
            if (cell2 is not null) return cell2;
        }

        // Decode the cell before acquiring the spinlock in order to minimize the
        // time that the lock is held.
        S2ShapeIndexCell cell = new();
        var decoder = EncodedCells.GetDecoder(index);
        if (!cell.Decode(NumShapeIds(), decoder))
        {
            return null;
        }
        // Recheck cell_decoded(i) once we hold the lock in case another thread
        // has decoded this cell in the meantime.
        if (Cells[index].IsValueCreated)
        {
            var cell2 = Cells[index].Value;
            if (cell2 is not null) return cell2;
        }

        // Update the cell, setting cells_[i] before cell_decoded(i).
        Cells[index] = new Lazy<S2ShapeIndexCell>(cell);
        if (CellCache.Count < MaxCellCacheSize())
        {
            CellCache.Add(index);
        }
        return cell;
    }

    private int MaxCellCacheSize()
    {
        // The cell cache is sized so that scanning decoded_cells_ in the destructor
        // costs about 30 cycles per decoded cell in the worst case.  (This overhead
        // is acceptable relative to the other costs of decoding each cell.)
        //
        // For example, if there are 65,536 cells then we won't need to scan
        // encoded_cells_ unless we decode at least (65536/2048) == 32 cells.  It
        // takes about 1 cycle per 64 cells to scan encoded_cells_, so that works
        // out to (65536/64) == 1024 cycles.  However this cost is amortized over
        // the 32 cells decoded, which works out to 32 cycles per cell.
        return CellIds.Count >> 11;
    }

    public override EnumeratorBase<S2ShapeIndexCell> GetNewEnumerator(InitialPosition pos) =>
        new Enumerator(this, pos);

    private new sealed class Enumerator : EnumeratorBase<S2ShapeIndexCell>, IEnumerator<S2ShapeIndexCell>
    {
        private readonly EncodedS2ShapeIndex _index;
        private readonly Func<int> NumCells;
        private int position;

        public Enumerator(EncodedS2ShapeIndex index, InitialPosition pos)
        {
            _index = index;
            NumCells = () => index.CellIds.Count; 
            position = pos == InitialPosition.BEGIN ? -1 : NumCells();
            Refresh();
        }

        public override S2ShapeIndexCell Current => Cell!;

        object IEnumerator.Current => Current;

        public override void Finish()
        {
            position = NumCells();
            Refresh();
        }

        public override void Dispose() { }

        public override bool MoveNext()
        {
            MyDebug.Assert(!Done());
            ++position;
            Refresh();
            return !Done();
        }
        public override bool MovePrevious()
        {
            if (position == 0) return false;
            --position;
            Refresh();
            return true;
        }
        public override void Reset() => position = -1;
        public override bool Done() => position >= NumCells();
        public override void SetPosition(int pos) => position = pos;

        public override bool Locate(S2Point target) => LocateImpl(this, target);
        public override S2CellRelation Locate(S2CellId target) => LocateImpl(this, target);
        private void Refresh()
        {
            if (position == NumCells())
            {
                SetFinished();
            }
            else
            {
                // It's faster to initialize the cell to nullptr even if it has already
                // been decoded, since algorithms frequently don't need it (i.e., based on
                // the S2CellId they might not need to look at the cell contents).
                SetState(_index.CellIds[position], null);
            }
        }

        public override void Seek(S2CellId target)
        {
            position = _index.CellIds.LowerBound(target);
            Refresh();
        }
        protected override S2ShapeIndexCell? GetCell() => _index.GetCell(position);
    }
}
