using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace S2Geometry
{
    using Options = MutableS2ShapeIndex.Options;
    using S2ShapeIndexIdCell = KeyData<S2CellId, S2ShapeIndexCell>;

    public sealed class EncodedS2ShapeIndex : S2ShapeIndex, IDisposable
    {
        // Creates an index that must be initialized by calling Init().
        public EncodedS2ShapeIndex() { }

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
        public bool Init(Decoder decoder, ShapeFactory shape_factory)
        {
            Minimize();
            if (!decoder.TryGetVarUInt64(out var max_edges_version)) return false;
            int version = (int)(max_edges_version & 3);
            if (version != MutableS2ShapeIndex.kCurrentEncodingVersionNumber)
            {
                return false;
            }
            Options_.MaxEdgesPerCell = (int)(max_edges_version >> 2);

            // AtomicShape is a subtype of atomic<S2Shape that changes the
            // default constructor value to kUndecodedShape().  This saves the effort of
            // initializing all the elements twice.
            shapes_ = new();
            shape_factory_ = (ShapeFactory)shape_factory.Clone();
            if (!cell_ids_.Init(decoder)) return false;

            // The cells_ elements are *uninitialized memory*.  Instead we have bit
            // vector (cells_decoded_) to indicate which elements of cells_ are valid.
            // This reduces constructor times by more than a factor of 50, since rather
            // than needing to initialize one 64-bit pointer per cell to zero, we only
            // need to initialize one bit per cell to zero.
            //
            // For very large S2ShapeIndexes the internal memset() call to initialize
            // cells_decoded_ still takes about 4 microseconds per million cells, but
            // this seems reasonable relative to other likely costs (I/O, etc).
            //
            // NOTE(ericv): DO NOT use make_unique<> here! make_unique<> allocates memory
            // using "new T[n]()", which initializes all elements of the array.  This
            // slows down some benchmarks by over 100x.
            //
            // cells_ = make_unique<std::atomic<S2ShapeIndexCell*>[]>(cell_ids_.size());
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            //                                NO NO NO
            cells_ = new Lazy<S2ShapeIndexCell>[cell_ids_.Count];

            return EncodedStringVector.Init(decoder, out encoded_cells_);
        }

        public void Dispose()
        {
            // Although Minimize() does slightly more than required for destruction
            // (i.e., it resets vector elements to their default values), this does not
            // affect benchmark times.
            Minimize();
        }

        public Options Options_ { get; }

        // The number of distinct shape ids in the index.  This equals the number of
        // shapes in the index provided that no shapes have ever been removed.
        // (Shape ids are not reused.)
        public override int NumShapeIds() { return shapes_.Count; }

        // Return a pointer to the shape with the given id, or null if the shape
        // has been removed from the index.
        public override S2Shape Shape(int id)
        {
            var shape = shapes_[id];
            if (shape != null) return shape;
            return GetShape(id);
        }

        // Minimizes memory usage by requesting that any data structures that can be
        // rebuilt should be discarded.  This method invalidates all iterators.
        //
        // Like all non-methods, this method is not thread-safe.
        public override void Minimize()
        {
            if (cells_ == null) return;  // Not initialized yet.

            foreach (var atomic_shape in shapes_)
            {
                if (atomic_shape.Value != null)
                    shapes_[atomic_shape.Key] = null;
            }

            if (cell_cache_.Count < MaxCellCacheSize())
            {
                // When only a tiny fraction of the cells are decoded, we keep track of
                // those cells in cell_cache_ to avoid the cost of scanning the
                // cells_decoded_ vector.  (The cost is only about 1 cycle per 64 cells,
                // but for a huge polygon with 1 million cells that's still 16000 cycles.)
                foreach (int pos in cell_cache_)
                    cells_[pos] = null;
            }
            else
            {
                for (var i = 0; i < cells_.Length; i++)
                {
                    if (cells_[i].IsValueCreated)
                    {
                        cells_[i] = null;
                    }
                }
            }
            cell_cache_.Clear();
        }

        public override (int pos, bool found) SeekCell(S2CellId target)
        {
            var pos = cell_ids_.LowerBound(target);

            if (pos < 0) return (~pos, false);
            return (pos, true);
        }

        // Returns the number of bytes currently occupied by the index (including any
        // unused space at the end of vectors, etc). It has the same thread safety
        // as the other "const" methods (see introduction).
        public override int SpaceUsed()
        {
            // TODO(ericv): Add SpaceUsed() method to S2Shape base class,and Include
            // memory owned by the allocated S2Shapes (here and in S2ShapeIndex).
            int size = Marshal.SizeOf(this);
            size += shapes_.Count * Marshal.SizeOf(typeof(S2Shape));
            size += cell_ids_.Count * Marshal.SizeOf(typeof(S2ShapeIndexCell));  // cells_
            size += cell_cache_.Count * sizeof(int);
            return size;
        }

        private S2Shape GetShape(int id)
        {
            // This method is called when a shape has not been decoded yet.
            var shape = shape_factory_[id];
            if (shape != null) shape.SetId(id);
            return shapes_[id];
        }

        public override S2CellId GetCellId(int index) => cell_ids_[index];
        public override S2ShapeIndexIdCell? GetIndexCell(int index)
        {
            if (cells_[index].IsValueCreated)
            {
                var cell2 = cells_[index].Value;
                if (cell2 != null) return new(GetCellId(index), cell2);
            }
            S2ShapeIndexCell cell = new();
            var decoder = encoded_cells_.GetDecoder(index);
            if (!cell.Decode(NumShapeIds(), decoder))
            {
                return null;
            }
            if (cell_cache_.Count < MaxCellCacheSize())
            {
                cell_cache_.Add(index);
            }
            return new(GetCellId(index), cell);
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
            return cell_ids_.Count >> 11;
        }

        public override IReversableEnumerator<S2ShapeIndexIdCell> GetNewEnumerator()
            => new EncodedS2ShapeIndexEnumerator(this);
        public override IEnumerable<S2ShapeIndexIdCell> GetNewEnumerable() 
        {
            var enumerator = GetNewEnumerator();
            while (enumerator.MoveNext())
                yield return enumerator.Current;
        }
        public override int GetEnumerableCount() => cell_ids_.Count;

        private ShapeFactory shape_factory_;

        // A vector containing all shapes in the index.  Initially all shapes are
        // set to kUndecodedShape(); as shapes are decoded, they are added to the
        // vector using atomic.compare_exchange_strong.
        private ConcurrentDictionary<int, S2Shape> shapes_;

        // A vector containing the S2CellIds of each cell in the index.
        private readonly EncodedS2CellIdVector cell_ids_ = new();

        // A vector containing the encoded contents of each cell in the index.
        private EncodedStringVector encoded_cells_;

        // A raw array containing the decoded contents of each cell in the index.
        // Initially all values are *uninitialized memory*.  The cells_decoded_
        // field below keeps track of which elements are present.
        private Lazy<S2ShapeIndexCell>[] cells_;

        // In order to minimize destructor time when very few cells of a large
        // S2ShapeIndex are needed, we keep track of the indices of the first few
        // cells to be decoded.  This lets us avoid scanning the cells_decoded_
        // vector when the number of cells decoded is very small.
        private readonly List<int> cell_cache_ = new();

        private class EncodedS2ShapeIndexEnumerator : IReversableEnumerator<S2ShapeIndexIdCell>
        {
            private readonly EncodedS2CellIdVector _cellIds;
            private int position;
            public EncodedS2ShapeIndexEnumerator(EncodedS2ShapeIndex index)
            { _cellIds = index.cell_ids_; position = -1; }

            public S2ShapeIndexIdCell Current => new S2ShapeIndexIdCell(_cellIds[position], null);

            object IEnumerator.Current => Current;

            public void Dispose() { }
            public bool MoveNext() => ++position >= 0 && position < _cellIds.Count;
            public bool MovePrevious() => --position >= 0 && position < _cellIds.Count;
            public void Reset() => position = -1;
            public bool Done() => position >= _cellIds.Count;
            public void SetPosition(int pos) => position = pos;
        }
    }
}
