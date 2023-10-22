// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The main documentation is with the class definition below.
// (Some helper classes are defined first.)
//
// S2ShapeIndex has two major subtypes:
//
//  - MutableS2ShapeIndex is for building new S2ShapeIndexes.  Indexes may
//    also be updated dynamically by inserting or deleting shapes.  Once an
//    index has been built it can be encoded compactly and later decoded as
//    either a MutableS2ShapeIndex or an EncodedS2ShapeIndex.
//
//  - EncodedS2ShapeIndex is an S2ShapeIndex implementation that works
//    directly with encoded data.  Rather than decoding everything in advance,
//    geometry is decoded incrementally (down to individual edges) as needed.
//    It can be initialized from a single block of data in nearly constant
//    time.  This saves large amounts of memory and is also much faster in the
//    common situation where geometric data is loaded from somewhere, decoded,
//    and then only a single operation is performed on it.  (Speedups for
//    large geometric objects can be over 1000x.)  It supports all
//    S2ShapeIndex operations including boolean operations, measuring
//    distances, etc.


// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The objects in the index are known as "shapes", and may consist of
// points, polylines, and/or polygons, possibly overlapping.  The index makes
// it very fast to answer queries such as finding nearby shapes, measuring
// distances, testing for intersection and containment, etc.
//
// Each object in the index implements the S2Shape interface.  An S2Shape is a
// collection of edges that optionally defines an interior.  The edges do not
// need to be connected, so for example an S2Shape can represent a polygon
// with multiple shells and/or holes, or a set of polylines, or a set of
// points.  All geometry within a single S2Shape must have the same dimension,
// so for example if you want to create an S2ShapeIndex containing a polyline
// and 10 points, then you will need at least two different S2Shape objects.
//
// There are two important types of S2ShapeIndex.  MutableS2ShapeIndex allows
// you to build an index incrementally by adding or removing shapes, whereas
// EncodedS2ShapeIndex works very efficiently with existing indexes by keeping
// the index data in encoded form (see introduction at the top of this file).
//
// Code that only needs read-only ("const") access to an index should use the
// S2ShapeIndex base class as the parameter type, so that it will work with
// any S2ShapeIndex subtype.  For example:
//
//   void DoSomething(S2ShapeIndex& index) {
//     ... works with MutableS2ShapeIndex or EncodedS2ShapeIndex ...
//   }
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
// Here is an example showing how to index a set of polygons and then
// determine which polygon(s) contain each of a set of query points:
//
//   void TestContainment(S2Point[] points,
//                        S2Polygon[] polygons) {
//     MutableS2ShapeIndex index;
//     foreach (var polygon in polygons) {
//       index.Add(absl.make_unique<S2Polygon.Shape>(polygon));
//     }
//     var query = index.MakeS2ContainsPointQuery();
//     foreach (var& point : points) {
//       for (S2Shape* shape in query.GetContainingShapes(point)) {
//         S2Polygon* polygon = polygons[shape.id()];
//         ... do something with (point, polygon) ...
//       }
//     }
//   }
//
// This example uses S2Polygon.Shape, which is one example of an S2Shape
// object.  S2Polyline and S2Loop also have nested Shape classes, and there are
// additional S2Shape types defined in *_shape.h.
//
// Internally, an S2ShapeIndex is essentially a map from S2CellIds to the set
// of shapes that intersect each S2CellId.  It is adaptively refined to ensure
// that no cell contains more than a small number of edges.
//
// In addition to implementing a shared set of virtual methods, all
// S2ShapeIndex subtypes define an Iterator type with the same API.  This
// makes it easy to convert code that uses a particular S2ShapeIndex subtype
// to instead use the abstract base class (or vice versa).  You can also
// choose to avoid the overhead of virtual method calls by making the
// S2ShapeIndex type a template argument, like this:
//
//   template <class IndexType>
//   void DoSomething(IndexType& index) {
//     for (typename IndexType.Iterator it(out index, S2ShapeIndex.InitialPosition.BEGIN);
//          !it.done(); it.Next()) {
//       ...
//     }
//   }
//
// The S2ShapeIndex subtypes provided by the S2 library are thread-compatible,
// meaning that const methods are const methods may be called concurrently
// from multiple threads, while non-const methods require exclusive access to
// the S2ShapeIndex.

namespace S2Geometry;

using System.Collections;
using System.Runtime.CompilerServices;
using static S2Geometry.S2ShapeIndex;

public abstract class S2ShapeIndex : IEnumerable<S2Shape>
{
    // Returns the number of distinct shape ids in the index.  This is the same
    // as the number of shapes provided that no shapes have ever been removed.
    // (Shape ids are never reused.)
    public abstract int NumShapeIds();

    // Returns a pointer to the shape with the given id, or null if the shape
    // has been removed from the index.
    public abstract S2Shape? Shape(int id);

    // Returns the number of bytes currently occupied by the index (including any
    // unused space at the end of vectors, etc).
    public abstract int SpaceUsed();

    // Minimizes memory usage by requesting that any data structures that can be
    // rebuilt should be discarded.  This method invalidates all iterators.
    //
    // Like all non-methods, this method is not thread-safe.
    public abstract void Minimize();

    // When passed to an Iterator constructor, specifies whether the iterator
    // should be positioned at the beginning of the index (BEGIN), the end of
    // the index (END), or arbitrarily (UNPOSITIONED).  By default iterators are
    // unpositioned, since this avoids an extra seek in this situation where one
    // of the seek methods (such as Locate) is immediately called.
    public enum InitialPosition { BEGIN, END, UNPOSITIONED };

    // ShapeFactory is an interface for decoding vectors of S2Shapes.  It allows
    // random access to the shapes in order to support lazy decoding.  See
    // s2shapeutil_coding.h for useful subtypes.
    public abstract class ShapeFactory : ICustomCloneable
    {
        // Returns the number of S2Shapes in the vector.
        public abstract int Count { get; }

        // Returns the S2Shape object corresponding to the given "shape_id".
        // Returns null if a shape cannot be decoded or a shape is missing
        // (e.g., because MutableS2ShapeIndex.Release() was called).
        public abstract S2Shape? this[int shape_id] { get; }

        // Returns a deep copy of this ShapeFactory.
        public abstract object CustomClone();
    }

    public virtual IEnumerator<S2Shape> GetEnumerator()
    {
        for (var i = 0; i < NumShapeIds(); i++)
        {
            yield return Shape(i)!;
        }
    }

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public abstract EnumeratorBase<S2ShapeIndexCell> GetNewEnumerator(InitialPosition pos);

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and
    // returns INDEXED.  Otherwise if T contains one or more (smaller) index
    // cells, it positions the iterator at the first such cell I and returns
    // SUBDIVIDED.  Otherwise it returns DISJOINT.
    public (S2CellRelation cellRelation, int pos) LocateCell(S2CellId target)
    {
        // Let T be the target, let I = cell_map_.GetLowerBound(T.RangeMin), and
        // let I' be the predecessor of I.  If T contains any index cells, then T
        // contains I.  Similarly, if T is contained by an index cell, then the
        // containing cell is either I or I'.  We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.

        var (pos, _) = SeekCell(target.RangeMin());
        var id = GetCellId(pos);
        if (id is not null)
        {
            if (id >= target && id.Value.RangeMin() <= target) return (S2CellRelation.INDEXED, pos);
            if (id <= target.RangeMax()) return (S2CellRelation.SUBDIVIDED, pos);
        }
        if (--pos >= 0)
        {
            id = GetCellId(pos);

            if (id is not null && id.Value.RangeMax() >= target)
                return (S2CellRelation.INDEXED, pos);
        }
        return (S2CellRelation.DISJOINT, 0);
    }

    // The default implementation of SearchPointPos(S2Point).
    public virtual (int pos, bool found) LocatePoint(S2Point target_point)
    {
        // Let I = cell_map_.lower_bound(T), where T is the leaf cell containing
        // "target_point".  Then if T is contained by an index cell, then the
        // containing cell is either I or I'.  We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.

        var target = new S2CellId(target_point);
        var (pos, _) = SeekCell(target);
        var cell1 = GetCellId(pos);
        if (cell1 is not null && cell1.Value.RangeMin() <= target) return (pos, true);
        var cell2 = GetCellId(pos - 1);
        if (cell2 is not null && cell2.Value.RangeMax() >= target) return (pos - 1, true);
        return (pos - 1, false);
    }

    // Returns 1) the zero-based index of item in the sorted Container,
    // if item is found;  otherwise, the index of the next element that is
    // larger than item or, if there is no larger element, List.Count.
    // 2) A boolean indicating if the item was found.
    public abstract (int pos, bool found) SeekCell(S2CellId target);

    public abstract S2CellId? GetCellId(int pos);
    public abstract S2ShapeIndexCell? GetCell(int pos);

    // A random access iterator that provides low-level access to the cells of
    // the index.  Cells are sorted in increasing order of S2CellId.
    public sealed class Enumerator(S2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED) : S2CellEnumerator<S2ShapeIndexCell>, IEnumerator<S2ShapeIndexCell>
    {
        private readonly EnumeratorBase<S2ShapeIndexCell> Enumerator_ = index.GetNewEnumerator(pos);

        // Returns the S2CellId of the current index cell.  If done() is true,
        // returns a value larger than any valid S2CellId (S2CellId.Sentinel()).
        public override S2CellId Id => Enumerator_.Id;

        // Returns a reference to the contents of the current index cell.
        // REQUIRES: !done()
        public S2ShapeIndexCell Cell => Enumerator_.Cell!;

        public override S2ShapeIndexCell Current => Enumerator_.Current;

        object IEnumerator.Current => Current!;

        // Positions the iterator at the next index cell.
        // REQUIRES: !done()
        public override bool MoveNext() => Enumerator_.MoveNext();

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        public override bool MovePrevious() => Enumerator_.MovePrevious();
        public override void Reset() => Enumerator_.Reset();
        public override bool Done() => Enumerator_.Done();

        // Positions the iterator past the last index cell.
        public override void Finish() => Enumerator_.Finish();
        public override void SetPosition(int position) => Enumerator_.SetPosition(position);
        public override void Dispose() { GC.SuppressFinalize(this); }

        // Positions the iterator at the first cell with id() >= target, or at the
        // end of the index if no such cell exists.
        public override void Seek(S2CellId target) => Enumerator_.Seek(target);

        // Positions the iterator at the cell containing target and returns true. If
        // no such cell exists, return false and leave the iterator in an undefined
        // (but valid) state.
        public override bool Locate(S2Point target) => LocateImpl(this, target);

        // Let T be the target S2CellId.  If T is contained by some index cell I
        // (including equality), this method positions the iterator at I and
        // returns INDEXED.  Otherwise if T contains one or more (smaller) index
        // cells, it positions the iterator at the first such cell I and returns
        // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
        // positioned arbitrarily.
        public override S2CellRelation Locate(S2CellId target) => LocateImpl(this, target);
    }

    // Each subtype of S2ShapeIndex should define an Iterator type derived
    // from the following base class.
    public abstract class EnumeratorBase<T> : S2CellEnumerator<T>
    {
        // Returns the S2CellId of the current index cell.  If done() is true,
        // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
        public override S2CellId Id { get => _Id; }
        private S2CellId _Id;

        // Returns the center point of the cell.
        // REQUIRES: !done()
        public S2Point Center
        {
            get
            {
                MyDebug.Assert(!Done());
                return _Id.ToPoint();
            }
        }

        // Returns a reference to the contents of the current index cell.
        // REQUIRES: !done()
        public virtual S2ShapeIndexCell? Cell
        {
            get
            {
                // Like other const methods, this method is thread-safe provided that it
                // does not overlap with calls to non-const methods.
                MyDebug.Assert(!Done());
                S2ShapeIndexCell? cell;
                lock (this)
                {
                    cell = CellInternal;
                }
                if (cell is null)
                {
                    cell = GetCell();
                    Cell = cell;
                }
                return cell;
            }
            private set
            {
                lock (this)
                {
                    CellInternal = value;
                }
            }
        }
        protected S2ShapeIndexCell? CellInternal { get; private set; }

        protected EnumeratorBase() { _Id = S2CellId.Sentinel; CellInternal = null; }

        // Returns true if the iterator is positioned past the last index cell.
        public override bool Done() => _Id == S2CellId.Sentinel;

        // Sets the iterator state.  "cell" typically points to the cell contents,
        // but may also be given as "nullptr" in order to implement decoding on
        // demand.  In that situation, the first that the client attempts to
        // access the cell contents, the GetCell() method is called and "cell_" is
        // updated in a thread-safe way.
        protected void SetState(S2CellId id, S2ShapeIndexCell? cell)
        {
            _Id = id;
            Cell = cell;
        }

        // Sets the iterator state so that done() is true.
        public void SetFinished()
        {
            SetState(S2CellId.Sentinel, null);
        }

        // This method is called to decode the contents of the current cell, if
        // set_state() was previously called with a nullptr "cell" argument.  This
        // allows decoding on demand for subtypes that keep the cell contents in
        // an encoded state.  It does not need to be implemented at all if
        // set_state() is always called with (cell != nullptr).
        //
        // REQUIRES: This method is thread-safe.
        // REQUIRES: Multiple calls to this method return the same value.
        protected abstract S2ShapeIndexCell? GetCell();
    }
}

public static class EnumeratorBaseExtensions
{
    public static IEnumerable GetEnumerator<T>(this EnumeratorBase<T> me)
    {
        while (me.MoveNext())
        {
            yield return me.Current;
        }
    }
}

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
public class S2ShapeIndexCell
{
    public S2ShapeIndexCell() { }

    // Returns the number of clipped shapes in this cell.
    public int NumClipped() => shapes_.Count;

    // Returns the clipped shape at the given index.  Shapes are kept sorted in
    // increasing order of shape id.
    //
    // REQUIRES: 0 <= i < num_clipped()
    public S2ClippedShape Clipped(int i) { return shapes_[i]; }

    // Returns a pointer to the clipped shape corresponding to the given shape,
    // or null if the shape does not intersect this cell.
    public S2ClippedShape? FindClipped(S2Shape shape) => FindClipped(shape.Id);
    public S2ClippedShape? FindClipped(int shape_id)
    {
        // Linear search is fine because the number of shapes per cell is typically
        // very small (most often 1), and is large only for pathological inputs
        // (e.g. very deeply nested loops).
        foreach (var s in shapes_)
        {
            if (s.ShapeId == shape_id) return s;
        }
        return null;
    }

    // Convenience method that returns the total number of edges in all clipped
    // shapes.
    public int NumEdges()
    {
        int n = 0;
        for (int i = 0; i < NumClipped(); ++i) n += Clipped(i).NumEdges;
        return n;
    }

    // Appends an encoded representation of the S2ShapeIndexCell to "encoder".
    // "num_shape_ids" should be set to index.num_shape_ids(); this information
    // allows the encoding to be more compact in some cases.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(int num_shape_ids, Encoder encoder)
    {
        // The encoding is designed to be especially compact in certain common
        // situations:
        //
        // 1. The S2ShapeIndex contains exactly one shape.
        //
        // 2. The S2ShapeIndex contains more than one shape, but a particular index
        //    cell contains only one shape (num_clipped == 1).
        //
        // 3. The edge ids for a given shape in a cell form a contiguous range.
        //
        // The details were optimized by constructing an S2ShapeIndex for each
        // feature in Google's geographic repository and measuring their total
        // encoded size.  The MutableS2ShapeIndex encoding (of which this function
        // is just one part) uses an average of 1.88 bytes per vertex for features
        // consisting of polygons or polylines.
        //
        // Note that this code does not bother handling num_shapes >= 2**28 or
        // num_edges >= 2**29.  This could be fixed using varint64 in a few more
        // places, but if a single cell contains this many shapes or edges then we
        // have bigger problems than just the encoding format :)
        if (num_shape_ids == 1)
        {
            // If the entire S2ShapeIndex contains just one shape, then we don't need
            // to encode any shape ids.  This is a very important and common case.
            MyDebug.Assert(NumClipped() == 1, "Index invariant: no empty cells.");
            S2ClippedShape clipped = Clipped(0);
            MyDebug.Assert(clipped.ShapeId == 0);
            int n = clipped.NumEdges;
            encoder.Ensure(Encoder.kVarintMax64 + n * Encoder.kVarintMax32);
            var ccc = clipped.ContainsCenter ? 2 : 0;
            if (n >= 2 && n <= 17 && clipped.Edges[n - 1] - clipped.Edges[0] == n - 1)
            {
                // The cell contains a contiguous range of edges (*most common case*).
                // If the starting edge id is small then we can encode the cell in one
                // byte.  (The n == 0 and n == 1 cases are encoded compactly below.)
                // This encoding uses a 1-bit tag because it is by far the most common.
                //
                // Encoding: bit 0: 0
                //           bit 1: contains_center
                //           bits 2-5: (num_edges - 2)
                //           bits 6+: edge_id
                encoder.PutVarInt64(clipped.Edges[0] << 6 | (n - 2) << 2 | ccc | 0);
            }
            else if (n == 1)
            {
                // The cell contains only one edge.  For edge ids up to 15, we can
                // encode the cell in a single byte.
                //
                // Encoding: bits 0-1: 1
                //           bit 2: contains_center
                //           bits 3+: edge_id
                encoder.PutVarInt64(clipped.Edges[0] << 3 | ccc | 1);
            }
            else
            {
                // General case (including n == 0, which is encoded compactly here).
                //
                // Encoding: bits 0-1: 3
                //           bit 2: contains_center
                //           bits 3+: num_edges
                encoder.PutVarInt64(n << 3 | ccc | 3);
                EncodeEdges(clipped, encoder);
            }
        }
        else
        {
            if (NumClipped() > 1)
            {
                // The cell contains more than one shape.  The tag for this encoding
                // must be distinguishable from the cases encoded below.  We can afford
                // to use a 3-bit tag because num_clipped is generally small.
                encoder.Ensure(Encoder.kVarintMax32);
                encoder.PutVarInt32((NumClipped() << 3) | 3);
            }
            // The shape ids are delta-encoded.
            int shape_id_base = 0;
            for (int j = 0; j < NumClipped(); ++j)
            {
                var clipped = Clipped(j);
                MyDebug.Assert(clipped.ShapeId >= shape_id_base);
                int shape_delta = clipped.ShapeId - shape_id_base;
                shape_id_base = clipped.ShapeId + 1;

                // Like the code above except that we also need to encode shape_id(s).
                // Because of this some choices are slightly different.
                int n = clipped.NumEdges;
                encoder.Ensure((n + 2) * Encoder.kVarintMax32);
                var ccc = clipped.ContainsCenter ? 2 : 0;
                if (n >= 1 && n <= 16 && clipped.Edges[n - 1] - clipped.Edges[0] == n - 1)
                {
                    // The clipped shape has a contiguous range of up to 16 edges.  This
                    // encoding uses a 1-bit tag because it is by far the most common.
                    //
                    // Encoding: bit 0: 0
                    //           bit 1: contains_center
                    //           bits 2+: edge_id
                    // Next value: bits 0-3: (num_edges - 1)
                    //             bits 4+: shape_delta
                    encoder.PutVarInt32(clipped.Edges[0] << 2 | ccc | 0);
                    encoder.PutVarInt32(shape_delta << 4 | (n - 1));
                }
                else if (n == 0)
                {
                    // Special encoding for clipped shapes with no edges.  Such shapes are
                    // common in polygon interiors.  This encoding uses a 3-bit tag in
                    // order to leave more bits available for the other encodings.
                    //
                    // NOTE(ericv): When num_clipped > 1, this tag could be 2 bits
                    // (because the tag used to indicate num_clipped > 1 can't appear).
                    // Alternatively, that tag can be considered reserved for future use.
                    //
                    // Encoding: bits 0-2: 7
                    //           bit 3: contains_center
                    //           bits 4+: shape_delta
                    encoder.PutVarInt32(shape_delta << 4 | ccc << 2 | 7);
                }
                else
                {
                    // General case.  This encoding uses a 2-bit tag, and the first value
                    // typically is encoded into one byte.
                    //
                    // Encoding: bits 0-1: 1
                    //           bit 2: contains_center
                    //           bits 3+: (num_edges - 1)
                    // Next value: shape_delta
                    encoder.PutVarInt32((n - 1) << 3 | ccc << 1 | 1);
                    encoder.PutVarInt32(shape_delta);
                    EncodeEdges(clipped, encoder);
                }
            }
        }
    }

    // Decodes an S2ShapeIndexCell, returning true on success.
    // "num_shape_ids" should be set to index.num_shape_ids().
    public bool Decode(int num_shape_ids, Decoder decoder)
    {
        // This function inverts the encodings documented above.
        if (num_shape_ids == 1)
        {
            // Entire S2ShapeIndex contains only one shape.
            if (!decoder.TryGetVarUInt64(out var header)) return false;
            if ((header & 1) == 0)
            {
                // The cell contains a contiguous range of edges.
                int num_edges = (int)(((header >> 2) & 15) + 2);
                S2ClippedShape clippedTmp = new(0 /*shape_id*/, num_edges, (header & 2) != 0);
                for (int i = 0, edge_id = (int)(header >> 6); i < num_edges; ++i)
                {
                    clippedTmp.Edges[i] = edge_id + i;
                }
                shapes_.Add(clippedTmp);
                return true;
            }
            if ((header & 2) == 0)
            {
                // The cell contains a single edge.
                S2ClippedShape clippedTmp = new(0 /*shape_id*/, 1 /*num_edges*/, (header & 4) != 0);
                clippedTmp.Edges[0] = (int)(header >> 3);
                shapes_.Add(clippedTmp);
                return true;
            }
            else
            {
                // The cell contains some other combination of edges.
                int num_edges = (int)(header >> 3);
                S2ClippedShape clippedTmp = new(0 /*shape_id*/, num_edges, (header & 4) != 0);
                var res = DecodeEdges(num_edges, clippedTmp, decoder);
                shapes_.Add(clippedTmp);
                return res;
            }
        }
        // S2ShapeIndex contains more than one shape.
        if (!decoder.TryGetVarUInt32(out var header32)) return false;
        int num_clipped = 1;
        if ((header32 & 7) == 3)
        {
            // This cell contains more than one shape.
            num_clipped = (int)(header32 >> 3);
            if (!decoder.TryGetVarUInt32(out header32)) return false;
        }
        int shape_id = 0;
        for (int j = 0; j < num_clipped; ++j, ++shape_id)
        {
            if (j > 0 && !decoder.TryGetVarUInt32(out header32)) return false;
            if ((header32 & 1) == 0)
            {
                // The clipped shape contains a contiguous range of edges.
                if (!decoder.TryGetVarUInt32(out var shape_id_count)) return false;
                shape_id += (int)(shape_id_count >> 4);
                int num_edges = (int)((shape_id_count & 15) + 1);
                S2ClippedShape clipped = new(shape_id, num_edges, (header32 & 2) != 0);
                for (int i = 0, edge_id = (int)(header32 >> 2); i < num_edges; ++i)
                {
                    clipped.Edges[i] = edge_id + i;
                }
                shapes_.Add(clipped);
            }
            else if ((header32 & 7) == 7)
            {
                // The clipped shape has no edges.
                shape_id += (int)(header32 >> 4);
                S2ClippedShape clipped = new(shape_id, 0, (header32 & 8) != 0);
                shapes_.Add(clipped);
            }
            else
            {
                // The clipped shape contains some other combination of edges.
                MyDebug.Assert((header32 & 3U) == 1U);
                if (!decoder.TryGetVarUInt32(out var shape_delta)) return false;
                shape_id += (int)shape_delta;
                int num_edges = (int)((header32 >> 3) + 1);
                S2ClippedShape clipped = new(shape_id, num_edges, (header32 & 4) != 0);
                if (!DecodeEdges(num_edges, clipped, decoder)) return false;
                shapes_.Add(clipped);
            }
        }
        return true;
    }

    private static void EncodeEdges(S2ClippedShape clipped, Encoder encoder)
    {
        // Each entry is an (edge_id, count) pair representing a contiguous range of
        // edges.  The edge ids are delta-encoded such that 0 represents the minimum
        // valid next edge id.
        //
        // Encoding: if bits 0-2 < 7: encodes (count - 1)
        //            - bits 3+: edge delta
        //           if bits 0-2 == 7:
        //            - bits 3+ encode (count - 8)
        //            - Next value is edge delta
        //
        // No count is encoded for the last edge (saving 3 bits).
        int edge_id_base = 0;
        int num_edges = clipped.NumEdges;
        for (int i = 0; i < num_edges; ++i)
        {
            int edge_id = clipped.Edges[i];
            MyDebug.Assert(edge_id >= edge_id_base);
            int delta = edge_id - edge_id_base;
            if (i + 1 == num_edges)
            {
                // This is the last edge; no need to encode an edge count.
                encoder.PutVarInt32(delta);
            }
            else
            {
                // Count the edges in this contiguous range.
                int count = 1;
                for (; i + 1 < num_edges && clipped.Edges[i + 1] == edge_id + count; ++i)
                {
                    ++count;
                }
                if (count < 8)
                {
                    // Count is encoded in low 3 bits of delta.
                    encoder.PutVarInt32(delta << 3 | (count - 1));
                }
                else
                {
                    // Count and delta are encoded separately.
                    encoder.PutVarInt32((count - 8) << 3 | 7);
                    encoder.PutVarInt32(delta);
                }
                edge_id_base = edge_id + count;
            }
        }
    }
    private static bool DecodeEdges(int num_edges, S2ClippedShape clipped, Decoder decoder)
    {
        // This function inverts the encodings documented above.
        Int32 edge_id = 0;
        for (int i = 0; i < num_edges;)
        {
            if (!decoder.TryGetVarUInt32(out var delta)) return false;
            if (i + 1 == num_edges)
            {
                // The last edge is encoded without an edge count.
                clipped.Edges[i++] = edge_id + (int)delta;
            }
            else
            {
                // Otherwise decode the count and edge delta.
                UInt32 count = (delta & 7) + 1;
                delta >>= 3;
                if (count == 8)
                {
                    count = delta + 8;
                    if (!decoder.TryGetVarUInt32(out delta)) return false;
                }
                edge_id += (int)delta;
                for (; count > 0; --count, ++i, ++edge_id)
                {
                    clipped.Edges[i] = edge_id;
                }
            }
        }
        return true;
    }

    // Allocate room for "n" additional clipped shapes in the cell, and return a
    // pointer to the first new clipped shape.  Expects that all new clipped
    // shapes will have a larger shape id than any current shape, and that shapes
    // will be added in increasing shape id order.
    public void AddCapacity(int n)
    {
        int size = shapes_.Count;
        shapes_.Capacity = size + n;
    }
    public void AddShape(S2ClippedShape s)
    {
        shapes_.Add(s);
    }

    private readonly List<S2ClippedShape> shapes_ = [];
}

// S2ClippedShape represents the part of a shape that intersects an S2Cell.
// It consists of the set of edge ids that intersect that cell, and a boolean
// indicating whether the center of the cell is inside the shape (for shapes
// that have an interior).
//
// Note that the edges themselves are not clipped; we always use the original
// edges for intersection tests so that the results will be the same as the
// original shape.
public class S2ClippedShape(Int32 shape_id, Int32 num_edges, bool containsCenter = false)
{
    #region Fields and Properties

    // The shape id of the clipped shape.
    public int ShapeId { get; } = shape_id;

    // Returns true if the center of the S2CellId is inside the shape.  Returns
    // false for shapes that do not have an interior.
    //
    // Set "contains_center_" to indicate whether this clipped shape contains the
    // center of the cell to which it belongs.
    public bool ContainsCenter { get; set; } = containsCenter;

    // The number of edges that intersect the S2CellId.
    public int NumEdges { get; } = num_edges;

    // If there are more than two edges, this field holds a pointer.
    // Otherwise it holds an array of edge ids.
    //
    // Returns the edge id of the given edge in this clipped shape.  Edges are
    // sorted in increasing order of edge id.
    //
    // REQUIRES: 0 <= i < num_edges()
    //
    // Set the i-th edge of this clipped shape to be the given edge of the
    // original shape.
    //
    // Owned by the containing S2ShapeIndexCell.
    public Int32[] Edges { get; } = new Int32[num_edges];

    #endregion
    #region Constructor

    #endregion

    // Returns true if the clipped shape contains the given edge id.
    public bool ContainsEdge(int id)
    {
        // Linear search is fast because the number of edges per shape is typically
        // very small (less than 10).
        for (int e = 0; e < NumEdges; ++e)
        {
            if (Edges[e] == id) return true;
        }
        return false;
    }
}
