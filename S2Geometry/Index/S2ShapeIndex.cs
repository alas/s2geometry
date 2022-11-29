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

public abstract class S2ShapeIndex
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

    // The possible relationships between a "target" cell and the cells of the
    // S2ShapeIndex.  If the target is an index cell or is contained by an index
    // cell, it is "INDEXED".  If the target is subdivided into one or more
    // index cells, it is "SUBDIVIDED".  Otherwise it is "DISJOINT".
    public enum CellRelation
    {
        INDEXED,       // Target is contained by an index cell
        SUBDIVIDED,    // Target is subdivided into one or more index cells
        DISJOINT       // Target does not intersect any index cells
    };

    // When passed to an Iterator constructor, specifies whether the iterator
    // should be positioned at the beginning of the index (BEGIN), the end of
    // the index (END), or arbitrarily (UNPOSITIONED).  By default iterators are
    // unpositioned, since this avoids an extra seek in this situation where one
    // of the seek methods (such as Locate) is immediately called.
    public enum InitialPosition { BEGIN, END, UNPOSITIONED };

    // A random access iterator that provides low-level access to the cells of
    // the index.  Cells are sorted in increasing order of S2CellId.
    public class Enumerator : IReversableEnumerator<S2ShapeIndexIdCell>
    {
        private readonly IReversableEnumerator<S2ShapeIndexIdCell> _it;

        // Constructs an iterator positioned as specified.  By default iterators
        // are unpositioned, since this avoids an extra seek in this situation
        // where one of the seek methods (such as Locate) is immediately called.
        //
        // If you want to position the iterator at the beginning, e.g. in order to
        // loop through the entire index, do this instead:
        //
        //   for (S2ShapeIndex.Iterator it(out index, S2ShapeIndex.InitialPosition.BEGIN);
        //        !it.done(); it.Next()) { ... }
        public Enumerator(S2ShapeIndex index)
        {
            _it = index.GetNewEnumerator();
        }

        // Returns the S2CellId of the current index cell.  If done() is true,
        // returns a value larger than any valid S2CellId (S2CellId.Sentinel()).
        public S2CellId Id => _it.Current.Item1;

        // Returns a reference to the contents of the current index cell.
        // REQUIRES: !done()
        public S2ShapeIndexCell Cell => _it.Current.Item2;

        public S2ShapeIndexIdCell Current => _it.Current;

        object IEnumerator.Current => Current;

        // Positions the iterator at the next index cell.
        // REQUIRES: !done()
        public bool MoveNext() => _it.MoveNext();

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        public bool MovePrevious() => _it.MovePrevious();
        public void Reset() => _it.Reset();
        public bool Done() => _it.Done();
        public void SetPosition(int position) => _it.SetPosition(position);
        public void Dispose() { GC.SuppressFinalize(this); }
    }

    public class ReversableEnumerator<T> : IReversableEnumerator<T>
    {
        private readonly IList<T> _arr;
        private int position;
        public ReversableEnumerator(IList<T> arr)
        { _arr = arr; position = -1; }

        public T Current => _arr[position];

        object IEnumerator.Current => Current;

        public void Dispose() { GC.SuppressFinalize(this); }
        public bool MoveNext() => ++position >= 0 && position < _arr.Count;
        public bool MovePrevious() => --position >= 0 && position < _arr.Count;
        public void Reset() => position = -1;
        public void SetPosition(int pos) => position = pos;
        public bool Done() => position >= _arr.Count;
    }

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
            yield return Shape(i);
        }
    }

    public abstract IReversableEnumerator<S2ShapeIndexIdCell> GetNewEnumerator();
    public abstract IEnumerable<S2ShapeIndexIdCell> GetNewEnumerable();
    public abstract int GetEnumerableCount();

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and
    // returns INDEXED.  Otherwise if T contains one or more (smaller) index
    // cells, it positions the iterator at the first such cell I and returns
    // SUBDIVIDED.  Otherwise it returns DISJOINT.
    public (CellRelation cellRelation, int pos) LocateCell(S2CellId target)
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
            if (id >= target && id.Value.RangeMin() <= target) return (CellRelation.INDEXED, pos);
            if (id <= target.RangeMax()) return (CellRelation.SUBDIVIDED, pos);
        }
        if (--pos >= 0)
        {
            id = GetCellId(pos);

            if (id != null && id.Value.RangeMax() >= target)
                return (CellRelation.INDEXED, pos);
        }
        return (CellRelation.DISJOINT, 0);
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
    public abstract S2ShapeIndexIdCell? GetIndexCell(int pos);
}

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
public class S2ShapeIndexCell
{
    public S2ShapeIndexCell() { }

    // Returns the number of clipped shapes in this cell.
    public int NumClipped() { return shapes_.Count; }

    // Returns the clipped shape at the given index.  Shapes are kept sorted in
    // increasing order of shape id.
    //
    // REQUIRES: 0 <= i < num_clipped()
    public S2ClippedShape Clipped(int i) { return shapes_[i]; }

    // Returns a pointer to the clipped shape corresponding to the given shape,
    // or null if the shape does not intersect this cell.
    public S2ClippedShape FindClipped(S2Shape shape)
    {
        return FindClipped(shape.Id);
    }
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
            System.Diagnostics.Debug.Assert(NumClipped() == 1);  // Index invariant: no empty cells.
            S2ClippedShape clipped = Clipped(0);
            System.Diagnostics.Debug.Assert(clipped.ShapeId == 0);
            int n = clipped.NumEdges;
            encoder.Ensure(Encoder.kVarintMax64 + n * Encoder.kVarintMax32);
            var ccc = clipped.ContainsCenter ? 2 : 0;
            if (n >= 2 && n <= 17 && clipped.Edge(n - 1) - clipped.Edge(0) == n - 1)
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
                encoder.PutVarInt64(clipped.Edge(0) << 6 | (n - 2) << 2 | ccc | 0);
            }
            else if (n == 1)
            {
                // The cell contains only one edge.  For edge ids up to 15, we can
                // encode the cell in a single byte.
                //
                // Encoding: bits 0-1: 1
                //           bit 2: contains_center
                //           bits 3+: edge_id
                encoder.PutVarInt64(clipped.Edge(0) << 3 | ccc | 1);
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
                System.Diagnostics.Debug.Assert(clipped.ShapeId >= shape_id_base);
                int shape_delta = clipped.ShapeId - shape_id_base;
                shape_id_base = clipped.ShapeId + 1;

                // Like the code above except that we also need to encode shape_id(s).
                // Because of this some choices are slightly different.
                int n = clipped.NumEdges;
                encoder.Ensure((n + 2) * Encoder.kVarintMax32);
                var ccc = clipped.ContainsCenter ? 2 : 0;
                if (n >= 1 && n <= 16 && clipped.Edge(n - 1) - clipped.Edge(0) == n - 1)
                {
                    // The clipped shape has a contiguous range of up to 16 edges.  This
                    // encoding uses a 1-bit tag because it is by far the most common.
                    //
                    // Encoding: bit 0: 0
                    //           bit 1: contains_center
                    //           bits 2+: edge_id
                    // Next value: bits 0-3: (num_edges - 1)
                    //             bits 4+: shape_delta
                    encoder.PutVarInt32(clipped.Edge(0) << 2 | ccc | 0);
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
            var clippedTmp = new S2ClippedShape();
            shapes_.Add(clippedTmp);
            if (!decoder.TryGetVarUInt64(out var header)) return false;
            if ((header & 1) == 0)
            {
                // The cell contains a contiguous range of edges.
                int num_edges = (int)(((header >> 2) & 15) + 2);
                clippedTmp.Init(0 /*shape_id*/, num_edges);
                clippedTmp.ContainsCenter = ((header & 2) != 0);
                for (int i = 0, edge_id = (int)(header >> 6); i < num_edges; ++i)
                {
                    clippedTmp.SetEdge(i, edge_id + i);
                }
                return true;
            }
            if ((header & 2) == 0)
            {
                // The cell contains a single edge.
                clippedTmp.Init(0 /*shape_id*/, 1 /*num_edges*/);
                clippedTmp.ContainsCenter = ((header & 4) != 0);
                clippedTmp.SetEdge(0, (int)(header >> 3));
                return true;
            }
            else
            {
                // The cell contains some other combination of edges.
                int num_edges = (int)(header >> 3);
                clippedTmp.Init(0 /*shape_id*/, num_edges);
                clippedTmp.ContainsCenter = ((header & 4) != 0);
                return DecodeEdges(num_edges, clippedTmp, decoder);
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
            var clipped = new S2ClippedShape();
            shapes_.Add(clipped);
            if (j > 0 && !decoder.TryGetVarUInt32(out header32)) return false;
            if ((header32 & 1) == 0)
            {
                // The clipped shape contains a contiguous range of edges.
                if (!decoder.TryGetVarUInt32(out var shape_id_count)) return false;
                shape_id += (int)(shape_id_count >> 4);
                int num_edges = (int)((shape_id_count & 15) + 1);
                clipped.Init(shape_id, num_edges);
                clipped.ContainsCenter = ((header32 & 2) != 0);
                for (int i = 0, edge_id = (int)(header32 >> 2); i < num_edges; ++i)
                {
                    clipped.SetEdge(i, edge_id + i);
                }
            }
            else if ((header32 & 7) == 7)
            {
                // The clipped shape has no edges.
                shape_id += (int)(header32 >> 4);
                clipped.Init(shape_id, 0);
                clipped.ContainsCenter = ((header32 & 8) != 0);
            }
            else
            {
                // The clipped shape contains some other combination of edges.
                System.Diagnostics.Debug.Assert((header32 & 3U) == 1U);
                if (!decoder.TryGetVarUInt32(out var shape_delta)) return false;
                shape_id += (int)shape_delta;
                int num_edges = (int)((header32 >> 3) + 1);
                clipped.Init(shape_id, num_edges);
                clipped.ContainsCenter = ((header32 & 4) != 0);
                if (!DecodeEdges(num_edges, clipped, decoder)) return false;
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
            int edge_id = clipped.Edge(i);
            System.Diagnostics.Debug.Assert(edge_id >= edge_id_base);
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
                for (; i + 1 < num_edges && clipped.Edge(i + 1) == edge_id + count; ++i)
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
                clipped.SetEdge(i++, edge_id + (int)delta);
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
                    clipped.SetEdge(i, edge_id);
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

    private readonly List<S2ClippedShape> shapes_ = new();
}

// S2ClippedShape represents the part of a shape that intersects an S2Cell.
// It consists of the set of edge ids that intersect that cell, and a boolean
// indicating whether the center of the cell is inside the shape (for shapes
// that have an interior).
//
// Note that the edges themselves are not clipped; we always use the original
// edges for intersection tests so that the results will be the same as the
// original shape.
public class S2ClippedShape
{
    public S2ClippedShape() { }
    public S2ClippedShape(Int32 shape_id, Int32 num_edges) => Init(shape_id, num_edges);

    // Initialize an S2ClippedShape to hold the given number of edges.
    public void Init(Int32 shape_id, Int32 num_edges)
    {
        ShapeId = shape_id;
        NumEdges = num_edges;
        ContainsCenter = false;
        edges_ = new Int32[num_edges];
    }

    // The shape id of the clipped shape.
    public int ShapeId { get; private set; }

    // Returns true if the center of the S2CellId is inside the shape.  Returns
    // false for shapes that do not have an interior.
    //
    // Set "contains_center_" to indicate whether this clipped shape contains the
    // center of the cell to which it belongs.
    public bool ContainsCenter { get; set; } = true;

    // The number of edges that intersect the S2CellId.
    public int NumEdges { get; private set; } = 31;

    // Returns the edge id of the given edge in this clipped shape.  Edges are
    // sorted in increasing order of edge id.
    //
    // REQUIRES: 0 <= i < num_edges()
    public int Edge(int i)
    {
        return edges_[i];
    }
    // Set the i-th edge of this clipped shape to be the given edge of the
    // original shape.
    public void SetEdge(int i, int edge)
    {
        edges_[i] = edge;
    }

    // Returns true if the clipped shape contains the given edge id.
    public bool ContainsEdge(int id)
    {
        // Linear search is fast because the number of edges per shape is typically
        // very small (less than 10).
        for (int e = 0; e < NumEdges; ++e)
        {
            if (Edge(e) == id) return true;
        }
        return false;
    }

    // If there are more than two edges, this field holds a pointer.
    // Otherwise it holds an array of edge ids.
    private Int32[] edges_;  // Owned by the containing S2ShapeIndexCell.
}
