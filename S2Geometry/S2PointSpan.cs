// S2PointSpan represents a view of an S2Point array.  It is used to pass
// vertex arrays to functions that don't care about the actual array type
// (e.g. std::vector<S2Point> or S2Point[]).
//
// NOTE: S2PointSpan has an implicit constructor from any container type with
// data() and size() methods (such as std::vector and std::array).  Therefore
// you can use such containers as arguments for any S2PointSpan parameter.
global using S2PointSpan = System.Collections.Generic.List<S2Geometry.Vector3<double>>;
global using S2PointLoopSpan = System.Collections.Generic.List<S2Geometry.Vector3<double>>;

namespace S2Geometry;


// Like S2PointSpan, except that operator[] maps index values in the range
// [n, 2*n-1] to the range [0, n-1] by subtracting n (where n == size()).
// In other words, two full copies of the vertex array are available.  (This
// is a compromise between convenience and efficiency, since computing the
// index modulo "n" is surprisingly expensive.)
//
// This property is useful for implementing algorithms where the elements of
// the span represent the vertices of a loop.
public static class S2PointLoopSpanExtensions
{
    // Like operator[], but allows index values in the range [0, 2*size()-1]
    // where each index i >= size() is mapped to i - size().
    public static S2Point GetPoint(this IList<S2Point> span, int i)
    {
        MyDebug.Assert(i >= 0);
        MyDebug.Assert(i < 2 * span.Count);

        int j = i - span.Count;
        return span[j < 0 ? i : j];
    }
}
