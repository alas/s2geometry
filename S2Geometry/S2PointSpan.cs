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
        Assert.True(i >= 0);
        Assert.True(i < 2 * span.Count);

        int j = i - span.Count;
        return span[j < 0 ? i : j];
    }
}
