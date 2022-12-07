namespace S2Geometry;

public class S2RegionUnionTests
{
    [Fact]
    internal void Test_S2RegionUnionTest_Basic()
    {
        S2RegionUnion ru_empty = new(new List<IS2Region>());
        Assert.Equal(0, ru_empty.Count());
        Assert.Equal(S2Cap.Empty, ru_empty.GetCapBound());
        Assert.Equal(S2LatLngRect.Empty, ru_empty.GetRectBound());
        var empty_clone = (S2RegionUnion)ru_empty.CustomClone();

        var two_point_region = new List<IS2Region>
        {
            new S2PointRegion(S2LatLng.FromDegrees(35, 40).ToPoint()),
            new S2PointRegion(S2LatLng.FromDegrees(-35, -40).ToPoint())
        };

        var two_points_orig = new S2RegionUnion(two_point_region);
        // two_point_region is in a valid, but unspecified, state.

        // Check that Clone() returns a deep copy.
        var two_points = (S2RegionUnion)two_points_orig.CustomClone();
        // The bounds below may not be exactly equal because the S2PointRegion
        // version converts each S2LatLng value to an S2Point and back.
        Assert.True(MakeLatLngRectOrDie("-35:-40,35:40")!.Value
            .ApproxEquals(two_points.GetRectBound()));

        S2Cell face0 = S2Cell.FromFace(0);
        Assert.True(two_points.MayIntersect(face0));
        Assert.False(two_points.Contains(face0));

        Assert.True(two_points.Contains(S2LatLng.FromDegrees(35, 40).ToPoint()));
        Assert.True(two_points.Contains(S2LatLng.FromDegrees(-35, -40).ToPoint()));
        Assert.False(two_points.Contains(S2LatLng.FromDegrees(0, 0).ToPoint()));

        // Check that we can Add() another region.
        var three_points = (S2RegionUnion)two_points.CustomClone();
        Assert.False(three_points.Contains(S2LatLng.FromDegrees(10, 10).ToPoint()));
        three_points.Add(new S2RegionUnion(new List<IS2Region> { new S2PointRegion(S2LatLng.FromDegrees(10, 10).ToPoint()) }));
        Assert.True(three_points.Contains(S2LatLng.FromDegrees(10, 10).ToPoint()));

        var coverer = new S2RegionCoverer();
        coverer.Options_.MaxCells = 1;
        coverer.GetCovering(two_points, out var covering);
        Assert.Single(covering);
        Assert.Equal(face0.Id, covering[0]);
    }
}
