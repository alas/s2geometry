namespace S2Geometry;

public class S2ProjectionsTests
{
    [Fact]
    internal void Test_PlateCarreeProjection_Interpolate()
    {
        _ = new PlateCarreeProjection(180);

        // Test that coordinates and/or arguments are not accidentally reversed.
        Assert.Equal(new R2Point(1.5, 6),
                  Projection.Interpolate(0.25, new R2Point(1, 5), new R2Point(3, 9)));

        // Test extrapolation.
        Assert.Equal(new R2Point(-3, 0), Projection.Interpolate(-2, new R2Point(1, 0), new R2Point(3, 0)));

        // Check that interpolation is exact at both endpoints.
        var a = new R2Point(1.234, -5.456e-20);
        var b = new R2Point(2.1234e-20, 7.456);
        Assert.Equal(a, Projection.Interpolate(0, a, b));
        Assert.Equal(b, Projection.Interpolate(1, a, b));
    }

    void TestProjectUnproject(Projection projection, R2Point px, S2Point x)
    {
        // The arguments are chosen such that projection is exact, but
        // unprojection may not be.
        Assert.Equal(px, projection.Project(x));
        Assert.True(S2.ApproxEquals(x, projection.Unproject(px)));
    }

    [Fact]
    internal void Test_PlateCarreeProjection_ProjectUnproject()
    {
        PlateCarreeProjection proj = new(180);
        TestProjectUnproject(proj, new R2Point(0, 0), new S2Point(1, 0, 0));
        TestProjectUnproject(proj, new R2Point(180, 0), new S2Point(-1, 0, 0));
        TestProjectUnproject(proj, new R2Point(90, 0), new S2Point(0, 1, 0));
        TestProjectUnproject(proj, new R2Point(-90, 0), new S2Point(0, -1, 0));
        TestProjectUnproject(proj, new R2Point(0, 90), new S2Point(0, 0, 1));
        TestProjectUnproject(proj, new R2Point(0, -90), new S2Point(0, 0, -1));
    }

    [Fact]
    internal void Test_MercatorProjection_ProjectUnproject()
    {
        MercatorProjection proj = new(180);
        double inf = double.PositiveInfinity;
        TestProjectUnproject(proj, new R2Point(0, 0), new S2Point(1, 0, 0));
        TestProjectUnproject(proj, new R2Point(180, 0), new S2Point(-1, 0, 0));
        TestProjectUnproject(proj, new R2Point(90, 0), new S2Point(0, 1, 0));
        TestProjectUnproject(proj, new R2Point(-90, 0), new S2Point(0, -1, 0));
        TestProjectUnproject(proj, new R2Point(0, inf), new S2Point(0, 0, 1));
        TestProjectUnproject(proj, new R2Point(0, -inf), new S2Point(0, 0, -1));

        // Test one arbitrary point as a sanity check.
        TestProjectUnproject(proj, new R2Point(0, 70.255578967830246), S2LatLng.FromRadians(1, 0).ToPoint());
    }
}
