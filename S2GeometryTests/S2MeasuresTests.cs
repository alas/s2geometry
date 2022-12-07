namespace S2Geometry;

public class S2MeasuresTests
{
    [Fact]
    internal void Test_S2_AngleMethods()
    {
        S2Point pz = new(0, 0, 1);
        S2Point p000 = new(1, 0, 0);
        S2Point p045 = new S2Point(1, 1, 0).Normalize();
        S2Point p090 = new(0, 1, 0);
        S2Point p180 = new(-1, 0, 0);

        Assert2.Near(S2.Angle(p000, pz, p045), S2.M_PI_4);
        Assert2.Near(S2.TurnAngle(p000, pz, p045), -3 * S2.M_PI_4);

        Assert2.Near(S2.Angle(p045, pz, p180), 3 * S2.M_PI_4);
        Assert2.Near(S2.TurnAngle(p045, pz, p180), -S2.M_PI_4);

        Assert2.Near(S2.Angle(p000, pz, p180), Math.PI);
        Assert2.Near(S2.TurnAngle(p000, pz, p180), 0);

        Assert2.Near(S2.Angle(pz, p000, p045), S2.M_PI_2);
        Assert2.Near(S2.TurnAngle(pz, p000, p045), S2.M_PI_2);

        Assert2.Near(S2.Angle(pz, p000, pz), 0);
        Assert2.Near(Math.Abs(S2.TurnAngle(pz, p000, pz)), Math.PI);
    }

    [Fact]
    internal void Test_S2_AreaMethods()
    {
        S2Point pz = new(0, 0, 1);
        S2Point p000 = new(1, 0, 0);
        S2Point p045 = new S2Point(1, 1, 0).Normalize();
        S2Point p090 = new(0, 1, 0);
        S2Point p180 = new(-1, 0, 0);

        Assert2.Near(S2.Area(p000, p090, pz), S2.M_PI_2);
        Assert2.Near(S2.Area(p045, pz, p180), 3 * S2.M_PI_4);

        // Make sure that Area() has good *relative* accuracy even for
        // very small areas.
        const double eps = 1e-10;
        S2Point pepsx = new S2Point(eps, 0, 1).Normalize();
        S2Point pepsy = new S2Point(0, eps, 1).Normalize();
        double expected1 = 0.5 * eps * eps;
        Assert2.Near(S2.Area(pepsx, pepsy, pz), expected1, 1e-14 * expected1);

        // Make sure that it can handle degenerate triangles.
        S2Point pr = new S2Point(0.257, -0.5723, 0.112).Normalize();
        S2Point pq = new S2Point(-0.747, 0.401, 0.2235).Normalize();
        Assert.Equal(0, S2.Area(pr, pr, pr));
        // The following test is not exact due to rounding error.
        Assert2.Near(S2.Area(pr, pq, pr), 0, S2.DoubleError);
        Assert.Equal(0, S2.Area(p000, p045, p090));

        double max_girard = 0;
        for (int i = 0; i < 10000; ++i)
        {
            S2Point p0 = S2Testing.RandomPoint();
            S2Point d1 = S2Testing.RandomPoint();
            S2Point d2 = S2Testing.RandomPoint();
            S2Point p1 = (p0 + S2.DoubleError * d1).Normalize();
            S2Point p2 = (p0 + S2.DoubleError * d2).Normalize();
            // The actual displacement can be as much as 1.2e-15 due to roundoff.
            // This yields a maximum triangle area of about 0.7e-30.
            Assert.True(S2.Area(p0, p1, p2) <= 0.7e-30);
            max_girard = Math.Max(max_girard, S2.GirardArea(p0, p1, p2));
        }
        // This check only passes if GirardArea() uses RobustCrossProd().
        Assert.True(max_girard <= 1e-14);

        // Try a very long and skinny triangle.
        S2Point p045eps = new S2Point(1, 1, eps).Normalize();
        double expected2 = 5.8578643762690495119753e-11;  // Mathematica.
        Assert2.Near(S2.Area(p000, p045eps, p090), expected2, 1e-9 * expected2);

        // Triangles with near-180 degree edges that sum to a quarter-sphere.
        const double eps2 = 1e-14;
        S2Point p000eps2 = new S2Point(1, 0.1 * eps2, eps2).Normalize();
        double quarter_area1 = S2.Area(p000eps2, p000, p045) +
                               S2.Area(p000eps2, p045, p180) +
                               S2.Area(p000eps2, p180, pz) +
                               S2.Area(p000eps2, pz, p000);
        Assert2.Near(quarter_area1, Math.PI);

        // Four other triangles that sum to a quarter-sphere.
        S2Point p045eps2 = new S2Point(1, 1, eps2).Normalize();
        double quarter_area2 = S2.Area(p045eps2, p000, p045) +
                               S2.Area(p045eps2, p045, p180) +
                               S2.Area(p045eps2, p180, pz) +
                               S2.Area(p045eps2, pz, p000);
        Assert2.Near(quarter_area2, Math.PI);

        // Compute the area of a hemisphere using four triangles with one near-180
        // degree edge and one near-degenerate edge.
        for (int i = 0; i < 100; ++i)
        {
            double lng = S2.M_2_PI * S2Testing.Random.RandDouble();
            S2Point p0 = S2LatLng.FromRadians(1e-20, lng).Normalized().ToPoint();
            S2Point p1 = S2LatLng.FromRadians(0, lng).Normalized().ToPoint();
            double p2_lng = lng + S2Testing.Random.RandDouble();
            S2Point p2 = S2LatLng.FromRadians(0, p2_lng).Normalized().ToPoint();
            S2Point p3 = S2LatLng.FromRadians(0, lng + Math.PI).Normalized().ToPoint();
            S2Point p4 = S2LatLng.FromRadians(0, lng + 5.0).Normalized().ToPoint();
            double area = (S2.Area(p0, p1, p2) + S2.Area(p0, p2, p3) +
                           S2.Area(p0, p3, p4) + S2.Area(p0, p4, p1));
            Assert2.Near(area, S2.M_2_PI, 2e-15);
        }

        // This tests a case where the triangle has zero area, but S2.Area()
        // computes (dmin > 0) due to rounding errors.
        Assert.Equal(0.0, S2.Area(S2LatLng.FromDegrees(-45, -170).ToPoint(),
                                S2LatLng.FromDegrees(45, -170).ToPoint(),
                                S2LatLng.FromDegrees(0, -170).ToPoint()));
    }
}
