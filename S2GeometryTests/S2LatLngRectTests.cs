// Most of the S2LatLngRect methods have trivial implementations that
// use the R1Interval and S1Interval classes, so most of the testing
// is done in those unit tests.

namespace S2Geometry;

public class S2LatLngRectTests
{
    [Fact]
    internal void Test_S2LatLngRect_EmptyAndFull()
    {
        // Test basic properties of empty and full rectangles.
        S2LatLngRect empty = S2LatLngRect.Empty;
        S2LatLngRect full = S2LatLngRect.Full;
        Assert.True(empty.IsValid());
        Assert.True(empty.IsEmpty());
        Assert.False(empty.IsPoint());
        Assert.True(full.IsValid());
        Assert.True(full.IsFull());
        Assert.False(full.IsPoint());
        // Check that the default S2LatLngRect is identical to Empty().
        S2LatLngRect default_empty = S2LatLngRect.Empty;
        Assert.True(default_empty.IsValid());
        Assert.True(default_empty.IsEmpty());
        Assert.Equal(empty.Lat, default_empty.Lat);
        Assert.Equal(empty.Lng, default_empty.Lng);
    }

    [Fact]
    internal void Test_S2LatLngRect_Accessors()
    {
        // Check various accessor methods.
        S2LatLngRect d1 = RectFromDegrees(-90, 0, -45, 180);
        Assert2.DoubleEqual(d1.LatLo().GetDegrees(), -90);
        Assert2.DoubleEqual(d1.LatHi().GetDegrees(), -45);
        Assert2.DoubleEqual(d1.LngLo().GetDegrees(), 0);
        Assert2.DoubleEqual(d1.LngHi().GetDegrees(), 180);
        Assert.Equal(d1.Lat, new(-S2.M_PI_2, -S2.M_PI_4));
        Assert.Equal(d1.Lng, new(0, Math.PI));
    }

    [Fact]
    internal void Test_S2LatLngRect_ApproxEquals()
    {
        // S1Interval and R1Interval have additional testing.

        Assert.True(S2LatLngRect.Empty.ApproxEquals(RectFromDegrees(1, 5, 1, 5)));
        Assert.True(RectFromDegrees(1, 5, 1, 5).ApproxEquals(S2LatLngRect.Empty));
        Assert.False(RectFromDegrees(1, 5, 1, 5).ApproxEquals(RectFromDegrees(2, 7, 2, 7)));

        // Test the max_error (double) parameter.
        Assert.True(RectFromDegrees(10, 10, 20, 20).ApproxEquals(RectFromDegrees(11, 11, 19, 19), S1Angle.FromDegrees(1.001)));
        Assert.False(RectFromDegrees(10, 10, 20, 20).ApproxEquals(RectFromDegrees(11, 11, 19, 19), S1Angle.FromDegrees(0.999)));

        // Test the max_error (S2LatLng) parameter.
        Assert.True(RectFromDegrees(0, 10, 20, 30).ApproxEquals(RectFromDegrees(-1, 8, 21, 32), S2LatLng.FromDegrees(1.001, 2.001)));
        Assert.False(RectFromDegrees(0, 10, 20, 30).ApproxEquals(RectFromDegrees(-1, 8, 21, 32), S2LatLng.FromDegrees(0.999, 1.999)));
    }

    [Fact]
    internal void Test_S2LatLngRect_FromCenterSize()
    {
        Assert.True(S2LatLngRect.FromCenterSize(S2LatLng.FromDegrees(80, 170), S2LatLng.FromDegrees(40, 60)).ApproxEquals(RectFromDegrees(60, 140, 90, -160)));
        Assert.True(S2LatLngRect.FromCenterSize(S2LatLng.FromDegrees(10, 40), S2LatLng.FromDegrees(210, 400)).IsFull());
        Assert.True(S2LatLngRect.FromCenterSize(S2LatLng.FromDegrees(-90, 180), S2LatLng.FromDegrees(20, 50)).ApproxEquals(RectFromDegrees(-90, 155, -80, -155)));
    }

    [Fact]
    internal void Test_S2LatLngRect_FromPoint()
    {
        S2LatLng p = S2LatLng.FromDegrees(23, 47);
        Assert.Equal(S2LatLngRect.FromPoint(p), new S2LatLngRect(p, p));
        Assert.True(S2LatLngRect.FromPoint(p).IsPoint());
    }

    [Fact]
    internal void Test_S2LatLngRect_FromPointPair()
    {
        Assert.Equal(S2LatLngRect.FromPointPair(S2LatLng.FromDegrees(-35, -140), S2LatLng.FromDegrees(15, 155)), RectFromDegrees(-35, 155, 15, -140));
        Assert.Equal(S2LatLngRect.FromPointPair(S2LatLng.FromDegrees(25, -70), S2LatLng.FromDegrees(-90, 80)), RectFromDegrees(-90, -70, 25, 80));
    }

    [Fact]
    internal void Test_S2LatLngRect_GetCenterSize()
    {
        S2LatLngRect r1 = new(new R1Interval(0, S2.M_PI_2), new S1Interval(-Math.PI, 0));
        Assert.Equal(r1.Center(), S2LatLng.FromRadians(S2.M_PI_4, -S2.M_PI_2));
        Assert.Equal(r1.Size(), S2LatLng.FromRadians(S2.M_PI_2, Math.PI));
        Assert.True(S2LatLngRect.Empty.Size().LatRadians < 0);
        Assert.True(S2LatLngRect.Empty.Size().LngRadians < 0);
    }

    [Fact]
    internal void Test_S2LatLngRect_GetVertex()
    {
        S2LatLngRect r1 = new(new R1Interval(0, S2.M_PI_2), new S1Interval(-Math.PI, 0));
        Assert.Equal(r1.Vertex(0), S2LatLng.FromRadians(0, Math.PI));
        Assert.Equal(r1.Vertex(1), S2LatLng.FromRadians(0, 0));
        Assert.Equal(r1.Vertex(2), S2LatLng.FromRadians(S2.M_PI_2, 0));
        Assert.Equal(r1.Vertex(3), S2LatLng.FromRadians(S2.M_PI_2, Math.PI));

        // Make sure that GetVertex() returns vertices in CCW order.
        for (int i = 0; i < 4; ++i)
        {
            double lat = S2.M_PI_4 * (i - 2);
            double lng = S2.M_PI_2 * (i - 2) + 0.2;
            S2LatLngRect r = new(
                new R1Interval(lat, lat + S2.M_PI_4),
                new S1Interval(
                    Math.IEEERemainder(lng, S2.M_2_PI),
                    Math.IEEERemainder(lng + S2.M_PI_2, S2.M_2_PI)));
            for (int k = 0; k < 4; ++k)
            {
                Assert.True(S2Pred.Sign(
                    r.Vertex(k - 1).ToPoint(),
                    r.Vertex(k).ToPoint(),
                    r.Vertex(k + 1).ToPoint()) > 0);
            }
        }
    }

    [Fact]
    internal void Test_S2LatLngRect_Contains()
    {
        // Contains(S2LatLng), InteriorContains(S2LatLng), Contains()
        S2LatLng eq_m180 = S2LatLng.FromRadians(0, -Math.PI);
        S2LatLng north_pole = S2LatLng.FromRadians(S2.M_PI_2, 0);
        S2LatLngRect r1 = new(eq_m180, north_pole);

        Assert.True(r1.Contains(S2LatLng.FromDegrees(30, -45)));
        Assert.True(r1.InteriorContains(S2LatLng.FromDegrees(30, -45)));
        Assert.False(r1.Contains(S2LatLng.FromDegrees(30, 45)));
        Assert.False(r1.InteriorContains(S2LatLng.FromDegrees(30, 45)));
        Assert.True(r1.Contains(eq_m180));
        Assert.False(r1.InteriorContains(eq_m180));
        Assert.True(r1.Contains(north_pole));
        Assert.False(r1.InteriorContains(north_pole));
        Assert.True(r1.Contains(new S2Point(0.5, -0.3, 0.1)));
        Assert.False(r1.Contains(new S2Point(0.5, 0.2, 0.1)));
    }

    [Fact]
    internal void Test_S2LatLngRect_IntervalOps()
    {
        // Contains(S2LatLngRect), InteriorContains(S2LatLngRect),
        // Intersects(), InteriorIntersects(), Union(), Intersection().
        //
        // Much more testing of these methods is done in s1interval_test
        // and r1interval_test.

        // Rectangle "r1" covers one-quarter of the sphere.
        S2LatLngRect r1 = RectFromDegrees(0, -180, 90, 0);

        // Test operations where one rectangle consists of a single point.
        S2LatLngRect r1_mid = RectFromDegrees(45, -90, 45, -90);
        TestIntervalOps(r1, r1_mid, "TTTT", r1, r1_mid);

        S2LatLngRect req_m180 = RectFromDegrees(0, -180, 0, -180);
        TestIntervalOps(r1, req_m180, "TFTF", r1, req_m180);

        S2LatLngRect rnorth_pole = RectFromDegrees(90, 0, 90, 0);
        TestIntervalOps(r1, rnorth_pole, "TFTF", r1, rnorth_pole);

        TestIntervalOps(r1, RectFromDegrees(-10, -1, 1, 20), "FFTT",
                        RectFromDegrees(-10, 180, 90, 20),
                        RectFromDegrees(0, -1, 1, 0));
        TestIntervalOps(r1, RectFromDegrees(-10, -1, 0, 20), "FFTF",
                        RectFromDegrees(-10, 180, 90, 20),
                        RectFromDegrees(0, -1, 0, 0));
        TestIntervalOps(r1, RectFromDegrees(-10, 0, 1, 20), "FFTF",
                        RectFromDegrees(-10, 180, 90, 20),
                        RectFromDegrees(0, 0, 1, 0));

        TestIntervalOps(RectFromDegrees(-15, -160, -15, -150),
                        RectFromDegrees(20, 145, 25, 155), "FFFF",
                        RectFromDegrees(-15, 145, 25, -150),
                        S2LatLngRect.Empty);
        TestIntervalOps(RectFromDegrees(70, -10, 90, -140),
                        RectFromDegrees(60, 175, 80, 5), "FFTT",
                        RectFromDegrees(60, -180, 90, 180),
                        RectFromDegrees(70, 175, 80, 5));

        // Check that the intersection of two rectangles that overlap in latitude
        // but not longitude is valid, and vice versa.
        TestIntervalOps(RectFromDegrees(12, 30, 60, 60),
                        RectFromDegrees(0, 0, 30, 18), "FFFF",
                        RectFromDegrees(0, 0, 60, 60), S2LatLngRect.Empty);
        TestIntervalOps(RectFromDegrees(0, 0, 18, 42),
                        RectFromDegrees(30, 12, 42, 60), "FFFF",
                        RectFromDegrees(0, 0, 42, 60), S2LatLngRect.Empty);
    }

    [Fact]
    internal void Test_BoundaryIntersects_EmptyRectangle()
    {
        S2LatLngRect rect = S2LatLngRect.Empty;
        var lo = rect.Lo().ToPoint();
        var hi = rect.Hi().ToPoint();
        Assert.False(rect.BoundaryIntersects(lo, lo));
        Assert.False(rect.BoundaryIntersects(lo, hi));
    }

    [Fact]
    internal void Test_BoundaryIntersects_FullRectangle()
    {
        S2LatLngRect rect = S2LatLngRect.Full;
        S2Point lo = rect.Lo().ToPoint();
        var hi = rect.Hi().ToPoint();
        Assert.False(rect.BoundaryIntersects(lo, lo));
        Assert.False(rect.BoundaryIntersects(lo, hi));
    }

    [Fact]
    internal void Test_BoundaryIntersects_SphericalLune()
    {
        // This rectangle only has two non-degenerate sides.
        S2LatLngRect rect = RectFromDegrees(-90, 100, 90, 120);
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("60:60"), MakePointOrDie("90:60")));
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("-60:110"), MakePointOrDie("60:110")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("-60:95"), MakePointOrDie("60:110")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("60:115"), MakePointOrDie("80:125")));
    }

    [Fact]
    internal void Test_BoundaryIntersects_NorthHemisphere()
    {
        // This rectangle only has only one non-degenerate side.
        S2LatLngRect rect = RectFromDegrees(0, -180, 90, 180);
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("60:-180"), MakePointOrDie("90:-180")));
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("60:-170"), MakePointOrDie("60:170")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("-10:-180"), MakePointOrDie("10:-180")));
    }

    [Fact]
    internal void Test_BoundaryIntersects_SouthHemisphere()
    {
        // This rectangle only has only one non-degenerate side.
        S2LatLngRect rect = RectFromDegrees(-90, -180, 0, 180);
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("-90:-180"), MakePointOrDie("-60:-180")));
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("-60:-170"), MakePointOrDie("-60:170")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("-10:-180"), MakePointOrDie("10:-180")));
    }

    [Fact]
    internal void Test_BoundaryIntersects_RectCrossingAntiMeridian()
    {
        S2LatLngRect rect = RectFromDegrees(20, 170, 40, -170);
        Assert.True(rect.Contains(MakePointOrDie("30:180")));

        // Check that crossings of all four sides are detected.
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("25:160"), MakePointOrDie("25:180")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("25:-160"), MakePointOrDie("25:-180")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("15:175"), MakePointOrDie("30:175")));
        Assert.True(rect.BoundaryIntersects(
            MakePointOrDie("45:175"), MakePointOrDie("30:175")));

        // Check that the edges on the opposite side of the sphere but at the same
        // latitude do not intersect the rectangle boundary.
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("25:-20"), MakePointOrDie("25:0")));
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("25:20"), MakePointOrDie("25:0")));
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("15:-5"), MakePointOrDie("30:-5")));
        Assert.False(rect.BoundaryIntersects(
            MakePointOrDie("45:-5"), MakePointOrDie("30:-5")));
    }

    [Fact]
    internal void Test_S2LatLngRect_AddPoint()
    {
        S2LatLngRect p = S2LatLngRect.Empty;
        p = p.AddPoint(S2LatLng.FromDegrees(0, 0));
        Assert.True(p.IsPoint());
        p = p.AddPoint(S2LatLng.FromRadians(0, -S2.M_PI_2));
        Assert.False(p.IsPoint());
        p = p.AddPoint(S2LatLng.FromRadians(S2.M_PI_4, -Math.PI));
        p = p.AddPoint(new S2Point(0, 0, 1));
        Assert.Equal(p, RectFromDegrees(0, -180, 90, 0));
    }

    [Fact]
    internal void Test_S2LatLngRect_Expanded()
    {
        Assert.True(RectFromDegrees(70, 150, 80, 170).
                    Expanded(S2LatLng.FromDegrees(20, 30)).
                    ApproxEquals(RectFromDegrees(50, 120, 90, -160)));
        Assert.True(S2LatLngRect.Empty.Expanded(S2LatLng.FromDegrees(20, 30)).
                    IsEmpty());
        Assert.True(S2LatLngRect.Full.Expanded(S2LatLng.FromDegrees(500, 500)).
                    IsFull());
        Assert.True(RectFromDegrees(-90, 170, 10, 20).
                    Expanded(S2LatLng.FromDegrees(30, 80)).
                    ApproxEquals(RectFromDegrees(-90, -180, 40, 180)));

        // Negative margins.
        Assert.True(RectFromDegrees(10, -50, 60, 70).
                    Expanded(S2LatLng.FromDegrees(-10, -10)).
                    ApproxEquals(RectFromDegrees(20, -40, 50, 60)));
        Assert.True(RectFromDegrees(-20, -180, 20, 180).
                    Expanded(S2LatLng.FromDegrees(-10, -10)).
                    ApproxEquals(RectFromDegrees(-10, -180, 10, 180)));
        Assert.True(RectFromDegrees(-20, -180, 20, 180).
                    Expanded(S2LatLng.FromDegrees(-30, -30)).IsEmpty());
        Assert.True(RectFromDegrees(-90, 10, 90, 11).
                    Expanded(S2LatLng.FromDegrees(-10, -10)).IsEmpty());
        Assert.True(RectFromDegrees(-90, 10, 90, 100).
                    Expanded(S2LatLng.FromDegrees(-10, -10)).
                    ApproxEquals(RectFromDegrees(-80, 20, 80, 90)));
        Assert.True(S2LatLngRect.Empty.Expanded(S2LatLng.FromDegrees(-50, -500)).
                    IsEmpty());
        Assert.True(S2LatLngRect.Full.Expanded(S2LatLng.FromDegrees(-50, -50)).
                    ApproxEquals(RectFromDegrees(-40, -180, 40, 180)));

        // Mixed margins.
        Assert.True(RectFromDegrees(10, -50, 60, 70).
                    Expanded(S2LatLng.FromDegrees(-10, 30)).
                    ApproxEquals(RectFromDegrees(20, -80, 50, 100)));
        Assert.True(RectFromDegrees(-20, -180, 20, 180).
                    Expanded(S2LatLng.FromDegrees(10, -500)).
                    ApproxEquals(RectFromDegrees(-30, -180, 30, 180)));
        Assert.True(RectFromDegrees(-90, -180, 80, 180).
                    Expanded(S2LatLng.FromDegrees(-30, 500)).
                    ApproxEquals(RectFromDegrees(-60, -180, 50, 180)));
        Assert.True(RectFromDegrees(-80, -100, 80, 150).
                    Expanded(S2LatLng.FromDegrees(30, -50)).
                    ApproxEquals(RectFromDegrees(-90, -50, 90, 100)));
        Assert.True(RectFromDegrees(0, -180, 50, 180).
                    Expanded(S2LatLng.FromDegrees(-30, 500)).IsEmpty());
        Assert.True(RectFromDegrees(-80, 10, 70, 20).
                    Expanded(S2LatLng.FromDegrees(30, -200)).IsEmpty());
        Assert.True(S2LatLngRect.Empty.Expanded(S2LatLng.FromDegrees(100, -100)).
                    IsEmpty());
        Assert.True(S2LatLngRect.Full.Expanded(S2LatLng.FromDegrees(100, -100)).
                    IsFull());

    }

    [Fact]
    internal void Test_S2LatLngRect_PolarClosure()
    {
        Assert.Equal(RectFromDegrees(-89, 0, 89, 1),
                  RectFromDegrees(-89, 0, 89, 1).PolarClosure());
        Assert.Equal(RectFromDegrees(-90, -180, -45, 180),
                  RectFromDegrees(-90, -30, -45, 100).PolarClosure());
        Assert.Equal(RectFromDegrees(89, -180, 90, 180),
                  RectFromDegrees(89, 145, 90, 146).PolarClosure());
        Assert.Equal(S2LatLngRect.Full,
                  RectFromDegrees(-90, -145, 90, -144).PolarClosure());
    }

    [Fact]
    internal void Test_ExpandedByDistance_PositiveDistance()
    {
        Assert.True(RectFromDegrees(0, 170, 0, -170).
                    ExpandedByDistance(S1Angle.FromDegrees(15)).ApproxEquals(
                        RectFromDegrees(-15, 155, 15, -155)));
        Assert.True(RectFromDegrees(60, 150, 80, 10).
                    ExpandedByDistance(S1Angle.FromDegrees(15)).ApproxEquals(
                        RectFromDegrees(45, -180, 90, 180)));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceNorthEast()
    {
        S2LatLngRect in_rect = RectFromDegrees(0.0, 0.0, 30.0, 90.0);
        S1Angle distance = S1Angle.FromDegrees(5.0);
        S2LatLngRect out_rect = in_rect.ExpandedByDistance(distance).ExpandedByDistance(-distance);
        Assert.True(out_rect.ApproxEquals(in_rect));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceSouthWest()
    {
        S2LatLngRect in_rect = RectFromDegrees(-30.0, -90.0, 0.0, 0.0);
        S1Angle distance = S1Angle.FromDegrees(5.0);

        S2LatLngRect out_rect =
            in_rect.ExpandedByDistance(distance).ExpandedByDistance(-distance);

        Assert.True(out_rect.ApproxEquals(in_rect));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLatWithNorthPole()
    {
        S2LatLngRect rect = RectFromDegrees(0.0, -90.0, 90.0, 180.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        Assert.True(rect.ApproxEquals(RectFromDegrees(5.0, 0.0, 85.0, 90.0)));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLatWithNorthPoleAndLngFull()
    {
        S2LatLngRect rect = RectFromDegrees(0.0, -180.0, 90.0, 180.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        Assert.True(rect.ApproxEquals(RectFromDegrees(5.0, -180.0, 90.0, 180.0)));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLatWithSouthPole()
    {
        S2LatLngRect rect = RectFromDegrees(-90.0, -90.0, 0.0, 180.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        Assert.True(rect.ApproxEquals(RectFromDegrees(-85.0, 0.0, -5.0, 90.0)));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLatWithSouthPoleAndLngFull()
    {
        S2LatLngRect rect = RectFromDegrees(-90.0, -180.0, 0.0, 180.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        Assert.True(rect.ApproxEquals(RectFromDegrees(-90.0, -180.0, -5.0, 180.0)));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLngFull()
    {
        S2LatLngRect rect = RectFromDegrees(0.0, -180.0, 30.0, 180.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        Assert.True(rect.ApproxEquals(RectFromDegrees(5.0, -180.0, 25.0, 180.0)));
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLatResultEmpty()
    {
        S2LatLngRect rect = RectFromDegrees(0.0, 0.0, 9.9, 90.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        Assert.True(rect.IsEmpty());
    }

    [Fact]
    internal void Test_ExpandedByDistance_NegativeDistanceLngResultEmpty()
    {
        S2LatLngRect rect = RectFromDegrees(0.0, 0.0, 30.0, 11.0)
            .ExpandedByDistance(-S1Angle.FromDegrees(5.0));

        // The cap center is at latitude 30 - 5 = 25 degrees. The length of the
        // latitude 25 degree line is 0.906 times the length of the equator. Thus the
        // cap whose radius is 5 degrees covers the rectangle whose latitude interval
        // is 11 degrees.
        Assert.True(rect.IsEmpty());
    }

    [Fact]
    internal void Test_S2LatLngRect_GetCapBound()
    {
        // Bounding cap at center is smaller:
        Assert.True(RectFromDegrees(-45, -45, 45, 45).GetCapBound().
                    ApproxEquals(S2Cap.FromCenterHeight(new S2Point(1, 0, 0), 0.5)));

        // Bounding cap at north pole is smaller:
        Assert.True(RectFromDegrees(88, -80, 89, 80).GetCapBound().
                    ApproxEquals(new S2Cap(new S2Point(0, 0, 1), S1Angle.FromDegrees(2))));

        // Longitude span > 180 degrees:
        Assert.True(RectFromDegrees(-30, -150, -10, 50).GetCapBound().
                    ApproxEquals(new S2Cap(new S2Point(0, 0, -1), S1Angle.FromDegrees(80))));

        // Ensure hemispheres are bounded conservatively.
        Assert.True(RectFromDegrees(-10, -100, 0, 100).GetCapBound().Radius >=
            S1ChordAngle.Right);
    }

    [Fact]
    internal void Test_S2LatLngRect_CellOps()
    {
        // Contains(S2Cell), MayIntersect(S2Cell), Intersects(S2Cell)

        // Special cases.
        TestCellOps(S2LatLngRect.Empty, S2Cell.FromFacePosLevel(3, 0, 0), 0);
        TestCellOps(S2LatLngRect.Full, S2Cell.FromFacePosLevel(2, 0, 0), 4);
        TestCellOps(S2LatLngRect.Full, S2Cell.FromFacePosLevel(5, 0, 25), 4);

        // This rectangle includes the first quadrant of face 0.  It's expanded
        // slightly because cell bounding rectangles are slightly conservative.
        S2LatLngRect r4 = RectFromDegrees(-45.1, -45.1, 0.1, 0.1);
        TestCellOps(r4, S2Cell.FromFacePosLevel(0, 0, 0), 3);
        TestCellOps(r4, S2Cell.FromFacePosLevel(0, 0, 1), 4);
        TestCellOps(r4, S2Cell.FromFacePosLevel(1, 0, 1), 0);

        // This rectangle intersects the first quadrant of face 0.
        S2LatLngRect r5 = RectFromDegrees(-10, -45, 10, 0);
        TestCellOps(r5, S2Cell.FromFacePosLevel(0, 0, 0), 3);
        TestCellOps(r5, S2Cell.FromFacePosLevel(0, 0, 1), 3);
        TestCellOps(r5, S2Cell.FromFacePosLevel(1, 0, 1), 0);

        // Rectangle consisting of a single point.
        TestCellOps(RectFromDegrees(4, 4, 4, 4), S2Cell.FromFace(0), 3);

        // Rectangles that intersect the bounding rectangle of a face
        // but not the face itself.
        TestCellOps(RectFromDegrees(41, -87, 42, -79), S2Cell.FromFace(2), 1);
        TestCellOps(RectFromDegrees(-41, 160, -40, -160), S2Cell.FromFace(5), 1);

        // This is the leaf cell at the top right hand corner of face 0.
        // It has two angles of 60 degrees and two of 120 degrees.
        S2Cell cell0tr = new(new S2Point(1 + 1e-12, 1, 1));
        _ = cell0tr.GetRectBound();
        S2LatLng v0 = new(cell0tr.VertexRaw(0));
        TestCellOps(RectFromDegrees(v0.Lat().GetDegrees() - 1e-8,
                                    v0.Lng().GetDegrees() - 1e-8,
                                    v0.Lat().GetDegrees() - 2e-10,
                                    v0.Lng().GetDegrees() + 1e-10), cell0tr, 1);

        // Rectangles that intersect a face but where no vertex of one region
        // is contained by the other region.  The first one passes through
        // a corner of one of the face cells.
        TestCellOps(RectFromDegrees(-37, -70, -36, -20), S2Cell.FromFace(5), 2);

        // These two intersect like a diamond and a square.
        S2Cell cell202 = S2Cell.FromFacePosLevel(2, 0, 2);
        S2LatLngRect bound202 = cell202.GetRectBound();
        TestCellOps(RectFromDegrees(bound202.Lo().Lat().GetDegrees() + 3,
                                    bound202.Lo().Lng().GetDegrees() + 3,
                                    bound202.Hi().Lat().GetDegrees() - 3,
                                    bound202.Hi().Lng().GetDegrees() - 3), cell202, 2);
    }

    [Fact]
    internal void Test_S2LatLngRect_EncodeDecode()
    {
        S2LatLngRect r = RectFromDegrees(-20, -80, 10, 20);
        Encoder encoder = new();
        r.Encode(encoder);
        var decoder = encoder.GetDecoder();
        var (success, decoded_rect) = S2LatLngRect.Decode(decoder);
        Assert.True(success);
        Assert.Equal(r, decoded_rect);
    }

    [Fact]
    internal void Test_S2LatLngRect_Area()
    {
        Assert.Equal(0.0, S2LatLngRect.Empty.Area());
        Assert2.DoubleEqual(S2.M_4_PI, S2LatLngRect.Full.Area());
        Assert2.DoubleEqual(S2.M_PI_2, RectFromDegrees(0, 0, 90, 90).Area());
    }

    [Fact]
    internal void Test_S2LatLngRect_GetCentroid()
    {

        // Empty and full rectangles.
        Assert.Equal(new S2Point(), S2LatLngRect.Empty.Centroid());
        Assert.True(S2LatLngRect.Full.Centroid().Norm() <= 1e-15);

        // Rectangles that cover the full longitude range.
        for (int i = 0; i < 100; ++i)
        {
            double lat1 = S2Testing.Random.UniformDouble(-S2.M_PI_2, S2.M_PI_2);
            double lat2 = S2Testing.Random.UniformDouble(-S2.M_PI_2, S2.M_PI_2);
            S2LatLngRect r = new(R1Interval.FromPointPair(lat1, lat2), S1Interval.Full);
            S2Point centroid = r.Centroid();
            Assert2.Near(0.5 * (Math.Sin(lat1) + Math.Sin(lat2)) * r.Area(), centroid.Z, S2.DoubleError);
            Assert.True(new R2Point(centroid.X, centroid.Y).GetNorm() <= 1e-15);
        }

        // Rectangles that cover the full latitude range.
        for (int i = 0; i < 100; ++i)
        {
            double lng1 = S2Testing.Random.UniformDouble(-Math.PI, Math.PI);
            double lng2 = S2Testing.Random.UniformDouble(-Math.PI, Math.PI);
            S2LatLngRect r = new(S2LatLngRect.FullLat,
                           S1Interval.FromPointPair(lng1, lng2));
            S2Point centroid = r.Centroid();
            Assert.True(Math.Abs(centroid.Z) <= 1e-15);
            Assert2.Near(r.Lng.GetCenter(), new S2LatLng(centroid).LngRadians, S2.DoubleError);
            double alpha = 0.5 * r.Lng.GetLength();
            // TODO(Alas): the next Assert fails sometimes
            Assert2.Near(0.25 * Math.PI * Math.Sin(alpha) / alpha * r.Area(),
                        new R2Point(centroid.X, centroid.Y).GetNorm(), S2.DoubleError);
        }

        // Finally, verify that when a rectangle is recursively split into pieces,
        // the centroids of the pieces add to give the centroid of their parent.
        // To make the code simpler we avoid rectangles that cross the 180 degree
        // line of longitude.
        TestCentroidSplitting(
            new S2LatLngRect(S2LatLngRect.FullLat, new S1Interval(-3.14, 3.14)),
            10 /*splits_left*/);
    }

    [Fact]
    internal void Test_S2LatLngRect_GetDistanceOverlapping()
    {
        // Check pairs of rectangles that overlap: (should all return 0):
        S2LatLngRect a = RectFromDegrees(0, 0, 2, 2);
        S2LatLngRect b = PointRectFromDegrees(0, 0);
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(a));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(b));
        Assert.Equal(S1Angle.FromRadians(0), b.GetDistance(b));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(S2LatLng.FromDegrees(0, 0)));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(RectFromDegrees(0, 1, 2, 3)));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(RectFromDegrees(0, 2, 2, 4)));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(RectFromDegrees(1, 0, 3, 2)));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(RectFromDegrees(2, 0, 4, 2)));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(RectFromDegrees(1, 1, 3, 3)));
        Assert.Equal(S1Angle.FromRadians(0), a.GetDistance(RectFromDegrees(2, 2, 4, 4)));
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDistanceRectVsPoint()
    {
        // Rect that spans 180.
        S2LatLngRect a = RectFromDegrees(-1, -1, 2, 1);
        VerifyGetDistance(a, PointRectFromDegrees(-2, -1));
        VerifyGetDistance(a, PointRectFromDegrees(1, 2));

        VerifyGetDistance(PointRectFromDegrees(-2, -1), a);
        VerifyGetDistance(PointRectFromDegrees(1, 2), a);

        VerifyGetRectPointDistance(a, S2LatLng.FromDegrees(-2, -1));
        VerifyGetRectPointDistance(a, S2LatLng.FromDegrees(1, 2));

        // Tests near the north pole.
        S2LatLngRect b = RectFromDegrees(86, 0, 88, 2);
        VerifyGetDistance(b, PointRectFromDegrees(87, 3));
        VerifyGetDistance(b, PointRectFromDegrees(87, -1));
        VerifyGetDistance(b, PointRectFromDegrees(89, 1));
        VerifyGetDistance(b, PointRectFromDegrees(89, 181));
        VerifyGetDistance(b, PointRectFromDegrees(85, 1));
        VerifyGetDistance(b, PointRectFromDegrees(85, 181));
        VerifyGetDistance(b, PointRectFromDegrees(90, 0));

        VerifyGetDistance(PointRectFromDegrees(87, 3), b);
        VerifyGetDistance(PointRectFromDegrees(87, -1), b);
        VerifyGetDistance(PointRectFromDegrees(89, 1), b);
        VerifyGetDistance(PointRectFromDegrees(89, 181), b);
        VerifyGetDistance(PointRectFromDegrees(85, 1), b);
        VerifyGetDistance(PointRectFromDegrees(85, 181), b);
        VerifyGetDistance(PointRectFromDegrees(90, 0), b);

        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(87, 3));
        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(87, -1));
        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(89, 1));
        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(89, 181));
        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(85, 1));
        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(85, 181));
        VerifyGetRectPointDistance(b, S2LatLng.FromDegrees(90, 0));

        // Rect that touches the north pole.
        S2LatLngRect c = RectFromDegrees(88, 0, 90, 2);
        VerifyGetDistance(c, PointRectFromDegrees(89, 3));
        VerifyGetDistance(c, PointRectFromDegrees(89, 90));
        VerifyGetDistance(c, PointRectFromDegrees(89, 181));
        VerifyGetDistance(PointRectFromDegrees(89, 3), c);
        VerifyGetDistance(PointRectFromDegrees(89, 90), c);
        VerifyGetDistance(PointRectFromDegrees(89, 181), c);
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDistanceRectVsRect()
    {
        // Rect that spans 180.
        S2LatLngRect a = RectFromDegrees(-1, -1, 2, 1);
        VerifyGetDistance(a, RectFromDegrees(0, 2, 1, 3));
        VerifyGetDistance(a, RectFromDegrees(-2, -3, -1, -2));

        // Tests near the south pole.
        S2LatLngRect b = RectFromDegrees(-87, 0, -85, 3);
        VerifyGetDistance(b, RectFromDegrees(-89, 1, -88, 2));
        VerifyGetDistance(b, RectFromDegrees(-84, 1, -83, 2));
        VerifyGetDistance(b, RectFromDegrees(-88, 90, -86, 91));
        VerifyGetDistance(b, RectFromDegrees(-84, -91, -83, -90));
        VerifyGetDistance(b, RectFromDegrees(-90, 181, -89, 182));
        VerifyGetDistance(b, RectFromDegrees(-84, 181, -83, 182));
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDistanceRandomPairs()
    {
        // Test random pairs.
        for (int i = 0; i < 10000; ++i)
        {
            S2LatLngRect a =
                S2LatLngRect.FromPointPair(new S2LatLng(S2Testing.RandomPoint()),
                                            new S2LatLng(S2Testing.RandomPoint()));
            S2LatLngRect b =
                S2LatLngRect.FromPointPair(new S2LatLng(S2Testing.RandomPoint()),
                                            new S2LatLng(S2Testing.RandomPoint()));
            VerifyGetDistance(a, b);


            S2LatLng c = new(S2Testing.RandomPoint());
            VerifyGetRectPointDistance(a, c);
            VerifyGetRectPointDistance(b, c);
        }
    }

    [Fact]
    internal void Test_S2LatLngRect_GetDirectedHausdorffDistanceRandomPairs()
    {
        // Test random pairs.
        int kIters = 1000;
        for (int i = 0; i < kIters; ++i)
        {
            S2LatLngRect a =
                S2LatLngRect.FromPointPair(new S2LatLng(S2Testing.RandomPoint()),
                                            new S2LatLng(S2Testing.RandomPoint()));
            S2LatLngRect b =
                S2LatLngRect.FromPointPair(new S2LatLng(S2Testing.RandomPoint()),
                                            new S2LatLng(S2Testing.RandomPoint()));
            // a and b are *minimum* bounding rectangles of two random points, in
            // particular, their Voronoi diagrams are always of the same topology. We
            // take the "complements" of a and b for more thorough testing.
            S2LatLngRect a2 = new(a.Lat, a.Lng.Complement());
            S2LatLngRect b2 = new(b.Lat, b.Lng.Complement());

            // Note that "a" and "b" come from the same distribution, so there is no
            // need to test pairs such as (b, a), (b, a2), etc.
            VerifyGetDirectedHausdorffDistance(a, b);
            VerifyGetDirectedHausdorffDistance(a, b2);
            VerifyGetDirectedHausdorffDistance(a2, b);
            VerifyGetDirectedHausdorffDistance(a2, b2);
        }
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDirectedHausdorffDistanceContained()
    {
        // Caller rect is contained in callee rect. Should return 0.
        S2LatLngRect a = RectFromDegrees(-10, 20, -5, 90);
        Assert.Equal(S1Angle.FromRadians(0),
                  a.GetDirectedHausdorffDistance(RectFromDegrees(-10, 20, -5, 90)));
        Assert.Equal(S1Angle.FromRadians(0),
                  a.GetDirectedHausdorffDistance(RectFromDegrees(-10, 19, -5, 91)));
        Assert.Equal(S1Angle.FromRadians(0),
                  a.GetDirectedHausdorffDistance(RectFromDegrees(-11, 20, -4, 90)));
        Assert.Equal(S1Angle.FromRadians(0),
                  a.GetDirectedHausdorffDistance(RectFromDegrees(-11, 19, -4, 91)));
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDirectHausdorffDistancePointToRect()
    {
        // The Hausdorff distance from a point to a rect should be the same as its
        // distance to the rect.
        S2LatLngRect a1 = PointRectFromDegrees(5, 8);
        S2LatLngRect a2 = PointRectFromDegrees(90, 10);  // north pole

        S2LatLngRect b = RectFromDegrees(-85, -50, -80, 10);
        Assert2.DoubleEqual(a1.GetDirectedHausdorffDistance(b).Radians,
                         a1.GetDistance(b).Radians);
        Assert2.DoubleEqual(a2.GetDirectedHausdorffDistance(b).Radians,
                         a2.GetDistance(b).Radians);

        b = RectFromDegrees(4, -10, 80, 10);
        Assert2.DoubleEqual(a1.GetDirectedHausdorffDistance(b).Radians,
                         a1.GetDistance(b).Radians);
        Assert2.DoubleEqual(a2.GetDirectedHausdorffDistance(b).Radians,
                         a2.GetDistance(b).Radians);

        b = RectFromDegrees(70, 170, 80, -170);
        Assert2.DoubleEqual(a1.GetDirectedHausdorffDistance(b).Radians,
                         a1.GetDistance(b).Radians);
        Assert2.DoubleEqual(a2.GetDirectedHausdorffDistance(b).Radians,
                         a2.GetDistance(b).Radians);
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDirectedHausdorffDistanceRectToPoint()
    {
        S2LatLngRect a = RectFromDegrees(1, -8, 10, 20);
        VerifyGetDirectedHausdorffDistance(a, PointRectFromDegrees(5, 8));
        VerifyGetDirectedHausdorffDistance(a, PointRectFromDegrees(-6, -100));
        // south pole
        VerifyGetDirectedHausdorffDistance(a, PointRectFromDegrees(-90, -20));
        // north pole
        VerifyGetDirectedHausdorffDistance(a, PointRectFromDegrees(90, 0));
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDirectedHausdorffDistanceRectToRectNearPole()
    {
        // Tests near south pole.
        S2LatLngRect a = RectFromDegrees(-87, 0, -85, 3);
        VerifyGetDirectedHausdorffDistance(a, RectFromDegrees(-89, 1, -88, 2));
        VerifyGetDirectedHausdorffDistance(a, RectFromDegrees(-84, 1, -83, 2));
        VerifyGetDirectedHausdorffDistance(a, RectFromDegrees(-88, 90, -86, 91));
        VerifyGetDirectedHausdorffDistance(a, RectFromDegrees(-84, -91, -83, -90));
        VerifyGetDirectedHausdorffDistance(a, RectFromDegrees(-90, 181, -89, 182));
        VerifyGetDirectedHausdorffDistance(a, RectFromDegrees(-84, 181, -83, 182));
    }
    [Fact]
    internal void Test_S2LatLngRect_GetDirectedHausdorffDistanceRectToRectDegenerateCases()
    {
        // Rectangles that contain poles.
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(0, 10, 90, 20), RectFromDegrees(-4, -10, 4, 0));
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(-4, -10, 4, 0), RectFromDegrees(0, 10, 90, 20));

        // Two rectangles share same or complement longitudinal intervals.
        S2LatLngRect a = RectFromDegrees(-50, -10, 50, 10);
        S2LatLngRect b = RectFromDegrees(30, -10, 60, 10);
        VerifyGetDirectedHausdorffDistance(a, b);
        S2LatLngRect c = new(a.Lat, a.Lng.Complement());
        VerifyGetDirectedHausdorffDistance(c, b);

        // rectangle a touches b_opposite_lng.
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(10, 170, 30, 180), RectFromDegrees(-50, -10, 50, 10));
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(10, -180, 30, -170), RectFromDegrees(-50, -10, 50, 10));

        // rectangle b's Voronoi diagram is degenerate (lng interval spans 180
        // degrees), and a touches the degenerate Voronoi vertex.
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(-30, 170, 30, 180), RectFromDegrees(-10, -90, 10, 90));
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(-30, -180, 30, -170), RectFromDegrees(-10, -90, 10, 90));

        // rectangle a touches a voronoi vertex of rectangle b.
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(-20, 105, 20, 110), RectFromDegrees(-30, 5, 30, 15));
        VerifyGetDirectedHausdorffDistance(
            RectFromDegrees(-20, 95, 20, 105), RectFromDegrees(-30, 5, 30, 15));
    }

    private static S2LatLngRect RectFromDegrees(double lat_lo, double lng_lo, double lat_hi, double lng_hi)
    {
        // Convenience method to construct a rectangle.  This method is
        // intentionally *not* in the S2LatLngRect interface because the
        // argument order is ambiguous, but hopefully it's not too confusing
        // within the context of this unit test.

        return new S2LatLngRect(S2LatLng.FromDegrees(lat_lo, lng_lo).Normalized(), S2LatLng.FromDegrees(lat_hi, lng_hi).Normalized());
    }

    private static void TestIntervalOps(S2LatLngRect x, S2LatLngRect y, string expected_relation, S2LatLngRect expected_union, S2LatLngRect expected_intersection)
    {
        // Test all of the interval operations on the given pair of intervals.
        // "expected_relation" is a sequence of "T" and "F" characters corresponding
        // to the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.

        Assert.Equal(x.Contains(y), expected_relation[0] == 'T');
        Assert.Equal(x.InteriorContains(y), expected_relation[1] == 'T');
        Assert.Equal(x.Intersects(y), expected_relation[2] == 'T');
        Assert.Equal(x.InteriorIntersects(y), expected_relation[3] == 'T');

        Assert.Equal(x.Contains(y), x.Union(y) == x);
        Assert.Equal(x.Intersects(y), !x.Intersection(y).IsEmpty());

        Assert.Equal(x.Union(y), expected_union);
        Assert.Equal(x.Intersection(y), expected_intersection);

        if (y.Size() == S2LatLng.FromRadians(0, 0))
        {
            S2LatLngRect r = x;
            r.AddPoint(y.Lo());
            Assert.Equal(r, expected_union);
        }
    }

    // This function assumes that GetDirectedHausdorffDistance() always returns
    // a distance from some point in a to b. So the function mainly tests whether
    // the returned distance is large enough, and only does a weak test on whether
    // it is small enough.
    private static void VerifyGetDirectedHausdorffDistance(S2LatLngRect a, S2LatLngRect b)
    {
        S1Angle hausdorff_distance = a.GetDirectedHausdorffDistance(b);

        const double kResolution = 0.1;
        S1Angle max_distance = S1Angle.Zero;

        int sample_size_on_lat =
            (int)(a.Lat.GetLength() / kResolution) + 1;
        int sample_size_on_lng =
            (int)(a.Lng.GetLength() / kResolution) + 1;
        double delta_on_lat = a.Lat.GetLength() / sample_size_on_lat;
        double delta_on_lng = a.Lng.GetLength() / sample_size_on_lng;

        double lng = a.Lng.Lo;
        for (int i = 0; i <= sample_size_on_lng; ++i, lng += delta_on_lng)
        {
            double lat = a.Lat.Lo;
            for (int j = 0; j <= sample_size_on_lat; ++j, lat += delta_on_lat)
            {
                S2LatLng latlng = S2LatLng.FromRadians(lat, lng).Normalized();
                S1Angle distance_to_b = b.GetDistance(latlng);

                if (distance_to_b >= max_distance)
                {
                    max_distance = distance_to_b;
                }
            }
        }

        Assert.True(max_distance.Radians <= hausdorff_distance.Radians + 1e-10);
        Assert.True(max_distance.Radians >= hausdorff_distance.Radians - kResolution);
    }

    // Recursively verify that when a rectangle is split into two pieces, the
    // centroids of the children sum to give the centroid of the parent.
    static void TestCentroidSplitting(S2LatLngRect r, int splits_left)
    {
        S2LatLngRect child0, child1;
        if (S2Testing.Random.OneIn(2))
        {
            double lat = S2Testing.Random.UniformDouble(r.Lat.Lo, r.Lat.Hi);
            child0 = new S2LatLngRect(new R1Interval(r.Lat.Lo, lat), r.Lng);
            child1 = new S2LatLngRect(new R1Interval(lat, r.Lat.Hi), r.Lng);
        }
        else
        {
            Assert.True(r.Lng.Lo <= r.Lng.Hi);
            double lng = S2Testing.Random.UniformDouble(r.Lng.Lo, r.Lng.Hi);
            child0 = new S2LatLngRect(r.Lat, new S1Interval(r.Lng.Lo, lng));
            child1 = new S2LatLngRect(r.Lat, new S1Interval(lng, r.Lng.Hi));
        }
        Assert.True((r.Centroid() - child0.Centroid() - child1.Centroid()).Norm() <= 1e-15);
        if (splits_left > 0)
        {
            TestCentroidSplitting(child0, splits_left - 1);
            TestCentroidSplitting(child1, splits_left - 1);
        }
    }

    private static void TestCellOps(S2LatLngRect r, S2Cell cell, int level)
    {
        // Test the relationship between the given rectangle and cell:
        // 0 == no intersection, 1 == MayIntersect, 2 == Intersects,
        // 3 == Vertex Containment, 4 == Contains

        bool vertex_contained = false;
        for (int i = 0; i < 4; ++i)
        {
            if (r.Contains(cell.VertexRaw(i)) ||
                (!r.IsEmpty() && cell.Contains(r.Vertex(i).ToPoint())))
                vertex_contained = true;
        }
        Assert.Equal(r.MayIntersect(cell), level >= 1);
        Assert.Equal(r.Intersects(cell), level >= 2);
        Assert.Equal(vertex_contained, level >= 3);
        Assert.Equal(r.Contains(cell), level >= 4);
    }
    // Returns the minimum distance from X to the latitude line segment defined by
    // the given latitude and longitude interval.
    private static S1Angle GetDistance(S2LatLng x, S1Angle lat, S1Interval interval)
    {
        Assert.True(x.IsValid());
        Assert.True(interval.IsValid());

        // Is X inside the longitude interval?
        if (interval.Contains(x.LngRadians))
            return S1Angle.FromRadians((x.Lat() - lat).Abs());

        // Return the distance to the closer endpoint.
        return new[] { x.GetDistance(new S2LatLng(lat, S1Angle.FromRadians(interval.Lo))),
                       x.GetDistance(new S2LatLng(lat, S1Angle.FromRadians(interval.Hi))) }.Min();
    }

    private static S1Angle BruteForceDistance(S2LatLngRect a, S2LatLngRect b)
    {
        if (a.Intersects(b))
            return S1Angle.FromRadians(0);

        // Compare every point in 'a' against every latitude edge and longitude edge
        // in 'b', and vice-versa, for a total of 16 point-vs-latitude-edge tests and
        // 16 point-vs-longitude-edge tests.
        var pnt_a = new S2LatLng[4];
        var pnt_b = new S2LatLng[4];
        pnt_a[0] = new S2LatLng(a.LatLo(), a.LngLo());
        pnt_a[1] = new S2LatLng(a.LatLo(), a.LngHi());
        pnt_a[2] = new S2LatLng(a.LatHi(), a.LngHi());
        pnt_a[3] = new S2LatLng(a.LatHi(), a.LngLo());
        pnt_b[0] = new S2LatLng(b.LatLo(), b.LngLo());
        pnt_b[1] = new S2LatLng(b.LatLo(), b.LngHi());
        pnt_b[2] = new S2LatLng(b.LatHi(), b.LngHi());
        pnt_b[3] = new S2LatLng(b.LatHi(), b.LngLo());

        // Make arrays containing the lo/hi latitudes and the lo/hi longitude edges.
        var lat_a = new S1Angle[2] { a.LatLo(), a.LatHi() };
        var lat_b = new S1Angle[2] { b.LatLo(), b.LatHi() };
        var lng_edge_a = new S2Point[][] { new S2Point[] { pnt_a[0].ToPoint(), pnt_a[3].ToPoint() },
                                new S2Point[]{ pnt_a[1].ToPoint(), pnt_a[2].ToPoint() } };
        var lng_edge_b = new S2Point[][] { new S2Point[] { pnt_b[0].ToPoint(), pnt_b[3].ToPoint() },
                                new S2Point[]{ pnt_b[1].ToPoint(), pnt_b[2].ToPoint() } };

        S1Angle min_distance = S1Angle.FromDegrees(180.0);
        for (int i = 0; i < 4; ++i)
        {
            // For each point in a and b.
            var current_a = pnt_a[i];
            var current_b = pnt_b[i];

            for (int j = 0; j < 2; ++j)
            {
                // Get distances to latitude and longitude edges.
                S1Angle a_to_lat = GetDistance(current_a, lat_b[j], b.Lng);
                S1Angle b_to_lat = GetDistance(current_b, lat_a[j], a.Lng);
                S1Angle a_to_lng = S2.GetDistance(current_a.ToPoint(), lng_edge_b[j][0], lng_edge_b[j][1]);
                S1Angle b_to_lng = S2.GetDistance(current_b.ToPoint(), lng_edge_a[j][0], lng_edge_a[j][1]);

                min_distance = new[] { min_distance, a_to_lat, b_to_lat, a_to_lng, b_to_lng }.Min();
            }
        }
        return min_distance;
    }

    private static S1Angle BruteForceRectPointDistance(S2LatLngRect a, S2LatLng b)
    {
        if (a.Contains(b))
        {
            return S1Angle.FromRadians(0);
        }

        S1Angle b_to_lo_lat = GetDistance(b, a.LatLo(), a.Lng);
        S1Angle b_to_hi_lat = GetDistance(b, a.LatHi(), a.Lng);
        S1Angle b_to_lo_lng = S2.GetDistance(b.ToPoint(), new S2LatLng(a.LatLo(), a.LngLo()).ToPoint(), new S2LatLng(a.LatHi(), a.LngLo()).ToPoint());
        S1Angle b_to_hi_lng = S2.GetDistance(b.ToPoint(), new S2LatLng(a.LatLo(), a.LngHi()).ToPoint(), new S2LatLng(a.LatHi(), a.LngHi()).ToPoint());
        return new[] { b_to_lo_lat, b_to_hi_lat, b_to_lo_lng, b_to_hi_lng }.Min();
    }

    // This method verifies a.GetDistance(b) by comparing its result against a
    // brute-force implementation. The correctness of the brute-force version is
    // much easier to verify by inspection.
    private static void VerifyGetDistance(S2LatLngRect a, S2LatLngRect b)
    {
        S1Angle distance1 = BruteForceDistance(a, b);
        S1Angle distance2 = a.GetDistance(b);
        Assert2.Near(distance1.Radians - distance2.Radians, 0, 1e-10);
    }

    private static S2LatLngRect PointRectFromDegrees(double lat, double lng)
    {
        return S2LatLngRect.FromPoint(
            S2LatLng.FromDegrees(lat, lng).Normalized());
    }

    // This method verifies a.GetDistance(b), where b is a S2LatLng, by comparing
    // its result against a.GetDistance(c), c being the point rectangle created
    // from b.
    private static void VerifyGetRectPointDistance(S2LatLngRect a, S2LatLng p)
    {
        S1Angle distance1 = BruteForceRectPointDistance(a, p.Normalized());
        S1Angle distance2 = a.GetDistance(p.Normalized());
        Assert2.Near(Math.Abs(distance1.Radians - distance2.Radians), 0, 1e-10);
    }
}
