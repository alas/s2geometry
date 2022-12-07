namespace S2Geometry;

public class S2EdgeTessellatorTests
{
    // The interpolation parameter actually used in the .cc file.
    private const double kBestFraction = 0.31215691082248312;
    private static readonly S1ChordAngle kMaxInterpolationError = new(S1Angle.FromRadians(1e-14));

    // When there are longitudes greater than 180 degrees due to wrapping, the
    // combination of projecting and unprojecting an S2Point can have slightly more
    // error than is allowed by ApproxEquals.
    private static readonly S1Angle kMaxProjError = S1Angle.FromRadians(2e-15);
    private readonly ITestOutputHelper _logger;

    internal S2EdgeTessellatorTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedNoTessellation() {
        PlateCarreeProjection proj=new(180);
        S2EdgeTessellator tess = new(proj, S1Angle.FromDegrees(0.01));
        List<R2Point> vertices = new();
        tess.AppendProjected(new S2Point(1, 0, 0), new S2Point(0, 1, 0), vertices);
        Assert.Equal(2, vertices.Count);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedNoTessellation() {
        PlateCarreeProjection proj = new(180);
        S2EdgeTessellator tess = new(proj, S1Angle.FromDegrees(0.01));
        List<S2Point> vertices = new();
        tess.AppendUnprojected(new R2Point(0, 30), new R2Point(0, 50), vertices);
        Assert.Equal(2, vertices.Count);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedWrapping() {
        // This tests that a projected edge that crosses the 180 degree meridian
        // goes the "short way" around the sphere.

        PlateCarreeProjection proj = new(180);
        S2EdgeTessellator tess = new(proj, S1Angle.FromDegrees(0.01));
        List<S2Point> vertices = new();
        tess.AppendUnprojected(new R2Point(-170, 0), new R2Point(170, 80), vertices);
        foreach (var v in vertices) {
            Assert.True(Math.Abs(S2LatLng.Longitude(v).GetDegrees()) >= 170);
        }
    }

    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedWrapping() {
        // This tests projecting a geodesic edge that crosses the 180 degree
        // meridian.  This results in a set of vertices that may be non-canonical
        // (i.e., absolute longitudes greater than 180 degrees) but that don't have
        // any sudden jumps in value, which is convenient for interpolating them.
        PlateCarreeProjection proj = new(180);
        S2EdgeTessellator tess = new(proj, S1Angle.FromDegrees(0.01));
        List<R2Point> vertices = new();
        tess.AppendProjected(S2LatLng.FromDegrees(0, -170).ToPoint(),
                             S2LatLng.FromDegrees(0, 170).ToPoint(), vertices);
        foreach (var v in vertices) {
            Assert.True(v.X <= -170);
        }
    }

    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedWrappingMultipleCrossings() {
        // Tests an edge chain that crosses the 180 degree meridian multiple times.
        // Note that due to coordinate wrapping, the last vertex of one edge may not
        // exactly match the first edge of the next edge after unprojection.
        PlateCarreeProjection proj = new(180);
        S2EdgeTessellator tess = new(proj, S1Angle.FromDegrees(0.01));
        List<S2Point> vertices = new();
        for (double lat = 1; lat <= 60; ++lat) {
            tess.AppendUnprojected(new R2Point(180 - 0.03 * lat, lat),
                                   new R2Point(-180 + 0.07 * lat, lat), vertices);
            tess.AppendUnprojected(new R2Point(-180 + 0.07 * lat, lat),
                                   new R2Point(180 - 0.03 * (lat + 1), lat + 1), vertices);
        }
        foreach (var v in vertices) {
            Assert.True(Math.Abs(S2LatLng.Longitude(v).GetDegrees()) >= 175);
        }
    }

    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedWrappingMultipleCrossings() {
        // The following loop crosses the 180 degree meridian four times (twice in
        // each direction).
        var loop = ParsePointsOrDie("0:160, 0:-40, 0:120, 0:-80, 10:120, " +
                                     "10:-40, 0:160");
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S1Angle.FromE7(1);
        S2EdgeTessellator tess = new(proj, tolerance);
        List<R2Point> vertices = new();
        for (int i = 0; i + 1 < loop.Count; ++i) {
            tess.AppendProjected(loop[i], loop[i + 1], vertices);
        }
        Assert.Equal(vertices.First(), vertices.Last());

        // Note that the R2Point coordinates are in (lng, lat) order.
        double min_lng = vertices[0].X;
        double max_lng = vertices[0].X;
        foreach (R2Point v in vertices) {
            min_lng = Math.Min(min_lng, v.X);
            max_lng = Math.Max(max_lng, v.X);
        }
        Assert.Equal(160, min_lng);
        Assert.Equal(640, max_lng);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_InfiniteRecursionBug() {
        PlateCarreeProjection proj = new(180);
        S1Angle kOneMicron = S1Angle.FromRadians(1e-6 / 6371.0);
        S2EdgeTessellator tess = new(proj, kOneMicron);
        List<R2Point> vertices = new();
        tess.AppendProjected(S2LatLng.FromDegrees(3, 21).ToPoint(),
                             S2LatLng.FromDegrees(1, -159).ToPoint(), vertices);
        Assert.Equal(36, vertices.Count);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedAccuracy() {
        MercatorProjection proj = new(180);
        S1Angle tolerance = S1Angle.FromDegrees(1e-5);
        R2Point pa = new(0, 0), pb = new(89.999999, 179);
        Stats stats = TestUnprojected(proj, tolerance, pa, pb, true);
        Assert.True(stats.Max <= 1.0);
    }

    // Repro case for b/110719057.
    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedAccuracyCrossEquator() {
        MercatorProjection proj = new(180);
        S1Angle tolerance = S1Angle.FromDegrees(1e-5);
        R2Point pa = new(-10, -10), pb = new(10, 10);
        Stats stats = TestUnprojected(proj, tolerance, pa, pb, true);
        Assert.True(stats.Max < 1.0);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedAccuracy() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance=S1Angle.FromE7(1);
        S2Point a = S2LatLng.FromDegrees(-89.999, -170).ToPoint();
        S2Point b = S2LatLng.FromDegrees(50, 100).ToPoint();
        Stats stats = TestProjected(proj, tolerance, a, b, true);
        Assert.True(stats.Max <= 1.0);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedAccuracyMidpointEquator() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S2Testing.MetersToAngle(1);
        R2Point a = new(80, 50), b = new(-80, -50);
        Stats stats = TestUnprojected(proj, tolerance, a, b, true);
        Assert.True(stats.Max <= 1.0);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedAccuracyMidpointEquator() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S2Testing.MetersToAngle(1);
        S2Point a = S2LatLng.FromDegrees(50, 80).ToPoint();
        S2Point b = S2LatLng.FromDegrees(-50, -80).ToPoint();
        Stats stats = TestProjected(proj, tolerance, a, b, true);
        Assert.True(stats.Max <= 1.0);
    }

    // Repro case for b/110719057.
    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedAccuracyCrossEquator() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S1Angle.FromE7(1);
        S2Point a = S2LatLng.FromDegrees(-20, -20).ToPoint();
        S2Point b = S2LatLng.FromDegrees(20, 20).ToPoint();
        Stats stats = TestProjected(proj, tolerance, a, b, true);
        Assert.True(stats.Max < 1.0);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedAccuracySeattleToNewYork() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S2Testing.MetersToAngle(1);
        S2Point seattle = S2LatLng.FromDegrees(47.6062, -122.3321).ToPoint();
        S2Point newyork = S2LatLng.FromDegrees(40.7128, -74.0059).ToPoint();
        Stats stats = TestProjected(proj, tolerance, seattle, newyork, true);
        Assert.True(stats.Max <= 1.0);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_MaxEdgeErrorPlateCarree() {
        PlateCarreeProjection proj = new(180);
        // Uncomment to test some nearby parameter values.
        // TestEdgeError(proj, 0.311);
        TestEdgeError(proj, kBestFraction);
        // TestEdgeError(proj, 0.313);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_MaxEdgeErrorMercator() {
        MercatorProjection proj = new(180);
        // Uncomment to test some nearby parameter values.
        // TestEdgeError(proj, 0.311);
        TestEdgeError(proj, kBestFraction);
        // TestEdgeError(proj, 0.313);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_RandomEdgesPlateCarree() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S2Testing.MetersToAngle(100);
        TestRandomEdges(proj, tolerance);
    }

    [Fact]
    internal void Test_S2EdgeTessellator_RandomEdgesMercator() {
        MercatorProjection proj = new(180);
        S1Angle tolerance = S2Testing.MetersToAngle(100);
        TestRandomEdges(proj, tolerance);
    }

    // TODO(ericv): Superceded by random edge tests above, remove?
    [Fact]
    internal void Test_S2EdgeTessellator_UnprojectedAccuracyRandomCheck() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S1Angle.FromDegrees(1e-3);
        int kIters = 5000;
#if DEBUG
        kIters = 250;
#endif
        for (int i = 0; i < kIters; ++i) {
            S2Testing.Random.Reset(i);
            double alat = S2Testing.Random.UniformDouble(-89.99, 89.99);
            double blat = S2Testing.Random.UniformDouble(-89.99, 89.99);
            double blon = S2Testing.Random.UniformDouble(0, 179);

            R2Point pa = new(0, alat), pb = new(blon, blat);
            Stats stats = TestUnprojected(proj, tolerance, pa, pb, false);
            Assert.True(stats.Max < 1.0);
        }
    }

    // XXX(ericv): Superceded by random edge tests above, remove?
    [Fact]
    internal void Test_S2EdgeTessellator_ProjectedAccuracyRandomCheck() {
        PlateCarreeProjection proj = new(180);
        S1Angle tolerance = S1Angle.FromDegrees(1e-3);
        int kIters = 5000;
#if DEBUG
        kIters = 250;
#endif

        for (int i = 0; i < kIters; ++i) {
            S2Testing.Random.Reset(i);
            double alat = S2Testing.Random.UniformDouble(-89.99, 89.99);
            double blat = S2Testing.Random.UniformDouble(-89.99, 89.99);
            double blon = S2Testing.Random.UniformDouble(-180, 180);

            S2Point a = S2LatLng.FromDegrees(alat, 0).ToPoint();
            S2Point b = S2LatLng.FromDegrees(blat, blon).ToPoint();
            Stats stats = TestProjected(proj, tolerance, a, b, false);
            Assert.True(stats.Max < 1.0);
        }
    }

    private static S1Angle GetMaxDistance(Projection proj,
                           R2Point px, S2Point x,
                           R2Point py, S2Point y,
                           DistType dist_type = DistType.GEOMETRIC)
    {
        // Step along the projected edge at a fine resolution and keep track of the
        // maximum distance of any point to the current geodesic edge.
        int kNumSteps = 100;
        S1ChordAngle max_dist = S1ChordAngle.Zero;
        for (double f = 0.5 / kNumSteps; f < 1.0; f += 1.0 / kNumSteps)
        {
            S1ChordAngle dist = S1ChordAngle.Infinity;
            S2Point p = proj.Unproject(Projection.Interpolate(f, px, py));
            if (dist_type == DistType.GEOMETRIC)
            {
                S2.UpdateMinDistance(p, x, y, ref dist);
            }
            else
            {
                Assert.True(dist_type == DistType.PARAMETRIC);
                dist = new S1ChordAngle(p, S2.Interpolate(f, x, y));
            }
            if (dist > max_dist) max_dist = dist;
        }
        // Ensure that the maximum distance estimate is a lower bound, not an upper
        // bound, since we only want to record a failure of the distance estimation
        // algorithm if the number it returns is definitely too small.
        return max_dist.PlusError(-S2.GetUpdateMinDistanceMaxError(max_dist)).ToAngle();
    }

    // Converts a projected edge to a sequence of geodesic edges and verifies that
    // the result satisfies the given tolerance.
    private Stats TestUnprojected(Projection proj, S1Angle tolerance, R2Point pa, R2Point pb_in, bool log_stats)
    {
        S2EdgeTessellator tess=new(proj, tolerance);
        List<S2Point> vertices = new();
        tess.AppendUnprojected(pa, pb_in, vertices);
        R2Point pb = proj.WrapDestination(pa, pb_in);
        Assert.True(new S1Angle(proj.Unproject(pa), vertices.First()) <= kMaxProjError);
        Assert.True(new S1Angle(proj.Unproject(pb), vertices.Last()) <= kMaxProjError);
        Stats stats = new();
        if (pa == pb)
        {
            Assert.Single(vertices);
            return stats;
        }
        // Precompute the normal to the projected edge.
        R2Point norm = R2Point.Normalize((pb - pa).GetOrtho());
        S2Point x = vertices[0];
        R2Point px = proj.Project(x);
        for (int i = 1; i < vertices.Count; ++i)
        {
            S2Point y = vertices[i];
            R2Point py = proj.WrapDestination(px, proj.Project(y));
            // Check that every vertex is on the projected edge.
            Assert.True((py - pa).DotProd(norm) <= 1e-14 * py.GetNorm());
            stats.Tally(GetMaxDistance(proj, px, x, py, y) / tolerance);
            x = y;
            px = py;
        }
        if (log_stats)
        {
            _logger.WriteLine($"{pa} to {pb}: {vertices.Count} vertices, {stats}");
        }
        return stats;
    }

    // Converts a geodesic edge to a sequence of projected edges and verifies that
    // the result satisfies the given tolerance.
    private Stats TestProjected(Projection proj, S1Angle tolerance, S2Point a, S2Point b, bool log_stats)
    {
        S2EdgeTessellator tess = new(proj, tolerance);
        List<R2Point> vertices = new();
        tess.AppendProjected(a, b, vertices);
        Assert.True(new S1Angle(a, proj.Unproject(vertices.First())) <= kMaxProjError);
        Assert.True(new S1Angle(b, proj.Unproject(vertices.Last())) <= kMaxProjError);
        Stats stats = new();
        if (a == b)
        {
            Assert.Single(vertices);
            return stats;
        }
        R2Point px = vertices[0];
        S2Point x = proj.Unproject(px);
        for (int i = 1; i < vertices.Count; ++i)
        {
            R2Point py = vertices[i];
            S2Point y = proj.Unproject(py);
            // Check that every vertex is on the geodesic edge.
            Assert.True(S2.IsDistanceLess(y, a, b, kMaxInterpolationError));
            stats.Tally(GetMaxDistance(proj, px, x, py, y) / tolerance);
            x = y;
            px = py;
        }
        if (log_stats)
        {
            _logger.WriteLine($"{vertices[0]} to {px}: {vertices.Count} vertices, {stats}");
        }
        return stats;
    }


    // Given a projection, this function repeatedly chooses a pair of edge
    // endpoints and measures the true distance between the geodesic and projected
    // edges that connect those two endpoints.  It then compares this to against
    // the distance measurement algorithm used by S2EdgeTessellator, which
    // consists of measuring the point-to-point distance between the edges at each
    // of two fractions "t" and "1-t", computing the maximum of those two
    // distances, and then scaling by the constant documented in the .cc file
    // (based on the idea that the distance between the edges as a function of the
    // interpolation fraction can be accurately modeled as a cubic polynomial).
    //
    // This function is used to (1) verify that the distance estimates are always
    // conservative, (2) verify the optimality of the interpolation fraction "t",
    // and (3) estimate the amount of overtessellation that occurs for various
    // types of edges (e.g., short vs. long edges, edges that follow lines of
    // latitude or longitude, etc).
    private void TestEdgeError(Projection proj, double t)
    {
        // Here we compute how much we need to scale the error measured at the
        // chosen interpolation fraction "t" in order to bound the error along the
        // entire edge, under the assumption that the error is a convex combination
        // of E1(x) and E2(x) (see comments in the .cc file).
        double x = 1 - 2 * t;
        double dlat = Math.Sin(0.5 * S2.M_PI_4 * (1 - x));
        double dlng = Math.Sin(S2.M_PI_4 * (1 - x));
        double dsin2 = dlat * dlat + dlng * dlng * Math.Sin(S2.M_PI_4 * x) * S2.M_SQRT1_2;
        double dsin2_max = 0.5 * (1 - S2.M_SQRT1_2);
        // Note that this is the reciprocal of the value used in the .cc file!
        double kScaleFactor = Math.Max((2 * Math.Sqrt(3) / 9) / (x * (1 - x * x)),
                                        Math.Asin(Math.Sqrt(dsin2_max)) / Math.Asin(Math.Sqrt(dsin2)));

        // Keep track of the average and maximum geometric and parametric errors.
        Stats stats_g = new(), stats_p = new();
        int kIters = 100000;
#if DEBUG
        kIters = 10000;
#endif

        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Testing.Random.Reset(iter);
            S2Point a = S2Testing.RandomPoint();
            S2Point b = S2Testing.RandomPoint();
            // Uncomment to test edges longer than 90 degrees.
            if (a.DotProd(b) < -1e-14) continue;
            // Uncomment to test edges than span more than 90 degrees longitude.
            // if (a[0] * b[0] + a[1] * b[1] < 0) continue;
            // Uncomment to only test edges of a certain length.
            // b = S2EdgeDistances.InterpolateAtDistance(S1Angle.FromRadians(1e-5), a, b);
            // Uncomment to only test edges that stay in one hemisphere.
            // if (a[2] * b[2] <= 0) continue;
            R2Point pa = proj.Project(a);
            R2Point pb = proj.WrapDestination(pa, proj.Project(b));
            S1Angle max_dist_g = GetMaxDistance(proj, pa, a, pb, b, DistType.GEOMETRIC);
            // Ignore edges where the error is too small.
            if (max_dist_g <= S2EdgeTessellator.kMinTolerance()) continue;
            S1Angle max_dist_p = GetMaxDistance(proj, pa, a, pb, b, DistType.PARAMETRIC);
            if (max_dist_p <= S2EdgeTessellator.kMinTolerance()) continue;

            // Compute the estimated error bound.
            S1Angle d1 = new(S2.Interpolate(t, a, b), proj.Unproject((1 - t) * pa + t * pb));
            S1Angle d2 = new(S2.Interpolate(1 - t, a, b), proj.Unproject(t * pa + (1 - t) * pb));
            S1Angle dist = kScaleFactor * new[] { S1Angle.FromRadians(1e-300), d1, d2 }.Max();

            // Compute the ratio of the true geometric/parametric errors to the
            // estimate error bound.
            double r_g = max_dist_g / dist;
            double r_p = max_dist_p / dist;

            // Our objective is to ensure that the geometric error ratio is at most 1.
            // (The parametric ratio is computed only for analysis purposes.)
            if (r_g > 0.99999)
            {
                // Log any edges where the ratio is exceeded (or nearly so).
                _logger.WriteLine($"{pa} to {pb}: ratio = {r_g}, dist = {max_dist_g.GetDegrees()}");
            }
            stats_g.Tally(r_g);
            stats_p.Tally(r_p);
        }
        _logger.WriteLine($"t = {t}, scale = {kScaleFactor}, G[{stats_g}], P[{stats_p}]");
        Assert.True(stats_g.Max <= kScaleFactor);
    }


    // Tessellates random edges using the given projection and tolerance, and
    // verifies that the expected criteria are satisfied.
    private void TestRandomEdges(Projection proj, S1Angle tolerance) {
        int kIters = 500;
#if DEBUG
        kIters = 50;
#endif
        double max_r2 = 0, max_s2 = 0;
        for (int iter = 0; iter < kIters; ++iter) {
            S2Testing.Random.Reset(iter);
            S2Point a = S2Testing.RandomPoint();
            S2Point b = S2Testing.RandomPoint();
            max_r2 = Math.Max(max_r2, TestProjected(proj, tolerance, a, b, false).Max);
            R2Point pa = proj.Project(a);
            R2Point pb = proj.Project(b);
            max_s2 = Math.Max(max_s2, TestUnprojected(proj, tolerance, pa, pb, false).Max);
        }
        _logger.WriteLine($"max_r2 = {max_r2}, max_s2 = {max_s2}");
        Assert.True(max_r2 <= 1.0);
        Assert.True(max_s2 <= 1.0);
    }

    private class Stats
    {
        private double sum;
        private int count;
        internal double Max { get; private set; }

        internal Stats()
        {
            Max = double.NegativeInfinity;
            sum = 0; count = 0;
        }

        internal void Tally(double v)
        {
            Assert.False(double.IsNaN(v));
            Max = Math.Max(v, Max);
            sum += v;
            count += 1;
        }

        internal double Avg => sum / count;

        public override string ToString() => $"avg = {Avg}, max = {Max}";
    };

    // Determines whether the distance between the two edges is measured
    // geometrically or parameterically (see algorithm description in .cc file).
    private enum DistType { GEOMETRIC, PARAMETRIC };
}
