namespace S2Geometry;

public class S2LoopMeasuresTests
{
    // Some standard loops to use in the tests (see descriptions below).
    private readonly S2PointLoopSpan full_;
    private readonly S2PointLoopSpan v_loop_;
    private readonly S2PointLoopSpan north_hemi_;
    private readonly S2PointLoopSpan north_hemi3_;
    private readonly S2PointLoopSpan west_hemi_;
    private readonly S2PointLoopSpan east_hemi_;
    private readonly S2PointLoopSpan candy_cane_;
    private readonly S2PointLoopSpan line_triangle_;
    private readonly S2PointLoopSpan skinny_chevron_;
    private readonly S2PointLoopSpan three_leaf_clover_;
    private readonly S2PointLoopSpan tessellated_loop_;

    public S2LoopMeasuresTests()
    {
        // The full loop is represented as a loop with no vertices.
        full_ = [];

        // A degenerate loop in the shape of a "V".
        v_loop_ = ParsePointsOrDie("5:1, 0:2, 5:3, 0:2");

        // The northern hemisphere, defined using two pairs of antipodal points.
        north_hemi_ = ParsePointsOrDie("0:-180, 0:-90, 0:0, 0:90");

        // The northern hemisphere, defined using three points 120 degrees apart.
        north_hemi3_ = ParsePointsOrDie("0:-180, 0:-60, 0:60");

        // The western hemisphere, defined using two pairs of antipodal points.
        west_hemi_ = ParsePointsOrDie("0:-180, -90:0, 0:0, 90:0");

        // The eastern hemisphere, defined using two pairs of antipodal points.
        east_hemi_ = ParsePointsOrDie("90:0, 0:0, -90:0, 0:-180");

        // A spiral stripe that slightly over-wraps the equator.
        candy_cane_ = ParsePointsOrDie(
            "-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70");

        // A completely degenerate triangle along the equator that Sign()
        // considers to be CCW.
        line_triangle_ = ParsePointsOrDie("0:1, 0:2, 0:3");

        // A nearly-degenerate CCW chevron near the equator with very long sides
        // (about 80 degrees).  Its area is less than 1e-640, which is too small
        // to represent in double precision.
        skinny_chevron_ = ParsePointsOrDie("0:0, -1e-320:80, 0:1e-320, 1e-320:80");

        // A loop where the same vertex appears three times.
        three_leaf_clover_ = ParsePointsOrDie(
            "0:0, -3:3, 3:3, 0:0, 3:0, 3:-3, 0:0, -3:-3, -3:0");

        // A loop with groups of 3 or more vertices in a straight line.
        tessellated_loop_ = ParsePointsOrDie(
            "10:34, 5:34, 0:34, -10:34, -10:36, -5:36, 0:36, 10:36");
    }

    [Fact]
    internal void Test_PruneDegeneracies_AllDegeneracies()
    {
        TestPruneDegeneracies("", "");
        TestPruneDegeneracies("a", "");
        TestPruneDegeneracies("aaaaa", "");
        TestPruneDegeneracies("ab", "");
        TestPruneDegeneracies("abb", "");
        TestPruneDegeneracies("aab", "");
        TestPruneDegeneracies("aba", "");
        TestPruneDegeneracies("abba", "");
        TestPruneDegeneracies("abcb", "");
        TestPruneDegeneracies("abcba", "");
        TestPruneDegeneracies("abcdcdedefedcbcdcb", "");
    }

    [Fact]
    internal void Test_PruneDegeneracies_SomeDegeneracies()
    {
        TestPruneDegeneracies("abc", "abc");
        TestPruneDegeneracies("abca", "abc");
        TestPruneDegeneracies("abcc", "abc");
        TestPruneDegeneracies("abccaa", "abc");
        TestPruneDegeneracies("aabbcc", "abc");
        TestPruneDegeneracies("abcdedca", "abc");
        TestPruneDegeneracies("abcbabcbcdc", "abc");
        TestPruneDegeneracies("xyzabcazy", "abc");
        TestPruneDegeneracies("xxyyzzaabbccaazzyyxx", "abc");
    }

    [Fact]
    internal void Test_GetCanonicalLoopOrder_AllDegeneracies()
    {
        TestCanonicalLoopOrder("", new S2.LoopOrder(0, 1));
        TestCanonicalLoopOrder("a", new S2.LoopOrder(0, 1));
        TestCanonicalLoopOrder("aaaaa", new S2.LoopOrder(0, 1));
        TestCanonicalLoopOrder("ba", new S2.LoopOrder(1, 1));
        TestCanonicalLoopOrder("bab", new S2.LoopOrder(1, 1));
        TestCanonicalLoopOrder("cbab", new S2.LoopOrder(2, 1));
        TestCanonicalLoopOrder("bacbcab", new S2.LoopOrder(8, -1));
    }

    [Fact]
    internal void Test_GetPerimeter_Empty()
    {
        Assert.Equal(S1Angle.Zero, S2.GetPerimeter([]));
    }

    [Fact]
    internal void Test_GetPerimeter_Octant()
    {
        var loop = ParsePointsOrDie("0:0, 0:90, 90:0");
        Assert2.DoubleEqual(3 * S2.M_PI_2, S2.GetPerimeter([.. loop]).Radians);
    }

    [Fact]
    internal void Test_GetPerimeter_MoreThanTwoPi()
    {
        // Make sure that GetPerimeter doesn't use S1ChordAngle, which can only
        // represent distances up to 2*Pi.
        var loop = ParsePointsOrDie("0:0, 0:90, 0:180, 90:0, 0:-90").ToArray();
        Assert2.DoubleEqual(5 * S2.M_PI_2, S2.GetPerimeter(loop).Radians);
    }

    [Fact]
    internal void Test_GetSignedArea_Underflow()
    {
        var loop = ParsePointsOrDie("0:0, 0:1e-88, 1e-88:1e-88, 1e-88:0");
        Assert.True(S2.GetSignedArea(loop) > 0);
    }

    [Fact]
    internal void Test_LoopTestBase_GetAreaConsistentWithCurvature()
    {
        TestAreaConsistentWithCurvature(full_);
        TestAreaConsistentWithCurvature(north_hemi_);
        TestAreaConsistentWithCurvature(north_hemi3_);
        TestAreaConsistentWithCurvature(west_hemi_);
        TestAreaConsistentWithCurvature(east_hemi_);
        TestAreaConsistentWithCurvature(candy_cane_);
        TestAreaConsistentWithCurvature(line_triangle_);
        TestAreaConsistentWithCurvature(skinny_chevron_);
        TestAreaConsistentWithCurvature(three_leaf_clover_);
        TestAreaConsistentWithCurvature(tessellated_loop_);
    }

    [Fact]
    internal void Test_LoopTestBase_GetSurfaceIntegralGreaterThan4Pi()
    {
        // This test demonstrates that even when GetSurfaceIntegral() returns a an
        // area greater than 4*Pi, GetSignedArea() still returns an accurate result.

        // GetSurfaceIntegral() returns an area > 4 * Pi for this loop.  (Note that
        // the result of GetSurfaceIntegral is only correct modulo 4 * Pi, and that
        // S2::GetSignedArea() automatically corrects for this.)
        S2PointLoopSpan loop1 = [
            new(1, 0, 0), new(0, 1, 1e-150), new S2Point(-1, -2, 0).Normalize(),
            new(-1, 0, 1e-50), new(0, 0, 1)
        ];
        Assert.True(new S2Loop(loop1).IsValid());
        Assert.True(S2.GetSurfaceIntegral(loop1, S2.SignedArea) > 4 * S2.M_PI + 0.1);
        TestAreaConsistentWithCurvature(loop1);
    }

    [Fact]
    internal void Test_LoopTestBase_GetAreaConsistentWithOrientation()
    {
        // Test that GetArea() returns an area near 0 for degenerate loops that
        // contain almost no points, and an area near 4*Pi for degenerate loops that
        // contain almost all points.

        const int kMaxVertices = 6;
        for (int i = 0; i < 50; ++i)
        {
            int num_vertices = 3 + S2Testing.Random.Uniform(kMaxVertices - 3 + 1);
            // Repeatedly choose N vertices that are exactly on the equator until we
            // find some that form a valid loop.
            S2PointLoopSpan loop = [];
            do
            {
                for (int i2 = 0; i2 < num_vertices; ++i2)
                {
                    // We limit longitude to the range [0, 90] to ensure that the loop is
                    // degenerate (as opposed to following the entire equator).
                    loop.Add(
                        S2LatLng.FromRadians(0, S2Testing.Random.RandDouble() * S2.M_PI_2).ToPoint());
                }
            } while (!new S2Loop(loop, S2Debug.DISABLE).IsValid());
            bool ccw = S2.IsNormalized(loop);
            // The error bound is sufficient for current tests but not guaranteed.
            _ = i + ": " + loop.ToDebugString();
            Assert2.Near(ccw ? 0 : S2.M_4_PI, S2.GetArea(loop), 1e-14);
            Assert.Equal(!ccw, new S2Loop(loop).Contains(new S2Point(0, 0, 1)));
        }
    }

    [Fact]
    internal void Test_LoopTestBase_GetAreaAccuracy()
    {
        // TODO(ericv): Test that GetArea() has an accuracy significantly better
        // than 1e-15 on loops whose area is small.
    }

    [Fact]
    internal void Test_LoopTestBase_GetAreaAndCentroid()
    {
        Assert.Equal(S2.M_4_PI, S2.GetArea(full_));
        Assert.Equal(S2Point.Empty, S2.GetCentroid(full_));

        Assert2.DoubleEqual(S2.GetArea(north_hemi_), S2.M_2_PI);
        Assert2.Near(S2.M_2_PI, S2.GetArea(east_hemi_), 1e-12);

        // Construct spherical caps of random height, and approximate their boundary
        // with closely spaces vertices.  Then check that the area and centroid are
        // correct.
        for (int iter = 0; iter < 50; ++iter)
        {
            // Choose a coordinate frame for the spherical cap.
            S2Testing.GetRandomFrame(out var x, out var y, out var z);

            // Given two points at latitude phi and whose longitudes differ by dtheta,
            // the geodesic between the two points has a maximum latitude of
            // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
            // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
            //
            // We want to position the vertices close enough together so that their
            // maximum distance from the boundary of the spherical cap is kMaxDist.
            // Thus we want Math.Abs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
            const double kMaxDist = 1e-6;
            double height = 2 * S2Testing.Random.RandDouble();
            double phi = Math.Asin(1 - height);
            double max_dtheta = 2 * Math.Acos(Math.Tan(Math.Abs(phi)) / Math.Tan(Math.Abs(phi) + kMaxDist));
            max_dtheta = Math.Min(Math.PI, max_dtheta);  // At least 3 vertices.

            S2PointLoopSpan loop = [];
            for (double theta = 0; theta < S2.M_2_PI;
                 theta += S2Testing.Random.RandDouble() * max_dtheta)
            {
                loop.Add(Math.Cos(theta) * Math.Cos(phi) * x +
                               Math.Sin(theta) * Math.Cos(phi) * y +
                               Math.Sin(phi) * z);
            }
            double area = S2.GetArea(loop);
            S2Point centroid = S2.GetCentroid(loop);
            double expected_area = S2.M_2_PI * height;
            Assert.True(Math.Abs(area - expected_area) <= S2.M_2_PI * kMaxDist);
            S2Point expected_centroid = expected_area * (1 - 0.5 * height) * z;
            Assert.True((centroid - expected_centroid).Norm() <= 2 * kMaxDist);
        }
    }

    [Fact]
    internal void Test_LoopTestBase_GetCurvature()
    {
        Assert.Equal(-S2.M_2_PI, S2.GetCurvature(full_));

        Assert.Equal(S2.M_2_PI, S2.GetCurvature(v_loop_));
        CheckCurvatureInvariants(v_loop_);

        // This curvature should be computed exactly.
        Assert.Equal(0, S2.GetCurvature(north_hemi3_));
        CheckCurvatureInvariants(north_hemi3_);

        Assert2.Near(0, S2.GetCurvature(west_hemi_), 1e-15);
        CheckCurvatureInvariants(west_hemi_);

        // We don't have an easy way to estimate the curvature of these loops, but
        // we can still check that the expected invariants hold.
        CheckCurvatureInvariants(candy_cane_);
        CheckCurvatureInvariants(three_leaf_clover_);

        Assert2.DoubleEqual(S2.M_2_PI, S2.GetCurvature(line_triangle_));
        CheckCurvatureInvariants(line_triangle_);

        Assert2.DoubleEqual(S2.M_2_PI, S2.GetCurvature(skinny_chevron_));
        CheckCurvatureInvariants(skinny_chevron_);

        // Build a narrow spiral loop starting at the north pole.  This is designed
        // to test that the error in GetCurvature is linear in the number of
        // vertices even when the partial sum of the curvatures gets very large.
        // The spiral consists of two "arms" defining opposite sides of the loop.
        // This is a pathological loop that contains many long parallel edges.
        int kArmPoints = 10000;    // Number of vertices in each "arm"
        double kArmRadius = 0.01;  // Radius of spiral.
        S2PointLoopSpan spiral = new(2 * kArmPoints)
        {
            [kArmPoints] = new S2Point(0, 0, 1)
        };
        for (int i = 0; i < kArmPoints; ++i)
        {
            double angle = S2.M_2_PI / 3 * i;
            double x = Math.Cos(angle);
            double y = Math.Sin(angle);
            double r1 = i * kArmRadius / kArmPoints;
            double r2 = (i + 1.5) * kArmRadius / kArmPoints;
            spiral[kArmPoints - i - 1] = new S2Point(r1 * x, r1 * y, 1).Normalize();
            spiral[kArmPoints + i] = new S2Point(r2 * x, r2 * y, 1).Normalize();
        }

        // Check that GetCurvature() is consistent with GetArea() to within the
        // error bound of the former.  We actually use a tiny fraction of the
        // worst-case error bound, since the worst case only happens when all the
        // roundoff errors happen in the same direction and this test is not
        // designed to achieve that.  The error in GetArea() can be ignored for the
        // purposes of this test since it is generally much smaller.
        Assert2.Near(
            S2.M_2_PI - S2.GetArea(spiral),
            S2.GetCurvature(spiral),
            0.01 * S2.GetCurvatureMaxError(spiral));
    }

    // Given a string where each character "ch" represents a vertex (such as
    // "abac"), returns a vector of S2Points of the form (ch, 0, 0).  Note that
    // these points are not unit length and therefore are not suitable for general
    // use, however they are useful for testing certain functions below.
    private static S2PointLoopSpan MakeTestLoop(string loop_str) =>
        S2Testing.StrPoints.StrToPoints(loop_str);

    // Given a loop whose vertices are represented as characters (such as "abcd" or
    // "abccb"), verify that S2.PruneDegeneracies() yields the loop "expected".
    private static void TestPruneDegeneracies(string input_str, string expected_str)
    {
        var input = MakeTestLoop(input_str);
        var actual_str = S2Testing.StrPoints.PointsToStr(S2.PruneDegeneracies(input));
        Assert.Equal(expected_str, actual_str);
    }

    // Given a loop whose vertices are represented as characters (such as "abcd" or
    // "abccb"), verify that S2.GetCanonicalLoopOrder returns the given result.
    private static void TestCanonicalLoopOrder(string input_str, S2.LoopOrder expected_order)
    {
        Assert.Equal(expected_order, S2.GetCanonicalLoopOrder(MakeTestLoop(input_str)));
    }

    private static void TestAreaConsistentWithCurvature(S2PointLoopSpan loop)
    {
        // Check that the area computed using GetArea() is consistent with the loop
        // curvature.  According to the Gauss-Bonnet theorem, the area of the loop
        // equals 2*Pi minus its curvature.
        double area = S2.GetArea(loop);
        double gauss_area = S2.M_2_PI - S2.GetCurvature(loop);
        // The error bound below is sufficient for current tests but not guaranteed.
        _ = loop.ToDebugString();
        Assert.True(Math.Abs(area - gauss_area) <= 1e-14);
    }

    private static void ExpectSameOrder(S2PointLoopSpan loop1, S2.LoopOrder order1, S2PointLoopSpan loop2, S2.LoopOrder order2)
    {
        Assert.Equal(loop1.Count, loop2.Count);
        int i1 = order1.First, i2 = order2.First;
        int dir1 = order1.Dir, dir2 = order2.Dir;
        for (int n = loop1.Count; --n >= 0;)
        {
            Assert.Equal(loop1[i1], loop2[i2]);
            i1 += dir1;
            i2 += dir2;
        }
    }

    // Check that the curvature is *identical* when the vertex order is
    // rotated, and that the sign is inverted when the vertices are reversed.
    private static void CheckCurvatureInvariants(S2PointLoopSpan loop_in)
    {
        S2.LoopOrder order_in = S2.GetCanonicalLoopOrder(loop_in);
        var loop = loop_in.ToList();
        double expected = S2.GetCurvature(loop);
        for (int i = 0; i < loop.Count; ++i)
        {
            loop.Reverse();
            Assert.Equal((expected == S2.M_2_PI) ? expected : -expected,
                      S2.GetCurvature(loop));
            ExpectSameOrder(loop_in, order_in, loop, S2.GetCanonicalLoopOrder(loop));
            loop.Reverse();
            loop.RotateInPlace(1);
            Assert.Equal(expected, S2.GetCurvature(loop));
            ExpectSameOrder(loop_in, order_in, loop, S2.GetCanonicalLoopOrder(loop));
        }
    }
}
