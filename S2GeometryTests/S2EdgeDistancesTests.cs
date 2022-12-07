namespace S2Geometry;

public class S2EdgeDistancesTests
{
    [Fact]
    internal void Test_S2_GetUpdateMinDistanceMaxError()
    {
        // Verify that the error is "reasonable" for a sampling of distances.
        CheckUpdateMinDistanceMaxError(0, 1.5e-15);
        CheckUpdateMinDistanceMaxError(1e-8, 1e-15);
        CheckUpdateMinDistanceMaxError(1e-5, 1e-15);
        CheckUpdateMinDistanceMaxError(0.05, 1e-15);
        CheckUpdateMinDistanceMaxError(S2.M_PI_2 - 1e-8, 2e-15);
        CheckUpdateMinDistanceMaxError(S2.M_PI_2, 2e-15);
        CheckUpdateMinDistanceMaxError(S2.M_PI_2 + 1e-8, 2e-15);
        CheckUpdateMinDistanceMaxError(Math.PI - 1e-5, 2e-10);
        CheckUpdateMinDistanceMaxError(Math.PI, 0);
    }

    [Fact]
    internal void Test_S2_GetUpdateMinInteriorDistanceMaxError()
    {
        // Check that the error bound returned by
        // GetUpdateMinInteriorDistanceMaxError() is large enough.
        for (int iter = 0; iter < 10000; ++iter)
        {
            S2Point a0 = S2Testing.RandomPoint();
            var lenRadians = Math.PI * Math.Pow(1e-20, S2Testing.Random.RandDouble());
            S1Angle len = S1Angle.FromRadians(lenRadians);
            if (S2Testing.Random.OneIn(4)) len = S1Angle.FromRadians(S2.M_PI) - len;
            S2Point a1 = S2.GetPointOnLine(a0, S2Testing.RandomPoint(), len);

            // TODO(ericv): The error bound holds for antipodal points, but the S2
            // predicates used to test the error do not support antipodal points yet.
            if (a1 == -a0) continue;
            S2Point n = S2.RobustCrossProd(a0, a1).Normalize();
            double f = Math.Pow(1e-20, S2Testing.Random.RandDouble());
            S2Point a = ((1 - f) * a0 + f * a1).Normalize();
            var rRadians = S2.M_PI_2 * Math.Pow(1e-20, S2Testing.Random.RandDouble());
            S1Angle r = S1Angle.FromRadians(rRadians);
            if (S2Testing.Random.OneIn(2)) r = S1Angle.FromRadians(S2.M_PI_2) - r;
            S2Point x = S2.GetPointOnLine(a, n, r);
            S1ChordAngle min_dist = S1ChordAngle.Infinity;
            if (!S2.UpdateMinInteriorDistance(x, a0, a1, ref min_dist))
            {
                --iter; continue;
            }
            double error = S2.GetUpdateMinDistanceMaxError(min_dist);
            Assert.True(S2Pred.CompareEdgeDistance(x, a0, a1, min_dist.PlusError(error)) <= 0);
            Assert.True(S2Pred.CompareEdgeDistance(x, a0, a1, min_dist.PlusError(-error)) >= 0);
        }
    }

    [Fact]
    internal void Test_S2_Distance()
    {
        CheckDistance(new(1, 0, 0), new(1, 0, 0), new(0, 1, 0), 0, new(1, 0, 0));
        CheckDistance(new(0, 1, 0), new(1, 0, 0), new(0, 1, 0), 0, new(0, 1, 0));
        CheckDistance(new(1, 3, 0), new(1, 0, 0), new(0, 1, 0), 0, new(1, 3, 0));
        CheckDistance(new(0, 0, 1), new(1, 0, 0), new(0, 1, 0), S2.M_PI_2, new(1, 0, 0));
        CheckDistance(new(0, 0, -1), new(1, 0, 0), new(0, 1, 0), S2.M_PI_2, new(1, 0, 0));
        CheckDistance(new(-1, -1, 0), new(1, 0, 0), new(0, 1, 0), 0.75 * Math.PI, new());

        CheckDistance(new(0, 1, 0), new(1, 0, 0), new(1, 1, 0), S2.M_PI_4, new(1, 1, 0));
        CheckDistance(new(0, -1, 0), new(1, 0, 0), new(1, 1, 0), S2.M_PI_2, new(1, 0, 0));

        CheckDistance(new(0, -1, 0), new(1, 0, 0), new(-1, 1, 0), S2.M_PI_2, new(1, 0, 0));
        CheckDistance(new(-1, -1, 0), new(1, 0, 0), new(-1, 1, 0), S2.M_PI_2, new(-1, 1, 0));

        CheckDistance(new(1, 1, 1), new(1, 0, 0), new(0, 1, 0), Math.Asin(Math.Sqrt(1.0 / 3)), new(1, 1, 0));
        CheckDistance(new(1, 1, -1), new(1, 0, 0), new(0, 1, 0), Math.Asin(Math.Sqrt(1.0 / 3)), new(1, 1, 0));

        CheckDistance(new(-1, 0, 0), new(1, 1, 0), new(1, 1, 0), 0.75 * Math.PI, new(1, 1, 0));
        CheckDistance(new(0, 0, -1), new(1, 1, 0), new(1, 1, 0), S2.M_PI_2, new(1, 1, 0));
        CheckDistance(new(-1, 0, 0), new(1, 0, 0), new(1, 0, 0), Math.PI, new(1, 0, 0));
    }

    [Fact]
    internal void Test_S2_UpdateMinInteriorDistanceLowerBoundOptimizationIsConservative()
    {
        // Verifies that AlwaysUpdateMinInteriorDistance() computes the lower bound
        // on the true distance conservatively.  (This test used to fail.)
        S2Point x = new(-0.017952729194524016, -0.30232422079175203, 0.95303607751077712);
        S2Point a = new(-0.017894725505830295, -0.30229974986194175, 0.95304493075220664);
        S2Point b = new(-0.017986591360900289, -0.30233851195954353, 0.95303090543659963);
        S1ChordAngle min_distance = S1ChordAngle.Infinity;
        Assert.True(S2.UpdateMinDistance(x, a, b, ref min_distance));
        min_distance = min_distance.Successor();
        Assert.True(S2.UpdateMinDistance(x, a, b, ref min_distance));
    }

    [Fact]
    internal void Test_S2_UpdateMinInteriorDistanceRejectionTestIsConservative()
    {
        // This test checks several representative cases where previously
        // UpdateMinInteriorDistance was failing to update the distance because a
        // rejection test was not being done conservatively.
        //
        // Note that all of the edges AB in this test are nearly antipodal.
        {
            S2Point x = new(1, -4.6547732744037044e-11, -5.6374428459823598e-89);
            S2Point a = new(1, -8.9031850507928352e-11, 0);
            S2Point b = new(-0.99999999999996347, 2.7030110029169596e-07,
                1.555092348806121e-99);
            var min_dist = S1ChordAngle.FromLength2(6.3897233584120815e-26);
            Assert.True(S2.UpdateMinInteriorDistance(x, a, b, ref min_dist));
        }
        {
            S2Point x = new(1, -4.7617930898495072e-13, 0);
            S2Point a = new(-1, -1.6065916409055676e-10, 0);
            S2Point b = new(1, 0, 9.9964883247706732e-35);
            var min_dist = S1ChordAngle.FromLength2(6.3897233584120815e-26);
            Assert.True(S2.UpdateMinInteriorDistance(x, a, b, ref min_dist));
        }
        {
            S2Point x = new(1, 0, 0);
            S2Point a = new(1, -8.4965026896454536e-11, 0);
            S2Point b = new(-0.99999999999966138, 8.2297529603339328e-07,
                9.6070344113320997e-21);
            var min_dist = S1ChordAngle.FromLength2(6.3897233584120815e-26);
            Assert.True(S2.UpdateMinInteriorDistance(x, a, b, ref min_dist));
        }
    }

    [Fact]
    internal void Test_S2_MaxDistance()
    {
        CheckMaxDistance(new(1, 0, 1), new(1, 0, 0), new(0, 1, 0), S2.M_PI_2);
        CheckMaxDistance(new(1, 0, -1), new(1, 0, 0), new(0, 1, 0), S2.M_PI_2);
        CheckMaxDistance(new(0, 1, 1), new(1, 0, 0), new(0, 1, 0), S2.M_PI_2);
        CheckMaxDistance(new(0, 1, -1), new(1, 0, 0), new(0, 1, 0), S2.M_PI_2);

        CheckMaxDistance(new(1, 1, 1), new(1, 0, 0), new(0, 1, 0), Math.Asin(Math.Sqrt(2.0 / 3)));
        CheckMaxDistance(new(1, 1, -1), new(1, 0, 0), new(0, 1, 0), Math.Asin(Math.Sqrt(2.0 / 3)));

        CheckMaxDistance(new(1, 0, 0), new(1, 1, 0), new(1, -1, 0), S2.M_PI_4);
        CheckMaxDistance(new(0, 1, 0), new(1, 1, 0), new(-1, 1, 0), S2.M_PI_4);
        CheckMaxDistance(new(0, 0, 1), new(0, 1, 1), new(0, -1, 1), S2.M_PI_4);

        CheckMaxDistance(new(0, 0, 1), new(1, 0, 0), new(1, 0, -1), 3 * S2.M_PI_4);
        CheckMaxDistance(new(0, 0, 1), new(1, 0, 0), new(1, 1, -S2.M_SQRT2), 3 * S2.M_PI_4);

        CheckMaxDistance(new(0, 0, 1), new(0, 0, -1), new(0, 0, -1), Math.PI);
    }

    private static void TestInterpolate(S2Point a, S2Point b, double t, S2Point expected)
    {
        a = a.Normalize();
        b = b.Normalize();
        expected = expected.Normalize();

        // We allow a bit more than the usual 1e-15 error tolerance because
        // interpolation uses trig functions.
        S1Angle kError = S1Angle.FromRadians(3e-15);
        Assert.True(new S1Angle(S2.Interpolate(a, b, t), expected) <= kError);

        // Now test the other interpolation functions.
        S1Angle r = t * new S1Angle(a, b);
        Assert.True(new S1Angle(S2.GetPointOnLine(a, b, r), expected) <= kError);
        if (a.DotProd(b) == 0)
        {  // Common in the test cases below.
            Assert.True(new S1Angle(S2.GetPointOnRay(a, b, r), expected) <= kError);
        }
        if (r.Radians >= 0 && r.Radians < 0.99 * S2.M_PI)
        {
            S1ChordAngle r_ca = new(r);
            Assert.True(new S1Angle(S2.GetPointOnLine(a, b, r_ca), expected) <= kError);
            if (a.DotProd(b) == 0)
            {
                Assert.True(new S1Angle(S2.GetPointOnRay(a, b, r_ca), expected) <= kError);
            }
        }
    }

    [Fact]
    internal void Test_S2_Interpolate()
    {
        // Choose test points designed to expose floating-point errors.
        S2Point p1 = new S2Point(0.1, 1e-30, 0.3).Normalize();
        S2Point p2 = new S2Point(-0.7, -0.55, -1e30).Normalize();

        // A zero-length edge.
        TestInterpolate(p1, p1, 0, p1);
        TestInterpolate(p1, p1, 1, p1);

        // Start, end, and middle of a medium-length edge.
        TestInterpolate(p1, p2, 0, p1);
        TestInterpolate(p1, p2, 1, p2);
        TestInterpolate(p1, p2, 0.5, 0.5 * (p1 + p2));

        // Test that interpolation is done using distances on the sphere rather than
        // linear distances.
        TestInterpolate(new S2Point(1, 0, 0), new S2Point(0, 1, 0), 1.0 / 3,
                        new S2Point(Math.Sqrt(3), 1, 0));
        TestInterpolate(new S2Point(1, 0, 0), new S2Point(0, 1, 0), 2.0 / 3,
                        new S2Point(1, Math.Sqrt(3), 0));

        // Test that interpolation is accurate on a long edge (but not so long that
        // the definition of the edge itself becomes too unstable).
        {
            double kLng = Math.PI - 1e-2;
            S2Point a = S2LatLng.FromRadians(0, 0).ToPoint();
            S2Point b = S2LatLng.FromRadians(0, kLng).ToPoint();
            for (double f = 0.4; f > 1e-15; f *= 0.1)
            {
                TestInterpolate(a, b, f, S2LatLng.FromRadians(0, f * kLng).ToPoint());
                TestInterpolate(a, b, 1 - f, S2LatLng.FromRadians(0, (1 - f) * kLng).ToPoint());
            }
        }

        // Test that interpolation on a 180 degree edge (antipodal endpoints) yields
        // a result with the correct distance from each endpoint.
        for (double t = 0; t <= 1; t += 0.125)
        {
            S2Point actual = S2.Interpolate(p1, -p1, t);
            Assert2.Near(new S1Angle(actual, p1).Radians, t * Math.PI, 3e-15);
        }
    }

    [Fact]
    internal void Test_S2_InterpolateCanExtrapolate()
    {
        S2Point i = new(1, 0, 0);
        S2Point j = new(0, 1, 0);
        // Initial vectors at 90 degrees.
        TestInterpolate(i, j, 0, new S2Point(1, 0, 0));
        TestInterpolate(i, j, 1, new S2Point(0, 1, 0));
        TestInterpolate(i, j, 1.5, new S2Point(-1, 1, 0));
        TestInterpolate(i, j, 2, new S2Point(-1, 0, 0));
        TestInterpolate(i, j, 3, new S2Point(0, -1, 0));
        TestInterpolate(i, j, 4, new S2Point(1, 0, 0));

        // Negative values of t.
        TestInterpolate(i, j, -1, new S2Point(0, -1, 0));
        TestInterpolate(i, j, -2, new S2Point(-1, 0, 0));
        TestInterpolate(i, j, -3, new S2Point(0, 1, 0));
        TestInterpolate(i, j, -4, new S2Point(1, 0, 0));

        // Initial vectors at 45 degrees.
        TestInterpolate(i, new S2Point(1, 1, 0), 2, new S2Point(0, 1, 0));
        TestInterpolate(i, new S2Point(1, 1, 0), 3, new S2Point(-1, 1, 0));
        TestInterpolate(i, new S2Point(1, 1, 0), 4, new S2Point(-1, 0, 0));

        // Initial vectors at 135 degrees.
        TestInterpolate(i, new S2Point(-1, 1, 0), 2, new S2Point(0, -1, 0));

        // Take a small fraction along the curve.
        S2Point p = S2.Interpolate(i, j, 0.001);
        // We should get back where we started.
        TestInterpolate(i, p, 1000, j);
    }

    [Fact]
    internal void Test_S2_RepeatedInterpolation()
    {
        // Check that points do not drift away from unit length when repeated
        // interpolations are done.
        for (int i = 0; i < 100; ++i)
        {
            S2Point a = S2Testing.RandomPoint();
            S2Point b = S2Testing.RandomPoint();
            for (int j = 0; j < 1000; ++j)
            {
                a = S2.Interpolate(a, b, 0.01);
            }
            Assert.True(a.IsUnitLength());
        }
    }

    [Fact]
    internal void Test_S2_EdgePairMinDistance()
    {
        // One edge is degenerate.
        CheckEdgePairMinDistance(new S2Point(1, 0, 1), new S2Point(1, 0, 1),
                                 new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 S2.M_PI_4, new S2Point(1, 0, 1), new S2Point(1, 0, 0));
        CheckEdgePairMinDistance(new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 new S2Point(1, 0, 1), new S2Point(1, 0, 1),
                                 S2.M_PI_4, new S2Point(1, 0, 0), new S2Point(1, 0, 1));

        // Both edges are degenerate.
        CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 new S2Point(0, 1, 0), new S2Point(0, 1, 0),
                                 S2.M_PI_2, new S2Point(1, 0, 0), new S2Point(0, 1, 0));

        // Both edges are degenerate and antipodal.
        CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 new S2Point(-1, 0, 0), new S2Point(-1, 0, 0),
                                 Math.PI, new S2Point(1, 0, 0), new S2Point(-1, 0, 0));

        // Two identical edges.
        CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 0, S2Point.Empty, S2Point.Empty);

        // Both edges are degenerate and identical.
        CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 0, new S2Point(1, 0, 0), new S2Point(1, 0, 0));

        // Edges that share exactly one vertex (all 4 possibilities).
        CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 new S2Point(0, 1, 0), new S2Point(0, 1, 1),
                                 0, new S2Point(0, 1, 0), new S2Point(0, 1, 0));
        CheckEdgePairMinDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0),
                                 new S2Point(0, 1, 0), new S2Point(0, 1, 1),
                                 0, new S2Point(0, 1, 0), new S2Point(0, 1, 0));
        CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 new S2Point(0, 1, 1), new S2Point(0, 1, 0),
                                 0, new S2Point(0, 1, 0), new S2Point(0, 1, 0));
        CheckEdgePairMinDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0),
                                 new S2Point(0, 1, 1), new S2Point(0, 1, 0),
                                 0, new S2Point(0, 1, 0), new S2Point(0, 1, 0));

        // Two edges whose interiors cross.
        CheckEdgePairMinDistance(new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 new S2Point(1, 0, -1), new S2Point(1, 0, 1),
                                 0, new S2Point(1, 0, 0), new S2Point(1, 0, 0));

        // The closest distance occurs between two edge endpoints, but more than one
        // endpoint pair is equally distant.
        CheckEdgePairMinDistance(new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 new S2Point(-1, 0, 0), new S2Point(-1, 0, 1),
                                 Math.Acos(-0.5), S2Point.Empty, new S2Point(-1, 0, 1));
        CheckEdgePairMinDistance(new S2Point(-1, 0, 0), new S2Point(-1, 0, 1),
                                 new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 Math.Acos(-0.5), new S2Point(-1, 0, 1), S2Point.Empty);
        CheckEdgePairMinDistance(new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 new S2Point(-1, 0, -1), new S2Point(-1, 0, 1),
                                 Math.Acos(-0.5), S2Point.Empty, S2Point.Empty);
    }

    [Fact]
    internal void Test_S2_EdgePairMaxDistance()
    {
        // Standard situation.  Same hemisphere, not degenerate.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 new S2Point(1, 1, 0), new S2Point(1, 1, 1),
                                 Math.Acos(1 / Math.Sqrt(3)));

        // One edge is degenerate.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 1), new S2Point(1, 0, 1),
                                 new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 Math.Acos(0.5));
        CheckEdgePairMaxDistance(new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                 new S2Point(1, 0, 1), new S2Point(1, 0, 1),
                                 Math.Acos(0.5));

        // Both edges are degenerate.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 new S2Point(0, 1, 0), new S2Point(0, 1, 0),
                                 S2.M_PI_2);

        // Both edges are degenerate and antipodal.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 new S2Point(-1, 0, 0), new S2Point(-1, 0, 0),
                                 Math.PI);

        // Two identical edges.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                 S2.M_PI_2);

        // Both edges are degenerate and identical.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                 0);

        // Antipodal reflection of one edge crosses the other edge.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 1), new S2Point(1, 0, -1),
                                 new S2Point(-1, -1, 0), new S2Point(-1, 1, 0),
                                 Math.PI);

        // One vertex of one edge touches the interior of the antipodal reflection
        // of the other edge.
        CheckEdgePairMaxDistance(new S2Point(1, 0, 1), new S2Point(1, 0, 0),
                                 new S2Point(-1, -1, 0), new S2Point(-1, 1, 0),
                                 Math.PI);
    }

    [Fact]
    internal void Test_S2_EdgeBNearEdgeA()
    {
        // Edge is near itself.
        Assert.True(IsEdgeBNearEdgeA("5:5, 10:-5", "5:5, 10:-5", 1e-6));

        // Edge is near its reverse
        Assert.True(IsEdgeBNearEdgeA("5:5, 10:-5", "10:-5, 5:5", 1e-6));

        // Short edge is near long edge.
        Assert.True(IsEdgeBNearEdgeA("10:0, -10:0", "2:1, -2:1", 1.0));

        // Long edges cannot be near shorter edges.
        Assert.False(IsEdgeBNearEdgeA("2:1, -2:1", "10:0, -10:0", 1.0));

        // Orthogonal crossing edges are not near each other...
        Assert.False(IsEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 1.0));

        // ... unless all points on B are within tolerance of A.
        Assert.True(IsEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 2.0));

        // Very long edges whose endpoints are close may have interior points that are
        // far apart.  An implementation that only considers the vertices of polylines
        // will incorrectly consider such edges as "close" when they are not.
        // Consider, for example, two consecutive lines of longitude.  As they
        // approach the poles, they become arbitrarily close together, but along the
        // equator they bow apart.
        Assert.False(IsEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 0.5));
        Assert.True(IsEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 1.5));

        // Make sure that the result is independent of the edge directions.
        Assert.True(IsEdgeBNearEdgeA("89:1, -89:1", "-89:2, 89:2", 1.5));

        // Cases where the point that achieves the maximum distance to A is the
        // interior point of B that is equidistant from the endpoints of A.  This
        // requires two long edges A and B whose endpoints are near each other but
        // where B intersects the perpendicular bisector of the endpoints of A in
        // the hemisphere opposite A's midpoint.  Furthermore these cases are
        // constructed so that the points where circ(A) is furthest from circ(B) do
        // not project onto the interior of B.
        Assert.False(IsEdgeBNearEdgeA("0:-100, 0:100", "5:-80, -5:80", 70.0));
        Assert.False(IsEdgeBNearEdgeA("0:-100, 0:100", "1:-35, 10:35", 70.0));

        // Make sure that the result is independent of the edge directions.
        Assert.False(IsEdgeBNearEdgeA("0:-100, 0:100", "5:80, -5:-80", 70.0));

        // The two arcs here are nearly as long as S2 edges can be (just shy of 180
        // degrees), and their endpoints are less than 1 degree apart.  Their
        // midpoints, however, are at opposite ends of the sphere along its equator.
        Assert.False(IsEdgeBNearEdgeA(
                         "0:-179.75, 0:-0.25", "0:179.75, 0:0.25", 1.0));

        // At the equator, the second arc here is 9.75 degrees from the first, and
        // closer at all other points.  However, the southern point of the second arc
        // (-1, 9.75) is too far from the first arc for the short-circuiting logic in
        // IsEdgeBNearEdgeA to apply.
        Assert.True(IsEdgeBNearEdgeA("40:0, -5:0", "39:0.975, -1:0.975", 1.0));

        // Same as above, but B's orientation is reversed, causing the angle between
        // the normal vectors of circ(B) and circ(A) to be (180-9.75) = 170.5 degrees,
        // not 9.75 degrees.  The greatest separation between the planes is still 9.75
        // degrees.
        Assert.True(IsEdgeBNearEdgeA("10:0, -10:0", "-.4:0.975, 0.4:0.975", 1.0));

        // A and B are on the same great circle, A and B partially overlap, but the
        // only part of B that does not overlap A is shorter than tolerance.
        Assert.True(IsEdgeBNearEdgeA("0:0, 1:0", "0.9:0, 1.1:0", 0.25));

        // A and B are on the same great circle, all points on B are close to A at its
        // second endpoint, (1,0).
        Assert.True(IsEdgeBNearEdgeA("0:0, 1:0", "1.1:0, 1.2:0", 0.25));

        // Same as above, but B's orientation is reversed.  This case is special
        // because the projection of the normal defining A onto the plane containing B
        // is the null vector, and must be handled by a special case.
        Assert.True(IsEdgeBNearEdgeA("0:0, 1:0", "1.2:0, 1.1:0", 0.25));
    }

    // Given a point X and an edge AB, check that the distance from X to AB is
    // "distance_radians" and the closest point on AB is "expected_closest".
    private static void CheckDistance(S2Point x, S2Point a, S2Point b, double distance_radians, S2Point expected_closest)
    {
        x = x.Normalize();
        a = a.Normalize();
        b = b.Normalize();
        expected_closest = expected_closest.Normalize();
        Assert2.Near(distance_radians, S2.GetDistance(x, a, b).Radians, S2.DoubleError);
        S2Point closest = S2.Project(x, a, b);
        Assert.True(S2Pred.CompareEdgeDistance(
            closest, a, b, new S1ChordAngle(S2.kProjectPerpendicularErrorS1Angle)) < 0);

        // If X is perpendicular to AB then there is nothing further we can expect.
        if (distance_radians != S2.M_PI_2)
        {
            if (expected_closest == new S2Point())
            {
                // This special value says that the result should be A or B.
                Assert.True(closest == a || closest == b);
            }
            else
            {
                Assert.True(S2.ApproxEquals(closest, expected_closest));
            }
        }
        S1ChordAngle min_distance = S1ChordAngle.Zero;
        Assert.False(S2.UpdateMinDistance(x, a, b, ref min_distance));
        min_distance = S1ChordAngle.Infinity;
        Assert.True(S2.UpdateMinDistance(x, a, b, ref min_distance));
        Assert2.Near(distance_radians, min_distance.ToAngle().Radians, S2.DoubleError);
    }

    private static void CheckMaxDistance(S2Point x, S2Point a, S2Point b, double distance_radians)
    {
        x = x.Normalize();
        a = a.Normalize();
        b = b.Normalize();

        S1ChordAngle max_distance = S1ChordAngle.Straight;
        Assert.False(S2.UpdateMaxDistance(x, a, b, ref max_distance));
        max_distance = S1ChordAngle.Negative;
        Assert.True(S2.UpdateMaxDistance(x, a, b, ref max_distance));
        Assert2.Near(distance_radians, max_distance.Radians(), S2.DoubleError);
    }

    // Chooses a random S2Point that is often near the intersection of one of the
    // coodinates planes or coordinate axes with the unit sphere.  (It is possible
    // to represent very small perturbations near such points.)
    private static S2Point ChoosePoint()
    {
        var x = S2Testing.RandomPoint().ToArray();
        for (int i = 0; i < 3; ++i)
        {
            if (S2Testing.Random.OneIn(3))
            {
                x[i] *= Math.Pow(1e-50, S2Testing.Random.RandDouble());
            }
        }
        return new S2Point(x).Normalize();
    }

    [Fact]
    internal void Test_S2_ProjectError()
    {
        for (int iter = 0; iter < 1000; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point a = ChoosePoint();
            S2Point b = ChoosePoint();
            S2Point n = S2.RobustCrossProd(a, b).Normalize();
            S2Point x = S2Testing.SamplePoint(new S2Cap(n, S1Angle.FromRadians(1e-15)));
            S2Point p = S2.Project(x, a, b);
            Assert.True(S2Pred.CompareEdgeDistance(
                p, a, b, new S1ChordAngle(S2.kProjectPerpendicularErrorS1Angle)) < 0);
        }
    }

    // Given two edges a0a1 and b0b1, check that the minimum distance between them
    // is "distance_radians", and that GetEdgePairClosestPoints() returns
    // "expected_a" and "expected_b" as the points that achieve this distance.
    // S2Point.Empty may be passed for "expected_a" or "expected_b" to indicate
    // that both endpoints of the corresponding edge are equally distant, and
    // therefore either one might be returned.
    //
    // Parameters are passed by value so that this function can normalize them.
    private static void CheckEdgePairMinDistance(S2Point a0, S2Point a1, S2Point b0, S2Point b1, double distance_radians, S2Point expected_a, S2Point expected_b)
    {
        a0 = a0.Normalize();
        a1 = a1.Normalize();
        b0 = b0.Normalize();
        b1 = b1.Normalize();
        expected_a = expected_a.Normalize();
        expected_b = expected_b.Normalize();
        var closest = S2.GetEdgePairClosestPoints(a0, a1, b0, b1);
        S2Point actual_a = closest.Item1;
        S2Point actual_b = closest.Item2;
        if (expected_a == S2Point.Empty)
        {
            // This special value says that the result should be a0 or a1.
            Assert.True(actual_a == a0 || actual_a == a1);
        }
        else
        {
            Assert.True(S2.ApproxEquals(expected_a, actual_a));
        }
        if (expected_b == S2Point.Empty)
        {
            // This special value says that the result should be b0 or b1.
            Assert.True(actual_b == b0 || actual_b == b1);
        }
        else
        {
            Assert.True(S2.ApproxEquals(expected_b, actual_b));
        }
        S1ChordAngle min_distance = S1ChordAngle.Zero;
        Assert.False(S2.UpdateEdgePairMinDistance(a0, a1, b0, b1, ref min_distance));
        min_distance = S1ChordAngle.Infinity;
        Assert.True(S2.UpdateEdgePairMinDistance(a0, a1, b0, b1, ref min_distance));
        Assert2.Near(distance_radians, min_distance.Radians(), S2.DoubleError);
    }

    // Given two edges a0a1 and b0b1, check that the maximum distance between them
    // is "distance_radians".  Parameters are passed by value so that this
    // function can normalize them.
    private static void CheckEdgePairMaxDistance(S2Point a0, S2Point a1, S2Point b0, S2Point b1, double distance_radians)
    {
        a0 = a0.Normalize();
        a1 = a1.Normalize();
        b0 = b0.Normalize();
        b1 = b1.Normalize();

        S1ChordAngle max_distance = S1ChordAngle.Straight;
        Assert.False(S2.UpdateEdgePairMaxDistance(a0, a1, b0, b1, ref max_distance));
        max_distance = S1ChordAngle.Negative;
        Assert.True(S2.UpdateEdgePairMaxDistance(a0, a1, b0, b1, ref max_distance));
        Assert2.Near(distance_radians, max_distance.Radians(), S2.DoubleError);
    }

    private static bool IsEdgeBNearEdgeA(string a_str, string b_str, double max_error_degrees)
    {
        var a = MakePolylineOrDie(a_str);
        Assert.Equal(2, a.NumVertices());
        var b = MakePolylineOrDie(b_str);
        Assert.Equal(2, b.NumVertices());
        return S2.IsEdgeBNearEdgeA(
            a.Vertex(0), a.Vertex(1),
            b.Vertex(0), b.Vertex(1),
            S1Angle.FromDegrees(max_error_degrees));
    }

    // Checks that the error returned by S2EdgeDistances.GetUpdateMinDistanceMaxError() for
    // the distance "input" (measured in radians) corresponds to a distance error
    // of less than "max_error" (measured in radians).
    //
    // The reason for the awkward phraseology above is that the value returned by
    // GetUpdateMinDistanceMaxError() is not a distance; it represents an error in
    // the *squared* distance.
    private static void CheckUpdateMinDistanceMaxError(double actual, double max_error)
    {
        S1ChordAngle ca = new(S1Angle.FromRadians(actual));
        S1Angle bound = ca.PlusError(S2.GetUpdateMinDistanceMaxError(ca)).ToAngle();
        Assert.True(bound.Radians - actual <= max_error);
    }

    [Fact]
    internal void Test_S2_GetPointToLeftS1Angle()
    {
        S2Point a = S2LatLng.FromDegrees(0, 0).ToPoint();
        S2Point b = S2LatLng.FromDegrees(0, 5).ToPoint();  // east
        S1Angle kDistance = S2Testing.MetersToAngle(10);

        S2Point c = S2.GetPointToLeft(a, b, kDistance);
        Assert2.Near(new S1Angle(a, c).Radians, kDistance.Radians, 1e-15);
        // CAB must be a right angle with C to the left of AB.
        Assert2.Near(S2.TurnAngle(c, a, b), S2.M_PI_2 /*radians*/, 1e-15);
    }

    [Fact]
    internal void Test_S2_GetPointToLeftS1ChordAngle()
    {
        S2Point a = S2LatLng.FromDegrees(0, 0).ToPoint();
        S2Point b = S2LatLng.FromDegrees(0, 5).ToPoint();  // east
        S1Angle kDistance = S2Testing.MetersToAngle(10);

        S2Point c = S2.GetPointToLeft(a, b, new S1ChordAngle(kDistance));
        Assert2.Near(new S1Angle(a, c).Radians, kDistance.Radians, 1e-15);
        // CAB must be a right angle with C to the left of AB.
        Assert2.Near(S2.TurnAngle(c, a, b), S2.M_PI_2 /*radians*/, 1e-15);
    }

    [Fact]
    internal void Test_S2_GetPointToRightS1Angle()
    {
        S2Point a = S2LatLng.FromDegrees(0, 0).ToPoint();
        S2Point b = S2LatLng.FromDegrees(0, 5).ToPoint();  // east
        S1Angle kDistance = S2Testing.MetersToAngle(10);

        S2Point c = S2.GetPointToRight(a, b, kDistance);
        Assert2.Near(new S1Angle(a, c).Radians, kDistance.Radians, 1e-15);
        // CAB must be a right angle with C to the right of AB.
        Assert2.Near(S2.TurnAngle(c, a, b), -S2.M_PI_2 /*radians*/, 1e-15);
    }

    [Fact]
    internal void Test_S2_GetPointToRightS1ChordAngle()
    {
        S2Point a = S2LatLng.FromDegrees(0, 0).ToPoint();
        S2Point b = S2LatLng.FromDegrees(0, 5).ToPoint();  // east
        S1Angle kDistance = S2Testing.MetersToAngle(10);

        S2Point c = S2.GetPointToRight(a, b, new S1ChordAngle(kDistance));
        Assert2.Near(new S1Angle(a, c).Radians, kDistance.Radians, 1e-15);
        // CAB must be a right angle with C to the right of AB.
        Assert2.Near(S2.TurnAngle(c, a, b), -S2.M_PI_2 /*radians*/, 1e-15);
    }
}
