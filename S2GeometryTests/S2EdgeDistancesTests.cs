using System;
using Xunit;

namespace S2Geometry
{
    public class S2EdgeDistancesTests
    {
        [Fact]
        public void Test_S2_GetUpdateMinDistanceMaxError() {
            // Verify that the error is "reasonable" for a sampling of distances.
            CheckUpdateMinDistanceMaxError(0, 1.5e-15);
            CheckUpdateMinDistanceMaxError(1e-8, 1e-15);
            CheckUpdateMinDistanceMaxError(1e-5, 1e-15);
            CheckUpdateMinDistanceMaxError(0.05, 1e-15);
            CheckUpdateMinDistanceMaxError(S2Constants.M_PI_2 - 1e-8, 2e-15);
            CheckUpdateMinDistanceMaxError(S2Constants.M_PI_2, 2e-15);
            CheckUpdateMinDistanceMaxError(S2Constants.M_PI_2 + 1e-8, 2e-15);
            CheckUpdateMinDistanceMaxError(Math.PI - 1e-5, 2e-10);
            CheckUpdateMinDistanceMaxError(Math.PI, 0);
        }

        [Fact]
        public void Test_S2_GetUpdateMinInteriorDistanceMaxError() {
            // Check that the error bound returned by
            // GetUpdateMinInteriorDistanceMaxError() is large enough.
            for (int iter = 0; iter < 10000; ++iter) {
                S2Point a0 = S2Testing.RandomPoint();
                S1Angle len = S1Angle.FromRadians(Math.PI * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
                S2Point a1 = S2EdgeDistances.InterpolateAtDistance(len, a0, S2Testing.RandomPoint());
                // TODO(ericv): If S2Pred.RobustCrossProd() is implemented, then we can
                // also test nearly-antipodal points here.  In theory the error bound can
                // be exceeded when the edge endpoints are antipodal to within 0.8e-13
                // radians, but the only examples found in testing require the endpoints
                // to be nearly-antipodal to within 1e-16 radians.
                S2Point n = S2PointUtil.RobustCrossProd(a0, a1).Normalized;
                double f = Math.Pow(1e-20, S2Testing.Random.RandDouble());
                S2Point a = ((1 - f) * a0 + f * a1).Normalized;
                S1Angle r = S1Angle.FromRadians(S2Constants.M_PI_2 * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
                if (S2Testing.Random.OneIn(2)) r = S1Angle.FromRadians(S2Constants.M_PI_2) - r;
                S2Point x = S2EdgeDistances.InterpolateAtDistance(r, a, n);
                S1ChordAngle min_dist = S1ChordAngle.Infinity;
                if (!S2EdgeDistances.UpdateMinInteriorDistance(x, a0, a1, ref min_dist)) {
                    --iter; continue;
                }
                double error = S2EdgeDistances.GetUpdateMinDistanceMaxError(min_dist);
                Assert.True(S2Pred.CompareEdgeDistance(x, a0, a1,
                                                      min_dist.PlusError(error)) <= 0);
                Assert.True(S2Pred.CompareEdgeDistance(x, a0, a1,
                                                      min_dist.PlusError(-error)) >= 0);
            }
        }

        [Fact]
        public void Test_S2_Distance() {
            CheckDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          0, new S2Point(1, 0, 0));
            CheckDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          0, new S2Point(0, 1, 0));
            CheckDistance(new S2Point(1, 3, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          0, new S2Point(1, 3, 0));
            CheckDistance(new S2Point(0, 0, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          S2Constants.M_PI_2, new S2Point(1, 0, 0));
            CheckDistance(new S2Point(0, 0, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          S2Constants.M_PI_2, new S2Point(1, 0, 0));
            CheckDistance(new S2Point(-1, -1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          0.75 * Math.PI, S2Point.Empty);

            CheckDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(1, 1, 0),
                          S2Constants.M_PI_4, new S2Point(1, 1, 0));
            CheckDistance(new S2Point(0, -1, 0), new S2Point(1, 0, 0), new S2Point(1, 1, 0),
                          S2Constants.M_PI_2, new S2Point(1, 0, 0));

            CheckDistance(new S2Point(0, -1, 0), new S2Point(1, 0, 0), new S2Point(-1, 1, 0),
                          S2Constants.M_PI_2, new S2Point(1, 0, 0));
            CheckDistance(new S2Point(-1, -1, 0), new S2Point(1, 0, 0), new S2Point(-1, 1, 0),
                          S2Constants.M_PI_2, new S2Point(-1, 1, 0));

            CheckDistance(new S2Point(1, 1, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          Math.Asin(Math.Sqrt(1.0/ 3)), new S2Point(1, 1, 0));
            CheckDistance(new S2Point(1, 1, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                          Math.Asin(Math.Sqrt(1.0/ 3)), new S2Point(1, 1, 0));

            CheckDistance(new S2Point(-1, 0, 0), new S2Point(1, 1, 0), new S2Point(1, 1, 0),
                          0.75 * Math.PI, new S2Point(1, 1, 0));
            CheckDistance(new S2Point(0, 0, -1), new S2Point(1, 1, 0), new S2Point(1, 1, 0),
                          S2Constants.M_PI_2, new S2Point(1, 1, 0));
            CheckDistance(new S2Point(-1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                          Math.PI, new S2Point(1, 0, 0));
        }

        [Fact]
        public void Test_S2_DistanceOptimizationIsConservative() {
            // Verifies that AlwaysUpdateMinInteriorDistance() computes the lower bound
            // on the true distance conservatively.  (This test used to fail.)
            S2Point x=new(-0.017952729194524016, -0.30232422079175203, 0.95303607751077712);
            S2Point a=new(-0.017894725505830295, -0.30229974986194175, 0.95304493075220664);
            S2Point b=new(-0.017986591360900289, -0.30233851195954353, 0.95303090543659963);
            S1ChordAngle min_distance = S1ChordAngle.Infinity;
            Assert.True(S2EdgeDistances.UpdateMinDistance(x, a, b, ref min_distance));
            min_distance = min_distance.Successor;
            Assert.True(S2EdgeDistances.UpdateMinDistance(x, a, b, ref min_distance));
        }

        [Fact]
        public void Test_S2_MaxDistance() {
            CheckMaxDistance(new S2Point(1, 0, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             S2Constants.M_PI_2);
            CheckMaxDistance(new S2Point(1, 0, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             S2Constants.M_PI_2);
            CheckMaxDistance(new S2Point(0, 1, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             S2Constants.M_PI_2);
            CheckMaxDistance(new S2Point(0, 1, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             S2Constants.M_PI_2);

            CheckMaxDistance(new S2Point(1, 1, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             Math.Asin(Math.Sqrt(2.0/ 3)));
            CheckMaxDistance(new S2Point(1, 1, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             Math.Asin(Math.Sqrt(2.0/ 3)));

            CheckMaxDistance(new S2Point(1, 0, 0), new S2Point(1, 1, 0), new S2Point(1, -1, 0),
                             S2Constants.M_PI_4);
            CheckMaxDistance(new S2Point(0, 1, 0), new S2Point(1, 1, 0), new S2Point(-1, 1, 0),
                             S2Constants.M_PI_4);
            CheckMaxDistance(new S2Point(0, 0, 1), new S2Point(0, 1, 1), new S2Point(0, -1, 1),
                             S2Constants.M_PI_4);

            CheckMaxDistance(new S2Point(0, 0, 1), new S2Point(1, 0, 0), new S2Point(1, 0, -1),
                             3 * S2Constants.M_PI_4);
            CheckMaxDistance(new S2Point(0, 0, 1), new S2Point(1, 0, 0), new S2Point(1, 1, -S2Constants.M_SQRT2),
                             3 * S2Constants.M_PI_4);

            CheckMaxDistance(new S2Point(0, 0, 1), new S2Point(0, 0, -1), new S2Point(0, 0, -1),
                             Math.PI);
        }

        [Fact]
        public void Test_S2_Interpolate() {
            // Choose test points designed to expose floating-point errors.
            S2Point p1 = new S2Point(0.1, 1e-30, 0.3).Normalized;
            S2Point p2 = new S2Point(-0.7, -0.55, -1e30).Normalized;

            // A zero-length edge.
            CheckInterpolate(0, p1, p1, p1);
            CheckInterpolate(1, p1, p1, p1);

            // Start, end, and middle of a medium-length edge.
            CheckInterpolate(0, p1, p2, p1);
            CheckInterpolate(1, p1, p2, p2);
            CheckInterpolate(0.5, p1, p2, 0.5 * (p1 + p2));

            // Test that interpolation is done using distances on the sphere rather than
            // linear distances.
            CheckInterpolate(1.0/ 3, new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             new S2Point(Math.Sqrt(3), 1, 0));
            CheckInterpolate(2.0/ 3, new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                             new S2Point(1, Math.Sqrt(3), 0));

            // Test that interpolation is accurate on a long edge (but not so long that
            // the definition of the edge itself becomes too unstable).
            {
                double kLng = Math.PI - 1e-2;
                S2Point a = S2LatLng.FromRadians(0, 0).ToPoint();
                S2Point b = S2LatLng.FromRadians(0, kLng).ToPoint();
                for (double f = 0.4; f > 1e-15; f *= 0.1) {
                    CheckInterpolate(f, a, b,
                                     S2LatLng.FromRadians(0, f * kLng).ToPoint());
                    CheckInterpolate(1 - f, a, b,
                                     S2LatLng.FromRadians(0, (1 - f) * kLng).ToPoint());
                }
            }

            // Test that interpolation on a 180 degree edge (antipodal endpoints) yields
            // a result with the correct distance from each endpoint.
            for (double t = 0; t <= 1; t += 0.125) {
                S2Point actual = S2EdgeDistances.Interpolate(t, p1, -p1);
                Assert2.Near(new S1Angle(actual, p1).Radians, t * Math.PI, 3e-15);
            }
        }

        [Fact]
        public void Test_S2_InterpolateCanExtrapolate() {
            S2Point i=new(1, 0, 0);
            S2Point j = new(0, 1, 0);
            // Initial vectors at 90 degrees.
            CheckInterpolate(0, i, j, new S2Point(1, 0, 0));
            CheckInterpolate(1, i, j, new S2Point(0, 1, 0));
            CheckInterpolate(1.5, i, j, new S2Point(-1, 1, 0));
            CheckInterpolate(2, i, j, new S2Point(-1, 0, 0));
            CheckInterpolate(3, i, j, new S2Point(0, -1, 0));
            CheckInterpolate(4, i, j, new S2Point(1, 0, 0));

            // Negative values of t.
            CheckInterpolate(-1, i, j, new S2Point(0, -1, 0));
            CheckInterpolate(-2, i, j, new S2Point(-1, 0, 0));
            CheckInterpolate(-3, i, j, new S2Point(0, 1, 0));
            CheckInterpolate(-4, i, j, new S2Point(1, 0, 0));

            // Initial vectors at 45 degrees.
            CheckInterpolate(2, i, new S2Point(1, 1, 0), new S2Point(0, 1, 0));
            CheckInterpolate(3, i, new S2Point(1, 1, 0), new S2Point(-1, 1, 0));
            CheckInterpolate(4, i, new S2Point(1, 1, 0), new S2Point(-1, 0, 0));

            // Initial vectors at 135 degrees.
            CheckInterpolate(2, i, new S2Point(-1, 1, 0), new S2Point(0, -1, 0));

            // Take a small fraction along the curve.
            S2Point p=S2EdgeDistances.Interpolate(0.001, i, j);
            // We should get back where we started.
            CheckInterpolate(1000, i, p, j);
        }

        [Fact]
        public void Test_S2_RepeatedInterpolation() {
            // Check that points do not drift away from unit length when repeated
            // interpolations are done.
            for (int i = 0; i < 100; ++i) {
                S2Point a = S2Testing.RandomPoint();
                S2Point b = S2Testing.RandomPoint();
                for (int j = 0; j < 1000; ++j) {
                    a = S2EdgeDistances.Interpolate(0.01, a, b);
                }
                Assert.True(a.IsUnitLength);
            }
        }

        [Fact]
        public void Test_S2_EdgePairMinDistance() {
            // One edge is degenerate.
            CheckEdgePairMinDistance(new S2Point(1, 0, 1), new S2Point(1, 0, 1),
                                     new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                     S2Constants.M_PI_4, new S2Point(1, 0, 1), new S2Point(1, 0, 0));
            CheckEdgePairMinDistance(new S2Point(1, -1, 0), new S2Point(1, 1, 0),
                                     new S2Point(1, 0, 1), new S2Point(1, 0, 1),
                                     S2Constants.M_PI_4, new S2Point(1, 0, 0), new S2Point(1, 0, 1));

            // Both edges are degenerate.
            CheckEdgePairMinDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                     new S2Point(0, 1, 0), new S2Point(0, 1, 0),
                                     S2Constants.M_PI_2, new S2Point(1, 0, 0), new S2Point(0, 1, 0));

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
        public void Test_S2_EdgePairMaxDistance() {
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
                                     S2Constants.M_PI_2);

            // Both edges are degenerate and antipodal.
            CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0),
                                     new S2Point(-1, 0, 0), new S2Point(-1, 0, 0),
                                     Math.PI);

            // Two identical edges.
            CheckEdgePairMaxDistance(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                     new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                                     S2Constants.M_PI_2);

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
        public void Test_S2_EdgeBNearEdgeA() {
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
            x = x.Normalized;
            a = a.Normalized;
            b = b.Normalized;
            expected_closest = expected_closest.Normalized;
            Assert2.Near(distance_radians, S2EdgeDistances.GetDistance(x, a, b).Radians, S2Constants.DoubleError);
            S2Point closest = S2EdgeDistances.Project(x, a, b);
            if (expected_closest == S2Point.Empty)
            {
                // This special value says that the result should be A or B.
                Assert.True(closest == a || closest == b);
            }
            else
            {
                Assert.True(S2PointUtil.ApproxEquals(closest, expected_closest));
            }
            S1ChordAngle min_distance = S1ChordAngle.Zero;
            Assert.False(S2EdgeDistances.UpdateMinDistance(x, a, b, ref min_distance));
            min_distance = S1ChordAngle.Infinity;
            Assert.True(S2EdgeDistances.UpdateMinDistance(x, a, b, ref min_distance));
            Assert2.Near(distance_radians, min_distance.ToAngle().Radians, S2Constants.DoubleError);
        }

        private static void CheckMaxDistance(S2Point x, S2Point a, S2Point b, double distance_radians)
        {
            x = x.Normalized;
            a = a.Normalized;
            b = b.Normalized;

            S1ChordAngle max_distance = S1ChordAngle.Straight;
            Assert.False(S2EdgeDistances.UpdateMaxDistance(x, a, b, ref max_distance));
            max_distance = S1ChordAngle.Negative;
            Assert.True(S2EdgeDistances.UpdateMaxDistance(x, a, b, ref max_distance));

            Assert2.Near(distance_radians, max_distance.Radians, S2Constants.DoubleError);
        }

        private static void CheckInterpolate(double t, S2Point a, S2Point b, S2Point expected)
        {
            a = a.Normalized;
            b = b.Normalized;
            expected = expected.Normalized;
            S2Point actual = S2EdgeDistances.Interpolate(t, a, b);

            // We allow a bit more than the usual 1e-15 error tolerance because
            // Interpolate() uses trig functions.
            Assert.True(S2PointUtil.ApproxEquals(expected, actual, S1Angle.FromRadians(3e-15)));
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
            a0 = a0.Normalized;
            a1 = a1.Normalized;
            b0 = b0.Normalized;
            b1 = b1.Normalized;
            expected_a = expected_a.Normalized;
            expected_b = expected_b.Normalized;
            var closest = S2EdgeDistances.GetEdgePairClosestPoints(a0, a1, b0, b1);
            S2Point actual_a = closest.Item1;
            S2Point actual_b = closest.Item2;
            if (expected_a == S2Point.Empty)
            {
                // This special value says that the result should be a0 or a1.
                Assert.True(actual_a == a0 || actual_a == a1);
            }
            else
            {
                Assert.True(S2PointUtil.ApproxEquals(expected_a, actual_a));
            }
            if (expected_b == S2Point.Empty)
            {
                // This special value says that the result should be b0 or b1.
                Assert.True(actual_b == b0 || actual_b == b1);
            }
            else
            {
                Assert.True(S2PointUtil.ApproxEquals(expected_b, actual_b));
            }
            S1ChordAngle min_distance = S1ChordAngle.Zero;
            Assert.False(S2EdgeDistances.UpdateEdgePairMinDistance(a0, a1, b0, b1, ref min_distance));
            min_distance = S1ChordAngle.Infinity;
            Assert.True(S2EdgeDistances.UpdateEdgePairMinDistance(a0, a1, b0, b1, ref min_distance));
            Assert2.Near(distance_radians, min_distance.Radians, S2Constants.DoubleError);
        }

        // Given two edges a0a1 and b0b1, check that the maximum distance between them
        // is "distance_radians".  Parameters are passed by value so that this
        // function can normalize them.
        private static void CheckEdgePairMaxDistance(S2Point a0, S2Point a1, S2Point b0, S2Point b1, double distance_radians)
        {
            a0 = a0.Normalized;
            a1 = a1.Normalized;
            b0 = b0.Normalized;
            b1 = b1.Normalized;

            S1ChordAngle max_distance = S1ChordAngle.Straight;
            Assert.False(S2EdgeDistances.UpdateEdgePairMaxDistance(a0, a1, b0, b1, ref max_distance));
            max_distance = S1ChordAngle.Negative;
            Assert.True(S2EdgeDistances.UpdateEdgePairMaxDistance(a0, a1, b0, b1, ref max_distance));
            Assert2.Near(distance_radians, max_distance.Radians, S2Constants.DoubleError);
        }

        private static bool IsEdgeBNearEdgeA(string a_str, string b_str, double max_error_degrees)
        {
            var a=S2TextFormat.MakePolylineOrDie(a_str);
            Assert.Equal(2, a.NumVertices);
            var b=S2TextFormat.MakePolylineOrDie(b_str);
            Assert.Equal(2, b.NumVertices);
            return S2EdgeDistances.IsEdgeBNearEdgeA(a.Vertex(0), a.Vertex(1),
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
            S1ChordAngle ca=new(S1Angle.FromRadians(actual));
            S1Angle bound = ca.PlusError(S2EdgeDistances.GetUpdateMinDistanceMaxError(ca)).ToAngle();
            Assert.True(bound.Radians - actual <= max_error);
        }
    }
}
