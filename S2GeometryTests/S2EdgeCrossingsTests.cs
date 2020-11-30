using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;
using Xunit.Abstractions;
using IntersectionMethod = S2Geometry.S2EdgeCrossings.Internal.IntersectionMethod;

namespace S2Geometry
{
    // CrossingSign, VertexCrossing, and EdgeOrVertexCrossing are tested in
    // S2EdgeCrosserTests.
    public class S2EdgeCrossingsTests
    {
        // The approximate maximum error in GetDistance() for small distances.
        private static readonly S1Angle kGetDistanceAbsError = S1Angle.FromRadians(3 * S2Constants.DoubleEpsilon);
        private readonly ITestOutputHelper _logger;

        public S2EdgeCrossingsTests(ITestOutputHelper logger) { _logger = logger; }

        [Fact]
        public void Test_S2EdgeUtil_IntersectionError() {
            // We repeatedly construct two edges that cross near a random point "p", and
            // measure the distance from the actual intersection point "x" to the
            // exact intersection point and also to the edges.

            var tally_ = S2EdgeCrossings.Internal.GetNewTally();
            S1Angle max_point_dist = S1Angle.Zero, max_edge_dist = S1Angle.Zero;
            for (int iter = 0; iter < 5000; ++iter) {
                // We construct two edges AB and CD that intersect near "p".  The angle
                // between AB and CD (expressed as a slope) is chosen randomly between
                // 1e-15 and 1e15 such that its logarithm is uniformly distributed.
                // Similarly, two edge lengths approximately between 1e-15 and 1 are
                // chosen.  The edge endpoints are chosen such that they are often very
                // close to the other edge (i.e., barely crossing).  Taken together this
                // ensures that we test both long and very short edges that intersect at
                // both large and very small angles.
                //
                // Sometimes the edges we generate will not actually cross, in which case
                // we simply try again.
                S2Testing.GetRandomFrame(out var p, out var d1, out var d2);
                double slope = 1e-15 * Math.Pow(1e30, S2Testing.Random.RandDouble());
                d2 = (d1 + slope * d2).Normalized;
                S2Point a, b, c, d;
                do {
                    double ab_len = Math.Pow(1e-15, S2Testing.Random.RandDouble());
                    double cd_len = Math.Pow(1e-15, S2Testing.Random.RandDouble());
                    double a_fraction = Math.Pow(1e-5, S2Testing.Random.RandDouble());
                    if (S2Testing.Random.OneIn(2)) a_fraction = 1 - a_fraction;
                    double c_fraction = Math.Pow(1e-5, S2Testing.Random.RandDouble());
                    if (S2Testing.Random.OneIn(2)) c_fraction = 1 - c_fraction;
                    a = (p - a_fraction * ab_len * d1).Normalized;
                    b = (p + (1 - a_fraction) * ab_len * d1).Normalized;
                    c = (p - c_fraction * cd_len * d2).Normalized;
                    d = (p + (1 - c_fraction) * cd_len * d2).Normalized;
                } while (S2EdgeCrossings.CrossingSign(a, b, c, d) <= 0);

                // Each constructed edge should be at most 1.5 * S2Constants.DoubleEpsilon away from the
                // original point P.
                Assert.True(S2EdgeDistances.GetDistance(p, a, b) <=
                      S1Angle.FromRadians(1.5 * S2Constants.DoubleEpsilon) + kGetDistanceAbsError);
                Assert.True(S2EdgeDistances.GetDistance(p, c, d) <=
                          S1Angle.FromRadians(1.5 * S2Constants.DoubleEpsilon) + kGetDistanceAbsError);

                // Verify that the expected intersection point is close to both edges and
                // also close to the original point P.  (It might not be very close to P
                // if the angle between the edges is very small.)
                S2Point expected = GetIntersectionExact(a, b, c, d);
                Assert.True(S2EdgeDistances.GetDistance(expected, a, b) <=
                          S1Angle.FromRadians(3 * S2Constants.DoubleEpsilon) + kGetDistanceAbsError);
                Assert.True(S2EdgeDistances.GetDistance(expected, c, d) <=
                          S1Angle.FromRadians(3 * S2Constants.DoubleEpsilon) + kGetDistanceAbsError);
                Assert.True(new S1Angle(expected, p) <=
                          S1Angle.FromRadians(3 * S2Constants.DoubleEpsilon / slope) +
                          S2EdgeCrossings.kIntersectionErrorS1Angle);

                // Now we actually test the GetIntersection() method.
                S2Point actual = S2EdgeCrossings.GetIntersection(a, b, c, d, tally_);
                S1Angle dist_ab = S2EdgeDistances.GetDistance(actual, a, b);
                S1Angle dist_cd = S2EdgeDistances.GetDistance(actual, c, d);
                Assert.True(dist_ab <= S2EdgeCrossings.kIntersectionErrorS1Angle + kGetDistanceAbsError);
                Assert.True(dist_cd <= S2EdgeCrossings.kIntersectionErrorS1Angle + kGetDistanceAbsError);
                max_edge_dist = new[] { max_edge_dist, dist_ab, dist_cd }.Max();
                S1Angle point_dist = new(expected, actual);
                Assert.True(point_dist <= S2EdgeCrossings.kIntersectionErrorS1Angle);
                max_point_dist = new[] { max_point_dist, point_dist }.Max();
            }
            PrintIntersectionStats(tally_);
        }

        [Fact]
        public void Test_S2EdgeUtil_GrazingIntersections() {
            // This test choose 5 points along a great circle (i.e., as collinear as
            // possible), and uses them to construct an edge AB and a triangle CDE such
            // that CD and CE both cross AB.  It then checks that the intersection
            // points returned by GetIntersection() have the correct relative ordering
            // along AB (to within kIntersectionError).
            var tally_ = S2EdgeCrossings.Internal.GetNewTally();
            for (int iter = 0; iter < 1000; ++iter) {
                S2Testing.GetRandomFrame(out var x, out var y, out _);
                S2Point a, b, c, d, e, ab;
                do {
                    a = ChooseSemicirclePoint(x, y);
                    b = ChooseSemicirclePoint(x, y);
                    c = ChooseSemicirclePoint(x, y);
                    d = ChooseSemicirclePoint(x, y);
                    e = ChooseSemicirclePoint(x, y);
                    ab = (a - b).CrossProd(a + b);
                } while (ab.Norm < 50 * S2Constants.DoubleEpsilon ||
                         S2EdgeCrossings.CrossingSign(a, b, c, d) <= 0 ||
                         S2EdgeCrossings.CrossingSign(a, b, c, e) <= 0);
                S2Point xcd = S2EdgeCrossings.GetIntersection(a, b, c, d, tally_);
                S2Point xce = S2EdgeCrossings.GetIntersection(a, b, c, e, tally_);
                // Essentially this says that if CDE and CAB have the same orientation,
                // then CD and CE should intersect along AB in that order.
                ab = ab.Normalized;
                if (new S1Angle(xcd, xce) > 2 * S2EdgeCrossings.kIntersectionErrorS1Angle) {
                    Assert.Equal(S2Pred.Sign(c, d, e) == S2Pred.Sign(c, a, b),
                              S2Pred.Sign(ab, xcd, xce) > 0);
                }
            }
            PrintIntersectionStats(tally_);
        }

        [Fact]
        public void Test_S2EdgeUtil_ExactIntersectionUnderflow() {
            // Tests that a correct intersection is computed even when two edges are
            // exactly collinear and the normals of both edges underflow in double
            // precision when normalized (see S2PointFromExact function for details).
            S2Point a0=new(1, 0, 0), a1 = new(1, 2e-300, 0);
            S2Point b0 = new(1, 1e-300, 0), b1 = new(1, 3e-300, 0);
            Assert.Equal(new S2Point(1, 1e-300, 0), S2EdgeCrossings.GetIntersection(a0, a1, b0, b1, null));
        }

        [Fact]
        public void Test_S2EdgeUtil_GetIntersectionInvariants() {
            // Test that the result of GetIntersection does not change when the edges
            // are swapped and/or reversed.  The number of iterations is high because it
            // is difficult to generate test cases that show that CompareEdges() is
            // necessary and correct, for example.
            int kIters = 50000;
#if DEBUG
            kIters = 5000;
#endif
            for (int iter = 0; iter < kIters; ++iter) {
                S2Point a, b, c, d;
                do {
                    // GetIntersectionStable() sorts the two edges by length, soruct
                    // edges (a,b) and (c,d) that cross and have exactly the same length.
                    // This can be done by swapping the "x" and "y" coordinates.
                    // [Swapping other coordinate pairs doesn't work because it changes the
                    // order of addition in Norm2() == (x**2 + y**2) + z**2.]
                    a = c = S2Testing.RandomPoint();
                    b = d = S2Testing.RandomPoint();
                    c = new S2Point(c.Y, c.X, c.Z);
                    d = new S2Point(d.Y, d.X, d.Z);
                } while (S2EdgeCrossings.CrossingSign(a, b, c, d) <= 0);
                Assert.Equal((a - b).Norm2, (c - d).Norm2);

                // Now verify that GetIntersection returns exactly the same result when
                // the edges are swapped and/or reversed.
                S2Point result = S2EdgeCrossings.GetIntersection(a, b, c, d, null);
                if (S2Testing.Random.OneIn(2)) { var tmp = a; a = b; b = tmp; }
                if (S2Testing.Random.OneIn(2)) { var tmp = c; c = d; d = tmp; }
                if (S2Testing.Random.OneIn(2)) { var tmp = a; b = d; d = tmp; }
                Assert.Equal(result, S2EdgeCrossings.GetIntersection(a, b, c, d, null));
            }
        }

        // This returns the true intersection point of two line segments (a0,a1) and
        // (b0,b1), with a relative error of at most S2Constants.DoubleEpsilon in each coordinate
        // (i.e., one ulp, or twice the double precision rounding error).
        private static S2Point GetIntersectionExact(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
        {
            S2Point x = S2EdgeCrossings.Internal.GetIntersectionExact(a0, a1, b0, b1);
            if (x.DotProd((a0 + a1) + (b0 + b1)) < 0) x = -x;
            return x;
        }

        // Chooses a point in the XY plane that is separated from X by at least 1e-15
        // (to avoid choosing too many duplicate points) and by at most Pi/2 - 1e-3
        // (to avoid nearly-diametric edges, since the test below is not sophisticated
        // enough to test such edges).
        private static S2Point ChooseSemicirclePoint(S2Point x, S2Point y)
        {
            double sign = (2 * S2Testing.Random.Uniform(2)) - 1;
            return (x + sign * 1e3 * Math.Pow(1e-18, S2Testing.Random.RandDouble()) * y).Normalized;
        }

        // This method records statistics about the intersection methods used by
        // GetIntersection().
        private void PrintIntersectionStats(Dictionary<IntersectionMethod, int> tally_)
        {
            int total = 0;
            var totals = S2EdgeCrossings.Internal.GetNewTally();
            foreach (var key in tally_.Keys)
            {
                total += tally_[key];
                totals[key] = total;
            }
            _logger.WriteLine($"{"Method":10} {"Successes":16} {"Attempts":16}  {"Rate":6}");
            foreach (var key in tally_.Keys)
            {
                var tal = tally_[key];
                if (tal == 0) continue;
                var name = S2EdgeCrossings.Internal.GetIntersectionMethodName(key);
                var tot = totals[key];
                var suc = 100.0 * tal / total;
                var att = 100.0 * tot / total;
                var rat = 100.0 * tal / tot;
                _logger.WriteLine($"{name:10} {tal:9} {suc:5.1f}% {tot:9} {att:5.1f}%  {rat:5.1f}%");
            }
            foreach (var key in tally_.Keys)
                tally_[key] = 0;
        }
    }
}
