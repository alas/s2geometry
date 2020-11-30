using System;
using Xunit;

namespace S2Geometry
{
    public class S2EdgeCrosserTests
    {
#if !DEBUG

        [Fact]
        public void Test_S2EdgeUtil_InvalidDefaultPoints() {
  // Check that default-constructed S2Point arguments don't cause crashes.
  S2Point point=new(0, 0, 0);
  TestCrossingSignInvalid(point, 0);
  TestEdgeOrVertexCrossingInvalid(point, false);
}

        [Fact]
        public void Test_S2EdgeUtil_InvalidNanPoints() {
  // Check that NaN S2Point arguments don't cause crashes.
  S2Point point=new(double.NaN, double.NaN, double.NaN);
  TestCrossingSignInvalid(point, -1);
  TestEdgeOrVertexCrossingInvalid(point, false);
}

// In non-debug builds, check that default-constructed and/or NaN S2Point
// arguments don't cause crashes, especially on the very first method call
// (since S2CopyingEdgeCrosser checks whether the first vertex of each edge is
// the same as the last vertex of the previous edged when deciding whether or
// not to call Restart).
private void TestCrossingSignInvalid(S2Point point, int expected) {
  S2EdgeCrosser crosser=new(point, point);
  Assert.Equal(expected, crosser.CrossingSign(point, point));
  S2CopyingEdgeCrosser crosser2=new(point, point);
  Assert.Equal(expected, crosser2.CrossingSign(point, point));
}

private void TestEdgeOrVertexCrossingInvalid(S2Point point, bool expected) {
  S2EdgeCrosser crosser=new(point, point);
  Assert.Equal(expected, crosser.EdgeOrVertexCrossing(point, point));
  S2CopyingEdgeCrosser crosser2=new(point, point);
  Assert.Equal(expected, crosser2.EdgeOrVertexCrossing(point, point));
}

#endif

        [Fact]
        public void Test_S2EdgeUtil_Crossings() {
            // The real tests of edge crossings are in s2{loop,polygon}_test,
            // but we do a few simple tests here.

            // Two regular edges that cross.
            TestCrossings(new S2Point(1, 2, 1), new S2Point(1, -3, 0.5),
                          new S2Point(1, -0.5, -3), new S2Point(0.1, 0.5, 3), 1, true);

            // Two regular edges that intersect antipodal points.
            TestCrossings(new S2Point(1, 2, 1), new S2Point(1, -3, 0.5),
                          new S2Point(-1, 0.5, 3), new S2Point(-0.1, -0.5, -3), -1, false);

            // Two edges on the same great circle that start at antipodal points.
            TestCrossings(new S2Point(0, 0, -1), new S2Point(0, 1, 0),
                          new S2Point(0, 0, 1), new S2Point(0, 1, 1), -1, false);

            // Two edges that cross where one vertex is S2PointUtil.Origin.
            TestCrossings(new S2Point(1, 0, 0), S2PointUtil.Origin,
                          new S2Point(1, -0.1, 1), new S2Point(1, 1, -0.1), 1, true);

            // Two edges that intersect antipodal points where one vertex is
            // S2PointUtil.Origin.
            TestCrossings(new S2Point(1, 0, 0), S2PointUtil.Origin,
                          new S2Point(-1, 0.1, -1), new S2Point(-1, -1, 0.1), -1, false);

            // Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
            // and edge CD is further CCW around (2,3,4) than AB.
            TestCrossings(new S2Point(2, 3, 4), new S2Point(-1, 2, 5),
                          new S2Point(7, -2, 3), new S2Point(2, 3, 4), 0, false);

            // Two edges that barely cross each other near the middle of one edge.  The
            // edge AB is approximately in the x=y plane, while CD is approximately
            // perpendicular to it and ends exactly at the x=y plane.
            TestCrossings(new S2Point(1, 1, 1), new S2Point(1, MathUtils.NextAfter(1, 0), -1),
                          new S2Point(11, -12, -1), new S2Point(10, 10, 1), 1, true);

            // In this version, the edges are separated by a distance of about 1e-15 (S2Constants.DoubleError).
            TestCrossings(new S2Point(1, 1, 1), new S2Point(1, MathUtils.NextAfter(1, 2), -1),
                          new S2Point(1, -1, 0), new S2Point(1, 1, 0), -1, false);

            // Two edges that barely cross each other near the end of both edges.  This
            // example cannot be handled using regular double-precision arithmetic due
            // to floating-point underflow.
            TestCrossings(new S2Point(0, 0, 1), new S2Point(2, -1e-323, 1),
                          new S2Point(1, -1, 1), new S2Point(1e-323, 0, 1), 1, true);

            // In this version, the edges are separated by a distance of about 1e-640.
            TestCrossings(new S2Point(0, 0, 1), new S2Point(2, 1e-323, 1),
                          new S2Point(1, -1, 1), new S2Point(1e-323, 0, 1), -1, false);

            // Two edges that barely cross each other near the middle of one edge.
            // Computing the exact determinant of some of the triangles in this test
            // requires more than 2000 bits of precision.
            TestCrossings(new S2Point(1, -1e-323, -1e-323), new S2Point(1e-323, 1, 1e-323),
                          new S2Point(1, -1, 1e-323), new S2Point(1, 1, 0), 1, true);

            // In this version, the edges are separated by a distance of about 1e-640.
            TestCrossings(new S2Point(1, 1e-323, -1e-323), new S2Point(-1e-323, 1, 1e-323),
                          new S2Point(1, -1, 1e-323), new S2Point(1, 1, 0), -1, false);
        }

        [Fact]
        public void Test_S2EdgeUtil_CollinearEdgesThatDontTouch() {
            int kIters = 500;
            for (int iter = 0; iter < kIters; ++iter) {
                S2Point a = S2Testing.RandomPoint();
                S2Point d = S2Testing.RandomPoint();
                S2Point b = S2EdgeDistances.Interpolate(0.05, a, d);
                S2Point c = S2EdgeDistances.Interpolate(0.95, a, d);
                Assert.True(0 > S2EdgeCrossings.CrossingSign(a, b, c, d));
                Assert.True(0 > S2EdgeCrossings.CrossingSign(a, b, c, d));
                S2EdgeCrosser crosser = new(a, b, c);
                Assert.True(0 > crosser.CrossingSign(d));
                Assert.True(0 > crosser.CrossingSign(c));
            }
        }

        [Fact]
        public void Test_S2EdgeUtil_CoincidentZeroLengthEdgesThatDontTouch() {
            // It is important that the edge primitives can handle vertices that exactly
            // exactly proportional to each other, i.e. that are not identical but are
            // nevertheless exactly coincident when projected onto the unit sphere.
            // There are various ways that such points can arise.  For example,
            // Normalize() itself is not idempotent: there exist distinct points A,B
            // such that Normalize(A) == B  and Normalize(B) == A.  Another issue is
            // that sometimes calls to Normalize() are skipped when the result of a
            // calculation "should" be unit length mathematically (e.g., when computing
            // the cross product of two orthonormal vectors).
            //
            // This test checks pairs of edges AB and CD where A,B,C,D are exactly
            // coincident on the sphere and the norms of A,B,C,D are monotonically
            // increasing.  Such edge pairs should never intersect.  (This is not
            // obvious, since it depends on the particular symbolic perturbations used
            // by S2Pred.Sign().  It would be better to replace this with a test that
            // says that the CCW results must be consistent with each other.)
            int kIters = 1000;
            for (int iter = 0; iter < kIters; ++iter) {
                // Construct a point P where every component is zero or a power of 2.
                var t = new double[3];
                for (int i = 0; i < 3; ++i) {
                    int binary_exp = S2Testing.Random.Skewed(11);
                    t[i] = (binary_exp > 1022) ? 0 : Math.Pow(2, -binary_exp);
                }
                // If all components were zero, try again.  Note that normalization may
                // convert a non-zero point into a zero one due to underflow (!)
                var p = new S2Point(t).Normalized;
                if (p == S2Point.Empty) { --iter; continue; }

                // Now every non-zero component should have exactly the same mantissa.
                // This implies that if we scale the point by an arbitrary factor, every
                // non-zero component will still have the same mantissa.  Scale the points
                // so that they are all distinct and are still very likely to satisfy
                // S2.IsUnitLength (which allows for a small amount of error in the norm).
                S2Point a = (1 - 3e-16) * p;
                S2Point b = (1 - 1e-16) * p;
                S2Point c = p;
                S2Point d = (1 + 2e-16) * p;
                if (!a.IsUnitLength || !d.IsUnitLength) {
                    --iter;
                    continue;
                }
                // Verify that the expected edges do not cross.
                Assert.True(0 > S2EdgeCrossings.CrossingSign(a, b, c, d));
                S2EdgeCrosser crosser = new(a, b, c);
                Assert.True(0 > crosser.CrossingSign(d));
                Assert.True(0 > crosser.CrossingSign(c));
            }
        }

        private static void TestCrossing(S2Point a, S2Point b,
                          S2Point c, S2Point d,
                          int robust, bool edge_or_vertex) {
            // Modify the expected result if two vertices from different edges match.
            if (a == c || a == d || b == c || b == d) robust = 0;
            Assert.Equal(robust, S2EdgeCrossings.CrossingSign(a, b, c, d));
            S2EdgeCrosser crosser = new(a, b, c);
            Assert.Equal(robust, crosser.CrossingSign(d));
            Assert.Equal(robust, crosser.CrossingSign(c));
            Assert.Equal(robust, crosser.CrossingSign(d, c));
            Assert.Equal(robust, crosser.CrossingSign(c, d));

            Assert.Equal(edge_or_vertex, S2EdgeCrossings.EdgeOrVertexCrossing(a, b, c, d));
            crosser.RestartAt(c);
            Assert.Equal(edge_or_vertex, crosser.EdgeOrVertexCrossing(d));
            Assert.Equal(edge_or_vertex, crosser.EdgeOrVertexCrossing(c));
            Assert.Equal(edge_or_vertex, crosser.EdgeOrVertexCrossing(d, c));
            Assert.Equal(edge_or_vertex, crosser.EdgeOrVertexCrossing(c, d));

            // Check that the crosser can be re-used.
            crosser.Init(c, d);
            crosser.RestartAt(a);
            Assert.Equal(robust, crosser.CrossingSign(b));
            Assert.Equal(robust, crosser.CrossingSign(a));

            // Now try all the same tests with CopyingEdgeCrosser.
            S2CopyingEdgeCrosser crosser2 = new(a, b, c);
            Assert.Equal(robust, crosser2.CrossingSign(d));
            Assert.Equal(robust, crosser2.CrossingSign(c));
            Assert.Equal(robust, crosser2.CrossingSign(d, c));
            Assert.Equal(robust, crosser2.CrossingSign(c, d));

            Assert.Equal(edge_or_vertex, S2EdgeCrossings.EdgeOrVertexCrossing(a, b, c, d));
            crosser2.RestartAt(c);
            Assert.Equal(edge_or_vertex, crosser2.EdgeOrVertexCrossing(d));
            Assert.Equal(edge_or_vertex, crosser2.EdgeOrVertexCrossing(c));
            Assert.Equal(edge_or_vertex, crosser2.EdgeOrVertexCrossing(d, c));
            Assert.Equal(edge_or_vertex, crosser2.EdgeOrVertexCrossing(c, d));

            // Check that the crosser can be re-used.
            crosser2.Init(c, d);
            crosser2.RestartAt(a);
            Assert.Equal(robust, crosser2.CrossingSign(b));
            Assert.Equal(robust, crosser2.CrossingSign(a));
        }

        private static void TestCrossings(S2Point a, S2Point b, S2Point c, S2Point d,
                           int robust, bool edge_or_vertex) {
            a = a.Normalized;
            b = b.Normalized;
            c = c.Normalized;
            d = d.Normalized;
            TestCrossing(a, b, c, d, robust, edge_or_vertex);
            TestCrossing(b, a, c, d, robust, edge_or_vertex);
            TestCrossing(a, b, d, c, robust, edge_or_vertex);
            TestCrossing(b, a, d, c, robust, edge_or_vertex);
            TestCrossing(a, a, c, d, -1, false);
            TestCrossing(a, b, c, c, -1, false);
            TestCrossing(a, a, c, c, -1, false);
            TestCrossing(a, b, a, b, 0, true);
            TestCrossing(c, d, a, b, robust, edge_or_vertex != (robust == 0));
        }
    }
}
