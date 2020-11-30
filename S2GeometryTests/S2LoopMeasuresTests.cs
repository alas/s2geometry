using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;
using static S2Geometry.S2TextFormat;

namespace S2Geometry
{
    public class S2LoopMeasuresTests
    {
        // Some standard loops to use in the tests (see descriptions below).
        private readonly S2Point[] full_;
        private readonly S2Point[] v_loop_;
        private readonly S2Point[] north_hemi_;
        private readonly S2Point[] north_hemi3_;
        private readonly S2Point[] west_hemi_;
        private readonly S2Point[] east_hemi_;
        private readonly S2Point[] candy_cane_;
        private readonly S2Point[] line_triangle_;
        private readonly S2Point[] skinny_chevron_;
        private readonly S2Point[] three_leaf_clover_;
        private readonly S2Point[] tessellated_loop_;

        public S2LoopMeasuresTests() {
            // The full loop is represented as a loop with no vertices.
            full_ = Array.Empty<S2Point>();

            // A degenerate loop in the shape of a "V".
            v_loop_ = ParsePointsOrDie("5:1, 0:2, 5:3, 0:2").ToArray();

            // The northern hemisphere, defined using two pairs of antipodal points.
            north_hemi_ = ParsePointsOrDie("0:-180, 0:-90, 0:0, 0:90").ToArray();

            // The northern hemisphere, defined using three points 120 degrees apart.
            north_hemi3_ = ParsePointsOrDie("0:-180, 0:-60, 0:60").ToArray();

            // The western hemisphere, defined using two pairs of antipodal points.
            west_hemi_ = ParsePointsOrDie("0:-180, -90:0, 0:0, 90:0").ToArray();

            // The eastern hemisphere, defined using two pairs of antipodal points.
            east_hemi_ = ParsePointsOrDie("90:0, 0:0, -90:0, 0:-180").ToArray();

            // A spiral stripe that slightly over-wraps the equator.
            candy_cane_ = ParsePointsOrDie(
                "-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70").ToArray();

            // A completely degenerate triangle along the equator that Sign()
            // considers to be CCW.
            line_triangle_ = ParsePointsOrDie("0:1, 0:2, 0:3").ToArray();

            // A nearly-degenerate CCW chevron near the equator with very long sides
            // (about 80 degrees).  Its area is less than 1e-640, which is too small
            // to represent in double precision.
            skinny_chevron_ = ParsePointsOrDie("0:0, -1e-320:80, 0:1e-320, 1e-320:80").ToArray();

            // A loop where the same vertex appears three times.
            three_leaf_clover_ = ParsePointsOrDie(
                "0:0, -3:3, 3:3, 0:0, 3:0, 3:-3, 0:0, -3:-3, -3:0").ToArray();

            // A loop with groups of 3 or more vertices in a straight line.
            tessellated_loop_ = ParsePointsOrDie(
                "10:34, 5:34, 0:34, -10:34, -10:36, -5:36, 0:36, 10:36").ToArray();
        }

        [Fact]
        public void Test_PruneDegeneracies_AllDegeneracies()
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
        public void Test_PruneDegeneracies_SomeDegeneracies()
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
        public void Test_GetCanonicalLoopOrder_AllDegeneracies()
        {
            TestCanonicalLoopOrder("", new S2LoopMeasures.LoopOrder(0, 1));
            TestCanonicalLoopOrder("a", new S2LoopMeasures.LoopOrder(0, 1));
            TestCanonicalLoopOrder("aaaaa", new S2LoopMeasures.LoopOrder(0, 1));
            TestCanonicalLoopOrder("ba", new S2LoopMeasures.LoopOrder(1, 1));
            TestCanonicalLoopOrder("bab", new S2LoopMeasures.LoopOrder(1, 1));
            TestCanonicalLoopOrder("cbab", new S2LoopMeasures.LoopOrder(2, 1));
            TestCanonicalLoopOrder("bacbcab", new S2LoopMeasures.LoopOrder(8, -1));
        }

        [Fact]
        public void Test_GetPerimeter_Empty() {
            Assert.Equal(S1Angle.Zero, S2LoopMeasures.GetPerimeter(Array.Empty<S2Point>()));
        }

        [Fact]
        public void Test_GetPerimeter_Octant() {
            var loop = ParsePointsOrDie("0:0, 0:90, 90:0");
            Assert2.Near(3 * S2Constants.M_PI_2, S2LoopMeasures.GetPerimeter(loop.ToArray()).Radians);
        }

        [Fact]
        public void Test_GetPerimeter_MoreThanTwoPi() {
            // Make sure that GetPerimeter doesn't use S1ChordAngle, which can only
            // represent distances up to 2*Pi.
            var loop = ParsePointsOrDie("0:0, 0:90, 0:180, 90:0, 0:-90").ToArray();
            Assert2.Near(5 * S2Constants.M_PI_2, S2LoopMeasures.GetPerimeter(loop).Radians);
        }

        [Fact]
        public void Test_LoopTestBase_GetAreaConsistentWithCurvature() {
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
        public void Test_LoopTestBase_GetAreaConsistentWithOrientation() {
            // Test that GetArea() returns an area near 0 for degenerate loops that
            // contain almost no points, and an area near 4*Pi for degenerate loops that
            // contain almost all points.

            const int kMaxVertices = 6;
            for (int i = 0; i < 50; ++i) {
                int num_vertices = 3 + S2Testing.Random.Uniform(kMaxVertices - 3 + 1);
                // Repeatedly choose N vertices that are exactly on the equator until we
                // find some that form a valid loop.
                var loopTmp = new List<S2Point>();
                do {
                    for (int i2 = 0; i2 < num_vertices; ++i2) {
                        // We limit longitude to the range [0, 90] to ensure that the loop is
                        // degenerate (as opposed to following the entire equator).
                        loopTmp.Add(
                            S2LatLng.FromRadians(0, S2Testing.Random.RandDouble() * S2Constants.M_PI_2).ToPoint());
                    }
                } while (!new S2Loop(loopTmp, false).IsValid);
                var loop = loopTmp.ToArray();
                bool ccw = S2LoopMeasures.IsNormalized(loop);
                // The error bound is sufficient for current tests but not guaranteed.
                _ = i + ": " + loop.ToDebugString();
                Assert2.Near(ccw ? 0 : S2Constants.M_4_PI, S2LoopMeasures.GetArea(loop), 1e-14);
                Assert.Equal(!ccw, new S2Loop(loop).Contains(new S2Point(0, 0, 1)));
            }
        }

        [Fact]
        public void Test_LoopTestBase_GetAreaAccuracy() {
            // TODO(ericv): Test that GetArea() has an accuracy significantly better
            // than 1e-15 on loops whose area is small.
        }

        [Fact]
        public void Test_LoopTestBase_GetAreaAndCentroid() {
            Assert.Equal(S2Constants.M_4_PI, S2LoopMeasures.GetArea(full_));
            Assert.Equal(S2Point.Empty, S2LoopMeasures.GetCentroid(full_));

            Assert2.Near(S2LoopMeasures.GetArea(north_hemi_), S2Constants.M_2_PI);
            Assert2.Near(S2Constants.M_2_PI, S2LoopMeasures.GetArea(east_hemi_), 1e-12);

            // Construct spherical caps of random height, and approximate their boundary
            // with closely spaces vertices.  Then check that the area and centroid are
            // correct.
            for (int iter = 0; iter < 50; ++iter) {
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

                var loopTmp = new List<S2Point>();
                for (double theta = 0; theta < S2Constants.M_2_PI;
                     theta += S2Testing.Random.RandDouble() * max_dtheta) {
                    loopTmp.Add(Math.Cos(theta) * Math.Cos(phi) * x +
                                   Math.Sin(theta) * Math.Cos(phi) * y +
                                   Math.Sin(phi) * z);
                }
                var loop = loopTmp.ToArray();
                double area = S2LoopMeasures.GetArea(loop);
                S2Point centroid = S2LoopMeasures.GetCentroid(loop);
                double expected_area = S2Constants.M_2_PI * height;
                Assert.True(Math.Abs(area - expected_area) <= S2Constants.M_2_PI * kMaxDist);
                S2Point expected_centroid = expected_area * (1 - 0.5 * height) * z;
                Assert.True((centroid - expected_centroid).Norm <= 2 * kMaxDist);
            }
        }

        [Fact]
        public void Test_LoopTestBase_GetCurvature() {
            Assert.Equal(-S2Constants.M_2_PI, S2LoopMeasures.GetCurvature(full_));

            Assert.Equal(S2Constants.M_2_PI, S2LoopMeasures.GetCurvature(v_loop_));
            CheckCurvatureInvariants(v_loop_);

            // This curvature should be computed exactly.
            Assert.Equal(0, S2LoopMeasures.GetCurvature(north_hemi3_));
            CheckCurvatureInvariants(north_hemi3_);

            Assert2.Near(0, S2LoopMeasures.GetCurvature(west_hemi_), 1e-15);
            CheckCurvatureInvariants(west_hemi_);

            // We don't have an easy way to estimate the curvature of these loops, but
            // we can still check that the expected invariants hold.
            CheckCurvatureInvariants(candy_cane_);
            CheckCurvatureInvariants(three_leaf_clover_);

            Assert2.Near(S2Constants.M_2_PI, S2LoopMeasures.GetCurvature(line_triangle_));
            CheckCurvatureInvariants(line_triangle_);

            Assert2.Near(S2Constants.M_2_PI, S2LoopMeasures.GetCurvature(skinny_chevron_));
            CheckCurvatureInvariants(skinny_chevron_);

            // Build a narrow spiral loop starting at the north pole.  This is designed
            // to test that the error in GetCurvature is linear in the number of
            // vertices even when the partial sum of the curvatures gets very large.
            // The spiral consists of two "arms" defining opposite sides of the loop.
            // This is a pathological loop that contains many long parallel edges.
            int kArmPoints = 10000;    // Number of vertices in each "arm"
            double kArmRadius = 0.01;  // Radius of spiral.
            var spiral = new S2Point[2 * kArmPoints];
            spiral[kArmPoints] = new S2Point(0, 0, 1);
            for (int i = 0; i < kArmPoints; ++i) {
                double angle = (S2Constants.M_2_PI / 3) * i;
                double x = Math.Cos(angle);
                double y = Math.Sin(angle);
                double r1 = i * kArmRadius / kArmPoints;
                double r2 = (i + 1.5) * kArmRadius / kArmPoints;
                spiral[kArmPoints - i - 1] = new S2Point(r1 * x, r1 * y, 1).Normalized;
                spiral[kArmPoints + i] = new S2Point(r2 * x, r2 * y, 1).Normalized;
            }

            // Check that GetCurvature() is consistent with GetArea() to within the
            // error bound of the former.  We actually use a tiny fraction of the
            // worst-case error bound, since the worst case only happens when all the
            // roundoff errors happen in the same direction and this test is not
            // designed to achieve that.  The error in GetArea() can be ignored for the
            // purposes of this test since it is generally much smaller.
            Assert2.Near(
                S2Constants.M_2_PI - S2LoopMeasures.GetArea(spiral),
                S2LoopMeasures.GetCurvature(spiral),
                0.01 * S2LoopMeasures.GetCurvatureMaxError(spiral));
        }

        // Given a string where each character "ch" represents a vertex (such as
        // "abac"), returns a vector of S2Points of the form (ch, 0, 0).  Note that
        // these points are not unit length and therefore are not suitable for general
        // use, however they are useful for testing certain functions below.
        private S2Point[] MakeTestLoop(string loop_str)
        {
            return S2Testing.StrPoints.StrToPoints(loop_str);
        }

        // Given a loop whose vertices are represented as characters (such as "abcd" or
        // "abccb"), verify that S2.PruneDegeneracies() yields the loop "expected".
        private void TestPruneDegeneracies(string input_str, string expected_str)
        {
            var input = MakeTestLoop(input_str);
            var actual_str = S2Testing.StrPoints.PointsToStr(S2LoopMeasures.PruneDegeneracies(input));
            Assert.Equal(expected_str, actual_str);
        }

        // Given a loop whose vertices are represented as characters (such as "abcd" or
        // "abccb"), verify that S2.GetCanonicalLoopOrder returns the given result.
        private void TestCanonicalLoopOrder(string input_str, S2LoopMeasures.LoopOrder expected_order) {
            Assert.Equal(expected_order, S2LoopMeasures.GetCanonicalLoopOrder(MakeTestLoop(input_str)));
        }

        private static void TestAreaConsistentWithCurvature(S2Point[] loop) {
            // Check that the area computed using GetArea() is consistent with the loop
            // curvature.  According to the Gauss-Bonnet theorem, the area of the loop
            // equals 2*Pi minus its curvature.
            double area = S2LoopMeasures.GetArea(loop);
            double gauss_area = S2Constants.M_2_PI - S2LoopMeasures.GetCurvature(loop);
            // The error bound below is sufficient for current tests but not guaranteed.
            _ = loop.ToDebugString();
            Assert.True(Math.Abs(area - gauss_area) <= 1e-14);
        }

        private static void ExpectSameOrder(S2Point[] loop1, S2LoopMeasures.LoopOrder order1, S2Point[] loop2, S2LoopMeasures.LoopOrder order2) {
            Assert.Equal(loop1.Length, loop2.Length);
            int i1 = order1.first, i2 = order2.first;
            int dir1 = order1.dir, dir2 = order2.dir;
            for (int n = loop1.Length; --n >= 0;) {
                Assert.Equal(loop1[i1], loop2[i2]);
                i1 += dir1;
                i2 += dir2;
            }
        }

        // Check that the curvature is *identical* when the vertex order is
        // rotated, and that the sign is inverted when the vertices are reversed.
        private static void CheckCurvatureInvariants(S2Point[] loop_in) {
            S2LoopMeasures.LoopOrder order_in = S2LoopMeasures.GetCanonicalLoopOrder(loop_in);
            var loop = loop_in;
            double expected = S2LoopMeasures.GetCurvature(loop);
            for (int i = 0; i < loop.Length; ++i) {
                loop.Reverse();
                Assert.Equal((expected == S2Constants.M_2_PI) ? expected : -expected,
                          S2LoopMeasures.GetCurvature(loop));
                ExpectSameOrder(loop_in, order_in, loop, S2LoopMeasures.GetCanonicalLoopOrder(loop));
                loop.Reverse();
                loop.RotateInPlace(1);
                Assert.Equal(expected, S2LoopMeasures.GetCurvature(loop));
                ExpectSameOrder(loop_in, order_in, loop, S2LoopMeasures.GetCanonicalLoopOrder(loop));
            }
        }
    }
}
