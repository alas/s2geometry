using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using Xunit;
using Xunit.Abstractions;

namespace S2Geometry
{
    public class S2PolylineTests
    {
        private readonly ITestOutputHelper _logger;

        public S2PolylineTests(ITestOutputHelper logger) { _logger = logger; }

        [Fact]
        public void Test_S2Polyline_Basic() {
            var vertices = Array.Empty<S2Point>();
            S2Polyline empty = new S2Polyline(vertices);
            Assert.Equal(S2LatLngRect.Empty, empty.GetRectBound());
            empty.Reverse();
            Assert.Equal(0, empty.NumVertices);

            S2LatLng[] latlngs = {
                S2LatLng.FromDegrees(0, 0),
                S2LatLng.FromDegrees(0, 90),
                S2LatLng.FromDegrees(0, 180)};
            S2Polyline semi_equator = new S2Polyline(latlngs);
            Assert.True(S2PointUtil.ApproxEquals(semi_equator.Interpolate(0.5),
                                         new S2Point(0, 1, 0)));
            semi_equator.Reverse();
            Assert.Equal(new S2Point(1, 0, 0), semi_equator.Vertex(2));
        }

        [Fact]
        public void Test_S2Polyline_MoveConstruct() {
            var line = MakePolyline("1:1, 4:4");
            S2Polyline moved = new S2Polyline(line);
            Assert.Equal(0, line.NumVertices);
            Assert.Equal(2, moved.NumVertices);
        }

        [Fact]
        public void Test_S2Polyline_MoveAssign() {
            var line = MakePolyline("1:1, 4:4");
            S2Polyline copied;
            copied = line;
            line = new S2Polyline();
            Assert.Equal(0, line.NumVertices);
            Assert.Equal(2, copied.NumVertices);
        }

        [Fact]
        public void Test_S2Polyline_GetLengthAndCentroid() {
            // Construct random great circles and divide them randomly into segments.
            // Then make sure that the length and centroid are correct.  Note that
            // because of the way the centroid is computed, it does not matter how
            // we split the great circle into segments.

            for (int i = 0; i < 100; ++i) {
                // Choose a coordinate frame for the great circle.
                S2Testing.GetRandomFrame(out var x, out var y, out _);

                var vertices = new List<S2Point>();
                for (double theta = 0; theta < S2Constants.M_2_PI;
                     theta += Math.Pow(S2Testing.Random.RandDouble(), 10)) {
                    S2Point p = Math.Cos(theta) * x + Math.Sin(theta) * y;
                    if (!vertices.Any() || p != vertices.Last())
                        vertices.Add(p);
                }
                // Close the circle.
                vertices.Add(vertices[0]);
                S2Polyline line = new S2Polyline(vertices.ToArray());
                S1Angle length = line.GetLength();
                Assert.True(Math.Abs(length.Radians - S2Constants.M_2_PI) <= 2e-14);
                S2Point centroid = line.GetCentroid();
                Assert.True(centroid.Norm <= 2e-14);
            }
        }

        [Fact]
        public void Test_S2Polyline_MayIntersect() {
            S2Point[] vertices = {new S2Point(1, -1.1, 0.8).Normalized,
                              new S2Point(1, -0.8, 1.1).Normalized};
            S2Polyline line = new S2Polyline(vertices);
            for (int face = 0; face < 6; ++face) {
                S2Cell cell = S2Cell.FromFace(face);
                Assert.Equal((face & 1) == 0, line.MayIntersect(cell));
            }
        }

        [Fact]
        public void Test_S2Polyline_Interpolate() {
            S2Point[] vertices = {new S2Point(1, 0, 0),
                              new S2Point(0, 1, 0),
                              new S2Point(0, 1, 1).Normalized,
                              new S2Point(0, 0, 1)};
            S2Polyline line = new S2Polyline(vertices);
            Assert.Equal(vertices[0], line.Interpolate(-0.1));
            Assert.True(S2PointUtil.ApproxEquals(line.Interpolate(0.1),
                new S2Point(1, Math.Tan(0.2 * S2Constants.M_PI_2), 0).Normalized));
            Assert.True(S2PointUtil.ApproxEquals(line.Interpolate(0.25),
                new S2Point(1, 1, 0).Normalized));
            Assert.Equal(vertices[1], line.Interpolate(0.5));
            Assert.True(S2PointUtil.ApproxEquals(vertices[2], line.Interpolate(0.75)));
            Assert.Equal(vertices[0], line.GetSuffix(-0.1, out var next_vertex));
            Assert.Equal(1, next_vertex);
            Assert.True(S2PointUtil.ApproxEquals(vertices[2],
                                         line.GetSuffix(0.75, out next_vertex)));
            Assert.Equal(3, next_vertex);
            Assert.Equal(vertices[3], line.GetSuffix(1.1, out next_vertex));
            Assert.Equal(4, next_vertex);

            // Check the case where the interpolation fraction is so close to 1 that
            // the interpolated point is identical to the last vertex.
            vertices = new[]{
                new S2Point(1, 1, 1).Normalized,
              new S2Point(1, 1, 1 + 1e-15).Normalized,
              new S2Point(1, 1, 1 + 2e-15).Normalized};
            S2Polyline short_line = new S2Polyline(vertices);
            Assert.Equal(vertices[2], short_line.GetSuffix(1.0 - 2e-16, out next_vertex));
            Assert.Equal(3, next_vertex);
        }

        [Fact]
        public void Test_S2Polyline_UnInterpolate() {
            var vertices = new List<S2Point>{ new S2Point(1, 0, 0) };
            S2Polyline point_line = new S2Polyline(vertices.ToArray());
            Assert2.Near(0.0, point_line.UnInterpolate(new S2Point(0, 1, 0), 1));

            vertices.Add(new S2Point(0, 1, 0));
            vertices.Add(new S2Point(0, 1, 1).Normalized);
            vertices.Add(new S2Point(0, 0, 1));
            S2Polyline line = new S2Polyline(vertices.ToArray());

            S2Point interpolated;
            interpolated = line.GetSuffix(-0.1, out var next_vertex);
            Assert2.Near(0.0, line.UnInterpolate(interpolated, next_vertex));
            interpolated = line.GetSuffix(0.0, out next_vertex);
            Assert2.Near(0.0, line.UnInterpolate(interpolated, next_vertex));
            interpolated = line.GetSuffix(0.5, out next_vertex);
            Assert2.Near(0.5, line.UnInterpolate(interpolated, next_vertex));
            interpolated = line.GetSuffix(0.75, out next_vertex);
            Assert2.Near(0.75, line.UnInterpolate(interpolated, next_vertex));
            interpolated = line.GetSuffix(1.1, out next_vertex);
            Assert2.Near(1.0, line.UnInterpolate(interpolated, next_vertex));

            // Check that the return value is clamped to 1.0.
            Assert2.Near(1.0, line.UnInterpolate(new S2Point(0, 1, 0), vertices.Count));
        }

        [Fact]
        public void Test_S2Polyline_Project() {
            S2LatLng[] latlngs = {
      S2LatLng.FromDegrees(0, 0), S2LatLng.FromDegrees(0, 1),
      S2LatLng.FromDegrees(0, 2), S2LatLng.FromDegrees(1, 2)};
            S2Polyline line = new S2Polyline(latlngs);

            Assert.True(S2PointUtil.ApproxEquals(line.Project(
                S2LatLng.FromDegrees(0.5, -0.5).ToPoint(), out var next_vertex),
                S2LatLng.FromDegrees(0, 0).ToPoint()));
            Assert.Equal(1, next_vertex);
            Assert.True(S2PointUtil.ApproxEquals(line.Project(
                S2LatLng.FromDegrees(0.5, 0.5).ToPoint(), out next_vertex),
                S2LatLng.FromDegrees(0, 0.5).ToPoint()));
            Assert.Equal(1, next_vertex);
            Assert.True(S2PointUtil.ApproxEquals(line.Project(
                S2LatLng.FromDegrees(0.5, 1).ToPoint(), out next_vertex),
                S2LatLng.FromDegrees(0, 1).ToPoint()));
            Assert.Equal(2, next_vertex);
            Assert.True(S2PointUtil.ApproxEquals(line.Project(
                S2LatLng.FromDegrees(-0.5, 2.5).ToPoint(), out next_vertex),
                S2LatLng.FromDegrees(0, 2).ToPoint()));
            Assert.Equal(3, next_vertex);
            Assert.True(S2PointUtil.ApproxEquals(line.Project(
                S2LatLng.FromDegrees(2, 2).ToPoint(), out next_vertex),
                S2LatLng.FromDegrees(1, 2).ToPoint()));
            Assert.Equal(4, next_vertex);
        }

        [Fact]
        public void Test_S2Polyline_IsOnRight() {
            S2LatLng[] latlngs = {
      S2LatLng.FromDegrees(0, 0), S2LatLng.FromDegrees(0, 1),
      S2LatLng.FromDegrees(0, 2), S2LatLng.FromDegrees(1, 2)};
            S2Polyline line = new S2Polyline(latlngs);

            Assert.True(line.IsOnRight(S2LatLng.FromDegrees(-0.5, 0.5).ToPoint()));
            Assert.False(line.IsOnRight(S2LatLng.FromDegrees(0.5, -0.5).ToPoint()));
            Assert.False(line.IsOnRight(S2LatLng.FromDegrees(0.5, 0.5).ToPoint()));
            Assert.False(line.IsOnRight(S2LatLng.FromDegrees(0.5, 1).ToPoint()));
            Assert.True(line.IsOnRight(S2LatLng.FromDegrees(-0.5, 2.5).ToPoint()));
            Assert.True(line.IsOnRight(S2LatLng.FromDegrees(1.5, 2.5).ToPoint()));

            // Explicitly test the case where the closest point is an interior vertex.
            latlngs = new[]{ S2LatLng.FromDegrees(0, 0), S2LatLng.FromDegrees(0, 1),
             S2LatLng.FromDegrees(-1, 0)};
            S2Polyline line2 = new S2Polyline(latlngs);

            // The points are chosen such that they are on different sides of the two
            // edges that the interior vertex is on.
            Assert.False(line2.IsOnRight(S2LatLng.FromDegrees(-0.5, 5).ToPoint()));
            Assert.False(line2.IsOnRight(S2LatLng.FromDegrees(5.5, 5).ToPoint()));
        }

        [Fact]
        public void Test_S2Polyline_IntersectsEmptyPolyline() {
            var line1 = MakePolyline("1:1, 4:4");
            S2Polyline empty_polyline = new S2Polyline();
            Assert.False(empty_polyline.Intersects(line1));
        }

        [Fact]
        public void Test_S2Polyline_IntersectsOnePointPolyline() {
            var line1 = MakePolyline("1:1, 4:4");
            var line2 = MakePolyline("1:1");
            Assert.False(line1.Intersects(line2));
        }

        [Fact]
        public void Test_S2Polyline_Intersects() {
            var line1 = MakePolyline("1:1, 4:4");
            var small_crossing = MakePolyline("1:2, 2:1");
            var small_noncrossing = MakePolyline("1:2, 2:3");
            var big_crossing = MakePolyline("1:2, 2:3, 4:3");

            Assert.True(line1.Intersects(small_crossing));
            Assert.False(line1.Intersects(small_noncrossing));
            Assert.True(line1.Intersects(big_crossing));
        }

        [Fact]
        public void Test_S2Polyline_IntersectsAtVertex() {
            var line1 = MakePolyline("1:1, 4:4, 4:6");
            var line2 = MakePolyline("1:1, 1:2");
            var line3 = MakePolyline("5:1, 4:4, 2:2");
            Assert.True(line1.Intersects(line2));
            Assert.True(line1.Intersects(line3));
        }

        [Fact]
        public void Test_S2Polyline_IntersectsVertexOnEdge() {
            var horizontal_left_to_right = MakePolyline("0:1, 0:3");
            var vertical_bottom_to_top = MakePolyline("-1:2, 0:2, 1:2");
            var horizontal_right_to_left = MakePolyline("0:3, 0:1");
            var vertical_top_to_bottom = MakePolyline("1:2, 0:2, -1:2");
            Assert.True(horizontal_left_to_right.Intersects(vertical_bottom_to_top));
            Assert.True(horizontal_left_to_right.Intersects(vertical_top_to_bottom));
            Assert.True(horizontal_right_to_left.Intersects(vertical_bottom_to_top));
            Assert.True(horizontal_right_to_left.Intersects(vertical_top_to_bottom));
        }

        [Fact]
        public void Test_S2Polyline_SpaceUsedEmptyPolyline() {
            var line = MakePolyline("");
            Assert.True(line.SpaceUsed() > 0);
        }

        [Fact]
        public void Test_S2Polyline_SpaceUsedNonEmptyPolyline() {
            var line = MakePolyline("1:1, 4:4, 4:6");
            Assert.True(line.SpaceUsed() > 3 * Marshal.SizeOf(typeof(S2Point)));
        }

        [Fact]
        public void Test_S2Polyline_SubsampleVerticesTrivialInputs() {
            // No vertices.
            CheckSubsample("", 1.0, "");
            // One vertex.
            CheckSubsample("0:1", 1.0, "0");
            // Two vertices.
            CheckSubsample("10:10, 11:11", 5.0, "0,1");
            // Three points on a straight line.
            // In theory, zero tolerance should work, but in practice there are floating
            // point errors.
            CheckSubsample("-1:0, 0:0, 1:0", 1e-15, "0,2");
            // Zero tolerance on a non-straight line.
            CheckSubsample("-1:0, 0:0, 1:1", 0.0, "0,1,2");
            // Negative tolerance should return all vertices.
            CheckSubsample("-1:0, 0:0, 1:1", -1.0, "0,1,2");
            // Non-zero tolerance with a straight line.
            CheckSubsample("0:1, 0:2, 0:3, 0:4, 0:5", 1.0, "0,4");

            // And finally, verify that we still do something reasonable if the client
            // passes in an invalid polyline with two or more adjacent vertices.
            CheckSubsample("0:1, 0:1, 0:1, 0:2", 0.0, "0,3", false);
        }

        [Fact]
        public void Test_S2Polyline_SubsampleVerticesSimpleExample() {
            var poly_str = "0:0, 0:1, -1:2, 0:3, 0:4, 1:4, 2:4.5, 3:4, 3.5:4, 4:4";
            CheckSubsample(poly_str, 3.0, "0,9");
            CheckSubsample(poly_str, 2.0, "0,6,9");
            CheckSubsample(poly_str, 0.9, "0,2,6,9");
            CheckSubsample(poly_str, 0.4, "0,1,2,3,4,6,9");
            CheckSubsample(poly_str, 0, "0,1,2,3,4,5,6,7,8,9");
        }

        [Fact]
        public void Test_S2Polyline_SubsampleVerticesGuarantees() {
            // Check that duplicate vertices are never generated.
            CheckSubsample("10:10, 12:12, 10:10", 5.0, "0");
            CheckSubsample("0:0, 1:1, 0:0, 0:120, 0:130", 5.0, "0,3,4");

            // Check that points are not collapsed if they would create a line segment
            // longer than 90 degrees, and also that the code handles original polyline
            // segments longer than 90 degrees.
            CheckSubsample("90:0, 50:180, 20:180, -20:180, -50:180, -90:0, 30:0, 90:0",
                           5.0, "0,2,4,5,6,7");

            // Check that the output polyline is parametrically equivalent and not just
            // geometrically equivalent, i.e. that backtracking is preserved.  The
            // algorithm achieves this by requiring that the points must be encountered
            // in increasing order of distance along each output segment, except for
            // points that are within "tolerance" of the first vertex of each segment.
            CheckSubsample("10:10, 10:20, 10:30, 10:15, 10:40", 5.0, "0,2,3,4");
            CheckSubsample("10:10, 10:20, 10:30, 10:10, 10:30, 10:40", 5.0, "0,2,3,5");
            CheckSubsample("10:10, 12:12, 9:9, 10:20, 10:30", 5.0, "0,4");
        }

        [Fact]
        public void Test_S2Polyline_ApproxEquals() {
            S1Angle degree = S1Angle.FromDegrees(1);

            // Close lines, differences within max_error.
            Assert.True(TestEquals("0:0, 0:10, 5:5",
                                   "0:0.1, -0.1:9.9, 5:5.2",
                                   0.5 * degree));

            // Close lines, differences outside max_error.
            Assert.False(TestEquals("0:0, 0:10, 5:5",
                                    "0:0.1, -0.1:9.9, 5:5.2",
                                    0.01 * degree));

            // Same line, but different number of vertices.
            Assert.False(TestEquals("0:0, 0:10, 0:20", "0:0, 0:20", 0.1 * degree));

            // Same vertices, in different order.
            Assert.False(TestEquals("0:0, 5:5, 0:10", "5:5, 0:10, 0:0", 0.1 * degree));
        }

        [Fact]
        public void Test_S2Polyline_EncodeDecode() {
            var polyline = MakePolyline("0:0, 0:10, 10:20, 20:30");
            Encoder encoder = new Encoder();
            polyline.Encode(encoder);
            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            var (success, decoded_polyline) = S2Polyline.DecodeStatic(decoder);
            Assert.True(success);
            Assert.True(decoded_polyline.ApproxEquals(polyline, S1Angle.Zero));
        }

        [Fact]
        public void Test_S2PolylineShape_Basic() {
            var polyline = MakePolyline("0:0, 1:0, 1:1, 2:1");
            var shape = new S2Polyline.Shape(polyline);
            Assert.Equal(polyline, shape.Polyline);
            Assert.Equal(3, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(3, shape.GetChain(0).Length);
            var edge2 = shape.GetEdge(2);
            Assert.Equal(S2LatLng.FromDegrees(1, 1).ToPoint(), edge2.V0);
            Assert.Equal(S2LatLng.FromDegrees(2, 1).ToPoint(), edge2.V1);
            Assert.Equal(1, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2PolylineShape_EmptyPolyline() {
            S2Polyline polyline = new S2Polyline();
            S2Polyline.Shape shape = new S2Polyline.Shape(polyline);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2PolylineOwningShape_Ownership() {
            // Debug mode builds will catch any memory leak below.
            var polyline = new S2Polyline(Array.Empty<S2Point>());
            _ = new S2Polyline.OwningShape(polyline);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_PolylineOverlapsSelf() {
            string pline = "1:1, 2:2, -1:10";
            TestNearlyCovers(pline, pline, 1e-10, true, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_PolylineDoesNotOverlapReverse() {
            TestNearlyCovers("1:1, 2:2, -1:10", "-1:10, 2:2, 1:1", 1e-10, false, false);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_PolylineOverlapsEquivalent() {
            // These two polylines trace the exact same polyline, but the second one uses
            // three points instead of two.
            TestNearlyCovers("1:1, 2:1", "1:1, 1.5:1, 2:1", 1e-10, true, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_ShortCoveredByLong() {
            // The second polyline is always within 0.001 degrees of the first polyline,
            // but the first polyline is too long to be covered by the second.
            TestNearlyCovers(
                "-5:1, 10:1, 10:5, 5:10", "9:1, 9.9995:1, 10.0005:5", 1e-3, false, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_PartialOverlapOnly() {
            // These two polylines partially overlap each other, but neither fully
            // overlaps the other.
            TestNearlyCovers("-5:1, 10:1", "0:1, 20:1", 1.0, false, false);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_ShortBacktracking() {
            // Two lines that backtrack a bit (less than 1.5 degrees) on different edges.
            // A simple greedy matching algorithm would fail on this example.
            string t1 = "0:0, 0:2, 0:1, 0:4, 0:5";
            string t2 = "0:0, 0:2, 0:4, 0:3, 0:5";
            TestNearlyCovers(t1, t2, 1.5, true, true);
            TestNearlyCovers(t1, t2, 0.5, false, false);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_LongBacktracking() {
            // Two arcs with opposite direction do not overlap if the shorter arc is
            // longer than max_error, but do if the shorter arc is shorter than max-error.
            TestNearlyCovers("5:1, -5:1", "1:1, 3:1", 1.0, false, false);
            TestNearlyCovers("5:1, -5:1", "1:1, 3:1", 2.5, false, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_IsResilientToDuplicatePoints() {
            // S2Polyines are not generally supposed to contain adjacent, identical
            // points, but it happens in practice.  We also set S2Debug.DISABLE so
            // debug binaries won't abort on such polylines.
            TestNearlyCovers("0:1, 0:2, 0:2, 0:3", "0:1, 0:1, 0:1, 0:3",
                             1e-10, true, true, false);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_CanChooseBetweenTwoPotentialStartingPoints() {
            // Can handle two possible starting points, only one of which leads to finding
            // a correct path.  In the first polyline, the edge from 0:1.1 to 0:0 and the
            // edge from 0:0.9 to 0:2 might be lucrative starting states for covering the
            // second polyline, because both edges are with the max_error of 1.5 degrees
            // from 0:10.  However, only the latter is actually effective.
            TestNearlyCovers("0:11, 0:0, 0:9, 0:20", "0:10, 0:15", 1.5, false, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_StraightAndWigglyPolylinesCoverEachOther() {
            TestNearlyCovers("40:1, 20:1",
                "39.9:0.9, 40:1.1, 30:1.15, 29:0.95, 28:1.1, 27:1.15, " +
                "26:1.05, 25:0.85, 24:1.1, 23:0.9, 20:0.99",
                0.2, true, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_MatchStartsAtLastVertex() {
            // The first polyline covers the second, but the matching segment starts at
            // the last vertex of the first polyline.
            TestNearlyCovers(
                "0:0, 0:2", "0:2, 0:3", 1.5, false, true);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_MatchStartsAtDuplicatedLastVertex() {
            TestNearlyCovers(
                "0:0, 0:2, 0:2, 0:2", "0:2, 0:3", 1.5, false, true, false);
        }

        [Fact]
        public void Test_S2PolylineCoveringTest_EmptyPolylines() {
            // We expect:
            //    anything.covers(empty) = true
            //    empty.covers(nonempty) = false
            TestNearlyCovers("0:1, 0:2", "", 0.0, false, true);
            TestNearlyCovers("", "", 0.0, true, true);
        }

        // Wraps S2TextFormat.MakePolyline in order to test Encode/Decode.
        private static S2Polyline MakePolyline(string str, bool assertIsValid = true)
        {
            var polyline = S2TextFormat.MakePolylineOrDie(str, assertIsValid);
            Encoder encoder = new Encoder();
            polyline.Encode(encoder);
            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            var (_, decoded_polyline) = S2Polyline.DecodeStatic(decoder);
            return decoded_polyline;
        }

        private static string JoinInts(int[] ints)
        {
            var result = new StringBuilder();
            int n = ints.Length;
            for (int i = 0; i + 1 < n; ++i)
            {
                result.Append(ints[i] + ",");
            }
            if (n > 0)
            {
                result.Append(ints[n - 1]);
            }
            return result.ToString();
        }

        private void CheckSubsample(string polyline_str, double tolerance_degrees, string expected_str, bool assertIsValid = true)
        {
            _logger.WriteLine($@"""{polyline_str}"", tolerance {tolerance_degrees}");
            var polyline = MakePolyline(polyline_str, assertIsValid);
            polyline.SubsampleVertices(S1Angle.FromDegrees(tolerance_degrees), out var indices);
            Assert.Equal(expected_str, JoinInts(indices));
        }

        private static bool TestEquals(string a_str, string b_str, S1Angle max_error)
        {
            var a = MakePolyline(a_str);
            var b = MakePolyline(b_str);
            return a.ApproxEquals(b, max_error);
        }

        private void TestNearlyCovers(string a_str, string b_str, double max_error_degrees, bool expect_b_covers_a, bool expect_a_covers_b, bool assertIsValid = true)
        {
            _logger.WriteLine($@"a=""{a_str}"", b=""{b_str}"", max error={max_error_degrees}");
            var a = MakePolyline(a_str, assertIsValid);
            var b = MakePolyline(b_str, assertIsValid);
            S1Angle max_error = S1Angle.FromDegrees(max_error_degrees);
            Assert.Equal(expect_b_covers_a, b.NearlyCovers(a, max_error));
            Assert.Equal(expect_a_covers_b, a.NearlyCovers(b, max_error));
        }
    }
}
