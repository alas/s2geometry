using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace S2Geometry
{
    public class S2TextFormatTests
    {
        private const int kIters = 10000;

        [Fact]
        public void Test_Degrees()
        {
            ExpectString("0:0", S2LatLng.FromDegrees(0, 0));
            ExpectString("1e-20:1e-30", S2LatLng.FromDegrees(1e-20, 1e-30));
        }

        [Fact]
        public void Test_E5()
        {
            for (var i = 0; i < kIters; i++)
            {
                var ll = S2LatLng.FromPoint(S2Testing.RandomPoint());
                var ll_e5 = S2LatLng.FromE5(ll.Lat.E5, ll.Lng.E5);
                ExpectMaxDigits(ll_e5, 5);
            }
        }

        [Fact]
        public void Test_E6()
        {
            for (var i = 0; i < kIters; i++)
            {
                var ll = S2LatLng.FromPoint(S2Testing.RandomPoint());
                var ll_e6 = S2LatLng.FromE6(ll.Lat.E6, ll.Lng.E6);
                ExpectMaxDigits(ll_e6, 6);
            }
        }

        [Fact]
        public void Test_E7()
        {
            ExpectMaxDigits(S2LatLng.FromDegrees(0, 0), 7);
            for (var i = 0; i < kIters; i++)
            {
                var ll = S2LatLng.FromPoint(S2Testing.RandomPoint());
                var ll_e7 = S2LatLng.FromE7(ll.Lat.E7, ll.Lng.E7);
                ExpectMaxDigits(ll_e7, 7);
            }
        }

        [Fact]
        public void Test_MinimalDigits_DoubleConstants()
        {
            // Verify that points specified as floating-point literals in degrees using
            // up to 10 digits after the decimal point are formatted with the minimal
            // number of digits.
            for (var i = 0; i < kIters; i++)
            {
                var max_digits = S2Testing.Random.Uniform(11);
                Int64 scale = (long)Math.Round(Math.Pow(10, max_digits));
                Int64 lat = (long)Math.Round(S2Testing.Random.UniformDouble(-90 * scale, 90 * scale));
                Int64 lng = (long)Math.Round(S2Testing.Random.UniformDouble(-180 * scale, 180 * scale));
                var ll = S2LatLng.FromDegrees(lat / (double)scale, lng / (double)scale);
                ExpectMaxDigits(ll, max_digits);
            }
        }

        [Fact]
        public void Test_UninitializedLoop()
        {
            var loop = new S2Loop(Array.Empty<S2Point>());
            Assert.Equal("", loop.ToDebugString());
        }

        [Fact]
        public void Test_EmptyLoop()
        {
            var empty = S2Loop.kEmpty;
            Assert.Equal("empty", empty.ToDebugString());
        }

        [Fact]
        public void Test_FullLoop()
        {
            var full = S2Loop.kFull;
            Assert.Equal("full", full.ToDebugString());
        }

        [Fact]
        public void Test_FullLoopSpan()
        {
            var points = Array.Empty<S2Point>();
            Assert.Equal("full", points.ToDebugStringLoop());
        }

        [Fact]
        public void Test_EmptyPolyline()
        {
            var polyline = new S2Polyline();
            Assert.Equal("", polyline.ToDebugString());
        }

        [Fact]
        public void Test_EmptyPointVector()
        {
            var points = Array.Empty<S2Point>();
            Assert.Equal("", points.ToDebugString());
        }

        [Fact]
        public void Test_EmptyPolygon()
        {
            var empty = new S2Polygon();
            Assert.Equal("empty", empty.ToDebugString());
        }

        [Fact]
        public void Test_FullPolygon()
        {
            var full = new S2Polygon(S2Loop.kFull);
            Assert.Equal("full", full.ToDebugString());
        }

        [Fact]
        public void Test_S2PolygonLoopSeparator()
        {
            // Shells and holes same direction.
            var loops = new[]{ "0:0, 0:5, 5:0", "1:1, 1:4, 4:1" };
            var polygon = S2TextFormat.MakePolygonOrDie(string.Join("; ", loops));
            Assert.Equal(string.Join(";\n", loops), polygon.ToDebugString());
            Assert.Equal(string.Join("; ", loops), polygon.ToDebugString("; "));
        }

        [Fact]
        public void Test_LaxPolygonLoopSeparator()
        {
            string kLoop1 = "0:0, 0:5, 5:0";
            string kLoop2 = "1:1, 4:1, 1:4";  // Interior on left of all loops.
            var polygon = S2TextFormat.MakeLaxPolygonOrDie(kLoop1 + "; " + kLoop2);
            Assert.Equal(kLoop1 + ";\n" + kLoop2, polygon.ToDebugString());
            Assert.Equal(kLoop1 + "; " + kLoop2, polygon.ToDebugString("; "));
        }

        [Fact]
        public void Test_MakeLaxPolygon_Empty()
        {
            // Verify that "" and "empty" both create empty polygons.
            var shape = S2TextFormat.MakeLaxPolygonOrDie("");
            Assert.Equal(0, shape.NumLoops);
            shape = S2TextFormat.MakeLaxPolygonOrDie("empty");
            Assert.Equal(0, shape.NumLoops);
        }

        [Fact]
        public void Test_MakeLaxPolygon_Full()
        {
            var shape = S2TextFormat.MakeLaxPolygonOrDie("full");
            Assert.Equal(1, shape.NumLoops);
            Assert.Equal(0, shape.NumLoopVertices(0));
        }

        [Fact]
        public void Test_MakeLaxPolygon_FullWithHole()
        {
            var shape = S2TextFormat.MakeLaxPolygonOrDie("full; 0:0");
            Assert.Equal(2, shape.NumLoops);
            Assert.Equal(0, shape.NumLoopVertices(0));
            Assert.Equal(1, shape.NumLoopVertices(1));
            Assert.Equal(1, shape.NumEdges);
        }

        [Fact]
        public void Test_S2ShapeIndex()
        {
            TestS2ShapeIndex("# #");
            TestS2ShapeIndex("0:0 # #");
            TestS2ShapeIndex("0:0 | 1:0 # #");
            TestS2ShapeIndex("0:0 | 1:0 # #");
            TestS2ShapeIndex("# 0:0, 0:0 #");
            TestS2ShapeIndex("# 0:0, 0:0 | 1:0, 2:0 #");
            TestS2ShapeIndex("# # 0:0");
            TestS2ShapeIndex("# # 0:0, 0:1");
            TestS2ShapeIndex("# # 0:0, 0:1, 1:0");
            TestS2ShapeIndex("# # 0:0, 0:1, 1:0; 2:2");
            TestS2ShapeIndex("# # full");
        }

        [Fact]
        public void Test_MakePoint_ValidInput()
        {
            Assert.True(S2TextFormat.MakePoint("-20:150", out var point));
            Assert.Equal(S2LatLng.FromDegrees(-20, 150).ToPoint(), point);
        }

        [Fact]
        public void Test_MakePoint_InvalidInput()
        {
            Assert.False(S2TextFormat.MakePoint("blah", out _));
        }

        [Fact]
        public void Test_ParseLatLngs_ValidInput()
        {
            var latlngs = new List<S2LatLng>();
            Assert.True(S2TextFormat.ParseLatLngs("-20:150, -20:151, -19:150", latlngs));
            Assert.Equal(3, latlngs.Count);
            Assert.Equal(latlngs[0], S2LatLng.FromDegrees(-20, 150));
            Assert.Equal(latlngs[1], S2LatLng.FromDegrees(-20, 151));
            Assert.Equal(latlngs[2], S2LatLng.FromDegrees(-19, 150));
        }

        [Fact]
        public void Test_ParseLatLngs_InvalidInput()
        {
            var latlngs = new List<S2LatLng>();
            Assert.False(S2TextFormat.ParseLatLngs("blah", latlngs));
        }

        [Fact]
        public void Test_ParsePoints_ValidInput()
        {
            var vertices = new List<S2Point>();
            Assert.True(S2TextFormat.ParsePoints("-20:150, -20:151, -19:150", vertices));
            Assert.Equal(3, vertices.Count);
            Assert.Equal(vertices[0], S2LatLng.FromDegrees(-20, 150).ToPoint());
            Assert.Equal(vertices[1], S2LatLng.FromDegrees(-20, 151).ToPoint());
            Assert.Equal(vertices[2], S2LatLng.FromDegrees(-19, 150).ToPoint());
        }

        [Fact]
        public void Test_ParsePoints_InvalidInput()
        {
            var vertices = new List<S2Point>();
            Assert.False(S2TextFormat.ParsePoints("blah",  vertices));
        }

        [Fact]
        public void Test_MakeLatLngRect_ValidInput()
        {
            Assert.True(S2TextFormat.MakeLatLngRect("-10:-10, 10:10", out var rect));
            Assert.Equal(rect, new S2LatLngRect(
                S2LatLng.FromDegrees(-10, -10),
                S2LatLng.FromDegrees(10, 10)));
        }

        [Fact]
        public void Test_MakeLatLngRect_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeLatLngRect("blah", out _));
        }

        [Fact]
        public void Test_MakeLatLng_ValidInput()
        {
            Assert.True(S2TextFormat.MakeLatLng("-12.3:45.6", out var latlng));
            Assert.Equal(latlng, S2LatLng.FromDegrees(-12.3, 45.6));
        }

        [Fact]
        public void Test_MakeLatLng_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeLatLng("blah", out _));
        }

        [Fact]
        public void Test_MakeCellId_ValidInput()
        {
            Assert.True(S2TextFormat.MakeCellId("3/", out var cellId));
            Assert.Equal(cellId, S2CellId.FromFace(3));
        }

        [Fact]
        public void Test_MakeCellId_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeCellId("blah", out _));
            Assert.False(S2TextFormat.MakeCellId("6/0", out _));
            Assert.False(S2TextFormat.MakeCellId("3/04", out _));
        }

        [Fact]
        public void Test_MakeCellUnion_ValidInput()
        {
            Assert.True(S2TextFormat.MakeCellUnion("1/3, 4/", out var cellUnion));
            var expected = new S2CellUnion(new List<S2CellId> {
                S2CellId.FromFace(1).Child(3), S2CellId.FromFace(4) });
            Assert.Equal(cellUnion, expected);
        }

        [Fact]
        public void Test_MakeCellUnion_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeCellUnion("abc", out _));
            Assert.False(S2TextFormat.MakeCellUnion("3/1 4/1", out _));
        }

        [Fact]
        public void Test_MakeLoop_ValidInput()
        {
            Assert.True(S2TextFormat.MakeLoop("-20:150, -20:151, -19:150", out var loop));
            var expected = new S2Loop(new[]
            {
                S2LatLng.FromDegrees(-20, 150).ToPoint(),
                S2LatLng.FromDegrees(-20, 151).ToPoint(),
                S2LatLng.FromDegrees(-19, 150).ToPoint(),
            });
            Assert.True(loop.BoundaryApproxEquals(expected));
        }

        [Fact]
        public void Test_MakeLoop_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeLoop("blah", out _));
        }

        [Fact]
        public void Test_MakePolyline_ValidInput()
        {
            Assert.True(S2TextFormat.MakePolyline("-20:150, -20:151, -19:150", out var polyline));
            var expected = new S2Polyline(new[] {
                S2LatLng.FromDegrees(-20, 150).ToPoint(),
                S2LatLng.FromDegrees(-20, 151).ToPoint(),
                S2LatLng.FromDegrees(-19, 150).ToPoint(),
            });
            Assert.True(polyline == expected);
        }

        [Fact]
        public void Test_MakePolyline_InvalidInput()
        {
            Assert.False(S2TextFormat.MakePolyline("blah", out _));
        }

        [Fact]
        public void Test_MakeLaxPolyline_ValidInput()
        {
            Assert.True(S2TextFormat.MakeLaxPolyline("-20:150, -20:151, -19:150", out var laxPolyline));

            // No easy equality check for LaxPolylines; check vertices instead.
            Assert.Equal(3, laxPolyline.NumVertices);
            Assert.True(new S2LatLng(laxPolyline.Vertex(0)).ApproxEquals(S2LatLng.FromDegrees(-20, 150)));
            Assert.True(new S2LatLng(laxPolyline.Vertex(1)).ApproxEquals(S2LatLng.FromDegrees(-20, 151)));
            Assert.True(new S2LatLng(laxPolyline.Vertex(2)).ApproxEquals(S2LatLng.FromDegrees(-19, 150)));
        }

        [Fact]
        public void Test_MakeLaxPolyline_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeLaxPolyline("blah", out _));
        }

        [Fact]
        public void Test_MakePolygon_ValidInput()
        {
            Assert.True(S2TextFormat.MakePolygon("-20:150, -20:151, -19:150", out var polygon));
            var vertices = new[]
            {
                S2LatLng.FromDegrees(-20, 150).ToPoint(),
                S2LatLng.FromDegrees(-20, 151).ToPoint(),
                S2LatLng.FromDegrees(-19, 150).ToPoint(),
            };
            var expected = new S2Polygon(new S2Loop(vertices));
            Assert.True(polygon == expected);
        }

        [Fact]
        public void Test_MakePolygon_InvalidInput()
        {
            Assert.False(S2TextFormat.MakePolygon("blah", out _));
        }

        [Fact]
        public void Test_MakePolygon_Empty()
        {
            // Verify that "" and "empty" both create empty polygons.
            Assert.True(S2TextFormat.MakePolygon("", out var polygon));
            Assert.True(polygon.IsEmpty);
            Assert.True(S2TextFormat.MakePolygon("empty", out polygon));
            Assert.True(polygon.IsEmpty);
        }

        [Fact]
        public void Test_MakePolygon_Full()
        {
            // Verify that "full" creates the full polygon.
            Assert.True(S2TextFormat.MakePolygon("full", out var polygon));
            Assert.True(polygon.IsFull);
        }

        [Fact]
        public void Test_MakeVerbatimPolygon_ValidInput()
        {
            Assert.True(S2TextFormat.MakeVerbatimPolygon("-20:150, -20:151, -19:150", out var polygon));
            var vertices = new[]
            {
                new [] {-20, 150 },
                new [] {-20, 151 },
                new [] {-19, 150 },
            }.Select(t => S2LatLng.FromDegrees(t[0], t[1]).ToPoint());
            var expected = new S2Polygon(new S2Loop(vertices));
            Assert.Equal(polygon, expected);
        }

        [Fact]
        public void Test_MakeVerbatimPolygon_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeVerbatimPolygon("blah", out _));
        }

        [Fact]
        public void Test_MakeLaxPolygon_ValidInput()
        {
            Assert.True(S2TextFormat.MakeLaxPolygon("-20:150, -20:151, -19:150", out var lax_polygon));
            // No easy equality check for LaxPolygons; check vertices instead.
            Assert.Equal(1, lax_polygon.NumLoops);
            Assert.Equal(3, lax_polygon.NumVertices);
            Assert.True(new S2LatLng(lax_polygon.LoopVertex(0, 0).First()).ApproxEquals(S2LatLng.FromDegrees(-20, 150)));
            Assert.True(new S2LatLng(lax_polygon.LoopVertex(0, 1).First()).ApproxEquals(S2LatLng.FromDegrees(-20, 151)));
            Assert.True(new S2LatLng(lax_polygon.LoopVertex(0, 2).First()).ApproxEquals(S2LatLng.FromDegrees(-19, 150)));
        }

        [Fact]
        public void Test_MakeLaxPolygon_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeLaxPolygon("blah", out _));
        }

        [Fact]
        public void Test_MakeIndex_ValidInput()
        {
            const string valid = "# 0:0, 0:0 | 1:0, 2:0 #";
            Assert.True(S2TextFormat.MakeIndex(valid, out var index));
            Assert.Equal(valid, index.ToDebugString());
        }

        [Fact]
        public void Test_MakeIndex_InvalidInput()
        {
            Assert.False(S2TextFormat.MakeIndex("# blah #", out _));
        }

        /// <summary>
        /// Verify that S2TextFormat.ToDebugString() formats the given lat/lng with at most
        /// "max_digits" after the decimal point and has no trailing zeros.
        /// </summary>
        private static void ExpectMaxDigits(S2LatLng ll, int maxDigits)
        {
            var result = ll.ToPoint().ToDebugString();
            var values = result.Split(':', StringSplitOptions.RemoveEmptyEntries);
            Assert.Equal(2, values.Length);
            foreach (var value in values)
            {
                var numDigits = 0;
                if (value.Contains('.'))
                {
                    numDigits = value.Split(".")[1].Length;
                    Assert.NotEqual('0', value.Last());
                }
                Assert.True(numDigits <= maxDigits);
            }
        }

        private static void ExpectString(string expected, S2LatLng ll)
        {
            Assert.Equal(expected, ll.ToPoint().ToDebugString());
        }

        private static void TestS2ShapeIndex(string str)
        {
            Assert.Equal(str, S2TextFormat.MakeIndexOrDie(str).ToDebugString());
        }
    }
}

