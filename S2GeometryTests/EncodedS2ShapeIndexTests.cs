using System.Collections.Generic;
using System;
using Xunit;
using Xunit.Abstractions;
using S2Geometry.S2BuilderUtil;
using Options = S2Geometry.S2BuilderUtil.S2PolylineLayer.Options;

namespace S2Geometry
{
    public class EncodedS2ShapeIndexTests
    {
        private readonly ITestOutputHelper _logger;

        public EncodedS2ShapeIndexTests(ITestOutputHelper logger) { _logger = logger; }

        [Fact]
        public void Test_EncodedS2ShapeIndex_Empty()
        {
            MutableS2ShapeIndex index = new();
            TestEncodedS2ShapeIndex<S2LaxPolylineShape>(index, 4); // S2LaxPolylineShape
        }

        [Fact]
        public void Test_EncodedS2ShapeIndex_OneEdge()
        {
            MutableS2ShapeIndex index = new();
            index.Add(S2TextFormat.MakeLaxPolylineOrDie("1:1, 2:2"));
            TestEncodedS2ShapeIndex<S2LaxPolylineShape>(index, 8); // S2LaxPolylineShape
        }

        [Fact]
        public void Test_EncodedS2ShapeIndex_RegularLoops()
        {
            (int num_edges, int expected_bytes)[] test_cases = {
    (4, 8),
    (8, 8),
    ( 16, 16),
    ( 64, 77),
    ( 256, 327),
    ( 4096, 8813),
    ( 65536, 168291),
  };
            foreach (var (num_edges, expected_bytes) in test_cases)
            {
                MutableS2ShapeIndex index = new();
                S2Testing.Random.Reset(num_edges);
                _logger.WriteLine($"num_edges = {num_edges}");
                S2Polygon polygon = new(S2Loop.MakeRegularLoop(new S2Point(3, 2, 1).Normalized,
                                                          S1Angle.FromDegrees(0.1),
                                                          num_edges));
                index.Add(new S2LaxPolygonShape(polygon));
                TestEncodedS2ShapeIndex<S2LaxPolygonShape>( // EncodedS2LaxPolygonShape
                    index, expected_bytes);
            }
        }

        [Fact]
        public void Test_EncodedS2ShapeIndex_OverlappingPointClouds()
        {
            (int num_shapes, int num_points_per_shape, int expected_bytes)[] test_cases = {
    (1, 50, 83),
    (2, 100, 583),
    (4, 100, 1383),
  };
            S2Cap cap = new(new S2Point(0.1, -0.4, 0.3).Normalized, S1Angle.FromDegrees(1));
            foreach (var (num_shapes, num_points_per_shape, expected_bytes) in test_cases)
            {
                MutableS2ShapeIndex index = new();
                S2Testing.Random.Reset(num_shapes);
                _logger.WriteLine($"num_shapes = {num_shapes}");
                for (int i = 0; i < num_shapes; ++i)
                {
                    List<S2Point> points = new();
                    for (int j = 0; j < num_points_per_shape; ++j)
                    {
                        points.Add(S2Testing.SamplePoint(cap));
                    }
                    index.Add(new S2PointVectorShape(points.ToArray()));
                }
                TestEncodedS2ShapeIndex<S2PointVectorShape>( // EncodedS2PointVectorShape
                    index, expected_bytes);
            }
        }

        [Fact]
        public void Test_EncodedS2ShapeIndex_OverlappingPolylines()
        {
            (int num_shapes, int num_shape_edges, int expected_bytes)[] test_cases = {
    (2, 50, 139),
    (10, 50, 777),
    (20, 50, 2219),
  };
            S2Cap cap = new(new S2Point(-0.2, -0.3, 0.4).Normalized, S1Angle.FromDegrees(0.1));
            foreach (var (num_shapes, num_shape_edges, expected_bytes) in test_cases)
            {
                S1Angle edge_len = 2 * cap.RadiusAngle / num_shape_edges;
                MutableS2ShapeIndex index = new();
                S2Testing.Random.Reset(num_shapes);
                _logger.WriteLine($"num_shapes = {num_shapes}");
                for (int i = 0; i < num_shapes; ++i)
                {
                    S2Point a = S2Testing.SamplePoint(cap), b = S2Testing.RandomPoint();
                    List<S2Point> vertices = new();
                    int n = num_shape_edges;
                    for (int j = 0; j <= n; ++j)
                    {
                        vertices.Add(S2EdgeDistances.InterpolateAtDistance(j * edge_len, a, b));
                    }
                    index.Add(new S2LaxPolylineShape(vertices.ToArray()));
                }
                TestEncodedS2ShapeIndex<S2LaxPolylineShape>( // S2LaxPolylineShape, EncodedS2LaxPolylineShape
                    index, expected_bytes);
            }
        }

        [Fact]
        public void Test_EncodedS2ShapeIndex_OverlappingLoops()
        {
            (int num_shapes, int max_edges_per_loop, int expected_bytes)[] test_cases = {
    (2, 250, 138),
    (5, 250, 1084),
    (25, 50, 3673),
  };
            S2Cap cap = new(new S2Point(-0.1, 0.25, 0.2).Normalized, S1Angle.FromDegrees(3));
            foreach (var (num_shapes, max_edges_per_loop, expected_bytes) in test_cases)
            {
                MutableS2ShapeIndex index = new();
                S2Testing.Random.Reset(num_shapes);
                _logger.WriteLine($"num_shapes = {num_shapes}");
                for (int i = 0; i < num_shapes; ++i)
                {
                    S2Point center = S2Testing.SamplePoint(cap);
                    double radius_fraction = S2Testing.Random.RandDouble();
                    // Scale the number of edges so that they are all about the same length
                    // (similar to modeling all geometry at a similar resolution).
                    int num_edges = Convert.ToInt32(Math.Max(3.0, max_edges_per_loop * radius_fraction));
                    S2Polygon polygon = new(S2Loop.MakeRegularLoop(
                        center, cap.RadiusAngle * radius_fraction, num_edges));
                    index.Add(new S2LaxPolygonShape(polygon));
                }
                TestEncodedS2ShapeIndex<S2LaxPolygonShape>( // S2LaxPolygonShape, EncodedS2LaxPolygonShape
                    index, expected_bytes);
            }
        }

        [Fact]
        public void Test_EncodedS2ShapeIndex_SnappedFractalPolylines()
        {
            MutableS2ShapeIndex index = new();
            S2Builder builder = new(new S2Builder.Options(new S2CellIdSnapFunction()));
            for (int i = 0; i < 5; ++i)
            {
                builder.StartLayer(new IndexedLaxPolylineLayer(index));
                S2Testing.Fractal fractal = new();
                fractal.SetLevelForApproxMaxEdges(3 * 256);
                var frame = S2PointUtil.GetFrame(S2LatLng.FromDegrees(10, i).ToPoint());
                var loop = fractal.MakeLoop(frame, S1Angle.FromDegrees(0.1));
                List<S2Point> vertices = new();
                S2Testing.AppendLoopVertices(loop, vertices);
                S2Polyline polyline = new(vertices.ToArray());
                builder.AddPolyline(polyline);
            }
            Assert.True(builder.Build(out _));
            TestEncodedS2ShapeIndex<S2LaxPolylineShape>( // EncodedS2LaxPolylineShape
                index, 8698);
        }

        private static bool DecodeHomegeneousShapeIndex<Shape>(EncodedS2ShapeIndex index, Decoder decoder)
            where Shape : S2Shape, IEncodeInit, new() =>
            index.Init(decoder, new S2ShapeUtilCoding.HomogeneousShapeFactory<Shape>(decoder));

        private static void TestEncodedS2ShapeIndex<Shape>(MutableS2ShapeIndex expected, int expected_bytes)
             where Shape : S2Shape, IEncodeInit, new()
        {
            Encoder encoder = new();
            S2ShapeUtilCoding.EncodeHomogeneousShapes(expected, encoder);
            int shapes_bytes = encoder.Length;
            expected.Encode(encoder);
            Assert.Equal(expected_bytes, encoder.Length - shapes_bytes);
            Decoder decoder = new(encoder.Buffer, 0, encoder.Length);
            EncodedS2ShapeIndex actual = new();
            Assert.True(DecodeHomegeneousShapeIndex<Shape>(actual, decoder));
            Assert.Equal(expected.Options_.MaxEdgesPerCell,
                      actual.Options_.MaxEdgesPerCell);
            S2ShapeTestsUtil.ExpectEqual(expected, actual);
        }

        // Like S2PolylineLayer, but converts the polyline to an S2LaxPolylineShape
        // and adds it to an S2ShapeIndex (if the polyline is non-empty).
        private class IndexedLaxPolylineLayer : S2Builder.Layer
        {
            private readonly MutableS2ShapeIndex index_;
            private readonly S2Polyline polyline_;
            private readonly S2PolylineLayer layer_;

            public IndexedLaxPolylineLayer(MutableS2ShapeIndex index, Options options = null)
            {
                index_ = index;
                polyline_ = new S2Polyline();
                layer_ = new(polyline_, options ?? new Options());
            }

            public override S2Builder.GraphOptions GraphOptions_()
            {
                return layer_.GraphOptions_();
            }

            public override void Build(S2Builder.Graph g, out S2Error error)
            {
                layer_.Build(g, out error);
                if (error.IsOk && polyline_.NumVertices > 0)
                {
                    index_.Add(new S2LaxPolylineShape(polyline_));
                }
            }
        }
    }
}
