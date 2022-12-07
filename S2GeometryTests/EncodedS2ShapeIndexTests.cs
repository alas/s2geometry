namespace S2Geometry;

using Options = S2PolylineLayer.Options;

public class EncodedS2ShapeIndexTests
{
    private readonly ITestOutputHelper _logger;

    internal EncodedS2ShapeIndexTests(ITestOutputHelper logger) { _logger = logger; }

    private static void TestEncodedS2ShapeIndex<Shape>(MutableS2ShapeIndex expected, int expected_bytes)
         where Shape : S2Shape, IInitEncoder<Shape>, new()
    {
        Encoder encoder = new();
        S2ShapeUtilCoding.EncodeHomogeneousShapes(expected, encoder);
        int shapes_bytes = encoder.Length();
        expected.Encode(encoder);
        Assert.Equal(expected_bytes, encoder.Length() - shapes_bytes);
        var decoder = encoder.Decoder();
        var (success, actual) = EncodedS2ShapeIndex.Factory(decoder,
            new S2ShapeUtilCoding.HomogeneousShapeFactory<Shape>(decoder));
        Assert.True(success);
        Assert.Equal(expected.Options_.MaxEdgesPerCell,
            actual!.Options_.MaxEdgesPerCell);
        S2ShapeUtil_Testing.ExpectEqual(expected, actual);
    }

    [Fact]
    internal void Test_EncodedS2ShapeIndex_Empty()
    {
        MutableS2ShapeIndex index = new();
        TestEncodedS2ShapeIndex<S2LaxPolylineShape>(index, 4); // S2LaxPolylineShape
    }

    [Fact]
    internal void Test_EncodedS2ShapeIndex_OneEdge()
    {
        MutableS2ShapeIndex index = new();
        index.Add(MakeLaxPolylineOrDie("1:1, 2:2"));
        TestEncodedS2ShapeIndex<S2LaxPolylineShape>(index, 8); // S2LaxPolylineShape
    }

    [Fact]
    internal void Test_EncodedS2ShapeIndex_RegularLoops()
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
            S2Polygon polygon = new(S2Loop.MakeRegularLoop(new S2Point(3, 2, 1).Normalize(),
                                                      S1Angle.FromDegrees(0.1),
                                                      num_edges));
            index.Add(new S2LaxPolygonShape(polygon));
            TestEncodedS2ShapeIndex<S2LaxPolygonShape>( // EncodedS2LaxPolygonShape
                index, expected_bytes);
        }
    }

    // TODO(b/232496949): This test relies on `random()` return values because
    // it tests an exact encoded byte size.  Either change it to accept a range
    // of sizes, or decode and check either the number of shapes, or possibly
    // the points themselves by resetting the RNG state.
    [Fact]
    internal void Test_EncodedS2ShapeIndex_OverlappingPointClouds()
    {
        (int num_shapes, int num_points_per_shape, int expected_bytes)[] test_cases = {
            (1, 50, 83),
            (2, 100, 583),
            (4, 100, 1383),
        };

        S2Cap cap = new(new S2Point(0.1, -0.4, 0.3).Normalize(), S1Angle.FromDegrees(1));
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

    // TODO(b/232496949): This test relies on `random()` return values.
    [Fact]
    internal void Test_EncodedS2ShapeIndex_OverlappingPolylines()
    {
        (int num_shapes, int num_shape_edges, int expected_bytes)[] test_cases = {
            (2, 50, 139),
            (10, 50, 777),
            (20, 50, 2219),
        };

        S2Cap cap = new(new S2Point(-0.2, -0.3, 0.4).Normalize(), S1Angle.FromDegrees(0.1));
        foreach (var (num_shapes, num_shape_edges, expected_bytes) in test_cases)
        {
            S1Angle edge_len = 2 * cap.RadiusAngle() / num_shape_edges;
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
                    vertices.Add(S2.GetPointOnLine(a, b, j * edge_len));
                }
                index.Add(new S2LaxPolylineShape(vertices.ToArray()));
            }
            TestEncodedS2ShapeIndex<S2LaxPolylineShape>( // S2LaxPolylineShape, EncodedS2LaxPolylineShape
                index, expected_bytes);
        }
    }

    // TODO(b/232496949): This test relies on `random()` return values.
    [Fact]
    internal void Test_EncodedS2ShapeIndex_OverlappingLoops()
    {
        (int num_shapes, int max_edges_per_loop, int expected_bytes)[] test_cases = {
            (2, 250, 138),
            (5, 250, 1084),
            (25, 50, 3673),
        };

        S2Cap cap = new(new S2Point(-0.1, 0.25, 0.2).Normalize(), S1Angle.FromDegrees(3));
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
                    center, cap.RadiusAngle() * radius_fraction, num_edges));
                index.Add(new S2LaxPolygonShape(polygon));
            }
            TestEncodedS2ShapeIndex<S2LaxPolygonShape>( // S2LaxPolygonShape, EncodedS2LaxPolygonShape
                index, expected_bytes);
        }
    }

    [Fact]
    internal void Test_EncodedS2ShapeIndex_SnappedFractalPolylines()
    {
        MutableS2ShapeIndex index = new();
        S2Builder builder = new(new S2Builder.Options(new S2CellIdSnapFunction()));
        for (int i = 0; i < 5; ++i)
        {
            builder.StartLayer(new IndexedLaxPolylineLayer(index));
            S2Testing.Fractal fractal = new();
            fractal.SetLevelForApproxMaxEdges(3 * 256);
            var frame = S2.GetFrame(S2LatLng.FromDegrees(10, i).ToPoint());
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

    [Fact]
    internal void Test_EncodedS2ShapeIndex_LazyDecode()
    {
        // Ensure that lazy decoding is thread-safe.  In other words, make sure that
        // nothing bad happens when multiple threads call "const" methods that cause
        // index and/or shape data to be decoded.
        LazyDecodeTest test = new();

        // The number of readers should be large enough so that it is likely that
        // several readers will be running at once (with a multiple-core CPU).
        const int kNumReaders = 8;
        const int kIters = 1000;
        test.Run(kNumReaders, kIters);
    }

    [Fact]
    internal void Test_EncodedS2ShapeIndex_JavaByteCompatibility()
    {
        MutableS2ShapeIndex expected = new();
        expected.Add(new S2Polyline.OwningShape(MakePolylineOrDie("0:0, 1:1")));
        expected.Add(new S2Polyline.OwningShape(MakePolylineOrDie("1:1, 2:2")));
        expected.Release(0);

        // bytes is the encoded data of an S2ShapeIndex with a null shape and a
        // polyline with one edge. It was derived by base-16 encoding the buffer of
        // an encoder to which expected was encoded.
        var bytes = Convert.FromHexString(
            "100036020102000000B4825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF913FFC7FB8" +
            "B805F6EF3F28516A6D8FDBA13F27DCF7C958DEA13F28C809010408020010");
        Decoder decoder = new(bytes, 0, bytes.Length);
        MutableS2ShapeIndex actual = new();
        Assert.True(actual.Init(decoder, S2ShapeUtilCoding.FullDecodeShapeFactory(decoder)));

        S2ShapeUtil_Testing.ExpectEqual(expected, actual);
    }

    // A test that repeatedly minimizes "index_" in one thread and then reads the
    // index_ concurrently from several other threads.  When all threads have
    // finished reading, the first thread minimizes the index again.
    //
    // Note that Minimize() is non-const and therefore does not need to be tested
    // concurrently with the const methods.
    internal class LazyDecodeTest : ReaderWriterTest
    {
        private readonly EncodedS2ShapeIndex index_;

        internal LazyDecodeTest()
        {
            // We generate one shape per dimension.  Each shape has vertices uniformly
            // distributed across the sphere, and the vertices for each dimension are
            // different.  Having fewer cells in the index is more likely to trigger
            // race conditions, and so shape 0 has 384 points, shape 1 is a polyline
            // with 96 vertices, and shape 2 is a polygon with 24 vertices.
            MutableS2ShapeIndex input = new();
            for (int dim = 0; dim < 3; ++dim)
            {
                int level = 3 - dim;  // See comments above.
                var verticesTmp = new List<S2Point>();
                for (var id = S2CellId.Begin(level); id != S2CellId.End(level); id = id.Next())
                {
                    verticesTmp.Add(id.ToPoint());
                }
                var vertices = verticesTmp.ToArray();
                switch (dim)
                {
                    case 0: input.Add(new S2PointVectorShape(vertices)); break;
                    case 1: input.Add(new S2LaxPolylineShape(vertices)); break;
                    default:
                        input.Add(new S2LaxPolygonShape(new[] { verticesTmp }.ToList()));
                        break;
                }
            }
            Encoder encoder = new();
            Assert.True(S2ShapeUtilCoding.CompactEncodeTaggedShapes(input, encoder));
            input.Encode(encoder);
            var decoder = encoder.Decoder();
            var (success, index) = EncodedS2ShapeIndex.Factory(decoder,
                S2ShapeUtilCoding.LazyDecodeShapeFactory(decoder));
            Assert.True(success);
            index_ = index!;
        }

        internal override void WriteOp()
        {
            index_.Minimize();
        }

        internal override void ReadOp()
        {
            S2ClosestEdgeQuery query = new(index_);
            for (int iter = 0; iter < 10; ++iter)
            {
                S2ClosestEdgeQuery.PointTarget target = new(S2Testing.RandomPoint());
                query.FindClosestEdge(target);
            }
        }
    }

    // Like S2PolylineLayer, but converts the polyline to an S2LaxPolylineShape
    // and adds it to an S2ShapeIndex (if the polyline is non-empty).
    private class IndexedLaxPolylineLayer : Layer
    {
        private readonly MutableS2ShapeIndex index_;
        private readonly S2Polyline polyline_;
        private readonly S2PolylineLayer layer_;

        internal IndexedLaxPolylineLayer(MutableS2ShapeIndex index, Options? options = null)
        {
            index_ = index;
            polyline_ = new();
            layer_ = new(polyline_, options ?? new Options());
        }

        public override GraphOptions GraphOptions_()
        {
            return layer_.GraphOptions_();
        }

        public override void Build(Graph g, out S2Error error)
        {
            layer_.Build(g, out error);
            if (error.IsOk() && polyline_.NumVertices() > 0)
            {
                index_.Add(new S2LaxPolylineShape(polyline_));
            }
        }
    }
}
