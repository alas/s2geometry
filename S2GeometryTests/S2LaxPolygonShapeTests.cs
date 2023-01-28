namespace S2Geometry;

using System.Runtime.InteropServices;
using Loop = List<S2Point>;

public class S2LaxPolygonShapeTests
{
    private readonly ITestOutputHelper _logger;

    public S2LaxPolygonShapeTests(ITestOutputHelper logger) { _logger = logger; }

    // Verifies that EncodedS2LaxPolygonShape behaves identically to
    // S2LaxPolygonShape. Also supports testing that the encoded form is identical
    // to the re-encoded form.

    private static void TestEncodedS2LaxPolygonShape(S2LaxPolygonShape original)
    {
        Encoder encoder=new();
        original.Encode(encoder, CodingHint.COMPACT);
        var decoder = encoder.GetDecoder();
        var (success, encoded) = EncodedS2LaxPolygonShape.Init(decoder);
        Assert.True(success);
        Assert.Equal(encoded!.NumLoops, original.NumLoops);
        Assert.Equal(encoded.NumVertices, original.NumVertices);
        Assert.Equal(encoded.NumEdges(), original.NumEdges());
        Assert.Equal(encoded.NumChains(), original.NumChains());
        Assert.Equal(encoded.Dimension(), original.Dimension());
        Assert.Equal(encoded.IsEmpty(), original.IsEmpty());
        Assert.Equal(encoded.IsFull(), original.IsFull());
        Assert.Equal(encoded.GetReferencePoint(), original.GetReferencePoint());
        for (int i = 0; i<original.NumLoops; ++i) {
            Assert.Equal(encoded.NumLoopVertices(i), original.NumLoopVertices(i));
            Assert.Equal(encoded.GetChain(i), original.GetChain(i));
            for (int j = 0; j<original.NumLoopVertices(i); ++j) {
                Assert.Equal(encoded.LoopVertex(i, j), original.LoopVertex(i, j));
                Assert.Equal(encoded.ChainEdge(i, j), original.ChainEdge(i, j));
            }
        }
        // Now test all the edges in a random order in order to exercise the cases
        // involving prev_loop_.
        var count = original.NumEdges();
        List<int> edge_ids=new(count);
        LinqUtils.Iota(edge_ids, 0, count);
        var mt = new PseudoRandom.MersenneTwister();
        edge_ids = edge_ids.Shuffle(mt.GenerateRandomN).ToList();
        foreach (var e in edge_ids)
        {
            Assert.Equal(encoded.GetChainPosition(e), original.GetChainPosition(e));
            Assert.Equal(encoded.GetEdge(e), original.GetEdge(e));
        }

        // Let's also test that the encoded form can be encoded, yielding the same
        // bytes as the originally encoded form.
        Encoder reencoder=new();
        encoded.Encode(reencoder, CodingHint.COMPACT);
        Assert.True(encoder == reencoder);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_EmptyPolygon()
    {
        _logger.WriteLine(
            $"(INFO) sizeof(S2LaxPolygonShape) == {SizeHelper.SizeOf(typeof(S2LaxPolygonShape))}");
        _logger.WriteLine(
            $"(INFO) sizeof(EncodedS2LaxPolygonShape) == {SizeHelper.SizeOf(typeof(EncodedS2LaxPolygonShape))}");

        S2LaxPolygonShape shape = new(new S2Polygon());
        Assert.Equal(0, shape.NumLoops);
        Assert.Equal(0, shape.NumVertices);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(2, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_Move()
    {
        // Construct a shape to use as the correct answer and a second identical shape
        // to be moved.
        List<Loop> loops = new()
        {
            new Loop(ParsePointsOrDie("0:0, 0:3, 3:3")),
            new Loop(ParsePointsOrDie("1:1, 2:2, 1:2"))
        };
        S2LaxPolygonShape correct=new(loops);
        S2LaxPolygonShape to_move=new(loops);
        Assert.Equal(correct, to_move);
        Assert.Equal(correct.Id, to_move.Id);

        // Test the move constructor.
        S2LaxPolygonShape move1=to_move;
        Assert.Equal(correct, move1);
        TestEncodedS2LaxPolygonShape(move1);
        Assert.Equal(correct.Id, move1.Id);
        Assert.Equal(loops.Count, move1.NumLoops);
        Assert.Equal(6, move1.NumVertices);
        for (int i = 0; i < loops.Count; ++i)
        {
            for (int j = 0; j < loops[i].Count; ++j)
            {
                Assert.Equal(loops[i][j], move1.LoopVertex(i, j));
            }
        }

        // Test the move-assignment operator.
        S2LaxPolygonShape move2;
        move2 = move1;
        Assert.Equal(correct, move2);
        TestEncodedS2LaxPolygonShape(move2);
        Assert.Equal(correct.Id, move2.Id);
        Assert.Equal(loops.Count, move2.NumLoops);
        Assert.Equal(6, move2.NumVertices);
        for (int i = 0; i < loops.Count; ++i)
        {
            for (int j = 0; j < loops[i].Count; ++j)
            {
                Assert.Equal(loops[i][j], move2.LoopVertex(i, j));
            }
        }
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_MoveFromShapeIndex()
    {
        // Construct an index containing shapes to be moved.
        List<Loop> loops = new()
        {
            new Loop(ParsePointsOrDie("0:0, 0:3, 3:3")),
            new Loop(ParsePointsOrDie("1:1, 2:2, 1:2"))
        };
        MutableS2ShapeIndex index=new();
        index.Add(new S2LaxPolygonShape(loops));
        index.Add(new S2LaxPolygonShape(loops));
        Assert.Equal(index.NumShapeIds(), 2);

        // Verify that the move constructor moves the id.
        var shape0 = (S2LaxPolygonShape)index.Shape(0)!;
        S2LaxPolygonShape moved_shape0 = shape0;
        Assert.Equal(moved_shape0.Id, 0);

        // Verify that the move-assignment operator moves the id.
        var shape1 = (S2LaxPolygonShape)index.Shape(1)!;
        S2LaxPolygonShape moved_shape1;
        moved_shape1 = shape1;
        Assert.Equal(moved_shape1.Id, 1);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_FullPolygon()
    {
        var shape = new S2LaxPolygonShape(new S2Polygon(MakeLoopOrDie("full")));
        Assert.Equal(1, shape.NumLoops);
        Assert.Equal(0, shape.NumVertices);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.True(shape.IsFull());
        Assert.True(shape.GetReferencePoint().Contained);
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_SingleVertexPolygon()
    {
        // S2Polygon doesn't support single-vertex loops, so we need to construct
        // the S2LaxPolygonShape directly.
        var loops = new List<List<S2Point>>
            {
                ParsePointsOrDie("0:0")
            };
        var shape = new S2LaxPolygonShape(loops);
        Assert.Equal(1, shape.NumLoops);
        Assert.Equal(1, shape.NumVertices);
        Assert.Equal(1, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(0, shape.GetChain(0).Start);
        Assert.Equal(1, shape.GetChain(0).Length);
        var edge = shape.GetEdge(0);
        Assert.Equal(loops[0][0], edge.V0);
        Assert.Equal(loops[0][0], edge.V1);
        Assert.True(edge == shape.ChainEdge(0, 0));
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_SingleLoopPolygon()
    {
        // Test S2Polygon constructor.
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        var shape = new S2LaxPolygonShape(new S2Polygon(new S2Loop(vertices)));
        Assert.Equal(1, shape.NumLoops);
        Assert.Equal(vertices.Count, shape.NumVertices);
        Assert.Equal(vertices.Count, shape.NumLoopVertices(0));
        Assert.Equal(vertices.Count, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(0, shape.GetChain(0).Start);
        Assert.Equal(vertices.Count, shape.GetChain(0).Length);
        for (int i = 0; i < vertices.Count; ++i)
        {
            Assert.Equal(vertices[i], shape.LoopVertex(0, i));
            var edge = shape.GetEdge(i);
            Assert.Equal(vertices[i], edge.V0);
            Assert.Equal(vertices[(i + 1) % vertices.Count], edge.V1);
            Assert.Equal(edge.V0, shape.ChainEdge(0, i).V0);
            Assert.Equal(edge.V1, shape.ChainEdge(0, i).V1);
        }
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.ContainsBruteForce(S2.Origin));
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_MultiLoopPolygon()
    {
        // Test S2Point[][] constructor.  Make sure that the loops are
        // oriented so that the interior of the polygon is always on the left.
        List<List<S2Point>> loops = new()
        {
            ParsePointsOrDie("0:0, 0:3, 3:3"),  // CCW
            ParsePointsOrDie("1:1, 2:2, 1:2")   // CW
        };
        S2LaxPolygonShape shape = new(loops);

        Assert.Equal(loops.Count, shape.NumLoops);
        int num_vertices = 0;
        Assert.Equal(loops.Count, shape.NumChains());
        for (int i = 0; i < loops.Count; ++i)
        {
            Assert.Equal(loops[i].Count, shape.NumLoopVertices(i));
            Assert.Equal(num_vertices, shape.GetChain(i).Start);
            Assert.Equal(loops[i].Count, shape.GetChain(i).Length);
            for (int j = 0; j < loops[i].Count; ++j)
            {
                Assert.Equal(loops[i][j], shape.LoopVertex(i, j));
                var edge = shape.GetEdge(num_vertices + j);
                Assert.Equal(loops[i][j], edge.V0);
                Assert.Equal(loops[i][(j + 1) % loops[i].Count], edge.V1);
            }
            num_vertices += loops[i].Count;
        }
        Assert.Equal(num_vertices, shape.NumVertices);
        Assert.Equal(num_vertices, shape.NumEdges());
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.ContainsBruteForce(S2.Origin));
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_MultiLoopS2Polygon()
    {
        // Verify that the orientation of loops representing holes is reversed when
        // converting from an S2Polygon to an S2LaxPolygonShape.
        var polygon = MakePolygonOrDie("0:0, 0:3, 3:3; 1:1, 1:2, 2:2");
        var shape = new S2LaxPolygonShape(polygon);
        for (int i = 0; i < polygon.NumLoops(); ++i)
        {
            var loop = polygon.Loop(i);
            for (int j = 0; j < loop.NumVertices; ++j)
            {
                Assert.Equal(loop.OrientedVertex(j),
                          shape.LoopVertex(i, j));
            }
        }
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_ManyLoopPolygon()
    {
        // Test a polygon with enough loops so that binary search is used to find
        // the loop containing a given edge.
        var loops = new List<List<S2Point>>();
        for (int i = 0; i < 100; ++i)
        {
            var center = S2LatLng.FromDegrees(0, i).ToPoint();
            var loop = S2Testing.MakeRegularPoints(
                center, S1Angle.FromDegrees(0.1),
                S2Testing.Random.Uniform(3));
            loops.Add(loop.ToList());
        }
        var shape = new S2LaxPolygonShape(loops);

        Assert.Equal(loops.Count, shape.NumLoops);
        int num_vertices = 0;
        Assert.Equal(loops.Count, shape.NumChains());
        for (int i = 0; i < loops.Count; ++i)
        {
            Assert.Equal(loops[i].Count, shape.NumLoopVertices(i));
            Assert.Equal(num_vertices, shape.GetChain(i).Start);
            Assert.Equal(loops[i].Count, shape.GetChain(i).Length);
            for (int j = 0; j < loops[i].Count; ++j)
            {
                Assert.Equal(loops[i][j], shape.LoopVertex(i, j));
                int e = num_vertices + j;
                Assert.Equal(shape.GetChainPosition(e), new S2Shape.ChainPosition(i, j));
                Assert.Equal(loops[i][j], shape.GetEdge(e).V0);
                Assert.Equal(loops[i][(j + 1) % loops[i].Count], shape.GetEdge(e).V1);
            }
            num_vertices += loops[i].Count;
        }
        Assert.Equal(num_vertices, shape.NumVertices);
        Assert.Equal(num_vertices, shape.NumEdges());

        // Now test all the edges in a random order in order to exercise the cases
        // involving prev_loop_.
        List<(int, int, int)> edges=new();
        for (int i = 0, e = 0; i < loops.Count; ++i)
        {
            for (int j = 0; j < loops[i].Count; ++j, ++e)
            {
                edges.Add((e, i, j));
            }
        }
        var mt = new PseudoRandom.MersenneTwister();
        edges = edges.Shuffle(mt.GenerateRandomN).ToList();
        // TODO(user,b/210097200): Use structured bindings when we require
        // C++17 in opensource.
        foreach (var (e, i, j) in edges)
        {
            Assert.Equal(shape.GetChainPosition(e), new S2Shape.ChainPosition(i, j));
            var v0 = loops[i][j];
            var v1 = loops[i][(j + 1) % loops[i].Count];
            Assert.Equal(shape.GetEdge(e), new S2Shape.Edge(v0, v1));
        }
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_DegenerateLoops()
    {
        var loops = new List<List<S2Point>>{
            ParsePointsOrDie("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
            ParsePointsOrDie("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
            ParsePointsOrDie("5:5, 6:6")
        };
        var shape = new S2LaxPolygonShape(loops);
        Assert.False(shape.GetReferencePoint().Contained);
        TestEncodedS2LaxPolygonShape(shape);
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_InvertedLoops()
    {
        var loops = new List<List<S2Point>>{
            ParsePointsOrDie("1:2, 1:1, 2:2"),
            ParsePointsOrDie("3:4, 3:3, 4:4")
        };
        var shape = new S2LaxPolygonShape(loops);
        Assert.True(shape.ContainsBruteForce(S2.Origin));
        TestEncodedS2LaxPolygonShape(shape);
    }

    private static void CompareS2LoopToShape(S2Loop loop, S2Shape shape)
    {
        MutableS2ShapeIndex index = new();
        index.Add(shape);
        S2Cap cap = loop.GetCapBound();
        var query = index.MakeS2ContainsPointQuery();
        for (int iter = 0; iter < 100; ++iter)
        {
            S2Point point = S2Testing.SamplePoint(cap);
            Assert.Equal(loop.Contains(point),
                query.ShapeContains(index.Shape(0), point));
        }
    }

    [Fact]
    internal void Test_S2LaxPolygonShape_CompareToS2Loop()
    {
        for (int iter = 0; iter < 100; ++iter)
        {
            var fractal = new S2Testing.Fractal();
            fractal.MaxLevel = (S2Testing.Random.Uniform(5));
            fractal.FractalDimension = (1 + S2Testing.Random.RandDouble());
            S2Point center = S2Testing.RandomPoint();
            var loop = fractal.MakeLoop(
                S2Testing.GetRandomFrameAt(center), S1Angle.FromDegrees(5));

            // Compare S2Loop to S2LaxLoopShape.
            CompareS2LoopToShape(loop, new S2LaxLoopShape(loop));

            // Compare S2Loop to S2LaxPolygonShape.
            var loops = new List<List<S2Point>> { loop.CloneVertices().ToList() };
            CompareS2LoopToShape(loop, new S2LaxPolygonShape(loops));
        }
    }
}
