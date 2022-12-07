namespace S2Geometry;

public class S2LaxPolylineShapeTests
{
    [Fact]
    internal void Test_S2LaxPolylineShape_NoVertices()
    {
        var vertices = Array.Empty<S2Point>();
        var shape = new S2LaxPolylineShape(vertices);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(1, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2LaxPolylineShape_OneVertex()
    {
        S2Point[] vertices = { new S2Point(1, 0, 0) };
        var shape = new S2LaxPolylineShape(vertices);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(1, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
    }

    [Fact]
    internal void Test_S2LaxPolylineShape_Move()
    {
        // Construct a shape to use as the correct answer and a second identical shape
        // to be moved.
        List<S2Point> vertices = ParsePointsOrDie("1:1, 4:4, 2:2, 3:3");
        S2LaxPolylineShape correct=new(vertices);
        S2LaxPolylineShape to_move=new(vertices);

        // Test the move constructor.
        S2LaxPolylineShape move1=to_move;
        Assert.Equal(correct, move1);
        Assert.Equal(correct.Id, move1.Id);
        Assert.Equal(vertices.Count, move1.NumVertices());
        for (int i = 0; i < move1.NumVertices(); ++i)
        {
            Assert.Equal(vertices[i], move1.Vertex(i));
        }
        // Test the move-assignment operator.
        S2LaxPolylineShape move2;
        move2 = move1;
        Assert.Equal(correct, move2);
        Assert.Equal(correct.Id, move2.Id);
        Assert.Equal(vertices.Count, move2.NumVertices());
        for (int i = 0; i < move2.NumVertices(); ++i)
        {
            Assert.Equal(vertices[i], move2.Vertex(i));
        }
    }

    [Fact]
    internal void Test_S2LaxPolylineShape_EdgeAccess()
    {
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1").ToArray();
        S2LaxPolylineShape shape = new(vertices);
        Assert.Equal(2, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(0, shape.GetChain(0).Start);
        Assert.Equal(2, shape.GetChain(0).Length);
        Assert.Equal(1, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        var edge0 = shape.GetEdge(0);
        Assert.Equal(vertices[0], edge0.V0);
        Assert.Equal(vertices[1], edge0.V1);
        var edge1 = shape.GetEdge(1);
        Assert.Equal(vertices[1], edge1.V0);
        Assert.Equal(vertices[2], edge1.V1);
    }

    [Fact]
    internal void Test_EncodedS2LaxPolylineShape_RoundtripEncoding()
    {
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1");
        S2LaxPolylineShape shape=new(vertices.ToArray());

        Encoder encoder=new();
        shape.Encode(encoder, CodingHint.COMPACT);
        var a_decoder = encoder.Decoder();
        var (a_success, a_shape) = EncodedS2LaxPolylineShape.Init(a_decoder);
        Assert.True(a_success);

        Encoder b_encoder=new();
        a_shape!.Encode(b_encoder, CodingHint.COMPACT);
        var b_decoder = b_encoder.Decoder();
        var (b_success, b_shape) = EncodedS2LaxPolylineShape.Init(b_decoder);
        Assert.True(b_success);
        S2ShapeUtil_Testing.ExpectEqual(shape, b_shape!);
    }
}
