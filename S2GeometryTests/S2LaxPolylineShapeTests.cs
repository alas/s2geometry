namespace S2Geometry;

public class S2LaxPolylineShapeTests
{
    [Fact]
    public void Test_S2LaxPolylineShape_NoVertices()
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
    public void Test_S2LaxPolylineShape_OneVertex()
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
    public void Test_S2LaxPolylineShape_MoveConstructor()
    {
        var original = MakeLaxPolylineOrDie("1:1, 4:4");
        S2LaxPolylineShape moved = new(original);
        Assert.Equal(0, original.NumVertices());
        Assert.Equal(2, moved.NumVertices());
    }

    [Fact]
    public void Test_S2LaxPolylineShape_MoveAssignmentOperator()
    {
        var original = MakeLaxPolylineOrDie("1:1, 4:4");
        S2LaxPolylineShape moved;
        moved = original;
        original = new();
        Assert.Equal(0, original.NumVertices());
        Assert.Equal(2, moved.NumVertices());
    }

    [Fact]
    public void Test_S2LaxPolylineShape_EdgeAccess()
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
    public void Test_EncodedS2LaxPolylineShape_RoundtripEncoding()
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
