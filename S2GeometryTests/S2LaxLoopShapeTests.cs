namespace S2Geometry;

public class S2LaxLoopShapeTests
{
    [Fact]
    internal void Test_S2LaxLoopShape_EmptyLoop()
    {
        // Test S2Loop constructor.
        var shape = new S2LaxLoopShape(S2Loop.KEmpty);
        Assert.Equal(0, shape.NumVertices);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(2, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2LaxLoopShape_Move()
    {
        // Construct a shape to use as the correct answer and a second identical shape
        // to be moved.
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        S2LaxLoopShape correct = new([.. vertices]);
        S2LaxLoopShape to_move = new([.. vertices]);

        // Test the move constructor.
        S2LaxLoopShape move1 = to_move;
        Assert.Equal(correct, move1);
        Assert.Equal(correct.Id, move1.Id);
        Assert.Equal(correct.NumVertices, move1.NumVertices);
        for (int i = 0; i < correct.NumVertices; ++i)
        {
            Assert.Equal(correct.Vertex(i), move1.Vertex(i));
        }

        // Test the move-assignment operator.
        S2LaxLoopShape move2 = new();
        move2 = move1;
        Assert.Equal(correct, move2);
        Assert.Equal(correct.Id, move2.Id);
        Assert.Equal(correct.NumVertices, move2.NumVertices);
        for (int i = 0; i < correct.NumVertices; ++i)
        {
            Assert.Equal(correct.Vertex(i), move2.Vertex(i));
        }
    }

    [Fact]
    internal void Test_S2LaxLoopShape_MoveFromShapeIndex()
    {
        // Construct an index containing shapes to be moved.
        MutableS2ShapeIndex index =
        [
            new S2LaxLoopShape([.. ParsePointsOrDie("0:0, 0:1, 1:1, 1:0")]),
            new S2LaxLoopShape([.. ParsePointsOrDie("0:0, 0:2, 2:2, 2:0")]),
        ];
        Assert.Equal(index.NumShapeIds(), 2);

        // Verify that the move constructor moves the id.
        var shape0 = (S2LaxLoopShape)index.Shape(0)!;
        S2LaxLoopShape moved_shape0 = shape0;
        Assert.Equal(moved_shape0.Id, 0);

        // Verify that the move-assignment operator moves the id.
        var shape1 = (S2LaxLoopShape)index.Shape(1)!;
        S2LaxLoopShape moved_shape1;
        moved_shape1 = shape1;
        Assert.Equal(moved_shape1.Id, 1);
    }

    [Fact]
    internal void Test_S2LaxLoopShape_NonEmptyLoop()
    {
        // Test S2Point[] constructor.
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        var shape = new S2LaxLoopShape([.. vertices]);
        Assert.Equal(vertices.Count, shape.NumVertices);
        Assert.Equal(vertices.Count, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(0, shape.GetChain(0).Start);
        Assert.Equal(vertices.Count, shape.GetChain(0).Length);
        for (int i = 0; i < vertices.Count; ++i)
        {
            Assert.Equal(vertices[i], shape.Vertex(i));
            var edge = shape.GetEdge(i);
            Assert.Equal(vertices[i], edge.V0);
            Assert.Equal(vertices[(i + 1) % vertices.Count], edge.V1);
        }
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2LaxClosedPolylineShape_NoInterior()
    {
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        var shape = new S2LaxClosedPolylineShape([.. vertices]);
        Assert.Equal(1, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2VertexIdLaxLoopShape_EmptyLoop()
    {
        S2VertexIdLaxLoopShape shape = new([], null);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumVertices);
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(2, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2VertexIdLaxLoopShape_Move()
    {
        // Construct a shape to use as the correct answer and a second identical shape
        // to be moved.
        var vertices = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        Int32[] vertex_ids = [0, 3, 2, 1];  // Inverted.
        S2VertexIdLaxLoopShape correct = new(vertex_ids, [.. vertices]);
        S2VertexIdLaxLoopShape to_move = new(vertex_ids, [.. vertices]);

        // Test the move constructor.
        S2VertexIdLaxLoopShape move1 = to_move;
        Assert.Equal(correct, move1);
        Assert.Equal(correct.Id, move1.Id);
        Assert.Equal(correct.NumVertices, move1.NumVertices);
        for (int i = 0; i < correct.NumVertices; ++i)
        {
            Assert.Equal(correct.Vertex(i), move1.Vertex(i));
        }

        // Test the move-assignment operator.
        S2VertexIdLaxLoopShape move2;
        move2 = move1;
        Assert.Equal(correct, move2);
        Assert.Equal(correct.Id, move2.Id);
        Assert.Equal(correct.NumVertices, move2.NumVertices);
        for (int i = 0; i < correct.NumVertices; ++i)
        {
            Assert.Equal(correct.Vertex(i), move2.Vertex(i));
        }
    }

    [Fact]
    internal void Test_S2VertexIdLaxLoopShape_MoveFromShapeIndex()
    {
        // Setup vertices and vertex ids.
        var vertices0 = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        Int32[] vertex_ids0 = [0, 3, 2, 1];
        var vertices1 = ParsePointsOrDie("0:0, 0:2, 2:2, 2:0");
        Int32[] vertex_ids1 = [0, 3, 2, 1];

        // Construct an index containing shapes to be moved.
        MutableS2ShapeIndex index =
        [
            new S2VertexIdLaxLoopShape(vertex_ids0, [.. vertices0]),
            new S2VertexIdLaxLoopShape(vertex_ids1, [.. vertices1]),
        ];
        Assert.Equal(index.NumShapeIds(), 2);

        // Verify that the move constructor moves the id.
        var shape0 = (S2VertexIdLaxLoopShape)index.Shape(0)!;
        S2VertexIdLaxLoopShape moved_shape0 = shape0;
        Assert.Equal(moved_shape0.Id, 0);

        // Verify that the move-assignment operator moves the id.
        var shape1 = (S2VertexIdLaxLoopShape)index.Shape(1)!;
        S2VertexIdLaxLoopShape moved_shape1;
        moved_shape1 = shape1;
        Assert.Equal(moved_shape1.Id, 1);
    }

    [Fact]
    internal void Test_S2VertexIdLaxLoopShape_InvertedLoop()
    {
        var vertex_array = ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
        var vertex_ids = new[] { 0, 3, 2, 1 };  // Inverted.
        var shape = new S2VertexIdLaxLoopShape(vertex_ids, [.. vertex_array]);
        Assert.Equal(4, shape.NumEdges());
        Assert.Equal(4, shape.NumVertices);
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(0, shape.GetChain(0).Start);
        Assert.Equal(4, shape.GetChain(0).Length);
        Assert.Equal(vertex_array[0], shape.Vertex(0));
        Assert.Equal(vertex_array[3], shape.Vertex(1));
        Assert.Equal(vertex_array[2], shape.Vertex(2));
        Assert.Equal(vertex_array[1], shape.Vertex(3));
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.True(shape.ContainsBruteForce(S2.Origin));
    }
}
