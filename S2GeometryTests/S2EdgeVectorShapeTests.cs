namespace S2Geometry;

public class S2EdgeVectorShapeTests
{
    [Fact]
    internal void Test_S2EdgeVectorShape_Empty()
    {
        S2EdgeVectorShape shape = new();
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(1, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2EdgeVectorShape_Move()
    {
        // Construct a shape to use as the correct answer and a second identical shape
        // to be moved.
        S2EdgeVectorShape correct=new();
        S2EdgeVectorShape to_move=new();
        S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
        const int kNumEdges = 100;
        for (int i = 0; i < kNumEdges; ++i)
        {
            var start_point = S2Testing.RandomPoint();
            var end_point = S2Testing.RandomPoint();
            correct.Add(start_point, end_point);
            to_move.Add(start_point, end_point);
        }

        // Test the move constructor.
        S2EdgeVectorShape move1 = to_move;
        Assert.Equal(correct, move1);
        Assert.Equal(correct.Id, move1.Id);

        // Test the move-assignment operator.
        S2EdgeVectorShape move2;
        move2 = move1;
        Assert.Equal(correct, move2);
        Assert.Equal(correct.Id, move2.Id);
    }

    [Fact]
    internal void Test_S2EdgeVectorShape_EdgeAccess()
    {
        S2EdgeVectorShape shape = new();
        S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
        int kNumEdges = 100;
        List<(S2Point, S2Point)> edges = new();
        for (int i = 0; i < kNumEdges; ++i)
        {
            S2Point a = S2Testing.RandomPoint();  // Control the evaluation order
            edges.Add((a, S2Testing.RandomPoint()));
            shape.Add(edges.Last().Item1, edges.Last().Item2);
        }
        Assert.Equal(kNumEdges, shape.NumEdges());
        Assert.Equal(kNumEdges, shape.NumChains());
        Assert.Equal(1, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        for (int i = 0; i < kNumEdges; ++i)
        {
            Assert.Equal(i, shape.GetChain(i).Start);
            Assert.Equal(1, shape.GetChain(i).Length);
            var edge = shape.GetEdge(i);
            Assert.Equal(edges[i].Item1, edge.V0);
            Assert.Equal(edges[i].Item2, edge.V1);
        }
    }

    [Fact]
    internal void Test_S2EdgeVectorShape_SingletonConstructor()
    {
        S2Point a = new(1, 0, 0), b = new(0, 1, 0);
        S2EdgeVectorShape shape = new(a, b);
        Assert.Equal(1, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        var edge = shape.GetEdge(0);
        Assert.Equal(a, edge.V0);
        Assert.Equal(b, edge.V1);
    }
}
