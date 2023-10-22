namespace S2Geometry;

public class S2PointVectorShapeTests
{
    [Fact]
    internal void Test_S2PointVectorShape_Empty()
    {
        S2Point[] points = [];
        S2PointVectorShape shape = new(points);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.Equal(0, shape.Dimension());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2PointVectorShape_ConstructionAndAccess()
    {
        const int kNumPoints = 100;
        S2Point[] points = new S2Point[kNumPoints];
        S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
        for (int i = 0; i < kNumPoints; ++i)
        {
            points[i] = S2Testing.RandomPoint();
        }
        S2PointVectorShape shape = new(points);

        Assert.Equal(kNumPoints, shape.NumEdges());
        Assert.Equal(kNumPoints, shape.NumChains());
        Assert.Equal(0, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        for (int i = 0; i < kNumPoints; ++i)
        {
            Assert.Equal(i, shape.GetChain(i).Start);
            Assert.Equal(1, shape.GetChain(i).Length);
            var edge = shape.GetEdge(i);
            S2Point pt = points[i];
            Assert.Equal(pt, edge.V0);
            Assert.Equal(pt, edge.V1);
        }
    }

    [Fact]
    internal void Test_S2PointVectorShape_Move()
    {
        // Construct a shape to use as the correct answer and a second identical shape
        // to be moved.
        List<S2Point> points=[];
        const int kNumPoints = 100;
        for (int i = 0; i < kNumPoints; ++i)
        {
            points.Add(S2Testing.RandomPoint());
        }
        S2PointVectorShape correct=new([.. points]);
        S2PointVectorShape to_move=new([.. points]);

        // Test the move constructor.
        S2PointVectorShape move1=to_move;
        Assert.Equal(correct, move1);
        Assert.Equal(correct.Id, move1.Id);

        // Test the move-assignment operator.
        S2PointVectorShape move2;
        move2 = move1;
        Assert.Equal(correct, move2);
        Assert.Equal(correct.Id, move2.Id);
    }
}
