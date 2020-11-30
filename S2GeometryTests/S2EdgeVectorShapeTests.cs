using Xunit;

namespace S2Geometry
{
    public class S2EdgeVectorShapeTests
    {
        [Fact]
        public void Test_S2EdgeVectorShape_Empty()
        {
            S2EdgeVectorShape shape = new();
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(1, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2EdgeVectorShape_EdgeAccess()
        {
            S2EdgeVectorShape shape = new();
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            int kNumEdges = 100;
            for (int i = 0; i < kNumEdges; ++i)
            {
                S2Point a = S2Testing.RandomPoint();  // Control the evaluation order
                shape.Add(a, S2Testing.RandomPoint());
            }
            Assert.Equal(kNumEdges, shape.NumEdges);
            Assert.Equal(kNumEdges, shape.NumChains());
            Assert.Equal(1, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int i = 0; i < kNumEdges; ++i)
            {
                Assert.Equal(i, shape.GetChain(i).Start);
                Assert.Equal(1, shape.GetChain(i).Length);
                var edge = shape.GetEdge(i);
                Assert.Equal(S2Testing.RandomPoint(), edge.V0);
                Assert.Equal(S2Testing.RandomPoint(), edge.V1);
            }
        }

        [Fact]
        public void Test_S2EdgeVectorShape_SingletonConstructor()
        {
            S2Point a = new(1, 0, 0), b = new(0, 1, 0);
            S2EdgeVectorShape shape = new(a, b);
            Assert.Equal(1, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            var edge = shape.GetEdge(0);
            Assert.Equal(a, edge.V0);
            Assert.Equal(b, edge.V1);
        }
    }
}
