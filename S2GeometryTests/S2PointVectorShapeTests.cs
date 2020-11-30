
using Xunit;

namespace S2Geometry
{
    public class S2PointVectorShapeTests
    {
        [Fact]
        public void Test_S2PointVectorShape_Empty() {
            S2Point[] points = System.Array.Empty<S2Point>();
            S2PointVectorShape shape = new S2PointVectorShape(points);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(0, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2PointVectorShape_ConstructionAndAccess() {
            const int kNumPoints = 100;
            S2Point[] points = new S2Point[kNumPoints];
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int i = 0; i < kNumPoints; ++i) {
                points[i] = S2Testing.RandomPoint();
            }
            S2PointVectorShape shape = new S2PointVectorShape(points);

            Assert.Equal(kNumPoints, shape.NumEdges);
            Assert.Equal(kNumPoints, shape.NumChains());
            Assert.Equal(0, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int i = 0; i < kNumPoints; ++i) {
                Assert.Equal(i, shape.GetChain(i).Start);
                Assert.Equal(1, shape.GetChain(i).Length);
                var edge = shape.GetEdge(i);
                S2Point pt = S2Testing.RandomPoint();
                Assert.Equal(pt, edge.V0);
                Assert.Equal(pt, edge.V1);
            }
        }
    }
}
