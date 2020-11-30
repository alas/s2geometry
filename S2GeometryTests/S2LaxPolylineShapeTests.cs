using Xunit;

namespace S2Geometry
{
    public class S2LaxPolylineShapeTests
    {
        [Fact]
        public void Test_S2LaxPolylineShape_NoVertices() {
            var vertices = System.Array.Empty<S2Point>();
            var shape = new S2LaxPolylineShape(vertices);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(1, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxPolylineShape_OneVertex() {
            S2Point[] vertices = { new S2Point(1, 0, 0) };
            var shape = new S2LaxPolylineShape(vertices);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(1, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
        }

        [Fact]
        public void Test_S2LaxPolylineShape_EdgeAccess() {
            var vertices = S2TextFormat.ParsePointsOrDie("0:0, 0:1, 1:1").ToArray();
            S2LaxPolylineShape shape = new S2LaxPolylineShape(vertices);
            Assert.Equal(2, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(2, shape.GetChain(0).Length);
            Assert.Equal(1, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            var edge0 = shape.GetEdge(0);
            Assert.Equal(vertices[0], edge0.V0);
            Assert.Equal(vertices[1], edge0.V1);
            var edge1 = shape.GetEdge(1);
            Assert.Equal(vertices[1], edge1.V0);
            Assert.Equal(vertices[2], edge1.V1);
        }
    }
}
