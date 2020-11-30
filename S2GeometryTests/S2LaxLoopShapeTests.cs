using S2Geometry.S2ShapeUtil;
using Xunit;

namespace S2Geometry
{
    public class S2LaxLoopShapeTests
    {
        [Fact]
        public void Test_S2LaxLoopShape_EmptyLoop() {
            // Test S2Loop constructor.
            var shape = new S2LaxLoopShape(S2Loop.kEmpty);
            Assert.Equal(0, shape.NumVertices);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(2, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxLoopShape_NonEmptyLoop() {
            // Test S2Point[] constructor.
            var vertices = S2TextFormat.ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
            var shape = new S2LaxLoopShape(vertices.ToArray());
            Assert.Equal(vertices.Count, shape.NumVertices);
            Assert.Equal(vertices.Count, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(vertices.Count, shape.GetChain(0).Length);
            for (int i = 0; i < vertices.Count; ++i) {
                Assert.Equal(vertices[i], shape.Vertex(i));
                var edge = shape.GetEdge(i);
                Assert.Equal(vertices[i], edge.V0);
                Assert.Equal(vertices[(i + 1) % vertices.Count], edge.V1);
            }
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxClosedPolylineShape_NoInterior() {
            var vertices = S2TextFormat.ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
            var shape = new S2LaxClosedPolylineShape(vertices.ToArray());
            Assert.Equal(1, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2VertexIdLaxLoopShape_EmptyLoop() {
            var shape = new S2VertexIdLaxLoopShape(System.Array.Empty<int>(), null);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumVertices);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(2, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2VertexIdLaxLoopShape_InvertedLoop() {
            var vertex_array = S2TextFormat.ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
            var vertex_ids = new []{ 0, 3, 2, 1 };  // Inverted.
            var shape = new S2VertexIdLaxLoopShape(vertex_ids, vertex_array.ToArray());
            Assert.Equal(4, shape.NumEdges);
            Assert.Equal(4, shape.NumVertices);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(4, shape.GetChain(0).Length);
            Assert.Equal(vertex_array[0], shape.Vertex(0));
            Assert.Equal(vertex_array[3], shape.Vertex(1));
            Assert.Equal(vertex_array[2], shape.Vertex(2));
            Assert.Equal(vertex_array[1], shape.Vertex(3));
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.True(shape.ContainsBruteForce(S2PointUtil.Origin));
        }
    }
}
