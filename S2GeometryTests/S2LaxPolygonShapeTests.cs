using S2Geometry.S2ShapeUtil;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace S2Geometry
{
    public class S2LaxPolygonShapeTests
    {
        [Fact]
        public void Test_S2LaxPolygonShape_EmptyPolygon() {
            var shape = new S2LaxPolygonShape(new S2Polygon());
            Assert.Equal(0, shape.NumLoops);
            Assert.Equal(0, shape.NumVertices);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(0, shape.NumChains());
            Assert.Equal(2, shape.Dimension());
            Assert.True(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxPolygonShape_FullPolygon() {
            var shape = new S2LaxPolygonShape(new S2Polygon(S2TextFormat.MakeLoopOrDie("full")));
            Assert.Equal(1, shape.NumLoops);
            Assert.Equal(0, shape.NumVertices);
            Assert.Equal(0, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.True(shape.IsFull);
            Assert.True(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxPolygonShape_SingleVertexPolygon() {
            // S2Polygon doesn't support single-vertex loops, so we need to construct
            // the S2LaxPolygonShape directly.
            var loops = new List<List<S2Point>>
            {
                S2TextFormat.ParsePointsOrDie("0:0")
            };
            var shape = new S2LaxPolygonShape(loops);
            Assert.Equal(1, shape.NumLoops);
            Assert.Equal(1, shape.NumVertices);
            Assert.Equal(1, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(1, shape.GetChain(0).Length);
            var edge = shape.GetEdge(0);
            Assert.Equal(loops[0][0], edge.V0);
            Assert.Equal(loops[0][0], edge.V1);
            Assert.True(edge == shape.ChainEdge(0, 0));
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxPolygonShape_SingleLoopPolygon() {
            // Test S2Polygon constructor.
            var vertices = S2TextFormat.ParsePointsOrDie("0:0, 0:1, 1:1, 1:0");
            var shape = new S2LaxPolygonShape(new S2Polygon(new S2Loop(vertices)));
            Assert.Equal(1, shape.NumLoops);
            Assert.Equal(vertices.Count, shape.NumVertices);
            Assert.Equal(vertices.Count, shape.NumLoopVertices(0));
            Assert.Equal(vertices.Count, shape.NumEdges);
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(vertices.Count, shape.GetChain(0).Length);
            for (int i = 0; i < vertices.Count; ++i) {
                Assert.Equal(vertices[i], shape.LoopVertex(0, i).First());
                var edge = shape.GetEdge(i);
                Assert.Equal(vertices[i], edge.V0);
                Assert.Equal(vertices[(i + 1) % vertices.Count], edge.V1);
                Assert.Equal(edge.V0, shape.ChainEdge(0, i).V0);
                Assert.Equal(edge.V1, shape.ChainEdge(0, i).V1);
            }
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.ContainsBruteForce(S2PointUtil.Origin));
        }

        [Fact]
        public void Test_S2LaxPolygonShape_MultiLoopPolygon() {
            // Test S2Point[][] constructor.  Make sure that the loops are
            // oriented so that the interior of the polygon is always on the left.
            var loops = new List<List<S2Point>>{
    S2TextFormat.ParsePointsOrDie("0:0, 0:3, 3:3"),  // CCW
    S2TextFormat.ParsePointsOrDie("1:1, 2:2, 1:2")   // CW
  };
            var shape = new S2LaxPolygonShape(loops);

            Assert.Equal(loops.Count, shape.NumLoops);
            int num_vertices = 0;
            Assert.Equal(loops.Count, shape.NumChains());
            for (int i = 0; i < loops.Count; ++i) {
                Assert.Equal(loops[i].Count, shape.NumLoopVertices(i));
                Assert.Equal(num_vertices, shape.GetChain(i).Start);
                Assert.Equal(loops[i].Count, shape.GetChain(i).Length);
                for (int j = 0; j < loops[i].Count; ++j) {
                    Assert.Equal(loops[i][j], shape.LoopVertex(i, j).First());
                    var edge = shape.GetEdge(num_vertices + j);
                    Assert.Equal(loops[i][j], edge.V0);
                    Assert.Equal(loops[i][(j + 1) % loops[i].Count], edge.V1);
                }
                num_vertices += loops[i].Count;
            }
            Assert.Equal(num_vertices, shape.NumVertices);
            Assert.Equal(num_vertices, shape.NumEdges);
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty);
            Assert.False(shape.IsFull);
            Assert.False(shape.ContainsBruteForce(S2PointUtil.Origin));
        }

        [Fact]
        public void Test_S2LaxPolygonShape_MultiLoopS2Polygon() {
            // Verify that the orientation of loops representing holes is reversed when
            // converting from an S2Polygon to an S2LaxPolygonShape.
            var polygon = S2TextFormat.MakePolygonOrDie("0:0, 0:3, 3:3; 1:1, 1:2, 2:2");
            var shape = new S2LaxPolygonShape(polygon);
            for (int i = 0; i < polygon.NumLoops(); ++i) {
                var loop = polygon.Loop(i);
                for (int j = 0; j < loop.NumVertices; ++j) {
                    Assert.Equal(loop.OrientedVertex(j),
                              shape.LoopVertex(i, j).First());
                }
            }
        }

        [Fact]
        public void Test_S2LaxPolygonShape_ManyLoopPolygon() {
            // Test a polygon with enough loops so that cumulative_vertices_ is used.
            var loops = new List<List<S2Point>>();
            for (int i = 0; i < 100; ++i) {
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
            for (int i = 0; i < loops.Count; ++i) {
                Assert.Equal(loops[i].Count, shape.NumLoopVertices(i));
                Assert.Equal(num_vertices, shape.GetChain(i).Start);
                Assert.Equal(loops[i].Count, shape.GetChain(i).Length);
                for (int j = 0; j < loops[i].Count; ++j) {
                    Assert.Equal(loops[i][j], shape.LoopVertex(i, j).First());
                    var edge = shape.GetEdge(num_vertices + j);
                    Assert.Equal(loops[i][j], edge.V0);
                    Assert.Equal(loops[i][(j + 1) % loops[i].Count], edge.V1);
                }
                num_vertices += loops[i].Count;
            }
            Assert.Equal(num_vertices, shape.NumVertices);
            Assert.Equal(num_vertices, shape.NumEdges);
        }

        [Fact]
        public void Test_S2LaxPolygonShape_DegenerateLoops() {
            var loops = new List<List<S2Point>>{
    S2TextFormat.ParsePointsOrDie("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
    S2TextFormat.ParsePointsOrDie("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
    S2TextFormat.ParsePointsOrDie("5:5, 6:6")
  };
            var shape = new S2LaxPolygonShape(loops);
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2LaxPolygonShape_InvertedLoops() {
            var loops = new List<List<S2Point>>{
    S2TextFormat.ParsePointsOrDie("1:2, 1:1, 2:2"),
    S2TextFormat.ParsePointsOrDie("3:4, 3:3, 4:4")
  };
            var shape = new S2LaxPolygonShape(loops);
            Assert.True(shape.ContainsBruteForce(S2PointUtil.Origin));
        }

        [Fact]
        public void Test_S2LaxPolygonShape_CompareToS2Loop()
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

        private static void CompareS2LoopToShape(S2Loop loop, S2Shape shape)
        {
            MutableS2ShapeIndex index = new();
            index.Add(shape);
            S2Cap cap = loop.GetCapBound();
            var query = index.MakeS2ContainsPointQuery();
            for (int iter = 0; iter < 100; ++iter) {
                S2Point point = S2Testing.SamplePoint(cap);
                Assert.Equal(loop.Contains(point),
                    query.ShapeContains(index.Shape(0), point));
            }
        }
    }
}
