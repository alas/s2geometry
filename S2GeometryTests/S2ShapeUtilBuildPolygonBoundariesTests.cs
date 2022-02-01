using System.Collections.Generic;
using Xunit;

namespace S2Geometry
{
    public class S2ShapeUtilBuildPolygonBoundariesTests
    {
        [Fact]
        public void Test_BuildPolygonBoundaries_NoComponents()
        {
            List<List<S2Shape>> faces = new ();
            List<List<S2Shape>> components = new ();
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Empty(faces);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_OneLoop()
        {
            TestLaxLoop a0 = new("0:0, 1:0, 0:1");  // Outer face
            TestLaxLoop a1 = new("0:0, 0:1, 1:0");
            var faces = new List<List<S2Shape>>();
            var components = new List<List<S2Shape>> { new List<S2Shape> { a0, a1 } };
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Equal(2, faces.Count);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_TwoLoopsSameComponent()
        {
            TestLaxLoop a0 = new("0:0, 1:0, 0:1");  // Outer face
            TestLaxLoop a1 = new("0:0, 0:1, 1:0");
            TestLaxLoop a2 = new("1:0, 0:1, 1:1");
            var faces = new List<List<S2Shape>>();
            var components = new List<List<S2Shape>> { new List<S2Shape> { a0, a1, a2 } };
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Equal(3, faces.Count);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_TwoNestedLoops()
        {
            TestLaxLoop a0 = new("0:0, 3:0, 0:3");  // Outer face
            TestLaxLoop a1 = new("0:0, 0:3, 3:0");
            TestLaxLoop b0 = new("1:1, 2:0, 0:2");  // Outer face
            TestLaxLoop b1 = new("1:1, 0:2, 2:0");
            var faces = new List<List<S2Shape>>();
            var components = new List<List<S2Shape>> { new List<S2Shape> { a0, a1 }, new List<S2Shape> { b0, b1 } };
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Equal(3, faces.Count);
            Assert.Equal((new S2Shape[] { b0, a1 }), faces[0]);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_TwoLoopsDifferentComponents()
        {
            TestLaxLoop a0 = new("0:0, 1:0, 0:1");  // Outer face
            TestLaxLoop a1 = new("0:0, 0:1, 1:0");
            TestLaxLoop b0 = new("0:2, 1:2, 0:3");  // Outer face
            TestLaxLoop b1 = new("0:2, 0:3, 1:2");
            var faces = new List<List<S2Shape>>();
            var components = new List<List<S2Shape>> { new List<S2Shape> { a0, a1 }, new List<S2Shape> { b0, b1 } };
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Equal(3, faces.Count);
            Assert.Equal((new S2Shape[] { a0, b0 }), faces[2]);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_OneDegenerateLoop()
        {
            TestLaxLoop a0 = new("0:0, 1:0, 0:0");
            var faces = new List<List<S2Shape>>();
            var components = new List<List<S2Shape>> { new List<S2Shape> { a0 } };
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Single(faces);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_TwoDegenerateLoops()
        {
            TestLaxLoop a0 = new("0:0, 1:0, 0:0");
            TestLaxLoop b0 = new("2:0, 3:0, 2:0");
            var faces = new List<List<S2Shape>>();
            var components = new List<List<S2Shape>> { new List<S2Shape> { a0, b0 } };
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Single(faces);
            Assert.Equal(2, faces[0].Count);
        }

        [Fact]
        public void Test_BuildPolygonBoundaries_ComplexTest1()
        {
            // Loops at index 0 are the outer (clockwise) loops.
            // Component "a" consists of 4 adjacent squares forming a larger square.
            TestLaxLoop a0 = new("0:0, 25:0, 50:0, 50:25, 50:50, 25:50, 0:50, 0:50");
            TestLaxLoop a1 = new("0:0, 0:25, 25:25, 25:0");
            TestLaxLoop a2 = new("0:25, 0:50, 25:50, 25:25");
            TestLaxLoop a3 = new("25:0, 25:25, 50:25, 50:0");
            TestLaxLoop a4 = new("25:25, 25:50, 50:50, 50:25");
            // Component "b" consists of a degenerate loop to the left of "a".
            TestLaxLoop b0 = new("0:-10, 10:-10");
            // Components "a1_a", "a1_b", and "a1_c" are located within "a1".
            TestLaxLoop a1_a0 = new("5:5, 20:5, 20:10, 5:10");
            TestLaxLoop a1_a1 = new("5:5, 5:10, 10:10, 10:5");
            TestLaxLoop a1_a2 = new("10:5, 10:10, 15:10, 15:5");
            TestLaxLoop a1_a3 = new("15:5, 15:10, 20:10, 20:5");
            TestLaxLoop a1_b0 = new("5:15, 20:15, 20:20, 5:20");
            TestLaxLoop a1_b1 = new("5:15, 5:20, 20:20, 20:15");
            TestLaxLoop a1_c0 = new("2:5, 2:10, 2:5");
            // Two components located inside "a1_a2" and "a1_a3".
            TestLaxLoop a1_a2_a0 = new("11:6, 14:6, 14:9, 11:9");
            TestLaxLoop a1_a2_a1 = new("11:6, 11:9, 14:9, 14:6");
            TestLaxLoop a1_a3_a0 = new("16:6, 19:9, 16:6");
            // Five component located inside "a3" and "a4".
            TestLaxLoop a3_a0 = new("30:5, 45:5, 45:20, 30:20");
            TestLaxLoop a3_a1 = new("30:5, 30:20, 45:20, 45:5");
            TestLaxLoop a4_a0 = new("30:30, 40:30, 30:30");
            TestLaxLoop a4_b0 = new("30:35, 40:35, 30:35");
            TestLaxLoop a4_c0 = new("30:40, 40:40, 30:40");
            TestLaxLoop a4_d0 = new("30:45, 40:45, 30:45");
            var components = new List<List<S2Shape>>
                {
                    new List<S2Shape>{a0, a1, a2, a3, a4},
                    new List<S2Shape>{b0},
                    new List<S2Shape>{a1_a0, a1_a1, a1_a2, a1_a3},
                    new List<S2Shape>{a1_b0, a1_b1},
                    new List<S2Shape>{a1_c0},
                    new List<S2Shape>{a1_a2_a0, a1_a2_a1},
                    new List<S2Shape>{a1_a3_a0},
                    new List<S2Shape>{a3_a0, a3_a1},
                    new List<S2Shape>{a4_a0},
                    new List<S2Shape>{a4_b0},
                    new List<S2Shape>{a4_c0},
                    new List<S2Shape>{a4_d0}
                };
            var expected_faces = new List<List<S2Shape>>
                {
                    new List<S2Shape>{a0, b0},
                    new List<S2Shape>{a1, a1_a0, a1_b0, a1_c0},
                    new List<S2Shape>{a1_a1},
                    new List<S2Shape>{a1_a2, a1_a2_a0},
                    new List<S2Shape>{a1_a2_a1},
                    new List<S2Shape>{a1_a3, a1_a3_a0},
                    new List<S2Shape>{a1_b1},
                    new List<S2Shape>{a2},
                    new List<S2Shape>{a3, a3_a0},
                    new List<S2Shape>{a3_a1},
                    new List<S2Shape>{a4, a4_a0, a4_b0, a4_c0, a4_d0}
                };
            var faces = new List<List<S2Shape>>();
            S2ShapeUtil.BuildPolygonBoundaries(components, faces);
            Assert.Equal(expected_faces.Count, faces.Count);
            SortFaces(expected_faces);
            SortFaces(faces);
            Assert.Equal(expected_faces, faces);
        }

        static void SortFaces(List<List<S2Shape>> faces)
        {
            foreach (var face in faces)
            {
                face.Sort();
            }
            faces.Sort();
        }

        public class TestLaxLoop : S2LaxLoopShape
        {
            public TestLaxLoop(string vertex_str)
            {
                var vertices = ParsePointsOrDie(vertex_str);
                Init(vertices.ToArray());
            }
        }
    }
}
