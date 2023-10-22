namespace S2Geometry;

public class S2ShapeUtilBuildPolygonBoundariesTests
{
    [Fact]
    internal void Test_BuildPolygonBoundaries_NoComponents()
    {
        List<List<S2Shape>> faces = [];
        List<List<S2Shape>> components = [];
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Empty(faces);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_OneLoop()
    {
        var a0 = BuildTestLaxLoop("0:0, 1:0, 0:1");  // Outer face
        var a1 = BuildTestLaxLoop("0:0, 0:1, 1:0");
        var faces = new List<List<S2Shape>>();
        var components = new List<List<S2Shape>> { new() { a0, a1 } };
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Equal(2, faces.Count);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_TwoLoopsSameComponent()
    {
        var a0 = BuildTestLaxLoop("0:0, 1:0, 0:1");  // Outer face
        var a1 = BuildTestLaxLoop("0:0, 0:1, 1:0");
        var a2 = BuildTestLaxLoop("1:0, 0:1, 1:1");
        var faces = new List<List<S2Shape>>();
        var components = new List<List<S2Shape>> { new() { a0, a1, a2 } };
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Equal(3, faces.Count);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_TwoNestedLoops()
    {
        var a0 = BuildTestLaxLoop("0:0, 3:0, 0:3");  // Outer face
        var a1 = BuildTestLaxLoop("0:0, 0:3, 3:0");
        var b0 = BuildTestLaxLoop("1:1, 2:0, 0:2");  // Outer face
        var b1 = BuildTestLaxLoop("1:1, 0:2, 2:0");
        var faces = new List<List<S2Shape>>();
        var components = new List<List<S2Shape>> { new() { a0, a1 }, new() { b0, b1 } };
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Equal(3, faces.Count);
        Assert.Equal(new S2Shape[] { b0, a1 }, faces[0]);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_TwoLoopsDifferentComponents()
    {
        var a0 = BuildTestLaxLoop("0:0, 1:0, 0:1");  // Outer face
        var a1 = BuildTestLaxLoop("0:0, 0:1, 1:0");
        var b0 = BuildTestLaxLoop("0:2, 1:2, 0:3");  // Outer face
        var b1 = BuildTestLaxLoop("0:2, 0:3, 1:2");
        var faces = new List<List<S2Shape>>();
        var components = new List<List<S2Shape>> { new() { a0, a1 }, new() { b0, b1 } };
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Equal(3, faces.Count);
        Assert.Equal(new S2Shape[] { a0, b0 }, faces[2]);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_OneDegenerateLoop()
    {
        var a0 = BuildTestLaxLoop("0:0, 1:0, 0:0");
        var faces = new List<List<S2Shape>>();
        var components = new List<List<S2Shape>> { new() { a0 } };
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Single(faces);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_TwoDegenerateLoops()
    {
        var a0 = BuildTestLaxLoop("0:0, 1:0, 0:0");
        var b0 = BuildTestLaxLoop("2:0, 3:0, 2:0");
        var faces = new List<List<S2Shape>>();
        var components = new List<List<S2Shape>> { new() { a0, b0 } };
        S2ShapeUtil.BuildPolygonBoundaries(components, faces);
        Assert.Single(faces);
        Assert.Equal(2, faces[0].Count);
    }

    [Fact]
    internal void Test_BuildPolygonBoundaries_ComplexTest1()
    {
        // Loops at index 0 are the outer (clockwise) loops.
        // Component "a" consists of 4 adjacent squares forming a larger square.
        var a0 = BuildTestLaxLoop("0:0, 25:0, 50:0, 50:25, 50:50, 25:50, 0:50, 0:50");
        var a1 = BuildTestLaxLoop("0:0, 0:25, 25:25, 25:0");
        var a2 = BuildTestLaxLoop("0:25, 0:50, 25:50, 25:25");
        var a3 = BuildTestLaxLoop("25:0, 25:25, 50:25, 50:0");
        var a4 = BuildTestLaxLoop("25:25, 25:50, 50:50, 50:25");
        // Component "b" consists of a degenerate loop to the left of "a".
        var b0 = BuildTestLaxLoop("0:-10, 10:-10");
        // Components "a1_a", "a1_b", and "a1_c" are located within "a1".
        var a1_a0 = BuildTestLaxLoop("5:5, 20:5, 20:10, 5:10");
        var a1_a1 = BuildTestLaxLoop("5:5, 5:10, 10:10, 10:5");
        var a1_a2 = BuildTestLaxLoop("10:5, 10:10, 15:10, 15:5");
        var a1_a3 = BuildTestLaxLoop("15:5, 15:10, 20:10, 20:5");
        var a1_b0 = BuildTestLaxLoop("5:15, 20:15, 20:20, 5:20");
        var a1_b1 = BuildTestLaxLoop("5:15, 5:20, 20:20, 20:15");
        var a1_c0 = BuildTestLaxLoop("2:5, 2:10, 2:5");
        // Two components located inside "a1_a2" and "a1_a3".
        var a1_a2_a0 = BuildTestLaxLoop("11:6, 14:6, 14:9, 11:9");
        var a1_a2_a1 = BuildTestLaxLoop("11:6, 11:9, 14:9, 14:6");
        var a1_a3_a0 = BuildTestLaxLoop("16:6, 19:9, 16:6");
        // Five component located inside "a3" and "a4".
        var a3_a0 = BuildTestLaxLoop("30:5, 45:5, 45:20, 30:20");
        var a3_a1 = BuildTestLaxLoop("30:5, 30:20, 45:20, 45:5");
        var a4_a0 = BuildTestLaxLoop("30:30, 40:30, 30:30");
        var a4_b0 = BuildTestLaxLoop("30:35, 40:35, 30:35");
        var a4_c0 = BuildTestLaxLoop("30:40, 40:40, 30:40");
        var a4_d0 = BuildTestLaxLoop("30:45, 40:45, 30:45");
        var components = new List<List<S2Shape>>
            {
                new() {a0, a1, a2, a3, a4},
                new() {b0},
                new() {a1_a0, a1_a1, a1_a2, a1_a3},
                new() {a1_b0, a1_b1},
                new() {a1_c0},
                new() {a1_a2_a0, a1_a2_a1},
                new() {a1_a3_a0},
                new() {a3_a0, a3_a1},
                new() {a4_a0},
                new() {a4_b0},
                new() {a4_c0},
                new() {a4_d0}
            };
        var expected_faces = new List<List<S2Shape>>
            {
                new() {a0, b0},
                new() {a1, a1_a0, a1_b0, a1_c0},
                new() {a1_a1},
                new() {a1_a2, a1_a2_a0},
                new() {a1_a2_a1},
                new() {a1_a3, a1_a3_a0},
                new() {a1_b1},
                new() {a2},
                new() {a3, a3_a0},
                new() {a3_a1},
                new() {a4, a4_a0, a4_b0, a4_c0, a4_d0}
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

    private static S2LaxLoopShape BuildTestLaxLoop(string vertex_str)
    {
        var vertices = ParsePointsOrDie(vertex_str);
        return new S2LaxLoopShape([.. vertices]);
    }
}
