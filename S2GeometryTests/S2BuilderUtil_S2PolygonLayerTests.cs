namespace S2Geometry;

using S2BuilderUtil;
using static S2Builder;

public class S2BuilderUtil_S2PolygonLayerTests
{
    [Fact]
    internal void Test_S2PolygonLayer_Empty()
    {
        TestS2PolygonUnchanged("");
    }

    [Fact]
    internal void Test_S2PolygonLayer_Full()
    {
        TestS2PolygonUnchanged("full");
    }

    [Fact]
    internal void Test_S2PolygonLayer_SmallLoop()
    {
        TestS2PolygonUnchanged("0:0, 0:1, 1:1");
    }

    [Fact]
    internal void Test_S2PolygonLayer_ThreeLoops()
    {
        // The second two loops are nested.
        TestS2PolygonUnchanged("0:1, 1:1, 0:0; " +
                               "3:3, 3:6, 6:6, 6:3; " +
                               "4:4, 4:5, 5:5, 5:4");
    }

    [Fact]
    internal void Test_S2PolygonLayer_PartialLoop()
    {
        TestS2PolygonError(["0:1, 2:3, 4:5"],
                           S2ErrorCode.BUILDER_EDGES_DO_NOT_FORM_LOOPS);
    }

    [Fact]
    internal void Test_S2PolygonLayer_InvalidPolygon()
    {
        TestS2PolygonError(["0:0, 0:10, 10:0, 10:10, 0:0"],
                           S2ErrorCode.LOOP_SELF_INTERSECTION);
    }

    [Fact]
    internal void Test_S2PolygonLayer_DuplicateInputEdges()
    {
        // Check that S2PolygonLayer can assemble polygons even when there are
        // duplicate edges (after sibling pairs are removed), and then report the
        // duplicate edges as an error.
        S2Builder builder = new(new Options());
        S2Polygon output = new();
        S2PolygonLayer.Options options = new()
        {
            Validate = true
        };
        builder.StartLayer(new S2PolygonLayer(output, options));
        builder.AddPolyline(MakePolylineOrDie(
            "0:0, 0:2, 2:2, 1:1, 0:2, 2:2, 2:0, 0:0"));
        Assert.False(builder.Build(out var error));
        Assert.Equal(S2ErrorCode.POLYGON_LOOPS_SHARE_EDGE, error.Code);
        Assert.Equal(2, output.NumLoops());
        S2Loop loop0 = MakeLoopOrDie("0:0, 0:2, 2:2, 2:0");
        S2Loop loop1 = MakeLoopOrDie("0:2, 2:2, 1:1");
        Assert.True(loop0 == output.Loop(0));
        Assert.True(loop1 == output.Loop(1));
    }

    [Fact]
    internal void Test_S2PolygonLayer_DirectedEdgeLabels()
    {
        TestEdgeLabels(EdgeType.DIRECTED);
    }

    [Fact]
    internal void Test_S2PolygonLayer_UndirectedEdgeLabels()
    {
        TestEdgeLabels(EdgeType.UNDIRECTED);
    }

    [Fact]
    internal void Test_S2PolygonLayer_LabelsRequestedButNotProvided()
    {
        // Tests the situation where labels are requested but none were provided.
        S2Builder builder=new(new S2Builder.Options());
        S2Polygon output=new();
        LabelSetIds label_set_ids=[];
        IdSetLexicon label_set_lexicon=new();
        builder.StartLayer(new S2PolygonLayer(
            output, label_set_ids, label_set_lexicon));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:1, 1:0, 0:0"));
        Assert.True(builder.Build(out S2Error error));
        Assert.Equal(label_set_ids.Count, 1);     // One loop.
        Assert.Equal(label_set_ids[0].Count, 3);  // Three edges.
        foreach (var label_set_id in label_set_ids[0])
        {
            Assert.Equal(label_set_id, IdSetLexicon.kEmptySetId);
        }
    }

    [Fact]
    internal void Test_S2PolygonLayer_ThreeLoopsIntoOne()
    {
        // Three loops (two shells and one hole) that combine into one.
        TestS2Polygon([
            "10:0, 0:0, 0:10, 5:10, 10:10, 10:5",
            "0:10, 0:15, 5:15, 5:10",
            "10:10, 5:10, 5:5, 10:5"],
            "10:5, 10:0, 0:0, 0:10, 0:15, 5:15, 5:10, 5:5");
    }

    [Fact]
    internal void Test_S2PolygonLayer_TrianglePyramid()
    {
        // A big CCW triangle containing 3 CW triangular holes.  The whole thing
        // looks like a pyramid of nine triangles.  The output consists of 6
        // positive triangles with no holes.
        TestS2Polygon([
  "0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1",
   "0:2, 1:1, 1:3",
   "0:4, 1:3, 1:5",
   "1:3, 2:2, 2:4"],
            "0:4, 0:6, 1:5; 2:4, 3:3, 2:2; 2:2, 1:1, 1:3; " +
            "1:1, 0:0, 0:2; 1:3, 0:2, 0:4; 1:3, 1:5, 2:4");
    }

    [Fact]
    internal void Test_S2PolygonLayer_ComplexNesting()
    {
        // A complex set of nested polygons, with the loops in random order and the
        // vertices in random cyclic order within each loop.  This test checks that
        // the order (after S2Polygon.InitNested is called) is preserved exactly,
        // whether directed or undirected edges are used.
        TestS2PolygonUnchanged(
            "47:15, 47:5, 5:5, 5:15; " +
            "35:12, 35:7, 27:7, 27:12; " +
            "1:50, 50:50, 50:1, 1:1; " +
            "42:22, 10:22, 10:25, 42:25; " +
            "47:30, 47:17, 5:17, 5:30; " +
            "7:27, 45:27, 45:20, 7:20; " +
            "37:7, 37:12, 45:12, 45:7; " +
            "47:47, 47:32, 5:32, 5:47; " +
            "50:60, 50:55, 1:55, 1:60; " +
            "25:7, 17:7, 17:12, 25:12; " +
            "7:7, 7:12, 15:12, 15:7");
    }

    [Fact]
    internal void Test_S2PolygonLayer_FiveLoopsTouchingAtOneCommonPoint()
    {
        // Five nested loops that touch at one common point.
        TestS2PolygonUnchanged("0:0, 0:10, 10:10, 10:0; " +
                               "0:0, 1:9, 9:9, 9:1; " +
                               "0:0, 2:8, 8:8, 8:2; " +
                               "0:0, 3:7, 7:7, 7:3; " +
                               "0:0, 4:6, 6:6, 6:4");
    }

    [Fact]
    internal void Test_S2PolygonLayer_FourNestedDiamondsTouchingAtTwoPointsPerPair()
    {
        // Four diamonds nested inside each other, where each diamond shares two
        // vertices with the diamond inside it and shares its other two vertices
        // with the diamond that contains it.  The resulting shape looks vaguely
        // like an eye made out of chevrons.
        TestS2Polygon([
  "0:10, -10:0, 0:-10, 10:0",
   "0:-20, -10:0, 0:20, 10:0",
   "0:-10, -5:0, 0:10, 5:0",
   "0:5, -5:0, 0:-5, 5:0"],
            "10:0, 0:10, -10:0, 0:20; " +
            "0:-20, -10:0, 0:-10, 10:0; " +
            "5:0, 0:-10, -5:0, 0:-5; " +
            "0:5, -5:0, 0:10, 5:0");
    }

    [Fact]
    internal void Test_S2PolygonLayer_SevenDiamondsTouchingAtOnePointPerPair()
    {
        // Seven diamonds nested within each other touching at one
        // point between each nested pair.
        TestS2PolygonUnchanged("0:-70, -70:0, 0:70, 70:0; " +
                               "0:-70, -60:0, 0:60, 60:0; " +
                               "0:-50, -60:0, 0:50, 50:0; " +
                               "0:-40, -40:0, 0:50, 40:0; " +
                               "0:-30, -30:0, 0:30, 40:0; " +
                               "0:-20, -20:0, 0:30, 20:0; " +
                               "0:-10, -20:0, 0:10, 10:0");
    }

    [Fact]
    internal void Test_IndexedS2PolygonLayer_AddsShape()
    {
        S2Builder builder = new(new Options());
        MutableS2ShapeIndex index = [];
        builder.StartLayer(new IndexedS2PolygonLayer(index));
        string polygon_str = "0:0, 0:10, 10:0";
        builder.AddPolygon(MakePolygonOrDie(polygon_str));
        Assert.True(builder.Build(out _));
        Assert.Equal(1, index.NumShapeIds());
        S2Polygon polygon = ((S2Polygon.Shape)index.Shape(0)!).Polygon;
        Assert.Equal(polygon_str, polygon.ToDebugString());
    }

    [Fact]
    internal void Test_IndexedS2PolygonLayer_IgnoresEmptyShape()
    {
        S2Builder builder = new(new Options());
        MutableS2ShapeIndex index = [];
        builder.StartLayer(new IndexedS2PolygonLayer(index));
        Assert.True(builder.Build(out var error));
        Assert.Equal(0, index.NumShapeIds());
    }

    private static void TestS2Polygon(string[] input_strs, string expected_str, EdgeType edge_type)
    {
        S2Builder builder = new(new Options());
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(
            output, new S2PolygonLayer.Options(edge_type)));
        bool is_full = false;
        foreach (var input_str in input_strs)
        {
            if (input_str == "full") is_full = true;
            builder.AddPolygon(MakeVerbatimPolygonOrDie(input_str));
        }
        builder.AddIsFullPolygonPredicate(IsFullPolygon(is_full));
        Assert.True(builder.Build(out _));
        // The input strings in tests may not be in normalized form, so we build an
        // S2Polygon and convert it back to a string.
        S2Polygon expected = MakePolygonOrDie(expected_str);
        Assert.Equal(expected.ToDebugString(),
                  output.ToDebugString());
    }

    private static void TestS2Polygon(string[] input_strs, string expected_str)
    {
        TestS2Polygon(input_strs, expected_str, EdgeType.DIRECTED);
        TestS2Polygon(input_strs, expected_str, EdgeType.UNDIRECTED);
    }

    private static void TestS2PolygonUnchanged(string input_str)
    {
        TestS2Polygon([input_str], input_str);
    }

    // Unlike the methods above, the input consists of a set of *polylines*.
    private static void TestS2PolygonError(string[] input_strs, S2ErrorCode expected_error, EdgeType edge_type)
    {
        S2Builder builder = new(new Options());
        S2Polygon output = new();
        S2PolygonLayer.Options options = new(edge_type)
        {
            Validate = true
        };
        builder.StartLayer(new S2PolygonLayer(output, options));
        foreach (var input_str in input_strs)
        {
            builder.AddPolyline(MakePolylineOrDie(input_str));
        }
        Assert.False(builder.Build(out var error));
        Assert.Equal(expected_error, error.Code);
    }

    private static void TestS2PolygonError(string[] input_strs, S2ErrorCode expected_error)
    {
        TestS2PolygonError(input_strs, expected_error, EdgeType.DIRECTED);
        TestS2PolygonError(input_strs, expected_error, EdgeType.UNDIRECTED);
    }


    private static void AddPolylineWithLabels(S2Polyline polyline, EdgeType edge_type,
        Int32 label_begin, S2Builder builder, EdgeLabelMap edge_label_map)
    {
        for (int i = 0; i + 1 < polyline.NumVertices(); ++i)
        {
            Int32 label = label_begin + i;
            builder.SetLabel(label);
            // With undirected edges, reverse the direction of every other input edge.
            int dir = edge_type == EdgeType.DIRECTED ? 1 : (i & 1);
            builder.AddEdge(polyline.Vertex(i + (1 - dir)), polyline.Vertex(i + dir));
            S2Point key = polyline.Vertex(i) + polyline.Vertex(i + 1);
            edge_label_map[key].Add(label);
        }
    }

    private static void TestEdgeLabels(EdgeType edge_type)
    {
        S2Builder builder = new(new Options());
        S2Polygon output = new();
        List<LabelSet> label_set_ids = [];
        IdSetLexicon label_set_lexicon = new();
        builder.StartLayer(new S2PolygonLayer(
            output, label_set_ids, label_set_lexicon,
            new S2PolygonLayer.Options(edge_type)));

        // We use a polygon consisting of 3 loops.  The loops are reordered and
        // some of the loops are inverted during S2Polygon construction.
        EdgeLabelMap edge_label_map = [];
        AddPolylineWithLabels(MakePolylineOrDie(
            "0:0, 9:1, 1:9, 0:0, 2:8, 8:2, 0:0, 0:10, 10:10, 10:0, 0:0"),
            edge_type, 0, builder, edge_label_map);
        Assert.True(builder.Build(out var error));
        int[] expected_loop_sizes = [4, 3, 3];
        Assert.Equal(expected_loop_sizes.Length, label_set_ids.Count);
        for (int i = 0; i < expected_loop_sizes.Length; ++i)
        {
            Assert.Equal(expected_loop_sizes[i], label_set_ids[i].Count);
            for (int j = 0; j < label_set_ids[i].Count; ++j)
            {
                S2Point key = output.Loop(i).Vertex(j) + output.Loop(i).Vertex(j + 1);
                var expected_labels = edge_label_map[key];
                Assert.Equal(expected_labels.Count,
                          label_set_lexicon.IdSet_(label_set_ids[i][j]).Count);
                Assert.Equal(expected_labels,
                    label_set_lexicon.IdSet_(label_set_ids[i][j]));
            }
        }
    }
}
