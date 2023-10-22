namespace S2Geometry;

using LabelSetIds = List<System.Int32>;

public class S2BuilderUtil_LaxPolylineLayerTests
{
    private static ITestOutputHelper _logger;

    public S2BuilderUtil_LaxPolylineLayerTests(ITestOutputHelper logger) { _logger = logger; }

    private static void TestS2LaxPolylineShape(
                    List<string> input_strs,
                    string expected_str, EdgeType edge_type,
                    Options? options = null)
    {
        _logger.WriteLine(edge_type == EdgeType.DIRECTED ? "DIRECTED" : "UNDIRECTED");
        S2Builder builder = new(options ?? new Options());
        S2LaxPolylineShape output = new();
        builder.StartLayer(new LaxPolylineLayer(
            output, new LaxPolylineLayer.Options(edge_type)));
        foreach (var input_str in input_strs)
        {
            builder.AddShape(MakeLaxPolylineOrDie(input_str));
        }

        Assert.True(builder.Build(out _));
        Assert.Equal(expected_str, S2TextFormat.ToDebugString(output));
    }

    // Convenience function that tests both directed and undirected edges.
    private static void TestS2LaxPolylineShape(
                List<string> input_strs, string expected_str,
                Options? options = null)
    {
        options ??= new Options();
        TestS2LaxPolylineShape(input_strs, expected_str, EdgeType.DIRECTED, options);
        TestS2LaxPolylineShape(input_strs, expected_str, EdgeType.UNDIRECTED, options);
    }

    private void TestS2LaxPolylineShapeUnchanged(string input_str)
    {
        TestS2LaxPolylineShape([input_str], input_str);
    }

    [Fact]
    internal void LaxPolylineLayer_NoEdges() {
        TestS2LaxPolylineShape([], "");
    }

    [Fact]
    internal void LaxPolylineLayer_OneEdge() {
        // Even with undirected edges, LaxPolylineLayer prefers to reconstruct edges
        // in their original direction.
        TestS2LaxPolylineShapeUnchanged("3:4, 1:1");
        TestS2LaxPolylineShapeUnchanged("1:1, 3:4");
    }

    [Fact]
    internal void LaxPolylineLayer_StraightLineWithBacktracking() {
        TestS2LaxPolylineShapeUnchanged(
            "0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0, 4:0");
    }

    [Fact]
    internal void LaxPolylineLayer_EarlyWalkTerminationWithEndLoop1() {
        // Test that the "early walk termination" code (which is needed by
        // S2LaxPolylineShapeVectorLayer in order to implement idempotency) does not
        // create two polylines when it is possible to assemble the edges into one.
        //
        // This example tests a code path where the early walk termination code
        // should not be triggered at all (but was at one point due to a bug).
        S2Builder.Options options = new()
        {
            SnapFunction = new IntLatLngSnapFunction(2)
        };
        TestS2LaxPolylineShape(["0:0, 0:2, 0:1"], "0:0, 0:1, 0:2, 0:1", options);
    }

    [Fact]
    internal void LaxPolylineLayer_EarlyWalkTerminationWithEndLoop2() {
        // This tests a different code path where the walk is terminated early
        // (yield a polyline with one edge), and then the walk is "maximimzed" by
        // appending a two-edge loop to the end.
        TestS2LaxPolylineShape(["0:0, 0:1", "0:2, 0:1", "0:1, 0:2"],
                         "0:0, 0:1, 0:2, 0:1");
    }

    [Fact]
    internal void LaxPolylineLayer_SimpleLoop() {
        TestS2LaxPolylineShapeUnchanged("0:0, 0:5, 5:5, 5:0, 0:0");
    }

    [Fact]
    internal void LaxPolylineLayer_ManyLoops() {
        // This polyline consists of many overlapping loops that keep returning to
        // the same starting vertex (2:2).  This tests whether the implementation is
        // able to assemble the polyline in the original order.
        TestS2LaxPolylineShapeUnchanged(
            "0:0, 2:2, 2:4, 2:2, 2:4, 4:4, 4:2, 2:2, 4:4, 4:2, 2:2, 2:0, 2:2, "+   
            "2:0, 4:0, 2:2, 4:2, 2:2, 0:2, 0:4, 2:2, 0:4, 0:2, 2:2, 0:4, 2:2, "+    
            "0:2, 2:2, 0:0, 0:2, 2:2, 0:0");
    }

    [Fact]
    internal void LaxPolylineLayer_UnorderedLoops() {
        // This test consists of 5 squares that touch diagonally, similar to the 5
        // white squares of a 3x3 chessboard.  The edges of these squares need to be
        // reordered to assemble them into a single unbroken polyline.
        TestS2LaxPolylineShape([
            "3:3, 3:2, 2:2, 2:3, 3:3",
      "1:0, 0:0, 0:1, 1:1, 1:0",
      "3:1, 3:0, 2:0, 2:1, 3:1",
      "1:3, 1:2, 0:2, 0:1, 1:3",
      "1:1, 1:2, 2:2, 2:1, 1:1",  // Central square
      ],
    "3:3, 3:2, 2:2, 2:1, 3:1, 3:0, 2:0, 2:1, 1:1, 1:0, 0:0, "+
    "0:1, 1:1, 1:2, 0:2, 0:1, 1:3, 1:2, 2:2, 2:3, 3:3");
    }

    [Fact]
    internal void LaxPolylineLayer_SplitEdges() {
        // Test reconstruction of a polyline where two edges have been split into
        // many pieces by crossing edges.  This example is particularly challenging
        // because (1) the edges form a loop, and (2) the first and last edges are
        // identical (but reversed).  This is designed to test the heuristics that
        // attempt to find the first edge of the input polyline.
        S2Builder.Options options = new()
        {
            SplitCrossingEdges = true,
            SnapFunction = new IntLatLngSnapFunction(7)
        };
        TestS2LaxPolylineShape(["0:10, 0:0, 1:0, -1:2, 1:4, -1:6, 1:8, -1:10, -5:0, 0:0, 0:10"],
      "0:10, 0:9, 0:7, 0:5, 0:3, 0:1, 0:0, 1:0, 0:1, -1:2, 0:3, 1:4, 0:5, "+
      "-1:6, 0:7, 1:8, 0:9, -1:10, -5:0, 0:0, 0:1, 0:3, 0:5, 0:7, 0:9, 0:10",
      options);
    }

    [Fact]
    internal void LaxPolylineLayer_SimpleEdgeLabels() {
        S2Builder builder=new(new Options());
        S2LaxPolylineShape output=new();
        LabelSetIds label_set_ids=[];
        IdSetLexicon label_set_lexicon=new();
        builder.StartLayer(new LaxPolylineLayer(
            output, label_set_ids, label_set_lexicon,
            new LaxPolylineLayer.Options(EdgeType.UNDIRECTED)));
        builder.SetLabel(5);
        builder.AddShape(MakeLaxPolylineOrDie("0:0, 0:1, 0:2"));
        builder.PushLabel(7);
        builder.AddShape(MakeLaxPolylineOrDie("0:3, 0:2"));
        builder.ClearLabels();
        builder.AddShape(MakeLaxPolylineOrDie("0:3, 0:4, 0:5"));
        builder.SetLabel(11);
        builder.AddShape(MakeLaxPolylineOrDie("0:6, 0:5"));
        Assert.True(builder.Build(out S2Error error));
        List<List<int>> expected = [new() { 5 }, new() { 5 }, new() { 5, 7 }, new() { }, new() { }, new() { 11 }];
        Assert.Equal(expected.Count, label_set_ids.Count);
        for (int i = 0; i < expected.Count; ++i)
        {
            Assert.Equal(expected[i].Count,
                      label_set_lexicon.IdSet_(label_set_ids[i]).Count);
            int j = 0;
            foreach (var label in label_set_lexicon.IdSet_(label_set_ids[i]))
            {
                Assert.Equal(expected[i][j++], label);
            }
        }
    }

    [Fact]
    internal void LaxPolylineLayer_AntipodalVertices() {
        S2Builder builder=new(new S2Builder.Options());
        S2LaxPolylineShape output=new();
        builder.StartLayer(new LaxPolylineLayer(output));
        builder.AddEdge(new S2Point(1, 0, 0), new S2Point(-1, 0, 0));
        Assert.True(builder.Build(out _));
        Assert.Equal(output.NumVertices(), 2);
        Assert.Equal(output.Vertex(0), new S2Point(1, 0, 0));
        Assert.Equal(output.Vertex(1), new S2Point(-1, 0, 0));
    }

    [Fact]
    internal void IndexedLaxPolylineLayer_AddsShape() {
        S2Builder builder=new(new Options());
        MutableS2ShapeIndex index=[];
        builder.StartLayer(new LaxPolylineLayer.IndexedLaxPolylineLayer(index));
        string polyline_str = "0:0, 0:10";
        builder.AddShape(MakeLaxPolylineOrDie(polyline_str));
        Assert.True(builder.Build(out S2Error error));
        Assert.Equal(1, index.NumShapeIds());
        S2LaxPolylineShape polyline = (S2LaxPolylineShape)index.Shape(0);
        Assert.Equal(polyline_str, S2TextFormat.ToDebugString(polyline));
    }

    [Fact]
    internal void IndexedLaxPolylineLayer_AddsEmptyShape() {
        S2Builder builder=new(new Options());
        MutableS2ShapeIndex index=[];
        builder.StartLayer(new LaxPolylineLayer.IndexedLaxPolylineLayer(index));
        Assert.True(builder.Build(out S2Error error));
        Assert.Equal(0, index.NumShapeIds());
    }
}
