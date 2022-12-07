namespace S2Geometry;

using S2BuilderUtil;
using static S2Builder;
using static S2Builder.GraphOptions;

public class S2BuilderUtil_ClosedSetNormalizerTests
{
    private readonly NormalizeTest normalizeTest = new();

    [Fact]
    internal void Test_NormalizeTest_EmptyGraphs()
    {
        normalizeTest.Run("# #", "# #");
    }

    [Fact]
    internal void Test_NormalizeTest_NonDegenerateInputs()
    {
        normalizeTest.Run("0:0 # 1:0, 1:1 | 1:2, 1:3 # 2:2, 2:3, 3:2",
            "0:0 # 1:0, 1:1 | 1:2, 1:3 # 2:2, 2:3, 3:2");
    }

    [Fact]
    internal void Test_NormalizeTest_PointShell()
    {
        normalizeTest.Run("# # 0:0", "0:0 # #");
    }

    [Fact]
    internal void Test_NormalizeTest_PointHole()
    {
        normalizeTest.Run("# # 0:0, 0:3, 3:0 | 1:1", "# # 0:0, 0:3, 3:0");
    }

    [Fact]
    internal void Test_NormalizeTest_PointPolyline()
    {
        // Verify that a single degenerate polyline edge is transformed into a
        // single point.  Note that since the polyline layer is undirected while the
        // point layer is not, this tests that the edge count is halved when the
        // edge is demoted.
        normalizeTest.Run("# 0:0, 0:0 #", "0:0 # #");
    }

    [Fact]
    internal void Test_NormalizeTest_SiblingPairShell()
    {
        normalizeTest.Run("# # 0:0, 1:0 ", "# 0:0, 1:0 #");
    }

    [Fact]
    internal void Test_NormalizeTest_SiblingPairHole()
    {
        normalizeTest.Run("# # 0:0, 0:3, 3:0; 0:0, 1:1",
            "# # 0:0, 0:3, 3:0");
    }

    [Fact]
    internal void Test_NormalizeTest_PointSuppressedByPolygonVertex()
    {
        normalizeTest.Run("0:0 | 0:1 | 1:0 # # 0:0, 0:1, 1:0",
            "# # 0:0, 0:1, 1:0");
        normalizeTest.SuppressLowerDimensions_ = false;
        normalizeTest.Run("0:0 | 0:1 | 1:0 # # 0:0, 0:1, 1:0",
            "0:0 | 0:1 | 1:0 # # 0:0, 0:1, 1:0");
    }

    [Fact]
    internal void Test_NormalizeTest_PointSuppressedByPolylineVertex()
    {
        normalizeTest.Run("0:0 | 0:1 # 0:0, 0:1 #", "# 0:0, 0:1 #");
        normalizeTest.SuppressLowerDimensions_ = false;
        normalizeTest.Run("0:0 | 0:1 # 0:0, 0:1 #", "0:0 | 0:1 # 0:0, 0:1 #");
    }

    [Fact]
    internal void Test_NormalizeTest_PointShellSuppressedByPolylineEdge()
    {
        // This tests the case where a single-point shell is demoted to a point
        // which is then suppressed by a matching polyline vertex.
        normalizeTest.Run("# 0:0, 1:0 # 0:0; 1:0", "# 0:0, 1:0 #");
        normalizeTest.SuppressLowerDimensions_ = false;
        normalizeTest.Run("# 0:0, 1:0 # 0:0; 1:0", "0:0 | 1:0 # 0:0, 1:0 #");
    }

    [Fact]
    internal void Test_NormalizeTest_PolylineEdgeSuppressedByPolygonEdge()
    {
        normalizeTest.Run("# 0:0, 0:1 # 0:0, 0:1, 1:0", "# # 0:0, 0:1, 1:0");
        normalizeTest.SuppressLowerDimensions_ = false;
        normalizeTest.Run("# 0:0, 0:1 # 0:0, 0:1, 1:0", "# 0:0, 0:1 # 0:0, 0:1, 1:0");
    }

    [Fact]
    internal void Test_NormalizeTest_PolylineEdgeSuppressedByReversePolygonEdge()
    {
        normalizeTest.GraphOptionsOut_[1].EdgeType_ = (EdgeType.DIRECTED);
        normalizeTest.Run("# 1:0, 0:0 # 0:0, 0:1, 1:0", "# # 0:0, 0:1, 1:0");
        normalizeTest.SuppressLowerDimensions_ = false;
        normalizeTest.Run("# 1:0, 0:0 # 0:0, 0:1, 1:0", "# 1:0, 0:0 # 0:0, 0:1, 1:0");
    }

    [Fact]
    internal void Test_NormalizeTest_DuplicateEdgeMerging()
    {
        // Verify that when DuplicateEdges.KEEP is specified, demoted edges are
        // added as new edges rather than being merged with existing ones.
        // (Note that NormalizeTest specifies DuplicateEdges.KEEP by default.)
        normalizeTest.Run("0:0 | 0:0 # 0:0, 0:0 | 0:1, 0:2 # 0:0; 0:1, 0:2",
            "0:0 | 0:0 | 0:0 | 0:0 # 0:1, 0:2 | 0:1, 0:2 #");
        // Now verify that the duplicate edges are merged if requested.
        normalizeTest.GraphOptionsOut_[0].DuplicateEdges_ = (DuplicateEdges.MERGE);
        normalizeTest.GraphOptionsOut_[1].DuplicateEdges_ = (DuplicateEdges.MERGE);
        normalizeTest.Run("0:0 | 0:0 # 0:0, 0:0 | 0:1, 0:2 # 0:0; 0:1, 0:2",
            "0:0 # 0:1, 0:2 #");
    }

    [Fact]
    internal void Test_ComputeUnion_MixedGeometry()
    {
        // Verifies that the code above works.  Features tested include:
        //  - Points and polylines in the interior of the other polygon are removed
        //  - Degenerate polygon shells are converted to points/polylines
        //  - Degenerate polygon holes are removed
        //  - Points coincident with polyline or polygon edges are removed
        //  - Polyline edges coincident with polygon edges are removed
        var a = MakeIndexOrDie(
            "0:0 | 10:10 | 20:20 # " +
            "0:0, 0:10 | 0:0, 10:0 | 15:15, 16:16 # " +
            "0:0, 0:10, 10:10, 10:0; 0:0, 1:1; 2:2; 10:10, 11:11; 12:12");
        var b = MakeIndexOrDie(
            "0:10 | 10:0 | 3:3 | 16:16 # " +
            "10:10, 0:10 | 10:10, 10:0 | 5:5, 6:6 # " +
            "19:19, 19:21, 21:21, 21:19");
        MutableS2ShapeIndex result = new();
        Assert.True(ComputeUnion(a, b, result, out _));
        Assert.Equal("12:12 # " +
                  "15:15, 16:16 | 10:10, 11:11 # " +
                  "0:0, 0:10, 10:10, 10:0; 19:19, 19:21, 21:21, 21:19",
                  result.ToDebugString());
    }


    // A test harness that sets default values for ClosedSetNormalizer.Options
    // and the S2Builder.GraphOptions for each of the three output layers.
    private class NormalizeTest
    {
        internal bool SuppressLowerDimensions_ { get; set; }
        internal List<GraphOptions> GraphOptionsOut_ { get; set; } = new();
        internal List<GraphClone> GraphClones { get; set; } = new();

        internal NormalizeTest()
        {
            SuppressLowerDimensions_ = true;
            // Set the default GraphOptions for building S2Points, S2Polylines, and
            // S2Polygons.  Tests can modify these options as necessary.  Most of the
            // defaults are KEEP so that we can verify edge counts in some cases.
            //
            // Polyline edges are undirected by default because (1) this case is
            // slightly more challenging and (2) it is expected to be common.
            GraphOptionsOut_.Add(  // Points
                new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP));
            GraphOptionsOut_.Add(  // Polylines
                new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP));
            GraphOptionsOut_.Add(  // Polygons
                new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP));
        }

        internal void Run(string input_str, string expected_str)
        {
            ClosedSetNormalizer.Options options = new()
            {
                SuppressLowerDimensions = (SuppressLowerDimensions_)
            };
            ClosedSetNormalizer normalizer = new(options, GraphOptionsOut_);

            S2Builder builder = new(new Options());
            List<Graph> input = new(), expected = new();
            AddLayers(input_str, normalizer.GraphOptions_, input, builder);
            AddLayers(expected_str, GraphOptionsOut_, expected, builder);
            // Populate the "input" and "expected" vectors.
            Assert.True(builder.Build(out _));

            var actual = normalizer.Run(input, out _);
            for (int dim = 0; dim < 3; ++dim)
            {
                Assert.True(expected[dim].Options == actual[dim].Options);
                Assert.Equal(S2TextFormat.ToDebugString(expected[dim]), S2TextFormat.ToDebugString(actual[dim]));
            }
        }

        private void AddLayers(string str, List<GraphOptions> graph_options, List<Graph> graphs_out, S2Builder builder)
        {
            var index = MakeIndexOrDie(str);
            for (int dim = 0; dim < 3; ++dim)
            {
                builder.StartLayer(new GraphAppendingLayer(graph_options[dim], graphs_out, GraphClones));
                foreach (S2Shape shape in index)
                {
                    if (shape.Dimension() != dim) continue;
                    int n = shape.NumEdges();
                    for (int e = 0; e < n; ++e)
                    {
                        S2Shape.Edge edge = shape.GetEdge(e);
                        builder.AddEdge(edge.V0, edge.V1);
                    }
                }
            }
        }

    }

    // If this code changes, please update the header file comments to match.
    private static bool ComputeUnion(S2ShapeIndex a, S2ShapeIndex b, MutableS2ShapeIndex index, out S2Error error)
    {
        S2PolylineVectorLayer.Options polyline_options = new();
        polyline_options.EdgeType_ = (EdgeType.UNDIRECTED);
        polyline_options.PolylineType_ = (Graph.PolylineType.WALK);
        polyline_options.DuplicateEdges_ = (DuplicateEdges.MERGE);
        var layers = new List<Layer>
            {
                new IndexedS2PointVectorLayer(index),
                new IndexedS2PolylineVectorLayer(index, polyline_options),
                new IndexedS2PolygonLayer(index),
            };
        S2BooleanOperation op = new(S2BooleanOperation.OpType.UNION, layers.NormalizeClosedSet());
        return op.Build(a, b, out error);
    }
}
