// Most of S2Builder.Graph is tested by the S2Builder.Layer implementations
// rather than here.

namespace S2Geometry;

using static S2Builder;
using static S2Builder.GraphOptions;
using static S2Builder.Graph;

public class S2BuilderGraphTests
{
    [Fact]
    internal void Test_Graph_LabelsRequestedButNotProvided()
    {
        // Tests the situation where labels are requested but none were provided.
        GraphOptions options=new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        List<S2Point> vertices=[new S2Point(1, 0, 0)];
        List<Edge> edges=[new(0, 0)];
        List<InputEdgeIdSetId> input_edge_id_set_ids=[0];
        IdSetLexicon input_edge_id_set_lexicon=new(), label_set_lexicon=new();
        List<LabelSetId> label_set_ids=[];  // Empty means no labels are present.
        Graph g = new(
            options, vertices, edges, input_edge_id_set_ids,
            input_edge_id_set_lexicon, label_set_ids, label_set_lexicon, null);
        Assert.True(g.LabelSetIds!.Count==0);
        Assert.Equal(g.LabelSetId(0), IdSetLexicon.kEmptySetId);
        Assert.Equal(g.Labels(0).Count, 0);  // Labels for input edge 0.
        LabelFetcher fetcher = new(g, EdgeType.DIRECTED);
        List<Label> labels=[];
        fetcher.Fetch(0, labels);         // Labels for graph edge 0.
        Assert.True(labels.Count==0);
    }

    [Fact]
    internal void Test_GetDirectedLoops_DegenerateEdges()
    {
        GraphClone gc = new();
        S2Builder builder = new(new Options());
        GraphOptions graph_options = new(
            EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP, SiblingPairs.KEEP);
        builder.StartLayer(new GraphCloningLayer(graph_options, gc));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        builder.AddShape(MakeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
        builder.AddShape(MakeLaxPolylineOrDie("0:3, 3:3, 0:3"));
        Assert.True(builder.Build(out _));
        Graph g = gc.Graph();
        DirectedComponent loops = [];
        Assert.True(g.GetDirectedLoops(LoopType.SIMPLE, loops, out _));
        Assert.Equal(3, loops.Count);
        Assert.Single(loops[0]);
        Assert.Equal(4, loops[1].Count);
        Assert.Equal(2, loops[2].Count);
    }

    [Fact]
    internal void Test_GetDirectedComponents_DegenerateEdges()
    {
        GraphClone gc = new();
        S2Builder builder = new(new Options());
        GraphOptions graph_options = new(
            EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP, SiblingPairs.CREATE);
        builder.StartLayer(new GraphCloningLayer(graph_options, gc));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        builder.AddShape(MakeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
        Assert.True(builder.Build(out _));
        Graph g = gc.Graph();
        List<DirectedComponent> components = [];
        Assert.True(g.GetDirectedComponents(DegenerateBoundaries.KEEP, components, out _));
        Assert.Equal(2, components.Count);
        Assert.Single(components[0]);
        Assert.Single(components[0][0]);
        Assert.Equal(2, components[1].Count);
        Assert.Equal(4, components[1][0].Count);
        Assert.Equal(4, components[1][1].Count);
    }

    [Fact]
    internal void Test_GetUndirectedComponents_DegenerateEdges()
    {
        GraphClone gc = new();
        S2Builder builder = new(new Options());
        GraphOptions graph_options = new(
            EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
        builder.StartLayer(new GraphCloningLayer(graph_options, gc));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        builder.AddShape(MakeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
        Assert.True(builder.Build(out _));
        Graph g = gc.Graph();
        List<UndirectedComponent> components = [];
        Assert.True(g.GetUndirectedComponents(LoopType.CIRCUIT, components, out _));
        // The result consists of two components, each with two complements.  Each
        // complement in this example has exactly one loop.  The loops in both
        // complements of the first component have 1 vertex, while the loops in both
        // complements of the second component have 4 vertices.
        Assert.Equal(2, components.Count);
        Assert.Single(components[0][0]);
        Assert.Single(components[0][0][0]);
        Assert.Single(components[0][1]);
        Assert.Single(components[0][1][0]);
        Assert.Single(components[1][0]);
        Assert.Equal(4, components[1][0][0].Count);
        Assert.Single(components[1][1]);
        Assert.Equal(4, components[1][1][0].Count);
    }

    [Fact]
    internal void Test_GetPolylines_UndirectedDegeneratePaths()
    {
        GraphClone gc = new();
        S2Builder builder = new(new Options());
        GraphOptions graph_options = new(
            EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
            DuplicateEdges.KEEP, SiblingPairs.KEEP);
        builder.StartLayer(new GraphCloningLayer(graph_options, gc));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        builder.AddShape(MakeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        Assert.True(builder.Build(out _));
        Graph g = gc.Graph();
        var polylines = g.GetPolylines(PolylineType.PATH);
        Assert.Equal(7, polylines.Count);
    }

    [Fact]
    internal void Test_GetPolylines_UndirectedDegenerateWalks()
    {
        GraphClone gc = new();
        S2Builder builder = new(new Options());
        GraphOptions graph_options = new(
            EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
            DuplicateEdges.KEEP, SiblingPairs.KEEP);
        builder.StartLayer(new GraphCloningLayer(graph_options, gc));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        builder.AddShape(MakeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
        builder.AddShape(MakeLaxPolylineOrDie("1:1, 1:1"));
        Assert.True(builder.Build(out _));
        Graph g = gc.Graph();
        var polylines = g.GetPolylines(PolylineType.WALK);
        Assert.Equal(2, polylines.Count);
        Assert.Equal(2, polylines[0].Count);
        Assert.Equal(5, polylines[1].Count);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        TestProcessEdges([new(0, 0), new(0, 0)], expected: [], options);
    }

    [Fact]
    internal void Test_ProcessEdges_KeepDuplicateDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        TestProcessEdges([new(0, 0), new(0, 0)], [new(0, 0), new(0, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_MergeDuplicateDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.MERGE, SiblingPairs.KEEP);
        TestProcessEdges([new(0, 0, [1]), new(0, 0, [2])],
            [new(0, 0, [1, 2])], options);
    }

    [Fact]
    internal void Test_ProcessEdges_MergeUndirectedDuplicateDegenerateEdges()
    {
        // Edge count should be reduced to 2 (i.e., one undirected edge), and all
        // labels should be merged.
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.MERGE, SiblingPairs.KEEP);
        TestProcessEdges([new(0, 0, [1]), new(0, 0), new(0, 0), new(0, 0, [2])],
               [new(0, 0, [1, 2]), new(0, 0, [1, 2])], options);
    }

    [Fact]
    internal void Test_ProcessEdges_ConvertedUndirectedDegenerateEdges()
    {
        // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
        // merges any edge labels.
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
        TestProcessEdges(
            [new(0, 0, [1]), new(0, 0), new(0, 0), new(0, 0, [2])],
            [new(0, 0, [1, 2]), new(0, 0, [1, 2])], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
    }

    [Fact]
    internal void Test_ProcessEdges_MergeConvertedUndirectedDuplicateDegenerateEdges()
    {
        // Like the test above, except that we also merge duplicates.
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.MERGE, SiblingPairs.REQUIRE);
        TestProcessEdges([new(0, 0, [1]), new(0, 0), new(0, 0), new(0, 0, [2])],
               [new(0, 0, [1, 2])], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessConnectedDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        // Test that degenerate edges are discarded if they are connnected to any
        // non-degenerate edges (whether they are incoming or outgoing, and whether
        // they are lexicographically before or after the degenerate edge).
        TestProcessEdges([new(0, 0), new(0, 1)], [new(0, 1)], options);
        TestProcessEdges([new(0, 0), new(1, 0)], [new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(1, 1)], [new(0, 1)], options);
        TestProcessEdges([new(1, 0), new(1, 1)], [new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessIsolatedDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        // Test that DISCARD_EXCESS merges any duplicate degenerate edges together.
        TestProcessEdges(
            [new(0, 0, [1]), new(0, 0, [2])],
            [new(0, 0, [1, 2])], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessUndirectedIsolatedDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        // Test that DISCARD_EXCESS merges any duplicate undirected degenerate edges
        // together.
        TestProcessEdges(
            [new(0, 0, [1]), new(0, 0), new(0, 0, [2]), new(0, 0)],
            [new(0, 0, [1, 2]), new(0, 0, [1, 2])], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessConvertedUndirectedIsolatedDegenerateEdges()
    {
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
                             DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
        // Test that DISCARD_EXCESS with SiblingPairs::REQUIRE merges any duplicate
        // edges together and converts the edges from UNDIRECTED to DIRECTED.
        TestProcessEdges(
            [new(0, 0, [1]), new(0, 0, [2]), new(0, 0, [3]), new(0, 0)],
            [new(0, 0, [1, 2, 3])], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
    }

    [Fact]
    internal void Test_ProcessEdges_SiblingPairsDiscardMergesDegenerateEdgeLabels()
    {
        // Test that when SiblingPairs.DISCARD or SiblingPairs.DISCARD_EXCESS
        // is specified, the edge labels of degenerate edges are merged together
        // (for consistency, since these options merge the labels of all
        // non-degenerate edges as well).
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.DISCARD);
        TestProcessEdges([new(0, 0, [1]), new(0, 0, [2]), new(0, 0, [3])],
               [new(0, 0, [1, 2, 3]), new(0, 0, [1, 2, 3]), new(0, 0, [1, 2, 3])],
               options);
        options.SiblingPairs_ = SiblingPairs.DISCARD_EXCESS;
        TestProcessEdges([new(0, 0, [1]), new(0, 0, [2]), new(0, 0, [3])],
               [new(0, 0, [1, 2, 3]), new(0, 0, [1, 2, 3]), new(0, 0, [1, 2, 3])],
               options);
    }

    [Fact]
    internal void Test_ProcessEdges_KeepSiblingPairs()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        TestProcessEdges([new(0, 1), new(1, 0)], [new(0, 1), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_MergeDuplicateSiblingPairs()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.MERGE, SiblingPairs.KEEP);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0)], [new(0, 1), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardSiblingPairs()
    {
        // Check that matched pairs are discarded, leaving behind any excess edges.
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.DISCARD);
        TestProcessEdges([new(0, 1), new(1, 0)], [], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)], [], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0)],
               [new(0, 1), new(0, 1)], options);
        TestProcessEdges([new(0, 1), new(1, 0), new(1, 0), new(1, 0)],
               [new(1, 0), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardSiblingPairsMergeDuplicates()
    {
        // Check that matched pairs are discarded, and then any remaining edges
        // are merged.
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.MERGE, SiblingPairs.DISCARD);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)], [], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0)], [new(0, 1)], options);
        TestProcessEdges([new(0, 1), new(1, 0), new(1, 0), new(1, 0)], [new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardUndirectedSiblingPairs()
    {
        // An undirected sibling pair consists of four edges, two in each direction
        // (see s2builder.h).  Since undirected edges always come in pairs, this
        // means that the result always consists of either 0 or 2 edges.
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.DISCARD);
        TestProcessEdges([new(0, 1), new(1, 0)], [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)], [], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0)],
                         [new(0, 1), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessSiblingPairs()
    {
        // Like SiblingPairs.DISCARD, except that one sibling pair is kept if the
        // result would otherwise be empty.
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
        TestProcessEdges([new(0, 1), new(1, 0)], [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)],
                         [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0)],
                         [new(0, 1), new(0, 1)], options);
        TestProcessEdges([new(0, 1), new(1, 0), new(1, 0), new(1, 0)],
                         [new(1, 0), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessSiblingPairsMergeDuplicates()
    {
        // Like SiblingPairs.DISCARD, except that one sibling pair is kept if the
        // result would otherwise be empty.
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.MERGE, SiblingPairs.DISCARD_EXCESS);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)],
                         [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0)], [new(0, 1)], options);
        TestProcessEdges([new(0, 1), new(1, 0), new(1, 0), new(1, 0)], [new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_DiscardExcessUndirectedSiblingPairs()
    {
        // Like SiblingPairs.DISCARD, except that one undirected sibling pair
        // (4 edges) is kept if the result would otherwise be empty.
        GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
        TestProcessEdges([new(0, 1), new(1, 0)], [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)],
                         [new(0, 1), new(0, 1), new(1, 0), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0)],
                         [new(0, 1), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_CreateSiblingPairs()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.CREATE);
        TestProcessEdges([new(0, 1)], [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1)],
               [new(0, 1), new(0, 1), new(1, 0), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_RequireSiblingPairs()
    {
        // Like SiblingPairs.CREATE, but generates an error.
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
        TestProcessEdges([new(0, 1), new(1, 0)], [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1)], [new(0, 1), new(1, 0)], options,
               S2ErrorCode.BUILDER_MISSING_EXPECTED_SIBLING_EDGES);
    }

    [Fact]
    internal void Test_ProcessEdges_CreateUndirectedSiblingPairs()
    {
        // An undirected sibling pair consists of 4 edges, but SiblingPairs.CREATE
        // also converts the graph to EdgeType.DIRECTED and cuts the number of
        // edges in half.
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.CREATE);
        TestProcessEdges([new(0, 1), new(1, 0)],
               [new(0, 1), new(1, 0)], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);

        options.EdgeType_ = EdgeType.UNDIRECTED;
        TestProcessEdges([new(0, 1), new(0, 1), new(1, 0), new(1, 0)],
               [new(0, 1), new(1, 0)], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);

        options.EdgeType_ = EdgeType.UNDIRECTED;
        TestProcessEdges([new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0)],
               [new(0, 1), new(0, 1), new(1, 0), new(1, 0)], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
    }

    [Fact]
    internal void Test_ProcessEdges_CreateSiblingPairsMergeDuplicates()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.MERGE, SiblingPairs.CREATE);
        TestProcessEdges([new(0, 1)], [new(0, 1), new(1, 0)], options);
        TestProcessEdges([new(0, 1), new(0, 1)], [new(0, 1), new(1, 0)], options);
    }

    [Fact]
    internal void Test_ProcessEdges_CreateUndirectedSiblingPairsMergeDuplicates()
    {
        GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.MERGE, SiblingPairs.CREATE);
        TestProcessEdges([new(0, 1), new(1, 0)],
                         [new(0, 1), new(1, 0)], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);

        options.EdgeType_ = EdgeType.UNDIRECTED;
        TestProcessEdges(
            [new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0)],
            [new(0, 1), new(1, 0)], options);
        Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
    }

    private static void TestProcessEdges(TestEdge[] input, TestEdge[] expected,
    GraphOptions options, S2ErrorCode expected_code = S2ErrorCode.OK)
    {
        List<Edge> edges = [];
        List<Int32> input_id_set_ids = [];
        IdSetLexicon id_set_lexicon = new();
        foreach (var e in input)
        {
            edges.Add(new Edge(e.Item1, e.Item2));
            if (e.InputIds is not null)
                input_id_set_ids.Add(id_set_lexicon.Add(e.InputIds));
        }
        Graph.ProcessEdges(options, edges, input_id_set_ids, id_set_lexicon, out var error);
        Assert.Equal(expected_code, error.Code);
        Assert.Equal(edges.Count, input_id_set_ids.Count);
        for (int i = 0; i < expected.Length; ++i)
        {
            var e = expected[i];
            Assert.True(i < edges.Count); // Not enough output edges
            Assert.Equal(new Edge(e.Item1, e.Item2), edges[i]); // $"(edge {i})";
            var id_set = id_set_lexicon.IdSet_(input_id_set_ids[i]);
            List<Int32> actual_ids = new(id_set);
            Assert.Equal(e.InputIds ?? [], actual_ids); // $"(edge {i})";
        }
        Assert.Equal(expected.Length, edges.Count); // "Too many output edges";
    }

    private void TestMakeSubgraph(
        Graph g,
        IdSetLexicon new_input_edge_id_set_lexicon,
        GraphOptions new_options,
        List<Edge> new_edges,
        List<InputEdgeIdSetId> new_input_edge_id_set_ids,
        GraphOptions expected_options,
        List<Edge> expected_edges,
        List<InputEdgeIdSetId> expected_input_edge_id_set_ids)
    {
        Graph new_g = g.MakeSubgraph(
            new_options, new_edges, new_input_edge_id_set_ids,
            new_input_edge_id_set_lexicon, null, out _)!;

        // Some parts of the graph should be the same.
        Assert.True(new_g.Vertices == g.Vertices);
        Assert.Equal(new_g.LabelSetIds, g.LabelSetIds);
        Assert.Equal(new_g.LabelSetLexicon, g.LabelSetLexicon);

        // Some parts of the graph should use the provided underlying storage.
        Assert.True(new_g.Edges == new_edges);
        Assert.True(new_g.InputEdgeIdSetIds == new_input_edge_id_set_ids);
        Assert.True(new_g.InputEdgeIdSetLexicon ==
          new_input_edge_id_set_lexicon);

        // The new graph should have the expected options.
        Assert.Equal(new_g.Options.EdgeType_, expected_options.EdgeType_);
        Assert.Equal(new_g.Options.DegenerateEdges_,
              expected_options.DegenerateEdges_);
        Assert.Equal(new_g.Options.DuplicateEdges_,
              expected_options.DuplicateEdges_);
        Assert.Equal(new_g.Options.SiblingPairs_, expected_options.SiblingPairs_);

        // The edges should be updated according to the requested options.
        Assert.Equal(new_g.Edges, expected_edges);
        Assert.Equal(new_g.InputEdgeIdSetIds, expected_input_edge_id_set_ids);
        Assert.Equal(new_g.InputEdgeIdSetLexicon,
              new_input_edge_id_set_lexicon);
    }

    [Fact]
    internal void Test_MakeSubgraph_UndirectedToUndirected()
{
    // Test that MakeSubgraph() doesn't transform edges into edge pairs when
    // creating an undirected subgraph of an undirected graph.
    GraphOptions options=new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                         DuplicateEdges.KEEP, SiblingPairs.KEEP);
    var vertices = ParsePointsOrDie("0:0, 0:1, 1:1");
    List<Edge> edges=[new(0, 0), new(0, 0), new(1, 2), new(2, 1)];
        List<InputEdgeIdSetId> input_edge_id_set_ids=[0, 0, 1, 1];
        List<LabelSetId> label_set_ids=[];
    IdSetLexicon input_edge_id_set_lexicon=new(), label_set_lexicon=new();
    Graph graph=new(
        options, vertices, edges, input_edge_id_set_ids,
        input_edge_id_set_lexicon, label_set_ids, label_set_lexicon, null);

    // Now create a subgraph with undirected edges but different options.
    GraphOptions new_options=new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
                             DuplicateEdges.KEEP, SiblingPairs.KEEP);
        List<Edge> expected_edges=[new(1, 2), new(2, 1)];
        List<InputEdgeIdSetId> expected_input_edge_id_set_ids=[1, 1];
    TestMakeSubgraph(
        graph, input_edge_id_set_lexicon,
        new_options, edges, input_edge_id_set_ids,
        new_options, expected_edges, expected_input_edge_id_set_ids);
    }

    [Fact]
    internal void Test_MakeSubgraph_DirectedToUndirected()
{
    // Test transforming directed edges into undirected edges (which doubles the
    // number of edges).
    GraphOptions options=new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                         DuplicateEdges.KEEP, SiblingPairs.KEEP);
    var vertices = ParsePointsOrDie("0:0, 0:1, 1:1");
        List<Edge> edges = [new(0, 0), new(0, 1), new(1, 2), new(1, 2), new(2, 1)];
        List<InputEdgeIdSetId> input_edge_id_set_ids = [1, 2, 3, 3, 3];
        List<LabelSetId> label_set_ids = [];
    IdSetLexicon input_edge_id_set_lexicon = new(), label_set_lexicon = new();
    Graph graph = new(
        options, vertices, edges, input_edge_id_set_ids,
        input_edge_id_set_lexicon, label_set_ids, label_set_lexicon, null);

    // Now create a subgraph with undirected edges and different options.
    GraphOptions new_options=new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                             DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
        List<Edge> expected_edges = [
            new(0, 0), new(0, 0),  // Undirected degenerate edge.
new (0, 1), new (1, 0),  // Undirected edge.
new(1, 2), new(2, 1)   // Undirected edge after discarding sibling pair.
    ];
    List<InputEdgeIdSetId> expected_input_edge_id_set_ids =
    [
        1, 1, 2, IdSetLexicon.kEmptySetId, 3, 3
    ];
    TestMakeSubgraph(
        graph, input_edge_id_set_lexicon,
        new_options, edges, input_edge_id_set_ids,
        new_options, expected_edges, expected_input_edge_id_set_ids);
}


private record TestEdge(int Item1, int Item2, List<int>? InputIds = null);
}
