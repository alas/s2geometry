// Most of S2Builder.Graph is tested by the S2Builder.Layer implementations
// rather than here.

using System;
using System.Collections.Generic;
using Xunit;

namespace S2Geometry
{
    using static S2Builder;
    using static S2Builder.GraphOptions;
    using static S2Builder.Graph;

    using Edge = KeyKey<int, int>;
    using DirectedComponent = List<List<int>>;
    using UndirectedComponent = Array2<List<List<int>>>;

    public class S2BuilderGraphTests
    {
        [Fact]
        public void Test_GetDirectedLoops_DegenerateEdges()
        {
            GraphClone gc = new();
            S2Builder builder = new(new S2Builder.Options());
            GraphOptions graph_options = new(
                EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                DuplicateEdges.KEEP, SiblingPairs.KEEP);
            builder.StartLayer(new GraphCloningLayer(graph_options, gc));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("0:3, 3:3, 0:3"));
            Assert.True(builder.Build(out _));
            Graph g = gc.Graph();
            DirectedComponent loops = new();
            Assert.True(g.GetDirectedLoops(LoopType.SIMPLE, loops, out _));
            Assert.Equal(3, loops.Count);
            Assert.Single(loops[0]);
            Assert.Equal(4, loops[1].Count);
            Assert.Equal(2, loops[2].Count);
        }

        [Fact]
        public void Test_GetDirectedComponents_DegenerateEdges()
        {
            GraphClone gc = new();
            S2Builder builder = new(new S2Builder.Options());
            GraphOptions graph_options = new(
                EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                DuplicateEdges.MERGE, SiblingPairs.CREATE);
            builder.StartLayer(new GraphCloningLayer(graph_options, gc));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
            Assert.True(builder.Build(out _));
            Graph g = gc.Graph();
            List<DirectedComponent> components = new();
            Assert.True(g.GetDirectedComponents(DegenerateBoundaries.KEEP, components, out _));
            Assert.Equal(2, components.Count);
            Assert.Single(components[0]);
            Assert.Single(components[0][0]);
            Assert.Equal(2, components[1].Count);
            Assert.Equal(4, components[1][0].Count);
            Assert.Equal(4, components[1][1].Count);
        }

        [Fact]
        public void Test_GetUndirectedComponents_DegenerateEdges()
        {
            GraphClone gc = new();
            S2Builder builder = new(new S2Builder.Options());
            GraphOptions graph_options = new(
                EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
                DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
            builder.StartLayer(new GraphCloningLayer(graph_options, gc));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
            Assert.True(builder.Build(out _));
            Graph g = gc.Graph();
            List<UndirectedComponent> components = new();
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
        public void Test_GetPolylines_UndirectedDegeneratePaths()
        {
            GraphClone gc = new();
            S2Builder builder = new(new S2Builder.Options());
            GraphOptions graph_options = new(
                EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                DuplicateEdges.KEEP, SiblingPairs.KEEP);
            builder.StartLayer(new GraphCloningLayer(graph_options, gc));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            Assert.True(builder.Build(out _));
            Graph g = gc.Graph();
            var polylines = g.GetPolylines(PolylineType.PATH);
            Assert.Equal(7, polylines.Count);
        }

        [Fact]
        public void Test_GetPolylines_UndirectedDegenerateWalks()
        {
            GraphClone gc = new();
            S2Builder builder = new(new S2Builder.Options());
            GraphOptions graph_options = new(
                EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                DuplicateEdges.KEEP, SiblingPairs.KEEP);
            builder.StartLayer(new GraphCloningLayer(graph_options, gc));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
            builder.AddShape(S2TextFormat.MakeLaxPolylineOrDie("1:1, 1:1"));
            Assert.True(builder.Build(out _));
            Graph g = gc.Graph();
            var polylines = g.GetPolylines(PolylineType.WALK);
            Assert.Equal(2, polylines.Count);
            Assert.Equal(2, polylines[0].Count);
            Assert.Equal(5, polylines[1].Count);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.KEEP);
            TestProcessEdges(new TestEdge[] { new(0, 0), new(0, 0) }, expected: Array.Empty<TestEdge>(), options);
        }

        [Fact]
        public void Test_ProcessEdges_KeepDuplicateDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                                 DuplicateEdges.KEEP, SiblingPairs.KEEP);
            TestProcessEdges(new TestEdge[] { new(0, 0), new(0, 0) }, new TestEdge[] { new(0, 0), new(0, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_MergeDuplicateDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                                 DuplicateEdges.MERGE, SiblingPairs.KEEP);
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int> { 1 }), new(0, 0, new List<int> { 2 }) },
                new TestEdge[] { new(0, 0, new List<int>{ 1, 2 }) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_MergeUndirectedDuplicateDegenerateEdges()
        {
            // Edge count should be reduced to 2 (i.e., one undirected edge), and all
            // labels should be merged.
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                                 DuplicateEdges.MERGE, SiblingPairs.KEEP);
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0), new(0, 0), new(0, 0, new List<int>{ 2 }) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1, 2 }), new(0, 0, new List<int>{ 1, 2 }) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_ConvertedUndirectedDegenerateEdges()
        {
            // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
            // merges any edge labels.
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                                 DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int> { 1 }), new(0, 0), new(0, 0), new(0, 0, new List<int> { 2 }) },
                   new TestEdge[] { new(0, 0, new List<int> { 1, 2 }), new(0, 0, new List<int> { 1, 2 }) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
        }

        [Fact]
        public void Test_ProcessEdges_MergeConvertedUndirectedDuplicateDegenerateEdges()
        {
            // Like the test above, except that we also merge duplicates.
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
                                 DuplicateEdges.MERGE, SiblingPairs.REQUIRE);
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0), new(0, 0), new(0, 0, new List<int>{ 2 }) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1, 2 }) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessConnectedDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                                 DuplicateEdges.KEEP, SiblingPairs.KEEP);
            // Test that degenerate edges are discarded if they are connnected to any
            // non-degenerate edges (whether they are incoming or outgoing, and whether
            // they are lexicographically before or after the degenerate edge).
            TestProcessEdges(new TestEdge[] { new(0, 0), new(0, 1) }, new TestEdge[] { new(0, 1) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 0), new(1, 0) }, new TestEdge[] { new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 1) }, new TestEdge[] { new(0, 1) }, options);
            TestProcessEdges(new TestEdge[] { new(1, 0), new(1, 1) }, new TestEdge[] { new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessIsolatedDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                                 DuplicateEdges.KEEP, SiblingPairs.KEEP);
            // Test that DISCARD_EXCESS does not merge any remaining duplicate
            // degenerate edges together.
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0, new List<int>{ 2 }) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0, new List<int>{ 2 }) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessUndirectedIsolatedDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
                                 DuplicateEdges.KEEP, SiblingPairs.KEEP);
            // Test that DISCARD_EXCESS does not merge any remaining duplicate
            // undirected degenerate edges together.
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0), new(0, 0, new List<int>{ 2 }), new(0, 0) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0), new(0, 0, new List<int>{ 2 }), new(0, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessConvertedUndirectedIsolatedDegenerateEdges()
        {
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
                                 DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
            // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
            // merges edge labels.
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0, new List<int>{ 2 }), new(0, 0, new List<int>{ 3 }), new(0, 0) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1, 2, 3 }), new(0, 0, new List<int>{ 1, 2, 3 }) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
        }

        [Fact]
        public void Test_ProcessEdges_SiblingPairsDiscardMergesDegenerateEdgeLabels()
        {
            // Test that when SiblingPairs.DISCARD or SiblingPairs.DISCARD_EXCESS
            // is specified, the edge labels of degenerate edges are merged together
            // (for consistency, since these options merge the labels of all
            // non-degenerate edges as well).
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                                 DuplicateEdges.KEEP, SiblingPairs.DISCARD);
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0, new List<int>{ 2 }), new(0, 0, new List<int>{ 3 }) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1, 2, 3 }), new(0, 0, new List<int>{ 1, 2, 3 }), new(0, 0, new List<int>{ 1, 2, 3 }) },
                   options);
            options.SiblingPairs_ = (SiblingPairs.DISCARD_EXCESS);
            TestProcessEdges(new TestEdge[] { new(0, 0, new List<int>{ 1 }), new(0, 0, new List<int>{ 2 }), new(0, 0, new List<int>{ 3 }) },
                   new TestEdge[] { new(0, 0, new List<int>{ 1, 2, 3 }), new(0, 0, new List<int>{ 1, 2, 3 }), new(0, 0, new List<int>{ 1, 2, 3 }) },
                   options);
        }

        [Fact]
        public void Test_ProcessEdges_KeepSiblingPairs()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.KEEP);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_MergeDuplicateSiblingPairs()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.MERGE, SiblingPairs.KEEP);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardSiblingPairs()
        {
            // Check that matched pairs are discarded, leaving behind any excess edges.
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.DISCARD);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) }, Array.Empty<TestEdge>(), options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) }, Array.Empty<TestEdge>(), options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0) },
                   new TestEdge[] { new(0, 1), new(0, 1) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0), new(1, 0), new(1, 0) },
                   new TestEdge[] { new(1, 0), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardSiblingPairsMergeDuplicates()
        {
            // Check that matched pairs are discarded, and then any remaining edges
            // are merged.
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.MERGE, SiblingPairs.DISCARD);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) }, Array.Empty<TestEdge>(), options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0), new(1, 0), new(1, 0) }, new TestEdge[] { new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardUndirectedSiblingPairs()
        {
            // An undirected sibling pair consists of four edges, two in each direction
            // (see s2builder.h).  Since undirected edges always come in pairs, this
            // means that the result always consists of either 0 or 2 edges.
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.DISCARD);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) }, Array.Empty<TestEdge>(), options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessSiblingPairs()
        {
            // Like SiblingPairs.DISCARD, except that one sibling pair is kept if the
            // result would otherwise be empty.
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(0, 1) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0), new(1, 0), new(1, 0) },
                             new TestEdge[] { new(1, 0), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessSiblingPairsMergeDuplicates()
        {
            // Like SiblingPairs.DISCARD, except that one sibling pair is kept if the
            // result would otherwise be empty.
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.MERGE, SiblingPairs.DISCARD_EXCESS);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0), new(1, 0), new(1, 0) }, new TestEdge[] { new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_DiscardExcessUndirectedSiblingPairs()
        {
            // Like SiblingPairs.DISCARD, except that one undirected sibling pair
            // (4 edges) is kept if the result would otherwise be empty.
            GraphOptions options = new(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_CreateSiblingPairs()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.CREATE);
            TestProcessEdges(new TestEdge[] { new(0, 1) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1) },
                   new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_RequireSiblingPairs()
        {
            // Like SiblingPairs.CREATE, but generates an error.
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1) }, new TestEdge[] { new(0, 1), new(1, 0) }, options,
                   S2ErrorCode.BUILDER_MISSING_EXPECTED_SIBLING_EDGES);
        }

        [Fact]
        public void Test_ProcessEdges_CreateUndirectedSiblingPairs()
        {
            // An undirected sibling pair consists of 4 edges, but SiblingPairs.CREATE
            // also converts the graph to EdgeType.DIRECTED and cuts the number of
            // edges in half.
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.KEEP, SiblingPairs.CREATE);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) },
                   new TestEdge[] { new(0, 1), new(1, 0) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);

            options.EdgeType_ = (EdgeType.UNDIRECTED);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) },
                   new TestEdge[] { new(0, 1), new(1, 0) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);

            options.EdgeType_ = (EdgeType.UNDIRECTED);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0) },
                   new TestEdge[] { new(0, 1), new(0, 1), new(1, 0), new(1, 0) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
        }

        [Fact]
        public void Test_ProcessEdges_CreateSiblingPairsMergeDuplicates()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.MERGE, SiblingPairs.CREATE);
            TestProcessEdges(new TestEdge[] { new(0, 1) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(0, 1) }, new TestEdge[] { new(0, 1), new(1, 0) }, options);
        }

        [Fact]
        public void Test_ProcessEdges_CreateUndirectedSiblingPairsMergeDuplicates()
        {
            GraphOptions options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                 DuplicateEdges.MERGE, SiblingPairs.CREATE);
            TestProcessEdges(new TestEdge[] { new(0, 1), new(1, 0) },
                             new TestEdge[] { new(0, 1), new(1, 0) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);

            options.EdgeType_ = (EdgeType.UNDIRECTED);
            TestProcessEdges(
                new TestEdge[] { new(0, 1), new(0, 1), new(0, 1), new(1, 0), new(1, 0), new(1, 0) },
                new TestEdge[] { new(0, 1), new(1, 0) }, options);
            Assert.Equal(EdgeType.DIRECTED, options.EdgeType_);
        }

        private static void TestProcessEdges(TestEdge[] input, TestEdge[] expected,
        GraphOptions options, S2ErrorCode expected_code = S2ErrorCode.OK)
        {
            List<Edge> edges = new();
            List<Int32> input_id_set_ids = new();
            IdSetLexicon id_set_lexicon = new();
            foreach (var e in input)
            {
                edges.Add(new Edge(e.Item1, e.Item2));
                if (e.InputIds != null)
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
                Assert.Equal(e.InputIds ?? new List<int>(), actual_ids); // $"(edge {i})";
            }
            Assert.Equal(expected.Length, edges.Count); // "Too many output edges";
        }

        private record TestEdge(int Item1, int Item2, List<int> InputIds = null);
    }
}
