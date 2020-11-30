using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace S2Geometry
{
    using static S2Builder;
    using static S2Builder.GraphOptions;
    
    using Edge = KeyKey<int, int>;

    public class S2BuilderUtil_TestingTests
    {
        [Fact]
        public void Test_GraphCloningLayer_MakeIndependentCopy()
        {
            // Also tests GraphClone.
            GraphClone gc = new();
            S2Builder builder = new(new S2Builder.Options());
            GraphOptions graph_options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
                                       DuplicateEdges.MERGE, SiblingPairs.KEEP);
            builder.StartLayer(new GraphCloningLayer(graph_options, gc));
            S2Point v0 = new(1, 0, 0), v1 = new(0, 1, 0);
            builder.SetLabel(14);
            builder.AddEdge(v0, v1);
            Assert.True(builder.Build(out _));
            Graph g = gc.Graph();
            Assert.True(graph_options == g.Options);
            Assert.Equal(new List<S2Point> { v0, v1 }, g.Vertices);
            Assert.Equal(new List<Edge> { new Edge(0, 1) }, g.Edges);
            Assert.Single(g.InputEdgeIds(0));
            Assert.Equal(0, g.InputEdgeIds(0).First());
            Graph.LabelFetcher fetcher = new(g, g.Options.EdgeType_);
            List<int> labels = new();
            fetcher.Fetch(0, labels);
            Assert.Equal(new List<int> { 14 }, labels);
            // S2Builder sets a default IsFullPolygonPredicate that returns an error.
            Assert.True(g.IsFullPolygonPredicate() != null);
        }

        [Fact]
        public void Test_GraphAppendingLayer_AppendsTwoGraphs()
        {
            List<Graph> graphs = new();
            List<GraphClone> clones = new();
            S2Builder builder = new(new S2Builder.Options());
            builder.StartLayer(new GraphAppendingLayer(
                new GraphOptions(), graphs, clones));
            S2Point v0 = new(1, 0, 0), v1 = new(0, 1, 0);
            builder.AddEdge(v0, v1);
            builder.StartLayer(new GraphAppendingLayer(
                new GraphOptions(), graphs, clones));
            builder.AddEdge(v1, v0);
            Assert.True(builder.Build(out _));
            Assert.Equal(2, graphs.Count);
            Assert.Equal(2, clones.Count);
            Assert.Equal(v0, graphs[0].Vertex(graphs[0].GetEdge(0).Item1));
            Assert.Equal(v1, graphs[1].Vertex(graphs[1].GetEdge(0).Item1));
        }
    }
}
