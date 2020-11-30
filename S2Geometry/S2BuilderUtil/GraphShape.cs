namespace S2Geometry.S2BuilderUtil
{
    using Graph = S2Builder.Graph;

    // An S2Shape representing the edges in an S2Builder.Graph.
    public sealed class GraphShape : S2Shape
    {
        public GraphShape(Graph g) { g_ = g; }
        public override int NumEdges => g_.NumEdges; public override Edge GetEdge(int e)
        {
            var g_edge = g_.GetEdge(e);
            return new Edge(g_.Vertex(g_edge.Item1), g_.Vertex(g_edge.Item2));
        }
        public override int Dimension() { return 1; }
        public override ReferencePoint GetReferencePoint()
        {
            return ReferencePoint.FromContained(false);
        }
        public override int NumChains() { return g_.NumEdges; }
        public override Chain GetChain(int i) { return new Chain(i, 1); }
        public override Edge ChainEdge(int i, int j)
        {
            Assert.True(j == 0);
            return GetEdge(i);
        }
        public override ChainPosition GetChainPosition(int e)
        {
            return new ChainPosition(e, 0);
        }

        private readonly Graph g_;
    }
}
