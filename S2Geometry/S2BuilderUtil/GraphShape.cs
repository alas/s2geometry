namespace S2Geometry.S2BuilderUtil;

using Graph = S2Builder.Graph;

// An S2Shape representing the edges in an S2Builder.Graph.
public sealed class GraphShape : S2Shape
{
    public GraphShape(Graph g) => g_ = g;

    public override int NumEdges() => g_.NumEdges;

    public override Edge GetEdge(int e)
    {
        var g_edge = g_.GetEdge(e);
        return new Edge(g_.Vertex(g_edge.ShapeId), g_.Vertex(g_edge.EdgeId));
    }
    public override int Dimension() => 1;
    public override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public override int NumChains() => g_.NumEdges;
    public override Chain GetChain(int i) => new(i, 1);
    public override Edge ChainEdge(int i, int j)
    {
        MyDebug.Assert(j == 0);
        return GetEdge(i);
    }
    public override ChainPosition GetChainPosition(int e) => new(e, 0);

    private readonly Graph g_;
}
