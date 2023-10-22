// A class that copies an Graph and owns the underlying data
// (unlike Graph, which is just a view).

namespace S2Geometry;

using EdgeVector = List<S2Shape.Edge>;

internal class GraphClone
{
    internal GraphClone() { }  // Must call Init().
    internal GraphClone(Graph g) { Init(g); }
    internal void Init(Graph g)
    {
        options_ = g.Options;
        vertices_ = g.Vertices;
        edges_ = g.Edges;
        input_edge_id_set_ids_ = g.InputEdgeIdSetIds;
        input_edge_id_set_lexicon_ = g.InputEdgeIdSetLexicon;
        label_set_ids_ = g.LabelSetIds;
        label_set_lexicon_ = g.LabelSetLexicon;
        is_full_polygon_predicate_ = g.IsFullPolygonPredicate();
        g_ = new Graph(
            options_, vertices_, edges_, input_edge_id_set_ids_,
            input_edge_id_set_lexicon_, label_set_ids_, label_set_lexicon_,
            is_full_polygon_predicate_);
    }
    internal Graph Graph() => g_;

    private GraphOptions options_;
    private List<S2Point> vertices_;
    private List<Edge> edges_;
    private List<Int32> input_edge_id_set_ids_;
    private IdSetLexicon input_edge_id_set_lexicon_;
    private List<Int32> label_set_ids_;
    private IdSetLexicon label_set_lexicon_;
    private IsFullPolygonPredicate is_full_polygon_predicate_;
    private Graph g_;
}

// A layer type that copies an Graph into a GraphClone object
// (which owns the underlying data, unlike Graph itself).
internal class GraphCloningLayer(GraphOptions graph_options, GraphClone gc) : Layer
{
    private readonly GraphClone gc_ = gc;

    public override GraphOptions GraphOptions_() => graph_options_;
    private readonly GraphOptions graph_options_ = graph_options;

    public override void Build(Graph g, out S2Error error)
    { error = S2Error.OK; gc_.Init(g); }
}

// A layer type that copies an Graph and appends it to a vector,
// and appends the corresponding GraphClone object (which owns the Graph data)
// to a separate vector.
internal class GraphAppendingLayer(GraphOptions graph_options, List<Graph> graphs, List<GraphClone> clones) : Layer
{
    private readonly List<Graph> graphs_ = graphs;
    private readonly List<GraphClone> clones_ = clones;

    public override GraphOptions GraphOptions_() => graph_options_;
    private readonly GraphOptions graph_options_ = graph_options;

    public override void Build(Graph g, out S2Error error)
    {
        error = S2Error.OK;
        clones_.Add(new GraphClone(g));
        graphs_.Add(clones_.Last().Graph());
    }
}

// A layer type that expects that the edges in the S2Builder::Graph passed to
// its Build() method should match the edges in the given S2ShapeIndex
// (including multiplicities).  This allows testing whether an algorithm
// produces a given multiset of edges without needing to specify a particular
// ordering of those edges.
internal class IndexMatchingLayer : Layer
{
    // Tests whether the edges passed to its Build() method match the edges in
    // the given S2ShapeIndex (including multiplicities).  If any differences
    // are found, sets "error" to a descriptive error message.
    //
    // If "dimension" is non-negative then only shapes of the given dimension
    // are used.  (This makes allows use with classes such as S2BooleanOperation
    // that output one S2Builder::Graph for each dimension.)
    internal IndexMatchingLayer(GraphOptions graph_options,
                              S2ShapeIndex index, int dimension = -1)
    {
        graph_options_ = graph_options; index_ = index; dimension_ = dimension;
    }

    // S2Builder interface:
    public override GraphOptions GraphOptions_()
    {
        return graph_options_;
    }

    public override void Build(Graph g, out S2Error error)
    {
        error = S2Error.OK;
        List<S2Shape.Edge> actual=[], expected=[];
        for (int e = 0; e < g.NumEdges; ++e)
        {
            var edge = g.GetEdge(e);
            actual.Add(new S2Shape.Edge(g.Vertex(edge.ShapeId),
                g.Vertex(edge.EdgeId)));
        }
        foreach (var shape in index_)
        {
            if (shape is null) continue;
            if (dimension_ >= 0 && shape.Dimension() != dimension_) continue;
            for (int e = shape.NumEdges(); --e >= 0;)
            {
                expected.Add(shape.GetEdge(e));
            }
        }
        actual.Sort();
        expected.Sort();

        // The edges are a multiset, so we can't use std::set_difference.
        List<S2Shape.Edge> missing=[], extra = [];
        var ai = actual.FirstOrDefault();
        var ei = expected.FirstOrDefault();
        var ai_i = 0;
        var ei_i = 0;
        var limit = actual.Count - 1;
        while (ai_i < actual.Count || ei_i < expected.Count)
        {
            if (ei_i == expected.Count || (ai_i != actual.Count && ai < ei))
            {
                extra.Add(ai); ai_i++;
            }
            else if (ai_i == limit || ei < ai)
            {
                missing.Add(ei); ei_i++;
            }
            else
            {
                ++ai_i;
                ++ei_i;
            }
            ai = actual[ai_i];
            ei = expected[ei_i];
        }
        if (missing.Count!=0 || extra.Count!=0)
        {
            // There may be errors in more than one dimension, so we append to the
            // existing error text.
            string label = "";
            if (dimension_ >= 0) label = $"Dimension {dimension_}: ";
            error = new S2Error(S2ErrorCode.FAILED_PRECONDITION,
                $"{label}Missing edges: {ToString(missing)} Extra edges: {ToString(extra)}\n");
        }
    }

    private static string ToString(EdgeVector edges)
    {
        string msg = "";
        foreach (var edge in edges) {
            var vertices = new[] { edge.V0, edge.V1 };
            if (!String.IsNullOrEmpty(msg)) msg += "; ";
            msg += S2TextFormat.ToDebugString(vertices);
        }
        return msg;
    }

    private readonly GraphOptions graph_options_;
    private readonly S2ShapeIndex index_;
    private readonly int dimension_;
}
