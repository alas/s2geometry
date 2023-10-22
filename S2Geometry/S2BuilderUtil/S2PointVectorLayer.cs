using static S2Geometry.S2Builder;

namespace S2Geometry.S2BuilderUtil;

using Options = S2PointVectorLayer.Options;
using DuplicateEdges = GraphOptions.DuplicateEdges;
using DegenerateEdges = GraphOptions.DegenerateEdges;
using SiblingPairs = GraphOptions.SiblingPairs;

// A layer type that collects degenerate edges as points.
// This layer expects all edges to be degenerate. In case of finding
// non-degenerate edges it sets S2Error but it still generates the
// output with degenerate edges.
public class S2PointVectorLayer(List<S2Point> points, LabelSet? label_set_ids, IdSetLexicon? label_set_lexicon, S2PointVectorLayer.Options? options = null) : Layer
{
    public class Options
    {
        public Options() { DuplicateEdges_ = DuplicateEdges.MERGE; }
        public Options(DuplicateEdges duplicate_edges) { DuplicateEdges_ = duplicate_edges; }

        // DEFAULT: DuplicateEdges.MERGE
        public DuplicateEdges DuplicateEdges_ { get; set; }
    }

    public S2PointVectorLayer(List<S2Point> points, Options? options = null) : this(points, null, null, options ?? new Options()) { }

    // Layer interface:
    public override GraphOptions GraphOptions_()
    {
        return new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP, options_.DuplicateEdges_, SiblingPairs.KEEP);
    }
    public override void Build(Graph g, out S2Error error)
    {
        error = S2Error.OK;
        Graph.LabelFetcher fetcher = new(g, EdgeType.DIRECTED);

        var labels = new InputEdgeLoop();  // Temporary storage for labels.
        for (var edge_id = 0; edge_id < g.Edges.Count; edge_id++)
        {
            var edge = g.GetEdge(edge_id);
            if (edge.ShapeId != edge.EdgeId)
            {
                error = new(S2ErrorCode.INVALID_ARGUMENT, "Found non-degenerate edges");
                continue;
            }
            points_.Add(g.Vertex(edge.ShapeId));
            if (label_set_ids_ is not null)
            {
                fetcher.Fetch(edge_id, labels);
                label_set_ids_.Add(label_set_lexicon_.Add(labels));
            }
        }
    }

    private readonly List<S2Point> points_ = points;
    private readonly LabelSet? label_set_ids_ = label_set_ids;
    private readonly IdSetLexicon? label_set_lexicon_ = label_set_lexicon;
    private readonly Options options_ = options ?? new Options();
}

// Like S2PointVectorLayer, but adds the points to a MutableS2ShapeIndex (if
// the point vector is non-empty).
public class IndexedS2PointVectorLayer : Layer
{
    public IndexedS2PointVectorLayer(MutableS2ShapeIndex index, Options? options = null)
    { index_ = index; layer_ = new S2PointVectorLayer(points_, options ?? new Options()); }

    public override GraphOptions GraphOptions_()
    {
        return layer_.GraphOptions_();
    }

    public override void Build(Graph g, out S2Error error)
    {
        layer_.Build(g, out error);
        if (error.IsOk() && points_.Count!=0)
        {
            index_.Add(new S2PointVectorShape([.. points_]));
        }
    }

    private readonly MutableS2ShapeIndex index_;
    private readonly List<S2Point> points_ = [];
    private readonly S2PointVectorLayer layer_;
}
