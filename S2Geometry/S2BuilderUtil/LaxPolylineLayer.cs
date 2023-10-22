// A layer type that assembles edges (directed or undirected) into an
// S2LaxPolylineShape.  Returns an error if the edges cannot be assembled into
// a single unbroken polyline.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of input edges whenever
// possible, so that if the input is a polyline and it is not modified by
// S2Builder, then the output will be the same polyline (even if the polyline
// backtracks on itself or forms a loop).
//
// LaxPolylineLayer does not support options such as discarding sibling pairs
// or merging duplicate edges because these options can split the polyline
// into several pieces.  TODO(ericv): Implement LaxPolylineVectorLayer.

namespace S2Geometry.S2BuilderUtil;

// Specifies that a polyline should be constructed using the given options,
// and that any labels attached to the input edges should be returned in
// "label_set_ids" and "label_set_lexicion".
//
// The labels associated with the edge "polyline.vertex({j, j+1})" can be
// retrieved as follows:
//
//   for (int32 label : label_set_lexicon.id_set(label_set_ids[j])) {...}
using LabelSetIds = List<System.Int32>;

using DegenerateEdges = S2Builder.GraphOptions.DegenerateEdges;
using DuplicateEdges = S2Builder.GraphOptions.DuplicateEdges;
using SiblingPairs = S2Builder.GraphOptions.SiblingPairs;
using PolylineType = S2Builder.Graph.PolylineType;

public class LaxPolylineLayer : S2Builder.Layer
{

    // Specifies that a polyline should be constructed using the given options.
    public LaxPolylineLayer(S2LaxPolylineShape polyline, Options? options = null)
        : this(polyline, null, null, options)
    {
    }
    public LaxPolylineLayer(S2LaxPolylineShape polyline, LabelSetIds? label_set_ids,
                   IdSetLexicon? label_set_lexicon, Options? options = null)
    {
        MyDebug.Assert((label_set_ids is null) == (label_set_lexicon is null));
        polyline_ = polyline;
        label_set_ids_ = label_set_ids;
        label_set_lexicon_ = label_set_lexicon;
        options_ = options ?? new Options();
    }

    // Layer interface:
    public override S2Builder.GraphOptions GraphOptions_()
    {
        return new S2Builder.GraphOptions(options_.EdgeType_, DegenerateEdges.KEEP,
                            DuplicateEdges.KEEP, SiblingPairs.KEEP);
    }

    public override void Build(S2Builder.Graph g, out S2Error error)
    {
        error = S2Error.OK;
        if (g.NumEdges == 0)
        {
            polyline_ = new(Array.Empty<S2Point>());
            return;
        }
        var edge_polylines =
            g.GetPolylines(PolylineType.WALK);
        if (edge_polylines.Count != 1)
        {
            error = new S2Error(S2ErrorCode.BUILDER_EDGES_DO_NOT_FORM_POLYLINE,
                        "Input edges cannot be assembled into polyline");
            return;
        }
        var edge_polyline = edge_polylines[0];
        List<S2Point> vertices = [];  // Temporary storage for vertices.
        vertices.Capacity = edge_polyline.Count;
        vertices.Add(g.Vertex(g.GetEdge(edge_polyline[0]).ShapeId));
        foreach (EdgeId e in edge_polyline)
        {
            vertices.Add(g.Vertex(g.GetEdge(e).EdgeId));
        }
        if (label_set_ids_ is not null)
        {
            S2Builder.Graph.LabelFetcher fetcher = new(g, options_.EdgeType_);
            List<Label> labels = [];  // Temporary storage for labels.
            label_set_ids_.Capacity = edge_polyline.Count;
            foreach (EdgeId e in edge_polyline)
            {
                fetcher.Fetch(e, labels);
                label_set_ids_.Add(label_set_lexicon_!.Add(labels));
            }
        }
        polyline_ = new(vertices);
    }

    private S2LaxPolylineShape? polyline_;
    private readonly LabelSetIds? label_set_ids_;
    private readonly IdSetLexicon? label_set_lexicon_;
    private readonly Options options_;

    public class Options
    {
        // Constructor that uses the default options (listed below).
        public Options()
        {
            EdgeType_ = S2Builder.EdgeType.DIRECTED;
        }

        // Constructor that specifies the edge type.
        public Options(S2Builder.EdgeType edge_type)
        {
            EdgeType_ = edge_type;
        }

        // Indicates whether the input edges provided to S2Builder are directed or
        // undirected.  Directed edges should be used whenever possible to avoid
        // ambiguity.
        //
        // DEFAULT: S2Builder::EdgeType::DIRECTED
        public S2Builder.EdgeType EdgeType_ { get; set; }
    }

    // Like LaxPolylineLayer, but adds the polyline to a MutableS2ShapeIndex (if
    // the polyline is non-empty).
    public class IndexedLaxPolylineLayer : S2Builder.Layer
    {
        public IndexedLaxPolylineLayer(MutableS2ShapeIndex index,
                                          Options? options = null)
        {
            index_ = index;
            polyline_ = new S2LaxPolylineShape();
            layer_ = new(polyline_, options ?? new Options());
        }

        public override S2Builder.GraphOptions GraphOptions_()
        {
            return layer_.GraphOptions_();
        }

        public override void Build(S2Builder.Graph g, out S2Error error)
        {
            layer_.Build(g, out error);
            if (error.IsOk() && polyline_.NumVertices() > 0)
            {
                index_.Add(polyline_);
            }
        }

        private readonly MutableS2ShapeIndex index_;
        private readonly S2LaxPolylineShape polyline_;
        private readonly LaxPolylineLayer layer_;
    }
}
