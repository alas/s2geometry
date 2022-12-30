// A layer type that assembles edges (directed or undirected) into an
// S2Polyline.  Returns an error if the edges cannot be assembled into a
// single unbroken polyline.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of input edges whenever
// possible, so that if the input is a polyline and it is not modified by
// S2Builder, then the output will be the same polyline (even if the polyline
// backtracks on itself or forms a loop).
//
// S2PolylineLayer does not support options such as discarding sibling pairs
// or merging duplicate edges because these options can split the polyline
// into several pieces.  Use S2PolylineVectorLayer if you need these features.

namespace S2Geometry.S2BuilderUtil;

using static S2Builder;
using static S2Builder.GraphOptions;

public class S2PolylineLayer : Layer
{
    public class Options
    {
        // Constructor that uses the default options (listed below).
        public Options()
        {
            EdgeType = EdgeType.DIRECTED; Validate = false;
        }

        // Constructor that specifies the edge type.
        public Options(EdgeType edge_type)
        {
            EdgeType = edge_type; Validate = false;
        }

        // Indicates whether the input edges provided to S2Builder are directed or
        // undirected.  Directed edges should be used whenever possible to avoid
        // ambiguity.
        //
        // DEFAULT: S2Builder.EdgeType.DIRECTED
        public EdgeType EdgeType { get; set; }

        // If true, calls FindValidationError() on the output polyline.  If any
        // error is found, it will be returned by S2Builder.Build().
        //
        // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
        // order to turn off the default error checking in debug builds.
        //
        // DEFAULT: false
        public bool Validate { get; set; }
    }

    // Specifies that a polyline should be constructed using the given options.
    public S2PolylineLayer(S2Polyline polyline, Options? options = null)
    {
        Init(polyline, null, null, options ?? new Options());
    }

    // Specifies that a polyline should be constructed using the given options,
    // and that any labels attached to the input edges should be returned in
    // "label_set_ids" and "label_set_lexicion".
    //
    // The labels associated with the edge "polyline.vertex({j, j+1})" can be
    // retrieved as follows:
    //
    //   for (Int32 label : label_set_lexicon.id_set(label_set_ids[j])) {...}
    public S2PolylineLayer(S2Polyline polyline, LabelSet label_set_ids, IdSetLexicon label_set_lexicon, Options? options = null)
    {
        Init(polyline, label_set_ids, label_set_lexicon, options ?? new Options());
    }

    // Layer interface:
    public override GraphOptions GraphOptions_()
    {
        // Remove edges that collapse to a single vertex, but keep duplicate and
        // sibling edges, since merging duplicates or discarding siblings can make
        // it impossible to assemble the edges into a single polyline.
        return new GraphOptions(options_.EdgeType, DegenerateEdges.DISCARD, DuplicateEdges.KEEP, SiblingPairs.KEEP);
    }
    public override void Build(Graph g, out S2Error error)
    {
        error = S2Error.OK;
        if (g.NumEdges == 0)
        {
            polyline_ = new S2Polyline(Array.Empty<S2Point>());
            return;
        }
        var edge_polylines = g.GetPolylines(Graph.PolylineType.WALK);
        if (edge_polylines.Count != 1)
        {
            error = new(S2ErrorCode.BUILDER_EDGES_DO_NOT_FORM_POLYLINE, "Input edges cannot be assembled into polyline");
            return;
        }
        var edge_polyline = edge_polylines[0];
        // Temporary storage for vertices.
        var vertices = new List<S2Point>(edge_polyline.Count)
            {
                g.Vertex(g.GetEdge(edge_polyline[0]).ShapeId)
            };
        foreach (var e in edge_polyline)
        {
            vertices.Add(g.Vertex(g.GetEdge(e).EdgeId));
        }
        if (label_set_ids_ is not null)
        {
            var fetcher = new Graph.LabelFetcher(g, options_.EdgeType);
            var labels = new List<Label>();  // Temporary storage for labels.
            label_set_ids_.Capacity = edge_polyline.Count;
            foreach (var e in edge_polyline)
            {
                fetcher.Fetch(e, labels);
                label_set_ids_.Add(label_set_lexicon_.Add(labels));
            }
        }
        polyline_ = new S2Polyline(vertices.ToArray());
        if (options_.Validate)
        {
            polyline_.FindValidationError(out error);
        }
    }

    private void Init(S2Polyline polyline, LabelSet? label_set_ids, IdSetLexicon? label_set_lexicon, Options options)
    {
        System.Diagnostics.Debug.Assert((label_set_ids is null) == (label_set_lexicon is null));
        polyline_ = polyline;
        label_set_ids_ = label_set_ids;
        label_set_lexicon_ = label_set_lexicon;
        options_ = options;

        if (options_.Validate)
        {
            polyline_.s2debug_override_ = S2Debug.DISABLE;
        }
    }

    private S2Polyline polyline_;
    private LabelSet? label_set_ids_;
    private IdSetLexicon? label_set_lexicon_;
    private Options options_;
}

// Like S2PolylineLayer, but adds the polyline to a MutableS2ShapeIndex (if the
// polyline is non-empty).
public class IndexedS2PolylineLayer : Layer
{

    public IndexedS2PolylineLayer(MutableS2ShapeIndex index, S2PolylineLayer.Options? options = null)
    {
        index_ = index; polyline_ = new S2Polyline();
        layer_ = new S2PolylineLayer(polyline_, options ?? new S2PolylineLayer.Options());
    }

    public override GraphOptions GraphOptions_()
    {
        return layer_.GraphOptions_();
    }

    public override void Build(Graph g, out S2Error error)
    {
        layer_.Build(g, out error);
        if (error.IsOk() && polyline_.NumVertices() > 0)
        {
            index_.Add(new S2Polyline.OwningShape(polyline_));
        }
    }

    private readonly MutableS2ShapeIndex index_;
    private readonly S2Polyline polyline_;
    private readonly S2PolylineLayer layer_;
}
