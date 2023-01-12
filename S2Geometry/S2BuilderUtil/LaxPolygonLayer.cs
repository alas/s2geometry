// Note that there are two supported output types for polygons: S2Polygon and
// S2LaxPolygonShape.  Use S2Polygon if you need the full range of operations
// that S2Polygon implements.  Use S2LaxPolygonShape if you want to represent
// polygons with zero-area degenerate regions, or if you need a type that has
// low memory overhead and fast initialization.  However, be aware that to
// convert from S2LaxPolygonShape to S2Polygon you will need to use S2Builder
// again.
//
// Similarly, there are two supported output formats for polygon meshes:
// S2PolygonMesh and S2LaxPolygonShapeVector.  Use S2PolygonMesh if you need
// to be able to determine which polygons are adjacent to each edge or vertex;
// otherwise use S2LaxPolygonShapeVector, which uses less memory and is faster
// to construct.

// A layer type that assembles edges (directed or undirected) into an
// S2LaxPolygonShape.  Returns an error if the edges cannot be assembled into
// loops.
//
// If the input edges are directed, they must be oriented such that the
// polygon interior is to the left of all edges.  Directed edges are always
// preferred (see S2Builder.EdgeType).
//
// LaxPolygonLayer is implemented such that if the input to S2Builder is a
// polygon and is not modified, then the output has the same cyclic ordering
// of loop vertices and the same loop ordering as the input polygon.
//
// If the given edge graph is degenerate (i.e., it consists entirely of
// degenerate edges and sibling pairs), then the IsFullPolygonPredicate
// associated with the edge graph is called to determine whether the output
// polygon should be empty (possibly with degenerate shells) or full (possibly
// with degenerate holes).  This predicate can be specified as part of the
// S2Builder input geometry.

namespace S2Geometry.S2BuilderUtil;

using static S2Builder;
using static S2Builder.GraphOptions;

public class LaxPolygonLayer : Layer
{
    public class Options
    {
        // Constructor that uses the default options (listed below).
        public Options() : this(EdgeType.DIRECTED) { }

        // Constructor that specifies the edge type.
        public Options(EdgeType edge_type)
        {
            EdgeType = edge_type;
            DegenerateBoundaries_ = DegenerateBoundaries.KEEP;
        }

        // Indicates whether the input edges provided to S2Builder are directed or
        // undirected.  Directed edges should be used whenever possible (see
        // S2Builder.EdgeType for details).
        //
        // If the input edges are directed, they should be oriented so that the
        // polygon interior is to the left of all edges.  This means that for a
        // polygon with holes, the outer loops ("shells") should be directed
        // counter-clockwise while the inner loops ("holes") should be directed
        // clockwise.  Note that S2Builder.AddPolygon() does this automatically.
        //
        // DEFAULT: S2Builder.EdgeType.DIRECTED
        public EdgeType EdgeType { get; set; }

        // Specifies whether degenerate boundaries should be discarded or kept.
        // (A degenerate boundary consists of either a sibling edge pair or an
        // edge from a vertex to itself.)  Optionally, degenerate boundaries may
        // be kept only if they represent shells, or only if they represent holes.
        //
        // This option is useful for normalizing polygons with various boundary
        // conditions.  For example, DISCARD_HOLES can be used to normalize closed
        // polygons (those that include their boundary), since degenerate holes do
        // not affect the set of points contained by such polygons.  Similarly,
        // DISCARD_SHELLS can be used to normalize polygons with open boundaries.
        // DISCARD is used to normalize polygons with semi-open boundaries (since
        // degenerate loops do not affect point containment in that case), and
        // finally KEEP is useful for working with any type of polygon where
        // degeneracies are assumed to contain an infinitesmal interior.  (This
        // last model is the most useful for working with simplified geometry,
        // since it maintains the closest fidelity to the original geometry.)
        //
        // DEFAULT: DegenerateBoundaries.KEEP
        public enum DegenerateBoundaries : byte
        {
            DISCARD, DISCARD_HOLES, DISCARD_SHELLS, KEEP
        }
        public DegenerateBoundaries DegenerateBoundaries_ { get; set; }
    }

    // Specifies that a polygon should be constructed using the given options.
    public LaxPolygonLayer(S2LaxPolygonShape polygon, Options? options = null)
    {
        Init(polygon, null, null, options ?? new Options());
    }

    // Specifies that a polygon should be constructed using the given options,
    // and that any labels attached to the input edges should be returned in
    // "label_set_ids" and "label_set_lexicion".
    //
    // The labels associated with the edge "polygon.chain_edge(i, j)"
    // can be retrieved as follows:
    //
    //   foreach (Int32 label in label_set_lexicon.id_set(label_set_ids[i][j])) {...}
    public LaxPolygonLayer(S2LaxPolygonShape polygon, LabelSetIds label_set_ids,
        IdSetLexicon label_set_lexicon, Options? options = null)
    {
        Init(polygon, label_set_ids, label_set_lexicon, options ?? new Options());
    }

    // Layer interface:
    public override GraphOptions GraphOptions_()
    {
        if (options_.DegenerateBoundaries_ == Options.DegenerateBoundaries.DISCARD)
        {
            // There should not be any duplicate edges, but if there are then we keep
            // them since this yields more comprehensible error messages.
            return new GraphOptions(options_.EdgeType, DegenerateEdges.DISCARD,
                                DuplicateEdges.KEEP, SiblingPairs.DISCARD);
        }
        else
        {
            // Keep at most one copy of each sibling pair and each isolated vertex.
            return new GraphOptions(options_.EdgeType, DegenerateEdges.DISCARD_EXCESS,
                                DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
        }
    }
    public override void Build(Graph g, out S2Error error)
    {
        if (label_set_ids_ is not null) label_set_ids_.Clear();
        if (g.Options.EdgeType_ == EdgeType.DIRECTED)
        {
            BuildDirected(g, out error);
        }
        else
        {
            error = new(S2ErrorCode.UNIMPLEMENTED, "Undirected edges not supported yet");
        }
    }

    private void Init(S2LaxPolygonShape polygon, LabelSetIds? label_set_ids, IdSetLexicon? label_set_lexicon, Options options)
    {
        MyDebug.Assert((label_set_ids is null) == (label_set_lexicon == null));
        polygon_ = polygon;
        label_set_ids_ = label_set_ids;
        label_set_lexicon_ = label_set_lexicon;
        options_ = options;
    }
    private static void AppendPolygonLoops(Graph g, List<EdgeLoop> edge_loops, List<List<S2Point>> loops)
    {
        foreach (var edge_loop in edge_loops)
        {
            var vertices = new List<S2Point>(edge_loop.Count);
            foreach (var edge_id in edge_loop)
            {
                vertices.Add(g.Vertex(g.GetEdge(edge_id).ShapeId));
            }
            loops.Add(vertices);
        }
    }
    private void AppendEdgeLabels(Graph g, List<EdgeLoop> edge_loops)
    {
        if (label_set_ids_ is null) return;

        List<Int32> labels = new();  // Temporary storage for labels.
        Graph.LabelFetcher fetcher = new(g, options_.EdgeType);
        foreach (var edge_loop in edge_loops)
        {
            var loop_label_set_ids = new List<Int32>(edge_loop.Count);
            foreach (var edge_id in edge_loop)
            {
                fetcher.Fetch(edge_id, labels);
                loop_label_set_ids.Add(label_set_lexicon_.Add(labels));
            }
            label_set_ids_.Add(loop_label_set_ids);
        }
    }
    private void BuildDirected(Graph g, out S2Error error)
    {
        // Some cases are implemented by constructing a new graph with certain
        // degenerate edges removed (overwriting "g").  "new_edges" is where the
        // edges for the new graph are stored.
        var new_edges = new List<Edge>();
        var new_input_edge_id_set_ids = new List<Int32>();
        var loops = new List<List<S2Point>>();
        var degenerate_boundaries = options_.DegenerateBoundaries_;
        if (degenerate_boundaries == Options.DegenerateBoundaries.DISCARD)
        {
            // This is the easiest case, since there are no degeneracies.
            if (g.NumEdges == 0) MaybeAddFullLoop(g, loops, out _);
        }
        else if (degenerate_boundaries == Options.DegenerateBoundaries.KEEP)
        {
            // S2LaxPolygonShape doesn't need to distinguish degenerate shells from
            // holes except when the entire graph is degenerate, in which case we need
            // to decide whether it represents an empty polygons possibly with
            // degenerate shells, or a full polygon possibly with degenerate holes.
            if (PolygonDegeneracy.IsFullyDegenerate(g))
            {
                MaybeAddFullLoop(g, loops, out _);
            }
        }
        else
        {
            // For DISCARD_SHELLS and DISCARD_HOLES we first determine whether any
            // degenerate loops of the given type exist, and if so we construct a new
            // graph with those edges removed (overwriting "g").
            bool discard_holes =
                (degenerate_boundaries == Options.DegenerateBoundaries.DISCARD_HOLES);
            var degeneracies = PolygonDegeneracy.FindPolygonDegeneracies(g, out error);
            if (!error.IsOk()) return;
            if (degeneracies.Count == g.NumEdges)
            {
                if (!degeneracies.Any())
                {
                    MaybeAddFullLoop(g, loops, out error);
                }
                else if (degeneracies[0].IsHole)
                {
                    loops.Add(new List<S2Point>());  // Full loop.
                }
            }
            var edges_to_discard = new List<Int32>();
            foreach (var degeneracy in degeneracies)
            {
                if (degeneracy.IsHole == discard_holes)
                {
                    edges_to_discard.Add((int)degeneracy.EdgeId);
                }
            }
            if (edges_to_discard.Any())
            {
                // Construct a new graph that discards the unwanted edges.
                edges_to_discard.Sort();
                DiscardEdges(g, edges_to_discard, new_edges, new_input_edge_id_set_ids);
                g = new Graph(g.Options, g.Vertices,
                          new_edges, new_input_edge_id_set_ids,
                          g.InputEdgeIdSetLexicon, g.LabelSetIds,
                          g.LabelSetLexicon, g.IsFullPolygonPredicate());
            }
        }
        var edge_loops = new List<EdgeLoop>();
        if (!g.GetDirectedLoops(Graph.LoopType.CIRCUIT, edge_loops, out error))
        {
            return;
        }
        AppendPolygonLoops(g, edge_loops, loops);
        AppendEdgeLabels(g, edge_loops);
        edge_loops.Clear();  // Release memory
        new_edges.Clear();
        new_input_edge_id_set_ids.Clear();
        polygon_.Init(loops);
    }


    // Returns all edges of "g" except for those identified by "edges_to_discard".
    private static void DiscardEdges(Graph g, List<Int32> edges_to_discard, List<Edge> new_edges, List<Int32> new_input_edge_id_set_ids)
    {
        MyDebug.Assert(edges_to_discard.IsSorted());
        new_edges.Clear();
        new_input_edge_id_set_ids.Clear();
        new_edges.Capacity = g.NumEdges;
        new_input_edge_id_set_ids.Capacity = g.NumEdges;
        var it = 0;
        for (int e = 0; e < g.NumEdges; ++e)
        {
            if (it < edges_to_discard.Count && e == edges_to_discard[it])
            {
                ++it;
            }
            else
            {
                new_edges.Add(g.GetEdge(e));
                new_input_edge_id_set_ids.Add(g.InputEdgeIdSetId(e));
            }
        }
        MyDebug.Assert(it == edges_to_discard.Count);
    }

    private static void MaybeAddFullLoop(Graph g, List<List<S2Point>> loops, out S2Error error)
    {
        if (g.IsFullPolygon(out error))
        {
            loops.Add(new List<S2Point>());  // Full loop.
        }
    }

    private S2LaxPolygonShape polygon_;
    private LabelSetIds? label_set_ids_;
    private IdSetLexicon? label_set_lexicon_;
    private Options options_;
}

// Like LaxPolygonLayer, but adds the polygon to a MutableS2ShapeIndex (if the
// polygon is non-empty).
public class IndexedLaxPolygonLayer : Layer
{
    public IndexedLaxPolygonLayer(MutableS2ShapeIndex index, LaxPolygonLayer.Options? options = null)
    {
        index_ = index; polygon_ = new S2LaxPolygonShape(); layer_ = new LaxPolygonLayer(polygon_, options ?? new LaxPolygonLayer.Options());
    }

    public override GraphOptions GraphOptions_()
    {
        return layer_.GraphOptions_();
    }

    public override void Build(Graph g, out S2Error error)
    {
        layer_.Build(g, out error);
        if (error.IsOk() && !polygon_.IsEmpty())
        {
            index_.Add(polygon_);
        }
    }

    private readonly MutableS2ShapeIndex index_;
    private readonly S2LaxPolygonShape polygon_;
    private readonly LaxPolygonLayer layer_;
}
