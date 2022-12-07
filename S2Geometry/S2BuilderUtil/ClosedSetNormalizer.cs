namespace S2Geometry.S2BuilderUtil;

using static S2Builder;
using static S2Builder.GraphOptions;

// The purpose of this class is to allow S2Builder.Layer implementations to
// remove polygon and polyline degeneracies by converting them to polylines or
// points.  Note that most clients should not use ClosedSetNormalizer itself,
// but should instead call NormalizeClosedSet defined below.
//
// A polyline degeneracy is a polyline consisting of a single degenerate edge.
// A polygon degeneracy is either a single-vertex loop (a degenerate edge from
// a vertex to itself) or a sibling edge pair (consisting of an edge and its
// corresponding reverse edge).  Polygon degeneracies are further classified
// as shells or holes depending on whether they are located in the exterior or
// interior of the polygon respectively.  For example, a single-vertex loop
// contained within a polygon shell would be classified as a hole.
//
// All objects are modeled as closed, i.e. polygons contain their boundaries
// and polylines contain their endpoints.  Note that under this model,
// degenerate polygon shells and holes need to be handled differently.
// Degenerate shells are converted to polylines or points, whereas degenerate
// holes do not affect the set of points contained by the polygon and are
// simply discarded.
//
// Specifically, given three S2Builder.Graphs (corresponding to points,
// polylines, and polygons), this class makes the following transformations:
//
//  - Polygon sibling edge pairs are either discarded (for holes) or converted
//    to a pair of polyline edges (for shells).
//  - Degenerate polygon edges are either discarded (for holes) or converted
//    to points (for shells).
//  - Degenerate polyline edges are converted to points.
//
// Optionally, this class further normalize the graphs by suppressing edges
// that are duplicates of higher-dimensional edges.  In other words:
//
//  - Polyline edges that coincide with polygon edges are discarded.
//  - Points that coincide with polyline or polygon vertices are discarded.
//
// (When edges are discarded, any labels attached to those edges are discarded
// as well.)
//
// This class takes three graphs as input and yields three graphs as output.
// However note that the output graphs are *not* independent objects; they may
// point to data in the input graphs or data owned by the ClosedSetNormalizer
// itself.  For this reason the input graphs and ClosedSetNormalizer must
// persist until the output graphs are no longer needed.
//
// Finally, note that although this class may be necessary in some situations
// (e.g., to implement the OGC Simple Features Access spec), in general the
// recommended approach to degeneracies is simply to keep them (by using a
// representation such as S2LaxPolygonShape or S2LaxPolylineShape).
// Keeping degeneracies has many advantages, such as not needing to deal with
// geometry of multiple dimensions, and being able to preserve polygon
// boundaries accurately (including degenerate holes).
public class ClosedSetNormalizer
{
    public class Options
    {
        public Options() => SuppressLowerDimensions = true;

        // If "suppress_lower_dimensions" is true, then the graphs are further
        // normalized by discarding lower-dimensional edges that coincide with
        // higher-dimensional edges.
        //
        // DEFAULT: true
        public bool SuppressLowerDimensions { get; set; }
    }

    // Constructs a ClosedSetNormalizer whose output will be three
    // S2Builder.Graphs with the given "graph_options_out".
    //
    // REQUIRES: graph_options_out.size() == 3
    // REQUIRES: graph_options_out[0].edge_type() == DIRECTED
    // REQUIRES: graph_options_out[1].sibling_pairs() != {CREATE, REQUIRE}
    // REQUIRES: graph_options_out[2].edge_type() == DIRECTED
    public ClosedSetNormalizer(Options options, List<GraphOptions> graph_options_out)
    {
        Options_ = options;
        graph_options_out_ = graph_options_out;
        GraphOptions_ = graph_options_out_;
        sentinel_ = new Edge(Int32.MaxValue, Int32.MaxValue);
        System.Diagnostics.Debug.Assert(graph_options_out_.Count == 3);
        System.Diagnostics.Debug.Assert(graph_options_out_[0].EdgeType_ == EdgeType.DIRECTED);
        System.Diagnostics.Debug.Assert(graph_options_out_[2].EdgeType_ == EdgeType.DIRECTED);

        // NOTE(ericv): Supporting these options would require some extra code in
        // order to handle undirected edges, and they are not useful for building
        // polylines anyway (they are intended for polygon meshes).
        System.Diagnostics.Debug.Assert(graph_options_out_[1].SiblingPairs_ != SiblingPairs.CREATE);
        System.Diagnostics.Debug.Assert(graph_options_out_[1].SiblingPairs_ != SiblingPairs.REQUIRE);

        // Set the GraphOptions for the input graphs to ensure that (1) they share a
        // common set of vertices, (2) degenerate edges are kept only if they are
        // isolated, and (3) multiple copies of siblings pairs are discarded.  (Note
        // that there may be multiple copies of isolated degenerate edges; clients
        // can eliminate them if desired using DuplicateEdges.MERGE.)
        for (var dim = 0; dim < 3; ++dim)
        {
            GraphOptions_[dim].AllowVertexFiltering = false;
        }
        GraphOptions_[1].DegenerateEdges_ = DegenerateEdges.DISCARD_EXCESS;
        GraphOptions_[2].DegenerateEdges_ = DegenerateEdges.DISCARD_EXCESS;
        GraphOptions_[2].SiblingPairs_ = SiblingPairs.DISCARD_EXCESS;
    }

    // Returns the ClosedSetNormalizer options.
    public Options Options_ { get; }

    // Returns the GraphOptions that should be used to construct the input
    // S2Builder.Graph of each dimension.
    public List<GraphOptions> GraphOptions_ { get; }

    // Normalizes the input graphs and returns a new set of graphs where
    // degeneracies have been discarded or converted to objects of lower
    // dimension.  input[d] is the graph representing edges of dimension "d".
    //
    // Note that the input graphs, their contents, and the ClosedSetNormalizer
    // itself must persist until the output of this class is no longer needed.
    // (To emphasize this requirement, a reference is returned.)
    public List<Graph> Run(List<Graph> g, out S2Error error)
    {
        // Ensure that the input graphs were built with our requested options.
        for (var dim = 0; dim < 3; ++dim)
        {
            System.Diagnostics.Debug.Assert(g[dim].Options == GraphOptions_[dim]);
        }

        if (Options_.SuppressLowerDimensions)
        {
            // Build the auxiliary data needed to suppress lower-dimensional edges.
            in_edges2_ = g[2].GetInEdgeIds();
            is_suppressed_.Capacity = g[0].Vertices.Count;
            for (var dim = 1; dim <= 2; ++dim)
            {
                for (var e = 0; e < g[dim].NumEdges; ++e)
                {
                    var edge = g[dim].GetEdge(e);
                    if (edge.ShapeId != edge.EdgeId)
                    {
                        is_suppressed_[edge.ShapeId] = true;
                        is_suppressed_[edge.EdgeId] = true;
                    }
                }
            }
        }

        // Compute the edges that belong in the output graphs.
        NormalizeEdges(g, out error);

        // If any edges were added or removed, we need to run Graph.ProcessEdges to
        // ensure that the edges satisfy the requested GraphOptions.  Note that
        // since edges are never added to dimension 2, we can use the edge count to
        // test whether any edges were removed.  If no edges were removed from
        // dimension 2, then no edges were added to dimension 1, and so we can again
        // use the edge count to test whether any edges were removed, etc.
        var modified = new bool[3];
        var any_modified = false;
        for (var dim = 2; dim >= 0; --dim)
        {
            if (new_edges_[dim].Count != g[dim].NumEdges) any_modified = true;
            modified[dim] = any_modified;
        }

        if (!any_modified)
        {
            for (var dim = 0; dim < 3; ++dim)
            {
                // Copy the graphs to ensure that they have the GraphOptions that were
                // originally requested.
                new_graphs_.Add(new Graph(
                    graph_options_out_[dim], g[dim].Vertices, g[dim].Edges,
                    g[dim].InputEdgeIdSetIds, g[dim].InputEdgeIdSetLexicon,
                    g[dim].LabelSetIds, g[dim].LabelSetLexicon,
                    g[dim].IsFullPolygonPredicate()));
            }
        }
        else
        {
            // Make a copy of input_edge_id_set_lexicon() so that ProcessEdges can
            // merge edges if necessary.
            new_input_edge_id_set_lexicon_ = g[0].InputEdgeIdSetLexicon;
            for (var dim = 0; dim < 3; ++dim)
            {
                if (modified[dim])
                {
                    Graph.ProcessEdges(graph_options_out_[dim], new_edges_[dim],
                        new_input_edge_ids_[dim],
                        new_input_edge_id_set_lexicon_, out error);
                }
                new_graphs_.Add(new Graph(
                    graph_options_out_[dim], g[dim].Vertices, new_edges_[dim],
                    new_input_edge_ids_[dim], new_input_edge_id_set_lexicon_,
                    g[dim].LabelSetIds, g[dim].LabelSetLexicon,
                    g[dim].IsFullPolygonPredicate()));
            }
        }

        return new_graphs_;
    }

    // Helper function that advances to the next edge in the given graph,
    // returning a sentinel value once all edges are exhausted.
    private Edge Advance(Graph g, ref Int32 e)
    {
        return (++e == g.NumEdges) ? sentinel_ : g.GetEdge(e);
    }

    // Helper function that advances to the next incoming edge in the given graph,
    // returning a sentinel value once all edges are exhausted.
    private Edge AdvanceIncoming(Graph g, List<Int32> in_edges, ref int i)
    {
        return ((in_edges[++i] == in_edges.Count) ? sentinel_ :
            Graph.Reverse(g.GetEdge(in_edges[i])));
    }
    private void NormalizeEdges(List<Graph> g, out S2Error error)
    {
        // Find the degenerate polygon edges and sibling pairs, and classify each
        // edge as belonging to either a shell or a hole.
        var degeneracies = PolygonDegeneracy.FindPolygonDegeneracies(g[2], out error);
        var degeneracy = 0;

        // Walk through the three edge vectors performing a merge join.  We also
        // maintain positions in two other auxiliary vectors: the vector of sorted
        // polygon degeneracies (degeneracies), and the vector of incoming polygon
        // edges (if we are suppressing lower-dimensional duplicate edges).
        var e0 = -1;
        var e1 = -1;
        var e2 = -1;  // Current position in g[dim].edges()
        var in_e2 = -1;  // Current position in in_edges2_
        var edge0 = Advance(g[0], ref e0);
        var edge1 = Advance(g[1], ref e1);
        var edge2 = Advance(g[2], ref e2);
        var in_edge2 = AdvanceIncoming(g[2], in_edges2_, ref in_e2);
        for (; ; )
        {
            if (edge2 <= edge1 && edge2 <= edge0)
            {
                if (edge2 == sentinel_) break;
                if (degeneracy == degeneracies.Count || degeneracies[degeneracy].EdgeId != e2)
                {
                    // Normal polygon edge (not part of a degeneracy).
                    AddEdge(2, g[2], e2);
                    while (Options_.SuppressLowerDimensions && edge1 == edge2)
                    {
                        edge1 = Advance(g[1], ref e1);
                    }
                }
                else if (!degeneracies[degeneracy++].IsHole)
                {
                    // Edge belongs to a degenerate shell.
                    if (edge2.ShapeId != edge2.EdgeId)
                    {
                        AddEdge(1, g[2], e2);
                        // Since this edge was demoted, make sure that it does not suppress
                        // any coincident polyline edge(s).
                        while (edge1 == edge2)
                        {
                            AddEdge(1, g[1], e1);
                            edge1 = Advance(g[1], ref e1);
                        }
                    }
                    else
                    {
                        // The test below is necessary because a single-vertex polygon shell
                        // can be discarded by a polyline edge incident to that vertex.
                        if (!IsSuppressed(edge2.ShapeId)) AddEdge(0, g[2], e2);
                    }
                }
                edge2 = Advance(g[2], ref e2);
            }
            else if (edge1 <= edge0)
            {
                if (edge1.ShapeId != edge1.EdgeId)
                {
                    // Non-degenerate polyline edge.  (Note that in_edges2_ is empty
                    // whenever "suppress_lower_dimensions" is false.)
                    while (in_edge2 < edge1)
                    {
                        in_edge2 = AdvanceIncoming(g[2], in_edges2_, ref in_e2);
                    }
                    if (edge1 != in_edge2) AddEdge(1, g[1], e1);
                }
                else
                {
                    // Degenerate polyline edge.
                    if (!IsSuppressed(edge1.ShapeId)) AddEdge(0, g[1], e1);
                    if (g[1].Options.EdgeType_ == EdgeType.UNDIRECTED) ++e1;
                }
                edge1 = Advance(g[1], ref e1);
            }
            else
            {
                // Input point.
                if (!IsSuppressed(edge0.ShapeId)) AddEdge(0, g[0], e0);
                edge0 = Advance(g[0], ref e0);
            }
        }
    }
    private void AddEdge(int new_dim, Graph g, Int32 e)
    {
        new_edges_[new_dim].Add(g.GetEdge(e));
        new_input_edge_ids_[new_dim].Add(g.InputEdgeIdSetId(e));
    }
    private bool IsSuppressed(Int32 v)
    {
        return Options_.SuppressLowerDimensions && is_suppressed_[v];
    }

    // Requested options for the output graphs.
    private readonly List<GraphOptions> graph_options_out_;

    // A sentinel value that compares larger than any valid edge.
    private readonly Edge sentinel_;
    // is_suppressed_[i] is true if vertex[i] belongs to a non-degenerate edge,
    // and therefore should be suppressed from the output graph for points.
    private readonly List<bool> is_suppressed_ = new();

    // A vector of incoming polygon edges sorted in lexicographic order.  This
    // is used to suppress directed polyline edges that match a polygon edge in
    // the reverse direction.
    private List<Int32> in_edges2_;

    // Output data.
    private readonly List<Graph> new_graphs_ = new();
    private readonly List<Edge>[] new_edges_ = new List<Edge>[3];
    private readonly List<Int32>[] new_input_edge_ids_ = new List<Int32>[3];
    private IdSetLexicon new_input_edge_id_set_lexicon_;
}


// This method implements the NormalizeClosedSet function.  The Create()
// method allocates a single object of this class whose ownership is shared
// (using shared_ptr) among the three returned S2Builder.Layers.  Here is how
// the process works:
//
//  - The returned layers are passed to a class (such as S2Builder or
//    S2BooleanOperation) that calls their Build methods.  We call these the
//    "input layers" because they provide the input to ClosedSetNormalizer.
//
//  - When Build() is called on the first two layers, pointers to the
//    corresponding Graph arguments are saved.
//
//  - When Build() is called on the third layer, ClosedSetNormalizer is used
//    to normalize the graphs, and then the Build() method of each of the
//    three output layers is called.
//
// TODO(ericv): Consider generalizing this technique as a public class.
public class NormalizeClosedSetImpl
{
    public static List<Layer> Create(List<Layer> output_layers, ClosedSetNormalizer.Options options)
    {
        var impl = new NormalizeClosedSetImpl(output_layers, options);
        var result = new List<Layer>();
        for (var dim = 0; dim < 3; ++dim)
        {
            result.Add(new DimensionLayer(dim, impl.normalizer_.GraphOptions_[dim], impl));
        }
        return result;
    }

    private NormalizeClosedSetImpl(List<Layer> output_layers, ClosedSetNormalizer.Options options)
    {
        output_layers_ = output_layers;
        normalizer_ = new ClosedSetNormalizer(options, new List<GraphOptions>
                {
                    output_layers_[0].GraphOptions_(),
                    output_layers_[1].GraphOptions_(),
                    output_layers_[2].GraphOptions_(),
                });
        graphs_ = new List<Graph?>() { null, null, null }; graphs_left_ = 3;
        System.Diagnostics.Debug.Assert(3 == output_layers_.Count);
    }

    private class DimensionLayer : Layer
    {
        public DimensionLayer(int dimension, GraphOptions graph_options, NormalizeClosedSetImpl impl)
        {
            dimension_ = dimension; graph_options_ = graph_options; impl_ = impl;
        }

        public override GraphOptions GraphOptions_() { return graph_options_; }
        private readonly GraphOptions graph_options_;

        public override void Build(Graph g, out S2Error error)
        {
            impl_.Build(dimension_, g, out error);
        }

        private readonly int dimension_;
        private readonly NormalizeClosedSetImpl impl_;
    }

    private void Build(int dimension, Graph g, out S2Error error)
    {
        error = S2Error.OK;
        // Errors are reported only on the last layer built.
        graphs_[dimension] = g;
        if (--graphs_left_ > 0) return;

        var output = normalizer_.Run(graphs_, out error);
        for (var dim = 0; dim < 3; ++dim)
        {
            output_layers_[dim].Build(output[dim], out error);
        }
    }

    private readonly List<Layer> output_layers_;
    private readonly ClosedSetNormalizer normalizer_;
    private readonly List<Graph?> graphs_;
    private int graphs_left_;
}

public static class NormalizeClosedSetX
{
    // Given a set of three output layers (one each for dimensions 0, 1, and 2),
    // returns a new set of layers that preprocess the input graphs using a
    // ClosedSetNormalizer with the given options.  This can be used to ensure
    // that the graphs passed to "output_layers" do not contain any polyline or
    // polygon degeneracies.
    //
    // Example showing how to compute the union of two S2ShapeIndexes containing
    // points, polylines, and/or polygons, and save the result as a collection of
    // S2Points, S2Polylines, and S2Polygons in another S2ShapeIndex (where
    // degeneracies have been normalized to objects of lower dimension, and
    // maximal polylines are constructed from undirected edges):
    //
    // bool ComputeUnion(S2ShapeIndex& a, S2ShapeIndex& b,
    //                   MutableS2ShapeIndex* index, S2Error error) {
    //   IndexedS2PolylineVectorLayer.Options polyline_options;
    //   polyline_options.set_edge_type(EdgeType.UNDIRECTED);
    //   polyline_options.set_polyline_type(Graph.PolylineType.WALK);
    //   polyline_options.set_duplicate_edges(DuplicateEdges.MERGE);
    //   LayerVector layers(3);
    //   layers[0] = new IndexedS2PointVectorLayer(index);
    //   layers[1] = new IndexedS2PolylineVectorLayer(
    //       index, polyline_options);
    //   layers[2] = new IndexedS2PolygonLayer(index);
    //   S2BooleanOperation op(S2BooleanOperation.OpType.UNION,
    //                         NormalizeClosedSet(move(layers)));
    //   return op.Build(a, b, error);
    // }
    public static List<Layer> NormalizeClosedSet(this List<Layer> output_layers, ClosedSetNormalizer.Options? options = null)
    {
        return NormalizeClosedSetImpl.Create(output_layers, options ?? new ClosedSetNormalizer.Options());
    }
}
