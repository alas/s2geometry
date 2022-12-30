// Given a set of possibly self-intersecting closed loops, this class computes
// a partitioning of the sphere into regions of constant winding number and
// returns the subset of those regions selected by a given winding rule.  This
// functionality can be used to implement N-way boolean polygon operations
// including union, intersection, and symmetric difference, as well as more
// exotic operations such as buffering and Minkowski sums.
//
// Recall that in the plane, the winding number of a closed curve around a
// point is an integer representing the number of times the curve travels
// counterclockwise around that point.  For example, points in the interior of
// a planar simple polygon with counterclockwise boundary have a winding
// number of +1, while points outside the polygon have a winding number of 0.
// If the boundary is clockwise then points in the interior have a winding
// number of -1.
//
// The interior of a complex self-intersecting closed boundary curve may be
// defined by choosing an appropriate winding rule.  For example, the interior
// could be defined as all regions whose winding number is positive, or all
// regions whose winding number is non-zero, or all regions whose winding
// number is odd.  Different winding rules are useful for implementing the
// various operations mentioned above (union, symmetric difference, etc).
//
// Unfortunately the concept of winding number on the sphere is not
// well-defined, due to the fact that a given closed curve does not partition
// the sphere into regions of constant winding number.  This is because the
// winding number changes when a point crosses either the given curve or that
// curve's reflection through the origin.
//
// Instead we engage in a slight abuse of terminology by modifying the concept
// of winding number on the sphere to be a relative one.  Given a set of
// closed curves on the sphere and the winding number of a reference point R,
// we define the winding number of every other point P by counting signed edge
// crossings.  When the edge RP crosses from the right side of a curve to its
// left the winding number increases by one, and when the edge RP crosses from
// the left side of a curve to its right the winding number decreases by one.
//
// This definition agrees with the classical one in the plane when R is taken
// to be the point at infinity with a winding number of zero.  It also agrees
// with the rule used by the S2 library to define polygon interiors, namely
// that the interior of a loop is the region to its left.  And most
// importantly, it satisfies the property that a closed curve partitions the
// sphere into regions of constant winding number.

namespace S2Geometry;

using EdgeType = S2Builder.EdgeType;
using Graph = S2Builder.Graph;
using GraphOptions = S2Builder.GraphOptions;

using DegenerateEdges = S2Builder.GraphOptions.DegenerateEdges;
using DuplicateEdges = S2Builder.GraphOptions.DuplicateEdges;
using SiblingPairs = S2Builder.GraphOptions.SiblingPairs;

public class S2WindingOperation
{
    private Options options_;
    private S2Builder builder_;
    private InputEdgeId ref_input_edge_id_;
    private int ref_winding_in_;
    private WindingRule rule_;

    // Default constructor; requires Init() to be called.
    public S2WindingOperation() { }

    // Convenience constructor that calls Init().
    public S2WindingOperation(S2Builder.Layer result_layer, Options? options = null)
    {
        Init(result_layer, options);
    }

    // Starts an operation that sends its output to the given S2Builder layer.
    // This method may be called more than once.
    public void Init(S2Builder.Layer result_layer, Options? options = null)
    {
        options_ = options ?? new();
        S2Builder.Options builder_options = new(options_.SnapFunction)
        {
            SplitCrossingEdges = true,
            MemoryTracker = options_.MemoryTracker,
        };
        builder_ = new(builder_options);
        builder_.StartLayer(new WindingLayer(this, result_layer));
    }

    // Adds a loop to the set of loops used to partition the sphere.  The given
    // loop may be degenerate or self-intersecting.
    public void AddLoop(S2PointLoopSpan loop)
    {
        builder_.AddLoop(loop);
    }

    // Given a reference point "ref_p" whose winding number is given to be
    // "ref_winding", snaps all the given input loops together using the given
    // snap_function() and then partitions the sphere into regions of constant
    // winding number.  As discussed above, the winding number increases by one
    // when any loop edge is crossed from right to left, and decreases by one
    // when any loop edge is crossed from left to right.  Each partition region
    // is classified as belonging to the interior of the result if and only if
    // its winding number matches the given rule (e.g. WindingRule::POSITIVE).
    //
    // The output consists of the boundary of this interior region plus possibly
    // certain degeneraices (as controlled by the include_degeneracies() option).
    // The boundary edges are sent to the S2Builder result layer specified in the
    // constructor, along with an appropriate IsFullPolygonPredicate that can be
    // used to distinguish whether the result is empty or full (even when
    // degeneracies are present).  Note that distingishing empty from full
    // results is a problem unique to spherical geometry.
    //
    // REQUIRES: error->ok() [an existing error will not be overwritten]
    public bool Build(S2Point ref_p, int ref_winding, WindingRule rule, out S2Error error)
    {
        // The reference point must be an S2Builder input vertex in order to
        // determine how its winding number is affected by snapping.  We record the
        // input edge id of this edge so that we can find it later.
        ref_input_edge_id_ = builder_.NumInputEdges();
        builder_.AddPoint(ref_p);
        ref_winding_in_ = ref_winding;
        rule_ = rule;
        return builder_.Build(out error);
    }

    public class Options
    {
        public Options()
        {
            SnapFunction = new S2BuilderUtil.IdentitySnapFunction(S1Angle.Zero);
        }

        // Convenience constructor that calls set_snap_function().
        public Options(S2BuilderUtil.SnapFunction snap_function)
        {
            SnapFunction = snap_function;
        }

        public Options(Options options)
        {
            SnapFunction = options.SnapFunction;
            IncludeDegeneracies = options.IncludeDegeneracies;
            MemoryTracker = options.MemoryTracker;
        }

        // Specifies the function used for snap rounding the output during the
        // call to Build().
        //
        // DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
        public S2BuilderUtil.SnapFunction SnapFunction { get; set; }

        // This option can be enabled to provide limited support for degeneracies
        // (i.e., sibling edge pairs and isolated vertices).  By default the output
        // does not include such edges because they do not bound any interior.  If
        // this option is true, then the output includes additional degeneracies as
        // follows:
        //
        //  - WindingRule::ODD outputs degeneracies whose multiplicity is odd.
        //
        //  - All other winding rules output degeneracies contained by regions
        //    whose winding number is zero.
        //
        // These rules are sufficient to implement the following useful operations:
        //
        //  - WindingRule::Odd can be used to compute the N-way symmetric
        //    difference of any combination of points, polylines, and polygons.
        //
        //  - WindingRule::POSITIVE can be used to compute the N-way union of any
        //    combination of points, polylines, and polygons.
        //
        // These statements only apply when closed boundaries are being used (see
        // S2BooleanOperation::{Polygon,Polyline}Model::CLOSED) and require the
        // client to convert points and polylines to degenerate loops and then back
        // again (e.g. using s2builderutil::ClosedSetNormalizer).  Note that the
        // handling of degeneracies is NOT sufficient for other polygon operations
        // (e.g. intersection) or other boundary models (e.g, semi-open).
        //
        // DEFAULT: false
        public bool IncludeDegeneracies { get; set; } = false;

        // Specifies that internal memory usage should be tracked using the given
        // S2MemoryTracker.  If a memory limit is specified and more more memory
        // than this is required then an error will be returned.  Example usage:
        //
        //   S2MemoryTracker tracker;
        //   tracker.set_limit(500 << 20);  // 500 MB
        //   S2WindingOperation::Options options;
        //   options.set_memory_tracker(&tracker);
        //   S2WindingOperation op{options};
        //   ...
        //   S2Error error;
        //   if (!op.Build(..., &error)) {
        //     if (error.code() == S2Error::RESOURCE_EXHAUSTED) {
        //       S2_LOG(ERROR) << error;  // Memory limit exceeded
        //     }
        //   }
        //
        // CAVEATS:
        //
        //  - Memory allocated by the output S2Builder layer is not tracked.
        //
        //  - While memory tracking is reasonably complete and accurate, it does
        //    not account for every last byte.  It is intended only for the
        //    purpose of preventing clients from running out of memory.
        //
        // DEFAULT: nullptr (memory tracking disabled)
        public S2MemoryTracker? MemoryTracker { get; set; } = null;
    }

    // Specifies the winding rule used to determine which regions belong to the
    // result.
    //
    // Note that additional winding rules may be created by adjusting the
    // winding number of the reference point.  For example, to select regions
    // whose winding number is at least 10, simply subtract 9 from ref_winding
    // and then use WindingRule::POSITIVE.  (This can be used to implement to
    // implement N-way polygon intersection).  Similarly, additional behaviors
    // can be obtained by reversing some of the input loops (e.g., this can be
    // used to compute one polygon minus the union of several other polygons).
    public enum WindingRule
    {
        POSITIVE,  // Winding number > 0
        NEGATIVE,  // Winding number < 0
        NON_ZERO,  // Winding number != 0
        ODD        // Winding number % 2 == 1
    }

    // The purpose of WindingOracle is to compute winding numbers with respect to
    // a set of input loops after snapping.  It is given the input edges (via
    // S2Builder), the output eges (an S2Builder::Graph), and the winding number
    // at a reference point with respect to the input edges ( "ref_input_edge_id"
    // and "ref_winding_in").  GetWindingNumber() can then be called to compute
    // the winding number at any point with respect to the snapped loops.  The
    // main algorithm uses this to compute the winding number at one arbitrary
    // vertex of each connected component of the S2Builder output graph.
    private class WindingOracle
    {
        public WindingOracle(InputEdgeId ref_input_edge_id, int ref_winding_in,
                    S2Builder builder, Graph g)
        {
            g_ = g;
            System.Diagnostics.Debug.Assert(g_.Options.EdgeType_ == EdgeType.DIRECTED);
            System.Diagnostics.Debug.Assert(g_.Options.DegenerateEdges_ == DegenerateEdges.KEEP);
            System.Diagnostics.Debug.Assert(g_.Options.DuplicateEdges_ == DuplicateEdges.KEEP);
            System.Diagnostics.Debug.Assert(g_.Options.SiblingPairs_ == SiblingPairs.KEEP);

            // Compute the winding number at the reference point after snapping (see
            // s2builderutil::GetSnappedWindingDelta).
            S2Point ref_in = builder.InputEdge(ref_input_edge_id).V0;
            VertexId ref_v = S2BuilderUtil.SnappedWindingDelta.FindFirstVertexId(ref_input_edge_id, g_);
            System.Diagnostics.Debug.Assert(ref_v >= 0);  // No errors are possible.
            CurrentRefPoint = g_.Vertex(ref_v);
            S2Error error;
            CurrentRefWinding = ref_winding_in + S2BuilderUtil.SnappedWindingDelta.GetSnappedWindingDelta(
                ref_in, ref_v, ie => true, builder, g_, out error);
            System.Diagnostics.Debug.Assert(error.IsOk());  // No errors are possible.

            // Winding numbers at other points are computed by counting signed edge
            // crossings.  If we need to do this many times, it is worthwhile to build an
            // index.  Note that although we initialize the index here, it is only built
            // the first time we use it (i.e., when brute_force_winding_tests_left_ < 0).
            //
            // As a small optimization, we set max_edges_per_cell() higher than its
            // default value.  This makes it faster to build the index but increases the
            // time per query.  This is a good tradeoff because the number of calls to
            // GetWindingNumber() is typically small relative to the number of graph
            // edges.  It also saves memory, which is important when the graph is very
            // large (e.g. because the input loops have many self-intersections).
            const int kMaxEdgesPerCell = 40;
            MutableS2ShapeIndex.Options options = new()
            {
                MaxEdgesPerCell = kMaxEdgesPerCell,
            };
            index_ = new(options)
            {
                MemoryTracker = builder.Options_.MemoryTracker,
            };
            index_.Add(new S2BuilderUtil.GraphShape(g_));
        }

        // Returns the winding number at the given point after snapping.
        // Returns the winding number at the given point "p".
        public int GetWindingNumber(S2Point p)
        {
            // Count signed edge crossings starting from the reference point, whose
            // winding number is known.  If we need to do this many times then we build
            // an S2ShapeIndex to speed up this process.
            S2EdgeCrosser crosser = new(CurrentRefPoint, p);
            int winding = CurrentRefWinding;
            if (--brute_force_winding_tests_left_ >= 0)
            {
                for (EdgeId e = 0; e < g_.NumEdges; ++e)
                {
                    winding += SignedCrossingDelta(crosser, e);
                }
            }
            else
            {
                S2CrossingEdgeQuery query = new(index_);
                foreach (var id in query.GetCandidates(CurrentRefPoint, p, index_.Shape(0)))
                {
                    winding += SignedCrossingDelta(crosser, id.EdgeId);
                }
            }
            // It turns out that GetWindingNumber() is called with a sequence of points
            // that are sorted in approximate S2CellId order.  This means that if we
            // update our reference point as we go along, the edges for which we need to
            // count crossings are much shorter on average.
            CurrentRefPoint = p;
            CurrentRefWinding = winding;
            return winding;
        }

        // Returns the current reference point.
        //
        // The current reference point.  Initially this is simply the snapped
        // position of the input reference point, but it is updated on each call to
        // GetWindingNumber (see below).
        public S2Point CurrentRefPoint { get; private set; }

        // Returns the winding number at the current reference point.
        //
        // The winding number at the current reference point.
        public int CurrentRefWinding { get; private set; }

        // Returns the change in winding number due to crossing the given graph edge.
        private int SignedCrossingDelta(S2EdgeCrosser crosser, EdgeId e) =>
            crosser.SignedEdgeOrVertexCrossing(
                g_.Vertex(g_.GetEdge(e).ShapeId),
                g_.Vertex(g_.GetEdge(e).EdgeId));

        private readonly Graph g_;

        // For each connected component of the S2Builder output graph, we determine
        // the winding number at one arbitrary vertex by counting edge crossings.
        // Initially we do this by brute force, but if there are too many connected
        // components then we build an S2ShapeIndex to speed up the process.
        //
        // Building an index takes about as long as 25 brute force queries.  The
        // classic competitive algorithm technique would be to do 25 brute force
        // queries and then build the index.  However in practice, inputs with one
        // connected component are common while inputs with 2-25 connected
        // components are rare.  Therefore we build the index as soon as we realize
        // that there is more than one component.
        //
        // Another idea is to count the connected components in advance, however it
        // turns out that this takes about 25% as long as building the index does.
        private int brute_force_winding_tests_left_ = 1;
        private readonly MutableS2ShapeIndex index_;  // Built only if needed.
    }

    // The actual winding number operation is implemented as an S2Builder layer.
    private class WindingLayer : S2Builder.Layer
    {
        public WindingLayer(S2WindingOperation op, S2Builder.Layer result_layer)
        {
            op_ = op;
            result_layer_ = result_layer;
            tracker_ = new(op.options_.MemoryTracker);
        }

        // Layer interface:
        public override GraphOptions GraphOptions_()
        {
            // The algorithm below has several steps with different graph requirements,
            // however the first step is to determine how snapping has affected the
            // winding number of the reference point.  This requires keeping all
            // degenerate edges and sibling pairs.  We also keep all duplicate edges
            // since this makes it easier to remove the reference edge.
            return new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                                DuplicateEdges.KEEP, SiblingPairs.KEEP);
        }

        public override void Build(Graph g, out S2Error error)
        {
            error = S2Error.OK;
            //if (!error.IsOk()) return;  // Abort if an error already exists.

            // The WindingOracle computes the winding number of any point on the sphere.
            // It requires knowing the winding number at a reference point which must be
            // an input vertex to S2Builder.  This is achieved by adding a degenerate
            // edge from the reference point to itself (ref_input_edge_id_).  Once we
            // have computed the change in winding number, we create a new graph with
            // this edge removed (since it should not be emitted to the result layer).
            WindingOracle oracle = new(op_.ref_input_edge_id_, op_.ref_winding_in_,
                                 op_.builder_, g);
            //System.Diagnostics.Debug.Assert(error.IsOk());  // No errors are possible.

            // Now that we have computed the change in winding number, we create a new
            // graph with the reference edge removed.  Note that S2MemoryTracker errors
            // are automatically copied to "error" by S2Builder.
            List<Edge> new_edges = new();
            List<InputEdgeIdSetId> new_input_edge_ids = new();
            if (!tracker_.AddSpace(new_edges, g.NumEdges - 1)) return;
            if (!tracker_.AddSpace(new_input_edge_ids, g.NumEdges - 1)) return;
            var new_input_edge_id_set_lexicon = g.InputEdgeIdSetLexicon;

            // Copy all of the edges except the reference edge.
            for (int e = 0; e < g.NumEdges; ++e)
            {
                if (g.InputEdgeIds(e).First() == op_.ref_input_edge_id_) continue;
                new_edges.Add(g.GetEdge(e));
                new_input_edge_ids.Add(g.InputEdgeIdSetId(e));
            }

            // Our goal is to assemble the given edges into loops that parition the
            // sphere.  In order to do this we merge duplicate edges and create sibling
            // edges so that every region can have its own directed boundary loop.
            //
            // Isolated degenerate edges and sibling pairs are preserved in order to
            // provide limited support for working with geometry of dimensions 0 and 1
            // (i.e., points and polylines).  Clients can simply convert these objects to
            // degenerate loops and then convert these degenerate loops back to
            // points/polylines in the output using s2builderutil::ClosedSetNormalizer.
            GraphOptions new_options = new(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                                     DuplicateEdges.MERGE, SiblingPairs.CREATE);
            Graph? new_graph = g.MakeSubgraph(
                new_options, new_edges, new_input_edge_ids,
                new_input_edge_id_set_lexicon, null, out error, tracker_);
            if (!error.IsOk()) return;

            // Now visit each connected component of the graph and assemble its edges
            // into loops.  For each loop we determine the winding number of the region
            // to its left, and if the winding number matches the given rule then we
            // add its edges to result_edges_.
            if (!ComputeBoundary(new_graph!, oracle, out error)) return;

            // Now we construct yet another S2Builder::Graph by starting with the edges
            // that bound the desired output region and then processing them according to
            // the client's requested GraphOptions.  (Note that ProcessEdges can change
            // these options in certain cases; see S2Builder::GraphOptions for details.)
            //
            // The IsFullPolygonPredicate below allows clients to distinguish full from
            // empty results (including cases where there are degeneracies).  Note that
            // we can use the winding number at any point on the sphere for this
            // purpose.
            S2Builder.IsFullPolygonPredicate is_full_polygon_predicate = (Graph g, out S2Error error) => {
                error = S2Error.OK;
                return MatchesRule(oracle.CurrentRefWinding);
            };
            var result_id_set_lexicon = new_graph.InputEdgeIdSetLexicon;
            var result_graph = new_graph.MakeSubgraph(
                result_layer_.GraphOptions_(), result_edges_, result_input_edge_ids_,
                result_id_set_lexicon, is_full_polygon_predicate, out error, tracker_);
            if (tracker_.Clear(new_edges) && tracker_.Clear(new_input_edge_ids))
            {
                result_layer_.Build(result_graph!, out error);
            }
        }

        private bool ComputeBoundary(Graph g, WindingOracle oracle, out S2Error error)
        {
            error = S2Error.OK;
            // We assemble the edges into loops using an algorithm similar to
            // S2Builder::Graph::GetDirectedComponents(), except that we also keep track
            // of winding numbers and output the relevant edges as we go along.
            //
            // The following accounts for sibling_map, left_turn_map, and edge_winding
            // (which have g.num_edges() elements each).
            var kTempUsage = 3 * sizeof(EdgeId) * g.NumEdges;
            if (!tracker_.Tally(kTempUsage)) return false;

            List<EdgeId> sibling_map = g.GetSiblingMap();
            List<EdgeId> left_turn_map;
            g.GetLeftTurnMap(sibling_map, out left_turn_map, out error);
            System.Diagnostics.Debug.Assert(error.IsOk());

            // A map from EdgeId to the winding number of the region it bounds.
            List<int> edge_winding = new(g.NumEdges);
            List<EdgeId> frontier = new();  // Unexplored sibling edges.
            for (EdgeId e_min = 0; e_min < g.NumEdges; ++e_min)
            {
                if (left_turn_map[e_min] < 0) continue;  // Already visited.

                // We have found a new connected component of the graph.  Each component
                // consists of a set of loops that partition the sphere.  We start by
                // choosing an arbitrary vertex "v0" and computing the winding number at
                // this vertex.  Recall that point containment is defined such that when a
                // set of loops partition the sphere, every point is contained by exactly
                // one loop.  Therefore the winding number at "v0" is the same as the
                // winding number of the adjacent loop that contains it.  We choose "e0" to
                // be an arbitrary edge of this loop (it is the incoming edge to "v0").
                VertexId v0 = g.GetEdge(e_min).EdgeId;
                EdgeId e0 = GetContainingLoopEdge(v0, e_min, g, left_turn_map, sibling_map);
                edge_winding[e0] = oracle.GetWindingNumber(g.Vertex(v0));

                // Now visit all loop edges in this connected component of the graph.
                // "frontier" is a stack of unexplored siblings of the edges visited far.
                if (!tracker_.AddSpace(frontier, 1)) return false;
                frontier.Add(e0);
                while (frontier.Any())
                {
                    EdgeId e = frontier.Last();
                    frontier.RemoveAt(frontier.Count - 1);
                    if (left_turn_map[e] < 0) continue;  // Already visited.

                    // Visit all edges of the loop starting at "e".
                    int winding = edge_winding[e];
                    for (EdgeId next; left_turn_map[e] >= 0; e = next)
                    {
                        // Count signed edge crossings to determine the winding number of
                        // the sibling region.  Input edges that snapped to "e" decrease the
                        // winding number by one (since we cross them from left to right),
                        // while input edges that snapped to its sibling edge increase the
                        // winding number by one (since we cross them from right to left).
                        EdgeId sibling = sibling_map[e];
                        int winding_minus = g.InputEdgeIds(e).Count;
                        int winding_plus = g.InputEdgeIds(sibling).Count;
                        int sibling_winding = winding - winding_minus + winding_plus;

                        // Output all edges that bound the region selected by the winding
                        // rule, plus certain degenerate edges.
                        if ((MatchesRule(winding) && !MatchesRule(sibling_winding)) ||
                            MatchesDegeneracy(winding, winding_minus, winding_plus))
                        {
                            OutputEdge(g, e);
                        }
                        next = left_turn_map[e];
                        left_turn_map[e] = -1;
                        // If the sibling hasn't been visited yet, add it to the frontier.
                        if (left_turn_map[sibling] >= 0)
                        {
                            edge_winding[sibling] = sibling_winding;
                            if (!tracker_.AddSpace(frontier, 1)) return false;
                            frontier.Add(sibling);
                        }
                    }
                }
            }
            tracker_.Untally(frontier);
            return tracker_.Tally(-kTempUsage);
        }

        // Given an incoming edge "start" to a vertex "v", returns an edge of the loop
        // that contains "v" with respect to the usual semi-open boundary rules.
        private static EdgeId GetContainingLoopEdge(VertexId v, EdgeId start, Graph g,
                                   List<EdgeId> left_turn_map, List<EdgeId> sibling_map)
        {
            // If the given edge is degenerate, this is an isolated vertex.
            System.Diagnostics.Debug.Assert(g.GetEdge(start).EdgeId == v);
            if (g.GetEdge(start).ShapeId == v) return start;

            // Otherwise search for the loop that contains "v".
            EdgeId e0 = start;
            for (; ; )
            {
                EdgeId e1 = left_turn_map[e0];
                System.Diagnostics.Debug.Assert(g.GetEdge(e0).EdgeId == v);
                System.Diagnostics.Debug.Assert(g.GetEdge(e1).ShapeId == v);

                // The first test below handles the case where there are only two edges
                // incident to this vertex (i.e., the vertex angle is 360 degrees).
                if (g.GetEdge(e0).ShapeId == g.GetEdge(e1).EdgeId ||
                    S2.AngleContainsVertex(g.Vertex(g.GetEdge(e0).ShapeId), g.Vertex(v),
                                            g.Vertex(g.GetEdge(e1).EdgeId)))
                {
                    return e0;
                }
                e0 = sibling_map[e1];
                System.Diagnostics.Debug.Assert(e0 != start);
            }
            throw new Exception("Shoul not be possible to get here");
        }

        private bool MatchesRule(int winding)
        {
            switch (op_.rule_)
            {
                case WindingRule.POSITIVE: return winding > 0;
                case WindingRule.NEGATIVE: return winding < 0;
                case WindingRule.NON_ZERO: return winding != 0;
                case WindingRule.ODD: return (winding & 1) != 0;
            }
            throw new Exception("Shoul not be possible to get here");
        }

        private bool MatchesDegeneracy(int winding, int winding_minus, int winding_plus)
        {
            if (!op_.options_.IncludeDegeneracies) return false;

            // A degeneracy is either a self-loop or a sibling pair where equal numbers
            // of input edges snapped to both edges of the pair.  The test below covers
            // both cases (because self-loops are their own sibling).
            if (winding_minus != winding_plus) return false;

            if (op_.rule_ == WindingRule.ODD)
            {
                // Any degeneracy whose multiplicity is odd should be part of the result
                // independent of the winding number of the region that contains it.
                // This rule allows computing symmetric differences of any combination of
                // points, polylines, and polygons (where the first two are represented as
                // degenerate loops).
                return (winding_plus & 1) != 0;
            }
            else
            {
                // For all other winding rules we output degeneracies only if they are
                // contained by a region of winding number zero.  Even though the interface
                // to this class does not provide enough information to allow consistent
                // handling of degeneracies in general, this rule is sufficient for several
                // important cases.  Specifically it allows computing N-way unions of any
                // mixture of points, polylines, and polygons by converting the points and
                // polylines to degenerate loops.  In this case all input loops are
                // degenerate or CCW, and the boundary of the result can be computed using
                // WindingRule::POSITIVE.  Since there are no clockwise loops, all
                // degeneracies contained by a region of winding number zero represent
                // degenerate shells and should be emitted.  (They can be converted back to
                // points/polylines using s2builderutil::ClosedSetNormalizer.)
                //
                // Similarly, this heuristic is sufficient to compute unions of points,
                // polylines, and polygons where all boundaries are clockwise (by using
                // WindingRule::NEGATIVE) or where all boundaries are of an unknown but
                // consistent oreientation (by using WindingRule::NON_ZERO).
                return winding == 0;
            }
        }

        // Adds an edge to the set of output edges.
        private void OutputEdge(Graph g, EdgeId e)
        {
            if (!tracker_.AddSpace(result_edges_, 1)) return;
            if (!tracker_.AddSpace(result_input_edge_ids_, 1)) return;
            result_edges_.Add(g.GetEdge(e));
            result_input_edge_ids_.Add(g.InputEdgeIdSetId(e));
        }

        // Constructor parameters.
        private readonly S2WindingOperation op_;
        private readonly S2Builder.Layer result_layer_;

        // The graph data that will be sent to result_layer_.
        private readonly List<Edge> result_edges_ = new();
        private readonly List<InputEdgeIdSetId> result_input_edge_ids_ = new();

        private readonly S2MemoryTracker.Client tracker_;
    }
}
