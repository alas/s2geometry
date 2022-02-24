// The following algorithm would not be necessary with planar geometry, since
// then winding numbers could be computed by starting at a point at infinity
// (whose winding number is zero) and counting signed edge crossings.  However
// points at infinity do not exist on the sphere.
//
// Instead we compute the change in winding number of a reference vertex R
// using only the set of edges incident to the snapped reference vertex R'.
// Essentially this involves looking at the set of input edges that snapped to
// R' and assembling them into edge chains.  These edge chains can be divided
// into two categories:
//
// (1) Edge chains that are entirely contained by the Voronoi region of R'.
//     This means that the input edges form a closed loop where every vertex
//     snaps to R'.  We can compute the change in winding number due to this
//     loop by simply choosing a point Z outside the Voronoi region of R' and
//     computing the winding numbers of R and R' relative to Z.
//
// (2) Edge chains that enter the Voronoi region of R' and later leave it.  In
//     this case the input chain has the form C = (A0, A1, ..., B0, B1) where
//     A0 and B1 are outside the Voronoi region of R' and all other vertices
//     snap to R'.  In the vicinity of R' this input chain snaps to a chain C'
//     of the form (A0', R', B1') where A0' is the second-last vertex in the
//     snapped edge chain for A0A1 and B1' is the second vertex in the snapped
//     edge chain for B0B1.  In principle we handle this similarly to the case
//     above (by finding a point Z whose change in winding number is known,
//     and then counting signed edge crossings along ZR with respect to C and
//     along ZR' with respect to C').  However the details are more
//     complicated and are described in GetSnappedWindingDelta().
//
// The total change in winding number is simply the sum of the changes in
// winding number due to each of these edge chains.

// A map that allows finding all the input edges that start at a given point.

namespace S2Geometry.S2BuilderUtil;

using static S2Builder;
using InputVertexEdgeMap = List<SnappedWindingDelta.InputVertexEdge>;

public class SnappedWindingDelta
{
    // The winding number returned when a usage error is detected.
    private const int kErrorResult = int.MaxValue;

    // A function that returns true if the given S2Builder input edge should be
    // ignored in the winding number calculation.  This means either that the edge
    // is not a loop edge (e.g., a non-closed polyline) or that this loop should
    // not affect the winding number.  This is useful for two purposes:
    //
    //  - To process graphs that contain polylines and points in addition to loops.
    //
    //  - To process graphs where the winding number is computed with respect to
    //    only a subset of the input loops.
    //
    // It can be default-constructed to indicate that no edges should be ignored.
    public delegate bool InputEdgeFilter(InputEdgeId inputEdgeId);

    // Given an S2Builder::Graph of output edges after snap rounding and a
    // reference vertex R, computes the change in winding number of R due to
    // snapping.  (See S2WindingOperation for an introduction to winding numbers on
    // the sphere.)  The return value can be added to the original winding number
    // of R to obtain the winding number of the corresponding snapped vertex R'.
    //
    // The algorithm requires that the S2Builder input edges consist entirely of
    // (possibly self-intersecting) closed loops.  If you need to process inputs
    // that include other types of geometry (e.g., non-closed polylines), you will
    // need to either (1) put them into a different S2Builder layer, (2) close the
    // polylines into loops (e.g. using GraphOptions::SiblingEdges::CREATE), or (3)
    // provide a suitable InputEdgeFilter (see above) so that the non-loop edges
    // can be ignored.
    //
    // The algorithm is designed to be robust for any input edge configuration and
    // snapping result.  However note that it cannot be used in conjunction with
    // edge chain simplification (S2Builder::Options::simplify_edge_chains).  It
    // also requires that S2Builder::GraphOptions be configured to keep all snapped
    // edges, even degenerate ones (see requirements below).
    //
    // "ref_in" is the reference vertex location before snapping.  It *must* be an
    // input vertex to S2Builder, however this is not checked.
    //
    // "ref_v" is the Graph::VertexId of the reference vertex after snapping.
    // (This can be found using the FindFirstVertexId() function below if desired.)
    //
    // "input_edge_filter" can optionally be used to ignore input edges that
    // should not affect the winding number calculation (such as polyline edges).
    // The value can be default-constructed (InputEdgeFilter{}) to use all edges.
    //
    // "builder" is the S2Builder that produced the given edge graph.  It is used
    // to map InputEdgeIds back to the original edge definitions, and also to
    // verify that no incompatible S2Builder::Options were used (see below).
    //
    // "g" is the S2Builder output graph of snapped edges.
    //
    // The only possible errors are usage errors, in which case "error" is set to
    // an appropriate error message and a very large value is returned.
    //
    // Example usage:
    //
    // This function is generally called from an S2Builder::Layer implementation.
    // We assume here that the reference vertex is the first vertex of the input
    // edge identified by "ref_input_edge_id_", and that its desired winding number
    // with respect to the input loops is "ref_winding_".
    //
    // using Graph = S2Builder::Graph;
    // class SomeLayer : public S2Builder::Layer {
    //  private:
    //   int ref_input_edge_id_;
    //   int ref_winding_;
    //   const S2Builder& builder_;
    //
    //  public:
    //   ...
    //   void Build(const Graph& g, S2Error* error) {
    //     // Find the positions of the reference vertex before and after snapping.
    //     S2Point ref_in = builder_.input_edge(ref_input_edge_id_).v0;
    //     Graph::VertexId ref_v =
    //         s2builderutil::FindFirstVertexId(ref_input_edge_id_, g);
    //     S2Point ref_out = g.vertex(ref_v);
    //
    //     // Compute the change in winding number due to snapping.
    //     S2Error error;
    //     ref_winding_ += s2builderutil::GetSnappedWindingDelta(
    //         ref_in, ref_v, InputEdgeFilter{}, builder_, g, error);
    //     S2_CHECK(error->ok());  // All errors are usage errors.
    //
    //     // Winding numbers of others points can now be found by counting signed
    //     // edge crossings (S2EdgeCrosser::SignedEdgeOrVertexCrossing) between
    //     // "ref_out" and the desired point.  Note that if DuplicateEdges::MERGE
    //     // or SiblingPairs::CREATE was used, each crossing has a multiplicity
    //     // equal to the number of non-filtered input edges that snapped to that
    //     // output edge.
    //   }
    // }
    //
    // REQUIRES: The input edges after filtering consist entirely of closed loops.
    //           (If DuplicateEdges::MERGE or SiblingPairs::CREATE was used,
    //           each graph edge has a multiplicity equal to the number of
    //           non-filtered input edges that snapped to it.)
    //
    // REQUIRES: g.options().edge_type() == DIRECTED
    // REQUIRES: g.options().degenerate_edges() == KEEP
    // REQUIRES: g.options().sibling_pairs() == {KEEP, REQUIRE, CREATE}
    // REQUIRES: builder.options().simplify_edge_chains() == false
    //
    // CAVEAT: The running time is proportional to the S2Builder::Graph size.  If
    //         you need to call this function many times on the same graph then
    //         use the alternate version below.  (Most clients only need to call
    //         GetSnappedWindingDelta() once per graph because the winding numbers
    //         of other points can be computed by counting signed edge crossings.)
    public static int GetSnappedWindingDelta(
        S2Point ref_in, VertexId ref_v,
        InputEdgeFilter? input_edge_filter, S2Builder builder,
        Graph g, out S2Error error)
    {
        return GetSnappedWindingDelta(
            ref_in, ref_v, GetIncidentEdgesBruteForce(ref_v, g),
            input_edge_filter, builder, g, out error);
    }

    // This version can be used when GetSnappedWindingDelta() needs to be called
    // many times on the same graph.  It is faster than the function above, but
    // less convenient to use because it requires the client to provide the set of
    // graph edges incident to the snapped reference vertex.  It runs in time
    // proportional to the size of this set.
    //
    // "incident_edges" is the set of incoming and outgoing graph edges incident
    // to ref_v.  (These edges can be found efficiently using Graph::VertexOutMap
    // and Graph::VertexInMap.)
    //
    // See the function above for the remaining parameters and requirements.
    public static int GetSnappedWindingDelta(
        S2Point ref_in, VertexId ref_v, List<EdgeId> incident_edges,
        InputEdgeFilter? input_edge_filter, S2Builder builder,
        Graph g, out S2Error error)
    {
        error = S2Error.OK;

        System.Diagnostics.Debug.Assert(!builder.Options_.SimplifyEdgeChains);
        System.Diagnostics.Debug.Assert(g.Options.EdgeType_ == S2Builder.EdgeType.DIRECTED);
        System.Diagnostics.Debug.Assert(g.Options.DegenerateEdges_ == GraphOptions.DegenerateEdges.KEEP);
        System.Diagnostics.Debug.Assert(
            new[]
            {
                GraphOptions.SiblingPairs.KEEP,
                GraphOptions.SiblingPairs.REQUIRE,
                GraphOptions.SiblingPairs.CREATE
            }.Contains(g.Options.SiblingPairs_));

        // First we group all the incident edges by input edge id, to handle the
        // problem that input edges can map to either one or two snapped edges.
        Dictionary<InputEdgeId, EdgeSnap> input_id_edge_map = new();
        foreach (EdgeId e in incident_edges)
        {
            var edge = g.GetEdge(e);
            foreach (var input_id in g.InputEdgeIds(e))
            {
                if (input_edge_filter != null && input_edge_filter(input_id)) continue;
                var snap = input_id_edge_map[input_id];
                snap.input = builder.input_edge(input_id);
                if (edge.ShapeId != ref_v) snap.v_in = edge.ShapeId;
                if (edge.EdgeId != ref_v) snap.v_out = edge.EdgeId;
            }
        }
        // Now we regroup the edges according to the reference vertex of the
        // corresponding input edge.  This makes it easier to assemble these edges
        // into (portions of) input edge loops.
        InputVertexEdgeMap input_vertex_edge_map = new();
        foreach (var entry in input_id_edge_map)
        {
            var snap = entry.Value;
            input_vertex_edge_map.AddSorted(new InputVertexEdge(snap.input.V0, snap));
        }

        // The position of the reference vertex after snapping.  In comments we will
        // refer to the reference vertex before and after snapping as R and R'.
        S2Point ref_out = g.Vertex(ref_v);

        // Now we repeatedly assemble edges into an edge chain and determine the
        // change in winding number due to snapping of that edge chain.  These
        // values are summed to compute the final winding number delta.
        //
        // An edge chain is either a closed loop of input vertices where every
        // vertex snapped to the reference vertex R', or a partial loop such that
        // all interior vertices snap to R' while the first and last vertex do not.
        // Note that the latter includes the case where the first and last input
        // vertices in the chain are identical but do not snap to R'.
        //
        // Essentially we compute the winding number of the unsnapped reference
        // vertex R with respect to the input chain and the winding number of the
        // snapped reference vertex R' with respect to the output chain, and then
        // subtract them.  In the case of open chains, we compute winding numbers as
        // if the chains had been closed in a way that preserves topology while
        // snapping (i.e., in a way that does not the cause chain to pass through
        // the reference vertex as it continuously deforms from the input to the
        // output).
        //
        // Any changes to this code should be validated by running the RandomLoops
        // unit test with at least 10 million iterations.
        int winding_delta = 0;
        while (input_vertex_edge_map.Any())
        {
            List<S2Point> chain_in = new(), chain_out = new();
            if (!BuildChain(ref_v, g, input_vertex_edge_map,
                            chain_in, chain_out, out error))
            {
                return kErrorResult;
            }
            if (chain_out.Count == 1)
            {
                // We have a closed chain C of input vertices such that every vertex
                // snaps to the reference vertex R'.  Therefore we can easily find a
                // point Z whose winding number is not affected by the snapping of C; it
                // just needs to be outside the Voronoi region of R'.  Since the snap
                // radius is at most 70 degrees, we can use a point 90 degrees away such
                // as S2::Ortho(R').
                //
                // We then compute the winding numbers of R and R' relative to Z.  We
                // compute the winding number of R by counting signed crossings of the
                // edge ZR, while the winding number of R' relative to Z is always zero
                // because the snapped chain collapses to a single point.
                System.Diagnostics.Debug.Assert(chain_out[0] == ref_out);         // Snaps to R'.
                System.Diagnostics.Debug.Assert(chain_in[0] == chain_in.Last());  // Chain is a loop.
                S2Point z = S2.Ortho(ref_out);
                winding_delta += 0 - GetEdgeWindingDelta(z, ref_in, chain_in);
            }
            else
            {
                // We have an input chain C = (A0, A1, ..., B0, B1) that snaps to a
                // chain C' = (A0', R', B1'), where A0 and B1 are outside the Voronoi
                // region of R' and all other input vertices snap to R'.  This includes
                // the case where A0 == B1 and also the case where the input chain
                // consists of only two vertices.  Note that technically the chain C
                // snaps to a supersequence of C', since A0A1 snaps to a chain whose
                // last two vertices are (A0', R') and B0B1 snaps to a chain whose first
                // two vertices are (R', B1').  This implies that A0 does not
                // necessarily snap to A0', and similarly for B1 and B1'.
                //
                // Note that A0 and B1 can be arbitrarily far away from R'.  This makes
                // it difficult (on the sphere) to construct a point Z whose winding
                // number is guaranteed not to be affected by snapping the edges A0A1
                // and B0B1.  Instead we construct two points Za and Zb such that Za is
                // guaranteed not be affected by the snapping of A0A1, Zb is guaranteed
                // not to be affected by the snapping of B0B1, and both points are
                // guaranteed not to be affected by the snapping of any other input
                // edges.  We can then compute the change in winding number of Zb by
                // considering only the single edge A0A1 that snaps to A0'R'.
                // Furthermore we can compute the latter by using Za as the reference
                // point, since its winding number is guaranteed not be affected by this
                // particular edge.
                //
                // Given the point Zb, whose change in winding number is now known, we
                // can compute the change in winding number of the reference vertex R.
                // We essentially want to count the signed crossings of ZbR' with respect
                // to C' and subtract the signed crossings of ZbR with respect to C,
                // which we will write as s(ZbR', C') - s(ZbR, C).
                //
                // However to do this we need to close both chains in a way that is
                // topologically consistent and does not affect the winding number of
                // Zb.  This can be achieved by counting the signed crossings of ZbR' by
                // following the two-edge path (Zb, R, R').  In other words, we compute
                // this as s(ZbR, C') + s(RR', C') - s(ZbR, C).  We can then compute
                // s(ZbR, C') - s(ZbR, C) by simply concatenating the vertices of C in
                // reverse order to C' to form a single closed loop.  The remaining term
                // s(RR', C') can be implemented as signed edge crossing tests, or more
                // directly by testing whether R is contained by the wedge C'.
                System.Diagnostics.Debug.Assert(chain_out.Count == 3);
                System.Diagnostics.Debug.Assert(chain_out[1] == ref_out);

                // Compute two points Za and Zb such that Za is not affected by the
                // snapping of any edge except possibly B0B1, and Zb is not affected by
                // the snapping of any edge except possibly A0A1.  Za and Zb are simply
                // the normals to the edges A0A1 and B0B1 respectively, with their sign
                // chosen to point away from the Voronoi site R'.  This ensures at least
                // 20 degrees of separation from all edges except the ones mentioned.
                S2Point za = S2.RobustCrossProd(chain_in[0], chain_in[1]).Normalize();
                S2Point zb = S2.RobustCrossProd(chain_in[^2], chain_in.Last())
                            .Normalize();
                if (za.DotProd(ref_out) > 0) za = -za;
                if (zb.DotProd(ref_out) > 0) zb = -zb;

                // We now want to determine the change in winding number of Zb due to
                // A0A1 snapping to A0'R'.  Conceptually we do this by closing these
                // two single-edge chains into loops L and L' and then computing
                // s(ZaZb, L') - s(ZaZb, L).  Recall that Za was constructed so as not
                // to be affected by the snapping of edge A0A1, however this is only
                // true provided that L can snap to L' without passing through Za.
                //
                // To achieve this we let L be the degenerate loop (A0, A1, A0), and L'
                // be the loop (A0', R', A1, A0, A0').  The only problem is that we need
                // to ensure that the edge A0A0' stays within 90 degrees of A0A1, since
                // otherwise when the latter edge snaps to the former it might pass
                // through Za.  (This problem arises because we only consider the last
                // two vertices (A0', R') that A0A1 snaps to.  If we used the full chain
                // of snapped vertices for A0A1 then L' would always stay within the
                // snap radius of this edge.)
                //
                // The simplest way to fix this problem is to insert a connecting vertex
                // Ac between A0 and A0'.  THis vertex acts as a proxy for the missing
                // snapped vertices, yielding the loop L' = (A0', R', A1, A0, Ac, A0').
                // The vertex Ac is located on the edge A0A1 and is at most 90 degrees
                // away from A0'.  This ensures that the chain (A0, Ac, A0') always
                // stays within the snap radius of the input edge A0A1.
                //
                // Similarly we insert a connecting vertex Bc between B0 and B1 to
                // ensure that the edge B1'B1 passes on the correct side of Zb.
                S2Point a0_connector = GetConnector(chain_in[1], chain_in[0],
                                                    chain_out[0]);
                S2Point b1_connector = GetConnector(chain_in[^2], chain_in.Last(),
                                                    chain_out[2]);

                // Compute the change in winding number for reference vertex Zb.  Note
                // that we must duplicate the first/last vertex of the loop since the
                // argument to GetEdgeWindingDelta() is a polyline.
                List<S2Point> chain_z = new()
                {
                    chain_out[0],
                    chain_out[1],
                    chain_in[1],
                    chain_in[0],
                    a0_connector,
                    chain_out[0]
                };
                winding_delta += GetEdgeWindingDelta(za, zb, chain_z);

                // Compute the change in winding number of ZbR due to snapping C to C'.
                // As before, conceptually we do this by closing these chains into loops
                // L and L' such that L snaps to L' without passing through Zb.  Again
                // this can be done by concatenating the vertices of C' with the
                // reversed vertices of C, along with the two extra connecting vertices
                // Ac and Bc to ensure that L and L' pass on the same side of Za and Zb.
                // This yields the loop (A0', R', B1', Bc, B1, B0, ..., A1, A0, Ac, A0').
                List<S2Point> chain_diff = chain_out;
                chain_diff.Add(b1_connector);
                chain_diff.AddRange(chain_in.AsEnumerable().Reverse());
                chain_diff.Add(a0_connector);
                chain_diff.Add(chain_out[0]);  // Close the loop.
                winding_delta += GetEdgeWindingDelta(zb, ref_in, chain_diff);

                // Compute the change in winding number of RR' with respect to C' only.
                // (This could also be computed using two calls to s2pred::OrderedCCW.)
                System.Diagnostics.Debug.Assert(chain_out[1] == ref_out);
                winding_delta += GetEdgeWindingDelta(ref_in, ref_out, chain_out);
            }
        }
        return winding_delta;
    }

    // Returns the set of incoming and outgoing edges incident to the given
    // vertex.  This method takes time linear in the size of the graph "g";
    // if you need to call this function many times then it is more efficient to
    // use Graph::VertexOutMap and Graph::VertexInMap instead.
    private static List<EdgeId> GetIncidentEdgesBruteForce(VertexId v, Graph g)
    {
        List<EdgeId> result = new();
        for (EdgeId e = 0; e < g.NumEdges; ++e)
        {
            if (g.GetEdge(e).ShapeId == v || g.GetEdge(e).EdgeId == v)
            {
                result.Add(e);
            }
        }
        return result;
    }

    private static bool BuildChain(
        VertexId ref_v, Graph g, InputVertexEdgeMap input_vertex_edge_map,
        List<S2Point> chain_in, List<S2Point> chain_out, out S2Error error)
    {
        error = S2Error.OK;
        System.Diagnostics.Debug.Assert(!chain_in.Any());
        System.Diagnostics.Debug.Assert(!chain_out.Any());

        // First look for an incoming edge to the reference vertex.  (This will be
        // the start of a chain that eventually leads to an outgoing edge.)
        {
            S2Point? snapKey = null;
            EdgeSnap? snapCopy = null;
            foreach (var it in input_vertex_edge_map)
            {
                var snap = it.EdgeSnap;
                if (snap.v_in >= 0)
                {
                    chain_out.Add(g.Vertex(snap.v_in));
                    snapKey = it.Point;
                    snapCopy = snap;
                    break;
                }
            }
            if (snapCopy == null)
            {
                snapKey = input_vertex_edge_map.First().Point;
                // Pick an arbitrary edge to start a closed loop.
                snapCopy = input_vertex_edge_map.First().EdgeSnap;
            }
            input_vertex_edge_map.Remove(new InputVertexEdge(snapKey!.Value, new()));

            chain_in.Add(snapCopy.Value.input.V0);
            chain_in.Add(snapCopy.Value.input.V1);
            chain_out.Add(g.Vertex(ref_v));
            if (snapCopy.Value.v_out >= 0)
            {
                // This input edge enters and immediately exits the Voronoi region.
                chain_out.Add(g.Vertex(snapCopy.Value.v_out));
                return true;
            }
        }

        // Now repeatedly add edges until the chain or loop is finished.
        while (chain_in.Last() != chain_in.First())
        {
            var value = new InputVertexEdge(chain_in.Last(), new());
            var firstIndex = input_vertex_edge_map.GetLowerBound(value);
            var first = input_vertex_edge_map[firstIndex];
            var secondIndex = input_vertex_edge_map.GetUpperBound(value);
            var second = input_vertex_edge_map[secondIndex];

            if (first == second)
            {
                error = new(S2ErrorCode.INVALID_ARGUMENT,
                    "Input edges (after filtering) do not form loops");
                return false;
            }
            EdgeSnap snap = first.EdgeSnap;
            input_vertex_edge_map.Remove(first);
            chain_in.Add(snap.input.V1);
            if (snap.v_out >= 0)
            {
                // The chain has exited the Voronoi region.
                chain_out.Add(g.Vertex(snap.v_out));
                break;
            }
        }
        return true;
    }

    // Returns the change in winding number along the edge AB with respect to the
    // given edge chain.  This is simply the sum of the signed edge crossings.
    private static int GetEdgeWindingDelta(S2Point a, S2Point b,
                            List<S2Point> chain)
    {
        System.Diagnostics.Debug.Assert(chain.Count > 0);

        int delta = 0;
        S2EdgeCrosser crosser = new(a, b, chain[0]);
        for (int i = 1; i < chain.Count; ++i)
        {
            delta += crosser.SignedEdgeOrVertexCrossing(chain[i]);
        }
        return delta;
    }

    // Given an input edge (B0, B1) that snaps to an edge chain (B0', B1', ...),
    // returns a connecting vertex "Bc" that can be used as a substitute for the
    // remaining snapped vertices "..." when computing winding numbers.  This
    // requires that (1) the edge (B1', Bc) does not intersect the Voronoi region
    // of B0', and (2) the edge chain (B0', B1', Bc, B1) always stays within the
    // snap radius of the input edge (B0, B1).
    private static S2Point GetConnector(S2Point b0, S2Point b1,
                         S2Point b1_snapped)
    {
        // If B1' within 90 degrees of B1, no connecting vertex is necessary.
        if (b1_snapped.DotProd(b1) >= 0) return b1;

        // Otherwise we use the point on (B0, B1) that is 90 degrees away from B1'.
        // This is sufficient to ensure conditions (1) and (2).
        S2Point x = S2.RobustCrossProd(b0, b1).CrossProd(b1_snapped).Normalize();
        return (x.DotProd(S2.Interpolate(b0, b1, 0.5)) >= 0) ? x : -x;
    }

    public static VertexId FindFirstVertexId(InputEdgeId input_edge_id, Graph g)
    {
        // A given input edge snaps to a chain of output edges.  To determine which
        // output vertex the source of the given input edge snaps to, we must find
        // the first edge in this chain.
        //
        // The search below takes linear time in the number of edges; it can be done
        // more efficiently if duplicate edges are not being merged and the mapping
        // returned by Graph::GetInputEdgeOrder() is available.  The algorithm would
        // also be much simpler if input_edge_id were known to be degenerate.
        Dictionary<VertexId, int> excess_degree_map = new();
        for (EdgeId e = 0; e < g.NumEdges; ++e)
        {
            var id_set = g.InputEdgeIds(e);
            foreach (InputEdgeId id in id_set)
            {
                if (id == input_edge_id)
                {
                    excess_degree_map[g.GetEdge(e).ShapeId] += 1;
                    excess_degree_map[g.GetEdge(e).EdgeId] -= 1;
                    break;
                }
            }
        }
        if (!excess_degree_map.Any()) return -1;  // Does not exist.

        // Look for the (unique) vertex whose excess degree is +1.
        foreach (var entry in excess_degree_map)
        {
            if (entry.Value == 1) return entry.Key;
        }
        // Otherwise "input_edge_id" must snap to a single degenerate edge.
        System.Diagnostics.Debug.Assert(excess_degree_map.Count == 1);
        return excess_degree_map.First().Key;
    }

    // An input edge may snap to zero, one, or two non-degenerate output edges
    // incident to the reference vertex, consisting of at most one incoming and
    // one outgoing edge.
    //
    // If v_in >= 0, an incoming edge to the reference vertex is present.
    // If v_out >= 0, an outgoing edge from the reference vertex is present.
    public struct EdgeSnap
    {
        public S2Shape.Edge input;
        public VertexId v_in;
        public VertexId v_out;

        public EdgeSnap(S2Shape.Edge input, int v_in = -1, int v_out = -1)
        {
            this.input = input;
            this.v_in = v_in;
            this.v_out = v_out;
        }
    }

    public readonly record struct InputVertexEdge(S2Point Point, EdgeSnap EdgeSnap)
        : IComparable<InputVertexEdge>
    {
        public int CompareTo(InputVertexEdge other) => Point.CompareTo(other.Point);

        public bool Equals(InputVertexEdge other) => other.Point.Equals(this.Point);

        public override int GetHashCode() => Point.GetHashCode();
    }
}
