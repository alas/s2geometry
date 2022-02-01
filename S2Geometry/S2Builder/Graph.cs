using System.Diagnostics.CodeAnalysis;
using System.Runtime.InteropServices;
using System.Text;

namespace S2Geometry;

using SiblingPairs = S2Builder.GraphOptions.SiblingPairs;
using DegenerateEdges = S2Builder.GraphOptions.DegenerateEdges;
using DuplicateEdges = S2Builder.GraphOptions.DuplicateEdges;

public partial class S2Builder
{
    // An S2Builder.Graph represents a collection of snapped edges that is passed
    // to a Layer for assembly.  (Example layers include polygons, polylines, and
    // polygon meshes.)  The Graph object does not own any of its underlying data;
    // it is simply a view of data that is stored elsewhere.  You will only
    // need this interface if you want to implement a new Layer subtype.
    //
    // The graph consists of vertices and directed edges.  Vertices are numbered
    // sequentially starting from zero.  An edge is represented as a pair of
    // vertex ids.  The edges are sorted in lexicographic order, therefore all of
    // the outgoing edges from a particular vertex form a contiguous range.
    //
    // S2Builder.Graph is movable and copyable.  Note that although this class
    // does not own the underlying vertex and edge data, S2Builder guarantees that
    // all Graph objects passed to S2Builder.Layer.Build() methods will remain
    // valid until all layers have been built.
    //
    // TODO(ericv): Consider pulling out the methods that are helper functions for
    // Layer implementations (such as GetDirectedLoops) into s2builderutil_graph.h.
    public class Graph
    {
        private readonly IsFullPolygonPredicate? is_full_polygon_predicate_;

        public GraphOptions Options { get; }

        /// <summary>
        /// Tthe number of vertices in the graph.
        /// </summary>
        public int NumVertices { get; }  // vertices_.Count requires division by 24.

        // Returns the entire set of vertices.
        public List<S2Point> Vertices { get; }

        // Returns the entire set of edges.
        public List<Edge> Edges { get; }

        // Low-level method that returns a vector where each element represents the
        // set of input edge ids that were snapped to a particular output edge.
        public EdgeIdSet InputEdgeIdSetIds { get; }

        // Returns a mapping from an int to a set of input edge ids.
        public IdSetLexicon InputEdgeIdSetLexicon { get; }

        // Low-level method that returns a vector where each element represents the
        // set of labels associated with a particular input edge.  Note that this
        // vector may be empty, which indicates that no labels are present.
        public LabelSet? LabelSetIds { get; }

        // Returns a mapping from a int to a set of labels.
        public IdSetLexicon? LabelSetLexicon { get; }

        // Defines a value larger than any valid InputEdgeId.
        private const int kMaxInputEdgeId = int.MaxValue;

        // The following value of InputEdgeId means that an edge does not
        // corresponds to any input edge.
        private const int kNoInputEdgeId = kMaxInputEdgeId - 1;

        #region Constructors

        // Note that most of the parameters are passed by reference and must
        // exist for the duration of the Graph object.  Notes on parameters:
        // "options":
        //    - the S2Builder.GraphOptions used to build the Graph.  In some cases these
        //      can be different than the options provided by the Layer.
        // "vertices":
        //   - a vector of S2Points indexed by VertexId.
        // "edges":
        //   - a vector of VertexId pairs (sorted in lexicographic order)
        //     indexed by int.
        // "input_edge_id_set_ids":
        //   - a vector indexed by int that allows access to the set of
        //     InputEdgeIds that were mapped to the given edge, by looking up the
        //     returned value (an int) in "input_edge_id_set_lexicon".
        // "input_edge_id_set_lexicon":
        //   - a class that maps an int to a set of InputEdgeIds.
        // "label_set_ids":
        //   - a vector indexed by InputEdgeId that allows access to the set of
        //     labels that were attached to the given input edge, by looking up the
        //     returned value (a LabelSetId) in the "label_set_lexicon".  This
        //     vector may be empty to indicate that no labels are present.
        // "label_set_lexicon":
        //   - a class that maps a LabelSetId to a set of S2Builder::Labels.  (Must
        //     be provided even if no labels are present.)
        // "is_full_polygon_predicate":
        //   - a predicate called to determine whether a graph consisting only of
        //     polygon degeneracies represents the empty polygon or the full polygon
        //     (see s2builder.h for details).
        public Graph(GraphOptions options, List<S2Point> vertices, List<Edge> edges,
            EdgeIdSet input_edge_id_set_ids, IdSetLexicon input_edge_id_set_lexicon,
            LabelSet? label_set_ids, IdSetLexicon? label_set_lexicon,
            IsFullPolygonPredicate? is_full_polygon_predicate)
        {
            Options = options; NumVertices = vertices.Count; Vertices = vertices;
            Edges = edges; InputEdgeIdSetIds = input_edge_id_set_ids;
            InputEdgeIdSetLexicon = input_edge_id_set_lexicon;
            LabelSetIds = label_set_ids;
            LabelSetLexicon = label_set_lexicon;
            is_full_polygon_predicate_ = is_full_polygon_predicate;

            Assert.True(edges.IsSorted());
            Assert.True(edges.Count == input_edge_id_set_ids.Count);
        }

        #endregion

        // Returns the vertex at the given index.
        public S2Point Vertex(int v) => Vertices[v];

        // Returns the total number of edges in the graph.
        public int NumEdges => Edges.Count;

        // Returns the endpoints of the given edge (as vertex indices).
        public Edge GetEdge(int e) => Edges[e];

        // Given an edge (src, dst), returns the reverse edge (dst, src).
        public static Edge Reverse(Edge e) => new(e.EdgeId, e.ShapeId);

        // Returns a vector of edge ids sorted in lexicographic order by
        // (destination, origin).  All of the incoming edges to each vertex form a
        // contiguous subrange of this ordering.
        public EdgeLoop GetInEdgeIds()
        {
            var in_edge_ids = new EdgeLoop(NumEdges);
            in_edge_ids.Iota(0, NumEdges);
            in_edge_ids.Sort((int ai, int bi) => StableLessThan(Reverse(GetEdge(ai)), Reverse(GetEdge(bi)), ai, bi));
            return in_edge_ids;
        }

        // Given a graph such that every directed edge has a sibling, returns a map
        // from EdgeId to the sibling EdgeId.  This method is identical to
        // GetInEdgeIds() except that (1) it requires edges to have siblings, and
        // (2) undirected degenerate edges are grouped together in pairs such that
        // one edge is the sibling of the other.  (The sibling of a directed
        // degenerate edge is itself.)  Handles duplicate edges correctly and is
        // also consistent with GetLeftTurnMap().
        //
        // REQUIRES: An option is chosen that guarantees sibling pairs:
        //     (options.sibling_pairs() == { REQUIRE, CREATE } ||
        //      options.edge_type() == UNDIRECTED)
        public EdgeLoop GetSiblingMap()
        {
            var in_edge_ids = GetInEdgeIds();
            MakeSiblingMap(in_edge_ids);
            return in_edge_ids;
        }

        // Like GetSiblingMap(), butructs the map starting from the vector of
        // incoming edge ids returned by GetInEdgeIds().  (This operation is a no-op
        // unless undirected degenerate edges are present, in which case such edges
        // are grouped together in pairs to satisfy the requirement that every edge
        // must have a sibling edge.)
        public void MakeSiblingMap(EdgeLoop in_edge_ids)
        {
            Assert.True(Options.SiblingPairs_ == SiblingPairs.REQUIRE ||
                   Options.SiblingPairs_ == SiblingPairs.CREATE ||
                   Options.EdgeType_ == EdgeType.UNDIRECTED);
            for (var e = 0; e < NumEdges; ++e)
            {
                Assert.True(GetEdge(e) == Reverse(GetEdge(in_edge_ids[e])));
            }
            if (Options.EdgeType_ == EdgeType.DIRECTED) return;
            if (Options.DegenerateEdges_ == DegenerateEdges.DISCARD) return;

            for (var e = 0; e < NumEdges; ++e)
            {
                var v = GetEdge(e).ShapeId;
                if (GetEdge(e).EdgeId == v)
                {
                    Assert.True(e + 1 < NumEdges);
                    Assert.True(GetEdge(e + 1).ShapeId == v);
                    Assert.True(GetEdge(e + 1).EdgeId == v);
                    Assert.True(in_edge_ids[e] == e);
                    Assert.True(in_edge_ids[e + 1] == e + 1);
                    in_edge_ids[e] = e + 1;
                    in_edge_ids[e + 1] = e;
                    ++e;
                }
            }
        }

        // Returns the set of input edge ids that were snapped to the given
        // edge.  ("Input edge ids" are assigned to input edges sequentially in
        // the order they are added to the builder.)  For example, if input
        // edges 2 and 17 were snapped to edge 12, then input_edge_ids(12)
        // returns a set containing the numbers 2 and 17.  Example usage:
        //
        //   for (var input_edge_id in g.input_edge_ids(e)) { ... }
        //
        // Please note the following:
        //
        //  - When edge chains are simplified, the simplified edge is assigned all
        //    the input edge ids associated with edges of the chain.
        //
        //  - Edges can also have multiple input edge ids due to edge merging
        //    (if DuplicateEdges.MERGE is specified).
        //
        //  - Siblings edges automatically created by EdgeType.UNDIRECTED or
        //    SiblingPairs.CREATE have an empty set of input edge ids.  (However
        //    you can use a LabelFetcher to retrieve the set of labels associated
        //    with both edges of a given sibling pair.)
        public List<int> InputEdgeIds(int e) => InputEdgeIdSetLexicon.IdSet_(InputEdgeIdSetIds[e]);

        // Low-level method that returns an integer representing the entire set of
        // input edge ids that were snapped to the given edge.  The elements of the
        // IdSet can be accessed using input_edge_id_set_lexicon().
        public int InputEdgeIdSetId(int e) => InputEdgeIdSetIds[e];

        // Returns the minimum input edge id that was snapped to this edge, or -1 if
        // no input edges were snapped (see SiblingPairs.CREATE).  This is
        // useful for layers that wish to preserve the input edge ordering as much
        // as possible (e.g., to ensure idempotency).
        public int MinInputEdgeId(int e)
        {
            var id_set = InputEdgeIds(e);
            return (id_set.Count == 0) ? kNoInputEdgeId : 0;
        }

        // Returns a vector containing the minimum input edge id for every edge.
        // If an edge has no input ids, kNoInputEdgeId is used.
        public InputEdgeLoop GetMinInputEdgeIds()
        {
            var min_input_ids = new InputEdgeLoop(NumEdges);
            for (var e = 0; e < NumEdges; ++e)
            {
                min_input_ids[e] = MinInputEdgeId(e);
            }
            return min_input_ids;
        }

        // Returns a vector of EdgeIds sorted by minimum input edge id.  This is an
        // approximation of the input edge ordering.
        public static EdgeLoop GetInputEdgeOrder(InputEdgeLoop min_input_edge_ids)
        {
            var order = new EdgeLoop(min_input_edge_ids.Count);
            order.Iota(0, min_input_edge_ids.Count);
            // Comparison function ensures sort is stable.
            order.Sort((int a, int b) => (min_input_edge_ids[a], a).CompareTo((min_input_edge_ids[b], b)));
            return order;
        }

        // Returns the set of labels associated with a given input edge.  Example:
        //   for (Label label : g.labels(input_edge_id)) { ... }
        // See also LabelFetcher, which returns the labels for a given graph edge.
        public LabelSet Labels(InputEdgeId e) => LabelSetLexicon.IdSet_(LabelSetId(e));

        // Low-level method that returns an integer representing the set of
        // labels associated with a given input edge.  The elements of
        // the IdSet can be accessed using label_set_lexicon().
        public int LabelSetId(int e) => !LabelSetIds.Any()
            ? IdSetLexicon.kEmptySetId
            : LabelSetIds[e];

        // Convenience method that calls is_full_polygon_predicate() to determine
        // whether a graph that consists only of polygon degeneracies represents the
        // empty polygon or the full polygon (see s2builder.h for details).
        public bool IsFullPolygon(out S2Error error) => is_full_polygon_predicate_(this, out error);

        // Returns a method that determines whether a graph that consists only of
        // polygon degeneracies represents the empty polygon or the full polygon
        // (see s2builder.h for details).
        public IsFullPolygonPredicate IsFullPolygonPredicate() => is_full_polygon_predicate_;

        // Returns a map "m" that maps each edge e=(v0,v1) to the following outgoing
        // edge around "v1" in clockwise order.  (This corresponds to making a "left
        // turn" at the vertex.)  By starting at a given edge and making only left
        // turns, you can construct a loop whose interior does not contain any edges
        // in the same connected component.
        //
        // If the incoming and outgoing edges around a vertex do not alternate
        // perfectly (e.g., there are two incoming edges in a row), then adjacent
        // (incoming, outgoing) pairs are repeatedly matched and removed.  This is
        // similar to finding matching parentheses in a string such as "(()())()".
        //
        // For sibling edge pairs, the incoming edge is assumed to immediately
        // follow the outgoing edge in clockwise order.  Thus a left turn is made
        // from an edge to its sibling only if there are no other outgoing edges.
        // With respect to the parentheses analogy, a sibling pair is ")(".
        // Similarly, if there are multiple copies of a sibling edge pair then the
        // duplicate incoming and outgoing edges are sorted in alternating order
        // (e.g., ")()(").
        //
        // Degenerate edges (edges from a vertex to itself) are treated as loops
        // consisting of a single edge.  This avoids the problem of deciding the
        // connectivity and ordering of such edges when they share a vertex with
        // other edges (possibly including other degenerate edges).
        //
        // If it is not possible to make a left turn from every input edge, this
        // method returns false and sets "error" appropriately.  In this situation
        // the left turn map is still valid except that any incoming edge where it
        // is not possible to make a left turn will have its entry set to -1.
        //
        // "in_edge_ids" should be equal to GetInEdgeIds() or GetSiblingMap().
        public bool GetLeftTurnMap(EdgeLoop in_edge_ids, out EdgeLoop left_turn_map, out S2Error error)
        {
            error = S2Error.OK;
            left_turn_map = new EdgeLoop();
            left_turn_map.Fill(-1, NumEdges);
            if (NumEdges == 0) return true;

            // Declare vectors outside the loop to avoid reallocating them each time.
            var v0_edges = new VertexEdgeList();
            var e_in = new EdgeLoop();
            var e_out = new EdgeLoop();

            // Walk through the two sorted arrays of edges (outgoing and incoming) and
            // gather all the edges incident to each vertex.  Then we sort those edges
            // and add an entry to the left turn map from each incoming edge to the
            // immediately following outgoing edge in clockwise order.
            int out_ = 0, in_ = 0;
            var out_edge = GetEdge(out_);
            var in_edge = GetEdge(in_edge_ids[in_]);
            var sentinel = new Edge(NumVertices, NumVertices);
            var min_edge = new[] { out_edge, Reverse(in_edge) }.Min();
            while (min_edge != sentinel)
            {
                // Gather all incoming and outgoing edges around vertex "v0".
                var v0 = min_edge.ShapeId;
                for (; min_edge.ShapeId == v0; min_edge = new[] { out_edge, Reverse(in_edge) }.Min())
                {
                    var v1 = min_edge.EdgeId;
                    // Count the number of copies of "min_edge" in each direction.
                    var out_begin = out_;
                    var in_begin = in_;
                    while (out_edge == min_edge)
                    {
                        out_edge = (++out_ == NumEdges) ? sentinel : GetEdge(out_);
                    }
                    while (Reverse(in_edge) == min_edge)
                    {
                        in_edge = (++in_ == NumEdges) ? sentinel : GetEdge(in_edge_ids[in_]);
                    }
                    if (v0 != v1)
                    {
                        AddVertexEdges(out_begin, out_, in_begin, in_, v1, v0_edges);
                    }
                    else
                    {
                        // Each degenerate edge becomes its own loop.
                        for (; in_begin < in_; ++in_begin)
                        {
                            left_turn_map[in_begin] = in_begin;
                        }
                    }
                }
                if (!v0_edges.Any()) continue;

                // Sort the edges in clockwise order around "v0".
                var min_endpoint = v0_edges.First().Endpoint;
                v0_edges.Sort(1, v0_edges.Count, new CcwComp(min_endpoint, Vertices, v0));
                // Match incoming with outgoing edges.  We do this by keeping a stack of
                // unmatched incoming edges.  We also keep a stack of outgoing edges with
                // no previous incoming edge, and match these at the end by wrapping
                // around circularly to the start of the edge ordering.
                foreach (var e in v0_edges)
                {
                    if (e.Incoming)
                    {
                        e_in.Add(in_edge_ids[e.Index]);
                    }
                    else if (e_in.Any())
                    {
                        left_turn_map[e_in.Last()] = e.Index;
                        e_in.RemoveAt(e_in.Count - 1);
                    }
                    else
                    {
                        e_out.Add(e.Index);  // Matched below.
                    }
                }
                // Pair up additional edges using the fact that the ordering is circular.
                e_out.Reverse();
                for (; e_out.Any() && e_in.Any(); e_out.RemoveAt(e_out.Count - 1), e_in.RemoveAt(e_in.Count - 1))
                {
                    left_turn_map[e_in.Last()] = e_out.Last();
                }
                // We only need to process unmatched incoming edges, since we are only
                // responsible for creating left turn map entries for those edges.
                if (e_in.Any() && error.IsOk())
                {
                    error = new(S2ErrorCode.BUILDER_EDGES_DO_NOT_FORM_LOOPS, "Given edges do not form loops (indegree != outdegree)");
                }
                e_in.Clear();
                e_out.Clear();
                v0_edges.Clear();
            }
            return error.IsOk();
        }

        // Given a set of duplicate outgoing edges (v0, v1) and a set of duplicate
        // incoming edges (v1, v0), this method assigns each edge an integer "rank" so
        // that the edges are sorted in a consistent order with respect to their
        // orderings around "v0" and "v1".  Usually there is just one edge, in which
        // case this is easy.  Sometimes there is one edge in each direction, in which
        // case the outgoing edge is always ordered before the incoming edge.
        //
        // In general, we allow any number of duplicate edges in each direction, in
        // which case outgoing edges are interleaved with incoming edges so as to
        // create as many degenerate (two-edge) loops as possible.  In order to get a
        // consistent ordering around "v0" and "v1", we move forwards through the list
        // of outgoing edges and backwards through the list of incoming edges.  If
        // there are more incoming edges, they go at the beginning of the ordering,
        // while if there are more outgoing edges then they go at the end.
        //
        // For example, suppose there are 2 edges "a,b" from "v0" to "v1", and 4 edges
        // "w,x,y,z" from "v1" to "v0".  Using lower/upper case letters to represent
        // incoming/outgoing edges, the clockwise ordering around v0 would be zyAxBw,
        // and the clockwise ordering around v1 would be WbXaYZ.  (Try making a
        // diagram with each edge as a separate arc.)
        private static void AddVertexEdges(int out_begin, int out_end, int in_begin,
            int in_end, int v1, VertexEdgeList v0_edges)
        {
            var rank = 0;
            // Any extra incoming edges go at the beginning of the ordering.
            while (in_end - in_begin > out_end - out_begin)
            {
                v0_edges.Add(new VertexEdge(true, --in_end, v1, rank++));
            }
            // Next we interleave as many outgoing and incoming edges as possible.
            while (in_end > in_begin)
            {
                v0_edges.Add(new VertexEdge(false, out_begin++, v1, rank++));
                v0_edges.Add(new VertexEdge(true, --in_end, v1, rank++));
            }
            // Any extra outgoing edges to at the end of the ordering.
            while (out_end > out_begin)
            {
                v0_edges.Add(new VertexEdge(false, out_begin++, v1, rank++));
            }
        }

        // Rotates the edges of "loop" if necessary so that the edge(s) with the
        // largest input edge ids are last.  This ensures that when an output loop
        // is equivalent to an input loop, their cyclic edge orders are the same.
        // "min_input_ids" is the output of GetMinInputEdgeIds().
        public static void CanonicalizeLoopOrder(InputEdgeLoop min_input_ids, EdgeLoop loop)
        {
            if (!loop.Any()) return;
            // Find the position of the element with the highest input edge id.  If
            // there are multiple such elements together (i.e., the edge was split
            // into several pieces by snapping it to several vertices), then we choose
            // the last such position in cyclic order (this attempts to preserve the
            // original loop order even when new vertices are added).  For example, if
            // the input edge id sequence is (7, 7, 4, 5, 6, 7) then we would rotate
            // it to obtain (4, 5, 6, 7, 7, 7).

            // The reason that we put the highest-numbered edge last, rather than the
            // lowest-numbered edge first, is that S2Loop.Invert() reverses the loop
            // edge order *except* for the last edge.  For example, the loop ABCD (with
            // edges AB, BC, CD, DA) becomes DCBA (with edges DC, CB, BA, AD).  Note
            // that the last edge is the same except for its direction (DA vs. AD).
            // This has the advantage that if an undirected loop is assembled with the
            // wrong orientation and later inverted (e.g. by S2Polygon.InitOriented),
            // we still end up preserving the original cyclic vertex order.
            var pos = 0;
            var saw_gap = false;
            for (var i = 1; i < loop.Count; ++i)
            {
                var cmp = min_input_ids[loop[i]] - min_input_ids[loop[pos]];
                if (cmp < 0)
                {
                    saw_gap = true;
                }
                else if (cmp > 0 || !saw_gap)
                {
                    pos = i;
                    saw_gap = false;
                }
            }
            if (++pos == loop.Count) pos = 0;  // Convert loop end to loop start.
            loop.RotateInPlace(pos);
        }

        // Sorts the given edge chains (i.e., loops or polylines) by the minimum
        // input edge id of each chains's first edge.  This ensures that when the
        // output consists of multiple loops or polylines, they are sorted in the
        // same order as they were provided in the input.
        public static void CanonicalizeVectorOrder(InputEdgeLoop min_input_ids, List<EdgeLoop> chains) =>
            // Comparison function ensures sort is stable.
            chains.Sort((EdgeLoop a, EdgeLoop b) => (min_input_ids[a[0]], a[0]).CompareTo((min_input_ids[b[0]], b[0])));

        // Builds loops from a set of directed edges, turning left at each vertex
        // until either a repeated vertex (for LoopType.SIMPLE) or a repeated edge
        // (for LoopType.CIRCUIT) is found.  (Use LoopType.SIMPLE if you intend to
        // construct an S2Loop.)
        //
        // Each loop is represented as a sequence of edges.  The edge ordering and
        // loop ordering are automatically canonicalized in order to preserve the
        // input ordering as much as possible.  Loops are non-crossing provided that
        // the graph contains no crossing edges.  If some edges cannot be turned
        // into loops, returns false and sets "error" appropriately.
        //
        // If any degenerate edges are present, then each such edge is treated as a
        // separate loop.  This is mainly useful in conjunction with
        // options.degenerate_edges() == DISCARD_EXCESS, in order to build polygons
        // that preserve degenerate geometry.
        //
        // REQUIRES: options.degenerate_edges() == {DISCARD, DISCARD_EXCESS}
        // REQUIRES: options.edge_type() == DIRECTED
        public bool GetDirectedLoops(LoopType loop_type, List<EdgeLoop> loops, out S2Error error)
        {
            Assert.True(Options.DegenerateEdges_ == DegenerateEdges.DISCARD ||
                   Options.DegenerateEdges_ == DegenerateEdges.DISCARD_EXCESS);
            Assert.True(Options.EdgeType_ == EdgeType.DIRECTED);

            if (!GetLeftTurnMap(GetInEdgeIds(), out var left_turn_map, out error)) return false;
            var min_input_ids = GetMinInputEdgeIds();

            // If we are breaking loops at repeated vertices, we maintain a map from
            // int to its position in "path".
            int[]? path_index = null;
            if (loop_type == LoopType.SIMPLE)
            {
                path_index = new int[NumVertices].Fill(-1);
            }

            // Visit edges in arbitrary order, and try to build a loop from each edge.
            var path = new EdgeLoop();
            for (var start = 0; start < NumEdges; ++start)
            {
                if (left_turn_map[start] < 0) continue;

                // Build a loop by making left turns at each vertex until we return to
                // "start".  We use "left_turn_map" to keep track of which edges have
                // already been visited by setting its entries to -1 as we go along.  If
                // we are building vertex cycles, then whenever we encounter a vertex that
                // is already part of the path, we "peel off" a loop by removing those
                // edges from the path so far.
                int next;
                for (var e = start; left_turn_map[e] >= 0; e = next)
                {
                    path.Add(e);
                    next = left_turn_map[e];
                    left_turn_map[e] = -1;
                    if (loop_type == LoopType.SIMPLE)
                    {
                        path_index[GetEdge(e).ShapeId] = path.Count - 1;
                        var loop_start = path_index[GetEdge(e).EdgeId];
                        if (loop_start < 0) continue;
                        // Peel off a loop from the path.
                        var loop = path.Skip(loop_start).ToList();
                        path.RemoveRange(loop_start, path.Count - loop_start);
                        foreach (var e2 in loop) path_index[GetEdge(e2).ShapeId] = -1;
                        CanonicalizeLoopOrder(min_input_ids, loop);
                        loops.Add(loop);
                    }
                }
                if (loop_type == LoopType.SIMPLE)
                {
                    Assert.True(!path.Any());  // Invariant.
                }
                else
                {
                    CanonicalizeLoopOrder(min_input_ids, path);
                    loops.Add(path);
                    path.Clear();
                }
            }
            CanonicalizeVectorOrder(min_input_ids, loops);
            return true;
        }
        public bool GetDirectedComponents(DegenerateBoundaries degenerate_boundaries, List<DirectedComponent> components, out S2Error error)
        {
            Assert.True(Options.DegenerateEdges_ == DegenerateEdges.DISCARD ||
                   (Options.DegenerateEdges_ == DegenerateEdges.DISCARD_EXCESS &&
                    degenerate_boundaries == DegenerateBoundaries.KEEP));
            Assert.True(Options.SiblingPairs_ == SiblingPairs.REQUIRE ||
                   Options.SiblingPairs_ == SiblingPairs.CREATE);
            Assert.True(Options.EdgeType_ == EdgeType.DIRECTED);  // Implied by above.

            var sibling_map = GetSiblingMap();
            if (!GetLeftTurnMap(sibling_map, out var left_turn_map, out error)) return false;
            var min_input_ids = GetMinInputEdgeIds();
            var frontier = new EdgeLoop();  // Unexplored sibling edges.

            // A map from EdgeId to the position of that edge in "path".  Only needed if
            // degenerate boundaries are being discarded.
            var path_index = new int[NumEdges];
            if (degenerate_boundaries == DegenerateBoundaries.DISCARD)
            {
                path_index.Fill(-1);
            }
            for (var start = 0; start < NumEdges; ++start)
            {
                if (left_turn_map[start] < 0) continue;  // Already used.

                // Build a connected component by keeping a stack of unexplored siblings
                // of the edges used so far.
                var component = new DirectedComponent();
                frontier.Add(start);
                while (frontier.Any())
                {
                    var e = frontier.Last();
                    frontier.RemoveAt(frontier.Count - 1);
                    if (left_turn_map[e] < 0) continue;  // Already used.

                    // Build a path by making left turns at each vertex until we complete a
                    // loop.  Whenever we encounter an edge that is a sibling of an edge
                    // that is already on the path, we "peel off" a loop consisting of any
                    // edges that were between these two edges.
                    var path = new EdgeLoop();
                    EdgeId next;
                    for (; left_turn_map[e] >= 0; e = next)
                    {
                        path.Add(e);
                        next = left_turn_map[e];
                        left_turn_map[e] = -1;
                        // If the sibling hasn't been visited yet, add it to the frontier.
                        var sibling = sibling_map[e];
                        if (left_turn_map[sibling] >= 0)
                        {
                            frontier.Add(sibling);
                        }
                        if (degenerate_boundaries == DegenerateBoundaries.DISCARD)
                        {
                            path_index[e] = path.Count - 1;
                            var sibling_index = path_index[sibling];
                            if (sibling_index < 0) continue;

                            // Common special case: the edge and its sibling are adjacent, in
                            // which case we can simply remove them from the path and continue.
                            if (sibling_index == path.Count - 2)
                            {
                                path.Capacity = sibling_index;
                                // We don't need to update "path_index" for these two edges
                                // because both edges of the sibling pair have now been used.
                                continue;
                            }
                            // Peel off a loop from the path.
                            var loop = path.Skip(sibling_index + 1).Take(path.Count - sibling_index - 2).ToList();
                            path.RemoveRange(sibling_index, path.Count - sibling_index);
                            // Mark the edges that are no longer part of the path.
                            foreach (var e2 in loop) path_index[e2] = -1;
                            CanonicalizeLoopOrder(min_input_ids, loop);
                            component.Add(loop);
                        }
                    }
                    // Mark the edges that are no longer part of the path.
                    if (degenerate_boundaries == DegenerateBoundaries.DISCARD)
                    {
                        foreach (var e2 in path) path_index[e2] = -1;
                    }
                    CanonicalizeLoopOrder(min_input_ids, path);
                    component.Add(path);
                }
                CanonicalizeVectorOrder(min_input_ids, component);
                components.Add(component);
            }
            // Sort the components to correspond to the input edge ordering.
            components.Sort((DirectedComponent a, DirectedComponent b) => min_input_ids[a[0][0]].CompareTo(min_input_ids[b[0][0]]));
            return true;
        }

        // Builds loops from a set of undirected edges, turning left at each vertex
        // until either a repeated vertex (for LoopType.SIMPLE) or a repeated edge
        // (for LoopType.CIRCUIT) is found.  The loops are further grouped into
        // "components" such that all the loops in a component are connected by
        // shared vertices.  Finally, the loops in each component are divided into
        // two "complements" such that every edge in one complement is the sibling
        // of an edge in the other complement.  This corresponds to the fact that
        // given any set of non-crossing undirected loops, there are exactly two
        // possible interpretations of the region that those loops represent (where
        // one possibility is the complement of the other).  This method does not
        // attempt to resolve this ambiguity, but instead returns both possibilities
        // for each connected component and lets the client choose among them.
        //
        // This method is used to build single polygons.  (Use GetDirectedComponents
        // to build polygon meshes, even when the input edges are undirected.)  To
        // convert the output of this method into a polygon, the client must choose
        // one complement from each component such that the entire set of loops is
        // oriented consistently (i.e., they define a region such that the interior
        // of the region is always on the left).  The non-chosen complements form
        // another set of loops that are also oriented consistently but represent
        // the complementary region on the sphere.  Finally, the client needs to
        // choose one of these two sets of loops based on heuristics (e.g., the area
        // of each region), since both sets of loops are equally valid
        // interpretations of the input.
        //
        // Each loop is represented as a sequence of edges.  The edge ordering and
        // loop ordering are automatically canonicalized in order to preserve the
        // input ordering as much as possible.  Loops are non-crossing provided that
        // the graph contains no crossing edges.  If some edges cannot be turned
        // into loops, returns false and sets "error" appropriately.
        //
        // REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
        // REQUIRES: options.edge_type() == UNDIRECTED
        // REQUIRES: options.siblings_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
        //           [since REQUIRE, CREATE convert the edge_type() to DIRECTED]
        public bool GetUndirectedComponents(LoopType loop_type, List<UndirectedComponent> components, out S2Error error)
        {
            Assert.True(Options.DegenerateEdges_ == DegenerateEdges.DISCARD ||
                   Options.DegenerateEdges_ == DegenerateEdges.DISCARD_EXCESS);
            Assert.True(Options.EdgeType_ == EdgeType.UNDIRECTED);

            var sibling_map = GetInEdgeIds();
            if (!GetLeftTurnMap(sibling_map, out var left_turn_map, out error)) return false;
            MakeSiblingMap(sibling_map);
            var min_input_ids = GetMinInputEdgeIds();

            // A stack of unexplored sibling edges.  Each sibling edge has a "slot"
            // (0 or 1) that indicates which of the two complements it belongs to.
            var frontier = new List<(int start, int slot)>();

            // If we are breaking loops at repeated vertices, we maintain a map from
            // int to its position in "path".
            int[]? path_index = null;
            if (loop_type == LoopType.SIMPLE) path_index = new int[NumVertices].Fill(-1);

            for (var min_start = 0; min_start < NumEdges; ++min_start)
            {
                if (left_turn_map[min_start] < 0) continue;  // Already used.

                // Build a connected component by keeping a stack of unexplored siblings
                // of the edges used so far.
                var component = new UndirectedComponent(() => new List<List<int>>());
                frontier.Add((min_start, 0));
                while (frontier.Any())
                {
                    var (start, slot) = frontier.Last();
                    frontier.RemoveAt(frontier.Count - 1);
                    if (left_turn_map[start] < 0) continue;  // Already used.

                    // Build a path by making left turns at each vertex until we return to
                    // "start".  We use "left_turn_map" to keep track of which edges have
                    // already been visited, and which complement they were assigned to, by
                    // setting its entries to negative values as we go along.
                    var path = new EdgeLoop();
                    int next;
                    for (var e = start; left_turn_map[e] >= 0; e = next)
                    {
                        path.Add(e);
                        next = left_turn_map[e];
                        left_turn_map[e] = MarkEdgeUsed(slot);
                        // If the sibling hasn't been visited yet, add it to the frontier.
                        var sibling = sibling_map[e];
                        if (left_turn_map[sibling] >= 0)
                        {
                            frontier.Add((sibling, 1 - slot));
                        }
                        else if (left_turn_map[sibling] != MarkEdgeUsed(1 - slot))
                        {
                            // Two siblings edges can only belong the same complement if the
                            // given undirected edges do not form loops.
                            error = new(S2ErrorCode.BUILDER_EDGES_DO_NOT_FORM_LOOPS, "Given undirected edges do not form loops");
                            return false;
                        }
                        if (loop_type == LoopType.SIMPLE)
                        {
                            // Whenever we encounter a vertex that is already part of the path,
                            // we "peel off" a loop by removing those edges from the path.
                            path_index[GetEdge(e).ShapeId] = path.Count - 1;
                            var loop_start = path_index[GetEdge(e).EdgeId];
                            if (loop_start < 0) continue;
                            var loop = path.Skip(loop_start).ToList();
                            path.RemoveRange(loop_start, path.Count - loop_start);
                            // Mark the vertices that are no longer part of the path.
                            foreach (var e2 in loop) path_index[GetEdge(e2).ShapeId] = -1;
                            CanonicalizeLoopOrder(min_input_ids, loop);
                            component[slot].Add(loop);
                        }
                    }
                    if (loop_type == LoopType.SIMPLE)
                    {
                        Assert.True(!path.Any());  // Invariant.
                    }
                    else
                    {
                        CanonicalizeLoopOrder(min_input_ids, path);
                        component[slot].Add(path);
                    }
                }
                CanonicalizeVectorOrder(min_input_ids, component[0]);
                CanonicalizeVectorOrder(min_input_ids, component[1]);
                // To save some work in S2PolygonLayer, we swap the two loop sets of the
                // component so that the loop set whose first loop most closely follows
                // the input edge ordering is first.  (If the input was a valid S2Polygon,
                // then this component will contain normalized loops.)
                if (min_input_ids[component[0][0][0]] > min_input_ids[component[1][0][0]])
                {
                    var tmp = component[0];
                    component[0] = component[1];
                    component[1] = tmp;
                }
                components.Add(component);
            }
            // Sort the components to correspond to the input edge ordering.
            components.Sort((UndirectedComponent a, UndirectedComponent b) => min_input_ids[a[0][0][0]].CompareTo(min_input_ids[b[0][0][0]]));
            return true;
        }

        // Encodes the index of one of the two complements of each component
        // (a.k.a. the "slot", either 0 or 1) as a negative EdgeId.
        private static int MarkEdgeUsed(int slot) => -1 - slot;

        // Builds polylines from a set of edges.  If "polyline_type" is PATH, then
        // only vertices of indegree and outdegree 1 (or degree 2 in the case of
        // undirected edges) will appear in the interior of polylines.  This
        // essentially generates one polyline for each edge chain in the graph.  If
        // "polyline_type" is WALK, then polylines may pass through the same vertex
        // or even the same edge multiple times (if duplicate edges are present),
        // and each polyline will be as long as possible.  This option is useful for
        // reconstructing a polyline that has been snapped to a lower resolution,
        // since snapping can cause edges to become identical.
        //
        // This method attempts to preserve the input edge ordering in order to
        // implement idempotency, even when there are repeated edges or loops.  This
        // is true whether directed or undirected edges are used.  Degenerate edges
        // are also handled appropriately.
        //
        // REQUIRES: options.sibling_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
        public List<EdgePolyline> GetPolylines(PolylineType polyline_type)
        {
            Assert.True(Options.SiblingPairs_ == SiblingPairs.DISCARD ||
                   Options.SiblingPairs_ == SiblingPairs.DISCARD_EXCESS ||
                   Options.SiblingPairs_ == SiblingPairs.KEEP);
            var builder = new PolylineBuilder(this);
            if (polyline_type == PolylineType.PATH)
            {
                return builder.BuildPaths();
            }
            else
            {
                return builder.BuildWalks();
            }
        }

        ////////////////////////////////////////////////////////////////////////
        //////////////// Helper Functions for Creating Graphs //////////////////

        // Given an unsorted collection of edges, transform them according to the
        // given set of S2Builder.GraphOptions.  This includes actions such as discarding
        // degenerate edges; merging duplicate edges; and canonicalizing sibling
        // edge pairs in several possible ways (e.g. discarding or creating them).
        // The output is suitable for passing to the Graph constructor.
        //
        // If options.edge_type() == EdgeType.UNDIRECTED, then all input edges
        // should already have been transformed into a pair of directed edges.
        //
        // "input_ids" is a vector of the same length as "edges" that indicates
        // which input edges were snapped to each edge.  This vector is also updated
        // appropriately as edges are discarded, merged, etc.
        //
        // Note that "options" may be modified by this method: in particular, if
        // edge_type() is UNDIRECTED and sibling_pairs() is CREATE or REQUIRE, then
        // half of the edges in each direction will be discarded and edge_type()
        // will be changed to DIRECTED (see S2Builder::GraphOptions::SiblingPairs).
        //
        // If "tracker" is provided then the memory usage of this method is tracked
        // and an error is returned if the specified memory limit would be exceeded.
        // This option requires that "new_edges" and "new_input_edge_id_set_ids" are
        // already being tracked, i.e. their current memory usage is reflected in
        // "tracker".  Note that "id_set_lexicon" typically uses a negligible amount
        // of memory and is not tracked.
        public static void ProcessEdges(GraphOptions options, List<Edge> edges,
            List<int> input_ids, IdSetLexicon id_set_lexicon,
            out S2Error error, S2MemoryTracker.Client tracker = null)
        {
            error = S2Error.OK;
            // Graph::EdgeProcessor uses 8 bytes per input edge (out_edges_ and
            // in_edges_) plus 12 bytes per output edge (new_edges_, new_input_ids_).
            // For simplicity we assume that num_input_edges == num_output_edges, since
            // Graph:EdgeProcessor does not increase the number of edges except possibly
            // in the case of SiblingPairs::CREATE (which we ignore).
            //
            //  vector<EdgeId> out_edges_;                 // Graph::EdgeProcessor
            //  vector<EdgeId> in_edges_;                  // Graph::EdgeProcessor
            //  vector<Edge> new_edges_;                   // Graph::EdgeProcessor
            //  vector<InputEdgeIdSetId> new_input_ids_;   // Graph::EdgeProcessor
            //
            // EdgeProcessor discards the "edges" and "input_ids" vectors and replaces
            // them with new vectors that could be larger or smaller.  To handle this
            // correctly, we untally these vectors now and retally them at the end.
            var kFinalPerEdge = Marshal.SizeOf(typeof(Edge)) + sizeof(InputEdgeIdSetId);
            var kTempPerEdge = kFinalPerEdge + 2 * sizeof(EdgeId);
            if (tracker != null)
            {
                tracker.TallyTemp(edges.Count * kTempPerEdge);
                tracker.Tally(-edges.Capacity * kFinalPerEdge);
            }
            if (tracker == null || tracker.ok())
            {
                var processor = new EdgeProcessor(options, edges, input_ids, id_set_lexicon);
                processor.Run(out error);
            }
            // Certain values of sibling_pairs() discard half of the edges and change
            // the edge_type() to DIRECTED (see the description of S2Builder.GraphOptions).
            if (options.SiblingPairs_ == SiblingPairs.REQUIRE ||
                options.SiblingPairs_ == SiblingPairs.CREATE)
            {
                options.EdgeType_ = EdgeType.DIRECTED;
            }
            if (tracker != null && !tracker.Tally(edges.Capacity * kFinalPerEdge))
            {
                error = tracker.error();
            }
        }

        // Given a set of vertices and edges, removes all vertices that do not have
        // any edges and returned the new, minimal set of vertices.  Also updates
        // each edge in "edges" to correspond to the new vertex numbering.  (Note
        // that this method does *not* merge duplicate vertices, it simply removes
        // vertices of degree zero.)
        //
        // The new vertex ordering is a subsequence of the original ordering,
        // therefore if the edges were lexicographically sorted before calling this
        // method then they will still be sorted after calling this method.
        //
        // The extra argument "tmp" points to temporary storage used by this method.
        // All calls to this method from a single thread can reuse the same
        // temporary storage.  It should initially point to an empty vector.  This
        // can make a big difference to efficiency when this method is called many
        // times (e.g. to extract the vertices for different layers), since the
        // incremental running time for each layer becomes O(edges.Count) rather
        // than O(vertices.Count + edges.Count).
        public static List<S2Point> FilterVertices(List<S2Point> vertices, List<Edge> edges, List<int> tmp)
        {
            // Gather the vertices that are actually used.
            var usedTmp = new SortedSet<int>();
            // used.Capacity = 2 * edges.Count;
            foreach (var e in edges)
            {
                usedTmp.Add(e.ShapeId);
                usedTmp.Add(e.EdgeId);
            }
            var used = usedTmp.ToList();

            // Build the list of new vertices, and generate a map from old vertex id to
            // new vertex id.
            var vmap = tmp;
            vmap.Capacity = vertices.Count;
            var new_vertices = new List<S2Point>(used.Count);
            for (var i = 0; i < used.Count; ++i)
            {
                new_vertices[i] = vertices[used[i]];
                vmap[used[i]] = i;
            }
            // Update the edges.
            for (var i = 0; i < edges.Count; i++)
            {
                var e = edges[i];
                edges[i] = new Edge(vmap[e.ShapeId], vmap[e.EdgeId]);
            }
            return new_vertices;
        }

        // A comparison function that allows stable sorting with sort (which is
        // fast but not stable).  It breaks ties between equal edges by comparing
        // their edge ids.
        public static int StableLessThan(Edge a, Edge b, int ai, int bi)
        {
            if (a.ShapeId.CompareTo(b.ShapeId) != 0) return a.ShapeId.CompareTo(b.ShapeId);
            if (a.EdgeId.CompareTo(b.EdgeId) != 0) return a.EdgeId.CompareTo(b.EdgeId);
            return ai.CompareTo(bi);  // Stable sort.
        }

        // Constructs a new graph with the given GraphOptions and containing the
        // given edges.  Each edge is associated with a (possibly empty) set of
        // input edges as specified by new_input_edge_id_set_ids (which must be the
        // same length as "new_edges") and the given IdSetLexicon (which allows
        // looking up the set of input edge ids associated with a graph edge).
        // Finally, the subgraph may also specify a new IsFullPolygonPredicate
        // (which is used to distinguish an empty polygon possibly with degenerate
        // shells from a full polygon possibly with degenerate holes).
        //
        // The output graph has the same set of vertices and edge labels as the
        // input graph (noting that edge labels are associated with *input* edges,
        // not graph edges).
        //
        // If new_options.edge_type() is UNDIRECTED then note the following:
        //
        //  - If this->options().edge_type() is DIRECTED then each input edge will
        //    be transformed into a pair of directed edges before calling
        //    ProcessEdges() above.
        //
        //  - If new_options.sibling_pairs() is CREATE or REQUIRE then ProcessEdges()
        //    will discard half of the edges in each direction and change edge_type()
        //    to DIRECTED (see S2Builder::GraphOptions::SiblingPairs).
        //
        // If "tracker" is provided then the memory usage of this method is tracked
        // and an error is returned if the specified memory limit would be exceeded.
        // This option requires that "new_edges" and "new_input_edge_id_set_ids" are
        // already being tracked, i.e. their current memory usage is reflected in
        // "tracker".  Note that "id_set_lexicon" typically uses a negligible amount
        // of memory and is not tracked.
        public Graph? MakeSubgraph(
            GraphOptions new_options, List<Edge> new_edges,
            List<InputEdgeIdSetId> new_input_edge_id_set_ids,
            IdSetLexicon new_input_edge_id_set_lexicon,
            IsFullPolygonPredicate is_full_polygon_predicate,
            out S2Error error, S2MemoryTracker.Client tracker = null)
        {
            if (Options.EdgeType_ == EdgeType.DIRECTED &&
      new_options.EdgeType_ == EdgeType.UNDIRECTED) {
    // Create a reversed edge for every edge.
    int n = new_edges.Count;
    if (tracker == null) {
      new_edges.Capacity = 2 * n;
        new_input_edge_id_set_ids.Capacity = 2 * n;
    } else if (!tracker.AddSpaceExact(new_edges, n) ||
               !tracker.AddSpaceExact(new_input_edge_id_set_ids, n)) {
      error = tracker.error();
      return null;
}
    for (int i = 0; i<n; ++i) {
      new_edges.Add(Graph.Reverse(new_edges[i]));
new_input_edge_id_set_ids.Add(IdSetLexicon.kEmptySetId);
    }
  }
  Graph.ProcessEdges(new_options, new_edges, new_input_edge_id_set_ids,
                      new_input_edge_id_set_lexicon, out error, tracker);
if (tracker != null && !tracker.ok()) return null;  // Graph would be invalid.
return new Graph(new_options, Vertices, new_edges, new_input_edge_id_set_ids,
             new_input_edge_id_set_lexicon, LabelSetIds,
             LabelSetLexicon, is_full_polygon_predicate);
}

        // Indicates whether loops should be simple cycles (no repeated vertices) or
        // circuits (which allow repeated vertices but not repeated edges).  In
        // terms of how the loops are built, this corresponds to closing off a loop
        // at the first repeated vertex vs. the first repeated edge.
        public enum LoopType { SIMPLE, CIRCUIT }

        // Indicates whether polylines should be "paths" (which don't allow
        // duplicate vertices, except possibly the first and last vertex) or
        // "walks" (which allow duplicate vertices and edges).
        public enum PolylineType { PATH, WALK }

        // Builds loops from a set of directed edges, turning left at each vertex
        // until a repeated edge is found (i.e., LoopType.CIRCUIT).  The loops are
        // further grouped into connected components, where each component consists
        // of one or more loops connected by shared vertices.
        //
        // This method is used to build polygon meshes from directed or undirected
        // input edges.  To convert the output of this method into a mesh, the
        // client must determine how the loops in different components are related
        // to each other: for example, several loops from different components may
        // bound the same region on the sphere, in which case all of those loops are
        // combined into a single polygon.  (See S2ShapeUtil.BuildPolygonBoundaries
        // and s2builderutil.LaxPolygonVectorLayer for details.)
        //
        // Note that loops may include both edges of a sibling pair.  When several
        // such edges are connected in a chain or a spanning tree, they form a
        // zero-area "filament".  The entire loop may be a filament (i.e., a
        // degenerate loop with an empty interior), or the loop may have have
        // non-empty interior with several filaments that extend inside it, or the
        // loop may consist of several "holes" connected by filaments.  These
        // filaments do not change the interior of any loop, so if you are only
        // interested in point containment then they can safely be removed by
        // setting the "degenerate_boundaries" parameter to DISCARD.  (They can't be
        // removed by setting (options.sibling_pairs() == DISCARD) because the two
        // siblings might belong to different polygons of the mesh.)  Note that you
        // can prevent multiple copies of sibling pairs by specifying
        // options.duplicate_edges() == MERGE.
        //
        // Each loop is represented as a sequence of edges.  The edge ordering and
        // loop ordering are automatically canonicalized in order to preserve the
        // input ordering as much as possible.  Loops are non-crossing provided that
        // the graph contains no crossing edges.  If some edges cannot be turned
        // into loops, returns false and sets "error" appropriately.
        //
        // REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
        //           (but requires DISCARD if degenerate_boundaries == DISCARD)
        // REQUIRES: options.sibling_pairs() == { REQUIRE, CREATE }
        //           [i.e., every edge must have a sibling edge]
        public enum DegenerateBoundaries { DISCARD, KEEP }

        // A struct for sorting the incoming and outgoing edges around a vertex "v0".
        public readonly struct VertexEdge
        {
            public VertexEdge(bool _incoming, int _index, int _endpoint, int _rank)
            {
                Incoming = _incoming; Index = _index;
                Endpoint = _endpoint; Rank = _rank;
            }
            public readonly bool Incoming;             // Is this an incoming edge to "v0"?
            public readonly int Index;       // Index of this edge in "edges_" or "in_edge_ids"
            public readonly int Endpoint;  // The other (not "v0") endpoint of this edge
            public readonly int Rank;                // Secondary key for edges with the same endpoint
        }

        public class VertexEdgeList : List<VertexEdge> { }

        // A class that maps vertices to their outgoing edge ids.  Example usage:
        //   VertexOutMap out(g);
        //   for (var e in out.edge_ids(v)) { ... }
        //   for (var edge in out.edges(v)) { ... }
        public class VertexOutMap
        {
            public VertexOutMap(Graph g) => Init(g);
            public void Init(Graph g)
            {
                edges_ = g.Edges;
                edge_begins_.Capacity = g.NumVertices + 1;
                var e = 0;
                for (var v = 0; v <= g.NumVertices; ++v)
                {
                    while (e < g.NumEdges && g.GetEdge(e).ShapeId < v) ++e;
                    edge_begins_.Add(e);
                }
            }

            public int Degree(int v) => EdgeIds(v).Count();
            public IEnumerable<Edge> Edges(int v) => edges_.Skip(edge_begins_[v]).Take(edge_begins_[v + 1]);
            public IEnumerable<int> EdgeIds(int v) => edge_begins_.Skip(v).Take(1);
            //return new VertexOutEdgeIds(edge_begins_[v], edge_begins_[v + 1]);

            // Return the edges (or edge ids) between a specific pair of vertices.
            public IEnumerable<Edge> Edges(int v0, int v1)
            {
                var val = new Edge(v0, v1);
                var lb = edges_.GetLowerBound(val, edge_begins_[v0], edge_begins_[v0 + 1]);
                var ub = edges_.GetUpperBound(val, edge_begins_[v0], edge_begins_[v0 + 1]);
                return edges_.Skip(lb).Take(ub - lb);
                /*var range = equal_range(edges_ + edge_begins_[v0],
                    edges_ + edge_begins_[v0 + 1],
                    new Edge(v0, v1));
                return new VertexOutEdges(range.Item1, range.Item2);*/
            }
            public List<int> EdgeIds(int v0, int v1)
            {
                var val = new Edge(v0, v1);
                var lb = edges_.GetLowerBound(val, edge_begins_[v0], edge_begins_[v0 + 1]);
                var ub = edges_.GetUpperBound(val, edge_begins_[v0], edge_begins_[v0 + 1]);
                return edge_begins_.Skip(lb).Take(ub - lb).ToList();
                /*var range = equal_range(edges_ + edge_begins_[v0],
                    edges_ + edge_begins_[v0 + 1],
                    new Edge(v0, v1));
                return new VertexOutEdges(range.Item1-edges_, range.Item2-edges_);*/
            }

            private List<Edge> edges_;
            private readonly EdgeLoop edge_begins_ = new();
        }

        // A class that maps vertices to their incoming edge ids.  Example usage:
        //   VertexInMap in(g);
        //   for (var e in in.edge_ids(v)) { ... }
        public class VertexInMap
        {
            public VertexInMap(Graph g) => Init(g);
            public void Init(Graph g)
            {
                InEdgeIds = g.GetInEdgeIds();
                in_edge_begins_.Capacity = g.NumVertices + 1;
                var e = 0;
                for (var v = 0; v <= g.NumVertices; ++v)
                {
                    while (e < g.NumEdges && g.GetEdge(InEdgeIds[e]).EdgeId < v) ++e;
                    in_edge_begins_.Add(e);
                }
            }

            public int Degree(int v) => EdgeIds(v).Count();
            public IEnumerable<int> EdgeIds(int v) => InEdgeIds.Skip(in_edge_begins_[v]).Take(in_edge_begins_[v + 1] - in_edge_begins_[v]);

            // Returns a sorted vector of all incoming edges (see GetInEdgeIds).
            public EdgeLoop InEdgeIds { get; private set; }

            private readonly EdgeLoop in_edge_begins_ = new();
        }

        // Convenience class to return the set of labels associated with a given
        // graph edge.  Note that due to snapping, one graph edge may correspond to
        // several different input edges and will have all of their labels.
        // This class is the preferred way to retrieve edge labels.
        //
        // The reason this is a class rather than a graph method is because for
        // undirected edges, we need to fetch the labels associated with both
        // siblings.  This is because only the original edge of the sibling pair has
        // labels; the automatically generated sibling edge does not.
        public class LabelFetcher
        {
            public LabelFetcher(Graph g, EdgeType edge_type) => Init(g, edge_type);

            // Prepares to fetch labels associated with the given edge type.  For
            // EdgeType.UNDIRECTED, labels associated with both edges of the sibling
            // pair will be returned.  "edge_type" is a parameter (rather than using
            // g.options().edge_type()) so that clients can explicitly control whether
            // labels from one or both siblings are returned.
            public void Init(Graph g, EdgeType edge_type)
            {
                g_ = g;
                edge_type_ = edge_type;
                if (edge_type == EdgeType.UNDIRECTED) sibling_map_ = g.GetSiblingMap();
            }

            // Returns the set of labels associated with edge "e" (and also the labels
            // associated with the sibling of "e" if edge_type() is UNDIRECTED).
            // Labels are sorted and duplicate labels are automatically removed.
            //
            // This method uses an output parameter rather than returning by value in
            // order to avoid allocating a new vector on every call to this method.
            public void Fetch(int e, List<int> labels)
            {
                labels.Clear();
                foreach (var input_edge_id in g_.InputEdgeIds(e))
                {
                    foreach (var label in g_.Labels(input_edge_id))
                    {
                        labels.AddSortedUnique(label);
                    }
                }
                if (edge_type_ == EdgeType.UNDIRECTED)
                {
                    foreach (var input_edge_id in g_.InputEdgeIds(sibling_map_[e]))
                    {
                        foreach (var label in g_.Labels(input_edge_id))
                        {
                            labels.AddSortedUnique(label);
                        }
                    }
                }
            }

            private Graph g_;
            private EdgeType edge_type_;
            private EdgeLoop sibling_map_;
        }

        private class CcwComp : IComparer<VertexEdge>
        {
            private readonly int min_endpoint_;
            private readonly List<S2Point> vertices_;
            private readonly int v0;

            public CcwComp(int min_endpoint, List<S2Point> vertices, int v0)
            {
                min_endpoint_ = min_endpoint;
                vertices_ = vertices;
                this.v0 = v0;
            }

            public int Compare([AllowNull] VertexEdge a, [AllowNull] VertexEdge b)
            {
                if (a.Endpoint == b.Endpoint) return a.Rank.CompareTo(b.Rank);
                if (a.Endpoint == min_endpoint_) return 1;
                if (b.Endpoint == min_endpoint_) return -1;
                var res = !S2Pred.OrderedCCW(vertices_[a.Endpoint], vertices_[b.Endpoint],
                    vertices_[min_endpoint_], vertices_[v0]);
                return res ? 1 : -1;
            }
        }

        private class EdgeProcessor
        {
            public EdgeProcessor(GraphOptions options, List<Edge> edges, List<int> input_ids, IdSetLexicon id_set_lexicon)
            {
                options_ = options; edges_ = edges;
                input_ids_ = input_ids; id_set_lexicon_ = id_set_lexicon;
                out_edges_ = new EdgeLoop(edges_.Count); in_edges_ = new EdgeLoop(edges_.Count);


                // Sort the outgoing and incoming edges in lexigraphic order.  We use a
                // stable sort to ensure that each undirected edge becomes a sibling pair,
                // even if there are multiple identical input edges.
                out_edges_.Iota(0, edges_.Count);
                out_edges_.Sort((int a, int b) => StableLessThan(edges_[a], edges_[b], a, b));
                in_edges_.Iota(0, edges_.Count);
                in_edges_.Sort((int a, int b) => StableLessThan(Reverse(edges_[a]), Reverse(edges_[b]), a, b));
                new_edges_.Fill(() => new(0, 0), edges_.Count);
                new_input_ids_.Fill(0, edges_.Count);
            }
            public void Run(out S2Error error)
            {
                error = S2Error.OK;
                var num_edges = edges_.Count;
                if (num_edges == 0)
                {
                    return;
                }

                // Walk through the two sorted arrays performing a merge join.  For each
                // edge, gather all the duplicate copies of the edge in both directions
                // (outgoing and incoming).  Then decide what to do based on "options_" and
                // how many copies of the edge there are in each direction.
                int out_ = 0, in_ = 0;
                var out_edge = edges_[out_edges_[out_]];
                var in_edge = edges_[in_edges_[in_]];
                var sentinel = new Edge(int.MaxValue, int.MaxValue);
                for (; ; )
                {
                    var edge = new Edge[] { out_edge, Reverse(in_edge) }.Min();
                    if (edge == sentinel) break;

                    int out_begin = out_, in_begin = in_;
                    while (out_edge == edge)
                    {
                        out_edge = (++out_ == num_edges) ? sentinel : edges_[out_edges_[out_]];
                    }
                    while (Reverse(in_edge) == edge)
                    {
                        in_edge = (++in_ == num_edges) ? sentinel : edges_[in_edges_[in_]];
                    }
                    var n_out = out_ - out_begin;
                    var n_in = in_ - in_begin;
                    if (edge.ShapeId == edge.EdgeId)
                    {
                        // This is a degenerate edge.
                        Assert.True(n_out == n_in);
                        if (options_.DegenerateEdges_ == DegenerateEdges.DISCARD)
                        {
                            continue;
                        }
                        if (options_.DegenerateEdges_ == DegenerateEdges.DISCARD_EXCESS &&
                            ((out_begin > 0 &&
                              edges_[out_edges_[out_begin - 1]].ShapeId == edge.ShapeId) ||
                             (out_ < num_edges && edges_[out_edges_[out_]].ShapeId == edge.ShapeId) ||
                             (in_begin > 0 &&
                              edges_[in_edges_[in_begin - 1]].EdgeId == edge.ShapeId) ||
                             (in_ < num_edges && edges_[in_edges_[in_]].EdgeId == edge.ShapeId)))
                        {
                            continue;  // There were non-degenerate incident edges, so discard.
                        }
                        // DegenerateEdges::DISCARD_EXCESS also merges degenerate edges.
                        bool merge =
                            (options_.DuplicateEdges_ == DuplicateEdges.MERGE ||
                             options_.DegenerateEdges_ == DegenerateEdges.DISCARD_EXCESS);
                        if (options_.EdgeType_ == EdgeType.UNDIRECTED &&
                            (options_.SiblingPairs_ == SiblingPairs.REQUIRE ||
                             options_.SiblingPairs_ == SiblingPairs.CREATE))
                        {
                            // When we have undirected edges and are guaranteed to have siblings,
                            // we cut the number of edges in half (see s2builder.h).
                            Assert.True(0 == (n_out & 1));  // Number of edges is always even.
                            AddEdges(merge ? 1 : (n_out / 2), edge, MergeInputIds(out_begin, out_));
                        }
                        else if (merge)
                        {
                            AddEdges(options_.EdgeType_ == EdgeType.UNDIRECTED ? 2 : 1,
                                     edge, MergeInputIds(out_begin, out_));
                        }
                        else if (options_.SiblingPairs_ == SiblingPairs.DISCARD ||
                                 options_.SiblingPairs_ == SiblingPairs.DISCARD_EXCESS)
                        {
                            // Any SiblingPair option that discards edges causes the labels of all
                            // duplicate edges to be merged together (see s2builder.h).
                            AddEdges(n_out, edge, MergeInputIds(out_begin, out_));
                        }
                        else
                        {
                            CopyEdges(out_begin, out_);
                        }
                    }
                    else if (options_.SiblingPairs_ == SiblingPairs.KEEP)
                    {
                        if (n_out > 1 && options_.DuplicateEdges_ == DuplicateEdges.MERGE)
                        {
                            AddEdge(edge, MergeInputIds(out_begin, out_));
                        }
                        else
                        {
                            CopyEdges(out_begin, out_);
                        }
                    }
                    else if (options_.SiblingPairs_ == SiblingPairs.DISCARD)
                    {
                        if (options_.EdgeType_ == EdgeType.DIRECTED)
                        {
                            // If n_out == n_in: balanced sibling pairs
                            // If n_out < n_in:  unbalanced siblings, in the form AB, BA, BA
                            // If n_out > n_in:  unbalanced siblings, in the form AB, AB, BA
                            if (n_out <= n_in) continue;
                            // Any option that discards edges causes the labels of all duplicate
                            // edges to be merged together (see s2builder.h).
                            AddEdges(options_.DuplicateEdges_ == DuplicateEdges.MERGE ?
                                     1 : (n_out - n_in), edge, MergeInputIds(out_begin, out_));
                        }
                        else
                        {
                            if ((n_out & 1) == 0) continue;
                            AddEdge(edge, MergeInputIds(out_begin, out_));
                        }
                    }
                    else if (options_.SiblingPairs_ == SiblingPairs.DISCARD_EXCESS)
                    {
                        if (options_.EdgeType_ == EdgeType.DIRECTED)
                        {
                            // See comments above.  The only difference is that if there are
                            // balanced sibling pairs, we want to keep one such pair.
                            if (n_out < n_in) continue;
                            AddEdges(options_.DuplicateEdges_ == DuplicateEdges.MERGE ?
                                     1 : Math.Max(1, n_out - n_in), edge, MergeInputIds(out_begin, out_));
                        }
                        else
                        {
                            AddEdges(((n_out & 1) != 0) ? 1 : 2, edge, MergeInputIds(out_begin, out_));
                        }
                    }
                    else
                    {
                        Assert.True(options_.SiblingPairs_ == SiblingPairs.REQUIRE ||
                               options_.SiblingPairs_ == SiblingPairs.CREATE);
                        if (error.IsOk() && options_.SiblingPairs_ == SiblingPairs.REQUIRE &&
                            (options_.EdgeType_ == EdgeType.DIRECTED ? (n_out != n_in)
                                                                        : ((n_out & 1) != 0)))
                        {
                            error = new(S2ErrorCode.BUILDER_MISSING_EXPECTED_SIBLING_EDGES,
                                        "Expected all input edges to have siblings, but some were missing");
                        }
                        if (options_.DuplicateEdges_ == DuplicateEdges.MERGE)
                        {
                            AddEdge(edge, MergeInputIds(out_begin, out_));
                        }
                        else if (options_.EdgeType_ == EdgeType.UNDIRECTED)
                        {
                            // Convert graph to use directed edges instead (see documentation of
                            // REQUIRE/CREATE for undirected edges).
                            AddEdges((n_out + 1) / 2, edge, MergeInputIds(out_begin, out_));
                        }
                        else
                        {
                            CopyEdges(out_begin, out_);
                            if (n_in > n_out)
                            {
                                // Automatically created edges have no input edge ids or labels.
                                AddEdges(n_in - n_out, edge, IdSetLexicon.kEmptySetId);
                            }
                        }
                    }
                }

                LinqUtils.Swap(edges_, new_edges_);
                LinqUtils.Swap(input_ids_, new_input_ids_);
                edges_.TrimExcess();
                input_ids_.TrimExcess();
            }

            private void AddEdge(Edge edge, int input_edge_id_set_id)
            {
                new_edges_.Add(edge);
                new_input_ids_.Add(input_edge_id_set_id);
            }
            private void AddEdges(int num_edges, Edge edge, int input_edge_id_set_id)
            {
                for (var i = 0; i < num_edges; ++i)
                {
                    AddEdge(edge, input_edge_id_set_id);
                }
            }
            private void CopyEdges(int out_begin, int out_end)
            {
                for (var i = out_begin; i < out_end; ++i)
                {
                    AddEdge(edges_[out_edges_[i]], input_ids_[out_edges_[i]]);
                }
            }
            private int MergeInputIds(int out_begin, int out_end)
            {
                if (out_end - out_begin == 1)
                {
                    return input_ids_[out_edges_[out_begin]];
                }
                tmp_ids_.Clear();
                for (var i = out_begin; i < out_end; ++i)
                {
                    foreach (var id in id_set_lexicon_.IdSet_(input_ids_[out_edges_[i]]))
                    {
                        tmp_ids_.Add(id);
                    }
                }
                return id_set_lexicon_.Add(tmp_ids_);
            }

            private readonly GraphOptions options_;
            private readonly List<Edge> edges_;
            private readonly List<int> input_ids_;
            private readonly IdSetLexicon id_set_lexicon_;
            private readonly EdgeLoop out_edges_;
            private readonly EdgeLoop in_edges_;

            private readonly List<Edge> new_edges_ = new();
            private readonly List<int> new_input_ids_ = new();

            private readonly InputEdgeLoop tmp_ids_ = new();
        }

        private class PolylineBuilder
        {
            public PolylineBuilder(Graph g)
            {
                g_ = g; in_ = new VertexInMap(g); out_ = new VertexOutMap(g);
                min_input_ids_ = g.GetMinInputEdgeIds();
                directed_ = g_.Options.EdgeType_ == EdgeType.DIRECTED;
                edges_left_ = g.NumEdges / (directed_ ? 1 : 2);
                used_ = new bool[g.NumEdges].Fill(false);

                if (!directed_)
                {
                    sibling_map_ = in_.InEdgeIds;
                    g.MakeSiblingMap(sibling_map_);
                }
            }
            public List<EdgePolyline> BuildPaths()
            {
                // First build polylines starting at all the vertices that cannot be in the
                // polyline interior (i.e., indegree != 1 or outdegree != 1 for directed
                // edges, or degree != 2 for undirected edges).  We consider the possible
                // starting edges in input edge id order so that we preserve the input path
                // direction even when undirected edges are used.  (Undirected edges are
                // represented by sibling pairs where only the edge in the input direction
                // is labeled with an input edge id.)
                var polylines = new List<EdgePolyline>();
                var edges = GetInputEdgeOrder(min_input_ids_);
                foreach (var e in edges)
                {
                    if (!used_[e] && !IsInterior(g_.GetEdge(e).ShapeId))
                    {
                        polylines.Add(BuildPath(e));
                    }
                }
                // If there are any edges left, they form non-intersecting loops.  We build
                // each loop and then canonicalize its edge order.  We consider candidate
                // starting edges in input edge id order in order to preserve the input
                // direction of undirected loops.  Even so, we still need to canonicalize
                // the edge order to ensure that when an input edge is split into an edge
                // chain, the loop does not start in the middle of such a chain.
                foreach (var e in edges)
                {
                    if (edges_left_ == 0) break;
                    if (used_[e]) continue;
                    var polyline = BuildPath(e);
                    CanonicalizeLoopOrder(min_input_ids_, polyline);
                    polylines.Add(polyline);
                }
                Assert.True(0 == edges_left_);

                // Sort the polylines to correspond to the input order (if possible).
                CanonicalizeVectorOrder(min_input_ids_, polylines);
                return polylines;
            }
            public List<EdgePolyline> BuildWalks()
            {
                // Note that some of this code is worst-case quadratic in the maximum vertex
                // degree.  This could be fixed with a few extra arrays, but it should not
                // be a problem in practice.

                // First, build polylines from all vertices where outdegree > indegree (or
                // for undirected edges, vertices whose degree is odd).  We consider the
                // possible starting edges in input edge id order, for idempotency in the
                // case where multiple input polylines share vertices or edges.
                var polylines = new List<EdgePolyline>();
                var edges = GetInputEdgeOrder(min_input_ids_);
                foreach (var e in edges)
                {
                    if (used_[e]) continue;
                    var v = g_.GetEdge(e).ShapeId;
                    var excess = ExcessDegree(v);
                    if (excess <= 0) continue;
                    excess -= excess_used_[v];
                    if (directed_ ? (excess <= 0) : (excess % 2 == 0)) continue;
                    ++excess_used_[v];
                    polylines.Add(BuildWalk(v));
                    --excess_used_[g_.GetEdge(polylines.Last().Last()).EdgeId];
                }
                // Now all vertices have outdegree == indegree (or even degree if undirected
                // edges are being used).  Therefore all remaining edges can be assembled
                // into loops.  We first try to expand the existing polylines if possible by
                // adding loops to them.
                if (edges_left_ > 0)
                {
                    foreach (var polyline in polylines)
                    {
                        MaximizeWalk(polyline);
                    }
                }
                // Finally, if there are still unused edges then we build loops.  If the
                // input is a polyline that forms a loop, then for idempotency we need to
                // start from the edge with minimum input edge id.  If the minimal input
                // edge was split into several edges, then we start from the first edge of
                // the chain.
                for (var i = 0; i < edges.Count && edges_left_ > 0; ++i)
                {
                    var e = edges[i];
                    if (used_[e]) continue;

                    // Determine whether the origin of this edge is the start of an edge
                    // chain.  To do this, we test whether (outdegree - indegree == 1) for the
                    // origin, considering only unused edges with the same minimum input edge
                    // id.  (Undirected edges have input edge ids in one direction only.)
                    var v = g_.GetEdge(e).ShapeId;
                    var id = min_input_ids_[e];
                    var excess = 0;
                    for (var j = i; j < edges.Count && min_input_ids_[edges[j]] == id; ++j)
                    {
                        var e2 = edges[j];
                        if (used_[e2]) continue;
                        if (g_.GetEdge(e2).ShapeId == v) ++excess;
                        if (g_.GetEdge(e2).EdgeId == v) --excess;
                    }
                    // It is also acceptable to start a polyline from any degenerate edge.
                    if (excess == 1 || g_.GetEdge(e).EdgeId == v)
                    {
                        var polyline = BuildWalk(v);
                        MaximizeWalk(polyline);
                        polylines.Add(polyline);
                    }
                }
                Assert.True(0 == edges_left_);

                // Sort the polylines to correspond to the input order (if possible).
                CanonicalizeVectorOrder(min_input_ids_, polylines);
                return polylines;
            }

            private bool IsInterior(int v)
            {
                if (directed_)
                {
                    return in_.Degree(v) == 1 && out_.Degree(v) == 1;
                }
                else
                {
                    return out_.Degree(v) == 2;
                }
            }
            private int ExcessDegree(int v) =>
                directed_ ? out_.Degree(v) - in_.Degree(v) : out_.Degree(v) % 2;

            private EdgePolyline BuildPath(int e)
            {
                // We simply follow edges until either we reach a vertex where there is a
                // choice about which way to go (where is_interior(v) is false), or we
                // return to the starting vertex (if the polyline is actually a loop).
                var polyline = new EdgePolyline();
                var start = g_.GetEdge(e).ShapeId;
                for (; ; )
                {
                    polyline.Add(e);
                    Assert.True(!used_[e]);
                    used_[e] = true;
                    if (!directed_) used_[sibling_map_[e]] = true;
                    --edges_left_;
                    var v = g_.GetEdge(e).EdgeId;
                    if (!IsInterior(v) || v == start) break;
                    if (directed_)
                    {
                        Assert.True(1 == out_.Degree(v));
                        e = out_.EdgeIds(v).First();
                    }
                    else
                    {
                        Assert.True(2 == out_.Degree(v));
                        foreach (var e2 in out_.EdgeIds(v)) if (!used_[e2]) e = e2;
                    }
                }
                return polyline;
            }
            private EdgePolyline BuildWalk(int v)
            {
                var polyline = new EdgePolyline();
                for (; ; )
                {
                    // Follow the edge with the smallest input edge id.
                    var best_edge = -1;
                    var best_out_id = int.MaxValue;
                    foreach (var e in out_.EdgeIds(v))
                    {
                        if (used_[e] || min_input_ids_[e] >= best_out_id) continue;
                        best_out_id = min_input_ids_[e];
                        best_edge = e;
                    }
                    if (best_edge < 0) return polyline;
                    // For idempotency when there are multiple input polylines, we stop the
                    // walk early if "best_edge" might be a continuation of a different
                    // incoming edge.
                    var excess = ExcessDegree(v) - excess_used_[v];
                    if (directed_ ? (excess < 0) : (excess % 2) == 1)
                    {
                        foreach (var e in in_.EdgeIds(v))
                        {
                            if (!used_[e] && min_input_ids_[e] <= best_out_id)
                            {
                                return polyline;
                            }
                        }
                    }
                    polyline.Add(best_edge);
                    used_[best_edge] = true;
                    if (!directed_) used_[sibling_map_[best_edge]] = true;
                    --edges_left_;
                    v = g_.GetEdge(best_edge).EdgeId;
                }
            }
            private void MaximizeWalk(EdgePolyline polyline)
            {
                // Examine all vertices of the polyline and check whether there are any
                // unused outgoing edges.  If so, then build a loop starting at that vertex
                // and insert it into the polyline.  (The walk is guaranteed to be a loop
                // because this method is only called when all vertices have equal numbers
                // of unused incoming and outgoing edges.)
                for (var i = 0; i <= polyline.Count; ++i)
                {
                    var v = i == 0 ? g_.GetEdge(polyline[i]).ShapeId
                                  : g_.GetEdge(polyline[i - 1]).EdgeId;
                    foreach (var e in out_.EdgeIds(v))
                    {
                        if (!used_[e])
                        {
                            var loop = BuildWalk(v);
                            Assert.True(v == g_.GetEdge(loop.Last()).EdgeId);
                            polyline.AddRange(loop);
                            Assert.True(used_[e]);  // All outgoing edges from "v" are now used.
                            break;
                        }
                    }
                }
            }

            private readonly Graph g_;
            private readonly VertexInMap in_;
            private readonly VertexOutMap out_;
            private readonly EdgeLoop sibling_map_;
            private readonly InputEdgeLoop min_input_ids_;
            private readonly bool directed_;
            private int edges_left_;
            private readonly bool[] used_;
            // A map of (outdegree(v) - indegree(v)) considering used edges only.
            private readonly Dictionary<int, int> excess_used_ = new(); // gtl.btree_map
        }
    }
}
