// A polygon degeneracy is either a degenerate edge (an edge from a vertex to
// itself) or a sibling edge pair (consisting of an edge and its corresponding
// reverse edge).  "is_hole" indicates whether the degeneracy corresponds to a
// polygon hole (as opposed to a polygon shell).
//
// Degeneracies are not allowed to coincide with any non-degenerate portions
// of the polygon's boundary (since that would make it impossible to classify
// the degeneracy as a shell or hole).  Specifically, degenerate edges must
// coincide only with other degenerate edges, and sibling pairs must coincide
// only with other sibling pairs.  (Below we require a slightly stronger
// condition, namely that sibling pairs cannot coincide with any other edges.)

namespace S2Geometry.S2BuilderUtil;

using System.Diagnostics.CodeAnalysis;
using static S2Builder;
using static S2Builder.GraphOptions;

public readonly record struct PolygonDegeneracy(
                uint EdgeId,
                bool IsHole) : IComparable<PolygonDegeneracy>
{
    public static readonly PolygonDegeneracy Zero = new(0, false);

    #region Constructors

    public PolygonDegeneracy(Int32 edgeId, bool isHole)
        : this((uint)edgeId, isHole) { }

    #endregion

    #region PolygonDegeneracy

    // Given a graph representing a polygon, finds all degenerate edges and
    // sibling pairs and classifies them as being either shells or holes.  The
    // result vector is sorted by edge id.
    //
    // REQUIRES: g.options().edge_type() == DIRECTED
    // REQUIRES: g.options().sibling_pairs() == DISCARD_EXCESS (or DISCARD)
    // REQUIRES: g.options().degenerate_edges() == DISCARD_EXCESS (or DISCARD)
    //
    // Usually callers will want to specify SiblingPairs.DISCARD_EXCESS and
    // DegenerateEdges.DISCARD_EXCESS in order to remove all redundant
    // degeneracies.  DISCARD is also allowed in case you want to keep only one
    // type of degeneracy (i.e., degenerate edges or sibling pairs).
    //
    // If the graph edges cannot be assembled into loops, the result is undefined.
    // (An error may or may not be returned.)
    public static List<PolygonDegeneracy> FindPolygonDegeneracies(Graph g, out S2Error error)
    {
        CheckGraphOptions(g);
        if (g.Options.DegenerateEdges_ == DegenerateEdges.DISCARD &&
            g.Options.SiblingPairs_ == SiblingPairs.DISCARD)
        {
            error = S2Error.OK;
            return new List<PolygonDegeneracy>();  // All degeneracies have already been discarded.
        }
        return new DegeneracyFinder(g).Run(out error);
    }

    // Given a graph representing a polygon, returns true the graph consists
    // entirely of degenerate edges and/or sibling pairs.  Such a graph represents
    // either the empty polygon together with a collection of degenerate shells,
    // or the full polygon together with a collection of degenerate holes.
    //
    // REQUIRES: g.options().edge_type() == DIRECTED
    // REQUIRES: g.options().sibling_pairs() == DISCARD_EXCESS (or DISCARD)
    // REQUIRES: g.options().degenerate_edges() == DISCARD_EXCESS (or DISCARD)
    public static bool IsFullyDegenerate(Graph g)
    {
        CheckGraphOptions(g);
        var edges = g.Edges;
        for (int e = 0; e < g.NumEdges; ++e)
        {
            var edge = edges[e];
            if (edge.ShapeId == edge.EdgeId) continue;
            if (edges.BinarySearch(Graph.Reverse(edge)) < 0)
            {
                return false;
            }
        }
        return true;
    }

    private static void CheckGraphOptions(Graph g)
    {
        System.Diagnostics.Debug.Assert(g.Options.EdgeType_ == EdgeType.DIRECTED);
        System.Diagnostics.Debug.Assert(g.Options.DegenerateEdges_ == DegenerateEdges.DISCARD ||
               g.Options.DegenerateEdges_ == DegenerateEdges.DISCARD_EXCESS);
        System.Diagnostics.Debug.Assert(g.Options.SiblingPairs_ == SiblingPairs.DISCARD ||
               g.Options.SiblingPairs_ == SiblingPairs.DISCARD_EXCESS);
    }

    #endregion

    #region IComparable

    public int CompareTo([AllowNull] PolygonDegeneracy other)
    {
        var c = EdgeId.CompareTo(other.EdgeId);
        if (c != 0) return c;

        return IsHole.CompareTo(other.IsHole);
    }

    public static bool operator <(PolygonDegeneracy x, PolygonDegeneracy y) => x.CompareTo(y) < 0;
    public static bool operator >(PolygonDegeneracy x, PolygonDegeneracy y) => x.CompareTo(y) > 0;
    public static bool operator <=(PolygonDegeneracy x, PolygonDegeneracy y) => x.CompareTo(y) <= 0;
    public static bool operator >=(PolygonDegeneracy x, PolygonDegeneracy y) => x.CompareTo(y) >= 0; 
    
    #endregion
}

// The algorithm builds a set of connected components containing all edges
// that form degeneracies.  The shell/hole status of each degeneracy is
// initially unknown, and is expressed relative to the root vertex: "is_hole"
// means that the degeneracy is a hole if and only if the root vertex turns
// out to be inside the polygon.
public struct Component
{
    // The root vertex from which this component was built.
    public Int32 Root;

    // +1 if "root" inside the polygon, -1 if outside, and 0 if unknown.
    public int RootSign;

    // The degeneracies found in this component.  "is_hole" is expressed
    // relative to the root vertex: the degeneracy is a hole iff the root vertex
    // turns out to be inside the polygon (i.e., root_sign > 0).
    public List<PolygonDegeneracy> Degeneracies;
}

// The actual implementation of FindPolygonDegeneracies.
public class DegeneracyFinder
{
    public DegeneracyFinder(Graph g)
    {
        g_ = g; in_ = new Graph.VertexInMap(g_); out_ = new Graph.VertexOutMap(g_);
    }

    public List<PolygonDegeneracy> Run(out S2Error error)
    {
        error = S2Error.OK;
        // Mark all degenerate edges and sibling pairs in the "is_edge_degeneracy_"
        // vector, and mark any vertices with unbalanced edges in the
        // "is_vertex_unbalanced_" vector.
        int num_degeneracies = ComputeDegeneracies();
        if (num_degeneracies == 0) return new List<PolygonDegeneracy>();

        // If all edges are degenerate, then use IsFullPolygon() to classify the
        // degeneracies (they are necessarily all the same type).
        if (num_degeneracies == g_.NumEdges)
        {
            bool is_hole = g_.IsFullPolygon(out error);
            var result = new List<PolygonDegeneracy>(g_.NumEdges);
            for (int e = 0; e < g_.NumEdges; ++e)
            {
                result[e] = new PolygonDegeneracy(e, is_hole);
            }
            return result;
        }

        // Otherwise repeatedly build components starting from an unvisited
        // degeneracy.  (This avoids building components that don't contain any
        // degeneracies.)  Each component records the "is_hole" status of each
        // degeneracy relative to the root vertex of that component.  If the
        // component contains any non-degenerate portions, then we also determine
        // whether the root vertex is contained by the component (root_sign).
        // In addition we keep track of the number of components that were
        // completely degenerate (to help us decide whether to build an index).
        var components = new List<Component>();
        Int32 known_vertex = -1;
        int known_vertex_sign = 0;
        int num_unknown_signs = 0;
        is_vertex_used_.Capacity = g_.NumVertices;
        for (int e = 0; e < g_.NumEdges; ++e)
        {
            if (is_edge_degeneracy_[e])
            {
                Int32 root = g_.GetEdge(e).ShapeId;
                if (is_vertex_used_[root]) continue;
                Component component = BuildComponent(root);
                if (component.RootSign == 0)
                {
                    ++num_unknown_signs;
                }
                else
                {
                    known_vertex = root;
                    known_vertex_sign = component.RootSign;
                }
                components.Add(component);
            }
        }

        // If some components have an unknown root_sign (i.e., it is unknown whether
        // the root vertex is contained by the polygon or not), we determine the
        // sign of those root vertices by counting crossings starting from a vertex
        // whose sign is known.  Depending on how many components we need to do this
        // for, it may be worthwhile to build an index first.
        if (num_unknown_signs > 0)
        {
            if (known_vertex_sign == 0)
            {
                known_vertex = FindUnbalancedVertex();
                known_vertex_sign = ContainsVertexSign(known_vertex);
            }
            const int kMaxUnindexedSignComputations = 25;  // Tuned using benchmarks.
            if (num_unknown_signs <= kMaxUnindexedSignComputations)
            {
                ComputeUnknownSignsBruteForce(known_vertex, known_vertex_sign, components);
            }
            else
            {
                ComputeUnknownSignsIndexed(known_vertex, known_vertex_sign, components);
            }
        }
        // Finally we convert the "is_hole" status of each degeneracy from a
        // relative value (compared to the component's root vertex) to an absolute
        // one, and sort all the degeneracies by EdgeId.
        return MergeDegeneracies(components);
    }

    // Methods are documented below.
    private int ComputeDegeneracies()
    {
        is_edge_degeneracy_.Capacity = g_.NumEdges;
        is_vertex_unbalanced_.Capacity = g_.NumVertices;
        int num_degeneracies = 0;
        var in_edge_ids = in_.InEdgeIds;
        int n = g_.NumEdges;
        for (int in_ = 0, out_ = 0; out_ < n; ++out_)
        {
            Edge out_edge = g_.GetEdge(out_);
            if (out_edge.ShapeId == out_edge.EdgeId)
            {
                is_edge_degeneracy_[out_] = true;
                ++num_degeneracies;
            }
            else
            {
                while (in_ < n && Graph.Reverse(g_.GetEdge(in_edge_ids[in_])) < out_edge)
                {
                    ++in_;
                }
                if (in_ < n && Graph.Reverse(g_.GetEdge(in_edge_ids[in_])) == out_edge)
                {
                    is_edge_degeneracy_[out_] = true;
                    ++num_degeneracies;
                }
                else
                {
                    // This edge does not have a sibling, which mean that we can determine
                    // whether either vertex is contained by the polygon (using semi-open
                    // boundaries) by examining only the edges incident to that vertex.
                    // We only mark the first vertex since there is no advantage to
                    // finding more than one unbalanced vertex per connected component.
                    is_vertex_unbalanced_[out_edge.ShapeId] = true;
                }
            }
        }
        return num_degeneracies;
    }

    // Build a connected component starting at the given root vertex.  The
    // information returned includes: the root vertex, whether the containment
    // status of the root vertex could be determined using only the edges in this
    // component, and a vector of the edges that belong to degeneracies along with
    // the shell/hole status of each such edge relative to the root vertex.
    private Component BuildComponent(int root)
    {
        var result = new Component { RootSign = 0, Degeneracies = new List<PolygonDegeneracy>() };
        result.Root = root;
        // We keep track of the frontier of unexplored vertices, and whether each
        // vertex is on the same side of the polygon boundary as the root vertex.
        var frontier = new List<(int, bool)> { (root, true) };
        is_vertex_used_[root] = true;
        while (frontier.Any())
        {
            var last = frontier.Last();
            Int32 v0 = last.Item1;
            bool v0_same_inside = last.Item2;  // Same as root vertex?
            frontier.RemoveAt(frontier.Count - 1);
            if (result.RootSign == 0 && is_vertex_unbalanced_[v0])
            {
                int v0_sign = ContainsVertexSign(v0);
                System.Diagnostics.Debug.Assert(v0_sign != 0);
                result.RootSign = v0_same_inside ? v0_sign : -v0_sign;
            }
            foreach (Int32 e in out_.EdgeIds(v0))
            {
                Int32 v1 = g_.GetEdge(e).EdgeId;
                bool same_inside = v0_same_inside ^ CrossingParity(v0, v1, false);
                if (is_edge_degeneracy_[e])
                {
                    result.Degeneracies.Add(new PolygonDegeneracy(e, same_inside));
                }
                if (is_vertex_used_[v1]) continue;
                same_inside ^= CrossingParity(v1, v0, true);
                frontier.Add((v1, same_inside));
                is_vertex_used_[v1] = true;
            }
        }
        return result;
    }

    // Counts the number of times that (v0, v1) crosses the edges incident to v0,
    // and returns the result modulo 2.  This is equivalent to calling
    // S2.VertexCrossing for the edges incident to v0, except that this
    // implementation is more efficient (since it doesn't need to determine which
    // two edge vertices are the same).
    //
    // If "include_same" is false, then the edge (v0, v1) and its sibling (v1, v0)
    // (if any) are excluded from the parity calculation.
    private bool CrossingParity(Int32 v0, Int32 v1, bool include_same)
    {
        int crossings = 0;
        S2Point p0 = g_.Vertex(v0);
        S2Point p1 = g_.Vertex(v1);
        S2Point p0_ref = S2.RefDir(p0);
        foreach (var edge in out_.Edges(v0))
        {
            if (edge.EdgeId == v1)
            {
                if (include_same) ++crossings;
            }
            else if (S2Pred.OrderedCCW(p0_ref, g_.Vertex(edge.EdgeId), p1, p0))
            {
                ++crossings;
            }
        }
        foreach (Int32 e in in_.EdgeIds(v0))
        {
            Edge edge = g_.GetEdge(e);
            if (edge.ShapeId == v1)
            {
                if (include_same) ++crossings;
            }
            else if (S2Pred.OrderedCCW(p0_ref, g_.Vertex(edge.ShapeId), p1, p0))
            {
                ++crossings;
            }
        }
        return (crossings & 1) != 0;
    }
    private Int32 FindUnbalancedVertex()
    {
        for (Int32 v = 0; v < g_.NumVertices; ++v)
        {
            if (is_vertex_unbalanced_[v]) return v;
        }
        throw new ApplicationException("Could not find previously marked unbalanced vertex");
    }
    private int ContainsVertexSign(Int32 v0)
    {
        var query = new S2ContainsVertexQuery(g_.Vertex(v0));
        foreach (var edge in out_.Edges(v0))
        {
            query.AddEdge(g_.Vertex(edge.EdgeId), 1);
        }
        foreach (Int32 e in in_.EdgeIds(v0))
        {
            query.AddEdge(g_.Vertex(g_.GetEdge(e).ShapeId), -1);
        }
        return query.ContainsSign();
    }
    // Determines any unknown signs of component root vertices by counting
    // crossings starting from a vertex whose sign is known.  This version simply
    // tests all edges for crossings.
    private void ComputeUnknownSignsBruteForce(Int32 known_vertex, int known_vertex_sign, List<Component> components)
    {
        var crosser = new S2EdgeCrosser();
        for (var i = 0; i < components.Count; i++)
        {
            var component = components[i];
            if (component.RootSign != 0) continue;
            bool inside = known_vertex_sign > 0;
            crosser.Init(g_.Vertex(known_vertex), g_.Vertex(component.Root));
            for (Int32 e = 0; e < g_.NumEdges; ++e)
            {
                if (is_edge_degeneracy_[e]) continue;
                var edge = g_.GetEdge(e);
                inside ^= crosser.EdgeOrVertexCrossing(g_.Vertex(edge.ShapeId), g_.Vertex(edge.EdgeId));
            }
            components[i] = new Component
            {
                Degeneracies = component.Degeneracies,
                Root = component.Root,
                RootSign = inside ? 1 : -1,
            };
        }
    }
    // Like ComputeUnknownSignsBruteForce, except that this method uses an index
    // to find the set of edges that cross a given edge.
    private void ComputeUnknownSignsIndexed(Int32 known_vertex, int known_vertex_sign, List<Component> components)
    {
        MutableS2ShapeIndex index = new();
        index.Add(new GraphShape(g_));
        var query = new S2CrossingEdgeQuery(index);
        var crossing_edges = new List<Edge>();
        var crosser = new S2EdgeCrosser();
        for (var i = 0; i < components.Count; i++)
        {
            var component = components[i];
            if (component.RootSign != 0) continue;
            var inside = known_vertex_sign > 0;
            crosser.Init(g_.Vertex(known_vertex), g_.Vertex(component.Root));
            query.GetCandidates(g_.Vertex(known_vertex), g_.Vertex(component.Root),
                               index.Shape(0), crossing_edges);
            foreach (var (_, edgeId) in crossing_edges)
            {
                int e = edgeId;
                if (is_edge_degeneracy_[e]) continue;
                inside ^= crosser.EdgeOrVertexCrossing(g_.Vertex(g_.GetEdge(e).ShapeId),
                                                       g_.Vertex(g_.GetEdge(e).EdgeId));
            }
            components[i] = new Component
            {
                Degeneracies = component.Degeneracies,
                Root = component.Root,
                RootSign = inside ? 1 : -1,
            };
        }
    }
    // Merges the degeneracies from all components together, and computes the
    // final "is_hole" status of each edge (since up to this point, the "is_hole"
    // value has been expressed relative to the root vertex of each component).
    private static List<PolygonDegeneracy> MergeDegeneracies(List<Component> components)
    {
        var result = new List<PolygonDegeneracy>();
        foreach (var component in components)
        {
            System.Diagnostics.Debug.Assert(component.RootSign != 0);
            bool invert = component.RootSign < 0;
            foreach (var d in component.Degeneracies)
            {
                result.Add(new PolygonDegeneracy((int)d.EdgeId, d.IsHole ^ invert));
            }
        }
        result.Sort();
        return result;
    }

    private readonly Graph g_;
    private readonly Graph.VertexInMap in_;
    private readonly Graph.VertexOutMap out_;
    private readonly List<bool> is_vertex_used_ = new();        // Has vertex been visited?
    private readonly List<bool> is_edge_degeneracy_ = new();    // Belongs to a degeneracy?
    private readonly List<bool> is_vertex_unbalanced_ = new();  // Has unbalanced sibling pairs?
}
