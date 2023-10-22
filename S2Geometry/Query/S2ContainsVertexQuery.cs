// This class determines whether a polygon contains one of its vertices given
// the edges incident to that vertex.  The result is +1 if the vertex is
// contained, -1 if it is not contained, and 0 if the incident edges consist
// of matched sibling pairs (in which case the result cannot be determined
// locally).
//
// Point containment is defined according to the "semi-open" boundary model
// (see S2VertexModel), which means that if several polygons tile the region
// around a vertex, then exactly one of those polygons contains that vertex.
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).

namespace S2Geometry;

public class S2ContainsVertexQuery(S2Point target)
{

    // Indicates that the polygon has an edge between "target" and "v" in the
    // given direction (+1 = outgoing, -1 = incoming, 0 = degenerate).
    public void AddEdge(S2Point v, int direction)
    {
        if (!edge_map_.TryAdd(v, direction))
            edge_map_[v] += direction;
    }

    // Returns +1 if the vertex is contained, -1 if it is not contained, and 0
    // if the incident edges consisted of matched sibling pairs.
    public int ContainsSign()
    {
        // Find the unmatched edge that is immediately clockwise from S2::RefDir(P)
        // but not equal to it.  The result is +1 iff this edge is outgoing.
        //
        // A loop with consecutive vertices A,B,C contains vertex B if and only if
        // the fixed vector R = S2::RefDir(B) is contained by the wedge ABC.  The
        // wedge is closed at A and open at C, i.e. the point B is inside the loop
        // if A = R but not if C = R.  This convention is required for compatibility
        // with S2::VertexCrossing.
        var reference_dir = S2.RefDir(target_);
        var best = new KeyValuePair<S2Point, int>(reference_dir, 0);
        foreach (var e in edge_map_)
        {
            MyDebug.Assert(Math.Abs(e.Value) <= 1);
            if (e.Value == 0) continue;  // This is a "matched" edge.
            if (S2Pred.OrderedCCW(reference_dir, best.Key, e.Key, target_))
            {
                best = e;
            }
        }
        return best.Value;
    }

    private readonly S2Point target_ = target;
    private readonly Dictionary<S2Point, int> edge_map_ = []; // absl::btree_map
}
