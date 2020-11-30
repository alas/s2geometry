using System;
using System.Collections.Generic;
namespace S2Geometry
{
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
    public class S2ContainsVertexQuery
    {
        // "target" is the vertex whose containment will be determined.
        public S2ContainsVertexQuery(S2Point target)
        {
            target_ = target;
            edge_map_ = new Dictionary<S2Point, int>();
        }

        // Indicates that the polygon has an edge between "target" and "v" in the
        // given direction (+1 = outgoing, -1 = incoming, 0 = degenerate).
        public void AddEdge(S2Point v, int direction)
        {
            if (edge_map_.ContainsKey(v))
                edge_map_[v] += direction;
            else
                edge_map_[v] = direction;
        }

        // Returns +1 if the vertex is contained, -1 if it is not contained, and 0
        // if the incident edges consisted of matched sibling pairs.
        public int ContainsSign() {
            // Find the unmatched edge that is immediately clockwise from S2PointUtil.Ortho(P).
            var reference_dir = S2PointUtil.Ortho(target_);
            var best = new KeyValuePair<S2Point, int>(reference_dir, 0);
            foreach (var e in edge_map_)
            {
                Assert.True(Math.Abs(e.Value) <= 1);
                if (e.Value == 0) continue;  // This is a "matched" edge.
                if (S2Pred.OrderedCCW(reference_dir, best.Key, e.Key, target_)) {
                    best = e;
                }
            }
            return best.Value;
        }

        private readonly S2Point target_;
        private readonly Dictionary<S2Point, int> edge_map_; // TODO use a btree
    }
}
