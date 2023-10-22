// S2HausdorffDistanceQuery is a helper class for computing discrete Hausdorff
// distances between two geometries. This can be useful for e.g. computing how
// far apart a highway and a semi-parallel frontage road ever get.
//
// Both geometries are provided as S2ShapeIndex objects. A S2ShapeIndex is a
// collection of points, polylines, and/or polygons. See s2shape_index.h for
// details.
//
// Discrete directed Hausdorff distance from target geometry A to source
// geometry B is defined as the maximum, over all vertices of A, of the closest
// edge distance from the vertex to geometry B.  It is called _discrete_ because
// the maximum is computed over all _vertices_ of A, rather than over other
// positions such as midpoints of the edges. The current implementation computes
// the discrete Hausdorff distance instead of exact Hausdorff distance because
// the latter incurs significantly larger computational expenses, while the
// former is suitable for most practical use cases.
//
// The undirected Hausdorff distance (usually referred to more simply as just
// Hausdorff distance) between geometries A and B is the maximum of the directed
// Hausdorff distances from A to B, and from B to A.
//
// The difference between directed and undirected Hausdorff distances can be
// illustrated by the following example. Let's say we have two polygonal
// geometries, one representing the continental contiguous United States, and
// the other Catalina island off the coast of California. The directed
// Hausdorff distance from the island to the continental US is going to be about
// 60 km - this is how far the furthest point on the island gets from the
// continent. At the same time, the directed Hausdorff distance from the
// continental US to Catalina island is several thousand kilometers - this is
// how far from the island one can get at the North-East corner of the US. The
// undirected Hausdorff distance between these two entities is the maximum of
// the two directed distances, that is a few thousand kilometers.
//
// For example, given two geometries A and B in the form of S2 shape indexes,
// the following code finds the directed Hausdorff distance from A to B and the
// [undirected] Hausdorff distance between A and B:
//
// bool ComputeHausdorffDistances(const S2ShapeIndex* A, const S2ShapeIndex* B,
//                                S1ChordAngle& directed_distance,
//                                S1ChordAngle& undirected_distance) {
//   S2HausdorffDistanceQuery query(S2HausdorffDistanceQuery::Options());
//
//   absl::optional<DirectedResult> directed_result =
//                  query.GetDirectedHausdorffDistance(A, B);
//   if (!directed_result) {
//     return false;
//   }
//   directed_distance = directed_result->distance();
//
//   absl::optional<Result> undirected_result =
//                  query.GetHausdorffDistance(A, B);
//   undirected_distance = undirected_result->distance();
//
//   return true;
// }
//
// For the definition of Hausdorff distance and other details see
// https://en.wikipedia.org/wiki/Hausdorff_distance.
//

namespace S2Geometry;

public class S2HausdorffDistanceQuery
{
    public Options Options_ { get; set; }

    public S2HausdorffDistanceQuery() => Options_ = new();

    public class Options
    {

        // The include_interiors flag (true by default) indicates that the distance
        // should be computed not just to polygon boundaries of the source index,
        // but to polygon interiors as well. For example, if target shape A is fully
        // contained inside the source shape B, and include_interiors is set to
        // true, then the directed Hausdorff distance from A to B is going to be
        // zero.
        public bool IncludeInteriors { get; set; } = true;
    }

    // DirectedResult stores the results of directed Hausdorff distance queries
    public class DirectedResult(S1ChordAngle distance, S2Point targetPoint)
    {
        // Returns the resulting directed Hausdorff distance value.
        public S1ChordAngle Distance { get; } = distance;

        // Returns the point on the target index on which the directed Hausdorff
        // distance is achieved.
        public S2Point TargetPoint { get; } = targetPoint;
    }

    // Result stores the output of [undirected] Hausdorff distance query. It
    // consists of two directed query results, forward and reverse.
    public class Result(DirectedResult target_to_source, DirectedResult source_to_target)
    {
        // Returns the const reference to the result for the target-to-source
        // directed Hausdorff distance call.
        public DirectedResult TargetToSource { get; } = target_to_source;

        // Returns the const reference to the result for the source-to-target
        // directed Hausdorff distance call.
        public DirectedResult SourceToTarget { get; } = source_to_target;

        // Returns the actual Hausdorff distance, which is the maximum of the two
        // directed query results.
        public S1ChordAngle GetDistance()
        {
            return S1ChordAngle.Max(TargetToSource.Distance,
                            SourceToTarget.Distance);
        }
    }

    // Compute directed Hausdorff distance from the target index to the source
    // index.  Returns nullopt iff at least one of the shape indexes is empty.
    //
    // Note that directed Hausdorff distance from geometry A (as target) to
    // geometry B (as source) is not (in general case) equal to that from B (as
    // target) to A (as source).
    public DirectedResult? GetDirectedResult(
          S2ShapeIndex target, S2ShapeIndex source)
    {
        S2ClosestEdgeQuery closest_edge_query=new(source);
        closest_edge_query.Options_.MaxResults = 1;
        closest_edge_query.Options_.IncludeInteriors = Options_.IncludeInteriors;
        S1ChordAngle max_distance = S1ChordAngle.Negative;
        S2Point source_point=new(), target_point=new();

        // This approximation of Haussdorff distance is based on computing closest
        // point distances from the _vertices_ of the target index to _edges_ of the
        // source index.  Hence we iterate over all shapes in the target index, then
        // over all chains in those shapes, then over all edges in those chains, and
        // then over the edges' vertices.
        foreach (var shape in target) {
            for (int chain_id = 0; chain_id < shape.NumChains(); ++chain_id)
            {
                var chain_length = shape.GetChain(chain_id).Length;
                // We include the first vertex (v0) of an edge only if this is the first
                // edge of a polyline. For point shapes (dim == 0) the first vertex is not
                // needed since it coincides with the second vertex (v1). For polygon
                // shapes (dim == 2) the first vertex is not needed since it coincides
                // with the second vertex of the previous edge (or of the last) edge.
                // TODO(b/212844787): Avoid loading vertices twice by using chain vertex
                // iterators.
                bool include_first_vertex = shape.Dimension() == 1;
                for (int offset = 0; offset < chain_length; ++offset)
                {
                    var edge = shape.ChainEdge(chain_id, offset);

                    if (include_first_vertex)
                    {
                        UpdateMaxDistance(edge.V0, closest_edge_query, ref max_distance,
                                          ref target_point, ref source_point);
                        include_first_vertex = false;
                    }
                    UpdateMaxDistance(edge.V1, closest_edge_query, ref max_distance,
                                      ref target_point, ref source_point);
                }
            }
        }

        if (max_distance.IsNegative())
        {
            return null;
        }
        else
        {
            return new DirectedResult(max_distance, target_point);
        }
    }

    // Same as the above method, only returns the actual distance, or
    // S1ChordAngle::Infinity() iff at least one of the shape indexes is empty.
    public S1ChordAngle GetDirectedDistance(S2ShapeIndex target,
                                   S2ShapeIndex source)
    {
        var directed_result = GetDirectedResult(target, source);
        return directed_result is not null ? directed_result.Distance : S1ChordAngle.Infinity;
    }

    // Compute the [undirected] Hausdorff distance between the target index
    // and the source index.  Returns nullopt iff at least one of the shape
    // indexes is empty.
    //
    // Note that the result of this query is symmetrical with respect to target
    // vs. source, i.e. if target and source indices are swapped, the
    // resulting Hausdorff distance remains unchanged.
    public Result? GetResult(S2ShapeIndex target, S2ShapeIndex source)
    {
        var target_to_source = GetDirectedResult(target, source);
        if (target_to_source is not null)
        {
            return new Result(target_to_source, GetDirectedResult(source, target)!);
        }
        else
        {
            return null;
        }
    }

    // Same as the above method, but only returns the maximum of forward and
    // reverse distances, or S1ChordAngle::Infinity() iff at least one of the
    // shape indexes is empty.
    public S1ChordAngle GetDistance(S2ShapeIndex target, S2ShapeIndex source)
    {
        var result = GetResult(target, source);
        return result is not null ? result.GetDistance() : S1ChordAngle.Infinity;
    }

    // This internally used function computes the closest edge distance from point
    // to the source index via closest_edge_query, and, if necessary, updates the
    // max_distance, the target_point and the source_point.
    private static void UpdateMaxDistance(S2Point point,
                    S2ClosestEdgeQuery closest_edge_query,
                    ref S1ChordAngle max_distance, ref S2Point target_point,
                    ref S2Point source_point) {
        // In case we already have a valid result, it can be used as the lower
        // bound estimate for the final Hausdorff distance. Therefore, if the
        // distance between the current target point and the last source point
        // does not exceed this lower bound, we can safely skip this target point, not
        // updating the maximim distance.
        if (!max_distance.IsNegative() &&
            S2Pred.CompareDistance(point, source_point, max_distance) <= 0)
        {
            return;
        }

      // Find the closest edge and the closest point in the source geometry
      // to the target point.
      S2ClosestEdgeQuery.PointTarget target=new(point);
        var closest_edge =
            closest_edge_query.FindClosestEdge(target);
      if (!closest_edge.IsEmpty() && max_distance<closest_edge.Distance) {
        max_distance = closest_edge.Distance;
        target_point = point;
        source_point = closest_edge_query.Project(point, closest_edge);
      }
    }
}
