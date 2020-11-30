using System;

namespace S2Geometry
{
    // S2MaxDistanceTarget represents a geometric object to which maximum distances
    // on the sphere are measured.
    //
    // Subtypes are defined below for measuring the distance to a point, an edge,
    // an S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
    public abstract class S2MaxDistanceTargets : S2DistanceTarget<S2MaxDistance> { }

    // An S2DistanceTarget subtype for computing the maximum distance to a point.
    public class S2MaxDistancePointTarget : S2MaxDistanceTargets
    {
        public S2MaxDistancePointTarget(S2Point point)
        {
            point_ = point;
        }

        // This method returns an S2Cap that bounds the antipode of the target.  (This
        // is the set of points whose S2MaxDistance to the target is
        // S2MaxDistance.Zero().)
        public sealed override S2Cap GetCapBound()
        {
            return new S2Cap(-point_, S1ChordAngle.Zero);
        }
        public sealed override bool UpdateMinDistance(S2Point p, ref S2MaxDistance min_dist)
        {
            return S2MaxDistance.UpdateMin(new S2MaxDistance(new S1ChordAngle(p, point_)), ref min_dist);
        }
        public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MaxDistance min_dist)
        {
            var dist = min_dist.Distance;
            if (S2EdgeDistances.UpdateMaxDistance(point_, v0, v1, ref dist))
            {
                S2MaxDistance.UpdateMin(new S2MaxDistance(dist), ref min_dist);
                return true;
            }
            return false;
        }
        public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MaxDistance min_dist)
        {
            return S2MaxDistance.UpdateMin(new S2MaxDistance(cell.GetMaxDistance(point_)), ref min_dist);
        }
        public sealed override bool VisitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor)
        {
            // For furthest points, we visit the polygons whose interior contains the
            // antipode of the target point.  (These are the polygons whose
            // S2MaxDistance to the target is S2MaxDistance.Zero().)
            return index.MakeS2ContainsPointQuery().VisitContainingShapes(
                -point_, (S2Shape shape) => visitor(shape, point_));
        }

        private readonly S2Point point_;
    }

    // An S2DistanceTarget subtype for computing the maximum distance to an edge.
    public class S2MaxDistanceEdgeTarget : S2MaxDistanceTargets
    {
        public S2MaxDistanceEdgeTarget(S2Point a, S2Point b)
        {
            a_ = a.Normalized;
            b_ = b.Normalized;
        }

        // This method returns an S2Cap that bounds the antipode of the target.  (This
        // is the set of points whose S2MaxDistance to the target is
        // S2MaxDistance.Zero().)
        public sealed override S2Cap GetCapBound()
        {
            // The following computes a radius equal to half the edge length in an
            // efficient and numerically stable way.
            double d2 = new S1ChordAngle(a_, b_).Length2;
            double r2 = (0.5 * d2) / (1 + Math.Sqrt(1 - 0.25 * d2));
            return new S2Cap(-(a_ + b_).Normalized, S1ChordAngle.FromLength2(r2));
        }
        public sealed override bool UpdateMinDistance(S2Point p, ref S2MaxDistance min_dist)
        {
            var dist = min_dist.Distance;
            if (S2EdgeDistances.UpdateMaxDistance(p, a_, b_, ref dist))
            {
                S2MaxDistance.UpdateMin(new S2MaxDistance(dist), ref min_dist);
                return true;
            }
            return false;
        }
        public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MaxDistance min_dist)
        {
            var dist = min_dist.Distance;
            if (S2EdgeDistances.UpdateEdgePairMaxDistance(a_, b_, v0, v1, ref dist))
            {
                S2MaxDistance.UpdateMin(new S2MaxDistance(dist), ref min_dist);
                return true;
            }
            return false;
        }
        public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MaxDistance min_dist)
        {
            return S2MaxDistance.UpdateMin(new S2MaxDistance(cell.GetMaxDistance(a_, b_)), ref min_dist);
        }
        public sealed override bool VisitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor)
        {
            // We only need to test one edge point.  That is because the method *must*
            // visit a polygon if it fully contains the target, and *is allowed* to
            // visit a polygon if it intersects the target.  If the tested vertex is not
            // contained, we know the full edge is not contained; if the tested vertex is
            // contained, then the edge either is fully contained (must be visited) or it
            // intersects (is allowed to be visited).  We visit the center of the edge so
            // that edge AB gives identical results to BA.
            var target = new S2MaxDistancePointTarget((a_ + b_).Normalized);
            return target.VisitContainingShapes(index, visitor);
        }

        private readonly S2Point a_, b_;
    }

    // An S2DistanceTarget subtype for computing the maximum distance to an S2Cell
    // (including the interior of the cell).
    public class S2MaxDistanceCellTarget : S2MaxDistanceTargets
    {
        public S2MaxDistanceCellTarget(S2Cell cell) { cell_ = cell; }

        // This method returns an S2Cap that bounds the antipode of the target.  (This
        // is the set of points whose S2MaxDistance to the target is
        // S2MaxDistance.Zero().)
        public sealed override S2Cap GetCapBound()
        {
            S2Cap cap = cell_.GetCapBound();
            return new S2Cap(-cap.Center, cap.Radius);
        }
        public sealed override bool UpdateMinDistance(S2Point p, ref S2MaxDistance min_dist)
        {
            return S2MaxDistance.UpdateMin(new S2MaxDistance(cell_.GetMaxDistance(p)), ref min_dist);
        }
        public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MaxDistance min_dist)
        {
            return S2MaxDistance.UpdateMin(new S2MaxDistance(cell_.GetMaxDistance(v0, v1)), ref min_dist);
        }
        public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MaxDistance min_dist)
        {
            return S2MaxDistance.UpdateMin(new S2MaxDistance(cell_.GetMaxDistance(cell)), ref min_dist);
        }
        public sealed override bool VisitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor)
        {
            // We only need to check one point here - cell center is simplest.
            // See comment at S2MaxDistanceEdgeTarget.VisitContainingShapes.
            var target = new S2MaxDistancePointTarget(cell_.GetCenter());
            return target.VisitContainingShapes(index, visitor);
        }

        private readonly S2Cell cell_;
    }

    // An S2DistanceTarget subtype for computing the maximum distance to an
    // S2ShapeIndex (a collection of points, polylines, and/or polygons).
    //
    // Note that ShapeIndexTarget has its own options:
    //
    //   include_interiors()
    //     - specifies that distances are measured to the boundary and interior
    //       of polygons in the S2ShapeIndex.  (If set to false, distance is
    //       measured to the polygon boundary only.)
    //       DEFAULT: true.
    //
    //   brute_force()
    //     - specifies that the distances should be computed by examining every
    //       edge in the S2ShapeIndex (for testing and debugging purposes).
    //       DEFAULT: false.
    //
    // These options are specified independently of the corresponding
    // S2FurthestEdgeQuery options.  For example, if include_interiors is true for
    // a ShapeIndexTarget but false for the S2FurthestEdgeQuery where the target
    // is used, then distances will be measured from the boundary of one
    // S2ShapeIndex to the boundary and interior of the other.
    //
    public class S2MaxDistanceShapeIndexTarget : S2MaxDistanceTargets
    {
        public S2MaxDistanceShapeIndexTarget(S2ShapeIndex index)
        {
            index_ = index;
            query_ = new S2FurthestEdgeQuery(index);
        }

        // Specifies that distance will be measured to the boundary and interior
        // of polygons in the S2ShapeIndex rather than to polygon boundaries only.
        //
        // DEFAULT: true
        public bool IncludeInteriors { get => query_.Options_.IncludeInteriors; set => query_.Options_.IncludeInteriors = value; }

        // Specifies that the distances should be computed by examining every edge
        // in the S2ShapeIndex (for testing and debugging purposes).
        //
        // DEFAULT: false
        public bool UseBruteForce { get => query_.Options_.UseBruteForce; set => query_.Options_.UseBruteForce = value; }

        internal override S1ChordAngle MaxError { set => query_.Options_.MaxError = (value); }

        // This method returns an S2Cap that bounds the antipode of the target.  (This
        // is the set of points whose S2MaxDistance to the target is
        // S2MaxDistance.Zero().)
        public sealed override S2Cap GetCapBound()
        {
            S2Cap cap = index_.MakeS2ShapeIndexRegion().GetCapBound();
            return new S2Cap(-cap.Center, cap.Radius);
        }
        public sealed override bool UpdateMinDistance(S2Point p, ref S2MaxDistance min_dist)
        {
            query_.Options_.MinDistance = (min_dist.Distance);
            var target = new S2FurthestEdgeQuery.PointTarget(p);
            var r = query_.FindFurthestEdge(target);
            if (r.ShapeId() < 0)
            {
                return false;
            }
            min_dist = new S2MaxDistance(r.Distance());
            return true;
        }
        public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MaxDistance min_dist)
        {
            query_.Options_.MinDistance = (min_dist.Distance);
            var target = new S2FurthestEdgeQuery.EdgeTarget(v0, v1);
            var r = query_.FindFurthestEdge(target);
            if (r.ShapeId() < 0) return false;
            min_dist = new S2MaxDistance(r.Distance());
            return true;
        }
        public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MaxDistance min_dist)
        {
            query_.Options_.MinDistance = (min_dist.Distance);
            var target = new S2FurthestEdgeQuery.CellTarget(cell);
            var r = query_.FindFurthestEdge(target);
            if (r.ShapeId() < 0) return false;
            min_dist = new S2MaxDistance(r.Distance());
            return true;
        }

        // For target types consisting of multiple connected components (such as
        // S2MaxDistanceShapeIndexTarget), this method should return the
        // polygons containing the antipodal reflection of *any* connected
        // component.  (It is sufficient to test containment of one vertex per
        // connected component, since the API allows us to also return any polygon
        // whose boundary has S2MaxDistance.Zero() to the target.)
        public sealed override bool VisitContainingShapes(S2ShapeIndex query_index, ShapeVisitor visitor)
        {
            // It is sufficient to find the set of chain starts in the target index
            // (i.e., one vertex per connected component of edges) that are contained by
            // the query index, except for one special case to handle full polygons.
            //
            // TODO(ericv): Do this by merge-joining the two S2ShapeIndexes, and share
            // the code with S2BooleanOperation.
            foreach (var shape in index_)
            {
                if (shape == null) continue;
                int num_chains = shape.NumChains();
                // Shapes that don't have any edges require a special case (below).
                var tested_point = false;
                for (int c = 0; c < num_chains; ++c)
                {
                    S2Shape.Chain chain = shape.GetChain(c);
                    if (chain.Length == 0) continue;
                    tested_point = true;
                    var target = new S2MaxDistancePointTarget(shape.ChainEdge(c, 0).V0);
                    if (!target.VisitContainingShapes(query_index, visitor))
                    {
                        return false;
                    }
                }
                if (!tested_point)
                {
                    // Special case to handle full polygons.
                    var ref_ = shape.GetReferencePoint();
                    if (!ref_.Contained) continue;
                    var target = new S2MaxDistancePointTarget(ref_.Point);
                    if (!target.VisitContainingShapes(query_index, visitor))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private readonly S2ShapeIndex index_;
        private readonly S2FurthestEdgeQuery query_;
    }
}
