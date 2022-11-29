// This file defines a collection of classes that are useful for computing
// minimum distances on the sphere.  Their purpose is to allow code to be
// shared among the various query classes that find nearby geometry, such as
// S2ClosestEdgeQuery, S2ClosestPointQuery, and S2ClosestCellQuery.


// S2MinDistance is a thin wrapper around S1ChordAngle that is used by classes
// such as S2ClosestEdgeQuery to compute minimum distances on the sphere (as
// opposed to maximum distances, ellipsoidal distances, etc).
//
// It implements the Distance concept defined by S2DistanceTarget (see
// s2distance_target.h for details).
global using S2MinDistance = S2Geometry.S1ChordAngle;

// S2MinDistanceTarget represents a geometric object to which distances are
// measured.  Specifically, it is used to compute minimum distances on the
// sphere (as opposed to maximum distances, ellipsoidal distances, etc).
//
// Subtypes are defined below for measuring the distance to a point, an edge,
// an S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
global using S2MinDistanceTarget = S2Geometry.S2DistanceTarget<S2Geometry.S1ChordAngle>;

namespace S2Geometry;

// An S2DistanceTarget subtype for computing the minimum distance to a point.
public class S2MinDistancePointTarget : S2MinDistanceTarget
{
    public S2MinDistancePointTarget(S2Point point)
    {
        point_ = point;
    }
    public sealed override S2Cap GetCapBound()
    {
        return new S2Cap(point_, S2MinDistance.Zero);
    }
    public sealed override bool UpdateMinDistance(S2Point p, ref S2MinDistance min_dist)
    {
        return S2MinDistance.UpdateMin(new S2MinDistance(p, point_), ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MinDistance min_dist)
    {
        return S2.UpdateMinDistance(point_, v0, v1, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MinDistance min_dist)
    {
        return S2MinDistance.UpdateMin(cell.Distance(point_), ref min_dist);
    }
    public sealed override bool VisitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor)
    {
        return index.MakeS2ContainsPointQuery().VisitContainingShapes(
            point_, (S2Shape shape) => visitor(shape, point_));
    }

    private readonly S2Point point_;
}

// An S2DistanceTarget subtype for computing the minimum distance to a edge.
public class S2MinDistanceEdgeTarget : S2MinDistanceTarget
{
    public S2MinDistanceEdgeTarget(S2Point a, S2Point b)
    {
        a_ = a; b_ = b;
    }
    public sealed override S2Cap GetCapBound()
    {
        // The following computes a radius equal to half the edge length in an
        // efficient and numerically stable way.
        double d2 = new S2MinDistance(a_, b_).Length2;
        double r2 = (0.5 * d2) / (1 + Math.Sqrt(1 - 0.25 * d2));
        return new S2Cap((a_ + b_).Normalize(), S2MinDistance.FromLength2(r2));
    }
    public sealed override bool UpdateMinDistance(S2Point p, ref S2MinDistance min_dist)
    {
        return S2.UpdateMinDistance(p, a_, b_, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MinDistance min_dist)
    {
        return S2.UpdateEdgePairMinDistance(a_, b_, v0, v1, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MinDistance min_dist)
    {
        return S2MinDistance.UpdateMin(cell.Distance(a_, b_), ref min_dist);
    }
    public sealed override bool VisitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor)
    {
        // We test the center of the edge in order to ensure that edge targets AB
        // and BA yield identical results (which is not guaranteed by the API but
        // users might expect).  Other options would be to test both endpoints, or
        // return different results for AB and BA in some cases.
        var target = new S2MinDistancePointTarget((a_ + b_).Normalize());
        return target.VisitContainingShapes(index, visitor);
    }

    private readonly S2Point a_;
    private readonly S2Point b_;
}

// An S2DistanceTarget subtype for computing the minimum distance to an S2Cell
// (including the interior of the cell).
public class S2MinDistanceCellTarget : S2MinDistanceTarget
{
    public S2MinDistanceCellTarget(S2Cell cell)
    {
        cell_ = cell;
    }
    public sealed override S2Cap GetCapBound()
    {
        return cell_.GetCapBound();
    }
    public sealed override bool UpdateMinDistance(S2Point p, ref S2MinDistance min_dist)
    {
        return S2MinDistance.UpdateMin(cell_.Distance(p), ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MinDistance min_dist)
    {
        return S2MinDistance.UpdateMin(cell_.Distance(v0, v1), ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MinDistance min_dist)
    {
        return S2MinDistance.UpdateMin(cell_.Distance(cell), ref min_dist);
    }
    public sealed override bool VisitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor)
    {
        // The simplest approach is simply to return the polygons that contain the
        // cell center.  Alternatively, if the index cell is smaller than the target
        // cell then we could return all polygons that are present in the
        // S2ShapeIndexCell, but since the index is built conservatively this may
        // include some polygons that don't quite intersect the cell.  So we would
        // either need to recheck for intersection more accurately, or weaken the
        // VisitContainingShapes contract so that it only guarantees approximate
        // intersection, neither of which seems like a good tradeoff.
        var target = new S2MinDistancePointTarget(cell_.Center());
        return target.VisitContainingShapes(index, visitor);
    }

    private readonly S2Cell cell_;
}

// An S2DistanceTarget subtype for computing the minimum distance to an
// S2CellUnion (including the interior of all cells).
public class S2MinDistanceCellUnionTarget : S2MinDistanceTarget
{
    public S2MinDistanceCellUnionTarget(S2CellUnion cell_union)
    {
        cell_union_ = cell_union;
        query_ = new S2ClosestCellQuery(index_);
        foreach (S2CellId cell_id in cell_union_)
        {
            index_.Add(cell_id, 0);
        }
        index_.Build();
    }

    // Specifies that the distances should be computed by examining every cell
    // in the S2CellIndex (for testing and debugging purposes).
    //
    // DEFAULT: false
    public bool UseBruteForce { get => query_.Options_.UseBruteForce; set => query_.Options_.UseBruteForce = value; }

    // Note that set_max_error() should not be called directly by clients; it is
    // used internally by the S2Closest*Query implementations.
    internal override S2MinDistance MaxError { set => query_.Options_.MaxError = (value); }

    public sealed override S2Cap GetCapBound()
    {
        return cell_union_.GetCapBound();
    }
    public sealed override bool UpdateMinDistance(S2Point p, ref S2MinDistance min_dist)
    {
        var target = new S2ClosestCellQuery.PointTarget(p);
        return UpdateMinDistance(target, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MinDistance min_dist)
    {
        var target = new S2ClosestCellQuery.EdgeTarget(v0, v1);
        return UpdateMinDistance(target, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MinDistance min_dist)
    {
        var target = new S2ClosestCellQuery.CellTarget(cell);
        return UpdateMinDistance(target, ref min_dist);
    }
    public sealed override bool VisitContainingShapes(S2ShapeIndex query_index, ShapeVisitor visitor)
    {
        foreach (var cell_id in cell_union_)
        {
            var target = new S2MinDistancePointTarget(cell_id.ToPoint());
            if (!target.VisitContainingShapes(query_index, visitor))
            {
                return false;
            }
        }
        return true;
    }

    private bool UpdateMinDistance(S2MinDistanceTarget target, ref S2MinDistance min_dist)
    {
        query_.Options_.MaxDistance = (min_dist);
        var r = query_.FindClosestCell(target);
        if (r.IsEmpty()) return false;
        min_dist = r.Distance;
        return true;
    }

    private readonly S2CellUnion cell_union_;
    private readonly S2CellIndex index_ = new();
    private readonly S2ClosestCellQuery query_;
}

// An S2DistanceTarget subtype for computing the minimum distance to an
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
// S2ClosestEdgeQuery options.  For example, if include_interiors is true for
// a ShapeIndexTarget but false for the S2ClosestEdgeQuery where the target
// is used, then distances will be measured from the boundary of one
// S2ShapeIndex to the boundary and interior of the other.
//
// Note that when the distance to a ShapeIndexTarget is zero because the
// target intersects the interior of the query index, you can find a point
// that achieves this zero distance by calling the VisitContainingShapes()
// method directly.  For example:
//
//   S2ClosestEdgeQuery.ShapeIndexTarget target(out target_index);
//   target.VisitContainingShapes(
//       query_index, [](S2Shape* containing_shape,
//                       S2Point target_point) {
//         ... do something with "target_point" ...
//         return false;  // Terminate search
//       }));
public class S2MinDistanceShapeIndexTarget : S2MinDistanceTarget
{
    public S2MinDistanceShapeIndexTarget(S2ShapeIndex index)
    {
        index_ = index;
        query_ = new S2ClosestEdgeQuery(index);
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

    // Note that set_max_error() should not be called directly by clients; it is
    // used internally by the S2Closest*Query implementations.
    internal override S2MinDistance MaxError { set => query_.Options_.MaxError = (value); }

    public sealed override S2Cap GetCapBound()
    {
        return index_.MakeS2ShapeIndexRegion().GetCapBound();
    }
    public sealed override bool UpdateMinDistance(S2Point p, ref S2MinDistance min_dist)
    {
        var target = new S2ClosestEdgeQuery.PointTarget(p);
        return UpdateMinDistance(target, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Point v0, S2Point v1, ref S2MinDistance min_dist)
    {
        var target = new S2ClosestEdgeQuery.EdgeTarget(v0, v1);
        return UpdateMinDistance(target, ref min_dist);
    }
    public sealed override bool UpdateMinDistance(S2Cell cell, ref S2MinDistance min_dist)
    {
        var target = new S2ClosestEdgeQuery.CellTarget(cell);
        return UpdateMinDistance(target, ref min_dist);
    }
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
            bool tested_point = false;
            for (int c = 0; c < num_chains; ++c)
            {
                var chain = shape.GetChain(c);
                if (chain.Length == 0) continue;
                tested_point = true;
                var v0 = shape.ChainEdge(c, 0).V0;
                var target = new S2MinDistancePointTarget(v0);
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
                var target = new S2MinDistancePointTarget(ref_.Point);
                if (!target.VisitContainingShapes(query_index, visitor))
                {
                    return false;
                }
            }
        }
        return true;
    }

    private bool UpdateMinDistance(S2MinDistanceTarget target, ref S2MinDistance min_dist)
    {
        query_.Options_.MaxDistance = (min_dist);
        var r = query_.FindClosestEdge(target);
        if (r.IsEmpty()) return false;
        min_dist = r.Distance;
        return true;
    }

    private readonly S2ShapeIndex index_;
    private readonly S2ClosestEdgeQuery query_;
}
