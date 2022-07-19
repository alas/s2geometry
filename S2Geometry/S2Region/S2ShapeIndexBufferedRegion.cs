// This class provides a way to expand an arbitrary collection of geometry by
// a fixed radius (an operation variously known as "buffering", "offsetting",
// or "Minkowski sum with a disc") in order to compute an S2CellId covering
// (see S2RegionCoverer).  The resulting covering contains all points within
// the given radius of any point in the original geometry.
//
// This class does not actually buffer the geometry; instead it implements the
// S2Region API by computing the distance from candidate S2CellIds to the
// original geometry.  If this distance is below the given radius then the
// S2CellId intersects the buffered geometry.  For example, if the original
// geometry consists of a single S2Point then the buffered geometry is exactly
// equivalent to an S2Cap with the given radius.  (Note that the region is not
// approximated as a polygonal loop.)
//
// Example usage:
//
// S2CellUnion GetBufferedCovering(S2ShapeIndex& index, S1Angle radius) {
//   S2RegionCoverer coverer;
//   coverer.Options().set_max_cells(20);
//   S2CellUnion covering;
//   S2ShapeIndexBufferedRegion region(out index, radius);
//   coverer.GetCovering(region, &covering);
//   return covering;
// }
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).

namespace S2Geometry;

public sealed class S2ShapeIndexBufferedRegion : IS2Region<S2ShapeIndexBufferedRegion>, IEquatable<S2ShapeIndexBufferedRegion>
{
    #region Fields, Constants

    public S1ChordAngle Radius { get; private set; }

    // In order to handle (radius_ == 0) corectly, we need to test whether
    // distances are less than or equal to "radius_".  This is done by testing
    // whether distances are less than radius_.Successor().
    private S1ChordAngle radius_successor_;

    private readonly S2ClosestEdgeQuery query_;  // This class is not thread-safe! 

    #endregion

    #region Constructors

    // Default constructor; requires Init() to be called.
    public S2ShapeIndexBufferedRegion() { }

    // Constructs a region representing all points within the given radius of
    // any point in the given S2ShapeIndex.
    public S2ShapeIndexBufferedRegion(S2ShapeIndex index, S1ChordAngle radius)
    {
        Radius = radius;
        radius_successor_ = radius.Successor();
        query_ = new S2ClosestEdgeQuery(index);
        query_.Options_.IncludeInteriors = true;
    }

    // Convenience constructor that accepts an S1Angle for the radius.
    // REQUIRES: radius >= S1Angle.Zero()
    public S2ShapeIndexBufferedRegion(S2ShapeIndex index, S1Angle radius)
      : this(index, new S1ChordAngle(radius)) { }

    #endregion

    #region S2ShapeIndexBufferedRegion

    public S2ShapeIndex Index() => query_.Index();

    #endregion

    #region ICustomCloneable

    /// <returns>
    /// returns a <b>shallow</b> copy; it does not make a copy of the
    /// underlying S2ShapeIndex.
    /// </returns>
    public object CustomClone() => new S2ShapeIndexBufferedRegion(Index(), Radius);

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        var orig_cap = Index().MakeS2ShapeIndexRegion().GetCapBound();
        return new S2Cap(orig_cap.Center, orig_cap.Radius + Radius);
    }
    public S2LatLngRect GetRectBound()
    {
        var orig_rect = Index().MakeS2ShapeIndexRegion().GetRectBound();
        return orig_rect.ExpandedByDistance(Radius.ToAngle());
    }

    // This method returns a small non-optimal covering that may include
    // duplicate or overlapping cells.  It should not be used directly.
    // Instead, use S2RegionCoverer.GetCovering or GetFastCovering.
    public void GetCellUnionBound(List<S2CellId> cellids)
    {
        // We start with a covering of the original S2ShapeIndex, and then expand it
        // by replacing each cell with a block of 4 cells whose union contains the
        // original cell buffered by the given radius.
        //
        // This increases the number of cells in the covering by a factor of 4 and
        // increases the covered area by a factor of 16, so it is not a very good
        // covering, but it is much better than always returning the 6 face cells.
        var orig_cellids = new List<S2CellId>();
        Index().MakeS2ShapeIndexRegion().GetCellUnionBound(orig_cellids);

        double radians = Radius.ToAngle().Radians;
        int max_level = S2.kMinWidth.GetLevelForMinValue(radians) - 1;
        if (max_level < 0)
        {
            S2Cap.Full.GetCellUnionBound(cellids);
            return;
        }
        cellids.Clear();
        foreach (S2CellId id in orig_cellids)
        {
            if (id.IsFace())
            {
                S2Cap.Full.GetCellUnionBound(cellids);
                return;
            }
            id.AppendVertexNeighbors(Math.Min(max_level, id.Level() - 1), cellids);
        }
    }

    // The implementation is approximate but conservative; it always returns
    // "false" if the cell is not contained by the buffered region, but it may
    // also return false in some cases where "cell" is in fact contained.
    public bool Contains(S2Cell cell)
    {
        // Return true if the buffered region is guaranteed to cover whole globe.
        if (radius_successor_ > S1ChordAngle.Straight) return true;

        // To implement this method perfectly would require computing the directed
        // Hausdorff distance, which is expensive (and not currently implemented).
        // However the following heuristic is almost as good in practice and much
        // cheaper to compute.

        // Return true if the unbuffered region contains this cell.
        if (Index().MakeS2ShapeIndexRegion().Contains(cell)) return true;

        // Otherwise approximate the cell by its bounding cap.
        //
        // NOTE(ericv): It would be slightly more accurate to first find the closest
        // point in the indexed geometry to the cell, and then measure the actual
        // maximum distance from that point to the cell (a poor man's Hausdorff
        // distance).  But based on actual tests this is not worthwhile.
        S2Cap cap = cell.GetCapBound();
        if (Radius < cap.Radius) return false;

        // Return true if the distance to the cell center plus the radius of the
        // cell's bounding cap is less than or equal to "radius_".
        var target = new S2ClosestEdgeQuery.PointTarget(cell.Center());
        return query_.IsDistanceLess(target, radius_successor_ - cap.Radius);
    }

    // Returns true if any buffered shape intersects "cell" (to within a very
    // small error margin).
    public bool MayIntersect(S2Cell cell)
    {
        // Return true if the distance is less than or equal to "radius_".
        var target = new S2ClosestEdgeQuery.CellTarget(cell);
        return query_.IsDistanceLess(target, radius_successor_);
    }

    // Returns true if the given point is contained by the buffered region,
    // i.e. if it is within the given radius of any original shape.
    public bool Contains(S2Point p)
    {
        var target = new S2ClosestEdgeQuery.PointTarget(p);
        // Return true if the distance is less than or equal to "radius_".
        return query_.IsDistanceLess(target, radius_successor_);
    }

    #endregion

    #region IEquatable

    public override bool Equals(object? other)
    {
        return other is S2ShapeIndexBufferedRegion reg && Equals(reg);
    }

    public bool Equals(S2ShapeIndexBufferedRegion? other)
    {
        return Radius.Equals(other.Radius);
    }

    public static bool operator ==(S2ShapeIndexBufferedRegion x, S2ShapeIndexBufferedRegion y)
        => Equals(x, y);

    public static bool operator !=(S2ShapeIndexBufferedRegion x, S2ShapeIndexBufferedRegion y)
        => !Equals(x, y);

    public override int GetHashCode()
    {
        return Radius.GetHashCode();
    }

    #endregion
}
