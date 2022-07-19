namespace S2Geometry;

using Base = S2ClosestCellQueryBase<S1ChordAngle>;
using S2MinDistanceTarget = S2DistanceTarget<S1ChordAngle>; // S1ChordAngle = S2MinDistance, S2MinDistanceTarget = Target

// S2ClosestCellQuery is a helper class for finding the closest cell(s) to a
// given point, edge, S2Cell, S2CellUnion, or geometry collection.  A typical
// use case would be to add a collection of S2Cell coverings to an S2CellIndex
// (representing a collection of original geometry), and then use
// S2ClosestCellQuery to find all coverings that are within a given distance
// of some target geometry (which could be represented exactly, or could also
// be a covering).  The distance to the original geometry corresponding to
// each covering could then be measured more precisely if desired.
//
// For example, here is how to find all cells that are closer than
// "distance_limit" to a given target point:
//
//   S2ClosestCellQuery query(out cell_index);
//   query.Options().set_max_distance(distance_limit);
//   S2ClosestCellQuery.PointTarget target(target_point);
//   foreach (var& result in query.FindClosestCells(out target)) {
//     // result.distance() is the distance to the target.
//     // result.cell_id() is the indexed S2CellId.
//     // result.label() is the integer label associated with the S2CellId.
//     DoSomething(target_point, result);
//   }
//
// You can find either the k closest cells, or all cells within a given
// radius, or both (i.e., the k closest cells up to a given maximum radius).
// By default *all* cells are returned, so you should always specify either
// max_results() or max_distance() or both.  You can also restrict the results
// to cells that intersect a given S2Region; for example:
//
//   S2LatLngRect rect(...);
//   query.Options().set_region(out rect);  // Does *not* take ownership.
//
// There is a FindClosestCell() convenience method that returns the closest
// cell.  However, if you only need to test whether the distance is above or
// below a given threshold (e.g., 10 km), it is typically much faster to use
// the IsDistanceLess() method instead.  Unlike FindClosestCell(), this method
// stops as soon as it can prove that the minimum distance is either above or
// below the threshold.  Example usage:
//
//   if (query.IsDistanceLess(out target, limit_distance)) ...
//
// To find the closest cells to a query edge rather than a point, use:
//
//   S2ClosestCellQuery.EdgeTarget target(v0, v1);
//   query.FindClosestCells(out target);
//
// Similarly you can find the closest cells to an S2Cell using an
// S2ClosestCellQuery.CellTarget, you can find the closest cells to an
// S2CellUnion using an S2ClosestCellQuery.CellUnionTarget, and you can find
// the closest cells to an arbitrary collection of points, polylines, and
// polygons by using an S2ClosestCellQuery.ShapeIndexTarget.
//
// The implementation is designed to be fast for both simple and complex
// geometric objects.
//
// See S2ClosestCellQueryBase for full documentation.
public class S2ClosestCellQuery
{
    #region MaxBruteForceIndexSize

    // The thresholds for using the brute force algorithm are generally tuned to
    // optimize IsDistanceLess (which compares the distance against a threshold)
    // rather than FindClosest (which actually computes the minimum distance).
    // This is because the former operation is (1) more common, (2) inherently
    // faster, and (3) closely related to finding all cells within a given
    // distance, which is also very common.

    // Break-even points:                   Point cloud      Cap coverings
    // BM_FindClosest                                18                 16
    // BM_IsDistanceLess                              8                  9
    private const int PointTarget_MaxBruteForceIndexSize = 9;

    // Break-even points:                   Point cloud      Cap coverings
    // BM_FindClosestToLongEdge                      14                 16
    // BM_IsDistanceLessToLongEdge                    5                  5
    private const int EdgeTarget_MaxBruteForceIndexSize = 5;

    // Break-even points:                   Point cloud      Cap coverings
    // BM_FindClosestToSmallCell                     12                 13
    // BM_IsDistanceLessToSmallCell                   6                  6
    //
    // Note that the primary use of CellTarget is to implement CellUnionTarget,
    // and therefore it is very important to optimize for the case where a
    // distance limit has been specified.
    private const int CellTarget_MaxBruteForceIndexSize = 6;

    // Break-even points:                   Point cloud      Cap coverings
    // BM_FindClosestToSmallCoarseCellUnion          12                 10
    // BM_IsDistanceLessToSmallCoarseCellUnion        7                  6
    private const int CellUnionTarget_MaxBruteForceIndexSize = 8;

    // Break-even points:                   Point cloud      Cap coverings
    // BM_FindClosestToSmallCoarseShapeIndex         10                  8
    // BM_IsDistanceLessToSmallCoarseShapeIndex       7                  6
    private const int ShapeIndexTarget_MaxBruteForceIndexSize = 7;

    #endregion

    // Options that control the set of cells returned.  Note that by default
    // *all* cells are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    //
    // See S2ClosestCellQueryBase.Options for the full set of options.
    public class Options : Base.Options
    {
        // Like set_max_distance(), except that cells whose distance is exactly
        // equal to "max_distance" are also returned.  Equivalent to calling
        // set_max_distance(max_distance.Successor()).
        public S1ChordAngle InclusiveMaxDistance { set => MaxDistance = (value.Successor()); }

        // Like set_inclusive_max_distance(), except that "max_distance" is also
        // increased by the maximum error in the distance calculation.  This
        // ensures that all cells whose true distance is less than or equal to
        // "max_distance" will be returned (along with some cells whose true
        // distance is slightly greater).
        //
        // Algorithms that need to do exact distance comparisons can use this
        // option to find a set of candidate cells that can then be filtered
        // further (e.g., using S2Pred.CompareDistance).
        public S1ChordAngle ConservativeMaxDistance
        {
            set => MaxDistance = value.PlusError(
                S2.GetUpdateMinDistanceMaxError(value)).Successor();
        }
    }

    // Target subtype that computes the closest distance to a point.
    public sealed class PointTarget : S2MinDistancePointTarget
    {
        public PointTarget(S2Point point) : base(point) { }
        public override int MaxBruteForceIndexSize => PointTarget_MaxBruteForceIndexSize;
    }

    // Target subtype that computes the closest distance to an edge.
    public sealed class EdgeTarget : S2MinDistanceEdgeTarget
    {
        public EdgeTarget(S2Point a, S2Point b) : base(a, b) { }
        public override int MaxBruteForceIndexSize => EdgeTarget_MaxBruteForceIndexSize;
    }

    // Target subtype that computes the closest distance to an S2Cell
    // (including the interior of the cell).
    public sealed class CellTarget : S2MinDistanceCellTarget
    {
        public CellTarget(S2Cell cell) : base(cell) { }
        public override int MaxBruteForceIndexSize => CellTarget_MaxBruteForceIndexSize;
    }

    // Target subtype that computes the closest distance to an S2CellUnion.
    public sealed class CellUnionTarget : S2MinDistanceCellUnionTarget
    {
        public CellUnionTarget(S2CellUnion cell_union) : base(cell_union) { }
        public override int MaxBruteForceIndexSize => CellUnionTarget_MaxBruteForceIndexSize;
    }

    // Target subtype that computes the closest distance to an S2ShapeIndex
    // (an arbitrary collection of points, polylines, and/or polygons).
    //
    // By default, distances are measured to the boundary and interior of
    // polygons in the S2ShapeIndex rather than to polygon boundaries only.
    // If you wish to change this behavior, you may call
    //
    //   target.set_include_interiors(false);
    //
    // (see S2MinDistanceShapeIndexTarget for details).
    public sealed class ShapeIndexTarget : S2MinDistanceShapeIndexTarget
    {
        public ShapeIndexTarget(S2ShapeIndex index) : base(index) { }
        public override int MaxBruteForceIndexSize => ShapeIndexTarget_MaxBruteForceIndexSize;
    }

    // Convenience constructor that calls Init().  Options may be specified here
    // or changed at any time using the Options() accessor method.
    //
    // REQUIRES: "index" must persist for the lifetime of this object.
    // REQUIRES: ReInit() must be called if "index" is modified.
    public S2ClosestCellQuery(S2CellIndex index, Options? options = null)
    {
        Options_ = options ?? new Options();
        base_.Init(index);
    }

    // Default constructor; requires Init() to be called.
    public S2ClosestCellQuery() { }

    // Reinitializes the query.  This method must be called if the underlying
    // S2CellIndex is modified (by calling Clear() and Build() again).
    public void ReInit()
    {
        base_.ReInit();
    }

    // Returns a reference to the underlying S2CellIndex.
    public S2CellIndex Index()
    {
        return base_.Index;
    }

    // Returns the query options.  Options can be modified between queries.
    public Options Options_ { get; private set; }

    // Returns the closest cells to the given target that satisfy the current
    // options.  This method may be called multiple times.
    public List<Base.Result> FindClosestCells(S2MinDistanceTarget target)
    {
        return base_.FindClosestCells(target, Options_);
    }

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    public void FindClosestCells(S2MinDistanceTarget target, List<Base.Result> results)
    {
        base_.FindClosestCells(target, Options_, results);
    }

    //////////////////////// Convenience Methods ////////////////////////

    // Returns the closest cell to the target.  If no cell satisfies the search
    // criteria, then the Result object will have distance == Infinity() and
    // IsEmpty == true.
    public Base.Result FindClosestCell(S2MinDistanceTarget target)
    {
        // System.Diagnostics.Debug.Assert(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
        Options tmp_options = Options_;
        tmp_options.MaxResults = (1);
        return base_.FindClosestCell(target, tmp_options);
    }

    // Returns the minimum distance to the target.  If the index or target is
    // empty, returns S1ChordAngle.Infinity.
    //
    // Use IsDistanceLess() if you only want to compare the distance against a
    // threshold value, since it is often much faster.
    public S1ChordAngle GetDistance(S2MinDistanceTarget target)
    {
        return FindClosestCell(target).Distance;
    }

    // Returns true if the distance to "target" is less than "limit".
    //
    // This method is usually much faster than GetDistance(), since it is much
    // less work to determine whether the minimum distance is above or below a
    // threshold than it is to calculate the actual minimum distance.
    public bool IsDistanceLess(S2MinDistanceTarget target, S1ChordAngle limit)
    {
        // System.Diagnostics.Debug.Assert(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
        Options tmp_options = Options_;
        tmp_options.MaxResults = (1);
        tmp_options.MaxDistance = (limit);
        tmp_options.MaxError = (S1ChordAngle.Straight);
        return !base_.FindClosestCell(target, tmp_options).IsEmpty();
    }

    // Like IsDistanceLess(), but also returns true if the distance to "target"
    // is exactly equal to "limit".
    public bool IsDistanceLessOrEqual(S2MinDistanceTarget target, S1ChordAngle limit)
    {
        // System.Diagnostics.Debug.Assert(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
        Options tmp_options = Options_;
        tmp_options.MaxResults = (1);
        tmp_options.InclusiveMaxDistance = (limit);
        tmp_options.MaxError = (S1ChordAngle.Straight);
        return !base_.FindClosestCell(target, tmp_options).IsEmpty();
    }

    // Like IsDistanceLessOrEqual(), except that "limit" is increased by the
    // maximum error in the distance calculation.  This ensures that this
    // function returns true whenever the true, exact distance is less than
    // or equal to "limit".
    public bool IsConservativeDistanceLessOrEqual(S2MinDistanceTarget target, S1ChordAngle limit)
    {
        // System.Diagnostics.Debug.Assert(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
        Options tmp_options = Options_;
        tmp_options.MaxResults = (1);
        tmp_options.ConservativeMaxDistance = (limit);
        tmp_options.MaxError = (S1ChordAngle.Straight);
        return !base_.FindClosestCell(target, tmp_options).IsEmpty();
    }

    private readonly Base base_ = new();
}
