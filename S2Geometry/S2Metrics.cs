// The following are various constants that describe the shapes and sizes of
// S2Cells (see s2coords.h and s2cell_id.h).  They are useful for deciding
// which cell level to use in order to satisfy a given condition (e.g. that
// cell vertices must be no further than "x" apart).  All of the raw constants
// are differential quantities; you can use the GetValue(level) method to
// compute the corresponding length or area on the unit sphere for cells at a
// given level.  The minimum and maximum bounds are valid for cells at all
// levels, but they may be somewhat conservative for very large cells
// (e.g. face cells).
//
// All of the values below were obtained by a combination of hand analysis and
// Mathematica.  In general, S2_TAN_PROJECTION produces the most uniform
// shapes and sizes of cells, S2_LINEAR_PROJECTION is considerably worse, and
// S2_QUADRATIC_PROJECTION is somewhere in between (but generally closer to
// the tangent projection than the linear one).
//
// Note that S2_LINEAR_PROJECTION can be useful for analysis even when another
// projection is being used, since it allows many cell metrics to be bounded
// in terms of (u,v) coordinates rather than (s,t) coordinates.  (With the
// linear projection, u = 2 * s - 1 and similarly for v.)  Similarly,
// S2_TAN_PROJECTION allows cell metrics to be bounded in terms of (u,v)
// coordinate changes when they are measured as distances on the unit sphere.

namespace S2Geometry;

public static partial class S2
{
    public abstract class Metric(double deriv, int dim)
    {

        // The "deriv" value of a metric is a derivative, and must be multiplied by
        // a length or area in (s,t)-space to get a useful value.
        public double Deriv { get; private set; } = deriv;

        public int Dim { get; private set; } = dim;

        // Return the value of a metric for cells at the given level. The value is
        // either a length or an area on the unit sphere, depending on the
        // particular metric.
        public double GetValue(int level) { return MathUtils.Ldexp(Deriv, -Dim * level); }

        // Return the level at which the metric has approximately the given value.
        // For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the level at which
        // the average cell edge length is approximately 0.1. The return value is
        // always a valid level.
        public int GetClosestLevel(double value)
        {
            return GetLevelForMaxValue((Dim == 1 ? S2.M_SQRT2 : 2) * value);
        }

        // Return the minimum level such that the metric is at most the given value,
        // or S2CellId::kMaxLevel if there is no such level. For example,
        // S2.kMaxDiag.GetLevelForMaxValue(0.1) returns the minimum level such
        // that all cell diagonal lengths are 0.1 or smaller.  The return value
        // is always a valid level.
        public int GetLevelForMaxValue(double value)
        {
            if (value <= 0) return S2.kMaxCellLevel;

            // This code is equivalent to computing a floating-point "level" value and
            // rounding up.  ilogb() returns the exponent corresponding to a fraction in
            // the range [1,2).
            int level = Math.ILogB(value / Deriv);
            level = Math.Max(0, Math.Min(S2.kMaxCellLevel, -(level >> (Dim - 1))));
            MyDebug.Assert(level == S2.kMaxCellLevel || GetValue(level) <= value);
            MyDebug.Assert(level == 0 || GetValue(level - 1) > value);
            return level;
        }

        // Return the maximum level such that the metric is at least the given value,
        // or 0 if there is no such level.  For example,
        // S2::kMinWidth.GetLevelForMinValue(0.1) returns the maximum level such that
        // all cells have a minimum width of 0.1 or larger.  The return value is
        // always a valid level.
        public int GetLevelForMinValue(double value)
        {
            if (value <= 0) return S2.kMaxCellLevel;

            // This code is equivalent to computing a floating-point "level" value and
            // rounding down.
            int level = Math.ILogB(Deriv / value);
            level = Math.Max(0, Math.Min(S2.kMaxCellLevel, level >> (Dim - 1)));
            MyDebug.Assert(level == 0 || GetValue(level) >= value);
            MyDebug.Assert(level == S2.kMaxCellLevel || GetValue(level + 1) < value);
            return level;
        }
    }

    public sealed class LengthMetric(double deriv) : Metric(deriv, 1)
    {
    }

    public sealed class AreaMetric(double deriv) : Metric(deriv, 2)
    {
    }

    // Each cell is bounded by four planes passing through its four edges and
    // the center of the sphere.  These metrics relate to the angle between each
    // pair of opposite bounding planes, or equivalently, between the planes
    // corresponding to two different s-values or two different t-values.  For
    // example, the maximum angle between opposite bounding planes for a cell at
    // level k is kMaxAngleSpan.GetValue(k), and the average angle span for all
    // cells at level k is approximately kAvgAngleSpan.GetValue(k).

    public static readonly LengthMetric kMinAngleSpan = new(
#if S2_LINEAR_PROJECTION
			1.0																			// 1.000
#elif S2_TAN_PROJECTION
			S2Constants.M_PI_2															// 1.570
#elif S2_QUADRATIC_PROJECTION
            4.0 / 3                                                                     // 1.333
#else
			0
#endif
        );

    public static readonly LengthMetric kMaxAngleSpan = new(
#if S2_LINEAR_PROJECTION
			2																			// 2.000
#elif S2_TAN_PROJECTION
			S2Constants.M_PI_2															// 1.570
#elif S2_QUADRATIC_PROJECTION
            1.704897179199218452                                                        // 1.704
#else
			0
#endif
        );

    // This is true for all projections.
    public static readonly LengthMetric kAvgAngleSpan = new(S2.M_PI_2);

    // The width of geometric figure is defined as the distance between two
    // parallel bounding lines in a given direction.  For cells, the minimum
    // width is always attained between two opposite edges, and the maximum
    // width is attained between two opposite vertices.  However, for our
    // purposes we redefine the width of a cell as the perpendicular distance
    // between a pair of opposite edges.  A cell therefore has two widths, one
    // in each direction.  The minimum width according to this definition agrees
    // with the classic geometric one, but the maximum width is different.  (The
    // maximum geometric width corresponds to kMaxDiag defined below.)
    //
    // For a cell at level k, the distance between opposite edges is at least
    // kMinWidth.GetValue(k) and at most kMaxWidth.GetValue(k).  The average
    // width in both directions for all cells at level k is approximately
    // kAvgWidth.GetValue(k).
    //
    // The width is useful for bounding the minimum or maximum distance from a
    // point on one edge of a cell to the closest point on the opposite edge.
    // For example, this is useful when "growing" regions by a fixed distance.
    //
    // Note that because S2Cells are not usually rectangles, the minimum width of
    // a cell is generally smaller than its minimum edge length.  (The interior
    // angles of an S2Cell range from 60 to 120 degrees.)

    public static readonly LengthMetric kMinWidth = new(
#if S2_LINEAR_PROJECTION
			Math.Sqrt(2.0 / 3)															// 0.816
#elif S2_TAN_PROJECTION
			Math.PI / (2 * S2Constants.M_SQRT2)											// 1.111
#elif S2_QUADRATIC_PROJECTION
            2 * S2.M_SQRT2 / 3                                                 // 0.943
#else
			0
#endif
        );

    // This is true for all projections.
    public static readonly LengthMetric kMaxWidth = kMaxAngleSpan;

    public static readonly LengthMetric kAvgWidth = new(
#if S2_LINEAR_PROJECTION
			1.411459345844456965														// 1.411
#elif S2_TAN_PROJECTION
			1.437318638925160885														// 1.437
#elif S2_QUADRATIC_PROJECTION
            1.434523672886099389                                                        // 1.434
#else
			0
#endif
        );

    // The minimum edge length of any cell at level k is at least
    // kMinEdge.GetValue(k), and the maximum is at most kMaxEdge.GetValue(k).
    // The average edge length is approximately kAvgEdge.GetValue(k).
    //
    // The edge length metrics can also be used to bound the minimum, maximum,
    // or average distance from the center of one cell to the center of one of
    // its edge neighbors.  In particular, it can be used to bound the distance
    // between adjacent cell centers along the space-filling Hilbert curve for
    // cells at any given level.

    public static readonly LengthMetric kMinEdge = new(
#if S2_LINEAR_PROJECTION
			2 * S2Constants.M_SQRT2 / 3													// 0.943
#elif S2_TAN_PROJECTION
			Math.PI / (2 * S2Constants.M_SQRT2)											// 1.111
#elif S2_QUADRATIC_PROJECTION
            2 * S2.M_SQRT2 / 3                                                 // 0.943
#else
			0
#endif
        );

    // This is true for all projections.
    public static readonly LengthMetric kMaxEdge = kMaxAngleSpan;

    public static readonly LengthMetric kAvgEdge = new(
#if S2_LINEAR_PROJECTION
			1.440034192955603643														// 0.943
#elif S2_TAN_PROJECTION
			1.461667032546739266														// 1.111
#elif S2_QUADRATIC_PROJECTION
            1.459213746386106062                                                        // 0.943
#else
			0
#endif
        );

    // The minimum diagonal length of any cell at level k is at least
    // kMinDiag.GetValue(k), and the maximum is at most kMaxDiag.GetValue(k).
    // The average diagonal length is approximately kAvgDiag.GetValue(k).
    //
    // The maximum diagonal also happens to be the maximum diameter of any cell,
    // and also the maximum geometric width (see the discussion above).  So for
    // example, the distance from an arbitrary point to the closest cell center
    // at a given level is at most half the maximum diagonal length.

    public static readonly LengthMetric kMinDiag = new(
#if S2_LINEAR_PROJECTION
			2 * S2Constants.M_SQRT2 / 3													// 0.943
#elif S2_TAN_PROJECTION
			Math.PI* S2Constants.M_SQRT2 / 3											// 1.481
#elif S2_QUADRATIC_PROJECTION
            8 * S2.M_SQRT2 / 9                                                 // 1.257
#else
			0
#endif
        );

    public static readonly LengthMetric kMaxDiag = new(
#if S2_LINEAR_PROJECTION
			2 * S2Constants.M_SQRT2														// 2.828
#elif S2_TAN_PROJECTION
			Math.PI* Math.Sqrt(2.0 / 3)													// 2.565
#elif S2_QUADRATIC_PROJECTION
            2.438654594434021032                                                        // 2.439
#else
			0
#endif
        );

    public static readonly LengthMetric kAvgDiag = new(
#if S2_LINEAR_PROJECTION
			2.031817866418812674														// 2.032
#elif S2_TAN_PROJECTION
			2.063623197195635753														// 2.064
#elif S2_QUADRATIC_PROJECTION
            2.060422738998471683                                                        // 2.060
#else
			0;
#endif
        );

    // The minimum area of any cell at level k is at least kMinArea.GetValue(k),
    // and the maximum is at most kMaxArea.GetValue(k).  The average area of all
    // cells at level k is exactly kAvgArea.GetValue(k).

    public static readonly AreaMetric kMinArea = new(
#if S2_LINEAR_PROJECTION
			4 / (3 * Math.Sqrt(3))														// 0.770
#elif S2_TAN_PROJECTION
			(Math.PI * Math.PI) / (4 * S2Constants.M_SQRT2)								// 1.745
#elif S2_QUADRATIC_PROJECTION
            8 * S2.M_SQRT2 / 9                                                 // 1.257
#else
			0
#endif
        );

    public static readonly AreaMetric kMaxArea = new(
#if S2_LINEAR_PROJECTION
			4																			// 4.000
#elif S2_TAN_PROJECTION
			Math.PI* Math.PI / 4														// 2.467
#elif S2_QUADRATIC_PROJECTION
            2.635799256963161491                                                        // 2.636
#else
			0
#endif
        );

    // This is true for all projections.
    public static readonly AreaMetric kAvgArea = new(S2.M_4_PI / 6);     // 2.094

    // This is the maximum edge aspect ratio over all cells at any level, where
    // the edge aspect ratio of a cell is defined as the ratio of its longest
    // edge length to its shortest edge length.

    public static readonly double kMaxEdgeAspect =
#if S2_LINEAR_PROJECTION
			S2Constants.M_SQRT2															// 1.414 
#elif S2_TAN_PROJECTION
			S2Constants.M_SQRT2															// 1.414
#elif S2_QUADRATIC_PROJECTION
            1.442615274452682920                                                        // 1.442
#else
			0																			// 0.000
#endif
        ;

    // This is the maximum diagonal aspect ratio over all cells at any level,
    // where the diagonal aspect ratio of a cell is defined as the ratio of its
    // longest diagonal length to its shortest diagonal length.
    //
    // This is true for all projections.
    public static readonly double kMaxDiagAspect = Math.Sqrt(3);                    // 1.732
}
