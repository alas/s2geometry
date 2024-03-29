// A SnapFunction that snaps every vertex to itself.  It should be used when
// vertices do not need to be snapped to a discrete set of locations (such as
// E7 lat/lngs), or when maximum accuracy is desired.
//
// If the given "snap_radius" is zero, then all input vertices are preserved
// exactly.  Otherwise, S2Builder merges nearby vertices to ensure that no
// vertex pair is closer than "snap_radius".  Furthermore, vertices are
// separated from non-incident edges by at least "min_edge_vertex_separation",
// equal to (0.5 * snap_radius).  For example, if the snap_radius is 1km, then
// vertices will be separated from non-incident edges by at least 500m.

namespace S2Geometry.S2BuilderUtil;

public class IdentitySnapFunction : SnapFunction
{
    // The default constructor uses a snap_radius of zero (i.e., no snapping).
    public IdentitySnapFunction() => snap_radius_ = S1Angle.Zero;

    // Convenience constructor that calls set_snap_radius().
    public IdentitySnapFunction(S1Angle snap_radius) => SnapRadius = snap_radius;

    // REQUIRES: snap_radius <= SnapFunction.kMaxSnapRadius()
    public override S1Angle SnapRadius
    {
        get => snap_radius_;
        set
        {
            MyDebug.Assert(value <= kMaxSnapRadius);
            snap_radius_ = value;
        }
    }
    private S1Angle snap_radius_;

    // For the identity snap function, all vertex pairs are separated by at
    // least snap_radius().
    //
    // Since SnapFunction does not move the input point, output vertices are
    // separated by the full snap_radius().
    public override S1Angle MinVertexSeparation() => snap_radius_;

    // For the identity snap function, edges are separated from all non-incident
    // vertices by at least 0.5 * snap_radius().
    //
    // In the worst case configuration, the edge-vertex separation is half of the
    // vertex separation.
    public override S1Angle MinEdgeVertexSeparation() => 0.5 * snap_radius_;

    public override S2Point SnapPoint(S2Point point) => point;

    public override object CustomClone() => new IdentitySnapFunction(snap_radius_);
}

// A SnapFunction that snaps vertices to S2CellId centers.  This can be useful
// if you want to encode your geometry compactly using S2Polygon.Encode(),
// for example.  You can snap to the centers of cells at any level.
//
// Every snap level has a corresponding minimum snap radius, which is simply
// the maximum distance that a vertex can move when snapped.  It is
// approximately equal to half of the maximum diagonal length for cells at the
// chosen level.  You can also set the snap radius to a larger value; for
// example, you could snap to the centers of leaf cells (1cm resolution) but
// set the snap_radius() to 10m.  This would result in significant extra
// simplification, without moving vertices unnecessarily (i.e., vertices that
// are at least 10m away from all other vertices will move by less than 1cm).
public class S2CellIdSnapFunction : SnapFunction
{
    // The default constructor snaps to S2Constants.kMaxCellLevel (i.e., the centers
    // of leaf cells), and uses the minimum allowable snap radius at that level.
    public S2CellIdSnapFunction()
    {
        Level = S2.kMaxCellLevel;
    }

    // Convenience constructor equivalent to calling set_level(level).
    public S2CellIdSnapFunction(int level)
    {
        Level = level;
    }

    private S2CellIdSnapFunction(int level, S1Angle snap_radius)
    {
        Level = level;
        SnapRadius = snap_radius;
    }

    // Snaps vertices to S2Cell centers at the given level.  As a side effect,
    // this method also resets "snap_radius" to the minimum value allowed at
    // this level:
    //
    //   set_snap_radius(MinSnapRadiusForLevel(level))
    //
    // This means that if you want to use a larger snap radius than the minimum,
    // you must call set_snap_radius() *after* calling set_level().
    public int Level
    {
        get => level_;
        set
        {
            MyDebug.Assert(value >= 0);
            MyDebug.Assert(value <= S2.kMaxCellLevel);
            level_ = value;
            SnapRadius = MinSnapRadiusForLevel(value);
        }
    }
    private int level_;

    // Defines the snap radius to be used (see s2builder.h).  The snap radius
    // must be at least the minimum value for the current level(), but larger
    // values can also be used (e.g., to simplify the geometry).
    //
    // REQUIRES: snap_radius >= MinSnapRadiusForLevel(level())
    // REQUIRES: snap_radius <= SnapFunction.kMaxSnapRadius()
    public override S1Angle SnapRadius
    {
        get => snap_radius_;
        set
        {
            MyDebug.Assert(value >= MinSnapRadiusForLevel(Level));
            MyDebug.Assert(value <= kMaxSnapRadius);
            snap_radius_ = value;
        }
    }
    private S1Angle snap_radius_;

    // Returns the minimum allowable snap radius for the given S2Cell level
    // (approximately equal to half of the maximum cell diagonal length).
    public static S1Angle MinSnapRadiusForLevel(int level)
    {
        // snap_radius() needs to be an upper bound on the true distance that a
        // point can move when snapped, taking into account numerical errors.
        //
        // The maximum error when converting from an S2Point to an S2CellId is
        // S2.kMaxDiag.deriv() * S2Constants.DoubleEpsilon.  The maximum error when converting an
        // S2CellId center back to an S2Point is 1.5 * S2Constants.DoubleEpsilon.  These add up to
        // just slightly less than 4 * S2Constants.DoubleEpsilon.
        return S1Angle.FromRadians(0.5 * S2.kMaxDiag.GetValue(level) + 4 * S2.DoubleEpsilon);
    }

    // Returns the minimum S2Cell level (i.e., largest S2Cells) such that
    // vertices will not move by more than "snap_radius".  This can be useful
    // when choosing an appropriate level to snap to.  The return value is
    // always a valid level (out of range values are silently clamped).
    //
    // If you want to choose the snap level based on a distance, and then use
    // the minimum possible snap radius for the chosen level, do this:
    //
    //   S2CellIdSnapFunction f(
    //       S2CellIdSnapFunction.LevelForMaxSnapRadius(distance));
    public static int LevelForMaxSnapRadius(S1Angle snap_radius)
    {
        // When choosing a level, we need to acount for the error bound of
        // 4 * S2Constants.DoubleEpsilon that is added by MinSnapRadiusForLevel().
        return S2.kMaxDiag.GetLevelForMaxValue(
            2 * (snap_radius.Radians - 4 * S2.DoubleEpsilon));
    }

    // For S2CellId snapping, the minimum separation between vertices depends on
    // level() and snap_radius().  It can vary between 0.5 * snap_radius()
    // and snap_radius().
    public override S1Angle MinVertexSeparation()
    {
        // We have three different bounds for the minimum vertex separation: one is
        // a constant bound, one is proportional to snap_radius, and one is equal to
        // snap_radius minus a constant.  These bounds give the best results for
        // small, medium, and large snap radii respectively.  We return the maximum
        // of the three bounds.
        //
        // 1. Constant bound: Vertices are always separated by at least
        //    kMinEdge(level), the minimum edge length for the chosen snap level.
        //
        // 2. Proportional bound: It can be shown that in the plane, the worst-case
        //    configuration has a vertex separation of 2 / Math.Sqrt(13) * snap_radius.
        //    This is verified in the unit test, except that on the sphere the ratio
        //    is slightly smaller at cell level 2 (0.54849 vs. 0.55470).  We reduce
        //    that value a bit more below to be conservative.
        //
        // 3. Best asymptotic bound: This bound bound is derived by observing we
        //    only select a new site when it is at least snap_radius() away from all
        //    existing sites, and the site can move by at most 0.5 * kMaxDiag(level)
        //    when snapped.
        S1Angle min_edge = S1Angle.FromRadians(S2.kMinEdge.GetValue(level_));
        S1Angle max_diag = S1Angle.FromRadians(S2.kMaxDiag.GetValue(level_));
        return S1Angle.Max(min_edge,
                   S1Angle.Max(0.548 * snap_radius_,  // 2 / Math.Sqrt(13) in the plane
                       snap_radius_ - 0.5 * max_diag));
    }

    // For S2CellId snapping, the minimum separation between edges and
    // non-incident vertices depends on level() and snap_radius().  It can
    // be as low as 0.219 * snap_radius(), but is typically 0.5 * snap_radius()
    // or more.
    public override S1Angle MinEdgeVertexSeparation()
    {
        // Similar to min_vertex_separation(), in this case we have four bounds: a
        // constant bound that holds only at the minimum snap radius, a constant
        // bound that holds for any snap radius, a bound that is proportional to
        // snap_radius, and a bound that approaches 0.5 * snap_radius asymptotically.
        //
        // 1. Constant bounds:
        //
        //    (a) At the minimum snap radius for a given level, it can be shown that
        //    vertices are separated from edges by at least 0.5 * kMinDiag(level) in
        //    the plane.  The unit test verifies this, except that on the sphere the
        //    worst case is slightly better: 0.5652980068 * kMinDiag(level).
        //
        //    (b) Otherwise, for arbitrary snap radii the worst-case configuration
        //    in the plane has an edge-vertex separation of Math.Sqrt(3/19) *
        //    kMinDiag(level), where Math.Sqrt(3/19) is about 0.3973597071.  The unit
        //    test verifies that the bound is slighty better on the sphere:
        //    0.3973595687 * kMinDiag(level).
        //
        // 2. Proportional bound: In the plane, the worst-case configuration has an
        //    edge-vertex separation of 2 * Math.Sqrt(3/247) * snap_radius, which is
        //    about 0.2204155075.  The unit test verifies this, except that on the
        //    sphere the bound is slightly worse for certain large S2Cells: the
        //    minimum ratio occurs at cell level 6, and is about 0.2196666953.
        //
        // 3. Best asymptotic bound: If snap_radius() is large compared to the
        //    minimum snap radius, then the best bound is achieved by 3 sites on a
        //    circular arc of radius "snap_radius", spaced "min_vertex_separation"
        //    apart.  An input edge passing just to one side of the center of the
        //    circle intersects the Voronoi regions of the two end sites but not the
        //    Voronoi region of the center site, and gives an edge separation of
        //    (min_vertex_separation ** 2) / (2 * snap_radius).  This bound
        //    approaches 0.5 * snap_radius for large snap radii, i.e.  the minimum
        //    edge-vertex separation approaches half of the minimum vertex
        //    separation as the snap radius becomes large compared to the cell size.

        S1Angle min_diag = S1Angle.FromRadians(S2.kMinDiag.GetValue(level_));
        if (SnapRadius == MinSnapRadiusForLevel(level_))
        {
            // This bound only holds when the minimum snap radius is being used.
            return 0.565 * min_diag;            // 0.500 in the plane
        }
        // Otherwise, these bounds hold for any snap_radius().
        S1Angle vertex_sep = MinVertexSeparation();
        return S1Angle.Max(0.397 * min_diag,          // Math.Sqrt(3 / 19) in the plane
                   S1Angle.Max(0.219 * snap_radius_,  // 2 * Math.Sqrt(3 / 247) in the plane
                       0.5 * (vertex_sep / snap_radius_) * vertex_sep));
    }

    public override S2Point SnapPoint(S2Point point) => new S2CellId(point).Parent(level_).ToPoint();

    public override object CustomClone() => new S2CellIdSnapFunction(level_, snap_radius_);
}

// A SnapFunction that snaps vertices to S2LatLng E5, E6, or E7 coordinates.
// These coordinates are expressed in degrees multiplied by a power of 10 and
// then rounded to the nearest integer.  For example, in E6 coordinates the
// point (23.12345651, -45.65432149) would become (23123457, -45654321).
//
// The main argument of the SnapFunction is the exponent for the power of 10
// that coordinates should be multipled by before rounding.  For example,
// IntLatLngSnapFunction(7) is a function that snaps to E7 coordinates.  The
// exponent can range from 0 to 10.
//
// Each exponent has a corresponding minimum snap radius, which is simply the
// maximum distance that a vertex can move when snapped.  It is approximately
// equal to 1/Math.Sqrt(2) times the nominal point spacing; for example, for
// snapping to E7 the minimum snap radius is (1e-7 / Math.Sqrt(2)) degrees.
// You can also set the snap radius to any value larger than this; this can
// result in significant extra simplification (similar to using a larger
// exponent) but does not move vertices unnecessarily.
public class IntLatLngSnapFunction : SnapFunction
{
    // The minum exponent supported for snapping.
    public const int kMinExponent = 0;

    // The maximum exponent supported for snapping.
    public const int kMaxExponent = 10;

    // The default constructor yields an invalid snap function.  You must set
    // the exponent explicitly before using it.
    public IntLatLngSnapFunction()
    {
        Exponent = -1; snap_radius_ = S1Angle.Zero; from_degrees_ = 0; to_degrees_ = 0;
    }

    // Convenience constructor equivalent to calling set_exponent(exponent).
    public IntLatLngSnapFunction(int exponent)
    {
        Exponent = exponent;
    }

    private IntLatLngSnapFunction(IntLatLngSnapFunction sf)
    {
        Exponent = sf.Exponent; snap_radius_ = sf.snap_radius_;
        from_degrees_ = sf.from_degrees_; to_degrees_ = sf.to_degrees_;
    }

    // Snaps vertices to points whose (lat, lng) coordinates are integers after
    // converting to degrees and multiplying by 10 raised to the given exponent.
    // For example, (exponent == 7) yields E7 coordinates.  As a side effect,
    // this method also resets "snap_radius" to the minimum value allowed for
    // this exponent:
    //
    //   set_snap_radius(MinSnapRadiusForExponent(exponent))
    //
    // This means that if you want to use a larger snap radius than the minimum,
    // you must call set_snap_radius() *after* calling set_exponent().
    //
    // REQUIRES: kMinExponent <= exponent <= kMaxExponent
    public int Exponent
    {
        get => _exponent;
        set
        {
            MyDebug.Assert(value >= kMinExponent);
            MyDebug.Assert(value <= kMaxExponent);
            _exponent = value;
            SnapRadius = MinSnapRadiusForExponent(value);

            // Precompute the scale factors needed for snapping.  Note that these
            // calculations need to exactly match the ones in s1angle.h to ensure
            // that the same S2Points are generated.
            double power = 1;
            for (int i = 0; i < value; ++i) power *= 10;
            from_degrees_ = power;
            to_degrees_ = 1 / power;
        }
    }
    private int _exponent;

    // Defines the snap radius to be used (see s2builder.h).  The snap radius
    // must be at least the minimum value for the current exponent(), but larger
    // values can also be used (e.g., to simplify the geometry).
    //
    // REQUIRES: snap_radius >= MinSnapRadiusForExponent(exponent())
    // REQUIRES: snap_radius <= SnapFunction.kMaxSnapRadius()
    public override S1Angle SnapRadius
    {
        get => snap_radius_; set
        {
            MyDebug.Assert(value >= MinSnapRadiusForExponent(Exponent));
            MyDebug.Assert(value <= kMaxSnapRadius);
            snap_radius_ = value;
        }
    }
    private S1Angle snap_radius_;

    // Returns the minimum allowable snap radius for the given exponent
    // (approximately equal to (Math.Pow(10, -exponent) / Math.Sqrt(2)) degrees).
    public static S1Angle MinSnapRadiusForExponent(int exponent)
    {
        // snap_radius() needs to be an upper bound on the true distance that a
        // point can move when snapped, taking into account numerical errors.
        //
        // The maximum errors in latitude and longitude can be bounded as
        // follows (as absolute errors in terms of S2Constants.DoubleEpsilon):
        //
        //                                      Latitude      Longitude
        // Convert to S2LatLng:                    1.000          1.000
        // Convert to degrees:                     1.032          2.063
        // Scale by 10**exp:                       0.786          1.571
        // Round to integer: 0.5 * S1Angle.Degrees(to_degrees_)
        // Scale by 10**(-exp):                    1.375          2.749
        // Convert to radians:                     1.252          1.503
        // ------------------------------------------------------------
        // Total (except for rounding)             5.445          8.886
        //
        // The maximum error when converting the S2LatLng back to an S2Point is
        //
        //   Math.Sqrt(2) * (maximum error in latitude or longitude) + 1.5 * S2Constants.DoubleEpsilon
        //
        // which works out to (9 * Math.Sqrt(2) + 1.5) * S2Constants.DoubleEpsilon radians.  Finally
        // we need to consider the effect of rounding to integer coordinates
        // (much larger than the errors above), which can change the position by
        // up to (Math.Sqrt(2) * 0.5 * to_degrees_) radians.
        double power = 1;
        for (int i = 0; i < exponent; ++i) power *= 10;
        return S1Angle.FromDegrees(S2.M_SQRT1_2 / power) +
                S1Angle.FromRadians((9 * S2.M_SQRT2 + 1.5) * S2.DoubleEpsilon);
    }

    // Returns the minimum exponent such that vertices will not move by more
    // than "snap_radius".  This can be useful when choosing an appropriate
    // exponent for snapping.  The return value is always a valid exponent
    // (out of range values are silently clamped).
    //
    // If you want to choose the exponent based on a distance, and then use
    // the minimum possible snap radius for that exponent, do this:
    //
    //   IntLatLngSnapFunction f(
    //       IntLatLngSnapFunction.ExponentForMaxSnapRadius(distance));
    public static int ExponentForMaxSnapRadius(S1Angle snap_radius)
    {
        // When choosing an exponent, we need to acount for the error bound of
        // (9 * Math.Sqrt(2) + 1.5) * S2Constants.DoubleEpsilon added by MinSnapRadiusForExponent().
        snap_radius -= S1Angle.FromRadians((9 * S2.M_SQRT2 + 1.5) * S2.DoubleEpsilon);
        snap_radius = S1Angle.Max(snap_radius, S1Angle.FromRadians(1e-30));
        double exponent = Math.Log10(S2.M_SQRT1_2 / snap_radius.GetDegrees());

        // There can be small errors in the calculation above, so to ensure that
        // this function is the inverse of MinSnapRadiusForExponent() we subtract a
        // small error tolerance.
        return Math.Max(kMinExponent,
                   Math.Min(kMaxExponent,
                       (int)Math.Ceiling(exponent - 2 * S2.DoubleEpsilon)));
    }

    // For IntLatLng snapping, the minimum separation between vertices depends on
    // exponent() and snap_radius().  It can vary between snap_radius()
    // and snap_radius().
    public override S1Angle MinVertexSeparation()
    {
        // We have two bounds for the minimum vertex separation: one is proportional
        // to snap_radius, and one is equal to snap_radius minus a constant.  These
        // bounds give the best results for small and large snap radii respectively.
        // We return the maximum of the two bounds.
        //
        // 1. Proportional bound: It can be shown that in the plane, the worst-case
        //    configuration has a vertex separation of (Math.Sqrt(2) / 3) * snap_radius.
        //    This is verified in the unit test, except that on the sphere the ratio
        //    is slightly smaller (0.471337 vs. 0.471404).  We reduce that value a
        //    bit more below to be conservative.
        //
        // 2. Best asymptotic bound: This bound bound is derived by observing we
        //    only select a new site when it is at least snap_radius() away from all
        //    existing sites, and snapping a vertex can move it by up to
        //    ((1 / Math.Sqrt(2)) * to_degrees_) degrees.
        return S1Angle.Max(0.471 * snap_radius_,        // Math.Sqrt(2) / 3 in the plane
            snap_radius_ - S1Angle.FromDegrees(S2.M_SQRT1_2 * to_degrees_));
    }

    // For IntLatLng snapping, the minimum separation between edges and
    // non-incident vertices depends on level() and snap_radius().  It can
    // be as low as 0.222 * snap_radius(), but is typically 0.39 * snap_radius()
    // or more.
    public override S1Angle MinEdgeVertexSeparation()
    {
        // Similar to min_vertex_separation(), in this case we have three bounds:
        // one is a constant bound, one is proportional to snap_radius, and one
        // approaches 0.5 * snap_radius asymptotically.
        //
        // 1. Constant bound: In the plane, the worst-case configuration has an
        //    edge-vertex separation of ((1 / Math.Sqrt(13)) * to_degrees_) degrees.
        //    The unit test verifies this, except that on the sphere the ratio is
        //    slightly lower when small exponents such as E1 are used
        //    (0.2772589 vs 0.2773501).
        //
        // 2. Proportional bound: In the plane, the worst-case configuration has an
        //    edge-vertex separation of (2 / 9) * snap_radius (0.222222222222).  The
        //    unit test verifies this, except that on the sphere the bound can be
        //    slightly worse with large exponents (e.g., E9) due to small numerical
        //    errors (0.222222126756717).
        //
        // 3. Best asymptotic bound: If snap_radius() is large compared to the
        //    minimum snap radius, then the best bound is achieved by 3 sites on a
        //    circular arc of radius "snap_radius", spaced "min_vertex_separation"
        //    apart (see S2CellIdSnapFunction.min_edge_vertex_separation).  This
        //    bound approaches 0.5 * snap_radius as the snap radius becomes large
        //    relative to the grid spacing.

        S1Angle vertex_sep = MinVertexSeparation();
        return S1Angle.Max(0.277 * S1Angle.FromDegrees(to_degrees_),  // 1/Math.Sqrt(13) in the plane
                   S1Angle.Max(0.222 * snap_radius_,               // 2/9 in the plane
                       0.5 * (vertex_sep / snap_radius_) * vertex_sep));
    }
    public override S2Point SnapPoint(S2Point point)
    {
        MyDebug.Assert(Exponent >= 0, "Make sure the snap function was initialized.");
        var input = new S2LatLng(point);
        var lat = Math.Round(input.LatDegrees() * from_degrees_);
        var lng = Math.Round(input.LngDegrees() * from_degrees_);
        return S2LatLng.FromDegrees(lat * to_degrees_, lng * to_degrees_).ToPoint();
    }
    public override object CustomClone() => new IntLatLngSnapFunction(this);

    private double from_degrees_, to_degrees_;
}

// A SnapFunction restricts the locations of the output vertices.  For
// example, there are predefined snap functions that require vertices to be
// located at S2CellId centers or at E5/E6/E7 coordinates.  The SnapFunction
// can also specify a minimum spacing between vertices (the "snap radius").
//
// A SnapFunction defines the following methods:
//
// 1. The SnapPoint() method, which snaps a point P to a nearby point (the
//    "candidate snap site").  Any point may be returned, including P
//    itself (this is the "identity snap function").
//
// 2. "snap_radius", the maximum distance that vertices can move when
//    snapped.  The snap_radius must be at least as large as the maximum
//    distance between P and SnapPoint(P) for any point P.
//
//    Note that the maximum distance that edge interiors can move when
//    snapped is slightly larger than "snap_radius", and is returned by the
//    function S2Builder::Options::max_edge_deviation() (see there for
//    details).
//
// 3. "min_vertex_separation", the guaranteed minimum distance between
//    vertices in the output.  This is generally a fraction of
//    "snap_radius" where the fraction depends on the snap function.
//
// 4. "min_edge_vertex_separation", the guaranteed minimum distance between
//    edges and non-incident vertices in the output.  This is generally a
//    fraction of "snap_radius" where the fraction depends on the snap
//    function.
//
// It is important to note that SnapPoint() does not define the actual
// mapping from input vertices to output vertices, since the points it
// returns (the candidate snap sites) are further filtered to ensure that
// they are separated by at least the snap radius.  For example, if you
// specify E7 coordinates (2cm resolution) and a snap radius of 10m, then a
// subset of points returned by SnapPoint will be chosen (the "snap sites"),
// and each input vertex will be mapped to the closest site.  Therefore you
// cannot assume that P is necessarily snapped to SnapPoint(P).
//
// S2Builder makes the following guarantees:
//
// 1. Every vertex is at a location returned by SnapPoint().
//
// 2. Vertices are within "snap_radius" of the corresponding input vertex.
//
// 3. Edges are within "max_edge_deviation" of the corresponding input edge
//    (a distance slightly larger than "snap_radius").
//
// 4. Vertices are separated by at least "min_vertex_separation"
//    (a fraction of "snap_radius" that depends on the snap function).
//
// 5. Edges and non-incident vertices are separated by at least
//    "min_edge_vertex_separation" (a fraction of "snap_radius").
//
// 6. Vertex and edge locations do not change unless one of the conditions
//    above is not already met (idempotency / stability).
//
// 7. The topology of the input geometry is preserved (up to the creation
//    of degeneracies).  This means that there exists a continuous
//    deformation from the input to the output such that no vertex
//    crosses an edge.
public abstract class SnapFunction : ICustomCloneable
{
    // The maximum supported snap radius (equivalent to about 7800km).
    //
    // The maximum snap radius is just large enough to support snapping to
    // S2CellId level 0.  It is equivalent to 7800km on the Earth's surface.
    //
    // This value can't be larger than 85.7 degrees without changing the code
    // related to min_edge_length_to_split_ca_, and increasing it to 90 degrees
    // or more would most likely require significant changes to the algorithm.
    public static readonly S1Angle kMaxSnapRadius = S1Angle.FromDegrees(70);

    // The maximum distance that vertices can move when snapped.  The snap
    // radius can be any value between zero and SnapFunction::kMaxSnapRadius().
    //
    // If the snap radius is zero, then vertices are snapped together only if
    // they are identical.  Edges will not be snapped to any vertices other
    // than their endpoints, even if there are vertices whose distance to the
    // edge is zero, unless split_crossing_edges() is true (see below).
    //
    // REQUIRES: snap_radius() <= SnapFunction.kMaxSnapRadius
    public abstract S1Angle SnapRadius { get; set; }

    // The guaranteed minimum distance between vertices in the output.
    // This is generally some fraction of "snap_radius".
    public abstract S1Angle MinVertexSeparation();

    // The guaranteed minimum spacing between edges and non-incident vertices
    // in the output.  This is generally some fraction of "snap_radius".
    public abstract S1Angle MinEdgeVertexSeparation();

    // Returns a candidate snap site for the given point.  The final vertex
    // locations are a subset of the snap sites returned by this function
    // (spaced at least "min_vertex_separation" apart).
    //
    // The only requirement is that SnapPoint(x) must return a point whose
    // distance from "x" is no greater than "snap_radius".
    public abstract S2Point SnapPoint(S2Point point);

    // Returns a deep copy of this SnapFunction.
    public abstract object CustomClone();
}
