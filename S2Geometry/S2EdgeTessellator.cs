// Given an edge in some 2D projection (e.g., Mercator), S2EdgeTessellator
// converts the edge into a chain of spherical geodesic edges such that the
// maximum distance between the original edge and the geodesic edge chain is
// at most "tolerance".  Similarly, it can convert a spherical geodesic edge
// into a chain of edges in a given 2D projection such that the maximum
// distance between the geodesic edge and the chain of projected edges is at
// most "tolerance".
// </summary>
// <remarks>
// Tessellation is implemented by subdividing the edge until the estimated
// maximum error is below the given tolerance.  Estimating error is a hard
// problem, especially when the only methods available are point evaluation of
// the projection and its inverse.  (These are the only methods that
// S2.Projection provides, which makes it easier and less error-prone to
// implement new projections.)
//
// One technique that significantly increases robustness is to treat the
// geodesic and projected edges as parametric curves rather than geometric ones.
// Given a spherical edge AB and a projection p:S2.R2, let f(t) be the
// normalized arc length parametrization of AB and let g(t) be the normalized
// arc length parameterization of the projected edge p(A)p(B).  (In other words,
// f(0)=A, f(1)=B, g(0)=p(A), g(1)=p(B).)  We now define the geometric error as
// the maximum distance from the point p^-1(g(t)) to the geodesic edge AB for
// any t in [0,1], where p^-1 denotes the inverse projection.  In other words,
// the geometric error is the maximum distance from any point on the projected
// edge (mapped back onto the sphere) to the geodesic edge AB.  On the other
// hand we define the parametric error as the maximum distance between the
// points f(t) and p^-1(g(t)) for any t in [0,1], i.e. the maximum distance
// (measured on the sphere) between the geodesic and projected points at the
// same interpolation fraction t.
//
// The easiest way to estimate the parametric error is to simply evaluate both
// edges at their midpoints and measure the distance between them (the "midpoint
// method").  This is very fast and works quite well for most edges, however it
// has one major drawback: it doesn't handle points of inflection (i.e., points
// where the curvature changes sign).  For example, edges in the Mercator and
// Plate Carree projections always curve towards the equator relative to the
// corresponding geodesic edge, so in these projections there is a point of
// inflection whenever the projected edge crosses the equator.  The worst case
// occurs when the edge endpoints have different longitudes but the same
// absolute latitude, since in that case the error is non-zero but the edges
// have exactly the same midpoint (on the equator).
//
// One solution to this problem is to split the input edges at all inflection
// points (i.e., along the equator in the case of the Mercator and Plate Carree
// projections).  However for general projections these inflection points can
// occur anywhere on the sphere (e.g., consider the Transverse Mercator
// projection).  This could be addressed by adding methods to the S2Projection
// interface to split edges at inflection points but this would make it harder
// and more error-prone to implement new projections.
//
// Another problem with this approach is that the midpoint method sometimes
// underestimates the true error even when edges do not cross the equator.  For
// the Plate Carree and Mercator projections, the midpoint method can
// underestimate the error by up to 3%.
//
// Both of these problems can be solved as follows.  We assume that the error
// can be modeled as a convex combination of two worst-case functions, one
// where the error is maximized at the edge midpoint and another where the
// error is *minimized* (i.e., zero) at the edge midpoint.  For example, we
// could choose these functions as:
//
//    E1(x) = 1 - x^2
//    E2(x) = x * (1 - x^2)
//
// where for convenience we use an interpolation parameter "x" in the range
// [-1, 1] rather than the original "t" in the range [0, 1].  Note that both
// error functions must have roots at x = {-1, 1} since the error must be zero
// at the edge endpoints.  E1 is simply a parabola whose maximum value is 1
// attained at x = 0, while E2 is a cubic with an additional root at x = 0,
// and whose maximum value is 2 * Math.Sqrt(3) / 9 attained at x = 1 / Math.Sqrt(3).
//
// Next, it is convenient to scale these functions so that the both have a
// maximum value of 1.  E1 already satisfies this requirement, and we simply
// redefine E2 as
//
//   E2(x) = x * (1 - x^2) / (2 * Math.Sqrt(3) / 9)
//
// Now define x0 to be the point where these two functions intersect, i.e. the
// point in the range (-1, 1) where E1(x0) = E2(x0).  This value has the very
// convenient property that if we evaluate the actual error E(x0), then the
// maximum error on the entire interval [-1, 1] is bounded by
//
//   E(x) <= E(x0) / E1(x0)
//
// since whether the error is modeled using E1 or E2, the resulting function
// has the same maximum value (namely E(x0) / E1(x0)).  If it is modeled as
// some other convex combination of E1 and E2, the maximum value can only
// decrease.
//
// Finally, since E2 is not symmetric about the y-axis, we must also allow for
// the possibility that the error is a convex combination of E1 and -E2.  This
// can be handled by evaluating the error at E(-x0) as well, and then
// computing the final error bound as
//
//   E(x) <= Math.Max(E(x0), E(-x0)) / E1(x0) .
//
// Effectively, this method is simply evaluating the error at two points about
// 1/3 and 2/3 of the way along the edges, and then scaling the maximum of
// these two errors by a constant factor.  Intuitively, the reason this works
// is that if the two edges cross somewhere in the interior, then at least one
// of these points will be far from the crossing.
//
// The actual algorithm implemented below has some additional refinements.
// First, edges longer than 90 degrees are always subdivided; this avoids
// various unusual situations that can happen with very long edges, and there
// is really no reason to avoid adding vertices to edges that are so long.
//
// Second, the error function E1 above needs to be modified to take into
// account spherical distortions.  (It turns out that spherical distortions are
// beneficial in the case of E2, i.e. they only make its error estimates
// slightly more conservative.)  To do this, we model E1 as the maximum error
// in a Plate Carree edge of length 90 degrees or less.  This turns out to be
// an edge from 45:-90 to 45:90 (in lat:lng format).  The corresponding error
// as a function of "x" in the range [-1, 1] can be computed as the distance
// between the the the Plate Caree edge point (45, 90 * x) and the geodesic
// edge point (90 - 45 * abs(x), 90 * sgn(x)).  Using the Haversine formula,
// the corresponding function E1 (normalized to have a maximum value of 1) is:
//
//   E1(x) =
//     asin(Math.Sqrt(sin(Pi / 8 * (1 - x)) ^ 2 +
//               sin(Pi / 4 * (1 - x)) ^ 2 * cos(Pi / 4) * sin(Pi / 4 * x))) /
//     asin(Math.Sqrt((1 - 1 / Math.Sqrt(2)) / 2))
//
// Note that this function does not need to be evaluated at runtime, it
// simply affects the calculation of the value x0 where E1(x0) = E2(x0)
// and the corresponding scaling factor C = 1 / E1(x0).
//
// ------------------------------------------------------------------
//
// In the case of the Mercator and Plate Carree projections this strategy
// produces a conservative upper bound (verified using 10 million random
// edges).  Furthermore the bound is nearly tight; the scaling constant is
// C = 1.19289, whereas the maximum observed value was 1.19254.
//
// Compared to the simpler midpoint evaluation method, this strategy requires
// more function evaluations (currently twice as many, but with a smarter
// tessellation algorithm it will only be 50% more).  It also results in a
// small amount of additional tessellation (about 1.5%) compared to the
// midpoint method, but this is due almost entirely to the fact that the
// midpoint method does not yield conservative error estimates.
//
// For random edges with a tolerance of 1 meter, the expected amount of
// overtessellation is as follows:
//
//                   Midpoint Method    Cubic Method
//   Plate Carree               1.8%            3.0%
//   Mercator                  15.8%           17.4%

namespace S2Geometry;

public class S2EdgeTessellator
{
    // Constructs an S2EdgeTessellator using the given projection and error
    // tolerance.  The projection object must be valid for the entire lifetime
    // of this object.  (Projections are typically declared once and reused.)
    //
    // Method            | Input                  | Output
    // ------------------|------------------------|-----------------------
    // AppendProjected   | S2 geodesics           | Planar projected edges
    // AppendUnprojected | Planar projected edges | S2 geodesics
    public S2EdgeTessellator(Projection projection, S1Angle tolerance)
    {
        proj_ = projection;
        if (tolerance < kMinTolerance())
            throw new ArgumentException("Tolerance too small");

        // Rather than scaling the error estimate as described above, instead we scale
        // the tolerance.  See algorithm description at the top of this file.
        scaled_tolerance_ = new S1ChordAngle(kScaleFactor * S1Angle.Max(tolerance, kMinTolerance()));
    }

    // Converts the spherical geodesic edge AB to a chain of planar edges in the
    // given projection and appends the corresponding vertices to "vertices".
    //
    // This method can be called multiple times with the same output vector to
    // convert an entire polyline or loop.  All vertices of the first edge are
    // appended, but the first vertex of each subsequent edge is omitted (and
    // must match the last vertex of the previous edge).
    //
    // If the given projection has one or more coordinate axes that "wrap", then
    // every vertex's coordinates will be as close as possible to the previous
    // vertex's coordinates.  Note that this may yield vertices whose
    // coordinates are outside the usual range.  For example, tessellating the
    // edge (0:170, 0:-170) (in lat:lng notation) yields (0:170, 0:190).
    public void AppendProjected(S2Point a, S2Point b, List<R2Point> vertices)
    {
        var pa = proj_.Project(a);
        if (!vertices.Any())
        {
            vertices.Add(pa);
        }
        else
        {
            pa = proj_.WrapDestination(vertices.Last(), pa);
            System.Diagnostics.Debug.Assert(vertices.Last() == pa); // Appended edges must form a chain
        }
        var pb = proj_.Project(b);
        AppendProjected(pa, a, pb, b, vertices);
    }

    // Converts the planar edge AB in the given projection to a chain of
    // spherical geodesic edges and appends the vertices to "vertices".
    //
    // This method can be called multiple times with the same output vector to
    // convert an entire polyline or loop.  All vertices of the first edge are
    // appended, but the first vertex of each subsequent edge is omitted (and is
    // required to match that last vertex of the previous edge).
    //
    // Note that to construct an S2Loop, you must call vertices.pop_back() at
    // the very end to eliminate the duplicate first and last vertex.  Note also
    // that if the given projection involves coordinate "wrapping" (e.g. across
    // the 180 degree meridian) then the first and last vertices may not be
    // exactly the same.
    public void AppendUnprojected(R2Point a, R2Point b, List<S2Point> vertices)
    {
        var pointA = proj_.Unproject(a);
        var pointB = proj_.Unproject(b);
        if (!vertices.Any())
        {
            vertices.Add(pointA);
        }
        else
        {
            // Note that coordinate wrapping can create a small amount of error.  For
            // example in the edge chain "0:-175, 0:179, 0:-177", the first edge is
            // transformed into "0:-175, 0:-181" while the second is transformed into
            // "0:179, 0:183".  The two coordinate pairs for the middle vertex
            // ("0:-181" and "0:179") may not yield exactly the same S2Point.
            System.Diagnostics.Debug.Assert(S2.ApproxEquals(vertices.Last(), pointA)); // Appended edges must form a chain
        }
        AppendUnprojected(a, pointA, b, pointB, vertices);
    }

    // Returns the minimum supported tolerance (which corresponds to a distance
    // less than one micrometer on the Earth's surface).
    public static S1Angle kMinTolerance()
    {
        // This distance is less than 1 micrometer on the Earth's surface, but is
        // still much larger than the expected projection and interpolation errors.
        return S1Angle.FromRadians(1e-13);
    }

    private S1ChordAngle EstimateMaxError(R2Point pa, S2Point a, R2Point pb, S2Point b)
    {
        // See the algorithm description at the top of this file.

        // We always tessellate edges longer than 90 degrees on the sphere, since the
        // approximation below is not robust enough to handle such edges.
        if (a.DotProd(b) < -1e-14) return S1ChordAngle.Infinity;

        const double t1 = kInterpolationFraction;
        const double t2 = 1 - kInterpolationFraction;
        S2Point mid1 = S2.Interpolate(a, b, t1);
        S2Point mid2 = S2.Interpolate(a, b, t2);
        S2Point pmid1 = proj_.Unproject(Projection.Interpolate(t1, pa, pb));
        S2Point pmid2 = proj_.Unproject(Projection.Interpolate(t2, pa, pb));
        return S1ChordAngle.Max(new S1ChordAngle(mid1, pmid1), new S1ChordAngle(mid2, pmid2));
    }

    // Like AppendProjected, but interpolates a projected edge and appends the
    // corresponding points on the sphere.
    private void AppendUnprojected(R2Point pa, S2Point a, R2Point pb, S2Point b, List<S2Point> vertices)
    {
        // See notes above regarding measuring the interpolation error.
        var pbTmp = proj_.WrapDestination(pa, pb);
        if (EstimateMaxError(pa, a, pbTmp, b) <= scaled_tolerance_)
        {
            vertices.Add(b);
        }
        else
        {
            R2Point pmid = Projection.Interpolate(0.5, pa, pbTmp);
            S2Point mid = proj_.Unproject(pmid);
            AppendUnprojected(pa, a, pmid, mid, vertices);
            AppendUnprojected(pmid, mid, pbTmp, b, vertices);
        }
    }

    // Given a geodesic edge AB, split the edge as necessary and append all
    // projected vertices except the first to "vertices".
    //
    // The maximum recursion depth is (Math.PI / kMinTolerance()) < 45, and the
    // frame size is small so stack overflow should not be an issue.
    private void AppendProjected(R2Point pa, S2Point a, R2Point pb, S2Point b, List<R2Point> vertices)
    {
        R2Point pbTmp = proj_.WrapDestination(pa, pb);
        if (EstimateMaxError(pa, a, pbTmp, b) <= scaled_tolerance_)
        {
            vertices.Add(pbTmp);
        }
        else
        {
            var mid = (a + b).Normalize();
            var pmid = proj_.WrapDestination(pa, proj_.Project(mid));
            AppendProjected(pa, a, pmid, mid, vertices);
            AppendProjected(pmid, mid, pbTmp, b, vertices);
        }
    }

    private readonly Projection proj_;

    // The given tolerance scaled by a constant fraction so that it can be
    // compared against the result returned by EstimateMaxError().
    private readonly S1ChordAngle scaled_tolerance_;

    // The interpolation fraction at which the two edges are evaluated in order to
    // measure the error between them.  (Edges are evaluated at two points measured
    // this fraction from either end.)  With respect to the algorithm description
    // above, this value is t0 = (1 - x0) / 2 in the range [0, 1] that corresponds
    // to x0 in the range [-1, 1] chosen such that E1(x0) == E2(x0).
    private const double kInterpolationFraction = 0.31215691082248312;

    // The following is the value of E1(x0) == E2(x0).
    private const double kScaleFactor = 0.83829992569888509;
}
