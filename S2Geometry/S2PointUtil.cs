// Uncomment the following line for testing purposes only.
// # define S2_TEST_ DEGENERACIES

// Defines additional operations for points on the unit sphere (in addition to
// the standard vector operations defined in "util/math/vector.h").

namespace S2Geometry;

public static partial class S2
{
    // Returns a unique "origin" on the sphere for operations that need a fixed
    // reference point.  In particular, this is the "point at infinity" used for
    // point-in-polygon testing (by counting the number of edge crossings).
    // To be clear, this is NOT (0,0,0), the origin of the coordinate system.
    public static readonly S2Point Origin =
#if S2_TEST_DEGENERACIES
            // This value makes polygon operations much slower, because it greatly
            // increases the number of degenerate cases that need to be handled using
            // S2Pred.ExpensiveSign().
            new(0, 0, 1);
#else
        // The origin should not be a point that is commonly used in edge tests in
        // order to avoid triggering code to handle degenerate cases.  (This rules
        // out the north and south poles.)  It should also not be on the boundary of
        // any low-level S2Cell for the same reason.
        //
        // The point chosen here is about 66km from the north pole towards the East
        // Siberian Sea.  See the unittest for more details.  It is written out
        // explicitly using floating-point literals because the optimizer doesn't
        // seem willing to evaluate Normalize() at compile time.
        new(-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195);
#endif

    // Returns true if two points are within the given distance of each other
    // (this is mainly useful for testing).  It is an error if either point is a
    // zero-length vector (default S2Point), but this is only checked in debug
    // mode.  In non-debug mode it will always return true.
    public static bool ApproxEquals(S2Point a, S2Point b) =>
        ApproxEquals(a, b, S1Angle.FromRadians(S2.DoubleError));
    public static bool ApproxEquals(S2Point a, S2Point b, S1Angle maxError)
    {
        MyDebug.Assert(a != S2Point.Empty);
        MyDebug.Assert(b != S2Point.Empty);
        return new S1Angle(a, b) <= maxError;
    }

    // Returns a unit-length vector that is orthogonal to "a".  Satisfies
    // Ortho(-a) = -Ortho(a) for all a.
    //
    // Note that S2Point_d also defines an "Ortho" method, but this one is
    // preferred for use in S2 code because it explicitly tries to avoid result
    // coordinates that are zero.  (This is a performance optimization that
    // reduces the amount of time spent in functions that handle degeneracies.)
    public static S2Point Ortho(S2Point a)
    {
#if S2_TEST_DEGENERACIES
            // S2Point.Ortho() always returns a point on the X-Y, Y-Z, or X-Z planes.
            // This leads to many more degenerate cases in polygon operations.
            return a.Ortho;
#else
        int k = a.LargestAbsComponent() - 1;
        if (k < 0) k = 2;
        double[] temp = [0.012, 0.0053, 0.00457];
        temp[k] = 1;
        S2Point s2 = new(temp);
        return a.CrossProd(s2).Normalize();
#endif
    }

    // Returns a unit-length vector used as the reference direction for deciding
    // whether a polygon with semi-open boundaries contains the given vertex "a"
    // (see S2ContainsVertexQuery).  The result is unit length and is guaranteed
    // to be different from the given point "a".
    public static S2Point RefDir(S2Point a) => S2.Ortho(a);

    // Rotates the given point about the given axis by the given angle.  "p" and
    // "axis" must be unit length; "angle" has no restrictions (e.g., it can be
    // positive, negative, greater than 360 degrees, etc).
    //
    // See also the closely related functions S2::GetPointOnRay() and
    // S2::GetPointOnLine(), which are declared in s2edge_distances.h.
    public static S2Point Rotate(S2Point p, S2Point axis, S1Angle angle)
    {
        MyDebug.Assert(p.IsUnitLength());
        MyDebug.Assert(axis.IsUnitLength());
        // Let M be the plane through P that is perpendicular to "axis", and let
        // "center" be the point where M intersects "axis".  We construct a
        // right-handed orthogonal frame (dx, dy, center) such that "dx" is the
        // vector from "center" to P, and "dy" has the same length as "dx".  The
        // result can then be expressed as (cos(angle)*dx + sin(angle)*dy + center).
        var center = p.DotProd(axis) * axis;
        var dx = p - center;
        var dy = axis.CrossProd(p);
        // Mathematically the result is unit length, but normalization is necessary
        // to ensure that numerical errors don't accumulate.
        return (Math.Cos(angle.Radians) * dx + Math.Sin(angle.Radians) * dy + center).Normalize();
    }

    // Extends the given point "z" on the unit sphere into a right-handed
    // coordinate frame of unit-length column vectors m = (x,y,z).  Note that the
    // vectors (x,y) are an orthonormal frame for the tangent space at "z", while
    // "z" itself is an orthonormal frame for the normal space at "z".
    public static S2PointS2Point GetFrame(S2Point z)
    {
        MyDebug.Assert(z.IsUnitLength());

        var ortho = Ortho(z);
        var prod = ortho.CrossProd(z); // Already unit-length.

        return new S2PointS2Point(
            prod[0], ortho[0], z[0],
            prod[1], ortho[1], z[1],
            prod[2], ortho[2], z[2]);
    }

    // Given an orthonormal basis "m" of column vectors and a point "p", returns
    // the coordinates of "p" with respect to the basis "m".  The resulting
    // point "q" satisfies the identity (m * q == p).
    public static S2Point ToFrame(S2PointS2Point m, S2Point p) =>
        // The inverse of an orthonormal matrix is its transpose.
        m.Transpose() * p;

    // Given an orthonormal basis "m" of column vectors and a point "q" with
    // respect to that basis, return the equivalent point "p" with respect to
    // the standard axis-aligned basis.  The result satisfies (p == m * q).
    public static S2Point FromFrame(S2PointS2Point m, S2Point q) =>
        m * q;
}
