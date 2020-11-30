// Uncomment the following line for testing purposes only.
// # define S2_TEST_ DEGENERACIES

using System;
namespace S2Geometry
{
    // Defines additional operations for points on the unit sphere (in addition to
    // the standard vector operations defined in "util/math/vector.h").
    public static class S2PointUtil
    {
        // Return a unique "origin" on the sphere for operations that need a fixed
        // reference point.  In particular, this is the "point at infinity" used for
        // point-in-polygon testing (by counting the number of edge crossings).
        public static readonly S2Point Origin = 
#if S2_TEST_DEGENERACIES
            // This value makes polygon operations much slower, because it greatly
            // increases the number of degenerate cases that need to be handled using
            // S2Pred.ExpensiveSign().
            new S2Point(0, 0, 1);
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
            new S2Point(-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195);
#endif

        // Return true if two points are within the given distance of each other (this
        // is mainly useful for testing). It is an error if either point is a
        // zero-length vector (default S2Point), but this is only checked in debug mode.
        // In non-debug mode it will always return true.
        public static bool ApproxEquals(S2Point a, S2Point b)
        {
            return ApproxEquals(a, b, S1Angle.FromRadians(S2Constants.DoubleError));
        }
        public static bool ApproxEquals(S2Point a, S2Point b, S1Angle maxError)
        {
            Assert.True(a != S2Point.Empty);
            Assert.True(b != S2Point.Empty);
            return new S1Angle(a, b) <= maxError;
        }

        // Return a unit-length vector that is orthogonal to "a".  Satisfies
        // Ortho(-a) = -Ortho(a) for all a.
        //
        // Note that Vector3_d also defines an "Ortho" method, but this one is
        // preferred for use in S2 code because it explicitly tries to avoid result
        // result coordinates that are zero.  (This is a performance optimization that
        // reduces the amount of time spent in functions which handle degeneracies.)
        public static S2Point Ortho(S2Point a)
        {
#if S2_TEST_DEGENERACIES
            // Vector3.Ortho() always returns a point on the X-Y, Y-Z, or X-Z planes.
            // This leads to many more degenerate cases in polygon operations.
            return a.Ortho;
#else
            int k = a.LargestAbsComponent - 1;
            if (k < 0) k = 2;
            var temp = new double[] { 0.012, 0.0053, 0.00457 };
            temp[k] = 1;
            var s2 = new S2Point(temp);
            return a.CrossProd(s2).Normalized;
#endif
        }

        // Return a vector "c" that is orthogonal to the given unit-length vectors
        // "a" and "b".  This function is similar to a.CrossProd(b) except that it
        // does a better job of ensuring orthogonality when "a" is nearly parallel
        // to "b", and it returns a non-zero result even when a == b or a == -b.
        //
        // It satisfies the following properties (RCP == RobustCrossProd):
        //
        //   (1) RCP(a,b) != 0 for all a, b
        //   (2) RCP(b,a) == -RCP(a,b) unless a == b or a == -b
        //   (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
        //   (4) RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
        //
        // The result is not guaranteed to be unit length.
        public static S2Point RobustCrossProd(S2Point a, S2Point b)
        {
            // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
            // approaches zero.  This leads to situations where a.CrossProd(b) is not
            // very orthogonal to "a" and/or "b".  We could fix this using Gram-Schmidt,
            // but we also want b.RobustCrossProd(a) == -a.RobustCrossProd(b).
            //
            // The easiest fix is to just compute the cross product of (b+a) and (b-a).
            // Mathematically, this cross product is exactly twice the cross product of
            // "a" and "b", but it has the numerical advantage that (b+a) and (b-a)
            // are always perpendicular (since "a" and "b" are unit length).  This
            // yields a result that is nearly orthogonal to both "a" and "b" even if
            // these two values differ only in the lowest bit of one component.

            Assert.True(a.IsUnitLength);
            Assert.True(b.IsUnitLength);
            S2Point x = (b + a).CrossProd(b - a);
            if (x != S2Point.Empty) return x;

            // The only result that makes sense mathematically is to return zero, but
            // we find it more convenient to return an arbitrary orthogonal vector.
            return Ortho(a);
        }

        // Rotate the given point about the given axis by the given angle.  "p" and
        // "axis" must be unit length; "angle" has no restrictions (e.g., it can be
        // positive, negative, greater than 360 degrees, etc).
        public static S2Point Rotate(S2Point p, S2Point axis, S1Angle angle)
        {
            Assert.True(p.IsUnitLength);
            Assert.True(axis.IsUnitLength);
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
            return (Math.Cos(angle.Radians) * dx + Math.Sin(angle.Radians) * dy + center).Normalized;
        }

        // Extend the given point "z" on the unit sphere into a right-handed
        // coordinate frame of unit-length column vectors m = (x,y,z).  Note that the
        // vectors (x,y) are an orthonormal frame for the tangent space at "z", while
        // "z" itself is an orthonormal frame for the normal space at "z".
        public static S2PointVector3 GetFrame(S2Point z)
        {
            Assert.True(z.IsUnitLength);

            var ortho = Ortho(z);
            return new S2PointVector3(
                ortho.CrossProd(z), // Already unit-length.
                ortho,
                z);
        }

        // Given an orthonormal basis "m" of column vectors and a point "p", return
        // the coordinates of "p" with respect to the basis "m".  The resulting
        // point "q" satisfies the identity (m * q == p).
        public static S2Point ToFrame(S2PointVector3 m, S2Point p)
        {
            // The inverse of an orthonormal matrix is its transpose.
            return m.Transpose() * p;
        }

        // Given an orthonormal basis "m" of column vectors and a point "q" with
        // respect to that basis, return the equivalent point "p" with respect to
        // the standard axis-aligned basis.  The result satisfies (p == m * q).
        public static S2Point FromFrame(S2PointVector3 m, S2Point q)
        {
            return m * q;
        }

        // Return true if the points A, B, C are strictly counterclockwise.  Return
        // false if the points are clockwise or collinear (i.e. if they are all
        // contained on some great circle).
        //
        // Due to numerical errors, situations may arise that are mathematically
        // impossible, e.g. ABC may be considered strictly CCW while BCA is not.
        // However, the implementation guarantees the following:
        //
        //   If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
        [Obsolete(@"This Method is Deprecated, Use ""1 == S2Pred.Sign(..."" instead.", false)]
        public static bool SimpleCCW(S2Point a, S2Point b, S2Point c)
        {
            // We compute the signed volume of the parallelepiped ABC.  The usual
            // formula for this is (AxB).C, but we compute it here using (CxA).B
            // in order to ensure that ABC and CBA are not both CCW.  This follows
            // from the following identities (which are true numerically, not just
            // mathematically):
            //
            //     (1) x.CrossProd(y) == -(y.CrossProd(x))
            //     (2) (-x).DotProd(y) == -(x.DotProd(y))

            return c.CrossProd(a).DotProd(b) > 0;
        }
    }
}
