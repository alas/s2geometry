using System;
namespace S2Geometry
{
    // Defines various angle and area measures on the sphere.  Also see
    // s2edge_distances.h, s2loop_measures.h, and s2polyline_measures.h.
    public static class S2Measures
    {
        // Return the interior angle at the vertex B in the triangle ABC.  The
        // return value is always in the range [0, Pi].  All points should be
        // normalized.  Ensures that Angle(a,b,c) == Angle(c,b,a) for all a,b,c.
        //
        // The angle is undefined if A or C is diametrically opposite from B, and
        // becomes numerically unstable as the length of edge AB or BC approaches
        // 180 degrees.
        public static double Angle(S2Point a, S2Point b, S2Point c)
        {
            // RobustCrossProd() is necessary to get good accuracy when two of the input
            // points are very close together.
            return S2PointUtil.RobustCrossProd(a, b).Angle(S2PointUtil.RobustCrossProd(c, b));
        }

        // Return the exterior angle at vertex B in the triangle ABC.  The return
        // value is positive if ABC is counterclockwise and negative otherwise.  If
        // you imagine a constant walking from A to B to C, this is the angle that the
        // constant turns at vertex B (positive = left = CCW, negative = right = CW).
        // This quantity is also known as the "geodesic curvature" at B.
        //
        // Ensures that TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all distinct
        // a,b,c. The result is undefined if (a == b || b == c), but is either
        // -Pi or Pi if (a == c).  All points should be normalized.
        public static double TurnAngle(S2Point a, S2Point b, S2Point c)
        {
            // We use RobustCrossProd() to get good accuracy when two points are very
            // close together, and Sign() to ensure that the sign is correct for
            // turns that are close to 180 degrees.
            //
            // Unfortunately we can't save RobustCrossProd(a, b) and pass it as the
            // optional 4th argument to Sign(), because Sign() requires a.CrossProd(b)
            // exactly (the robust version differs in magnitude).
            double angle = S2PointUtil.RobustCrossProd(a, b).Angle(S2PointUtil.RobustCrossProd(b, c));

            // Don't return Sign() * angle because it is legal to have (a == c).
            return (S2Pred.Sign(a, b, c) > 0) ? angle : -angle;
        }

        // Return the area of triangle ABC.  This method combines two different
        // algorithms to get accurate results for both large and small triangles.
        // The maximum error is about 5e-15 (about 0.25 square meters on the Earth's
        // surface), the same as GirardArea() below, but unlike that method it is
        // also accurate for small triangles.  Example: when the true area is 100
        // square meters, Area() yields an error about 1 trillion times smaller than
        // GirardArea().
        //
        // All points should be unit length, and no two points should be antipodal.
        // The area is always positive.
        public static double Area(S2Point a, S2Point b, S2Point c)
        {
            Assert.True(a.IsUnitLength);
            Assert.True(b.IsUnitLength);
            Assert.True(c.IsUnitLength);
            // This method is based on l'Huilier's theorem,
            //
            //   tan(E/4) = Math.Sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
            //
            // where E is the spherical excess of the triangle (i.e. its area),
            //       a, b, c, are the side lengths, and
            //       s is the semiperimeter (a + b + c) / 2 .
            //
            // The only significant source of error using l'Huilier's method is the
            // cancellation error of the terms (s-a), (s-b), (s-c).  This leads to a
            // *relative* error of about 1e-16 * s / Math.Min(s-a, s-b, s-c).  This compares
            // to a relative error of about 1e-15 / E using Girard's formula, where E is
            // the true area of the triangle.  Girard's formula can be even worse than
            // this for very small triangles, e.g. a triangle with a true area of 1e-30
            // might evaluate to 1e-5.
            //
            // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
            // dmin = Math.Min(s-a, s-b, s-c).  This basically includes all triangles
            // except for extremely long and skinny ones.
            //
            // Since we don't know E, we would like a conservative upper bound on
            // the triangle area in terms of s and dmin.  It's possible to show that
            // E <= k1 * s * Math.Sqrt(s * dmin), where k1 = 2*Math.Sqrt(3)/Pi (about 1).
            // Using this, it's easy to show that we should always use l'Huilier's
            // method if dmin >= k2 * s^5, where k2 is about 1e-2.  Furthermore,
            // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
            // k3 is about 0.1.  Since the best case error using Girard's formula
            // is about 1e-15, this means that we shouldn't even consider it unless
            // s >= 3e-4 or so.
            //
            // TODO(ericv): Implement rigorous error bounds (analysis already done).
            double sa = b.Angle(c);
            double sb = c.Angle(a);
            double sc = a.Angle(b);
            double s = 0.5 * (sa + sb + sc);
            if (s >= 3e-4)
            {
                // Consider whether Girard's formula might be more accurate.
                double s2 = s * s;
                double dmin = s - Math.Max(sa, Math.Max(sb, sc));
                if (dmin < 1e-2 * s * s2 * s2)
                {
                    // This triangle is skinny enough to consider using Girard's formula.
                    // We increase the area by the approximate maximum error in the Girard
                    // calculation in order to ensure that this test is conservative.
                    double area = GirardArea(a, b, c);
                    if (dmin < s * (0.1 * (area + 5e-15))) return area;
                }
            }
            // Use l'Huilier's formula.
            return 4 * Math.Atan(Math.Sqrt(Math.Max(0.0, Math.Tan(0.5 * s) * Math.Tan(0.5 * (s - sa)) *
                                     Math.Tan(0.5 * (s - sb)) * Math.Tan(0.5 * (s - sc)))));
        }

        // Return the area of the triangle computed using Girard's formula.  All
        // points should be unit length, and no two points should be antipodal.
        //
        // This method is about twice as fast as Area() but has poor relative
        // accuracy for small triangles.  The maximum error is about 5e-15 (about
        // 0.25 square meters on the Earth's surface) and the average error is about
        // 1e-15.  These bounds apply to triangles of any size, even as the maximum
        // edge length of the triangle approaches 180 degrees.  But note that for
        // such triangles, tiny perturbations of the input points can change the
        // true mathematical area dramatically.
        public static double GirardArea(S2Point a, S2Point b, S2Point c)
        {
            // This is equivalent to the usual Girard's formula but is slightly more
            // accurate, faster to compute, and handles a == b == c without a special
            // case.  RobustCrossProd() is necessary to get good accuracy when two of
            // the input points are very close together.

            var ab = S2PointUtil.RobustCrossProd(a, b);
            var bc = S2PointUtil.RobustCrossProd(b, c);
            var ac = S2PointUtil.RobustCrossProd(a, c);
            return Math.Max(0.0, ab.Angle(ac) - ab.Angle(bc) + bc.Angle(ac));
        }

        // Like Area(), but returns a positive value for counterclockwise triangles
        // and a negative value otherwise.
        public static double SignedArea(S2Point a, S2Point b, S2Point c)
        {
            return S2Pred.Sign(a, b, c) * Area(a, b, c);
        }
    }
}
