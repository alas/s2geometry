using System;
namespace S2Geometry
{
    // Defines a collection of functions for computing the distance to an edge,
    // interpolating along an edge, projecting points onto edges, etc.
    public static class S2EdgeDistances
    {
        /////////////////////////////////////////////////////////////////////////////
        ///////////////            (point, edge) functions            ///////////////

        // Returns the minimum distance from X to any point on the edge AB.  All
        // arguments should be unit length.  The result is very accurate for small
        // distances but may have some numerical error if the distance is large
        // (approximately Pi/2 or greater).  The case A == B is handled correctly.
        //
        // If you want to compare a distance against a fixed threshold, e.g.
        //    if (S2EdgeDistances.GetDistance(x, a, b) < limit)
        // then it is significantly faster to use UpdateMinDistance() below.
        public static S1Angle GetDistance(S2Point x, S2Point a, S2Point b)
        {
            var min_dist = S1ChordAngle.Zero;
            AlwaysUpdateMinDistance(x, a, b, ref min_dist, true);
            return min_dist.ToAngle();
        }

        // Returns true if the distance from X to the edge AB is less than "limit".
        // (Specify limit.Successor() for "less than or equal to".)  This method is
        // significantly faster than GetDistance().  If you want to compare against a
        // fixed S1Angle, you should convert it to an S1ChordAngle once and save the
        // value, since this step is relatively expensive.
        //
        // See S2Pred.CompareEdgeDistance() for an exact version of this predicate.
        public static bool IsDistanceLess(S2Point x, S2Point a, S2Point b, S1ChordAngle limit)
        {
            return UpdateMinDistance(x, a, b, ref limit);
        }

        // If the distance from X to the edge AB is less than "min_dist", this
        // method updates "min_dist" and returns true.  Otherwise it returns false.
        // The case A == B is handled correctly.
        //
        // Use this method when you want to compute many distances and keep track of
        // the minimum.  It is significantly faster than using GetDistance(),
        // because (1) using S1ChordAngle is much faster than S1Angle, and (2) it
        // can save a lot of work by not actually computing the distance when it is
        // obviously larger than the current minimum.
        public static bool UpdateMinDistance(S2Point x, S2Point a, S2Point b, ref S1ChordAngle min_dist)
        {
            return AlwaysUpdateMinDistance(x, a, b, ref min_dist, false);
        }

        // If the maximum distance from X to the edge AB is greater than "max_dist",
        // this method updates "max_dist" and returns true.  Otherwise it returns false.
        // The case A == B is handled correctly.
        public static bool UpdateMaxDistance(S2Point x, S2Point a, S2Point b, ref S1ChordAngle max_dist)
        {
            var dist = S1ChordAngle.Max(new S1ChordAngle(x, a), new S1ChordAngle(x, b));
            if (dist > S1ChordAngle.Right)
            {
                AlwaysUpdateMinDistance(-x, a, b, ref dist, true);
                dist = S1ChordAngle.Straight - dist;
            }
            if (max_dist < dist)
            {
                max_dist = dist;
                return true;
            }

            return false;
        }

        // Returns the maximum error in the result of UpdateMinDistance (and
        // associated functions such as UpdateMinInteriorDistance, IsDistanceLess,
        // etc), assuming that all input points are normalized to within the bounds
        // guaranteed by S2Point.Normalize().  The error can be added or subtracted
        // from an S1ChordAngle "x" using x.PlusError(error).
        //
        // Note that accuracy goes down as the distance approaches 0 degrees or 180
        // degrees (for different reasons).  Near 0 degrees the error is acceptable
        // for all practical purposes (about 1.2e-15 radians ~= 8 nanometers).  For
        // exactly antipodal points the maximum error is quite high (0.5 meters), but
        // this error drops rapidly as the points move away from antipodality
        // (approximately 1 millimeter for points that are 50 meters from antipodal,
        // and 1 micrometer for points that are 50km from antipodal).
        //
        // TODO(ericv): Currently the error bound does not hold for edges whose
        // endpoints are antipodal to within about 1e-15 radians (less than 1 micron).
        // This could be fixed by extending S2.RobustCrossProd to use higher
        // precision when necessary.
        public static double GetUpdateMinDistanceMaxError(S1ChordAngle dist)
        {
            // There are two cases for the maximum error in UpdateMinDistance(),
            // depending on whether the closest point is interior to the edge.
            return Math.Max(GetUpdateMinInteriorDistanceMaxError(dist),
                       dist.GetS2PointConstructorMaxError());
        }

        // Returns true if the minimum distance from X to the edge AB is attained at
        // an interior point of AB (i.e., not an endpoint), and that distance is less
        // than "limit".  (Specify limit.Successor() for "less than or equal to".)
        public static bool IsInteriorDistanceLess(S2Point x, S2Point a, S2Point b, S1ChordAngle limit)
        {
            return UpdateMinInteriorDistance(x, a, b, ref limit);
        }

        // If the minimum distance from X to AB is attained at an interior point of AB
        // (i.e., not an endpoint), and that distance is less than "min_dist", then
        // this method updates "min_dist" and returns true.  Otherwise returns false.
        public static bool UpdateMinInteriorDistance(S2Point x, S2Point a, S2Point b, ref S1ChordAngle min_dist)
        {
            double xa2 = (x - a).Norm2, xb2 = (x - b).Norm2;
            return AlwaysUpdateMinInteriorDistance(x, a, b, xa2, xb2, ref min_dist, false);
        }

        // Returns the point along the edge AB that is closest to the point X.
        // The fractional distance of this point along the edge AB can be obtained
        // using GetDistanceFraction() above.  Requires that all vectors have
        // unit length.
        public static S2Point Project(S2Point x, S2Point a, S2Point b)
        {
            return Project(x, a, b, S2PointUtil.RobustCrossProd(a, b));
        }

        // A slightly more efficient version of Project() where the cross product of
        // the two endpoints has been precomputed.  The cross product does not need to
        // be normalized, but should be computed using S2.RobustCrossProd() for the
        // most accurate results.  Requires that x, a, and b have unit length.
        public static S2Point Project(S2Point x, S2Point a, S2Point b, S2Point a_cross_b)
        {
            Assert.True(a.IsUnitLength);
            Assert.True(b.IsUnitLength);
            Assert.True(x.IsUnitLength);

            // Find the closest point to X along the great circle through AB.
            S2Point p = x - (x.DotProd(a_cross_b) / a_cross_b.Norm2) * a_cross_b;

            // If this point is on the edge AB, then it's the closest point.
            //if (S2PointUtil.SimpleCCW(a_cross_b, a, p) && S2PointUtil.SimpleCCW(p, b, a_cross_b))
            if (S2Pred.Sign(a_cross_b, a, p) == 1 && S2Pred.Sign(p, b, a_cross_b) == 1)
            {
                return p.Normalized;
            }
            // Otherwise, the closest point is either A or B.
            return ((x - a).Norm2 <= (x - b).Norm2) ? a : b;
        }

        /////////////////////////////////////////////////////////////////////////////
        ///////////////         (point along edge) functions          ///////////////


        // Given a point X and an edge AB, returns the distance ratio AX / (AX + BX).
        // If X happens to be on the line segment AB, this is the fraction "t" such
        // that X == Interpolate(t, A, B).  Requires that A and B are distinct.
        public static double GetDistanceFraction(S2Point x, S2Point a0, S2Point a1)
        {
            Assert.True(a0 != a1);
            var d0 = x.Angle(a0);
            var d1 = x.Angle(a1);
            return d0 / (d0 + d1);
        }

        // Returns the point X along the line segment AB whose distance from A is the
        // given fraction "t" of the distance AB.  Does NOT require that "t" be
        // between 0 and 1.  Note that all distances are measured on the surface of
        // the sphere, so this is more complicated than just computing (1-t)*a + t*b
        // and normalizing the result.
        public static S2Point Interpolate(double t, S2Point a, S2Point b)
        {
            if (t == 0) return a;
            if (t == 1) return b;
            var ab = new S1Angle(a, b);
            return InterpolateAtDistance(t * ab, a, b);
        }

        // Like Interpolate(), except that the parameter "ax" represents the desired
        // distance from A to the result X rather than a fraction between 0 and 1.
        public static S2Point InterpolateAtDistance(S1Angle ax_angle, S2Point a, S2Point b)
        {
            var ax = ax_angle.Radians;

            Assert.True(a.IsUnitLength);
            Assert.True(b.IsUnitLength);

            // Use RobustCrossProd() to compute the tangent vector at A towards B.  The
            // result is always perpendicular to A, even if A=B or A=-B, but it is not
            // necessarily unit length.  (We effectively normalize it below.)
            var normal = S2PointUtil.RobustCrossProd(a, b);
            var tangent = normal.CrossProd(a);
            Assert.True(tangent != S2Point.Empty);

            // Now compute the appropriate linear combination of A and "tangent".  With
            // infinite precision the result would always be unit length, but we
            // normalize it anyway to ensure that the error is within acceptable bounds.
            // (Otherwise errors can build up when the result of one interpolation is
            // fed into another interpolation.)
            return (Math.Cos(ax) * a + (Math.Sin(ax) / tangent.Norm) * tangent).Normalized;
        }

        /////////////////////////////////////////////////////////////////////////////
        ///////////////            (edge, edge) functions             ///////////////


        // Like UpdateMinDistance(), but computes the minimum distance between the
        // given pair of edges.  (If the two edges cross, the distance is zero.)
        // The cases a0 == a1 and b0 == b1 are handled correctly.
        public static bool UpdateEdgePairMinDistance(S2Point a0, S2Point a1, S2Point b0, S2Point b1, ref S1ChordAngle min_dist)
        {
            if (min_dist == S1ChordAngle.Zero)
            {
                return false;
            }
            if (S2EdgeCrossings.CrossingSign(a0, a1, b0, b1) > 0)
            {
                min_dist = S1ChordAngle.Zero;
                return true;
            }
            // Otherwise, the minimum distance is achieved at an endpoint of at least
            // one of the two edges.  We use "|" rather than "||" below to ensure that
            // all four possibilities are always checked.
            //
            // The calculation below computes each of the six vertex-vertex distances
            // twice (this could be optimized).
            return (UpdateMinDistance(a0, b0, b1, ref min_dist) |
                    UpdateMinDistance(a1, b0, b1, ref min_dist) |
                    UpdateMinDistance(b0, a0, a1, ref min_dist) |
                    UpdateMinDistance(b1, a0, a1, ref min_dist));
        }

        // As above, but for maximum distances.  If one edge crosses the antipodal
        // reflection of the other, the distance is Pi.
        public static bool UpdateEdgePairMaxDistance(S2Point a0, S2Point a1, S2Point b0, S2Point b1, ref S1ChordAngle max_dist)
        {
            if (max_dist == S1ChordAngle.Straight)
            {
                return false;
            }
            if (S2EdgeCrossings.CrossingSign(a0, a1, -b0, -b1) > 0)
            {
                max_dist = S1ChordAngle.Straight;
                return true;
            }
            // Otherwise, the maximum distance is achieved at an endpoint of at least
            // one of the two edges.  We use "|" rather than "||" below to ensure that
            // all four possibilities are always checked.
            //
            // The calculation below computes each of the six vertex-vertex distances
            // twice (this could be optimized).
            return (UpdateMaxDistance(a0, b0, b1, ref max_dist) |
                    UpdateMaxDistance(a1, b0, b1, ref max_dist) |
                    UpdateMaxDistance(b0, a0, a1, ref max_dist) |
                    UpdateMaxDistance(b1, a0, a1, ref max_dist));
        }

        // Returns the pair of points (a, b) that achieves the minimum distance
        // between edges a0a1 and b0b1, where "a" is a point on a0a1 and "b" is a
        // point on b0b1.  If the two edges intersect, "a" and "b" are both equal to
        // the intersection point.  Handles a0 == a1 and b0 == b1 correctly.
        public static (S2Point, S2Point) GetEdgePairClosestPoints(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
        {
            if (S2EdgeCrossings.CrossingSign(a0, a1, b0, b1) > 0)
            {
                S2Point x = S2EdgeCrossings.GetIntersection(a0, a1, b0, b1, null);
                return (x, x);
            }
            // We save some work by first determining which vertex/edge pair achieves
            // the minimum distance, and then computing the closest point on that edge.
            var min_dist = S1ChordAngle.Zero;
            AlwaysUpdateMinDistance(a0, b0, b1, ref min_dist, true);

            var closest_vertex = 0;
            if (UpdateMinDistance(a1, b0, b1, ref min_dist)) { closest_vertex = 1; }
            if (UpdateMinDistance(b0, a0, a1, ref min_dist)) { closest_vertex = 2; }
            if (UpdateMinDistance(b1, a0, a1, ref min_dist)) { closest_vertex = 3; }
            return closest_vertex switch
            {
                0 => (a0, Project(a0, b0, b1)),
                1 => (a1, Project(a1, b0, b1)),
                2 => (Project(b0, a0, a1), b0),
                3 => (Project(b1, a0, a1), b1),
                _ => throw new ApplicationException("Unreached (to suppress Android compiler warning)"),
            };
        }

        // Returns true if every point on edge B=b0b1 is no further than "tolerance"
        // from some point on edge A=a0a1.  Equivalently, returns true if the directed
        // Hausdorff distance from B to A is no more than "tolerance".
        // Requires that tolerance is less than 90 degrees.
        public static bool IsEdgeBNearEdgeA(S2Point a0, S2Point a1, S2Point b0, S2Point b1, S1Angle tolerance)
        {
            Assert.True(tolerance.Radians < S2Constants.M_PI_2);
            Assert.True(tolerance.Radians > 0);
            // The point on edge B=b0b1 furthest from edge A=a0a1 is either b0, b1, or
            // some interior point on B.  If it is an interior point on B, then it must be
            // one of the two points where the great circle containing B (circ(B)) is
            // furthest from the great circle containing A (circ(A)).  At these points,
            // the distance between circ(B) and circ(A) is the angle between the planes
            // containing them.

            S2Point a_ortho = S2PointUtil.RobustCrossProd(a0, a1).Normalized;
            S2Point a_nearest_b0 = Project(b0, a0, a1, a_ortho);
            S2Point a_nearest_b1 = Project(b1, a0, a1, a_ortho);
            // If a_nearest_b0 and a_nearest_b1 have opposite orientation from a0 and a1,
            // we invert a_ortho so that it points in the same direction as a_nearest_b0 x
            // a_nearest_b1.  This helps us handle the case where A and B are oppositely
            // oriented but otherwise might be near each other.  We check orientation and
            // invert rather than computing a_nearest_b0 x a_nearest_b1 because those two
            // points might be equal, and have an unhelpful cross product.
            if (S2Pred.Sign(a_ortho, a_nearest_b0, a_nearest_b1) < 0) a_ortho *= -1;

            // To check if all points on B are within tolerance of A, we first check to
            // see if the endpoints of B are near A.  If they are not, B is not near A.
            var b0_distance = new S1Angle(b0, a_nearest_b0);
            var b1_distance = new S1Angle(b1, a_nearest_b1);
            if (b0_distance > tolerance || b1_distance > tolerance)
                return false;

            // If b0 and b1 are both within tolerance of A, we check to see if the angle
            // between the planes containing B and A is greater than tolerance.  If it is
            // not, no point on B can be further than tolerance from A (recall that we
            // already know that b0 and b1 are close to A, and S2Edges are all shorter
            // than 180 degrees).  The angle between the planes containing circ(A) and
            // circ(B) is the angle between their normal vectors.
            var b_ortho = S2PointUtil.RobustCrossProd(b0, b1).Normalized;
            var planar_angle = new S1Angle(a_ortho, b_ortho);
            if (planar_angle <= tolerance)
                return true;


            // As planar_angle approaches Math.PI, the projection of a_ortho onto the plane
            // of B approaches the null vector, and normalizing it is numerically
            // unstable.  This makes it unreliable or impossible to identify pairs of
            // points where circ(A) is furthest from circ(B).  At this point in the
            // algorithm, this can only occur for two reasons:
            //
            //  1.) b0 and b1 are closest to A at distinct endpoints of A, in which case
            //      the opposite orientation of a_ortho and b_ortho means that A and B are
            //      in opposite hemispheres and hence not close to each other.
            //
            //  2.) b0 and b1 are closest to A at the same endpoint of A, in which case
            //      the orientation of a_ortho was chosen arbitrarily to be that of a0
            //      cross a1.  B must be shorter than 2*tolerance and all points in B are
            //      close to one endpoint of A, and hence to A.
            //
            // The logic applies when planar_angle is robustly greater than S2Constants.M_PI_2, but
            // may be more computationally expensive than the logic beyond, so we choose a
            // value close to Math.PI.
            if (planar_angle >= S1Angle.FromRadians(Math.PI - 0.01))
            {
                return (new S1Angle(b0, a0) < new S1Angle(b0, a1)) ==
                    (new S1Angle(b1, a0) < new S1Angle(b1, a1));
            }

            // Finally, if either of the two points on circ(B) where circ(B) is furthest
            // from circ(A) lie on edge B, edge B is not near edge A.
            //
            // The normalized projection of a_ortho onto the plane of circ(B) is one of
            // the two points along circ(B) where it is furthest from circ(A).  The other
            // is -1 times the normalized projection.
            S2Point furthest = (a_ortho - a_ortho.DotProd(b_ortho) * b_ortho).Normalized;
            Assert.True(furthest.IsUnitLength);
            S2Point furthest_inv = -1 * furthest;

            // A point p lies on B if you can proceed from b_ortho to b0 to p to b1 and
            // back to b_ortho without ever turning right.  We test this for furthest and
            // furthest_inv, and return true if neither point lies on B.
            return !((S2Pred.Sign(b_ortho, b0, furthest) > 0 &&
                      S2Pred.Sign(furthest, b1, b_ortho) > 0) ||
                     (S2Pred.Sign(b_ortho, b0, furthest_inv) > 0 &&
                      S2Pred.Sign(furthest_inv, b1, b_ortho) > 0));
        }

        // Returns the maximum error in the result of UpdateMinInteriorDistance,
        // assuming that all input points are normalized to within the bounds
        // guaranteed by S2Point.Normalize().  The error can be added or subtracted
        // from an S1ChordAngle "x" using x.PlusError(error).
        public static double GetUpdateMinInteriorDistanceMaxError(S1ChordAngle dist)
        {
            // If a point is more than 90 degrees from an edge, then the minimum
            // distance is always to one of the endpoints, not to the edge interior.
            if (dist >= S1ChordAngle.Right) return 0.0;

            // This bound includes all source of error, assuming that the input points
            // are normalized to within the bounds guaranteed to S2Point.Normalize().
            // "a" and "b" are components of chord length that are perpendicular and
            // parallel to the plane containing the edge respectively.
            double b = Math.Min(1.0, 0.5 * dist.Length2);
            double a = Math.Sqrt(b * (2 - b));
            return ((2.5 + 2 * Math.Sqrt(3) + 8.5 * a) * a +
                    (2 + 2 * Math.Sqrt(3) / 3 + 6.5 * (1 - b)) * b +
                    (23 + 16 / Math.Sqrt(3)) * S2Constants.DoubleEpsilon) * S2Constants.DoubleEpsilon;
        }

        // If the minimum distance from X to AB is attained at an interior point of AB
        // (i.e., not an endpoint), and that distance is less than "min_dist" or
        // "always_update" is true, then update "min_dist" and return true.  Otherwise
        // return false.
        //
        // The "Always" in the function name refers to the template argument, i.e.
        // AlwaysUpdateMinInteriorDistance<true> always updates the given distance,
        // while AlwaysUpdateMinInteriorDistance<false> does not.  This optimization
        // increases the speed of GetDistance() by about 10% without creating code
        // duplication.
        private static bool AlwaysUpdateMinInteriorDistance(S2Point x, S2Point a, S2Point b, double xa2, double xb2, ref S1ChordAngle min_dist, bool always_update)
        {
            Assert.True(x.IsUnitLength && a.IsUnitLength && b.IsUnitLength);
            Assert.True(xa2 == (x - a).Norm2);
            Assert.True(xb2 == (x - b).Norm2);

            // The closest point on AB could either be one of the two vertices (the
            // "vertex case") or in the interior (the "interior case").  Let C = A x B.
            // If X is in the spherical wedge extending from A to B around the axis
            // through C, then we are in the interior case.  Otherwise we are in the
            // vertex case.
            //
            // Check whether we might be in the interior case.  For this to be true, XAB
            // and XBA must both be acute angles.  Checking this condition exactly is
            // expensive, so instead we consider the planar triangle ABX (which passes
            // through the sphere's interior).  The planar angles XAB and XBA are always
            // less than the corresponding spherical angles, so if we are in the
            // interior case then both of these angles must be acute.
            //
            // We check this by computing the squared edge lengths of the planar
            // triangle ABX, and testing acuteness using the law of cosines:
            //
            //             Math.Max(XA^2, XB^2) < Math.Min(XA^2, XB^2) + AB^2
            //
            if (Math.Max(xa2, xb2) >= Math.Min(xa2, xb2) + (a - b).Norm2)
            {
                return false;
            }
            // The minimum distance might be to a point on the edge interior.  Let R
            // be closest point to X that lies on the great circle through AB.  Rather
            // than computing the geodesic distance along the surface of the sphere,
            // instead we compute the "chord length" through the sphere's interior.
            // If the squared chord length exceeds min_dist.Length2 then we can
            // return "false" immediately.
            //
            // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
            // is the point X projected onto the plane through the great circle AB.
            // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
            // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
            // is faster and the corresponding distance on the Earth's surface is
            // accurate to within 1% for distances up to about 1800km.
            S2Point c = S2PointUtil.RobustCrossProd(a, b);
            double c2 = c.Norm2;
            double x_dot_c = x.DotProd(c);
            double x_dot_c2 = x_dot_c * x_dot_c;
            if (!always_update && x_dot_c2 > c2 * min_dist.Length2)
            {
                // The closest point on the great circle AB is too far away.  We need to
                // test this using ">" rather than ">=" because the actual minimum bound
                // on the distance is (x_dot_c2 / c2), which can be rounded differently
                // than the (more efficient) multiplicative test above.
                return false;
            }
            // Otherwise we do the exact, more expensive test for the interior case.
            // This test is very likely to succeed because of the conservative planar
            // test we did initially.
            S2Point cx = c.CrossProd(x);
            if (a.DotProd(cx) >= 0 || b.DotProd(cx) <= 0)
            {
                return false;
            }
            // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
            // This calculation has good accuracy for all chord lengths since it
            // is based on both the dot product and cross product (rather than
            // deriving one from the other).  However, note that the chord length
            // representation itself loses accuracy as the angle approaches Pi.
            double qr = 1 - Math.Sqrt(cx.Norm2 / c2);
            double dist2 = (x_dot_c2 / c2) + (qr * qr);
            if (!always_update && dist2 >= min_dist.Length2)
            {
                return false;
            }
            min_dist = S1ChordAngle.FromLength2(dist2);
            return true;
        }

        // This function computes the distance from a point X to a line segment AB.
        // If the distance is less than "min_dist" or "always_update" is true, it
        // updates "min_dist" and returns true.  Otherwise it returns false.
        //
        // The "Always" in the function name refers to the template argument, i.e.
        // AlwaysUpdateMinDistance<true> always updates the given distance, while
        // AlwaysUpdateMinDistance<false> does not.  This optimization increases the
        // speed of GetDistance() by about 10% without creating code duplication.
        private static bool AlwaysUpdateMinDistance(S2Point x, S2Point a, S2Point b, ref S1ChordAngle min_dist, bool always_update)
        {
            Assert.True(x.IsUnitLength && a.IsUnitLength && b.IsUnitLength);

            double xa2 = (x - a).Norm2, xb2 = (x - b).Norm2;
            if (AlwaysUpdateMinInteriorDistance(x, a, b, xa2, xb2, ref min_dist, always_update))
            {
                return true;  // Minimum distance is attained along the edge interior.
            }
            // Otherwise the minimum distance is to one of the endpoints.
            double dist2 = Math.Min(xa2, xb2);
            if (!always_update && dist2 >= min_dist.Length2)
            {
                return false;
            }
            min_dist = S1ChordAngle.FromLength2(dist2);
            return true;
        }
    }
}
