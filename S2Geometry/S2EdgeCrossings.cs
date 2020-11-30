using System;
using System.Collections.Generic;
using System.Linq;

namespace S2Geometry
{
    // Defines functions related to determining whether two geodesic edges cross
    // and for computing intersection points.
    //
    // The predicates CrossingSign(), VertexCrossing(), and EdgeOrVertexCrossing()
    // are robust, meaning that they produce correct, consistent results even in
    // pathological cases.  See s2predicates.h for additional robust predicates.
    //
    // See also S2EdgeCrosser (which efficiently tests an edge against a sequence
    // of other edges) and S2CrossingEdgeQuery (which uses an index to speed up
    // the process).
    //
    // All error bounds in this file are expressed in terms of the maximum
    // rounding error for a floating-point type.  The rounding error is half of
    // the numeric_limits<T>.epsilon() / SPred2.DBL_ERR / S2Constants.DoubleEpsilon value.
    public static class S2EdgeCrossings
    {
        // This function determines whether the edge AB intersects the edge CD.
        // Returns +1 if AB crosses CD at a point that is interior to both edges.
        // Returns  0 if any two vertices from different edges are the same.
        // Returns -1 otherwise.
        //
        // Note that if an edge is degenerate (A == B or C == D), the return value
        // is 0 if two vertices from different edges are the same and -1 otherwise.
        //
        // Properties of CrossingSign:
        //
        //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
        //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
        //  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
        //  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
        //
        // This function implements an exact, consistent perturbation model such
        // that no three points are ever considered to be collinear.  This means
        // that even if you have 4 points A, B, C, D that lie exactly in a line
        // (say, around the equator), C and D will be treated as being slightly to
        // one side or the other of AB.  This is done in a way such that the
        // results are always consistent (see S2Pred.Sign).
        //
        // Note that if you want to check an edge against a collection of other edges,
        // it is much more efficient to use an S2EdgeCrosser (see s2edge_crosser.h).
        public static int CrossingSign(S2Point a, S2Point b, S2Point c, S2Point d)
        {
            var crosser = new S2EdgeCrosser(a, b, c);
            return crosser.CrossingSign(d);
        }

        // Given two edges AB and CD where at least two vertices are identical
        // (i.e. CrossingSign(a,b,c,d) == 0), this function defines whether the
        // two edges "cross" in a such a way that point-in-polygon containment tests
        // can be implemented by counting the number of edge crossings.  The basic
        // rule is that a "crossing" occurs if AB is encountered after CD during a
        // CCW sweep around the shared vertex starting from a fixed reference point.
        //
        // Note that according to this rule, if AB crosses CD then in general CD
        // does not cross AB.  However, this leads to the correct result when
        // counting polygon edge crossings.  For example, suppose that A,B,C are
        // three consecutive vertices of a CCW polygon.  If we now consider the edge
        // crossings of a segment BP as P sweeps around B, the crossing number
        // changes parity exactly when BP crosses BA or BC.
        //
        // Useful properties of VertexCrossing (VC):
        //
        //  (1) VC(a,a,c,d) == VC(a,b,c,c) == false
        //  (2) VC(a,b,a,b) == VC(a,b,b,a) == true
        //  (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
        //  (3) If exactly one of a,b equals one of c,d, then exactly one of
        //      VC(a,b,c,d) and VC(c,d,a,b) is true
        //
        // It is an error to call this method with 4 distinct vertices.
        public static bool VertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d)
        {
            // If A == B or C == D there is no intersection.  We need to check this
            // case first in case 3 or more input points are identical.
            if (a == b || c == d) return false;

            // If any other pair of vertices is equal, there is a crossing if and only
            // if OrderedCCW() indicates that the edge AB is further CCW around the
            // shared vertex O (either A or B) than the edge CD, starting from an
            // arbitrary fixed reference point.
            //
            // Optimization: if AB=CD or AB=DC, we can avoid most of the calculations.
            if (a == c) return (b == d) || S2Pred.OrderedCCW(S2PointUtil.Ortho(a), d, b, a);
            if (b == d) return S2Pred.OrderedCCW(S2PointUtil.Ortho(b), c, a, b);

            if (a == d) return (b == c) || S2Pred.OrderedCCW(S2PointUtil.Ortho(a), c, b, a);
            if (b == c) return S2Pred.OrderedCCW(S2PointUtil.Ortho(b), d, a, b);

            throw new ArgumentException("VertexCrossing called with 4 distinct vertices");
        }

        // A convenience function that calls CrossingSign() to handle cases
        // where all four vertices are distinct, and VertexCrossing() to handle
        // cases where two or more vertices are the same.  This defines a crossing
        // function such that point-in-polygon containment tests can be implemented
        // by simply counting edge crossings.
        public static bool EdgeOrVertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d)
        {
            int crossing = CrossingSign(a, b, c, d);
            if (crossing < 0) return false;
            if (crossing > 0) return true;
            return VertexCrossing(a, b, c, d);
        }

        // Given two edges AB and CD such that CrossingSign(A, B, C, D) > 0, returns
        // their intersection point.  Useful properties of GetIntersection (GI):
        //
        //  (1) GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
        //  (2) GI(c,d,a,b) == GI(a,b,c,d)
        //
        // The returned intersection point X is guaranteed to be very close to the
        // true intersection point of AB and CD, even if the edges intersect at a
        // very small angle.  See "kIntersectionError" below for details.
        public static S2Point GetIntersection(S2Point a0, S2Point a1, S2Point b0, S2Point b1, Dictionary<Internal.IntersectionMethod, int> intersectionMethodTally)
        {
            Assert.True(CrossingSign(a0, a1, b0, b1) > 0);

            // It is difficult to compute the intersection point of two edges accurately
            // when the angle between the edges is very small.  Previously we handled
            // this by only guaranteeing that the returned intersection point is within
            // kIntersectionError of each edge.  However, this means that when the edges
            // cross at a very small angle, the computed result may be very far from the
            // true intersection point.
            //
            // Instead this function now guarantees that the result is always within
            // kIntersectionError of the true intersection.  This requires using more
            // sophisticated techniques and in some cases extended precision.
            //
            // Three different techniques are implemented, but only two are used:
            //
            //  - GetIntersectionSimple() computes the intersection point using
            //    numerically stable cross products in "double" precision.
            //
            //  - GetIntersectionStable() computes the intersection point using
            //    projection and interpolation, taking care to minimize cancellation
            //    error.  This method has 2 floating point precision versions.
            //
            //  - GetIntersectionExact() computes the intersection point using exact
            //    arithmetic and converts the final result back to an S2Point.
            //
            // We don't actually use the first method (GetIntersectionSimple) because it
            // turns out that GetIntersectionStable() is twice as fast and also much
            // more accurate (even in double precision).  The "double" version
            // (only available on Intel platforms) uses 80-bit precision and is about
            // twice as slow.  The exact arithmetic version is about 100x slower.
            //
            // So our strategy is to first call GetIntersectionStable() in double
            // precision; if that doesn't work and this platform supports "double",
            // then we try again in "double"; if that doesn't work then we fall
            // back to exact arithmetic.

            const bool kUseSimpleMethod = false;
            const bool kHasLongDouble = false; // S2Pred.DBL_ERR < S2Pred.LD_ERR;
            Internal.IntersectionMethod method;
            if (kUseSimpleMethod && GetIntersectionSimple(a0, a1, b0, b1, out var result))
            {
                method = Internal.IntersectionMethod.SIMPLE;
            }
            else if (kUseSimpleMethod && kHasLongDouble && GetIntersectionSimpleLD(a0, a1, b0, b1, out result))
            {
                method = Internal.IntersectionMethod.SIMPLE_LD;
            }
            else if (GetIntersectionStable(a0, a1, b0, b1, out result))
            {
                method = Internal.IntersectionMethod.STABLE;
            }
            else if (kHasLongDouble && GetIntersectionStableLD(a0, a1, b0, b1, out result))
            {
                method = Internal.IntersectionMethod.STABLE_LD;
            }
            else
            {
                result = Internal.GetIntersectionExact(a0, a1, b0, b1);
                method = Internal.IntersectionMethod.EXACT;
            }
            
            if (intersectionMethodTally != null)
            {
                intersectionMethodTally[method]++;
            }

            // Make sure the intersection point is on the correct side of the sphere.
            // Since all vertices are unit length, and edges are less than 180 degrees,
            // (a0 + a1) and (b0 + b1) both have positive dot product with the
            // intersection point.  We use the sum of all vertices to make sure that the
            // result is unchanged when the edges are swapped or reversed.
            if (result.DotProd((a0 + a1) + (b0 + b1)) < 0) result = -result;

            // Make sure that the intersection point lies on both edges.
            Assert.True(ApproximatelyOrdered(a0, result, a1, kIntersectionError));
            Assert.True(ApproximatelyOrdered(b0, result, b1, kIntersectionError));

            return result;
        }

#pragma warning disable IDE1006 // Estilos de nombres
        // kIntersectionError is an upper bound on the distance from the intersection
        // point returned by GetIntersection() to the true intersection point.
        //
        // kIntersectionError can be set somewhat arbitrarily, because the algorithm
        // uses more precision when necessary in order to achieve the specified error.
        // The only strict requirement is that kIntersectionError >= 2 * DBL_ERR
        // radians.  However, using a larger error tolerance makes the algorithm more
        // efficient because it reduces the number of cases where exact arithmetic is
        // needed.
        public const double kIntersectionError = 8 * S2Pred.DBL_ERR;
        public static readonly S1Angle kIntersectionErrorS1Angle = S1Angle.FromRadians(kIntersectionError);

        // This value can be used as the S2Builder snap_radius() to ensure that edges
        // that have been displaced by up to kIntersectionError are merged back
        // together again.  For example this can happen when geometry is intersected
        // with a set of tiles and then unioned.  It is equal to twice the
        // intersection error because input edges might have been displaced in
        // opposite directions.
        public const double kIntersectionMergeRadius = 2 * kIntersectionError;
        public static readonly S1Angle kIntersectionMergeRadiusS1Angle = S1Angle.FromRadians(kIntersectionMergeRadius);

        // The following functions are not part of the public API.  Currently they are
        // only used internally for testing purposes.
        public static class Internal
        {
            // The maximum error in the method above.
            public static readonly S1Angle kIntersectionExactError = S1Angle.FromRadians(2 * S2Pred.DBL_ERR);
#pragma warning restore IDE1006 // Estilos de nombres

            public static Dictionary<IntersectionMethod, int> GetNewTally()
            {
                return new Dictionary<IntersectionMethod, int>
                {
                    { IntersectionMethod.SIMPLE, 0 },
                    { IntersectionMethod.SIMPLE_LD, 0 },
                    { IntersectionMethod.STABLE, 0 },
                    { IntersectionMethod.STABLE_LD, 0 },
                    { IntersectionMethod.EXACT, 0 },
                };
            }

            public static string GetIntersectionMethodName(IntersectionMethod method)
            {
                return method switch
                {
                    IntersectionMethod.SIMPLE => "Simple",
                    IntersectionMethod.SIMPLE_LD => "Simple_ld",
                    IntersectionMethod.STABLE => "Stable",
                    IntersectionMethod.STABLE_LD => "Stable_ld",
                    IntersectionMethod.EXACT => "Exact",
                    _ => "Unknown",
                };
            }

            // Compute the intersection point of (a0, a1) and (b0, b1) using exact
            // arithmetic.  Note that the result is not exact because it is rounded to
            // double precision.  Also, the intersection point is not guaranteed to have
            // the correct sign (i.e., the return value may need to be negated).
            //
            // Returns the intersection point of two edges computed using exact arithmetic
            // and rounded to the nearest representable S2Point.
            public static S2Point GetIntersectionExact(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
            {
                // Since we are using exact arithmetic, we don't need to worry about
                // numerical stability.
                var a0_xf = a0.ToExact();
                var a1_xf = a1.ToExact();
                var b0_xf = b0.ToExact();
                var b1_xf = b1.ToExact();
                var a_norm_xf = a0_xf.CrossProd(a1_xf);
                var b_norm_xf = b0_xf.CrossProd(b1_xf);
                var x_xf = a_norm_xf.CrossProd(b_norm_xf);

                // The final Normalize() call is done in double precision, which creates a
                // directional error of up to 2 * DBL_ERR.  (ToDouble() and "Normalized" each
                // contribute up to DBL_ERR of directional error.)
                S2Point x = S2Point.FromExact(x_xf);

                if (x == S2Point.Empty)
                {
                    // The two edges are exactly collinear, but we still consider them to be
                    // "crossing" because of simulation of simplicity.  Out of the four
                    // endpoints, exactly two lie in the interior of the other edge.  Of
                    // those two we return the one that is lexicographically smallest.
                    x = new S2Point(10, 10, 10);  // Greater than any valid S2Point
                    S2Point a_norm = S2Point.FromExact(a_norm_xf);
                    S2Point b_norm = S2Point.FromExact(b_norm_xf);
                    if (a_norm == S2Point.Empty || b_norm == S2Point.Empty)
                    {
                        // TODO(ericv): To support antipodal edges properly, we would need to
                        // add an S2Pred.CrossProd() function that computes the cross product
                        // using simulation of simplicity and rounds the result to the nearest
                        // floating-point representation.
                        throw new ApplicationException("Exactly antipodal edges not supported by GetIntersection");
                    }
                    if (S2Pred.OrderedCCW(b0, a0, b1, b_norm) && a0 < x) x = a0;
                    if (S2Pred.OrderedCCW(b0, a1, b1, b_norm) && a1 < x) x = a1;
                    if (S2Pred.OrderedCCW(a0, b0, a1, a_norm) && b0 < x) x = b0;
                    if (S2Pred.OrderedCCW(a0, b1, a1, a_norm) && b1 < x) x = b1;
                }
                Assert.True(x.IsUnitLength);
                return x;
            }

            // The list of intersection methods implemented by GetIntersection().
            public enum IntersectionMethod
            {
                None,
                SIMPLE,
                SIMPLE_LD,
                STABLE,
                STABLE_LD,
                EXACT,
                NUM_METHODS,
            }
        }

        // Computes the cross product of "x" and "y", normalizes it to be unit length,
        // and stores the result in "result".  Also returns the length of the cross
        // product before normalization, which is useful for estimating the amount of
        // error in the result.  For numerical stability, "x" and "y" should both be
        // approximately unit length.
        private static double RobustNormalWithLength(S2Point x, S2Point y, out S2Point result)
        {
            // This computes 2 * (x.CrossProd(y)), but has much better numerical
            // stability when "x" and "y" are unit length.
            var tmp = (x - y).CrossProd(x + y);
            var length = tmp.Norm;
            if (length != 0) {
                result = (1 / length) * tmp;
            }
            else
            {
                result = S2Point.Empty;
            }
            return 0.5 * length;  // Since tmp == 2 * (x.CrossProd(y))
        }

        // If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
        // to within an error of at most kIntersectionError by this function, then set
        // "result" to the intersection point and return true.
        //
        // The intersection point is not guaranteed to have the correct sign
        // (i.e., it may be either "result" or "-result").
        private static bool GetIntersectionSimple(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
        {
            // The code below computes the intersection point as
            //
            //    (a0.CrossProd(a1)).CrossProd(b0.CrossProd(b1))
            //
            // except that it has better numerical stability and also computes a
            // guaranteed error bound.
            //
            // Each cross product is computed as (X-Y).CrossProd(X+Y) using unit-length
            // input vectors, which eliminates most of the cancellation error.  However
            // the error in the direction of the cross product can still become large if
            // the two points are extremely close together.  We can show that as long as
            // the length of the cross product is at least (16 * Math.Sqrt(3) + 24) * DBL_ERR
            // (about 6e-15), then the directional error is at most 5 * T_ERR (about
            // 3e-19 when T == "double").  (DBL_ERR appears in the first formula
            // because the inputs are assumed to be normalized in double precision
            // rather than in the given type T.)
            //
            // The third cross product is different because its inputs already have some
            // error.  Letting "result_len" be the length of the cross product, it can
            // be shown that the error is at most
            //
            //   (2 + 2 * Math.Sqrt(3) + 12 / result_len) * T_ERR
            //
            // We want this error to be at most kIntersectionError, which is true as
            // long as "result_len" is at least kMinResultLen defined below.

            // On some platforms "double" is the same as "double", and on these
            // platforms this method always returns false (e.g. ARM, Win32).  Rather
            // than testing this directly, instead we look at kMinResultLen since this
            // is a direct measure of whether "double" has sufficient accuracy to
            // be useful.  If kMinResultLen > 0.5, it means that this method will fail
            // even for edges that meet at an angle of 30 degrees.  (On Intel platforms
            // kMinResultLen corresponds to an intersection angle of about 0.04
            // degrees.)
            Assert.True(kMinResultLen <= 0.5);

            result = S2Point.Empty;
            return RobustNormalWithLength(a0, a1, out var a_norm) >= kMinNormalLength &&
                   RobustNormalWithLength(b0, b1, out var b_norm) >= kMinNormalLength &&
                   RobustNormalWithLength(a_norm, b_norm, out result) >= kMinResultLen;
        }

        private static readonly double kMinNormalLength = (16 * Math.Sqrt(3) + 24) * S2Pred.DBL_ERR;
        private static readonly double kMinResultLen = 12 / (kIntersectionError / S2Pred.TT_ERR - (2 + 2 * Math.Sqrt(3)));

        private static bool GetIntersectionSimpleLD(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
        {
            if (GetIntersectionSimple(a0.ToLD(), a1.ToLD(), b0.ToLD(), b1.ToLD(), out var result_ld))
            {
                result = S2Point.FromLD(result_ld);
                return true;
            }
            else
            {
                result = S2Point.Empty;
            }
            return false;
        }

        // Given a point X and a vector "a_norm" (not necessarily unit length),
        // compute x.DotProd(a_norm) and return a bound on the error in the result.
        // The remaining parameters allow this dot product to be computed more
        // accurately and efficiently.  They include the length of "a_norm"
        // ("a_norm_len") and the edge endpoints "a0" and "a1".

        private static double GetProjection(S2Point x, S2Point a_norm, double a_norm_len, S2Point a0, S2Point a1, out double error)
        {
            // The error in the dot product is proportional to the lengths of the input
            // vectors, so rather than using "x" itself (a unit-length vector) we use
            // the vectors from "x" to the closer of the two edge endpoints.  This
            // typically reduces the error by a huge factor.
            S2Point x0 = x - a0;
            S2Point x1 = x - a1;
            var x0_dist2 = x0.Norm2;
            var x1_dist2 = x1.Norm2;

            // If both distances are the same, we need to be careful to choose one
            // endpoint deterministically so that the result does not change if the
            // order of the endpoints is reversed.
            double dist, result;
            if (x0_dist2 < x1_dist2 || (x0_dist2 == x1_dist2 && x0 < x1))
            {
                dist = Math.Sqrt(x0_dist2);
                result = x0.DotProd(a_norm);
            }
            else
            {
                dist = Math.Sqrt(x1_dist2);
                result = x1.DotProd(a_norm);
            }
            // This calculation bounds the error from all sources: the computation of
            // the normal, the subtraction of one endpoint, and the dot product itself.
            // (DBL_ERR appears because the input points are assumed to be normalized in
            // double precision rather than in the given type T.)
            //
            // For reference, the bounds that went into this calculation are:
            // ||N'-N|| <= ((1 + 2 * Math.Sqrt(3))||N|| + 32 * Math.Sqrt(3) * DBL_ERR) * T_ERR
            // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * T_ERR
            // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * T_ERR
            error = (((3.5 + 2 * Math.Sqrt(3)) * a_norm_len + 32 * Math.Sqrt(3) * S2Pred.DBL_ERR) * dist + 1.5 * Math.Abs(result)) * S2Pred.TT_ERR;
            return result;
        }

        // Helper function for GetIntersectionStable().  It expects that the edges
        // (a0,a1) and (b0,b1) have been sorted so that the first edge is longer.

        private static bool GetIntersectionStableSorted(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
        {
            Assert.True((a1 - a0).Norm2 >= (b1 - b0).Norm2);
            result = S2Point.Empty;

            // Compute the normal of the plane through (a0, a1) in a stable way.
            S2Point a_norm = (a0 - a1).CrossProd(a0 + a1);
            var a_norm_len = a_norm.Norm;
            var b_len = (b1 - b0).Norm;

            // Compute the projection (i.e., signed distance) of b0 and b1 onto the
            // plane through (a0, a1).  Distances are scaled by the length of a_norm.
            var b0_dist = GetProjection(b0, a_norm, a_norm_len, a0, a1, out var b0_error);
            var b1_dist = GetProjection(b1, a_norm, a_norm_len, a0, a1, out var b1_error);

            // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is
            // |b0_dist - b1_dist|.  Note that b0_dist and b1_dist generally have
            // opposite signs because b0 and b1 are on opposite sides of (a0, a1).  The
            // code below finds the intersection point by interpolating along the edge
            // (b0, b1) to a fractional distance of b0_dist / (b0_dist - b1_dist).
            //
            // It can be shown that the maximum error in the interpolation fraction is
            //
            //     (b0_dist * b1_error - b1_dist * b0_error) /
            //        (dist_sum * (dist_sum - error_sum))
            //
            // We save ourselves some work by scaling the result and the error bound by
            // "dist_sum", since the result is normalized to be unit length anyway.
            var dist_sum = Math.Abs(b0_dist - b1_dist);
            var error_sum = b0_error + b1_error;
            if (dist_sum <= error_sum)
            {
                return false;  // Error is unbounded in this case.
            }
            S2Point x = b0_dist * b1 - b1_dist * b0;
            var error = b_len * Math.Abs(b0_dist * b1_error - b1_dist * b0_error) /
                (dist_sum - error_sum) + 2 * S2Pred.TT_ERR * dist_sum;

            // Finally we normalize the result, compute the corresponding error, and
            // check whether the total error is acceptable.
            var x_len2 = x.Norm2;
            if (x_len2 < S2Constants.DoubleMinNorm)
            {
                // If x.Norm2 is less than the minimum normalized value of T, x_len might
                // lose precision and the result might fail to satisfy S2.IsUnitLength().
                // TODO(ericv): Implement S2.RobustNormalize().
                return false;
            }
            var x_len = Math.Sqrt(x_len2);
            var kMaxError = kIntersectionError;
            if (error > (kMaxError - S2Pred.TT_ERR) * x_len)
            {
                return false;
            }
            result = (1 / x_len) * x;
            return true;
        }

        // Returns whether (a0,a1) is less than (b0,b1) with respect to a total
        // ordering on edges that is invariant under edge reversals.
        private static bool CompareEdges(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
        {
            var a = new[] { a0, a1 }.Min();
            if (b0 >= b1) { var tmp = b0; b0 = b1; b1 = tmp; }
            return a < b0 || (a == b0 && b0 < b1);
        }

        // If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
        // to within an error of at most kIntersectionError by this function, then set
        // "result" to the intersection point and return true.
        //
        // The intersection point is not guaranteed to have the correct sign
        // (i.e., it may be either "result" or "-result").
        private static bool GetIntersectionStable(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
        {
            // Sort the two edges so that (a0,a1) is longer, breaking ties in a
            // deterministic way that does not depend on the ordering of the endpoints.
            // This is desirable for two reasons:
            //  - So that the result doesn't change when edges are swapped or reversed.
            //  - It reduces error, since the first edge is used to compute the edge
            //    normal (where a longer edge means less error), and the second edge
            //    is used for interpolation (where a shorter edge means less error).
            double a_len2 = (a1 - a0).Norm2;
            double b_len2 = (b1 - b0).Norm2;
            if (a_len2 < b_len2 || (a_len2 == b_len2 && CompareEdges(a0, a1, b0, b1)))
            {
                return GetIntersectionStableSorted(b0, b1, a0, a1, out result);
            }
            else
            {
                return GetIntersectionStableSorted(a0, a1, b0, b1, out result);
            }
        }

        private static bool GetIntersectionStableLD(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
        {
            if (GetIntersectionStable(a0.ToLD(), a1.ToLD(), b0.ToLD(), b1.ToLD(), out var result_ld))
            {
                result = S2Point.FromLD(result_ld);
                return true;
            }
            else
            {
                result = S2Point.Empty;
                return false;
            }
        }

#pragma warning disable IDE0051 // Quitar miembros privados no utilizados
        private static S2Point S2PointFromExact(S2Point xf)
#pragma warning restore IDE0051 // Quitar miembros privados no utilizados
        {
            // If all components of "x" have absolute value less than about 1e-154,
            // then x.Norm2 is zero in double precision due to underflow.  Therefore
            // we need to scale "x" by an appropriate power of 2 before the conversion.
            var x = new S2Point(xf[0].ToDouble(), xf[1].ToDouble(), xf[2].ToDouble());
            if (x.Norm2 > 0) return x.Normalized;

            // Scale so that the largest component magnitude is in the range [0.5, 1).
            int exp = S2Constants.kMinExp - 1;
            for (int i = 0; i < 3; ++i)
            {
                if (double.IsNormal(xf[i]))
                    exp = (int)Math.Max(exp, Math.Exp(xf[i]));
            }
            if (exp < S2Constants.kMinExp)
            {
                return S2Point.Empty;
            }
            return new S2Point(
                MathUtils.Ldexp(xf[0], -exp).ToDouble(),
                MathUtils.Ldexp(xf[1], -exp).ToDouble(),
                MathUtils.Ldexp(xf[2], -exp).ToDouble()).Normalized;
        }

        // Given three points "a", "x", "b", returns true if these three points occur
        // in the given order along the edge (a,b) to within the given tolerance.
        // More precisely, either "x" must be within "tolerance" of "a" or "b", or
        // when "x" is projected onto the great circle through "a" and "b" it must lie
        // along the edge (a,b) (i.e., the shortest path from "a" to "b").
        private static bool ApproximatelyOrdered(S2Point a, S2Point x, S2Point b, double tolerance)
        {
            if ((x - a).Norm2 <= tolerance * tolerance) return true;
            if ((x - b).Norm2 <= tolerance * tolerance) return true;
            return S2Pred.OrderedCCW(a, x, b, S2PointUtil.RobustCrossProd(a, b).Normalized);
        }
    }
}
