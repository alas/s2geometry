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

namespace S2Geometry;

public static partial class S2
{
    // kIntersectionError is an upper bound on the distance from the intersection
    // point returned by GetIntersection() to the true intersection point.
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

    // kRobustCrossProdError is an upper bound on the angle between the vector
    // returned by RobustCrossProd(a, b) and the true cross product of "a" and "b".
    // Note that cases where "a" and "b" are exactly proportional but not equal
    // (e.g. a = -b or a = (1+epsilon)*b) are handled using symbolic perturbations
    // in order to ensure that the result is non-zero and consistent with S2::Sign.
    public const double kRobustCrossProdError = 6 * S2Pred.DBL_ERR;
    public static readonly S1Angle kRobustCrossProdErrorS1Angle = S1Angle.FromRadians(kRobustCrossProdError);

    // kRobustCrossProdError can be set somewhat arbitrarily because the algorithm
    // uses more precision as needed in order to achieve the specified error.  The
    // only strict requirement is that kRobustCrossProdError >= DBL_ERR, since
    // this is the minimum error even when using exact arithmetic.  We set the
    // error somewhat larger than this so that virtually all cases can be handled
    // using ordinary double-precision arithmetic.
    private static readonly bool StaticAssert01 = (kRobustCrossProdError == 6 * S2Pred.DBL_ERR) ? true : throw new Exception("update comment");

    // kIntersectionError can also be set somewhat arbitrarily (see above) except
    // that in this case the error using exact arithmetic is up to 2 * DBL_ERR,
    // and the error limit is set to 8 * DBL_ERR so that virtually all cases can
    // be handled using ordinary double-precision arithmetic.
    private static readonly bool StaticAssert02 = (kIntersectionError == 8 * S2Pred.DBL_ERR) ? true : throw new Exception("update comment");

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

    // Returns true if the angle ABC contains its vertex B.  Containment is
    // defined such that if several polygons tile the region around a vertex, then
    // exactly one of those polygons contains that vertex.  Returns false for
    // degenerate angles of the form ABA.
    //
    // Note that this method is not sufficient to determine vertex containment in
    // polygons with duplicate vertices (such as the polygon ABCADE).  Use
    // S2ContainsVertexQuery for such polygons.  S2::AngleContainsVertex(a, b, c)
    // is equivalent to using S2ContainsVertexQuery as follows:
    //
    //    S2ContainsVertexQuery query(b);
    //    query.AddEdge(a, -1);  // incoming
    //    query.AddEdge(c, 1);   // outgoing
    //    return query.ContainsSign() > 0;
    //
    // Useful properties of AngleContainsVertex:
    //
    //  (1) AngleContainsVertex(a,b,a) == false
    //  (2) AngleContainsVertex(a,b,c) == !AngleContainsVertex(c,b,a) unless a == c
    //  (3) Given vertices v_1 ... v_k ordered cyclically CCW around vertex b,
    //      AngleContainsVertex(v_{i+1}, b, v_i) is true for exactly one value of i.
    //
    // REQUIRES: a != b && b != c
    public static bool AngleContainsVertex(S2Point a, S2Point b, S2Point c)
    {
        // A loop with consecutive vertices A, B, C contains vertex B if and only if
        // the fixed vector R = S2::RefDir(B) is contained by the wedge ABC.  The
        // wedge is closed at A and open at C, i.e. the point B is inside the loop
        // if A = R but not if C = R.
        //
        // Note that the test below is written so as to get correct results when the
        // angle ABC is degenerate.  If A = C or C = R it returns false, and
        // otherwise if A = R it returns true.
        Assert.True(a != b && b != c);
        return !S2Pred.OrderedCCW(S2.RefDir(b), c, a, b);
    }

    // Given two edges AB and CD where at least two vertices are identical
    // (i.e. CrossingSign(a,b,c,d) == 0), this function defines whether the
    // two edges "cross" in such a way that point-in-polygon containment tests can
    // be implemented by counting the number of edge crossings.  The basic rule is
    // that a "crossing" occurs if AB is encountered after CD during a CCW sweep
    // around the shared vertex starting from a fixed reference point.
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
        if (a == c) return (b == d) || S2Pred.OrderedCCW(S2.RefDir(a), d, b, a);
        if (b == d) return S2Pred.OrderedCCW(S2.RefDir(b), c, a, b);

        if (a == d) return (b == c) || S2Pred.OrderedCCW(S2.RefDir(a), c, b, a);
        if (b == c) return S2Pred.OrderedCCW(S2.RefDir(b), d, a, b);

        throw new ArgumentException("VertexCrossing called with 4 distinct vertices");
    }

    // Like VertexCrossing() but returns -1 if AB crosses CD from left to right,
    // +1 if AB crosses CD from right to left, and 0 otherwise.  This implies that
    // if CD bounds some region according to the "interior is on the left" rule,
    // this function returns -1 when AB exits the region and +1 when AB enters.
    //
    // This is a helper method that allows computing the change in winding number
    // from point A to point B by summing the signed edge crossings of AB with the
    // edges of the loop(s) used to define the winding number.
    //
    // Useful properties of SignedVertexCrossing (SVC):
    //
    //  (1) SVC(a,a,c,d) == SVC(a,b,c,c) == 0
    //  (2) SVC(a,b,a,b) == +1
    //  (3) SVC(a,b,b,a) == -1
    //  (6) SVC(a,b,c,d) == -SVC(a,b,d,c) == -SVC(b,a,c,d) == SVC(b,a,d,c)
    //  (3) If exactly one of a,b equals one of c,d, then exactly one of
    //      SVC(a,b,c,d) and SVC(c,d,a,b) is non-zero
    //
    // It is an error to call this method with 4 distinct vertices.
    public static int SignedVertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d)
    {
        if (a == b || c == d) return 0;

        // See VertexCrossing.  The sign of the crossing is +1 if both edges are
        // outgoing or both edges are incoming with respect to the common vertex
        // and -1 otherwise.
        if (a == c)
        {
            return ((b == d) || S2Pred.OrderedCCW(S2.RefDir(a), d, b, a)) ? 1 : 0;
        }
        if (b == d) return S2Pred.OrderedCCW(S2.RefDir(b), c, a, b) ? 1 : 0;

        if (a == d)
        {
            return ((b == c) || S2Pred.OrderedCCW(S2.RefDir(a), c, b, a)) ? -1 : 0;
        }
        if (b == c) return S2Pred.OrderedCCW(S2.RefDir(b), d, a, b) ? -1 : 0;

        System.Diagnostics.Debug.WriteLine("(DFATAL) " + "SignedVertexCrossing called with 4 distinct vertices");
        return 0;
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
    public static S2Point GetIntersection(S2Point a0, S2Point a1, S2Point b0, S2Point b1, Dictionary<Internal.IntersectionMethod, int>? intersectionMethodTally)
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
        // (only available on some platforms) uses 80-bit precision and is about
        // twice as slow.  The exact arithmetic version is about 100x slower.
        //
        // So our strategy is to first call GetIntersectionStable() in double
        // precision; if that doesn't work and this platform supports "double",
        // then we try again in "double"; if that doesn't work then we fall
        // back to exact arithmetic.

        const bool kUseSimpleMethod = false;
        Internal.IntersectionMethod method;
        if (kUseSimpleMethod && GetIntersectionSimple(a0, a1, b0, b1, out var result))
        {
            method = Internal.IntersectionMethod.SIMPLE;
        }
        else if (kUseSimpleMethod && S2Pred.kHasLongDouble && GetIntersectionSimpleLD(a0, a1, b0, b1, out result))
        {
            method = Internal.IntersectionMethod.SIMPLE_LD;
        }
        else if (GetIntersectionStable(a0, a1, b0, b1, out result))
        {
            method = Internal.IntersectionMethod.STABLE;
        }
        else if (S2Pred.kHasLongDouble && GetIntersectionStableLD(a0, a1, b0, b1, out result))
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

        // Make sure that the intersection point lies on both edges.
        Assert.True(ApproximatelyOrdered(a0, result, a1, kIntersectionError));
        Assert.True(ApproximatelyOrdered(b0, result, b1, kIntersectionError));

        return result;
    }

    // The following functions are not part of the public API.  Currently they are
    // only used internally for testing purposes.
    public static class Internal
    {
        // The maximum error in the method above.
        // NOTE(Alas): The method is "ExactCrossProd"
        public static readonly S1Angle kExactCrossProdError = S1Angle.FromRadians(S2Pred.DBL_ERR);

        // The maximum error in the method above.
        // NOTE(Alas): The method is "GetIntersectionExact"
        public static readonly S1Angle kIntersectionExactError = S1Angle.FromRadians(2 * S2Pred.DBL_ERR);

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

        // Evaluates the cross product of unit-length vectors "a" and "b" in a
        // numerically stable way, returning true if the error in the result is
        // guaranteed to be at most kRobustCrossProdError.
        public static bool GetStableCrossProd(S2Point a, S2Point b, out S2Point result)
        {
            // We compute the cross product (a - b) x (a + b).  Mathematically this is
            // exactly twice the cross product of "a" and "b", but it has the numerical
            // advantage that (a - b) and (a + b) are nearly perpendicular (since "a" and
            // "b" are unit length).  This yields a result that is nearly orthogonal to
            // both "a" and "b" even if these two values differ only very slightly.
            //
            // The maximum directional error in radians when this calculation is done in
            // precision T (where T is a floating-point type) is:
            //
            //   (1 + 2 * sqrt(3) + 32 * sqrt(3) * DBL_ERR / ||N||) * T_ERR
            //
            // where ||N|| is the norm of the result.  To keep this error to at most
            // kRobustCrossProdError, assuming this value is much less than 1, we need
            //
            //   (1 + 2 * sqrt(3) + 32 * sqrt(3) * DBL_ERR / ||N||) * T_ERR <= kErr
            //
            //   ||N|| >= 32 * sqrt(3) * DBL_ERR / (kErr / T_ERR - (1 + 2 * sqrt(3)))
            //
            // From this you can see that in order for this calculation to ever succeed in
            // double precision, we must have kErr > (1 + 2 * sqrt(3)) * DBL_ERR, which is
            // about 4.46 * DBL_ERR.  We actually set kRobustCrossProdError == 6 * DBL_ERR
            // (== 3 * DBL_EPSILON) in order to minimize the number of cases where higher
            // precision is needed; in particular, higher precision is only necessary when
            // "a" and "b" are closer than about 18 * DBL_ERR == 9 * DBL_EPSILON.
            // (80-bit precision can handle inputs as close as 2.5 * LDBL_EPSILON.)
            var T_ERR = S2Pred.DBL_ERR;
            var kMinNorm =
                (32 * S2Pred.kSqrt3 * S2Pred.DBL_ERR) /
                (kRobustCrossProdError / T_ERR - (1 + 2 * S2Pred.kSqrt3));

            result = (a - b).CrossProd(a + b);
            return result.Norm2() >= kMinNorm * kMinNorm;
        }

        // Explicitly instantiate this function so that we can use it in tests without
        // putting its definition in a header file.
        //public bool GetStableCrossProd<double>(S2Point, S2Point, S2Point);
        //public bool GetStableCrossProd<long double>(Vector3_ld, Vector3_ld, Vector3_ld);

        // Returns the cross product of two points computed using exact arithmetic and
        // then symbolic perturbations if necessary, rounded to double-precision and
        // scaled so that the result can be normalized to an S2Point with loss of
        // precision due to floating-point underflow.
        //
        // REQUIRES: a != b (this case should be handled before calling this function)
        public static S2Point ExactCrossProd(S2Point a, S2Point b)
        {
            Assert.True(a != b);
            S2Point result_xf = a.ToExact().CrossProd(b.ToExact());
            if (!S2Pred.IsZero(result_xf))
            {
                return NormalizableFromExact(result_xf);
            }
            // SymbolicCrossProd() requires that a < b.
            if (a < b)
            {
                return EnsureNormalizable(SymbolicCrossProd(a, b));
            }
            else
            {
                return -EnsureNormalizable(SymbolicCrossProd(b, a));
            }
        }

        // Returns the cross product of two points using symbolic perturbations, rounded
        // to double-precision and scaled so that the result can be normalized to an
        // S2Point with loss of precision due to floating-point underflow.
        //
        // REQUIRES: a != b
        // REQUIRES: a and b are linearly dependent
        public static S2Point SymbolicCrossProd(S2Point a, S2Point b)
        {
            Assert.True(a != b);
            // SymbolicCrossProdSorted() requires that a < b.
            if (a < b)
            {
                return EnsureNormalizable(SymbolicCrossProdSorted(a, b));
            }
            else
            {
                return -EnsureNormalizable(SymbolicCrossProdSorted(b, a));
            }
        }

        // Compute the intersection point of (a0, a1) and (b0, b1) using exact
        // arithmetic.  Note that the result is not exact because it is rounded to
        // double precision.
        //
        // Returns the intersection point of two edges computed using exact arithmetic
        // and rounded to the nearest representable S2Point.
        public static S2Point GetIntersectionExact(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
        {
            // Since we are using exact arithmetic, we don't need to worry about
            // numerical stability.
            var a_norm_xf = a0.ToExact().CrossProd(a1.ToExact());
            var b_norm_xf = b0.ToExact().CrossProd(b1.ToExact());
            var x_xf = a_norm_xf.CrossProd(b_norm_xf);

            // The final Normalize() call is done in double precision, which creates a
            // directional error of up to 2 * DBL_ERR.  (NormalizableFromExact() and
            // Normalize() each contribute up to DBL_ERR of directional error.)
            if (!S2Pred.IsZero(x_xf))
            {
                // Make sure that we return the intersection point rather than its antipode.
                return S2Pred.Sign(a0, a1, b1) * ToS2Point(x_xf);
            }

            // The two edges are exactly collinear, but we still consider them to be
            // "crossing" because of simulation of simplicity.  The most principled way to
            // handle this situation is to use symbolic perturbations, similar to what
            // S2::RobustCrossProd and s2pred::Sign do.  This is certainly possible, but
            // it turns out that there are approximately 18 cases to consider (compared to
            // the 4 cases for RobustCrossProd and 13 for s2pred::Sign).
            //
            // For now we use a heuristic that simply chooses a plausible intersection
            // point.  Out of the four endpoints, exactly two lie in the interior of the
            // other edge.  Of those two we return the one that is lexicographically
            // smallest.
            var a_norm = ToS2Point(a_norm_xf);
            var b_norm = ToS2Point(b_norm_xf);
            if (a_norm == new S2Point(0, 0, 0)) a_norm = SymbolicCrossProd(a0, a1);
            if (b_norm == new S2Point(0, 0, 0)) b_norm = SymbolicCrossProd(b0, b1);

            S2Point x=new(10, 10, 10);  // Greater than any valid S2Point.
            if (S2Pred.OrderedCCW(b0, a0, b1, b_norm) && a0 < x) x = a0;
            if (S2Pred.OrderedCCW(b0, a1, b1, b_norm) && a1 < x) x = a1;
            if (S2Pred.OrderedCCW(a0, b0, a1, a_norm) && b0 < x) x = b0;
            if (S2Pred.OrderedCCW(a0, b1, a1, a_norm) && b1 < x) x = b1;

            Assert.True(x.IsUnitLength());
            return x;
        }

        // The following field is used exclusively by s2edge_crossings_test.cc to
        // measure how often each intersection method is used by GetIntersection().
        // If non-nullptr, then it points to an array of integers indexed by an
        // IntersectionMethod enum value.  Each call to GetIntersection() increments
        // the array entry corresponding to the intersection method that was used.
        //extern int* intersection_method_tally_;

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

        // The following classes are used as template arguments to S2EdgeCrosserBase in
        // order to create two versions, namely S2EdgeCrosser itself (which takes
        // (const S2Point *) arguments and requires points to be stored persistently in
        // memory) and S2CopyingEdgeCrosser (which takes (const S2Point&) arguments and
        // makes its own copies of all points).
        //
        // These classes are intended to be used like pointers, i.e. operator* and
        // operator-> return the S2Point value.  They also define an implicit
        // conversion operator that returns the underlying representation as either a
        // const pointer or const reference.

        //class S2Point_PointerRep {...};

        //class S2Point_ValueRep {...};
    }

    // Returns a vector whose direction is guaranteed to be very close to the exact
    // mathematical cross product of the given unit-length vectors "a" and "b", but
    // whose magnitude is arbitrary.  Unlike a.CrossProd(b), this statement is true
    // even when "a" and "b" are very nearly parallel (i.e., a ~= b or a ~= -b).
    // Specifically, the direction of the result vector differs from the exact cross
    // product by at most kRobustCrossProdError radians (see below).
    //
    // When a == -b exactly, the result is consistent with the symbolic perturbation
    // model used by S2::Sign (see s2predicates.h).  In other words, even antipodal
    // point pairs have a consistent and well-defined edge between them.  (In fact
    // this is true for any pair of distinct points whose vectors are parallel.)
    //
    // When a == b exactly, an arbitrary vector orthogonal to "a" is returned.
    // [From a strict mathematical viewpoint it would be better to return (0, 0, 0),
    // but this behavior helps to avoid special cases in client code.]
    //
    // This function has the following properties (RCP == RobustCrossProd):
    //
    //   (1) RCP(a, b) != 0 for all a, b
    //   (2) RCP(b, a) == -RCP(a, b) unless a == b
    //   (3) RCP(-a, b) == -RCP(a, b) unless a and b are exactly proportional
    //   (4) RCP(a, -b) == -RCP(a, b) unless a and b are exactly proportional
    //
    // Note that if you want the result to be unit-length, you must call Normalize()
    // explicitly.  (The result is always scaled such that Normalize() can be called
    // without precision loss due to floating-point underflow.)
    public static S2Point RobustCrossProd(S2Point a, S2Point b)
    {
        Assert.True(a.IsUnitLength());
        Assert.True(b.IsUnitLength());

        // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
        // approaches zero.  This leads to situations where a.CrossProd(b) is not
        // very orthogonal to "a" and/or "b".  To solve this problem robustly requires
        // falling back to extended precision, arbitrary precision, and even symbolic
        // perturbations to handle the case when "a" and "b" are exactly
        // proportional, e.g. a == -b (see s2predicates.cc for details).
        S2Point result;
        if (Internal.GetStableCrossProd(a, b, out result))
        {
            return result;
        }
        // Handle the (a == b) case now, before doing expensive arithmetic.  The only
        // result that makes sense mathematically is to return zero, but it turns out
        // to reduce the number of special cases in client code if we instead return
        // an arbitrary orthogonal vector.
        if (a == b)
        {
            return Ortho(a);
        }
        // Next we try using "long double" precision (if available).
        S2Point result_ld;
        if (S2Pred.kHasLongDouble && Internal.GetStableCrossProd(a.ToLD(), b.ToLD(), out result_ld))
        {
            return result_ld;
        }
        // Otherwise we fall back to exact arithmetic, then symbolic perturbations.
        return Internal.ExactCrossProd(a, b);
    }

    // Returns the cross product of "a" and "b" after symbolic perturbations.
    // (These perturbations only affect the result if "a" and "b" are exactly
    // collinear, e.g. if a == -b or a == (1+eps) * b.)  The result may not be
    // normalizable (i.e., EnsureNormalizable() should be called on the result).
    public static S2Point SymbolicCrossProdSorted(S2Point a, S2Point b)
    {
        Assert.True(a < b);
        Assert.True(S2Pred.IsZero(a.ToExact().CrossProd(b.ToExact())));

        // The following code uses the same symbolic perturbation model as S2::Sign.
        // The particular sequence of tests below was obtained using Mathematica
        // (although it would be easy to do it by hand for this simple case).
        //
        // Just like the function SymbolicallyPerturbedSign() in s2predicates.cc,
        // every input coordinate x[i] is assigned a symbolic perturbation dx[i].  We
        // then compute the cross product
        //
        //     (a + da).CrossProd(b + db) .
        //
        // The result is a polynomial in the perturbation symbols.  For example if we
        // did this in one dimension, the result would be
        //
        //     a * b + b * da + a * db + da * db
        //
        // where "a" and "b" have numerical values and "da" and "db" are symbols.
        // In 3 dimensions the result is similar except that the coefficients are
        // 3-vectors rather than scalars.
        //
        // Every possible S2Point has its own symbolic perturbation in each coordinate
        // (i.e., there are about 3 * 2**192 symbols).  The magnitudes of the
        // perturbations are chosen such that if x < y lexicographically, the
        // perturbations for "y" are much smaller than the perturbations for "x".
        // Similarly, the perturbations for the coordinates of a given point x are
        // chosen such that dx[0] is much smaller than dx[1] which is much smaller
        // than dx[2].  Putting this together with fact the inputs to this function
        // have been sorted so that a < b lexicographically, this tells us that
        //
        //     da[2] > da[1] > da[0] > db[2] > db[1] > db[0]
        //
        // where each perturbation is so much smaller than the previous one that we
        // don't even need to consider it unless the coefficients of all previous
        // perturbations are zero.  In fact, each succeeding perturbation is so small
        // that we don't need to consider it unless the coefficient of all products of
        // the previous perturbations are zero.  For example, we don't need to
        // consider the coefficient of db[1] unless the coefficient of db[2]*da[0] is
        // zero.
        //
        // The follow code simply enumerates the coefficients of the perturbations
        // (and products of perturbations) that appear in the cross product above, in
        // order of decreasing perturbation magnitude.  The first non-zero
        // coefficient determines the result.  The easiest way to enumerate the
        // coefficients in the correct order is to pretend that each perturbation is
        // some tiny value "eps" raised to a power of two:
        //
        // eps**    1      2      4      8     16     32
        //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]
        //
        // Essentially we can then just count in binary and test the corresponding
        // subset of perturbations at each step.  So for example, we must test the
        // coefficient of db[2]*da[0] before db[1] because eps**12 > eps**16.

        if (b[0] != 0 || b[1] != 0)
        {           // da[2]
            return new S2Point(-b[1], b[0], 0);
        }
        if (b[2] != 0)
        {                        // da[1]
            return new S2Point(b[2], 0, 0);         // Note that b[0] == 0.
        }

        // None of the remaining cases can occur in practice, because we can only get
        // to this point if b = (0, 0, 0).  Nevertheless, even (0, 0, 0) has a
        // well-defined direction under the symbolic perturbation model.
        Assert.True(b[1] == 0 && b[2] == 0);        // da[0] coefficients (always zero)

        if (a[0] != 0 || a[1] != 0)
        {          // db[2]
            return new S2Point(a[1], -a[0], 0);
        }

        // The following coefficient is always non-zero, so we can stop here.
        //
        // It may seem strange that we are returning (1, 0, 0) as the cross product
        // without even looking at the sign of a[2].  (Wouldn't you expect
        // (0, 0, -1) x (0, 0, 0) and (0, 0, 1) x (0, 0, 0) to point in opposite
        // directions?)  It's worth pointing out that in this function there is *no
        // relationship whatsoever* between the vectors "a" and "-a", because the
        // perturbations applied to these vectors may be entirely different.  This is
        // why the identity "RobustCrossProd(-a, b) == -RobustCrossProd(a, b)" does
        // not hold whenever "a" and "b" are linearly dependent (i.e., proportional).
        // [As it happens the two cross products above actually do point in opposite
        // directions, but for example (1, 1, 1) x (2, 2, 2) = (-2, 2, 0) and
        // (-1, -1, -1) x (2, 2, 2) = (-2, 2, 0) do not.]
        return new S2Point(1, 0, 0);                   // db[2] * da[1]
    }

    // Returns true if the given vector's magnitude is large enough such that the
    // angle to another vector of the same magnitude can be measured using Angle()
    // without loss of precision due to floating-point underflow.  (This requirement
    // is also sufficient to ensure that Normalize() can be called without risk of
    // precision loss.)
    public static bool IsNormalizable(S2Point p)
    {
        // Let ab = RobustCrossProd(a, b) and cd = RobustCrossProd(cd).  In order for
        // ab.Angle(cd) to not lose precision, the squared magnitudes of ab and cd
        // must each be at least 2**-484.  This ensures that the sum of the squared
        // magnitudes of ab.CrossProd(cd) and ab.DotProd(cd) is at least 2**-968,
        // which ensures that any denormalized terms in these two calculations do
        // not affect the accuracy of the result (since all denormalized numbers are
        // smaller than 2**-1022, which is less than DBL_ERR * 2**-968).
        //
        // The fastest way to ensure this is to test whether the largest component of
        // the result has a magnitude of at least 2**-242.
        return Math.Max(Math.Abs(p[0]), Math.Max(Math.Abs(p[1]), Math.Abs(p[2]))) >= MathUtils.Ldexp(1, -242);
    }

    // Scales a 3-vector as necessary to ensure that the result can be normalized
    // without loss of precision due to floating-point underflow.
    //
    // REQUIRES: p != (0, 0, 0)
    public static S2Point EnsureNormalizable(S2Point p)
    {
        Assert.True(p != new S2Point(0, 0, 0));
        if (!IsNormalizable(p))
        {
            // We can't just scale by a fixed factor because the smallest representable
            // double is 2**-1074, so if we multiplied by 2**(1074 - 242) then the
            // result might be so large that we couldn't square it without overflow.
            //
            // Note that we must scale by a power of two to avoid rounding errors,
            // and that the calculation of "pmax" is free because IsNormalizable()
            // is inline.  The code below scales "p" such that the largest component is
            // in the range [1, 2).
            double p_max = Math.Max(Math.Abs(p[0]), Math.Max(Math.Abs(p[1]), Math.Abs(p[2])));

            // The expression below avoids signed overflow for any value of ilogb().
            return MathUtils.Ldexp(2, -1 - Math.ILogB(p_max)) * p;
        }
        return p;
    }

    // Converts an ExactFloat vector to a double-precision vector, scaling the
    // result as necessary to ensure that the result can be normalized without loss
    // of precision due to floating-point underflow.  (This method doesn't actually
    // call Normalize() since that would create additional error in situations
    // where normalization is not necessary.)
    public static S2Point NormalizableFromExact(S2Point xf)
    {
        S2Point x=new(xf[0].ToDouble(), xf[1].ToDouble(), xf[2].ToDouble());
        if (IsNormalizable(x))
        {
            return x;
        }
        // Scale so that the largest component magnitude is in the range [0.5, 1).
        // Note that the exponents involved could be much smaller than those
        // representable by an IEEE double precision float.
        int exp = ExactFloat.kMinExp - 1;
        for (int i = 0; i < 3; ++i)
        {
            if (xf[i].IsNormal()) exp = Math.Max(exp, xf[i].Exp());
        }
        if (exp < ExactFloat.kMinExp)
        {
            return new S2Point(0, 0, 0);  // The exact result is (0, 0, 0).
        }
        return new S2Point(
            MathUtils.Ldexp(xf[0], -exp).ToDouble(),
            MathUtils.Ldexp(xf[1], -exp).ToDouble(),
            MathUtils.Ldexp(xf[2], -exp).ToDouble());
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
        var length = tmp.Norm();
        if (length != 0)
        {
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
    private static bool GetIntersectionSimple(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
    {
        result = S2Point.Empty;

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
        // the length of the cross product is at least (16 * sqrt(3) + 24) * DBL_ERR
        // (about 6e-15), then the directional error is at most 5 * T_ERR (about
        // 3e-19 when T == "double").  (DBL_ERR appears in the first formula
        // because the inputs are assumed to be normalized in double precision
        // rather than in the given type T.)
        //
        // The third cross product is different because its inputs already have some
        // error.  Letting "result_len" be the length of the cross product, it can
        // be shown that the error is at most
        //
        //   (2 + 2 * sqrt(3) + 12 / result_len) * T_ERR
        //
        // We want this error to be at most kIntersectionError, which is true as
        // long as "result_len" is at least kMinResultLen defined below.

        // On some platforms "double" is the same as "double", and on these
        // platforms this method always returns false (e.g. ARM, Win32).  Rather
        // than testing this directly, instead we look at kMinResultLen since this
        // is a direct measure of whether "double" has sufficient accuracy to
        // be useful.  If kMinResultLen >= 0.5, it means that this method will fail
        // even for edges that meet at an angle of 30 degrees.  (On Intel platforms
        // kMinResultLen corresponds to an intersection angle of about 0.04
        // degrees.)
        Assert.True(kMinResultLen <= 0.5);
        if (kMinResultLen >= 0.5) return false;

        if (RobustNormalWithLength(a0, a1, out var a_norm) >= kMinNormalLength &&
            RobustNormalWithLength(b0, b1, out var b_norm) >= kMinNormalLength &&
            RobustNormalWithLength(a_norm, b_norm, out result) >= kMinResultLen)
        {
            // Make sure that we return the intersection point rather than its antipode.
            result *= (a_norm.DotProd(b1 - b0) < 0) ? -1 : 1;
            return true;
        }
        return false;
    }

    private static readonly double kMinNormalLength = (16 * S2Pred.kSqrt3 + 24) * S2Pred.DBL_ERR;
    private static readonly double kMinResultLen = 12 / (kIntersectionError / S2Pred.TT_ERR - (2 + 2 * S2Pred.kSqrt3));

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
        var x0_dist2 = x0.Norm2();
        var x1_dist2 = x1.Norm2();

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
        // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * DBL_ERR) * T_ERR
        // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * T_ERR
        // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * T_ERR
        error = (((3.5 + 2 * S2Pred.kSqrt3) * a_norm_len + 32 * S2Pred.kSqrt3 * S2Pred.DBL_ERR) * dist + 1.5 * Math.Abs(result)) * S2Pred.TT_ERR;
        return result;
    }

    // Helper function for GetIntersectionStable().  It expects that the edges
    // (a0,a1) and (b0,b1) have been sorted so that the first edge is longer.

    private static bool GetIntersectionStableSorted(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
    {
        Assert.True((a1 - a0).Norm2() >= (b1 - b0).Norm2());
        result = S2Point.Empty;

        // Compute the normal of the plane through (a0, a1) in a stable way.
        S2Point a_norm = (a0 - a1).CrossProd(a0 + a1);
        var a_norm_len = a_norm.Norm();
        var b_len = (b1 - b0).Norm();

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
        //
        // Make sure that we return the intersection point rather than its antipode.
        // It is sufficient to ensure that (b0_dist - b1_dist) is non-negative.
        if (b0_dist < b1_dist)
        {
            b0_dist = -b0_dist;
            b1_dist = -b1_dist;
        }
        var dist_sum = b0_dist - b1_dist;
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
        var x_len2 = x.Norm2();
        if (x_len2 < S2.DoubleMinNorm)
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
    private static bool GetIntersectionStable(S2Point a0, S2Point a1, S2Point b0, S2Point b1, out S2Point result)
    {
        // Sort the two edges so that (a0,a1) is longer, breaking ties in a
        // deterministic way that does not depend on the ordering of the endpoints.
        // This is desirable for two reasons:
        //  - So that the result doesn't change when edges are swapped or reversed.
        //  - It reduces error, since the first edge is used to compute the edge
        //    normal (where a longer edge means less error), and the second edge
        //    is used for interpolation (where a shorter edge means less error).
        double a_len2 = (a1 - a0).Norm2();
        double b_len2 = (b1 - b0).Norm2();
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

    private static S2Point ToS2Point(S2Point xf)
    {
        return NormalizableFromExact(xf).Normalize();
    }

    // Given three points "a", "x", "b", returns true if these three points occur
    // in the given order along the edge (a,b) to within the given tolerance.
    // More precisely, either "x" must be within "tolerance" of "a" or "b", or
    // when "x" is projected onto the great circle through "a" and "b" it must lie
    // along the edge (a,b) (i.e., the shortest path from "a" to "b").
    private static bool ApproximatelyOrdered(S2Point a, S2Point x, S2Point b, double tolerance)
    {
        if ((x - a).Norm2() <= tolerance * tolerance) return true;
        if ((x - b).Norm2() <= tolerance * tolerance) return true;
        return S2Pred.OrderedCCW(a, x, b, S2.RobustCrossProd(a, b).Normalize());
    }
}
