// This class contains various predicates that are guaranteed to produce
// correct, consistent results.  They are also relatively efficient.  This is
// achieved by computing conservative error bounds and falling back to high
// precision or even exact arithmetic when the result is uncertain.  Such
// predicates are useful in implementing robust algorithms.
// 
// s2edge_crossings.h contains the following exact predicates that test for
// edge crossings.  (Note that usually you should use S2EdgeCrosser, which
// implements them in a much more efficient way.)
//
//   int CrossingSign(const S2Point& a0, const S2Point& a1,
//                    const S2Point& b0, const S2Point& b1);
//
//   bool EdgeOrVertexCrossing(const S2Point& a0, const S2Point& a1,
//                             const S2Point& b0, const S2Point& b1);
//
// It also contains the following functions, which compute their result to
// within a guaranteed tolerance and are consistent with the predicates
// defined here (including using symbolic perturbations when necessary):
//
//   S2Point RobustCrossProd(const S2Point& a, const S2Point& b);
//
//   S2Point GetIntersection(const S2Point& a, const S2Point& b,
//                           const S2Point& c, const S2Point& d);
// 
// TODO(ericv): Add InCircleSign() (the Voronoi/Delaunay predicate).
// (This is trickier than the usual textbook implementations because we want
// to model S2Points as lying exactly on the mathematical unit sphere.)

namespace S2Geometry;

public static class S2Pred
{
    /// <summary>
    /// Returns +1 if the points A, B, C are counterclockwise, -1 if the points
    /// are clockwise, and 0 if any two points are the same.  This function is
    /// essentially like taking the sign of the determinant of ABC, except that
    /// it has additional logic to make sure that the above properties hold even
    /// when the three points are coplanar, and to deal with the limitations of
    /// floating-point arithmetic.
    ///
    /// Sign satisfies the following conditions:
    ///
    ///  (1) Sign(a,b,c) == 0 if and only if a == b, b == c, or c == a
    ///  (2) Sign(b,c,a) == Sign(a,b,c) for all a,b,c
    ///  (3) Sign(c,b,a) == -Sign(a,b,c) for all a,b,c
    ///
    /// In other words:
    ///
    ///  (1) The result is zero if and only if two points are the same.
    ///  (2) Rotating the order of the arguments does not affect the result.
    ///  (3) Exchanging any two arguments inverts the result.
    ///
    /// On the other hand, note that it is not true in general that
    /// Sign(-a,b,c) == -Sign(a,b,c), or any similar identities
    /// involving antipodal points.
    /// </summary>
    public static int Sign(S2Point a, S2Point b, S2Point c)
    {
        // We don't need RobustCrossProd() here because Sign() does its own
        // error estimation and calls ExpensiveSign() if there is any uncertainty
        // about the result.
        return Sign(a, b, c, a.CrossProd(b));
    }

    // Given 4 points on the unit sphere, return true if the edges OA, OB, and
    // OC are encountered in that order while sweeping CCW around the point O.
    // You can think of this as testing whether A <= B <= C with respect to the
    // CCW ordering around O that starts at A, or equivalently, whether B is
    // contained in the range of angles (inclusive) that starts at A and extends
    // CCW to C.  Properties:
    //
    //  (1) If OrderedCCW(a,b,c,o) && OrderedCCW(b,a,c,o), then a == b
    //  (2) If OrderedCCW(a,b,c,o) && OrderedCCW(a,c,b,o), then b == c
    //  (3) If OrderedCCW(a,b,c,o) && OrderedCCW(c,b,a,o), then a == b == c
    //  (4) If a == b or b == c, then OrderedCCW(a,b,c,o) is true
    //  (5) Otherwise if a == c, then OrderedCCW(a,b,c,o) is false
    //
    // REQUIRES: a != o && b != o && c != o
    public static bool OrderedCCW(S2Point a, S2Point b, S2Point c, S2Point o)
    {
        MyDebug.Assert(a != o && b != o && c != o);

        // The last inequality below is ">" rather than ">=" so that we return true
        // if A == B or B == C, and otherwise false if A == C.  Recall that
        // Sign(x,y,z) == -Sign(z,y,x) for all x,y,z.

        var sum = 0;
        if (Sign(b, o, a) >= 0) ++sum;
        if (Sign(c, o, b) >= 0) ++sum;
        if (Sign(a, o, c) > 0) ++sum;
        return sum >= 2;
    }

    // Returns -1, 0, or +1 according to whether AX < BX, A == B, or AX > BX
    // respectively.  Distances are measured with respect to the positions of X,
    // A, and B as though they were reprojected to lie exactly on the surface of
    // the unit sphere.  Furthermore, this method uses symbolic perturbations to
    // ensure that the result is non-zero whenever A != B, even when AX == BX
    // exactly, or even when A and B project to the same point on the sphere.
    // Such results are guaranteed to be self-consistent, i.e. if AB < BC and
    // BC < AC, then AB < AC.
    public static int CompareDistances(S2Point x, S2Point a, S2Point b)
    {
        // We start by comparing distances using dot products (i.e., cosine of the
        // angle), because (1) this is the cheapest technique, and (2) it is valid
        // over the entire range of possible angles.  (We can only use the sin^2
        // technique if both angles are less than 90 degrees or both angles are
        // greater than 90 degrees.)
        int sign = TriageCompareCosDistances(x, a, b);
        if (sign != 0) return sign;

        // Optimization for (a == b) to avoid falling back to exact arithmetic.
        if (a == b) return 0;

        // It is much better numerically to compare distances using cos(angle) if
        // the distances are near 90 degrees and sin^2(angle) if the distances are
        // near 0 or 180 degrees.  We only need to check one of the two angles when
        // making this decision because the fact that the test above failed means
        // that angles "a" and "b" are very close together.
        double cos_ax = a.DotProd(x);
        if (cos_ax > S2.M_SQRT1_2)
        {
            // Angles < 45 degrees.
            sign = CompareSin2Distances(x, a, b);
        }
        else if (cos_ax < -S2.M_SQRT1_2)
        {
            // Angles > 135 degrees.  sin^2(angle) is decreasing in this range.
            sign = -CompareSin2Distances(x, a, b);
        }
        else if (kHasLongDouble)
        {
            // We've already tried double precision, so continue with "double".
            sign = TriageCompareCosDistances(x.ToLD(), a.ToLD(), b.ToLD());
        }
        if (sign != 0) return sign;
        sign = ExactCompareDistances(x.ToExact(), a.ToExact(), b.ToExact());
        if (sign != 0) return sign;
        return SymbolicCompareDistances(x, a, b);
    }

    // Returns -1, 0, or +1 according to whether the distance XY is less than,
    // equal to, or greater than "r" respectively.  Distances are measured with
    // respect the positions of all points as though they are projected to lie
    // exactly on the surface of the unit sphere.
    public static int CompareDistance(S2Point x, S2Point y, S1ChordAngle r)
    {
        // As with CompareDistances(), we start by comparing dot products because
        // the sin^2 method is only valid when the distance XY and the limit "r" are
        // both less than 90 degrees.
        int sign = TriageCompareCosDistance(x, y, r.Length2);
        if (sign != 0) return sign;

        // Optimization for (x == y) to avoid falling back to exact arithmetic.
        if (r.Length2 == 0 && x == y) return 0;

        // Unlike with CompareDistances(), it's not worth using the sin^2 method
        // when the distance limit is near 180 degrees because the S1ChordAngle
        // representation itself has has a rounding error of up to 2e-8 radians for
        // distances near 180 degrees.
        if (r < K45Degrees)
        {
            sign = TriageCompareSin2Distance(x, y, r.Length2);
            if (kHasLongDouble && sign == 0)
            {
                sign = TriageCompareSin2Distance(x.ToLD(), y.ToLD(), r.Length2.ToLD());
            }
        }
        else if (kHasLongDouble)
        {
            sign = TriageCompareCosDistance(x.ToLD(), y.ToLD(), r.Length2.ToLD());
        }
        if (sign != 0) return sign;
        return ExactCompareDistance(x.ToExact(), y.ToExact(), (decimal)r.Length2);
    }

    // Returns -1, 0, or +1 according to whether the distance from the point X to
    // the edge A is less than, equal to, or greater than "r" respectively.
    // Distances are measured with respect the positions of all points as though
    // they were projected to lie exactly on the surface of the unit sphere.
    //
    // REQUIRES: A0 and A1 do not project to antipodal points (e.g., A0 == -A1).
    //           This requires that (A0 != C * A1) for anyant C < 0.
    //
    // NOTE(ericv): All of the predicates defined here could be extended to handle
    // edges consisting of antipodal points by implementing additional symbolic
    // perturbation logic (similar to Sign) in order to rigorously define the
    // direction of such edges.
    public static int CompareEdgeDistance(S2Point x, S2Point a0, S2Point a1, S1ChordAngle r)
    {
        // Check that the edge does not consist of antipodal points.  (This catches
        // the most common case -- the full test is in ExactCompareEdgeDistance.)
        MyDebug.Assert(a0 != -a1);

        int sign = TriageCompareEdgeDistance(x, a0, a1, r.Length2);
        if (sign != 0) return sign;

        // Optimization for the case where the edge is degenerate.
        if (a0 == a1) return CompareDistance(x, a0, r);

        if (kHasLongDouble)
        {
            sign = TriageCompareEdgeDistance(x.ToLD(), a0.ToLD(), a1.ToLD(), r.Length2.ToLD());
            if (sign != 0) return sign;
        }
        return ExactCompareEdgeDistance(x, a0, a1, r);
    }

    // Returns -1, 0, or +1 according to whether the distance from edge A edge B
    // is less than, equal to, or greater than "r" respectively.  Distances are
    // measured with respect the positions of all points as though they were
    // projected to lie exactly on the surface of the unit sphere.
    //
    // REQUIRES: A0 and A1 do not project to antipodal points (e.g., A0 == -A1).
    // REQUIRES: B0 and B1 do not project to antipodal points (e.g., B0 == -B1).
    public static int CompareEdgePairDistance(S2Point a0, S2Point a1, S2Point b0, S2Point b1, S1ChordAngle r)
    {
        // The following logic mimics S2::UpdateEdgePairMinDistance().

        // If the edges cross or share an endpoint, the minimum distance is zero.
        if (S2.CrossingSign(a0, a1, b0, b1) >= 0)
        {
            return (r.Length2 > 0) ? -1 : (r.Length2 < 0) ? 1 : 0;
        }
        // Otherwise, the minimum distance is achieved at an endpoint of at least
        // one of the two edges.
        return new[]
        {
            CompareEdgeDistance(a0, b0, b1, r),
            CompareEdgeDistance(a1, b0, b1, r),
            CompareEdgeDistance(b0, a0, a1, r),
            CompareEdgeDistance(b1, a0, a1, r)
        }.Min();
    }

    // Returns -1, 0, or +1 according to whether the normal of edge A has
    // negative, zero, or positive dot product with the normal of edge B.  This
    // essentially measures whether the edges A and B are closer to proceeding in
    // the same direction or in opposite directions around the sphere.
    //
    // This method returns an exact result, i.e. the result is zero if and only if
    // the two edges are exactly perpendicular or at least one edge is degenerate.
    // (i.e., both edge endpoints project to the same point on the sphere).
    //
    // CAVEAT: This method does not use symbolic perturbations.  Therefore it can
    // return zero even when A0 != A1 and B0 != B1, e.g. if (A0 == C * A1) exactly
    // for someant C > 0 (which is possible even when both points are
    // considered "normalized").
    //
    // REQUIRES: Neither edge can consist of antipodal points (e.g., A0 == -A1)
    //           (see comments in CompareEdgeDistance).
    public static int CompareEdgeDirections(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
    {
        // Check that no edge consists of antipodal points.  (This catches the most
        // common case -- the full test is in ExactCompareEdgeDirections.)
        MyDebug.Assert(a0 != -a1);
        MyDebug.Assert(b0 != -b1);

        int sign = TriageCompareEdgeDirections(a0, a1, b0, b1);
        if (sign != 0) return sign;

        // Optimization for the case where either edge is degenerate.
        if (a0 == a1 || b0 == b1) return 0;

        if (kHasLongDouble)
        {
            sign = TriageCompareEdgeDirections(a0.ToLD(), a1.ToLD(), b0.ToLD(), b1.ToLD());
            if (sign != 0) return sign;
        }
        return ExactCompareEdgeDirections(a0.ToExact(), a1.ToExact(), b0.ToExact(), b1.ToExact());
    }

    // Returns Sign(X0, X1, Z) where Z is the circumcenter of triangle ABC.
    // The return value is +1 if Z is to the left of edge X, and -1 if Z is to the
    // right of edge X.  The return value is zero if A == B, B == C, or C == A
    // (exactly), and also if X0 and X1 project to identical points on the sphere
    // (e.g., X0 == X1).
    //
    // The result is determined with respect to the positions of all points as
    // though they were projected to lie exactly on the surface of the unit
    // sphere.  Furthermore this method uses symbolic perturbations to compute a
    // consistent non-zero result even when Z lies exactly on edge X.
    //
    // REQUIRES: X0 and X1 do not project to antipodal points (e.g., X0 == -X1)
    //           (see comments in CompareEdgeDistance).
    public static int EdgeCircumcenterSign(S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c)
    {
        // Check that the edge does not consist of antipodal points.  (This catches
        // the most common case -- the full test is in ExactEdgeCircumcenterSign.)
        MyDebug.Assert(x0 != -x1);

        int abc_sign = Sign(a, b, c);
        int sign = TriageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign);
        if (sign != 0) return sign;

        // Optimization for the cases that are going to return zero anyway, in order
        // to avoid falling back to exact arithmetic.
        if (x0 == x1 || a == b || b == c || c == a) return 0;

        if (kHasLongDouble)
        {
            sign = TriageEdgeCircumcenterSign(
                x0.ToLD(), x1.ToLD(), a.ToLD(), b.ToLD(), c.ToLD(), abc_sign);
            if (sign != 0) return sign;
        }
        sign = ExactEdgeCircumcenterSign(x0.ToExact(), x1.ToExact(), a.ToExact(), b.ToExact(), c.ToExact(), abc_sign);
        if (sign != 0) return sign;

        // Unlike the other methods, SymbolicEdgeCircumcenterSign does not depend
        // on the sign of triangle ABC.
        return SymbolicEdgeCircumcenterSign(x0, x1, a, b, c);
    }

    // This is a specialized method that is used to compute the intersection of an
    // edge X with the Voronoi diagram of a set of points, where each Voronoi
    // region is intersected with a disc of fixed radius "r".
    //
    // Given two sites A and B and an edge (X0, X1) such that d(A,X0) < d(B,X0)
    // and both sites are within the given distance "r" of edge X, this method
    // intersects the Voronoi region of each site with a disc of radius r and
    // determines whether either region has an empty intersection with edge X.  It
    // returns FIRST if site A has an empty intersection, SECOND if site B has an
    // empty intersection, NEITHER if neither site has an empty intersection, or
    // UNCERTAIN if A == B exactly.  Note that it is not possible for both
    // intersections to be empty because of the requirement that both sites are
    // within distance r of edge X.  (For example, the only reason that Voronoi
    // region A can have an empty intersection with X is that site B is closer to
    // all points on X that are within radius r of site A.)
    //
    // The result is determined with respect to the positions of all points as
    // though they were projected to lie exactly on the surface of the unit
    // sphere.  Furthermore this method uses symbolic perturbations to compute a
    // consistent non-zero result even when A and B lie on opposite sides of X
    // such that the Voronoi edge between them exactly coincides with edge X, or
    // when A and B are distinct but project to the same point on the sphere
    // (i.e., they are linearly dependent).
    //
    // REQUIRES: r < S1ChordAngle.Right (90 degrees)
    // REQUIRES: S2Pred.CompareDistances(x0, a, b) < 0
    // REQUIRES: S2Pred.CompareEdgeDistance(a, x0, x1, r) <= 0
    // REQUIRES: S2Pred.CompareEdgeDistance(b, x0, x1, r) <= 0
    // REQUIRES: X0 and X1 do not project to antipodal points (e.g., X0 == -X1)
    //           (see comments in CompareEdgeDistance).
    public static Excluded GetVoronoiSiteExclusion(S2Point a, S2Point b, S2Point x0, S2Point x1, S1ChordAngle r)
    {
        MyDebug.Assert(r < S1ChordAngle.Right);
        MyDebug.Assert(CompareDistances(x0, a, b) < 0, "(implies a != b)");
        MyDebug.Assert(CompareEdgeDistance(a, x0, x1, r) <= 0);
        MyDebug.Assert(CompareEdgeDistance(b, x0, x1, r) <= 0);
        // Check that the edge does not consist of antipodal points.  (This catches
        // the most common case -- the full test is in ExactVoronoiSiteExclusion.)
        MyDebug.Assert(x0 != -x1);

        // If one site is closer than the other to both endpoints of X, then it is
        // closer to every point on X.  Note that this also handles the case where A
        // and B are equidistant from every point on X (i.e., X is the perpendicular
        // bisector of AB), because CompareDistances uses symbolic perturbations to
        // ensure that either A or B is considered closer (in a consistent way).
        // This also ensures that the choice of A or B does not depend on the
        // direction of X.
        if (CompareDistances(x1, a, b) < 0)
        {
            return Excluded.SECOND;  // Site A is closer to every point on X.
        }

        Excluded result = TriageVoronoiSiteExclusion(a, b, x0, x1, r.Length2);
        if (result != Excluded.UNCERTAIN) return result;

        if (kHasLongDouble)
        {
            result = TriageVoronoiSiteExclusion(a.ToLD(), b.ToLD(), x0.ToLD(), x1.ToLD(), r.Length2.ToLD());
            if (result != Excluded.UNCERTAIN) return result;
        }

        return ExactVoronoiSiteExclusion(a.ToExact(), b.ToExact(), x0.ToExact(), x1.ToExact(), (decimal)r.Length2);
    }

    public enum Excluded
    {
        None, //invalid value
        FIRST,
        SECOND,
        NEITHER,
        UNCERTAIN,
    }

    #region Low-Level Methods

    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use, since they require the client to do
    // all the actual crossing tests.

    /// <summary>
    /// A more efficient version of Sign that allows the precomputed
    /// cross-product of A and B to be specified.  (Unlike the 3 argument
    /// version this method is also inlined.)  Note that "a_cross_b" must be
    /// computed using CrossProd rather than S2::RobustCrossProd.
    ///
    /// REQUIRES: a_cross_b == a.CrossProd(b)
    /// </summary>
    public static int Sign(S2Point a, S2Point b, S2Point c, S2Point a_cross_b)
    {
        var sign = TriageSign(a, b, c, a_cross_b);
        if (sign == 0) sign = ExpensiveSign(a, b, c);
        return sign;
    }

    /// <summary>
    /// This version of Sign returns +1 if the points are definitely CCW, -1 if
    /// they are definitely CW, and 0 if two points are identical or the result
    /// is uncertain.  Uncertain cases can be resolved, if desired, by calling
    /// ExpensiveSign.
    ///
    /// The purpose of this method is to allow additional cheap tests to be done,
    /// where possible, in order to avoid calling ExpensiveSign unnecessarily.
    ///
    /// REQUIRES: a_cross_b == a.CrossProd(b)
    /// </summary>
    public static int TriageSign(S2Point a, S2Point b, S2Point c, S2Point a_cross_b)
    {
        // kMaxDetError is the maximum error in computing (AxB).C where all vectors
        // are unit length.  Using standard inequalities, it can be shown that
        //
        //  fl(AxB) = AxB + D where |D| <= (|AxB| + (2/sqrt(3))*|A|*|B|) * e
        //
        // where "fl()" denotes a calculation done in floating-point arithmetic,
        // |x| denotes either absolute value or the L2-norm as appropriate, and
        // e = 0.5*S2Constants.DoubleEpsilon.  Similarly,
        //
        //  fl(B.C) = B.C + d where |d| <= (1.5*|B.C| + 1.5*|B|*|C|) * e .
        //
        // Applying these bounds to the unit-length vectors A,B,C and neglecting
        // relative error (which does not affect the sign of the result), we get
        //
        //  fl((AxB).C) = (AxB).C + d where |d| <= (2.5 + 2/sqrt(3)) * e
        //
        // which is about 3.6548 * e, or 1.8274 * S2Constants.DoubleEpsilon.
        double kMaxDetError = 1.8274 * S2.DoubleEpsilon;
        MyDebug.Assert(a.IsUnitLength());
        MyDebug.Assert(b.IsUnitLength());
        MyDebug.Assert(c.IsUnitLength());
        MyDebug.Assert(a_cross_b == a.CrossProd(b));
        double det = a_cross_b.DotProd(c);

#if s2debug
        // Double-check borderline cases in debug mode.
        MyDebug.Assert(Math.Abs(det) <= kMaxDetError || Math.Abs(det) >= 100 * kMaxDetError || det * ExpensiveSign(a, b, c) > 0);
#endif

        if (det > kMaxDetError) return 1;
        if (det < -kMaxDetError) return -1;
        return 0;
    }

    // This function is invoked by Sign() if the sign of the determinant is
    // uncertain.  It always returns a non-zero result unless two of the input
    // points are the same.  It uses a combination of multiple-precision
    // arithmetic and symbolic perturbations to ensure that its results are
    // always self-consistent (cf. Simulation of Simplicity, Edelsbrunner and
    // Muecke).  The basic idea is to assign an infinitesimal symbolic
    // perturbation to every possible S2Point such that no three S2Points are
    // collinear and no four S2Points are coplanar.  These perturbations are so
    // small that they do not affect the sign of any determinant that was
    // non-zero before the perturbations.  If "perturb" is false, then instead
    // the exact sign of the unperturbed input points is returned, which can be
    // zero even when all three points are distinct.
    //
    // Unlike Sign(), this method does not require the input points to be
    // normalized.
    //
    // ExpensiveSign() uses arbitrary-precision arithmetic and the "simulation of
    // simplicity" technique in order to be completely robust (i.e., to return
    // consistent results for all possible inputs).
    public static int ExpensiveSign(S2Point a, S2Point b, S2Point c, bool perturb = true)
    {
        // Return zero if and only if two points are the same.  This ensures (1).
        if (a == b || b == c || c == a) return 0;

        // Next we try recomputing the determinant still using floating-point
        // arithmetic but in a more precise way.  This is more expensive than the
        // simple calculation done by TriageSign(), but it is still *much* cheaper
        // than using arbitrary-precision arithmetic.  This optimization is able to
        // compute the correct determinant sign in virtually all cases except when
        // the three points are truly collinear (e.g., three points on the equator).
        var det_sign = StableSign(a, b, c);
        if (det_sign != 0) return det_sign;
        return ExactSign(a, b, c, perturb);
    }

    #endregion

    // All error bounds in this file are expressed in terms of the maximum
    // rounding error for a floating-point type.  The rounding error is half of
    // the numeric_limits<T>.epsilon() value.
    // constexpr double DBL_ERR = rounding_epsilon<double>();
    // constexpr double LD_ERR = rounding_epsilon<double>();
    public const double TT_ERR = S2.DoubleEpsilon;
    public const double DBL_ERR = S2.DoubleEpsilon;
    public const double LD_ERR = S2.DoubleEpsilon;
    public const bool kHasLongDouble = LD_ERR < DBL_ERR;

    // Define sqrt(3) as a constant so that we can use it with constexpr.
    // Unfortunately we can't use M_SQRT3 because some client libraries define
    // this symbol without first checking whether it already exists.
    public const double kSqrt3 = 1.7320508075688772935274463415058;

    // A predefined S1ChordAngle representing (approximately) 45 degrees.
    private static S1ChordAngle K45Degrees { get; } = S1ChordAngle.FromLength2(2 - S2.M_SQRT2);

    // Efficiently tests whether an ExactFloat vector is (0, 0, 0).
    public static bool IsZero(S2Point a) =>
        double.Sign(a[0]) == 0 && double.Sign(a[1]) == 0 && double.Sign(a[2]) == 0;
    public static bool IsZero(Vector3<ExactFloat> a) =>
        ExactFloat.Sign(a[0]) == 0 && ExactFloat.Sign(a[1]) == 0 && ExactFloat.Sign(a[2]) == 0;

    // Compute the determinant in a numerically stable way.  Unlike TriageSign(),
    // this method can usually compute the correct determinant sign even when all
    // three points are as collinear as possible.  For example if three points are
    // spaced 1km apart along a random line on the Earth's surface using the
    // nearest representable points, there is only a 0.4% chance that this method
    // will not be able to find the determinant sign.  The probability of failure
    // decreases as the points get closer together; if the collinear points are
    // 1 meter apart, the failure rate drops to 0.0004%.
    //
    // This method could be extended to also handle nearly-antipodal points (and
    // in fact an earlier version of this code did exactly that), but antipodal
    // points are rare in practice so it seems better to simply fall back to
    // exact arithmetic in that case.
    public static int StableSign(S2Point a, S2Point b, S2Point c)
    {
        var ab = b - a;
        var bc = c - b;
        var ca = a - c;
        var ab2 = ab.Norm2();
        var bc2 = bc.Norm2();
        var ca2 = ca.Norm2();

        // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
        // cyclically permuted if necessary so that AB is the longest edge.  (This
        // minimizes the magnitude of cross product.)  At the same time we also
        // compute the maximum error in the determinant.  Using a similar technique
        // to the one used for kMaxDetError, the error is at most
        //
        //   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
        //
        // where e = 0.5 * S2Constants.DoubleEpsilon.  If the determinant magnitude is larger than
        // this value then we know its sign with certainty.
        var kDetErrorMultiplier = 3.2321 * S2.DoubleEpsilon;  // see above
        double det, max_error;
        if (ab2 >= bc2 && ab2 >= ca2)
        {
            // AB is the longest edge, so compute (A-C)x(B-C).C.
            det = -ca.CrossProd(bc).DotProd(c);
            max_error = kDetErrorMultiplier * Math.Sqrt(ca2 * bc2);
        }
        else if (bc2 >= ca2)
        {
            // BC is the longest edge, so compute (B-A)x(C-A).A.
            det = -ab.CrossProd(ca).DotProd(a);
            max_error = kDetErrorMultiplier * Math.Sqrt(ab2 * ca2);
        }
        else
        {
            // CA is the longest edge, so compute (C-B)x(A-B).B.
            det = -bc.CrossProd(ab).DotProd(b);
            max_error = kDetErrorMultiplier * Math.Sqrt(bc2 * ab2);
        }
        // Errors smaller than this value may not be accurate due to underflow.
        var kMinNoUnderflowError = kDetErrorMultiplier * Math.Sqrt(double.MinValue);
        if (max_error < kMinNoUnderflowError) return 0;

        return (Math.Abs(det) <= max_error) ? 0 : (det > 0) ? 1 : -1;
    }

    // The following function returns the sign of the determinant of three points
    // A, B, C under a model where every possible S2Point is slightly perturbed by
    // a unique infinitesmal amount such that no three perturbed points are
    // collinear and no four points are coplanar.  The perturbations are so small
    // that they do not change the sign of any determinant that was non-zero
    // before the perturbations, and therefore can be safely ignored unless the
    // determinant of three points is exactly zero (using multiple-precision
    // arithmetic).
    //
    // Since the symbolic perturbation of a given point is fixed (i.e., the
    // perturbation is the same for all calls to this method and does not depend
    // on the other two arguments), the results of this method are always
    // self-consistent.  It will never return results that would correspond to an
    // "impossible" configuration of non-degenerate points.
    //
    // Requirements:
    //   The 3x3 determinant of A, B, C must be exactly zero.
    //   The points must be distinct, with A < B < C in lexicographic order.
    //
    // Returns:
    //   +1 or -1 according to the sign of the determinant after the symbolic
    // perturbations are taken into account.
    //
    // Reference:
    //   "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM Transactions on
    //   Graphics, 1990).
    //
    private static int SymbolicallyPerturbedSign(S2Point a, S2Point b, S2Point c, S2Point b_cross_c)
    {
        // This method requires that the points are sorted in lexicographically
        // increasing order.  This is because every possible S2Point has its own
        // symbolic perturbation such that if A < B then the symbolic perturbation
        // for A is much larger than the perturbation for B.
        //
        // Alternatively, we could sort the points in this method and keep track of
        // the sign of the permutation, but it is more efficient to do this before
        // converting the inputs to the multi-precision representation, and this
        // also lets us re-use the result of the cross product B x C.
        MyDebug.Assert(a < b && b < c);

        // Every input coordinate x[i] is assigned a symbolic perturbation dx[i].
        // We then compute the sign of the determinant of the perturbed points,
        // i.e.
        //               | a[0]+da[0]  a[1]+da[1]  a[2]+da[2] |
        //               | b[0]+db[0]  b[1]+db[1]  b[2]+db[2] |
        //               | c[0]+dc[0]  c[1]+dc[1]  c[2]+dc[2] |
        //
        // The perturbations are chosen such that
        //
        //   da[2] > da[1] > da[0] > db[2] > db[1] > db[0] > dc[2] > dc[1] > dc[0]
        //
        // where each perturbation is so much smaller than the previous one that we
        // don't even need to consider it unless the coefficients of all previous
        // perturbations are zero.  In fact, it is so small that we don't need to
        // consider it unless the coefficient of all products of the previous
        // perturbations are zero.  For example, we don't need to consider the
        // coefficient of db[1] unless the coefficient of db[2]*da[0] is zero.
        //
        // The follow code simply enumerates the coefficients of the perturbations
        // (and products of perturbations) that appear in the determinant above, in
        // order of decreasing perturbation magnitude.  The first non-zero
        // coefficient determines the sign of the result.  The easiest way to
        // enumerate the coefficients in the correct order is to pretend that each
        // perturbation is some tiny value "eps" raised to a power of two:
        //
        // eps**    1      2      4      8     16     32     64     128    256
        //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]  dc[2]  dc[1]  dc[0]
        //
        // Essentially we can then just count in binary and test the corresponding
        // subset of perturbations at each step.  So for example, we must test the
        // coefficient of db[2]*da[0] before db[1] because eps**12 > eps**16.
        //
        // Of course, not all products of these perturbations appear in the
        // determinant above, since the determinant only contains the products of
        // elements in distinct rows and columns.  Thus we don't need to consider
        // da[2]*da[1], db[1]*da[1], etc.  Furthermore, sometimes different pairs of
        // perturbations have the same coefficient in the determinant; for example,
        // da[1]*db[0] and db[1]*da[0] have the same coefficient (c[2]).  Therefore
        // we only need to test this coefficient the first time we encounter it in
        // the binary order above (which will be db[1]*da[0]).
        //
        // The sequence of tests below also appears in Table 4-ii of the paper
        // referenced above, if you just want to look it up, with the following
        // translations: [a,b,c] . [i,j,k] and [0,1,2] . [1,2,3].  Also note that
        // some of the signs are different because the opposite cross product is
        // used (e.g., B x C rather than C x B).

        int det_sign = double.Sign(b_cross_c[2]);            // da[2]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(b_cross_c[1]);                // da[1]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(b_cross_c[0]);                // da[0]
        if (det_sign != 0) return det_sign;

        det_sign = double.Sign(c[0] * a[1] - c[1] * a[0]);     // db[2]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(c[0]);                        // db[2] * da[1]
        if (det_sign != 0) return det_sign;
        det_sign = -double.Sign(c[1]);                     // db[2] * da[0]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(c[2] * a[0] - c[0] * a[2]);     // db[1]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(c[2]);                        // db[1] * da[0]
        if (det_sign != 0) return det_sign;
        // The following test is listed in the paper, but it is redundant because
        // the previous tests guarantee that C == (0, 0, 0).
        MyDebug.Assert(0 == double.Sign(c[1] * a[2] - c[2] * a[1]), "db[0]");

        det_sign = double.Sign(a[0] * b[1] - a[1] * b[0]);     // dc[2]
        if (det_sign != 0) return det_sign;
        det_sign = -double.Sign(b[0]);                     // dc[2] * da[1]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(b[1]);                        // dc[2] * da[0]
        if (det_sign != 0) return det_sign;
        det_sign = double.Sign(a[0]);                        // dc[2] * db[1]
        if (det_sign != 0) return det_sign;
        return 1;                                     // dc[2] * db[1] * da[0]
    }

    /// <summary>
    /// Compute the determinant using exact arithmetic and/or symbolic
    /// permutations.  Requires that the three points are distinct.
    /// </summary>
    public static int ExactSign(S2Point a, S2Point b, S2Point c, bool perturb)
    {
        MyDebug.Assert(a != b && b != c && c != a);

        // Sort the three points in lexicographic order, keeping track of the sign
        // of the permutation.  (Each exchange inverts the sign of the determinant.)
        int perm_sign = 1;
        S2Point buff;
        if (a > b) { buff = a; a = b; b = buff; perm_sign = -perm_sign; }
        if (b > c) { buff = b; b = c; c = buff; perm_sign = -perm_sign; }
        if (a > b) { buff = a; a = b; b = buff; perm_sign = -perm_sign; }
        MyDebug.Assert(a < b && b < c);

        // ~Construct multiple-precision versions of the sorted points and~ compute
        // their exact 3x3 determinant.
        var b_cross_c = b.CrossProd(c);
        var det = a.DotProd(b_cross_c);

        // The precision of ExactFloat is high enough that the result should always
        // be exact (no rounding was performed).
        MyDebug.Assert(!double.IsNaN(det));
        //Assert.True(det.prec() < det.max_prec());

        // If the exact determinant is non-zero, we're done.
        int det_sign = double.Sign(det);
        if (det_sign == 0 && perturb)
        {
            // Otherwise, we need to resort to symbolic perturbations to resolve the
            // sign of the determinant.
            det_sign = SymbolicallyPerturbedSign(a, b, c, b_cross_c);
            MyDebug.Assert(0 != det_sign);
        }
        return perm_sign * det_sign;
    }

    // Returns cos(XY), and sets "error" to the maximum error in the result.
    // REQUIRES: "x" and "y" satisfy S2.IsNormalized().
    private static double GetCosDistance(S2Point x, S2Point y, out double error)
    {
        double c = x.DotProd(y);
        error = 9.5 * DBL_ERR * Math.Abs(c) + 1.5 * DBL_ERR;
        return c;
    }

#pragma warning disable IDE0051 // Remove unused private members
    private static double GetCosDistance_HighPrecision(S2Point x, S2Point y, out double error)
#pragma warning restore IDE0051 // Remove unused private members
    {
        // With "double" precision it is worthwhile to compensate for length
        // errors in "x" and "y", since they are only unit length to within the
        // precision of "double".  (This would also reduce the errorant
        // slightly in the method above but is not worth the additional effort.)
        double c = x.DotProd(y) / Math.Sqrt(x.Norm2() * y.Norm2());
        error = 7 * LD_ERR * Math.Abs(c) + 1.5 * LD_ERR;
        return c;
    }

    // Returns sin**2(XY), where XY is the angle between X and Y, and sets "error"
    // to the maximum error in the result.
    //
    // REQUIRES: "x" and "y" satisfy S2.IsNormalized().
    private static double GetSin2Distance(S2Point x, S2Point y, out double error)
    {
        // The (x-y).CrossProd(x+y) trick eliminates almost all of error due to "x"
        // and "y" being not quite unit length.  This method is extremely accurate
        // for small distances; the *relative* error in the result is O(DBL_ERR) for
        // distances as small as DBL_ERR.
        S2Point n = (x - y).CrossProd(x + y);
        double d2 = 0.25 * n.Norm2();
        error = (21 + 4 * Math.Sqrt(3)) * DBL_ERR * d2 +
                  32 * Math.Sqrt(3) * DBL_ERR * DBL_ERR * Math.Sqrt(d2) +
                  768 * DBL_ERR * DBL_ERR * DBL_ERR * DBL_ERR;
        return d2;
    }

#pragma warning disable IDE0051 // Remove unused private members
    private static double GetSin2Distance_HighPrecision(S2Point x, S2Point y, out double error)
#pragma warning restore IDE0051 // Remove unused private members
    {
        // In "double" precision it is worthwhile to compensate for length
        // errors in "x" and "y", since they are only unit length to within the
        // precision of "double".  Otherwise the "d2" error coefficient below would
        // be (16 * DBL_ERR + (5 + 4 * sqrt(3)) * LD_ERR), which is much larger.
        // (Dividing by the squared norms of "x" and "y" would also reduce the error
        // constant slightly in the double-precision version, but this is not worth
        // the additional effort.)
        S2Point n = (x - y).CrossProd(x + y);
        double d2 = 0.25 * n.Norm2() / (x.Norm2() * y.Norm2());
        error = (13 + 4 * Math.Sqrt(3)) * LD_ERR * d2 +
                  32 * Math.Sqrt(3) * DBL_ERR * LD_ERR * Math.Sqrt(d2) +
                  768 * DBL_ERR * DBL_ERR * LD_ERR * LD_ERR;
        return d2;
    }

    private static int TriageCompareCosDistances(S2Point x, S2Point a, S2Point b)
    {
        double cos_ax = GetCosDistance(a, x, out var cos_ax_error);
        double cos_bx = GetCosDistance(b, x, out var cos_bx_error);
        double diff = cos_ax - cos_bx;
        double error = cos_ax_error + cos_bx_error;
        return (diff > error) ? -1 : (diff < -error) ? 1 : 0;
    }

    private static int TriageCompareSin2Distances(S2Point x, S2Point a, S2Point b)
    {
        double sin2_ax = GetSin2Distance(a, x, out var sin2_ax_error);
        double sin2_bx = GetSin2Distance(b, x, out var sin2_bx_error);
        double diff = sin2_ax - sin2_bx;
        double error = sin2_ax_error + sin2_bx_error;
        return (diff > error) ? 1 : (diff < -error) ? -1 : 0;
    }

    private static int ExactCompareDistances(Vector3<ExactFloat> x, Vector3<ExactFloat> a, Vector3<ExactFloat> b)
    {
        // This code produces the same result as though all points were reprojected
        // to lie exactly on the surface of the unit sphere.  It is based on testing
        // whether x.DotProd(a.Normalize()) < x.DotProd(b.Normalize()), reformulated
        // so that it can be evaluated using exact arithmetic.
        var cos_ax = x.DotProd(a);
        var cos_bx = x.DotProd(b);
        // If the two values have different signs, we need to handle that case now
        // before squaring them below.
        int a_sign = ExactFloat.Sign(cos_ax), b_sign = ExactFloat.Sign(cos_bx);
        if (a_sign != b_sign)
        {
            return (a_sign > b_sign) ? -1 : 1;  // If cos(AX) > cos(BX), then AX < BX.
        }
        var cmp = cos_bx * cos_bx * a.Norm2() - cos_ax * cos_ax * b.Norm2();
        return a_sign * ExactFloat.Sign(cmp);
    }

    // Given three points such that AX == BX (exactly), returns -1, 0, or +1
    // according whether AX < BX, AX == BX, or AX > BX after symbolic
    // perturbations are taken into account.
#pragma warning disable IDE0060 // Remove unused parameter
    private static int SymbolicCompareDistances(S2Point x, S2Point a, S2Point b)
#pragma warning restore IDE0060 // Remove unused parameter
    {
        // Our symbolic perturbation strategy is based on the following model.
        // Similar to "simulation of simplicity", we assign a perturbation to every
        // point such that if A < B, then the symbolic perturbation for A is much,
        // much larger than the symbolic perturbation for B.  We imagine that
        // rather than projecting every point to lie exactly on the unit sphere,
        // instead each point is positioned on its own tiny pedestal that raises it
        // just off the surface of the unit sphere.  This means that the distance AX
        // is actually the true distance AX plus the (symbolic) heights of the
        // pedestals for A and X.  The pedestals are infinitesmally thin, so they do
        // not affect distance measurements except at the two endpoints.  If several
        // points project to exactly the same point on the unit sphere, we imagine
        // that they are placed on separate pedestals placed close together, where
        // the distance between pedestals is much, much less than the height of any
        // pedestal.  (There are a finite number of S2Points, and therefore a finite
        // number of pedestals, so this is possible.)
        //
        // If A < B, then A is on a higher pedestal than B, and therefore AX > BX.
        return (a < b) ? 1 : (a > b) ? -1 : 0;
    }

    private static int CompareSin2Distances(S2Point x, S2Point a, S2Point b)
    {
        int sign = TriageCompareSin2Distances(x, a, b);
        if (kHasLongDouble && sign == 0)
        {
            sign = TriageCompareSin2Distances(x.ToLD(), a.ToLD(), b.ToLD());
        }
        return sign;
    }

    private static int TriageCompareCosDistance(S2Point x, S2Point y, double r2)
    {
        double cos_xy = GetCosDistance(x, y, out var cos_xy_error);
        double cos_r = 1 - 0.5 * r2;
        double cos_r_error = 2 * TT_ERR * cos_r;
        double diff = cos_xy - cos_r;
        double error = cos_xy_error + cos_r_error;
        return (diff > error) ? -1 : (diff < -error) ? 1 : 0;
    }

    private static int TriageCompareSin2Distance(S2Point x, S2Point y, double r2)
    {
        MyDebug.Assert(r2 < 2.0, "Only valid for distance limits < 90 degrees.");

        double sin2_xy = GetSin2Distance(x, y, out var sin2_xy_error);
        double sin2_r = r2 * (1 - 0.25 * r2);
        double sin2_r_error = 3 * TT_ERR * sin2_r;
        double diff = sin2_xy - sin2_r;
        double error = sin2_xy_error + sin2_r_error;
        return (diff > error) ? 1 : (diff < -error) ? -1 : 0;
    }

    private static int ExactCompareDistance(Vector3<ExactFloat> x, Vector3<ExactFloat> y, ExactFloat r2) =>
        ExactCompareDistance(x, y, r2.Value);
    private static int ExactCompareDistance(Vector3<ExactFloat> x, Vector3<ExactFloat> y, decimal r2)
    {
        // This code produces the same result as though all points were reprojected
        // to lie exactly on the surface of the unit sphere.  It is based on
        // comparing the cosine of the angle XY (when both points are projected to
        // lie exactly on the sphere) to the given threshold.
        var cos_xy = x.DotProd(y).Value;
        var cos_r = 1 - 0.5M * r2;
        // If the two values have different signs, we need to handle that case now
        // before squaring them below.
        int xy_sign = MathM.Sign(cos_xy), r_sign = MathM.Sign(cos_r);
        if (xy_sign != r_sign)
        {
            return (xy_sign > r_sign) ? -1 : 1;  // If cos(XY) > cos(r), then XY < r.
        }
        var cmp = cos_r * cos_r * x.Norm2().Value * y.Norm2().Value - cos_xy * cos_xy;
        return xy_sign * MathM.Sign(cmp);
    }

    // Helper function that compares the distance XY against the squared chord
    // distance "r2" using the given precision "T".

    private static int TriageCompareDistance(S2Point x, S2Point y, double r2)
    {
        // The Sin2 method is much more accurate for small distances, but it is only
        // valid when the actual distance and the distance limit are both less than
        // 90 degrees.  So we always start with the Cos method.
        int sign = TriageCompareCosDistance(x, y, r2);
        if (sign == 0 && r2 < K45Degrees.Length2)
        {
            sign = TriageCompareSin2Distance(x, y, r2);
        }
        return sign;
    }

    // Helper function that returns "a0" or "a1", whichever is closer to "x".
    // Also returns the squared distance from the returned point to "x" in "ax2".

    private static S2Point GetClosestVertex(S2Point x, S2Point a0, S2Point a1, out double ax2)
    {
        double a0x2 = (a0 - x).Norm2();
        double a1x2 = (a1 - x).Norm2();
        if (a0x2 < a1x2 || (a0x2 == a1x2 && a0 < a1))
        {
            ax2 = a0x2;
            return a0;
        }
        else
        {
            ax2 = a1x2;
            return a1;
        }
    }

    // Helper function that returns -1, 0, or +1 according to whether the distance
    // from "x" to the great circle through (a0, a1) is less than, equal to, or
    // greater than the given squared chord length "r2".  This method computes the
    // squared sines of the distances involved, which is more accurate when the
    // distances are small (less than 45 degrees).
    //
    // The remaining parameters are functions of (a0, a1) and are passed in
    // because they have already been computed: n = (a0 - a1) x (a0 + a1),
    // n1 = n.Norm, and n2 = n.Norm2.

    private static int TriageCompareLineSin2Distance(S2Point x, S2Point a0, S2Point a1, double r2, S2Point n, double n1, double n2)
    {
        // The minimum distance is to a point on the edge interior.  Since the true
        // distance to the edge is always less than 90 degrees, we can return
        // immediately if the limit is 90 degrees or larger.
        if (r2 >= 2.0) return -1;  // distance < limit

        // Otherwise we compute sin^2(distance to edge) to get the best accuracy
        // when the distance limit is small (e.g., S2EdgeCrossings.kIntersectionError).
        double n2sin2_r = n2 * r2 * (1 - 0.25 * r2);
        double n2sin2_r_error = 6 * TT_ERR * n2sin2_r;
        double ax2, xDn = (x - GetClosestVertex(x, a0, a1, out ax2)).DotProd(n);
        double xDn2 = xDn * xDn;
        double c1 = ((3.5 + 2 * Math.Sqrt(3)) * n1 + 32 * Math.Sqrt(3) * DBL_ERR) * TT_ERR * Math.Sqrt(ax2);
        double xDn2_error = 4 * TT_ERR * xDn2 + (2 * Math.Abs(xDn) + c1) * c1;

        /*
        // If we are using extended precision, then it is worthwhile to recompute
        // the length of X more accurately.  Otherwise we use the fact that X is
        // guaranteed to be unit length to with a tolerance of 4 * DBL_ERR.
        if (T_ERR < DBL_ERR)
        {
            n2sin2_r *= x.Norm2;
            n2sin2_r_error += 4 * T_ERR * n2sin2_r;
        }
        else*/
        {
            n2sin2_r_error += 8 * DBL_ERR * n2sin2_r;
        }
        double diff = xDn2 - n2sin2_r;
        double error = xDn2_error + n2sin2_r_error;
        return (diff > error) ? 1 : (diff < -error) ? -1 : 0;
    }

    // Like TriageCompareLineSin2Distance, but this method computes the squared
    // cosines of the distances involved.  It is more accurate when the distances
    // are large (greater than 45 degrees).

#pragma warning disable IDE0060 // Quitar el parámetro no utilizado
    private static int TriageCompareLineCos2Distance(S2Point x, S2Point a0, S2Point a1, double r2, S2Point n, double n1, double n2)
#pragma warning restore IDE0060 // Quitar el parámetro no utilizado
    {
        // The minimum distance is to a point on the edge interior.  Since the true
        // distance to the edge is always less than 90 degrees, we can return
        // immediately if the limit is 90 degrees or larger.
        if (r2 >= 2.0) return -1;  // distance < limit

        // Otherwise we compute cos^2(distance to edge).
        double cos_r = 1 - 0.5 * r2;
        double n2cos2_r = n2 * cos_r * cos_r;
        double n2cos2_r_error = 7 * TT_ERR * n2cos2_r;

        // The length of M = X.CrossProd(N) is the cosine of the distance.
        double m2 = x.CrossProd(n).Norm2();
        double m1 = Math.Sqrt(m2);
        double m1_error = ((1 + 8 / Math.Sqrt(3)) * n1 + 32 * Math.Sqrt(3) * DBL_ERR) * TT_ERR;
        double m2_error = 3 * TT_ERR * m2 + (2 * m1 + m1_error) * m1_error;

        // If we are using extended precision, then it is worthwhile to recompute
        // the length of X more accurately.  Otherwise we use the fact that X is
        // guaranteed to be unit length to within a tolerance of 4 * DBL_ERR.
        /*if (T_ERR < DBL_ERR)
        {
            n2cos2_r *= x.Norm2;
            n2cos2_r_error += 4 * T_ERR * n2cos2_r;
        }
        else*/
        {
            n2cos2_r_error += 8 * DBL_ERR * n2cos2_r;
        }
        double diff = m2 - n2cos2_r;
        double error = m2_error + n2cos2_r_error;
        return (diff > error) ? -1 : (diff < -error) ? 1 : 0;
    }

    private static int TriageCompareLineDistance(S2Point x, S2Point a0, S2Point a1, double r2, S2Point n, double n1, double n2)
    {
        if (r2 < K45Degrees.Length2)
        {
            return TriageCompareLineSin2Distance(x, a0, a1, r2, n, n1, n2);
        }
        else
        {
            return TriageCompareLineCos2Distance(x, a0, a1, r2, n, n1, n2);
        }
    }

    private static int TriageCompareEdgeDistance(S2Point x, S2Point a0, S2Point a1, double r2)
    {
        // First we need to decide whether the closest point is an edge endpoint or
        // somewhere in the interior.  To determine this we compute a plane
        // perpendicular to (a0, a1) that passes through X.  Letting M be the normal
        // to this plane, the closest point is in the edge interior if and only if
        // a0.M < 0 and a1.M > 0.  Note that we can use "<" rather than "<=" because
        // if a0.M or a1.M is zero exactly then it doesn't matter which code path we
        // follow (since the distance to an endpoint and the distance to the edge
        // interior are exactly the same in this case).
        S2Point n = (a0 - a1).CrossProd(a0 + a1);
        S2Point m = n.CrossProd(x);
        // For better accuracy when the edge (a0,a1) is very short, we subtract "x"
        // before computing the dot products with M.
        S2Point a0_dir = a0 - x;
        S2Point a1_dir = a1 - x;
        double a0_sign = a0_dir.DotProd(m);
        double a1_sign = a1_dir.DotProd(m);
        double n2 = n.Norm2();
        double n1 = Math.Sqrt(n2);
        double n1_error = ((3.5 + 8 / Math.Sqrt(3)) * n1 + 32 * Math.Sqrt(3) * DBL_ERR) * TT_ERR;
        double a0_sign_error = n1_error * a0_dir.Norm();
        double a1_sign_error = n1_error * a1_dir.Norm();
        if (a0_sign < a0_sign_error && a1_sign > -a1_sign_error)
        {
            if (a0_sign > -a0_sign_error || a1_sign < a1_sign_error)
            {
                // It is uncertain whether minimum distance is to an edge vertex or to
                // the edge interior.  We compute both distances and check whether they
                // yield the same result; otherwise the result is uncertain.
                int vertex_sign = Math.Min(TriageCompareDistance(x, a0, r2),
                    TriageCompareDistance(x, a1, r2));
                int line_sign = TriageCompareLineDistance(x, a0, a1, r2, n, n1, n2);
                return (vertex_sign == line_sign) ? line_sign : 0;
            }
            // The minimum distance is to the edge interior.
            return TriageCompareLineDistance(x, a0, a1, r2, n, n1, n2);
        }
        // The minimum distance is to an edge endpoint.
        return Math.Min(TriageCompareDistance(x, a0, r2),
            TriageCompareDistance(x, a1, r2));
    }

    // REQUIRES: the closest point to "x" is in the interior of edge (a0, a1).
    private static int ExactCompareLineDistance(S2Point x, S2Point a0, S2Point a1, double r2)
    {
        // Since we are given that the closest point is in the edge interior, the
        // true distance is always less than 90 degrees (which corresponds to a
        // squared chord length of 2.0).
        if (r2 >= 2.0) return -1;  // distance < limit

        // Otherwise compute the edge normal
        S2Point n = a0.CrossProd(a1);
        double sin_d = x.DotProd(n);
        double sin2_r = r2 * (1 - 0.25 * r2);
        double cmp = sin_d * sin_d - sin2_r * x.Norm2() * n.Norm2();
        return double.Sign(cmp);
    }

    private static int ExactCompareEdgeDistance(S2Point x, S2Point a0, S2Point a1, S1ChordAngle r)
    {
        // Even if previous calculations were uncertain, we might not need to do
        // *all* the calculations in exact arithmetic here.  For example it may be
        // easy to determine whether "x" is closer to an endpoint or the edge
        // interior.  The only calculation where we always use exact arithmetic is
        // when measuring the distance to the extended line (great circle) through
        // "a0" and "a1", since it is virtually certain that the previous floating
        // point calculations failed in that case.
        //
        // CompareEdgeDirections requires that no edge has antipodal endpoints,
        // therefore we need to handle the cases a0 == -x, a1 == -x separately.
        if (a0 != -x && a1 != -x &&
                CompareEdgeDirections(a0, a1, a0, x) > 0 &&
                CompareEdgeDirections(a0, a1, x, a1) > 0)
        {
            // The closest point to "x" is along the interior of the edge.
            return ExactCompareLineDistance(x, a0, a1, r.Length2);
        }
        else
        {
            // The closest point to "x" is one of the edge endpoints.
            return Math.Min(CompareDistance(x, a0, r), CompareDistance(x, a1, r));
        }
    }

    private static int TriageCompareEdgeDirections(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
    {
        var na = (a0 - a1).CrossProd(a0 + a1);
        var nb = (b0 - b1).CrossProd(b0 + b1);
        var na_len = na.Norm();
        var nb_len = nb.Norm();
        var cos_ab = na.DotProd(nb);
        var cos_ab_error = ((5 + 4 * Math.Sqrt(3)) * na_len * nb_len +
                          32 * Math.Sqrt(3) * DBL_ERR * (na_len + nb_len)) * TT_ERR;
        return (cos_ab > cos_ab_error) ? 1 : (cos_ab < -cos_ab_error) ? -1 : 0;
    }

    private static bool ArePointsLinearlyDependent(Vector3<ExactFloat> x, Vector3<ExactFloat> y)
    {
        return IsZero(x.CrossProd(y));
    }

    private static bool ArePointsAntipodal(Vector3<ExactFloat> x, Vector3<ExactFloat> y)
    {
        return ArePointsLinearlyDependent(x, y) && ExactFloat.Sign(x.DotProd(y)) < 0;
    }

    private static int ExactCompareEdgeDirections(Vector3<ExactFloat> a0, Vector3<ExactFloat> a1, Vector3<ExactFloat> b0, Vector3<ExactFloat> b1)
    {
        MyDebug.Assert(!ArePointsAntipodal(a0, a1));
        MyDebug.Assert(!ArePointsAntipodal(b0, b1));
        return ExactFloat.Sign(a0.CrossProd(a1).DotProd(b0.CrossProd(b1)));
    }

    // If triangle ABC has positive sign, returns its circumcenter.  If ABC has
    // negative sign, returns the negated circumcenter.

    private static S2Point GetCircumcenter(S2Point a, S2Point b, S2Point c, out double error)
    {
        // We compute the circumcenter using the intersection of the perpendicular
        // bisectors of AB and BC.  The formula is essentially
        //
        //    Z = ((A x B) x (A + B)) x ((B x C) x (B + C)),
        //
        // except that we compute the cross product (A x B) as (A - B) x (A + B)
        // (and similarly for B x C) since this is much more stable when the inputs
        // are unit vectors.
        S2Point ab_diff = a - b, ab_sum = a + b;
        S2Point bc_diff = b - c, bc_sum = b + c;
        S2Point nab = ab_diff.CrossProd(ab_sum);
        var nab_len = nab.Norm();
        var ab_len = ab_diff.Norm();
        S2Point nbc = bc_diff.CrossProd(bc_sum);
        var nbc_len = nbc.Norm();
        var bc_len = bc_diff.Norm();
        S2Point mab = nab.CrossProd(ab_sum);
        S2Point mbc = nbc.CrossProd(bc_sum);
        error = ((16 + 24 * Math.Sqrt(3)) * TT_ERR +
                      8 * DBL_ERR * (ab_len + bc_len)) * nab_len * nbc_len +
                     128 * Math.Sqrt(3) * DBL_ERR * TT_ERR * (nab_len + nbc_len) +
                     3 * 4096 * DBL_ERR * DBL_ERR * TT_ERR * TT_ERR;
        return mab.CrossProd(mbc);
    }

    private static int TriageEdgeCircumcenterSign(S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c, int abc_sign)
    {
        // Compute the circumcenter Z of triangle ABC, and then test which side of
        // edge X it lies on.
        S2Point z = GetCircumcenter(a, b, c, out var z_error);
        S2Point nx = (x0 - x1).CrossProd(x0 + x1);
        // If the sign of triangle ABC is negative, then we have computed -Z and the
        // result should be negated.
        double result = abc_sign * nx.DotProd(z);

        double z_len = z.Norm();
        double nx_len = nx.Norm();
        double nx_error = ((1 + 2 * Math.Sqrt(3)) * nx_len + 32 * Math.Sqrt(3) * DBL_ERR) * TT_ERR;
        double result_error = (3 * TT_ERR * nx_len + nx_error) * z_len + z_error * nx_len;
        return (result > result_error) ? 1 : (result < -result_error) ? -1 : 0;
    }

    private static int ExactEdgeCircumcenterSign(Vector3<ExactFloat> x0, Vector3<ExactFloat> x1, Vector3<ExactFloat> a, Vector3<ExactFloat> b, Vector3<ExactFloat> c, int abc_sign)
    {
        // Return zero if the edge X is degenerate.  (Also see the comments in
        // SymbolicEdgeCircumcenterSign.)
        if (ArePointsLinearlyDependent(x0, x1))
        {
            MyDebug.Assert(x0.DotProd(x1).Value > 0M, "Antipodal edges not allowed.");
            return 0;
        }
        // The simplest predicate for testing whether the sign is positive is
        //
        // (1)  (X0 x X1) . (|C|(A x B) + |A|(B x C) + |B|(C x A)) > 0
        //
        // where |A| denotes A.Norm and the expression after the "." represents
        // the circumcenter of triangle ABC.  (This predicate is terrible from a
        // numerical accuracy point of view, but that doesn't matter since we are
        // going to use exact arithmetic.)  This predicate also assumes that
        // triangle ABC is CCW (positive sign); we correct for that below.
        //
        // The only problem with evaluating this inequality is that computing |A|,
        // |B| and |C| requires square roots.  To avoid this problem we use the
        // standard technique of rearranging the inequality to isolate at least one
        // square root and then squaring both sides.  We need to repeat this process
        // twice in order to eliminate all the square roots, which leads to a
        // polynomial predicate of degree 20 in the input arguments.
        //
        // Rearranging (1) we get
        //
        //      (X0 x X1) . (|C|(A x B) + |A|(B x C)) > |B|(X0 x X1) . (A x C)
        //
        // Before squaring we need to check the sign of each side.  If the signs are
        // different then we know the result without squaring, and if the signs are
        // both negative then after squaring both sides we need to invert the
        // result.  Define
        //
        //      dAB = (X0 x X1) . (A x B)
        //      dBC = (X0 x X1) . (B x C)
        //      dCA = (X0 x X1) . (C x A)
        //
        // Then we can now write the inequality above as
        //
        // (2)  |C| dAB + |A| dBC > -|B| dCA
        //
        // The RHS of (2) is positive if dCA < 0, and the LHS of (2) is positive if
        // (|C| dAB + |A| dBC) > 0.  Since the LHS has square roots, we need to
        // eliminate them using the same process.  Rewriting the LHS as
        //
        // (3)  |C| dAB > -|A| dBC
        //
        // we again need to check the signs of both sides.  Let's start with that.
        // We also precompute the following values because they are used repeatedly
        // when squaring various expressions below:
        //
        //     abc2 = |A|^2 dBC^2
        //     bca2 = |B|^2 dCA^2
        //     cab2 = |C|^2 dAB^2
        var nx = x0.CrossProd(x1);
        var dab = nx.DotProd(a.CrossProd(b));
        var dbc = nx.DotProd(b.CrossProd(c));
        var dca = nx.DotProd(c.CrossProd(a));
        var abc2 = a.Norm2() * (dbc * dbc);
        var bca2 = b.Norm2() * (dca * dca);
        var cab2 = c.Norm2() * (dab * dab);

        // If the two sides of (3) have different signs (including the case where
        // one side is zero) then we know the result.  Also, if both sides are zero
        // then we know the result.  The following logic encodes this.
        int lhs3_sgn = ExactFloat.Sign(dab), rhs3_sgn = -ExactFloat.Sign(dbc);
        int lhs2_sgn = Math.Max(-1, Math.Min(1, lhs3_sgn - rhs3_sgn));
        if (lhs2_sgn == 0 && lhs3_sgn != 0)
        {
            // Both sides of (3) have the same non-zero sign, so square both sides.
            // If both sides were negative then invert the result.
            lhs2_sgn = ExactFloat.Sign(cab2 - abc2) * lhs3_sgn;
        }
        // Now if the two sides of (2) have different signs then we know the result
        // of this entire function.
        int rhs2_sgn = -ExactFloat.Sign(dca);
        int result = Math.Max(-1, Math.Min(1, lhs2_sgn - rhs2_sgn));
        if (result == 0 && lhs2_sgn != 0)
        {
            // Both sides of (2) have the same non-zero sign, so square both sides.
            // (If both sides were negative then we invert the result below.)
            // This gives
            //
            //        |C|^2 dAB^2 + |A|^2 dBC^2 + 2 |A| |C| dAB dBC > |B|^2 dCA^2
            //
            // This expression still has square roots (|A| and |C|), so we rewrite as
            //
            // (4)    2 |A| |C| dAB dBC > |B|^2 dCA^2 - |C|^2 dAB^2 - |A|^2 dBC^2 .
            //
            // Again, if the two sides have different signs then we know the result.
            int lhs4_sgn = ExactFloat.Sign(dab) * ExactFloat.Sign(dbc);
            var rhs4 = bca2 - cab2 - abc2;
            result = Math.Max(-1, Math.Min(1, lhs4_sgn - ExactFloat.Sign(rhs4)));
            if (result == 0 && lhs4_sgn != 0)
            {
                // Both sides of (4) have the same non-zero sign, so square both sides.
                // If both sides were negative then invert the result.
                result = MathM.Sign(4 * abc2.Value * cab2.Value - rhs4.Value * rhs4.Value) * lhs4_sgn;
            }
            // Correct the sign if both sides of (2) were negative.
            result *= lhs2_sgn;
        }
        // If the sign of triangle ABC is negative, then we have computed -Z and the
        // result should be negated.
        return abc_sign * result;
    }

    // Like Sign, except this method does not use symbolic perturbations when
    // the input points are exactly coplanar with the origin (i.e., linearly
    // dependent).  Clients should never use this method, but it is useful here in
    // order to implement the combined pedestal/axis-aligned perturbation scheme
    // used by some methods (such as EdgeCircumcenterSign).
    private static int UnperturbedSign(S2Point a, S2Point b, S2Point c)
    {
        int sign = TriageSign(a, b, c, a.CrossProd(b));
        if (sign == 0) sign = ExpensiveSign(a, b, c, false /*perturb*/);
        return sign;
    }

    // Given arguments such that ExactEdgeCircumcenterSign(x0, x1, a, b, c) == 0,
    // returns the value of Sign(X0, X1, Z) (where Z is the circumcenter of
    // triangle ABC) after symbolic perturbations are taken into account.  The
    // result is zero only if X0 == X1, A == B, B == C, or C == A.  (It is nonzero
    // if these pairs are exactly proportional to each other but not equal.)
    private static int SymbolicEdgeCircumcenterSign(S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c)
    {
        // We use the same perturbation strategy as SymbolicCompareDistances.  Note
        // that pedestal perturbations of X0 and X1 do not affect the result,
        // because Sign(X0, X1, Z) does not change when its arguments are scaled
        // by a positive factor.  Therefore we only need to consider A, B, C.
        // Suppose that A is the smallest lexicographically and therefore has the
        // largest perturbation.  This has the effect of perturbing the circumcenter
        // of ABC slightly towards A, and since the circumcenter Z was previously
        // exactly collinear with edge X, this implies that after the perturbation
        // Sign(X0, X1, Z) == UnperturbedSign(X0, X1, A).  (We want the result
        // to be zero if X0, X1, and A are linearly dependent, rather than using
        // symbolic perturbations, because these perturbations are defined to be
        // much, much smaller than the pedestal perturbation of B and C that are
        // considered below.)
        //
        // If A is also exactly collinear with edge X, then we move on to the next
        // smallest point lexicographically out of {B, C}.  It is easy to see that
        // as long as A, B, C are all distinct, one of these three Sign calls
        // will be nonzero, because if A, B, C are all distinct and collinear with
        // edge X then their circumcenter Z coincides with the normal of X, and
        // therefore Sign(X0, X1, Z) is nonzero.
        //
        // This function could be extended to handle the case where X0 and X1 are
        // linearly dependent as follows.  First, suppose that every point has both
        // a pedestal peturbation as described above, and also the three
        // axis-aligned perturbations described in the "Simulation of Simplicity"
        // paper, where all pedestal perturbations are defined to be much, much
        // larger than any axis-aligned perturbation.  Note that since pedestal
        // perturbations have no effect on Sign, we can use this model for *all*
        // the S2 predicates, which ensures that all the various predicates are
        // fully consistent with each other.
        //
        // With this model, the strategy described above yields the correct result
        // unless X0 and X1 are exactly linearly dependent.  When that happens, then
        // no perturbation (pedestal or axis-aligned) of A,B,C affects the result,
        // and no pedestal perturbation of X0 or X1 affects the result, therefore we
        // need to consider the smallest axis-aligned perturbation of X0 or X1.  The
        // first perturbation that makes X0 and X1 linearly independent yields the
        // result.  Supposing that X0 < X1, this is the perturbation of X0[2] unless
        // both points are multiples of [0, 0, 1], in which case it is the
        // perturbation of X0[1].  The sign test can be implemented by computing the
        // perturbed cross product of X0 and X1 and taking the dot product with the
        // exact value of Z.  For example if X0[2] is perturbed, the perturbed cross
        // product is proportional to (0, 0, 1) x X1 = (-X1[1], x1[0], 0).  Note
        // that if the dot product with Z is exactly zero, then it is still
        // necessary to fall back to pedestal perturbations of A, B, C, but one of
        // these perturbations is now guaranteed to succeed.

        // If any two triangle vertices are equal, the result is zero.
        if (a == b || b == c || c == a) return 0;

        // Sort A, B, C in lexicographic order.
        S2Point buff;
        if (b < a) { buff = a; a = b; b = buff; }
        if (c < b) { buff = c; c = b; b = buff; }
        if (b < a) { buff = b; b = a; a = buff; }

        // Now consider the perturbations in decreasing order of size.
        int sign = UnperturbedSign(x0, x1, a);
        if (sign != 0) return sign;
        sign = UnperturbedSign(x0, x1, b);
        if (sign != 0) return sign;
        return UnperturbedSign(x0, x1, c);
    }

    private static Excluded TriageVoronoiSiteExclusion(S2Point a, S2Point b, S2Point x0, S2Point x1, double r2)
    {
        // Define the "coverage disc" of a site S to be the disc centered at S with
        // radius r (i.e., squared chord angle length r2).  Similarly, define the
        // "coverage interval" of S along an edge X to be the intersection of X with
        // the coverage disc of S.  The coverage interval can be represented as the
        // point at the center of the interval and an angle that measures the
        // semi-width or "radius" of the interval.
        //
        // To test whether site A excludes site B along the input edge X, we test
        // whether the coverage interval of A contains the coverage interval of B.
        // Let "ra" and "rb" be the radii (semi-widths) of the two intervals, and
        // let "d" be the angle between their center points.  Then "a" properly
        // contains "b" if (ra - rb > d), and "b" contains "a" if (rb - ra > d).
        // Note that only one of these conditions can be true.  Therefore we can
        // determine whether one site excludes the other by checking whether
        //
        // (1)   |rb - ra| > d
        //
        // and use the sign of (rb - ra) to determine which site is excluded.
        //
        // The actual code is based on the following.  Let A1 and B1 be the unit
        // vectors A and B scaled by cos(r) (these points are inside the sphere).
        // The planes perpendicular to OA1 and OA2 cut off two discs of radius r
        // around A and B.  Now consider the two lines (inside the sphere) where
        // these planes intersect the plane containing the input edge X, and let A2
        // and B2 be the points on these lines that are closest to A and B.  The
        // coverage intervals of A and B can be represented as an interval along
        // each of these lines, centered at A2 and B2.  Let P1 and P2 be the
        // endpoints of the coverage interval for A, and let Q1 and Q2 be the
        // endpoints of the coverage interval for B.  We can view each coverage
        // interval as either a chord through the sphere's interior, or as a segment
        // of the original edge X (by projecting the chord onto the sphere's
        // surface).
        //
        // To check whether B's interval is contained by A's interval, we test
        // whether both endpoints of B's interval (Q1 and Q2) are contained by A's
        // interval.  E.g., we could test whether Qi.DotProd(A2) > A2.Norm2.
        //
        // However rather than constructing the actual points A1, A2, and so on, it
        // turns out to be more efficient to compute the sines and cosines
        // ("components") of the various angles and then use trigonometric
        // identities.  Predicate (1) can be expressed as
        //
        //      |sin(rb - ra)| > sin(d)
        //
        // provided that |d| <= Pi/2 (which must be checked), and then expanded to
        //
        // (2)  |sin(rb) cos(ra) - sin(ra) cos(rb)| > sin(d) .
        //
        // The components of the various angles can be expressed using dot and cross
        // products based on the construction above:
        //
        //   sin(ra) = sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2) / |aXn|
        //   cos(ra) = cos(r) |a| |n| / |aXn|
        //   sin(rb) = sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) / |bXn|
        //   cos(rb) = cos(r) |b| |n| / |bXn|
        //   sin(d)  = (aXb).n |n| / (|aXn| |bXn|)
        //   cos(d)  = (aXn).(bXn) / (|aXn| |bXn|)
        //
        // Also, the squared chord length r2 is equal to 4 * sin^2(r / 2), which
        // yields the following relationships:
        //
        //   sin(r)  = sqrt(r2 (1 - r2 / 4))
        //   cos(r)  = 1 - r2 / 2
        //
        // We then scale both sides of (2) by |aXn| |bXn| / |n| (in order to
        // minimize the number of calculations and to avoid divisions), which gives:
        //
        //    cos(r) ||a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) -
        //            |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)| > (aXb).n
        //
        // Furthermore we can substitute |a| = |b| = 1 (as long as this is taken
        // into account in the error bounds), yielding
        //
        // (3)   cos(r) |sqrt(sin^2(r) |n|^2 - |b.n|^2) -
        //               sqrt(sin^2(r) |n|^2 - |a.n|^2)| > (aXb).n
        //
        // The code below is more complicated than this because many expressions
        // have been modified for better numerical stability.  For example, dot
        // products between unit vectors are computed using (x - y).DotProd(x + y),
        // and the dot product between a point P and the normal N of an edge X is
        // measured using (P - Xi).DotProd(N) where Xi is the endpoint of X that is
        // closer to P.

        S2Point n = (x0 - x1).CrossProd(x0 + x1);  // 2 * x0.CrossProd(x1)
        var n2 = n.Norm2();
        var n1 = Math.Sqrt(n2);
        // This factor is used in the error terms of dot products with "n" below.
        var Dn_error = ((3.5 + 2 * Math.Sqrt(3)) * n1 + 32 * Math.Sqrt(3) * DBL_ERR) * TT_ERR;

        var cos_r = 1 - 0.5 * r2;
        var sin2_r = r2 * (1 - 0.25 * r2);
        var n2sin2_r = n2 * sin2_r;

        // "ra" and "rb" denote sin(ra) and sin(rb) after the scaling above.
        var aDn = (a - GetClosestVertex(a, x0, x1, out var ax2)).DotProd(n);
        var aDn2 = aDn * aDn;
        var aDn_error = Dn_error * Math.Sqrt(ax2);
        var ra2 = n2sin2_r - aDn2;
        var ra2_error = (8 * DBL_ERR + 4 * TT_ERR) * aDn2 +
            (2 * Math.Abs(aDn) + aDn_error) * aDn_error + 6 * TT_ERR * n2sin2_r;
        // This is the minimum possible value of ra2, which is used to bound the
        // derivative of sqrt(ra2) in computing ra_error below.
        var min_ra2 = ra2 - ra2_error;
        if (min_ra2 < 0) return Excluded.UNCERTAIN;
        var ra = Math.Sqrt(ra2);
        // Includes the ra2 subtraction error above.
        var ra_error = 1.5 * TT_ERR * ra + 0.5 * ra2_error / Math.Sqrt(min_ra2);

        var bDn = (b - GetClosestVertex(b, x0, x1, out var bx2)).DotProd(n);
        var bDn2 = bDn * bDn;
        var bDn_error = Dn_error * Math.Sqrt(bx2);
        var rb2 = n2sin2_r - bDn2;
        var rb2_error = (8 * DBL_ERR + 4 * TT_ERR) * bDn2 +
            (2 * Math.Abs(bDn) + bDn_error) * bDn_error + 6 * TT_ERR * n2sin2_r;
        var min_rb2 = rb2 - rb2_error;
        if (min_rb2 < 0) return Excluded.UNCERTAIN;
        var rb = Math.Sqrt(rb2);
        // Includes the rb2 subtraction error above.
        var rb_error = 1.5 * TT_ERR * rb + 0.5 * rb2_error / Math.Sqrt(min_rb2);

        // The sign of LHS(3) before taking the absolute value determines which site
        // may be excluded by the other.  If it is positive then A may be excluded,
        // and if it is negative then B may be excluded.
        var lhs3 = cos_r * (rb - ra);
        var abs_lhs3 = Math.Abs(lhs3);
        var lhs3_error = cos_r * (ra_error + rb_error) + 3 * TT_ERR * abs_lhs3;

        // Now we evaluate the RHS of (3), which is proportional to sin(d).
        S2Point aXb = (a - b).CrossProd(a + b);  // 2 * a.CrossProd(b)
        var aXb1 = aXb.Norm();
        var sin_d = 0.5 * aXb.DotProd(n);
        var sin_d_error = (4 * DBL_ERR + (2.5 + 2 * Math.Sqrt(3)) * TT_ERR) * aXb1 * n1 +
            16 * Math.Sqrt(3) * DBL_ERR * TT_ERR * (aXb1 + n1);

        // If LHS(3) is definitely less than RHS(3), neither site excludes the other.
        var result = abs_lhs3 - sin_d;
        var result_error = lhs3_error + sin_d_error;
        if (result < -result_error) return Excluded.NEITHER;

        // d < 0 means that when AB is projected onto the great circle through X0X1,
        // it proceeds in the opposite direction as X0X1.  Normally we have d > 0
        // since GetVoronoiSiteExclusion requires d(A,X0) < d(B,X0) (corresponding to
        // the fact that sites are processed in order of increasing distance from X0).
        //
        // However when edge X is long and/or the snap radius "r" is large, there
        // are two cases where where d < 0 can occur:
        //
        // 1. d(A,X0) > Pi/2 and d(B,X1) < Pi/2 or the symmetric case (swap < and >).
        //    This can only happen when d(X0,X1) + r > Pi/2.  Note that {A,B} may
        //    project to the interior of edge X or beyond its endpoints.  In this
        //    case A is kept iff d(A,X0) < Pi/2, and B is kept otherwise.  Note that
        //    the sign of (rb - ra) is not sufficient to determine which point is
        //    kept in the situation where d(X0,X1) + r > Pi.
        //
        // 2. A is beyond endpoint X0, B is beyond endpoint X1, and AB wraps around
        //    the sphere in the opposite direction from edge X.  This case can only
        //    happen when d(X0,X1) + r > Pi.  Here each site is closest to one
        //    endpoint of X and neither excludes the other.
        //
        // The algorithm that handles both cases is to keep A if d(A,X0) < Pi/2 and
        // to keep B if d(B,X1) < Pi/2.  (One of these conditions is always true.)
        if (sin_d < -sin_d_error)
        {
            // Site A is kept if ca < 0 and excluded if ca > 0.
            var r90 = S1ChordAngle.Right.Length2;
            int ca = TriageCompareCosDistance(a, x0, r90);
            int cb = TriageCompareCosDistance(b, x1, r90);
            if (ca < 0 && cb < 0) return Excluded.NEITHER;
            if (ca <= 0 && cb <= 0) return Excluded.UNCERTAIN;  // One or both kept?
                                                                 // Since either ca or cb is 1, we know the result even if the distance
                                                                 // comparison for the other site was uncertain.
            MyDebug.Assert(ca <= 0 || cb <= 0);
            return (ca > 0) ? Excluded.FIRST : Excluded.SECOND;
        }
        if (sin_d <= sin_d_error) return Excluded.UNCERTAIN;

        // Otherwise, before proceeding further we need to check that |d| <= Pi/2.
        // In fact, |d| < Pi/2 is enough because of the requirement that r < Pi/2.
        // The following expression represents cos(d) after scaling; it is
        // equivalent to (aXn).(bXn) but has about 30% less error.
        var cos_d = a.DotProd(b) * n2 - aDn * bDn;
        var cos_d_error =
            ((8 * DBL_ERR + 5 * TT_ERR) * Math.Abs(aDn) + aDn_error) * Math.Abs(bDn) +
            (Math.Abs(aDn) + aDn_error) * bDn_error + (8 * DBL_ERR + 8 * TT_ERR) * n2;
        if (cos_d <= -cos_d_error) return Excluded.NEITHER;

        // Potential optimization: if the sign of cos(d) is uncertain, then instead
        // we could check whether cos(d) >= cos(r).  Unfortunately this is fairly
        // expensive since it requires computing denominator |aXn||bXn| of cos(d)
        // and the associated error bounds.  In any case this case is relatively
        // rare so it seems better to punt.
        if (cos_d < cos_d_error) return Excluded.UNCERTAIN;

        // Now we can finish checking the results of predicate (3).
        if (result <= result_error) return Excluded.UNCERTAIN;
        MyDebug.Assert(abs_lhs3 > lhs3_error);
        return (lhs3 > 0) ? Excluded.FIRST : Excluded.SECOND;
    }

    private static Excluded ExactVoronoiSiteExclusion(Vector3<ExactFloat> a, Vector3<ExactFloat> b, Vector3<ExactFloat> x0, Vector3<ExactFloat> x1, ExactFloat r2) =>
        ExactVoronoiSiteExclusion(a, b, x0, x1, r2.Value);
    private static Excluded ExactVoronoiSiteExclusion(Vector3<ExactFloat> a, Vector3<ExactFloat> b, Vector3<ExactFloat> x0, Vector3<ExactFloat> x1, decimal r2)
    {
        MyDebug.Assert(!ArePointsAntipodal(x0, x1));

        // Recall that one site excludes the other if
        //
        // (1)  |sin(rb - ra)| > sin(d)
        //
        // and that the sign of (rb - ra) determines which site is excluded (see the
        // comments in TriageVoronoiSiteExclusion).  To evaluate this using exact
        // arithmetic, we expand this to the same predicate as before:
        //
        // (2)    cos(r) ||a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) -
        //                |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)| > (aXb).n
        //
        // We also need to verify that d <= Pi/2, which is implemented by checking
        // that sin(d) >= 0 and cos(d) >= 0.
        //
        // To eliminate the square roots we use the standard technique of
        // rearranging the inequality to isolate at least one square root and then
        // squaring both sides.  We need to repeat this process twice in order to
        // eliminate all the square roots, which leads to a polynomial predicate of
        // degree 20 in the input arguments (i.e., degree 4 in each of "a", "b",
        // "x0", "x1", and "r2").


        // Before squaring we need to check the sign of each side.  If the RHS of
        // (2) is negative (corresponding to sin(d) < 0), then we need to apply the
        // logic in TriageVoronoiSiteExclusion.
        var n = x0.CrossProd(x1);
        var rhs2 = a.CrossProd(b).DotProd(n);
        int rhs2_sgn = ExactFloat.Sign(rhs2);
        if (rhs2_sgn < 0)
        {
            // This is the d < 0 case.  See comments in TriageVoronoiSiteExclusion.
            var r90 = S1ChordAngle.Right.Length2;
            int ca = ExactCompareDistance(a, x0, (decimal)r90);
            int cb = ExactCompareDistance(b, x1, (decimal)r90);
            if (ca < 0 && cb < 0) return Excluded.NEITHER;
            MyDebug.Assert(ca != 0 && cb != 0, "This is guaranteed since d < 0.");
            MyDebug.Assert(ca < 0 || cb < 0, "At least one site must be kept.");
            return (ca > 0) ? Excluded.FIRST : Excluded.SECOND;
        }

        // We also check that cos(d) >= 0.  Given what else we need to compute, it
        // is cheaper use the identity (aXn).(bXn) = (a.b) |n|^2 - (a.n)(b.n) .
        var n2 = n.Norm2();
        var aDn = a.DotProd(n);
        var bDn = b.DotProd(n);
        var cos_d = a.DotProd(b) * n2 - aDn * bDn;
        if (ExactFloat.Sign(cos_d) < 0) return Excluded.NEITHER;

        // Otherwise we continue evaluating the LHS of (2), defining
        //    sa = |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)
        //    sb = |a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) .
        // The sign of the LHS of (2) (before taking the absolute value) determines
        // which coverage interval is larger and therefore which site is potentially
        // being excluded.
        var a2 = a.Norm2();
        var b2 = b.Norm2();
        var n2sin2_r = r2 * (1 - 0.25M * r2) * n2.Value;
        var sa2 = b2.Value * (n2sin2_r * a2.Value - aDn.Value * aDn.Value);
        var sb2 = a2.Value * (n2sin2_r * b2.Value - bDn.Value * bDn.Value);
        int lhs2_sgn = MathM.Sign(sb2 - sa2);

        if (lhs2_sgn == 0)
        {
            // If the RHS of (2) is zero as well (i.e., d == 0) then both sites are
            // equidistant from every point on edge X.  This case requires symbolic
            // perturbations, but it should already have been handled in
            // GetVoronoiSiteExclusion() (see the call to CompareDistances).
            MyDebug.Assert(rhs2_sgn > 0);
            return Excluded.NEITHER;
        }
        // Next we square both sides of (2), yielding
        //
        //      cos^2(r) (sb^2 + sa^2 - 2 sa sb) > (aXb.n)^2
        //
        // which can be rearranged to give
        //
        // (3)  cos^2(r) (sb^2 + sa^2) - (aXb.n)^2 > 2 cos^2(r) sa sb .
        //
        // The RHS of (3) is always non-negative, but we still need to check the
        // sign of the LHS.
        var cos_r = 1 - 0.5M * r2;
        var cos2_r = cos_r * cos_r;
        var lhs3 = cos2_r * (sa2 + sb2) - (rhs2 * rhs2).Value;
        if (decimal.Sign(lhs3) < 0) return Excluded.NEITHER;

        // Otherwise we square both sides of (3) to obtain:
        //
        // (4)  LHS(3)^2  >  4 cos^4(r) sa^2 sb^2
        var lhs4 = lhs3 * lhs3;
        var rhs4 = 4 * cos2_r * cos2_r * sa2 * sb2;
        int result = decimal.Sign(lhs4 - rhs4);
        if (result < 0) return Excluded.NEITHER;
        if (result == 0)
        {
            // We have |rb - ra| = d and d > 0.  This implies that one coverage
            // interval contains the other, but not properly: the two intervals share
            // a common endpoint.  The distance from each site to that point is
            // exactly "r", therefore we need to use symbolic perturbations.  Recall
            // that site A is considered closer to an equidistant point if and only if
            // A > B.  Therefore if (rb > ra && A > B) or (ra > rb && B > A) then each
            // site is closer to at least one point and neither site is excluded.
            //
            // Ideally this logic would be in a separate SymbolicVoronoiSiteExclusion
            // method for better testing, but this is not convenient because it needs
            // lhs_sgn (which requires exact arithmetic to compute).
            if ((lhs2_sgn > 0) == (a > b)) return Excluded.NEITHER;
        }
        // At this point we know that one of the two sites is excluded.  The sign of
        // the LHS of (2) (before the absolute value) determines which one.
        return (lhs2_sgn > 0) ? Excluded.FIRST : Excluded.SECOND;
    }

#if s2debug

    // The following functions are not part of the public API.  Currently they are
    // only used internally for testing purposes.

    public static int StableSign_Test(S2Point a, S2Point b, S2Point c) =>
        StableSign(a, b, c);

    public static int ExactSign_Test(S2Point a, S2Point b, S2Point c, bool perturb) =>
        ExactSign(a, b, c, perturb);

    public static int SymbolicallyPerturbedSign_Test(S2Point a, S2Point b, S2Point c, S2Point b_cross_c) =>
        SymbolicallyPerturbedSign(a, b, c, b_cross_c);

    public static int TriageCompareCosDistances_Test(S2Point x, S2Point a, S2Point b) =>
        TriageCompareCosDistances(x, a, b);

    public static int TriageCompareSin2Distances_Test(S2Point x, S2Point a, S2Point b) =>
        TriageCompareSin2Distances(x, a, b);

    public static int ExactCompareDistances_Test(Vector3<ExactFloat> x, Vector3<ExactFloat> a, Vector3<ExactFloat> b) =>
        ExactCompareDistances(x, a, b);

    public static int SymbolicCompareDistances_Test(S2Point x, S2Point a, S2Point b) =>
        SymbolicCompareDistances(x, a, b);

    public static int TriageCompareSin2Distance_Test(S2Point x, S2Point y, double r2) =>
        TriageCompareSin2Distance(x, y, r2);

    public static int TriageCompareCosDistance_Test(S2Point x, S2Point y, double r2) =>
        TriageCompareCosDistance(x, y, r2);

    public static int ExactCompareDistance_Test(Vector3<ExactFloat> x, Vector3<ExactFloat> y, ExactFloat r2) =>
        ExactCompareDistance_Test(x, y, r2.Value);
    public static int ExactCompareDistance_Test(Vector3<ExactFloat> x, Vector3<ExactFloat> y, decimal r2) =>
        ExactCompareDistance(x, y, r2);

    public static int TriageCompareEdgeDistance_Test(S2Point x, S2Point a0, S2Point a1, double r2) =>
        TriageCompareEdgeDistance(x, a0, a1, r2);

    public static int ExactCompareEdgeDistance_Test(S2Point x, S2Point a0, S2Point a1, S1ChordAngle r) =>
        ExactCompareEdgeDistance(x, a0, a1, r);

    public static int TriageCompareEdgeDirections_Test(S2Point a0, S2Point a1, S2Point b0, S2Point b1) =>
        TriageCompareEdgeDirections(a0, a1, b0, b1);

    public static int ExactCompareEdgeDirections_Test(Vector3<ExactFloat> a0, Vector3<ExactFloat> a1, Vector3<ExactFloat> b0, Vector3<ExactFloat> b1) =>
        ExactCompareEdgeDirections(a0, a1, b0, b1);

    public static int TriageEdgeCircumcenterSign_Test(S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c, int abc_sign) =>
        TriageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign);

    public static int ExactEdgeCircumcenterSign_Test(Vector3<ExactFloat> x0, Vector3<ExactFloat> x1, Vector3<ExactFloat> a, Vector3<ExactFloat> b, Vector3<ExactFloat> c, int abc_sign) =>
        ExactEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign);

    public static int SymbolicEdgeCircumcenterSign_Test(S2Point x0, S2Point x1, S2Point a_arg, S2Point b_arg, S2Point c_arg) =>
        SymbolicEdgeCircumcenterSign(x0, x1, a_arg, b_arg, c_arg);

    public static Excluded TriageVoronoiSiteExclusion_Test(S2Point a, S2Point b, S2Point x0, S2Point x1, double r2) =>
        TriageVoronoiSiteExclusion(a, b, x0, x1, r2);

    public static Excluded ExactVoronoiSiteExclusion_Test(Vector3<ExactFloat> a, Vector3<ExactFloat> b, Vector3<ExactFloat> x0, Vector3<ExactFloat> x1, ExactFloat r2) =>
        ExactVoronoiSiteExclusion_Test(a, b, x0, x1, r2.Value);
    public static Excluded ExactVoronoiSiteExclusion_Test(Vector3<ExactFloat> a, Vector3<ExactFloat> b, Vector3<ExactFloat> x0, Vector3<ExactFloat> x1, decimal r2) =>
        ExactVoronoiSiteExclusion(a, b, x0, x1, r2);

#endif
}
