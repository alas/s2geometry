namespace S2Geometry;

using System.Text;

public class S2PredicatesTests
{
    private const int consistency_iters = 5000; // Number of iterations for precision consistency tests
    private static readonly string[] kPrecisionNames = new[] { "double", "double", "exact", "symbolic" };

    // If `sizeof(long double) == sizeof(double)`, then we will never do
    // calculations with `long double` and instead fall back to exact.
    private const Precision kLongDoublePrecision = S2Pred.kHasLongDouble ? Precision.LONG_DOUBLE : Precision.EXACT;

    [Fact]
    internal void Test_epsilon_for_digits_recursion()
    {
        Assert.Equal(1.0, EpsilonForDigits(0));
        Assert.Equal(MathUtils.Ldexp(1.0, -24), EpsilonForDigits(24));
        Assert.Equal(MathUtils.Ldexp(1.0, -53), EpsilonForDigits(53));
        Assert.Equal(MathUtils.Ldexp(1.0, -64), EpsilonForDigits(64));
        Assert.Equal(MathUtils.Ldexp(1.0, -106), EpsilonForDigits(106));
        Assert.Equal(MathUtils.Ldexp(1.0, -113), EpsilonForDigits(113));
    }

    [Fact]
    internal void Test_Sign_CollinearPoints()
    {
        // The following points happen to be *exactly collinear* along a line that it
        // approximate tangent to the surface of the unit sphere.  In fact, C is the
        // exact midpoint of the line segment AB.  All of these points are close
        // enough to unit length to satisfy S2.IsUnitLength().
        S2Point a = new(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
        S2Point b = new(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
        S2Point c = new(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
        Assert.Equal(c - a, b - c);
        Assert.NotEqual(0, S2Pred.Sign(a, b, c));
        Assert.Equal(S2Pred.Sign(a, b, c), S2Pred.Sign(b, c, a));
        Assert.Equal(S2Pred.Sign(a, b, c), -S2Pred.Sign(c, b, a));

        // The points "x1" and "x2" are exactly proportional, i.e. they both lie
        // on a common line through the origin.  Both points are considered to be
        // normalized, and in fact they both satisfy (x == x.Normalize()).
        // Therefore the triangle (x1, x2, -x1) consists of three distinct points
        // that all lie on a common line through the origin.
        S2Point x1 = new(0.99999999999999989, 1.4901161193847655e-08, 0);
        S2Point x2 = new(1, 1.4901161193847656e-08, 0);
        Assert.Equal(x1, x1.Normalize());
        Assert.Equal(x2, x2.Normalize());
        Assert.NotEqual(0, S2Pred.Sign(x1, x2, -x1));
        Assert.Equal(S2Pred.Sign(x1, x2, -x1), S2Pred.Sign(x2, -x1, x1));
        Assert.Equal(S2Pred.Sign(x1, x2, -x1), -S2Pred.Sign(-x1, x2, x1));

        // Here are two more points that are distinct, exactly proportional, and
        // that satisfy (x == x.Normalize()).
        S2Point x3 = new S2Point(1, 1, 1).Normalize();
        S2Point x4 = 0.99999999999999989 * x3;
        Assert.Equal(x3, x3.Normalize());
        Assert.Equal(x4, x4.Normalize());
        Assert.NotEqual(x3, x4);
        Assert.NotEqual(0, S2Pred.Sign(x3, x4, -x3));

        // The following two points demonstrate that Normalize() is not idempotent,
        // i.e. y0.Normalize() != y0.Normalize().Normalize().  Both points satisfy
        // S2.IsNormalized(), though, and the two points are exactly proportional.
        S2Point y0 = new(1, 1, 0);
        S2Point y1 = y0.Normalize();
        S2Point y2 = y1.Normalize();
        Assert.NotEqual(y1, y2);
        Assert.Equal(y2, y2.Normalize());
        Assert.NotEqual(0, S2Pred.Sign(y1, y2, -y1));
        Assert.Equal(S2Pred.Sign(y1, y2, -y1), S2Pred.Sign(y2, -y1, y1));
        Assert.Equal(S2Pred.Sign(y1, y2, -y1), -S2Pred.Sign(-y1, y2, y1));
    }

    [Fact]
    internal void Test_SignTest_StressTest()
    {
        // The run time of this test is *cubic* in the parameter below.
        const int kNumPointsPerCircle = 17;

        // This test is randomized, so it is beneficial to run it several times.
        for (int iter = 0; iter < 3; ++iter)
        {
            // The most difficult great circles are the ones in the X-Y, Y-Z, and X-Z
            // planes, for two reasons.  First, when one or more coordinates are close
            // to zero then the perturbations can be much smaller, since floating
            // point numbers are spaced much more closely together near zero.  (This
            // tests the handling of things like underflow.)  The second reason is
            // that most of the cases of SymbolicallyPerturbedSign() can only be
            // reached when one or more input point coordinates are zero.
            TestGreatCircle(new(1, 0, 0), new(0, 1, 0), kNumPointsPerCircle);
            TestGreatCircle(new(1, 0, 0), new(0, 0, 1), kNumPointsPerCircle);
            TestGreatCircle(new(0, -1, 0), new(0, 0, 1), kNumPointsPerCircle);

            // This tests a great circle where at least some points have X, Y, and Z
            // coordinates with exactly the same mantissa.  One useful property of
            // such points is that when they are scaled (e.g. multiplying by 1+eps),
            // all such points are exactly collinear with the origin.
            TestGreatCircle(new(1 << 25, 1, -8), new(-4, -(1 << 20), 1), kNumPointsPerCircle);
        }
    }

    [Fact]
    internal void Test_StableSignTest_FailureRate()
    {
        // Verify that StableSign() is able to handle most cases where the three
        // points are as collinear as possible.  (For reference, TriageSign() fails
        // virtually 100% of the time on this test.)
        //
        // Note that the failure rate *decreases* as the points get closer together,
        // and the decrease is approximately linear.  For example, the failure rate
        // is 0.4% for collinear points spaced 1km apart, but only 0.0004% for
        // collinear points spaced 1 meter apart.

        Assert.True(StableSignTest.GetFailureRate(1.0) < 0.01);  //  1km spacing: <  1% (actual 0.4%)
        Assert.True(StableSignTest.GetFailureRate(10.0) < 0.1);  // 10km spacing: < 10% (actual 4%)
    }

    [Fact]
    internal void Test_Sign_SymbolicPerturbationCodeCoverage()
    {
        // The purpose of this test is simply to get code coverage of
        // SymbolicallyPerturbedSign().  Let M_1, M_2, ... be the sequence of
        // submatrices whose determinant sign is tested by that function.  Then the
        // i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
        //
        //    det(M) = 0
        //    det(M_j) = 0 for j < i
        //    det(M_i) != 0
        //    A < B < C in lexicographic order.
        //
        // I checked that reversing the sign of any of the "return" statements in
        // SymbolicallyPerturbedSign() will cause this test to fail.

        // det(M_1) = b0*c1 - b1*c0
        CheckSymbolicSign(1, new(-3, -1, 0), new(-2, 1, 0), new(1, -2, 0));

        // det(M_2) = b2*c0 - b0*c2
        CheckSymbolicSign(1, new(-6, 3, 3), new(-4, 2, -1), new(-2, 1, 4));

        // det(M_3) = b1*c2 - b2*c1
        CheckSymbolicSign(1, new(0, -1, -1), new(0, 1, -2), new(0, 2, 1));
        // From this point onward, B or C must be zero, or B is proportional to C.

        // det(M_4) = c0*a1 - c1*a0
        CheckSymbolicSign(1, new(-1, 2, 7), new(2, 1, -4), new(4, 2, -8));

        // det(M_5) = c0
        CheckSymbolicSign(1, new(-4, -2, 7), new(2, 1, -4), new(4, 2, -8));

        // det(M_6) = -c1
        CheckSymbolicSign(1, new(0, -5, 7), new(0, -4, 8), new(0, -2, 4));

        // det(M_7) = c2*a0 - c0*a2
        CheckSymbolicSign(1, new(-5, -2, 7), new(0, 0, -2), new(0, 0, -1));

        // det(M_8) = c2
        CheckSymbolicSign(1, new(0, -2, 7), new(0, 0, 1), new(0, 0, 2));
        // From this point onward, C must be zero.

        // det(M_9) = a0*b1 - a1*b0
        CheckSymbolicSign(1, new(-3, 1, 7), new(-1, -4, 1), S2Point.Empty);

        // det(M_10) = -b0
        CheckSymbolicSign(1, new(-6, -4, 7), new(-3, -2, 1), S2Point.Empty);

        // det(M_11) = b1
        CheckSymbolicSign(-1, new(0, -4, 7), new(0, -2, 1), S2Point.Empty);

        // det(M_12) = a0
        CheckSymbolicSign(-1, new(-1, -4, 5), new(0, 0, -3), S2Point.Empty);

        // det(M_13) = 1
        CheckSymbolicSign(1, new(0, -4, 5), new(0, 0, -5), S2Point.Empty);
    }

    [Fact]
    internal void Test_CompareDistances_Coverage()
    {
        // This test attempts to exercise all the code paths in all precisions.

        // Test TriageCompareSin2Distances.
        TestCompareDistances(new(1, 1, 1), new(1, 1 - 1e-15, 1), new(1, 1, 1 + 2e-15), -1, Precision.DOUBLE, Sin2Distances);
        TestCompareDistances(new(1, 1, 0), new(1, 1 - 1e-15, 1e-21), new(1, 1 - 1e-15, 0), 1, Precision.DOUBLE, Sin2Distances);
        TestCompareDistances(new(2, 0, 0), new(2, -1, 0), new(2, 1, 1e-8), -1, kLongDoublePrecision, Sin2Distances);
        TestCompareDistances(new(2, 0, 0), new(2, -1, 0), new(2, 1, 1e-100), -1, Precision.EXACT, Sin2Distances);
        TestCompareDistances(new(1, 0, 0), new(1, -1, 0), new(1, 1, 0), 1, Precision.SYMBOLIC, Sin2Distances);
        TestCompareDistances(new(1, 0, 0), new(1, 0, 0), new(1, 0, 0), 0, Precision.SYMBOLIC, Sin2Distances);

        // Test TriageCompareCosDistances.
        TestCompareDistances(new(1, 1, 1), new(1, -1, 0), new(-1, 1, 3e-15), 1, Precision.DOUBLE, CosDistances);
        TestCompareDistances(new(1, 0, 0), new(1, 1e-30, 0), new(-1, 1e-40, 0), -1, Precision.DOUBLE, CosDistances);
        TestCompareDistances(new(1, 1, 1), new(1, -1, 0), new(-1, 1, 3e-18), 1, kLongDoublePrecision, CosDistances);
        TestCompareDistances(new(1, 1, 1), new(1, -1, 0), new(-1, 1, 1e-100), 1, Precision.EXACT, CosDistances);
        TestCompareDistances(new(1, 1, 1), new(1, -1, 0), new(-1, 1, 0), -1, Precision.SYMBOLIC, CosDistances);
        TestCompareDistances(new(1, 1, 1), new(1, -1, 0), new(1, -1, 0), 0, Precision.SYMBOLIC, CosDistances);

        // Test TriageCompareSin2Distances using distances greater than 90 degrees.
        TestCompareDistances(new(1, 1, 0), new(-1, -1 + 1e-15, 0), new(-1, -1, 0), -1, Precision.DOUBLE, MinusSin2Distances);
        TestCompareDistances(new(-1, -1, 0), new(1, 1 - 1e-15, 0), new(1, 1 - 1e-15, 1e-21), 1, Precision.DOUBLE, MinusSin2Distances);
        TestCompareDistances(new(-1, -1, 0), new(2, 1, 0), new(2, 1, 1e-8), 1, kLongDoublePrecision, MinusSin2Distances);
        TestCompareDistances(new(-1, -1, 0), new(2, 1, 0), new(2, 1, 1e-30), 1, Precision.EXACT, MinusSin2Distances);
        TestCompareDistances(new(-1, -1, 0), new(2, 1, 0), new(1, 2, 0), -1, Precision.SYMBOLIC, MinusSin2Distances);
    }

    [Fact]
    internal void Test_CompareDistances_Consistency()
    {
        // This test chooses random point pairs that are nearly equidistant from a
        // target point, and then checks that the answer given by a method at one
        // level of precision is consistent with the answer given at the next higher
        // level of precision.
        //
        // The way the .cc file is structured, we can only do comparisons using a
        // specific precision if we also choose the specific distance calculation
        // method.  The code below checks that the Cos, Sin2, and MinusSin2 methods
        // are consistent across their entire valid range of inputs, and also
        // simulates the logic in CompareDistance that chooses which method to use
        // in order to gather statistics about how often each precision is needed.
        // (These statistics are only useful for coverage purposes, not benchmarks,
        // since the input points are chosen to be pathological worst cases.)
        TestCompareDistancesConsistency(new(1, 0, 0), new(0, -1, 0), new(0, 1, 0), CosDistances);
        var sin2_stats = new PrecisionStats();
        var cos_stats = new PrecisionStats();
        var minus_sin2_stats = new PrecisionStats();
        for (int iter = 0; iter < consistency_iters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point x = ChoosePoint();
            S2Point dir = ChoosePoint();
            S1Angle r = S1Angle.FromRadians(S2.M_PI_2 * Math.Pow(1e-30, S2Testing.Random.RandDouble()));
            if (S2Testing.Random.OneIn(2)) r = S1Angle.FromRadians(S2.M_PI_2) - r;
            if (S2Testing.Random.OneIn(2)) r = S1Angle.FromRadians(S2.M_PI_2) + r;
            S2Point a = S2.GetPointOnLine(x, dir, r);
            S2Point b = S2.GetPointOnLine(x, -dir, r);
            Precision prec = TestCompareDistancesConsistency(x, a, b, CosDistances);
            if (r.GetDegrees() >= 45 && r.GetDegrees() <= 135) cos_stats.Tally(prec);
            // The Sin2 method is only valid if both distances are less than 90
            // degrees, and similarly for the MinusSin2 method.  (In the actual
            // implementation these methods are only used if both distances are less
            // than 45 degrees or greater than 135 degrees respectively.)
            if (r.Radians < S2.M_PI_2 - 1e-14)
            {
                prec = TestCompareDistancesConsistency(x, a, b, Sin2Distances);
                if (r.GetDegrees() < 45)
                {
                    // Don't skew the statistics by recording degenerate inputs.
                    if (a == b)
                    {
                        Assert.Equal(Precision.SYMBOLIC, prec);
                    }
                    else
                    {
                        sin2_stats.Tally(prec);
                    }
                }
            }
            else if (r.Radians > S2.M_PI_2 + 1e-14)
            {
                prec = TestCompareDistancesConsistency(x, a, b, MinusSin2Distances);
                if (r.GetDegrees() > 135) minus_sin2_stats.Tally(prec);
            }
        }
    }

    [Fact]
    internal void Test_CompareDistance_Coverage()
    {
        // Test TriageCompareSin2Distance.
        TestCompareDistance(new(1, 1, 1), new(1, 1 - 1e-15, 1),
            S1ChordAngle.FromRadians(1e-15), -1, Precision.DOUBLE, Sin2Distance);
        TestCompareDistance(new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.FromRadians(S2.M_PI_4), -1, kLongDoublePrecision, Sin2Distance);
        TestCompareDistance(new(1, 1e-40, 0), new(1 + S2.DoubleEpsilon, 1e-40, 0),
            S1ChordAngle.FromRadians(0.9 * S2.DoubleEpsilon * 1e-40), 1, Precision.EXACT, Sin2Distance);
        TestCompareDistance(new(1, 1e-40, 0), new(1 + S2.DoubleEpsilon, 1e-40, 0),
            S1ChordAngle.FromRadians(1.1 * S2.DoubleEpsilon * 1e-40), -1, Precision.EXACT, Sin2Distance);
        TestCompareDistance(new(1, 0, 0), new(1 + S2.DoubleEpsilon, 0, 0),
            S1ChordAngle.Zero, 0, Precision.EXACT, Sin2Distance);

        // Test TriageCompareCosDistance.
        TestCompareDistance(new(1, 0, 0), new(1, 1e-8, 0),
            S1ChordAngle.FromRadians(1e-7), -1, Precision.DOUBLE, CosDistance);
        TestCompareDistance(new(1, 0, 0), new(-1, 1e-8, 0),
            S1ChordAngle.FromRadians(Math.PI - 1e-7), 1, Precision.DOUBLE, CosDistance);
        TestCompareDistance(new(1, 1, 0), new(1, -1 - 2 * S2.DoubleEpsilon, 0),
            S1ChordAngle.Right, 1, Precision.DOUBLE, CosDistance);
        TestCompareDistance(new(1, 1, 0), new(1, -1 - S2.DoubleEpsilon, 0),
            S1ChordAngle.Right, 1, kLongDoublePrecision, CosDistance);
        TestCompareDistance(new(1, 1, 0), new(1, -1, 1e-30),
            S1ChordAngle.Right, 0, Precision.EXACT, CosDistance);
        // The angle between these two points is exactly 60 degrees.
        TestCompareDistance(new(1, 1, 0), new(0, 1, 1),
            S1ChordAngle.FromLength2(1), 0, Precision.EXACT, CosDistance);
    }

    [Fact]
    internal void Test_CompareDistance_Consistency()
    {
        // This test chooses random inputs such that the distance between points X
        // and Y is very close to the threshold distance "r".  It then checks that
        // the answer given by a method at one level of precision is consistent with
        // the answer given at the next higher level of precision.  See also the
        // comments in the CompareDistances consistency test.
        var sin2_stats = new PrecisionStats();
        var cos_stats = new PrecisionStats();
        for (int iter = 0; iter < consistency_iters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point x = ChoosePoint();
            S2Point dir = ChoosePoint();
            S1Angle r = S1Angle.FromRadians(S2.M_PI_2 * Math.Pow(1e-30, S2Testing.Random.RandDouble()));
            if (S2Testing.Random.OneIn(2)) r = S1Angle.FromRadians(S2.M_PI_2) - r;
            if (S2Testing.Random.OneIn(5)) r = S1Angle.FromRadians(S2.M_PI_2) + r;
            S2Point y = S2.GetPointOnLine(x, dir, r);
            Precision prec = TestCompareDistanceConsistency(x, y, new S1ChordAngle(r), CosDistance);
            if (r.GetDegrees() >= 45) cos_stats.Tally(prec);
            if (r.Radians < S2.M_PI_2 - 1e-14)
            {
                prec = TestCompareDistanceConsistency(x, y, new S1ChordAngle(r), Sin2Distance);
                if (r.GetDegrees() < 45) sin2_stats.Tally(prec);
            }
        }
    }

    [Fact]
    internal void Test_CompareEdgeDistance_Coverage()
    {
        // Test TriageCompareLineSin2Distance.
        TestCompareEdgeDistance(
            new(1, 1e-10, 1e-15), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.FromRadians(1e-15 + S2.DoubleEpsilon), -1, Precision.DOUBLE);
        TestCompareEdgeDistance(
            new(1, 1, 1e-15), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.FromRadians(1e-15 + S2.DoubleEpsilon), -1, kLongDoublePrecision);
        TestCompareEdgeDistance(
            new(1, 1, 1e-40), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.FromRadians(1e-40), -1, Precision.EXACT);
        TestCompareEdgeDistance(
            new(1, 1, 0), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.Zero, 0, Precision.EXACT);

        // Test TriageCompareLineCos2Distance.
        TestCompareEdgeDistance(
            new(1e-15, 0, 1), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.FromRadians(S2.M_PI_2 - 1e-15 - 3 * S2.DoubleEpsilon),
            1, Precision.DOUBLE);
        TestCompareEdgeDistance(
            new(1e-15, 0, 1), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.FromRadians(S2.M_PI_2 - 1e-15 - S2.DoubleEpsilon),
            1, kLongDoublePrecision);
        TestCompareEdgeDistance(
            new(1e-40, 0, 1), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.Right, -1, Precision.EXACT);
        TestCompareEdgeDistance(
            new(0, 0, 1), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.Right, 0, Precision.EXACT);

        // Test cases where the closest point is an edge endpoint.
        TestCompareEdgeDistance(
            new(1e-15, -1, 0), new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.Right, -1, Precision.DOUBLE);
        TestCompareEdgeDistance(
            new(-1, -1, 1), new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.Right, 1, Precision.DOUBLE);
        TestCompareEdgeDistance(
            new(1e-18, -1, 0), new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.Right, -1, kLongDoublePrecision);
        TestCompareEdgeDistance(
            new(1e-100, -1, 0), new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.Right, -1, Precision.EXACT);
        TestCompareEdgeDistance(
            new(0, -1, 0), new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.Right, 0, Precision.EXACT);

        // Test cases where x == -a0 or x == -a1.
        TestCompareEdgeDistance(
            new(-1, 0, 0), new(1, 0, 0), new(1, 1, 0),
            S1ChordAngle.Right, 1, Precision.DOUBLE);
        TestCompareEdgeDistance(
            new(-1, 0, 0), new(1, 0, 0), new(1e-18, 1, 0),
            S1ChordAngle.Right, 1, kLongDoublePrecision);
        TestCompareEdgeDistance(
            new(-1, 0, 0), new(1, 0, 0), new(1e-100, 1, 0),
            S1ChordAngle.Right, 1, Precision.EXACT);
        TestCompareEdgeDistance(
            new(0, -1, 0), new(1, 0, 0), new(0, 1, 0),
            S1ChordAngle.Right, 0, Precision.EXACT);
    }

    [Fact]
    internal void Test_CompareEdgeDistance_Consistency()
    {
        // This test chooses random inputs such that the distance between "x" and
        // the line (a0, a1) is very close to the threshold distance "r".  It then
        // checks that the answer given by a method at one level of precision is
        // consistent with the answer given at the next higher level of precision.
        // See also the comments in the CompareDistances consistency test.
        PrecisionStats stats = new();
        for (int iter = 0; iter < consistency_iters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point a0 = ChoosePoint();
            S1Angle len = S1Angle.FromRadians(Math.PI * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
            S2Point a1 = S2.GetPointOnLine(a0, ChoosePoint(), len);
            if (S2Testing.Random.OneIn(2)) a1 = -a1;
            if (a0 == -a1) continue;  // Not allowed by API.
            S2Point n = S2.RobustCrossProd(a0, a1).Normalize();
            double f = Math.Pow(1e-20, S2Testing.Random.RandDouble());
            S2Point a = ((1 - f) * a0 + f * a1).Normalize();
            S1Angle r = S1Angle.FromRadians(S2.M_PI_2 * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
            if (S2Testing.Random.OneIn(2)) r = S1Angle.FromRadians(S2.M_PI_2) - r;
            S2Point x = S2.GetPointOnLine(a, n, r);
            if (S2Testing.Random.OneIn(5))
            {
                // Replace "x" with a random point that is closest to an edge endpoint.
                do
                {
                    x = ChoosePoint();
                } while (S2Pred.CompareEdgeDirections(a0, x, a0, a1) > 0 &&
                         S2Pred.CompareEdgeDirections(x, a1, a0, a1) > 0);
                r = new[] { new S1Angle(x, a0), new S1Angle(x, a1) }.Min();
            }
            Precision prec = TestCompareEdgeDistanceConsistency(x, a0, a1, new S1ChordAngle(r));
            stats.Tally(prec);
        }

        // Checks that the result at one level of precision is consistent with the
        // result at the next higher level of precision.  Returns the minimum
        // precision that yielded a non-zero result.
        static Precision TestCompareEdgeDistanceConsistency(S2Point x, S2Point a0, S2Point a1, S1ChordAngle r)
        {
            int dbl_sign = S2Pred.TriageCompareEdgeDistance_Test(x, a0, a1, r.Length2);
            int ld_sign = S2Pred.TriageCompareEdgeDistance_Test(x.ToLD(), a0.ToLD(), a1.ToLD(), r.Length2.ToLD());
            int exact_sign = S2Pred.ExactCompareEdgeDistance_Test(x, a0, a1, r);
            Assert.Equal(exact_sign, S2Pred.CompareEdgeDistance(x, a0, a1, r));
            if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);
            if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
            return (ld_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.LONG_DOUBLE : Precision.DOUBLE;
        }
    }

    [Fact]
    internal void Test_CompareEdgeDirections_Coverage()
    {
        TestCompareEdgeDirections(new(1, 0, 0), new(1, 1, 0),
                                  new(1, -1, 0), new(1, 0, 0),
                                  1, Precision.DOUBLE);
        TestCompareEdgeDirections(new(1, 0, 1.5e-15), new(1, 1, 0),
                                  new(0, -1, 0), new(0, 0, 1),
                                  1, Precision.DOUBLE);
        TestCompareEdgeDirections(new(1, 0, 1e-18), new(1, 1, 0),
                                  new(0, -1, 0), new(0, 0, 1),
                                  1, kLongDoublePrecision);
        TestCompareEdgeDirections(new(1, 0, 1e-50), new(1, 1, 0),
                                  new(0, -1, 0), new(0, 0, 1),
                                  1, Precision.EXACT);
        TestCompareEdgeDirections(new(1, 0, 0), new(1, 1, 0),
                                  new(0, -1, 0), new(0, 0, 1),
                                  0, Precision.EXACT);
    }

    [Fact]
    internal void Test_CompareEdgeDirections_Consistency()
    {
        // This test chooses random pairs of edges that are nearly perpendicular,
        // then checks that the answer given by a method at one level of precision
        // is consistent with the answer given at the next higher level of
        // precision.  See also the comments in the CompareDistances test.
        PrecisionStats stats = new();
        for (int iter = 0; iter < consistency_iters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point a0 = ChoosePoint();
            S1Angle a_len = S1Angle.FromRadians(Math.PI * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
            S2Point a1 = S2.GetPointOnLine(a0, ChoosePoint(), a_len);
            S2Point a_norm = S2.RobustCrossProd(a0, a1).Normalize();
            S2Point b0 = ChoosePoint();
            S1Angle b_len = S1Angle.FromRadians(Math.PI * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
            S2Point b1 = S2.GetPointOnLine(b0, a_norm, b_len);
            if (a0 == -a1 || b0 == -b1) continue;  // Not allowed by API.
            Precision prec = TestCompareEdgeDirectionsConsistency(a0, a1, b0, b1);
            // Don't skew the statistics by recording degenerate inputs.
            if (a0 == a1 || b0 == b1)
            {
                Assert.Equal(Precision.EXACT, prec);
            }
            else
            {
                stats.Tally(prec);
            }
        }

        // Checks that the result at one level of precision is consistent with the
        // result at the next higher level of precision.  Returns the minimum
        // precision that yielded a non-zero result.
        static Precision TestCompareEdgeDirectionsConsistency(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
        {
            int dbl_sign = S2Pred.TriageCompareEdgeDirections_Test(a0, a1, b0, b1);
            int ld_sign = S2Pred.TriageCompareEdgeDirections_Test(a0.ToLD(), a1.ToLD(), b0.ToLD(), b1.ToLD());
            int exact_sign = S2Pred.ExactCompareEdgeDirections_Test(a0.ToExact(), a1.ToExact(), b0.ToExact(), b1.ToExact());
            Assert.Equal(exact_sign, S2Pred.CompareEdgeDirections(a0, a1, b0, b1));
            if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);
            if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
            return (ld_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.LONG_DOUBLE : Precision.DOUBLE;
        }
    }

    [Fact]
    internal void Test_EdgeCircumcenterSign_Coverage()
    {
        TestEdgeCircumcenterSign(
            new(1, 0, 0), new(1, 1, 0),
            new(0, 0, 1), new(1, 0, 1), new(0, 1, 1),
            1, Precision.DOUBLE);
        TestEdgeCircumcenterSign(
            new(1, 0, 0), new(1, 1, 0),
            new(0, 0, -1), new(1, 0, -1), new(0, 1, -1),
            -1, Precision.DOUBLE);
        TestEdgeCircumcenterSign(
            new(1, -1, 0), new(1, 1, 0),
            new(1, -1e-5, 1), new(1, 1e-5, -1), new(1, 1 - 1e-5, 1e-5),
            -1, Precision.DOUBLE);
        TestEdgeCircumcenterSign(
            new(1, -1, 0), new(1, 1, 0),
            new(1, -1e-5, 1), new(1, 1e-5, -1), new(1, 1 - 1e-9, 1e-5),
            -1, kLongDoublePrecision);
        TestEdgeCircumcenterSign(
            new(1, -1, 0), new(1, 1, 0),
            new(1, -1e-5, 1), new(1, 1e-5, -1), new(1, 1 - 1e-15, 1e-5),
            -1, Precision.EXACT);
        TestEdgeCircumcenterSign(
            new(1, -1, 0), new(1, 1, 0),
            new(1, -1e-5, 1), new(1, 1e-5, -1), new(1, 1, 1e-5),
            1, Precision.SYMBOLIC);

        // This test falls back to the second symbolic perturbation:
        TestEdgeCircumcenterSign(
            new(1, -1, 0), new(1, 1, 0),
            new(0, -1, 0), new(0, 0, -1), new(0, 0, 1),
            -1, Precision.SYMBOLIC);

        // This test falls back to the third symbolic perturbation:
        TestEdgeCircumcenterSign(
            new(0, -1, 1), new(0, 1, 1),
            new(0, 1, 0), new(0, -1, 0), new(1, 0, 0),
            -1, Precision.SYMBOLIC);
    }

    [Fact]
    internal void Test_CompareEdgePairDistance_Coverage()
    {
        // Since CompareEdgePairDistance() is implemented using other predicates, we
        // only test to verify that those predicates are being used correctly.
        S2Point x=new(1, 0, 0), y=new(0, 1, 0), z=new(0, 0, 1);
        S2Point a=new(1, 1e-100, 1e-99), b=new(1, 1e-100, -1e-99);

        // Test cases where the edges have an interior crossing.
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, y, a, b, S1ChordAngle.Zero), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, y, a, b, S1ChordAngle.FromRadians(1)), -1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, y, a, b, S1ChordAngle.FromRadians(-1)), 1);

        // Test cases where the edges share an endpoint.
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, y, x, z, S1ChordAngle.FromRadians(0)), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, y, z, x, S1ChordAngle.FromRadians(0)), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(y, x, x, z, S1ChordAngle.FromRadians(0)), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(y, x, z, x, S1ChordAngle.FromRadians(0)), 0);

        // Test cases where one edge is degenerate.
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, x, x, y, S1ChordAngle.FromRadians(0)), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, y, x, x, S1ChordAngle.FromRadians(0)), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, x, y, z, S1ChordAngle.FromRadians(1)), 1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(y, z, x, x, S1ChordAngle.FromRadians(1)), 1);

        // Test cases where both edges are degenerate.
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, x, x, x, S1ChordAngle.FromRadians(0)), 0);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, x, y, y, S1ChordAngle.FromRadians(1)), 1);

        // Test cases where the minimum distance is non-zero and is achieved at each
        // of the four edge endpoints.
        S1ChordAngle kHi = S1ChordAngle.FromRadians(1e-100 + 1e-115);
        S1ChordAngle kLo = S1ChordAngle.FromRadians(1e-100 - 1e-115);
        Assert.Equal(S2Pred.CompareEdgePairDistance(a, y, x, z, kHi), -1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(a, y, x, z, kLo), 1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(y, a, x, z, kHi), -1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(y, a, x, z, kLo), 1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, z, a, y, kHi), -1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, z, a, y, kLo), 1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, z, y, a, kHi), -1);
        Assert.Equal(S2Pred.CompareEdgePairDistance(x, z, y, a, kLo), 1);
    }

    [Fact]
    internal void Test_EdgeCircumcenterSign_Consistency()
    {
        // This test chooses random a random edge X, then chooses a random point Z
        // on the great circle through X, and finally choose three points A, B, C
        // that are nearly equidistant from X.  It then checks that the answer given
        // by a method at one level of precision is consistent with the answer given
        // at the next higher level of precision.
        PrecisionStats stats = new();
        for (int iter = 0; iter < consistency_iters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point x0 = ChoosePoint();
            S2Point x1 = ChoosePoint();
            if (x0 == -x1) continue;  // Not allowed by API.
            double c0 = (S2Testing.Random.OneIn(2) ? -1 : 1) * Math.Pow(1e-20, S2Testing.Random.RandDouble());
            double c1 = (S2Testing.Random.OneIn(2) ? -1 : 1) * Math.Pow(1e-20, S2Testing.Random.RandDouble());
            S2Point z = (c0 * x0 + c1 * x1).Normalize();
            S1Angle r = S1Angle.FromRadians(Math.PI * Math.Pow(1e-30, S2Testing.Random.RandDouble()));
            S2Point a = S2.GetPointOnLine(z, ChoosePoint(), r);
            S2Point b = S2.GetPointOnLine(z, ChoosePoint(), r);
            S2Point c = S2.GetPointOnLine(z, ChoosePoint(), r);
            Precision prec = TestEdgeCircumcenterSignConsistency(x0, x1, a, b, c);
            // Don't skew the statistics by recording degenerate inputs.
            if (x0 == x1)
            {
                // This precision would be SYMBOLIC if we handled this degeneracy.
                Assert.Equal(Precision.EXACT, prec);
            }
            else if (a == b || b == c || c == a)
            {
                Assert.Equal(Precision.SYMBOLIC, prec);
            }
            else
            {
                stats.Tally(prec);
            }
        }

        // Checks that the result at one level of precision is consistent with the
        // result at the next higher level of precision.  Returns the minimum
        // precision that yielded a non-zero result.
        static Precision TestEdgeCircumcenterSignConsistency(S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c)
        {
            int abc_sign = S2Pred.Sign(a, b, c);
            int dbl_sign = S2Pred.TriageEdgeCircumcenterSign_Test(x0, x1, a, b, c, abc_sign);
            int ld_sign = S2Pred.TriageEdgeCircumcenterSign_Test(
                x0.ToLD(), x1.ToLD(), a.ToLD(), b.ToLD(), c.ToLD(), abc_sign);
            int exact_sign = S2Pred.ExactEdgeCircumcenterSign_Test(
                x0.ToExact(), x1.ToExact(), a.ToExact(), b.ToExact(), c.ToExact(), abc_sign);
            if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);
            if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
            if (exact_sign != 0)
            {
                Assert.Equal(exact_sign, S2Pred.EdgeCircumcenterSign(x0, x1, a, b, c));
                return (ld_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.LONG_DOUBLE : Precision.DOUBLE;
            }
            else
            {
                // Unlike the other methods, SymbolicEdgeCircumcenterSign has the
                // precondition that the exact sign must be zero.
                int symbolic_sign = S2Pred.SymbolicEdgeCircumcenterSign_Test(x0, x1, a, b, c);
                Assert.Equal(symbolic_sign, S2Pred.EdgeCircumcenterSign(x0, x1, a, b, c));
                return Precision.SYMBOLIC;
            }
        }
    }

    [Fact]
    internal void Test_VoronoiSiteExclusion_Coverage()
    {
        // Both sites are closest to edge endpoint X0.
        TestVoronoiSiteExclusion(
            new(1, -1e-5, 0), new(1, -2e-5, 0),
            new(1, 0, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1e-3),
            S2Pred.Excluded.SECOND, Precision.DOUBLE);

        // Both sites are closest to edge endpoint X1.
        TestVoronoiSiteExclusion(
            new(1, 1, 1e-30), new(1, 1, -1e-20),
            new(1, 0, 0), new (1, 1, 0), S1ChordAngle.FromRadians(1e-10),
            S2Pred.Excluded.SECOND, Precision.DOUBLE);

        // Test cases where neither site is excluded.
        TestVoronoiSiteExclusion(
            new(1, -1e-10, 1e-5), new(1, 1e-10, -1e-5),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1e-4),
            S2Pred.Excluded.NEITHER, Precision.DOUBLE);
        TestVoronoiSiteExclusion(
            new(1, -1e-10, 1e-5), new(1, 1e-10, -1e-5),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1e-5),
            S2Pred.Excluded.NEITHER, kLongDoublePrecision);
        TestVoronoiSiteExclusion(
            new(1, -1e-17, 1e-5), new(1, 1e-17, -1e-5),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1e-4),
            S2Pred.Excluded.NEITHER, kLongDoublePrecision);
        TestVoronoiSiteExclusion(
            new(1, -1e-20, 1e-5), new(1, 1e-20, -1e-5),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1e-5),
            S2Pred.Excluded.NEITHER, Precision.EXACT);

        // Test cases where the first site is excluded.  (Tests where the second
        // site is excluded are constructed by TestVoronoiSiteExclusion.)
        TestVoronoiSiteExclusion(
          new(1, -1e-6, 1.0049999999e-5), new(1, 0, -1e-5),
          new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1.005e-5),
          S2Pred.Excluded.FIRST, Precision.DOUBLE);
        TestVoronoiSiteExclusion(
            new(1, -1.00105e-6, 1.0049999999e-5), new(1, 0, -1e-5),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1.005e-5),
            S2Pred.Excluded.FIRST, kLongDoublePrecision);
        TestVoronoiSiteExclusion(
            new(1, -1e-6, 1.005e-5), new(1, 0, -1e-5),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1.005e-5),
            S2Pred.Excluded.FIRST, kLongDoublePrecision);
        TestVoronoiSiteExclusion(
            new(1, -1e-31, 1.005e-30), new(1, 0, -1e-30),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1.005e-30),
            S2Pred.Excluded.FIRST, Precision.EXACT);
        TestVoronoiSiteExclusion(
            new(1, -1e-31, 1.005e-30), new(1, 0, -1e-30),
            new(1, -1, 0), new(1, 1, 0), S1ChordAngle.FromRadians(1.005e-30),
            S2Pred.Excluded.FIRST, Precision.EXACT);

        // Test cases for the (d < 0) portion of the algorithm (see .cc file).  In
        // all of these cases A is closer to X0, B is closer to X1, and AB goes in
        // the opposite direction as edge X when projected onto it (since this is
        // what d < 0 means).

        // 1. Cases that require Pi/2 < d(X0,X1) + r < Pi.  Only one site is kept.
        //
        //    - A and B project to the interior of X.
        TestVoronoiSiteExclusion(
            new(1, -1e-5, 1e-4), new(1, -1.00000001e-5, 0),
            new(-1, -1, 0), new(1, 0, 0), S1ChordAngle.FromRadians(1),
            S2Pred.Excluded.FIRST, Precision.DOUBLE);
        //    - A and B project to opposite sides of X1.
        TestVoronoiSiteExclusion(
            new(1, 1e-10, 0.1), new(1, -1e-10, 1e-8),
            new(-1, -1, 0), new(1, 0, 0), S1ChordAngle.FromRadians(1),
            S2Pred.Excluded.FIRST, Precision.DOUBLE);
        //    - A and B both project to points past X1, and B is closer to the great
        //      circle through edge X.
        TestVoronoiSiteExclusion(
            new(1, 2e-10, 0.1), new(1, 1e-10, 0),
            new(-1, -1, 0), new(1, 0, 0), S1ChordAngle.FromRadians(1),
            S2Pred.Excluded.FIRST, Precision.DOUBLE);
        //    - Like the test above, but A is closer to the great circle through X.
        TestVoronoiSiteExclusion(
            new(1, 1.1, 0), new(1, 1.01, 0.01),
            new(-1, -1, 0), new(1, 0, 0), S1ChordAngle.FromRadians(1),
            S2Pred.Excluded.FIRST, Precision.DOUBLE);

        // 2. Cases that require d(X0,X1) + r > Pi and where only one site is kept.
        //
        //    - B is closer to edge X (in fact it's right on the edge), but when A
        //      and B are projected onto the great circle through X they are more
        //      than 90 degrees apart.  This case requires that the sin(d) < 0 case
        //      in the algorithm is handled *before* the cos(d) < 0 case.
        TestVoronoiSiteExclusion(
            new(1, 1.1, 0), new(1, -1, 0),
            new(-1, 0, 0), new(1, -1e-10, 0), S1ChordAngle.FromDegrees(70),
            S2Pred.Excluded.FIRST, Precision.DOUBLE);

        // 3. Cases that require d(X0,X1) + r > Pi and where both sites are kept.
        //
        //    - A projects to a point past X0, B projects to a point past X1,
        //      neither site should be excluded, and A is closer to the great circle
        //      through edge X.
        TestVoronoiSiteExclusion(
            new(-1, 0.1, 0.001), new(1, 1.1, 0),
            new(-1, -1, 0), new(1, 0, 0), S1ChordAngle.FromRadians(1),
            S2Pred.Excluded.NEITHER, Precision.DOUBLE);
        //    - Like the above, but B is closer to the great circle through edge X.
        TestVoronoiSiteExclusion(
            new(-1, 0.1, 0), new(1, 1.1, 0.001),
            new(-1, -1, 0), new(1, 0, 0), S1ChordAngle.FromRadians(1),
            S2Pred.Excluded.NEITHER, Precision.DOUBLE);

        // These two sites are exactly 60 degrees away from the point (1, 1, 0),
        // which is the midpoint of edge X.  This case requires symbolic
        // perturbations to resolve correctly.  Site A is closer to every point in
        // its coverage interval except for (1, 1, 0), but site B is considered
        // closer to that point symbolically.
        TestVoronoiSiteExclusion(
            new(0, 1, 1), new(1, 0, 1),
            new(0, 1, 1), new(1, 0, -1), S1ChordAngle.FromLength2(1),
            S2Pred.Excluded.NEITHER, Precision.EXACT);

        // This test is similar except that site A is considered closer to the
        // equidistant point (-1, 1, 0), and therefore site B is excluded.
        TestVoronoiSiteExclusion(
            new(0, 1, 1), new(-1, 0, 1),
            new(0, 1, 1), new(-1, 0, -1), S1ChordAngle.FromLength2(1),
            S2Pred.Excluded.SECOND, Precision.EXACT);
    }

    [Fact]
    internal void Test_VoronoiSiteExclusion_Consistency()
    {
        // This test chooses random a random edge X, a random point P on that edge,
        // and a random threshold distance "r".  It then choose two sites A and B
        // whose distance to P is almost exactly "r".  This ensures that the
        // coverage intervals for A and B will (almost) share a common endpoint.  It
        // then checks that the answer given by a method at one level of precision
        // is consistent with the answer given at higher levels of precision.
        PrecisionStats stats = new();
        for (int iter = 0; iter < consistency_iters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point x0 = ChoosePoint();
            S2Point x1 = ChoosePoint();
            if (x0 == -x1) continue;  // Not allowed by API.
            double f = Math.Pow(1e-20, S2Testing.Random.RandDouble());
            S2Point p = ((1 - f) * x0 + f * x1).Normalize();
            S1Angle r1 = S1Angle.FromRadians(S2.M_PI_2 * Math.Pow(1e-20, S2Testing.Random.RandDouble()));
            S2Point a = S2.GetPointOnLine(p, ChoosePoint(), r1);
            S2Point b = S2.GetPointOnLine(p, ChoosePoint(), r1);
            // Check that the other API requirements are met.
            S1ChordAngle r = new(r1);
            if (S2Pred.CompareEdgeDistance(a, x0, x1, r) > 0) continue;
            if (S2Pred.CompareEdgeDistance(b, x0, x1, r) > 0) continue;
            if (S2Pred.CompareDistances(x0, a, b) > 0) { (b, a) = (a, b); }
            if (a == b) continue;

            Precision prec = TestVoronoiSiteExclusionConsistency(a, b, x0, x1, r);
            // Don't skew the statistics by recording degenerate inputs.
            if (x0 == x1)
            {
                Assert.Equal(Precision.DOUBLE, prec);
            }
            else
            {
                stats.Tally(prec);
            }
        }

        // Checks that the result at one level of precision is consistent with the
        // result at the next higher level of precision.  Returns the minimum
        // precision that yielded a non-zero result.
        static Precision TestVoronoiSiteExclusionConsistency(S2Point a, S2Point b, S2Point x0, S2Point x1, S1ChordAngle r)
        {

            // The internal methods require this (see TestVoronoiSiteExclusion).
            if (S2Pred.CompareDistances(x1, a, b) < 0) return Precision.DOUBLE;

            S2Pred.Excluded dbl_result = S2Pred.TriageVoronoiSiteExclusion_Test(a, b, x0, x1, r.Length2);
            S2Pred.Excluded ld_result = S2Pred.TriageVoronoiSiteExclusion_Test(
                a.ToLD(), b.ToLD(), x0.ToLD(), x1.ToLD(), r.Length2.ToLD());
            S2Pred.Excluded exact_result = S2Pred.ExactVoronoiSiteExclusion_Test(
                a.ToExact(), b.ToExact(), x0.ToExact(), x1.ToExact(), r.Length2);
            Assert.Equal(exact_result, S2Pred.GetVoronoiSiteExclusion(a, b, x0, x1, r));

            Assert.NotEqual(S2Pred.Excluded.UNCERTAIN, exact_result);
            if (ld_result == S2Pred.Excluded.UNCERTAIN)
            {
                Assert.Equal(S2Pred.Excluded.UNCERTAIN, dbl_result);
                return Precision.EXACT;
            }
            Assert.Equal(exact_result, ld_result);
            if (dbl_result == S2Pred.Excluded.UNCERTAIN)
            {
                return Precision.LONG_DOUBLE;
            }
            Assert.Equal(exact_result, dbl_result);
            return Precision.DOUBLE;
        }
    }

    // Verifies that VoronoiSiteExclusion(a, b, x0, x1, r) == expected_result, and
    // furthermore checks that the minimum required precision is "expected_prec".
    private static void TestVoronoiSiteExclusion(S2Point a, S2Point b, S2Point x0, S2Point x1, S1ChordAngle r, S2Pred.Excluded expected_result, Precision expected_prec)
    {
        const S2Pred.Excluded UNCERTAIN = S2Pred.Excluded.UNCERTAIN;

        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!a.IsUnitLength()) a = a.Normalize();
        if (!b.IsUnitLength()) b = b.Normalize();
        if (!x0.IsUnitLength()) x0 = x0.Normalize();
        if (!x1.IsUnitLength()) x1 = x1.Normalize();

        // The internal methods (Triage, Exact, etc) require that site A is closer
        // to X0 and site B is closer to X1.  GetVoronoiSiteExclusion has special
        // code to handle the case where this is not true.  We need to duplicate
        // that code here.  Essentially, since the API requires site A to be closer
        // than site B to X0, then if site A is also closer to X1 then site B must
        // be excluded.
        if (S2Pred.CompareDistances(x1, a, b) < 0)
        {
            Assert.Equal(S2Pred.Excluded.SECOND, expected_result);
            // We don't know what precision was used by CompareDistances(), but we
            // arbitrarily require the test to specify it as DOUBLE.
            Assert.Equal(Precision.DOUBLE, expected_prec);
        }
        else
        {
            S2Pred.Excluded dbl_result = S2Pred.TriageVoronoiSiteExclusion_Test(a, b, x0, x1, r.Length2);
            S2Pred.Excluded ld_result = S2Pred.TriageVoronoiSiteExclusion_Test(
                a.ToLD(), b.ToLD(), x0.ToLD(), x1.ToLD(), r.Length2.ToLD());
            S2Pred.Excluded exact_result = S2Pred.ExactVoronoiSiteExclusion_Test(
                a.ToExact(), b.ToExact(), x0.ToExact(), x1.ToExact(), r.Length2);

            // Check that the results are correct (if not UNCERTAIN), and also that if
            // dbl_result is not UNCERTAIN then so is ld_result, etc.
            Assert.Equal(expected_result, exact_result);
            if (ld_result != UNCERTAIN) Assert.Equal(exact_result, ld_result);
            if (dbl_result != UNCERTAIN) Assert.Equal(ld_result, dbl_result);

            Precision actual_prec = (dbl_result != UNCERTAIN ? Precision.DOUBLE :
                                     ld_result != UNCERTAIN ? Precision.LONG_DOUBLE : Precision.EXACT);
            Assert.Equal(expected_prec, actual_prec);
        }
        // Make sure that the top-level function returns the expected result.
        Assert.Equal(expected_result, S2Pred.GetVoronoiSiteExclusion(a, b, x0, x1, r));

        // If site B is closer to X1, then the same site should be excluded (if any)
        // when we swap the sites and the edge direction.
        S2Pred.Excluded swapped_result =
            expected_result == S2Pred.Excluded.FIRST ? S2Pred.Excluded.SECOND :
            expected_result == S2Pred.Excluded.SECOND ? S2Pred.Excluded.FIRST : expected_result;
        if (S2Pred.CompareDistances(x1, b, a) < 0)
        {
            Assert.Equal(swapped_result, S2Pred.GetVoronoiSiteExclusion(b, a, x1, x0, r));
        }
    }

    // Verifies that EdgeCircumcenterSign(x0, x1, a, b, c) == expected_sign, and
    // furthermore checks that the minimum required precision is "expected_prec".
    private static void TestEdgeCircumcenterSign(S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c, int expected_sign, Precision expected_prec)
    {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!x0.IsUnitLength()) x0 = x0.Normalize();
        if (!x1.IsUnitLength()) x1 = x1.Normalize();
        if (!a.IsUnitLength()) a = a.Normalize();
        if (!b.IsUnitLength()) b = b.Normalize();
        if (!c.IsUnitLength()) c = c.Normalize();

        int abc_sign = S2Pred.Sign(a, b, c);
        int dbl_sign = S2Pred.TriageEdgeCircumcenterSign_Test(x0, x1, a, b, c, abc_sign);
        int ld_sign = S2Pred.TriageEdgeCircumcenterSign_Test(
            x0.ToLD(), x1.ToLD(), a.ToLD(), b.ToLD(), c.ToLD(), abc_sign);
        int exact_sign = S2Pred.ExactEdgeCircumcenterSign_Test(
            x0.ToExact(), x1.ToExact(), a.ToExact(), b.ToExact(), c.ToExact(), abc_sign);
        int actual_sign = (exact_sign != 0 ? exact_sign :
                           S2Pred.SymbolicEdgeCircumcenterSign_Test(x0, x1, a, b, c));

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        Assert.Equal(expected_sign, actual_sign);
        if (exact_sign != 0) Assert.Equal(exact_sign, actual_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);

        Precision actual_prec = ((dbl_sign != 0) ? Precision.DOUBLE :
                                 (ld_sign != 0) ? Precision.LONG_DOUBLE :
                                 (exact_sign != 0) ? Precision.EXACT : Precision.SYMBOLIC);
        Assert.Equal(expected_prec, actual_prec);

        // Make sure that the top-level function returns the expected result.
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, a, b, c));

        // Check various identities involving swapping or negating arguments.
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, a, c, b));
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, b, a, c));
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, b, c, a));
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, c, a, b));
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, c, b, a));
        Assert.Equal(-expected_sign, S2Pred.EdgeCircumcenterSign(x1, x0, a, b, c));
        Assert.Equal(expected_sign, S2Pred.EdgeCircumcenterSign(-x0, -x1, a, b, c));
        if (actual_sign == exact_sign)
        {
            // Negating the input points may not preserve the result when symbolic
            // perturbations are used, since -X is not an exact multiple of X.
            Assert.Equal(-expected_sign, S2Pred.EdgeCircumcenterSign(x0, x1, -a, -b, -c));
        }
    }

    // Verifies that CompareEdgeDirections(a0, a1, b0, b1) == expected_sign, and
    // furthermore checks that the minimum required precision is "expected_prec".
    private static void TestCompareEdgeDirections(S2Point a0, S2Point a1, S2Point b0, S2Point b1, int expected_sign, Precision expected_prec)
    {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!a0.IsUnitLength()) a0 = a0.Normalize();
        if (!a1.IsUnitLength()) a1 = a1.Normalize();
        if (!b0.IsUnitLength()) b0 = b0.Normalize();
        if (!b1.IsUnitLength()) b1 = b1.Normalize();

        int dbl_sign = S2Pred.TriageCompareEdgeDirections_Test(a0, a1, b0, b1);
        int ld_sign = S2Pred.TriageCompareEdgeDirections_Test(a0.ToLD(), a1.ToLD(), b0.ToLD(), b1.ToLD());
        int exact_sign = S2Pred.ExactCompareEdgeDirections_Test(a0.ToExact(), a1.ToExact(), b0.ToExact(), b1.ToExact());

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        Assert.Equal(expected_sign, exact_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);

        Precision actual_prec = ((dbl_sign != 0) ? Precision.DOUBLE : (ld_sign != 0) ? Precision.LONG_DOUBLE : Precision.EXACT);
        Assert.Equal(expected_prec, actual_prec);

        // Make sure that the top-level function returns the expected result.
        Assert.Equal(expected_sign, S2Pred.CompareEdgeDirections(a0, a1, b0, b1));

        // Check various identities involving swapping or negating arguments.
        Assert.Equal(expected_sign, S2Pred.CompareEdgeDirections(b0, b1, a0, a1));
        Assert.Equal(expected_sign, S2Pred.CompareEdgeDirections(-a0, -a1, b0, b1));
        Assert.Equal(expected_sign, S2Pred.CompareEdgeDirections(a0, a1, -b0, -b1));
        Assert.Equal(-expected_sign, S2Pred.CompareEdgeDirections(a1, a0, b0, b1));
        Assert.Equal(-expected_sign, S2Pred.CompareEdgeDirections(a0, a1, b1, b0));
        Assert.Equal(-expected_sign, S2Pred.CompareEdgeDirections(-a0, a1, b0, b1));
        Assert.Equal(-expected_sign, S2Pred.CompareEdgeDirections(a0, -a1, b0, b1));
        Assert.Equal(-expected_sign, S2Pred.CompareEdgeDirections(a0, a1, -b0, b1));
        Assert.Equal(-expected_sign, S2Pred.CompareEdgeDirections(a0, a1, b0, -b1));
    }

    // Verifies that CompareEdgeDistance(x, a0, a1, r) == expected_sign, and
    // furthermore checks that the minimum required precision is "expected_prec".
    private static void TestCompareEdgeDistance(S2Point x, S2Point a0, S2Point a1, S1ChordAngle r, int expected_sign, Precision expected_prec)
    {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!x.IsUnitLength()) x = x.Normalize();
        if (!a0.IsUnitLength()) a0 = a0.Normalize();
        if (!a1.IsUnitLength()) a1 = a1.Normalize();

        int dbl_sign = S2Pred.TriageCompareEdgeDistance_Test(x, a0, a1, r.Length2);
        int ld_sign = S2Pred.TriageCompareEdgeDistance_Test(x.ToLD(), a0.ToLD(), a1.ToLD(), r.Length2.ToLD());
        int exact_sign = S2Pred.ExactCompareEdgeDistance_Test(x, a0, a1, r);

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        Assert.Equal(expected_sign, exact_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);

        Precision actual_prec = ((dbl_sign != 0) ? Precision.DOUBLE : (ld_sign != 0) ? Precision.LONG_DOUBLE : Precision.EXACT);
        Assert.Equal(expected_prec, actual_prec);

        // Make sure that the top-level function returns the expected result.
        Assert.Equal(expected_sign, S2Pred.CompareEdgeDistance(x, a0, a1, r));
    }

    // Checks that the result at one level of precision is consistent with the
    // result at the next higher level of precision.  Returns the minimum
    // precision that yielded a non-zero result.
    private static Precision TestCompareDistanceConsistency(S2Point x, S2Point y, S1ChordAngle r, Triage triage)
    {
        int dbl_sign = triage(x, y, r);
        int ld_sign = triage(x.ToLD(), y.ToLD(), r);
        int exact_sign = S2Pred.ExactCompareDistance_Test(x.ToExact(), y.ToExact(), r.Length2);
        Assert.Equal(exact_sign, S2Pred.CompareDistance(x, y, r));
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        return (ld_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.LONG_DOUBLE : Precision.DOUBLE;
    }

    // Verifies that CompareDistance(x, y, r) == expected_sign, and furthermore
    // checks that the minimum required precision is "expected_prec" when the
    // distance calculation method defined by CompareDistanceWrapper is used.
    private static void TestCompareDistance(S2Point x, S2Point y, S1ChordAngle r, int expected_sign, Precision expected_prec, Triage triage)
    {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!x.IsUnitLength()) x = x.Normalize();
        if (!y.IsUnitLength()) y = y.Normalize();

        int dbl_sign = triage(x, y, r);
        int ld_sign = triage(x.ToLD(), y.ToLD(), r);
        int exact_sign = S2Pred.ExactCompareDistance_Test(x.ToExact(), y.ToExact(), r.Length2);

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        Assert.Equal(expected_sign, exact_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);

        Precision actual_prec = ((dbl_sign != 0) ? Precision.DOUBLE : (ld_sign != 0) ? Precision.LONG_DOUBLE : Precision.EXACT);
        Assert.Equal(expected_prec, actual_prec);

        // Make sure that the top-level function returns the expected result.
        Assert.Equal(expected_sign, S2Pred.CompareDistance(x, y, r));

        // Mathematically, if d(X, Y) < r then d(-X, Y) > (Pi - r).  Unfortunately
        // there can be rounding errors when computing the supplementary distance,
        // so to ensure the two distances are exactly supplementary we need to do
        // the following.
        S1ChordAngle r_supp = S1ChordAngle.Straight - r;
        r = S1ChordAngle.Straight - r_supp;
        Assert.Equal(-S2Pred.CompareDistance(x, y, r), S2Pred.CompareDistance(-x, y, r_supp));
    }

    // Checks that the result at one level of precision is consistent with the
    // result at the next higher level of precision.  Returns the minimum
    // precision that yielded a non-zero result.
    private static Precision TestCompareDistancesConsistency(S2Point x, S2Point a, S2Point b, TriagePoints triage)
    {
        int dbl_sign = triage(x, a, b);
        int ld_sign = triage(x.ToLD(), a.ToLD(), b.ToLD());
        int exact_sign = S2Pred.ExactCompareDistances_Test(
            x.ToExact(), a.ToExact(), b.ToExact());
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        if (exact_sign != 0)
        {
            Assert.Equal(exact_sign, S2Pred.CompareDistances(x, a, b));
            return (ld_sign == 0) ? Precision.EXACT : (dbl_sign == 0)
                ? Precision.LONG_DOUBLE : Precision.DOUBLE;
        }
        else
        {
            // Unlike the other methods, SymbolicCompareDistances has the
            // precondition that the exact sign must be zero.
            int symbolic_sign = S2Pred.SymbolicCompareDistances_Test(x, a, b);
            Assert.Equal(symbolic_sign, S2Pred.CompareDistances(x, a, b));
            return Precision.SYMBOLIC;
        }
    }

    // Given 3 points A, B, C that are exactly coplanar with the origin and where
    // A < B < C in lexicographic order, verify that ABC is counterclockwise (if
    // expected == 1) or clockwise (if expected == -1) using ExpensiveSign().
    //
    // This method is intended specifically for checking the cases where
    // symbolic perturbations are needed to break ties.
    private static void CheckSymbolicSign(int expected, S2Point a, S2Point b, S2Point c)
    {
        Assert.True(a < b);
        Assert.True(b < c);
        Assert.Equal(0, a.DotProd(b.CrossProd(c)));

        // Use ASSERT rather than EXPECT to suppress spurious error messages.
        Assert.Equal(expected, S2Pred.ExpensiveSign(a, b, c));
        Assert.Equal(expected, S2Pred.ExpensiveSign(b, c, a));
        Assert.Equal(expected, S2Pred.ExpensiveSign(c, a, b));
        Assert.Equal(-expected, S2Pred.ExpensiveSign(c, b, a));
        Assert.Equal(-expected, S2Pred.ExpensiveSign(b, a, c));
        Assert.Equal(-expected, S2Pred.ExpensiveSign(a, c, b));
    }

    // Returns 2 ** (-digits).  This could be implemented using "ldexp" except
    // that ldexp is notexpr in C++11.
    private static double EpsilonForDigits(int digits)
    {
        return (digits < 64 ? (1.0 / (1UL << digits))
            : (EpsilonForDigits(digits - 63) / (1UL << 63)));
    }

    // Verifies that CompareDistances(x, a, b) == expected_sign, and furthermore
    // checks that the minimum required precision is "expected_prec" when the
    // distance calculation method defined by CompareDistancesWrapper is used.
    private static void TestCompareDistances(S2Point x, S2Point a, S2Point b, int expected_sign, Precision expected_prec, TriagePoints triage)
    {
        // Don't normalize the arguments unless necessary (to allow testing points
        // that differ only in magnitude).
        if (!x.IsUnitLength()) x = x.Normalize();
        if (!a.IsUnitLength()) a = a.Normalize();
        if (!b.IsUnitLength()) b = b.Normalize();

        int dbl_sign = triage(x, a, b);
        int ld_sign = triage(x.ToLD(), a.ToLD(), b.ToLD());
        int exact_sign = S2Pred.ExactCompareDistances_Test(x.ToExact(), a.ToExact(), b.ToExact());
        int actual_sign = (exact_sign != 0 ? exact_sign : S2Pred.SymbolicCompareDistances_Test(x, a, b));

        // Check that the signs are correct (if non-zero), and also that if dbl_sign
        // is non-zero then so is ld_sign, etc.
        Assert.Equal(expected_sign, actual_sign);
        if (exact_sign != 0) Assert.Equal(exact_sign, actual_sign);
        if (ld_sign != 0) Assert.Equal(exact_sign, ld_sign);
        if (dbl_sign != 0) Assert.Equal(ld_sign, dbl_sign);

        Precision actual_prec = ((dbl_sign != 0) ? Precision.DOUBLE :
                                 (ld_sign != 0) ? Precision.LONG_DOUBLE :
                                 (exact_sign != 0) ? Precision.EXACT : Precision.SYMBOLIC);
        Assert.Equal(expected_prec, actual_prec);

        // Make sure that the top-level function returns the expected result.
        Assert.Equal(expected_sign, S2Pred.CompareDistances(x, a, b));

        // Check that reversing the arguments negates the result.
        Assert.Equal(-expected_sign, S2Pred.CompareDistances(x, b, a));
    }

    #region Triage Delegates

    // Delegates for testing the various distance calculation methods.

    internal delegate int Triage(S2Point x, S2Point y, S1ChordAngle r);

    internal static readonly Triage Sin2Distance = (S2Point x, S2Point y, S1ChordAngle r)
        => S2Pred.TriageCompareSin2Distance_Test(x, y, r.Length2);

    internal static readonly Triage CosDistance = (S2Point x, S2Point y, S1ChordAngle r)
        => S2Pred.TriageCompareCosDistance_Test(x, y, r.Length2);

    // The following helper classes allow us to test the various distance
    // calculation methods using a common test framework.
    internal delegate int TriagePoints(S2Point x, S2Point y, S2Point b);
    internal static readonly TriagePoints Sin2Distances = (S2Point x, S2Point a, S2Point b)
        => S2Pred.TriageCompareSin2Distances_Test(x, a, b);

    internal static readonly TriagePoints CosDistances = (S2Point x, S2Point a, S2Point b)
        => S2Pred.TriageCompareCosDistances_Test(x, a, b);

    // Compares distances greater than 90 degrees using sin^2(distance).
    internal static readonly TriagePoints MinusSin2Distances = (S2Point x, S2Point a, S2Point b)
        => -S2Pred.TriageCompareSin2Distances_Test(-x, a, b);

    #endregion

    // Construct approximately "n" points near the great circle through A and B,
    // then sort them and test whether they are sorted.
    private static void TestGreatCircle(S2Point a, S2Point b, int n)
    {
        a = a.Normalize();
        b = b.Normalize();
        List<S2Point> points = new(){ a, b };
        while (points.Count < n)
        {
            AddDegeneracy(points);
        }
        // Remove any (0, 0, 0) points that were accidentically created, then sort
        // the points and remove duplicates.
        points.Remove(new S2Point(0, 0, 0));
        points = new SortedSet<S2Point>(points).ToList();
        Assert.True(points.Count >= n / 2);

        SortAndTest(points, a);
        SortAndTest(points, b);
        foreach (var origin in points)
        {
            SortAndTest(points, origin);
        }

        // Add zero or more (but usually one) point that is likely to trigger
        // Sign() degeneracies among the given points.
        static void AddDegeneracy(List<S2Point> points)
        {
            var aIndex = S2Testing.Random.Uniform(points.Count);
            var bIndex = S2Testing.Random.Uniform(points.Count);
            S2Point a = points[aIndex];
            S2Point b = points[bIndex];
            int coord = S2Testing.Random.Uniform(3);
            switch (S2Testing.Random.Uniform(8))
            {
                case 0:
                    // Add a random point (not uniformly distributed) along the great
                    // circle AB.
                    points.Add((S2Testing.Random.UniformDouble(-1, 1) * a +
                                  S2Testing.Random.UniformDouble(-1, 1) * b).Normalize());
                    break;
                case 1:
                    // Perturb one coordinate by the minimum amount possible.
                    a = a.SetAxis(coord, MathUtils.NextAfter(a[coord], S2Testing.Random.OneIn(2) ? 2 : -2));
                    points.Add(a.Normalize());
                    break;
                case 2:
                    // Perturb one coordinate by up to 1e-15.
                    a = a.SetAxis(coord, 1e-15 * S2Testing.Random.UniformDouble(-1, 1));
                    points.Add(a.Normalize());
                    break;
                case 3:
                    // Scale a point just enough so that it is different while still being
                    // considered normalized.
                    a *= S2Testing.Random.OneIn(2) ? (1 + 2e-16) : (1 - 1e-16);
                    if (a.IsUnitLength()) points.Add(a);
                    break;
                case 4:
                    {
                        // Add the intersection point of AB with X=0, Y=0, or Z=0.
                        S2Point dir = new(0, 0, 0);
                        dir = dir.SetAxis(coord, S2Testing.Random.OneIn(2) ? 1 : -1);
                        var norm = S2.RobustCrossProd(a, b).Normalize();
                        if (norm.Norm2() > 0)
                        {
                            points.Add(S2.RobustCrossProd(dir, norm).Normalize());
                        }
                        break;
                    }
                case 5:
                    // Add two closely spaced points along the tangent at A to the great
                    // circle through AB.
                    AddTangentPoints(a, b, points);
                    break;
                case 6:
                    // Add two closely spaced points along the tangent at A to the great
                    // circle through A and the X-axis.
                    AddTangentPoints(a, new S2Point(1, 0, 0), points);
                    break;
                case 7:
                    // Add the negative of a point.
                    points.Add(-a);
                    break;
            }
        }
    }

    // Sort the points around the given origin, and then do some consistency
    // checks to verify that they are actually sorted.
    private static void SortAndTest(List<S2Point> points, S2Point origin)
    {
        var sorted = new List<S2Point>();
        SortCCW(points, origin, sorted);
        TestCCW(sorted, origin);

        // Given a set of points with no duplicates, first remove "origin" from
        // "points" (if it exists) and then sort the remaining points in CCW order
        // around "origin" putting the result in "sorted".
        static void SortCCW(List<S2Point> points, S2Point origin, List<S2Point> sorted)
        {
            // Make a copy of the points with "origin" removed.
            sorted.Clear();
            sorted.AddRange(points);
            sorted.Remove(origin);

            // Sort the points CCW around the origin starting at (*sorted)[0].
            var less = new SignTest.LessCCW(origin, sorted[0]);
            sorted.Sort(less);
        }

        // Test exhaustively whether the points in "sorted" are sorted circularly
        // CCW around "origin".
        static void TestCCW(List<S2Point> sorted, S2Point origin)
        {
            int n = sorted.Count;
            int total_num_ccw = 0;
            int last_num_ccw = CountCCW(sorted, origin, n - 1);
            for (int start = 0; start < n; ++start)
            {
                int num_ccw = CountCCW(sorted, origin, start);
                // Each iteration we increase the start index by 1, therefore the number
                // of CCW triangles should decrease by at most 1.
                Assert.True(num_ccw >= last_num_ccw - 1);
                total_num_ccw += num_ccw;
                last_num_ccw = num_ccw;
            }
            // We have tested all triangles of the form OAB.  Exactly half of these
            // should be CCW.
            Assert.Equal(n * (n - 1) / 2, total_num_ccw);

            // Given a set of points sorted circularly CCW around "origin", and the
            // index "start" of a point A, count the number of CCW triangles OAB over
            // all sorted points B not equal to A.  Also check that the results of the
            // CCW tests are consistent with the hypothesis that the points are sorted.
            static int CountCCW(List<S2Point> sorted, S2Point origin, int start)
            {
                int num_ccw = 0;
                int last_sign = 1;
                int n = sorted.Count;
                for (int j = 1; j < n; ++j)
                {
                    int sign = S2Pred.Sign(origin, sorted[start], sorted[(start + j) % n]);
                    Assert.NotEqual(0, sign);
                    if (sign > 0) ++num_ccw;

                    // Since the points are sorted around the origin, we expect to see a
                    // (possibly empty) sequence of CCW triangles followed by a (possibly
                    // empty) sequence of CW triangles.
                    Assert.False(sign > 0 && last_sign < 0);
                    last_sign = sign;
                }
                return num_ccw;
            }
        }
    }

    // Add two points A1 and A2 that are slightly offset from A along the
    // tangent toward B, and such that A, A1, and A2 are exactly collinear
    // (i.e. even with infinite-precision arithmetic).
    private static void AddTangentPoints(S2Point a, S2Point b, List<S2Point> points)
    {
        S2Point dir = S2.RobustCrossProd(a, b).CrossProd(a).Normalize();
        if (dir == new S2Point(0, 0, 0)) return;
        for (; ; )
        {
            S2Point delta = 1e-15 * S2Testing.Random.RandDouble() * dir;
            if ((a + delta) != a && (a + delta) - a == a - (a - delta) &&
                (a + delta).IsUnitLength() && (a - delta).IsUnitLength())
            {
                points.Add(a + delta);
                points.Add(a - delta);
                return;
            }
        }
    }

    // Chooses a random S2Point that is often near the intersection of one of the
    // coodinates planes or coordinate axes with the unit sphere.  (It is possible
    // to represent very small perturbations near such points.)
    private static S2Point ChoosePoint()
    {
        S2Point p = S2Testing.RandomPoint();
        var x = new double[] { p.X, p.Y, p.Z };
        for (int i = 0; i < 3; ++i)
        {
            if (S2Testing.Random.OneIn(3))
            {
                x[i] *= Math.Pow(1e-50, S2Testing.Random.RandDouble());
            }
        }
        return new S2Point(x).Normalize();
    }

    [Fact]
    internal void Test_Sign_StableSignUnderflow()
    {
        // Verify that StableSign returns zero (indicating that the result is
        // uncertain) when its error calculation underflows.
        S2Point a = new(1, 1.9535722048627587e-90, 7.4882501322554515e-80);
        S2Point b = new(1, 9.6702373087191359e-127, 3.706704857169321e-116);
        S2Point c = new(1, 3.8163353663361477e-142, 1.4628419538608985e-131);

        Assert.Equal(S2Pred.StableSign(a, b, c), 0);
        Assert.Equal(S2Pred.ExactSign(a, b, c, true), 1);
        Assert.Equal(S2Pred.Sign(a, b, c), 1);
    }

    // This test repeatedly constructs some number of points that are on or nearly
    // on a given great circle.  Then it chooses one of these points as the
    // "origin" and sorts the other points in CCW order around it.  Of course,
    // since the origin is on the same great circle as the points being sorted,
    // nearly all of these tests are degenerate.  It then does various consistency
    // checks to verify that the points are indeed sorted in CCW order.
    //
    // It is easier to think about what this test is doing if you imagine that the
    // points are in general position rather than on a great circle.
    internal class SignTest
    {
        // The following method is used to sort a collection of points in CCW order
        // around a given origin.  It returns true if A comes before B in the CCW
        // ordering (starting at an arbitrary fixed direction).
        internal class LessCCW : IComparer<S2Point>
        {
            private readonly S2Point origin_;
            private readonly S2Point start_;

            internal LessCCW(S2Point origin, S2Point start)
            { origin_ = origin; start_ = start; }

            public int Compare(S2Point a, S2Point b)
            {
                return S2Pred.OrderedCCW(start_, b, a, origin_) ? -1 : 1;
            }
        }
    }

    internal static class StableSignTest
    {
        // Estimate the probability that S2.StableSign() will not be able to compute
        // the determinant sign of a triangle A, B, C consisting of three points
        // that are as collinear as possible and spaced the given distance apart.
        internal static double GetFailureRate(double km)
        {
            int kIters = 1000;
            int failure_count = 0;
            double m = Math.Tan(S2Testing.KmToAngle(km).Radians);
            for (int iter = 0; iter < kIters; ++iter)
            {
                S2Testing.GetRandomFrame(out var a, out var x, out _);
                S2Point b = (a - m * x).Normalize();
                S2Point c = (a + m * x).Normalize();
                int sign = S2Pred.StableSign_Test(a, b, c);
                if (sign != 0)
                {
                    Assert.Equal(S2Pred.ExactSign_Test(a, b, c, true), sign);
                }
                else
                {
                    ++failure_count;
                }
            }
            return ((double)failure_count) / kIters;
        }
    }

    internal enum Precision { DOUBLE, LONG_DOUBLE, EXACT, SYMBOLIC, NUM_PRECISIONS };

    // A helper class that keeps track of how often each precision was used and
    // generates a string for logging purposes.
    internal class PrecisionStats
    {
        private readonly int[] counts_ = new int[4].Fill(0); // NUM_PRECISIONS

        internal PrecisionStats() { }

        internal void Tally(Precision precision) { ++counts_[(int)precision]; }

        public override string ToString()
        {
            StringBuilder sb = new();
            int total = 0;
            for (int i = 0; i < 4; ++i) // NUM_PRECISIONS
            {
                sb.AppendFormat($"{kPrecisionNames[i]}={counts_[i]:6d}, ");
                total += counts_[i];
            }
            sb.AppendFormat($"total={total:6d}");
            return sb.ToString();
        }
    }
}
