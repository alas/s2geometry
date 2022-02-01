namespace S2Geometry;

using static S2.Internal;

// CrossingSign, VertexCrossing, and EdgeOrVertexCrossing are tested in
// S2EdgeCrosserTests.
public class S2EdgeCrossingsTests
{
    // The approximate maximum error in GetDistance() for small distances.
    private static readonly S1Angle kGetDistanceAbsError = S1Angle.FromRadians(3 * S2.DoubleEpsilon);
    private readonly ITestOutputHelper _logger;
    public enum Precision { DOUBLE, LONG_DOUBLE, EXACT, SYMBOLIC, NUM_PRECISIONS }
    private static readonly string[] kPrecisionNames = new[] { "double", "long double", "exact", "symbolic" };
    public S2EdgeCrossingsTests(ITestOutputHelper logger) { _logger = logger; }

    // A helper class that keeps track of how often each precision was used and
    // generates a string for logging purposes.
    public class PrecisionStats
    {
        public PrecisionStats()
        {
            for (int i = 0; i < counts_.Length; i++)
                counts_[i] = 0;
        }
        public void Tally(Precision precision) { ++counts_[(int)precision]; }
        public override string ToString()
        {
            System.Text.StringBuilder result = new();
            int total = 0;
            for (int i = 0; i < (int)Precision.NUM_PRECISIONS; ++i)
            {
                result.AppendFormat("{0}={1:D6}, ", kPrecisionNames[i], counts_[i]);
                total += counts_[i];
            }
            result.AppendFormat("total={0:D6}", total);
            return result.ToString();
        }
        private readonly int[] counts_ = new int[(int)Precision.NUM_PRECISIONS];
    }

    // Checks that the RobustCrossProd result at one level of precision is
    // consistent with the result at the next higher level of precision.  Also
    // verifies that the results satisfy the identities documented in the header
    // file, including consistency with S2::Sign.  Returns the minimum precision
    // that yielded a non-zero result.
    private static Precision TestRobustCrossProdError(S2Point a, S2Point b)
    {
        // The following error bound includes the call to RobustCrossProd(), the
        // call to ExactCrossProd() (which has a small amount of error due to rounding
        // to a representable result), and two calls to Normalize().  Altogether this
        // adds up to 9 * DBL_ERR == 4.5 * DBL_EPSILON.
        //
        // The maximum error observed during testing is less than 1/3 of this value.
        S1ChordAngle kMaxError = new(S2.kRobustCrossProdErrorS1Angle +
                               kExactCrossProdError +
                         S1Angle.FromRadians(2 * S2Pred.DBL_ERR));

        S2Point result = S2.RobustCrossProd(a, b).Normalize();

        // Test that RobustCrossProd is consistent with S2::Sign.  We don't just check
        // that S2::Sign(a, b, result) == 1 since that would only test whether the two
        // functions agree on the direction of the cross product to within 90 degrees.
        // Instead we check that the sign of RobustCrossProd(a, b).DotProd(c) matches
        // S2::Sign(a, b, c) unless "c" is within kRobustCrossProdError of the plane
        // containing "a" and "b".  We test this by generating 2 pairs of vectors that
        // straddle the plane containing "a" and "b", where the two pairs are spaced
        // 90 degrees apart around the origin.
        S2Point offset = S2.kRobustCrossProdError * result;
        S2Point a90 = result.CrossProd(a);
        Assert.Equal(S2Pred.Sign(a, b, result), 1);
        Assert.True(result.DotProd(a + offset) > 0);
        Assert.True(result.DotProd(a - offset) < 0);
        Assert.True(result.DotProd(a90 + offset) > 0);
        Assert.True(result.DotProd(a90 - offset) < 0);

        // The following identities are true unless a, b are linearly dependent.
        bool have_exact = !S2Pred.IsZero(a.ToExact().CrossProd(b.ToExact()));
        if (have_exact)
        {
            Assert.Equal(S2.RobustCrossProd(-a, b).Normalize(), -result);
            Assert.Equal(S2.RobustCrossProd(a, -b).Normalize(), -result);
        }
        S2Point result_exact;
        if (a == b)
        {
            // This case isn't handled by ExactCrossProd() because it is supposed to be
            // handled long before we get to that point.
            result_exact = S2.Ortho(a).Normalize();
        }
        else
        {
            // Note that ExactCrossProd returns an exact *or* symbolic result.
            result_exact = ExactCrossProd(a, b).Normalize();

            // The following identity is true whenever a != b.
            Assert.Equal(S2.RobustCrossProd(b, a).Normalize(), -result);
        }

        S2Point tmp_result_ld;
        bool have_ld = GetStableCrossProd(a.ToLD(), b.ToLD(), out tmp_result_ld);

        S2Point result_dbl;
        bool have_dbl = GetStableCrossProd(a, b, out result_dbl);
        if (have_dbl)
        {
            result_dbl = result_dbl.Normalize();
            Assert.True(have_ld);
            Assert.Equal(result_dbl, result);
            Assert.True(S2Pred.CompareDistance(result_dbl, result_exact, kMaxError) < 0);
            return Precision.DOUBLE;
        }
        else if (have_ld)
        {
            S2Point result_ld = tmp_result_ld.Normalize();
            Assert.True(have_exact);
            Assert.Equal(result_ld, result);
            Assert.True(S2Pred.CompareDistance(result_ld, result_exact, kMaxError) < 0);
            return Precision.LONG_DOUBLE;
        }
        else
        {
            Assert.Equal(result_exact, result);
            return have_exact ? Precision.EXACT : Precision.SYMBOLIC;
        }
    }

    private static void TestRobustCrossProd(S2Point a, S2Point b, S2Point expected_result, Precision expected_prec)
    {
        Assert.Equal(S2Pred.Sign(a, b, expected_result), 1);
        Assert.Equal(S2.RobustCrossProd(a, b).Normalize(), expected_result);
        Assert.Equal(TestRobustCrossProdError(a, b), expected_prec);
    }

    [Fact]
    public void Test_S2_RobustCrossProdCoverage()
    {
        // Note that RobustCrossProd() returns a non-unit length result.  In "double"
        // and "long double" precision the magnitude is twice the usual cross product.
        // In exact precision it is equal to the usual cross product unless the
        // magnitude is so small that calling Normalize() would lose precision, in
        // which case it is scaled up appropriately.  When symbolic perturbations are
        // used, the result magnitude is arbitrary.  For this reason we only check
        // the result after calling Normalize().
        TestRobustCrossProd(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
                        new S2Point(0, 0, 1), Precision.DOUBLE);
        TestRobustCrossProd(new S2Point(20 * S2Pred.DBL_ERR, 1, 0), new S2Point(0, 1, 0),
                            new S2Point(0, 0, 1), Precision.DOUBLE);
        TestRobustCrossProd(new S2Point(16 * S2Pred.DBL_ERR, 1, 0), new S2Point(0, 1, 0),
                            new S2Point(0, 0, 1), Precision.LONG_DOUBLE);

        // The following bounds hold for both 80-bit and "double double" precision.
        TestRobustCrossProd(new S2Point(5 * S2Pred.LD_ERR, 1, 0), new S2Point(0, 1, 0),
                            new S2Point(0, 0, 1), Precision.LONG_DOUBLE);
        TestRobustCrossProd(new S2Point(4 * S2Pred.LD_ERR, 1, 0), new S2Point(0, 1, 0),
                            new S2Point(0, 0, 1), Precision.EXACT);

        // Test that exact results are scaled up when they would be too small.
        TestRobustCrossProd(new S2Point(5e-324, 1, 0), new S2Point(0, 1, 0),
                            new S2Point(0, 0, 1), Precision.EXACT);


        // And even when the exact cross product underflows in double precision.
        TestRobustCrossProd(new S2Point(5e-324, 1, 0), new S2Point(5e-324, 1 - S2Pred.DBL_ERR, 0),
                            new S2Point(0, 0, -1), Precision.EXACT);

        // Test symbolic results.
        TestRobustCrossProd(new S2Point(1, 0, 0), new S2Point(1 + S2.DoubleEpsilon, 0, 0),
                            new S2Point(0, 1, 0), Precision.SYMBOLIC);
        TestRobustCrossProd(new S2Point(0, 1 + S2.DoubleEpsilon, 0), new S2Point(0, 1, 0),
                            new S2Point(1, 0, 0), Precision.SYMBOLIC);
        TestRobustCrossProd(new S2Point(0, 0, 1), new S2Point(0, 0, -1),
                            new S2Point(-1, 0, 0), Precision.SYMBOLIC);

        // Finally, test some symbolic perturbation cases that can't happen in
        // practice but that are implemented for completeness.
        Assert.Equal(SymbolicCrossProd(new S2Point(-1, 0, 0), new S2Point(0, 0, 0)),
                  new S2Point(0, 1, 0));
        Assert.Equal(SymbolicCrossProd(new S2Point(0, 0, 0), new S2Point(0, -1, 0)),
                  new S2Point(1, 0, 0));
        Assert.Equal(SymbolicCrossProd(new S2Point(0, 0, 0), new S2Point(0, 0, -1)),
                  new S2Point(-1, 0, 0));
    }

    [Fact]
    public void Test_S2_SymbolicCrossProdConsistentWithSign()
    {
        // Tests that RobustCrossProd() is consistent with s2pred::Sign() when
        // symbolic perturbations are necessary.  This implies that the two vectors
        // A and B must be linearly dependent.  We test all possible orderings of
        // the components of vector A = (x, y, z) and all possible orderings of the
        // two vectors A and B.
        foreach (double x in new[] { -1, 0, 1 })
        {
            foreach (double y in new[] { -1, 0, 1 })
            {
                foreach (double z in new[] { -1, 0, 1 })
                {
                    var a = new S2Point(x, y, z).Normalize();
                    if (a == new S2Point()) continue;
                    foreach (double scale in new[] { -1.0, 1 - S2Pred.DBL_ERR, 1 + 2 * S2Pred.DBL_ERR })
                    {
                        S2Point b = scale * a;
                        Assert.True(b.IsUnitLength());
                        Assert.True(S2Pred.Sign(a, b, S2.RobustCrossProd(a, b).Normalize()) > 0);
                    }
                }
            }
        }
    }

    [Fact]
    public void Test_S2_RobustCrossProdMagnitude()
    {
        // Test that angles can be measured between vectors returned by
        // RobustCrossProd without loss of precision due to underflow.  The following
        // magnitude (1e-100) is demonstrates that it is not sufficient to simply
        // ensure that the result of RobustCrossProd() is normalizable, since in this
        // case Normalize() could be called on the result of both cross products
        // without precision loss but it is still not possible to measure the angle
        // between them without scaling the result.
        Assert.Equal(S2.RobustCrossProd(new S2Point(1, 0, 0),
                                  new S2Point(1, 1e-100, 0)).
              Angle(S2.RobustCrossProd(new S2Point(1, 0, 0),
                                        new S2Point(1, 0, 1e-100))),
              S2.M_PI_2);

        // Verify that this is true even when symbolic perturbations are necessary.
        Assert.Equal(S2.RobustCrossProd(new S2Point(-1e-100, 0, 1),
                                      new S2Point(1e-100, 0, -1)).
                  Angle(S2.RobustCrossProd(new S2Point(0, -1e-100, 1),
                                            new S2Point(0, 1e-100, -1))),
                  S2.M_PI_2);
    }

    // Chooses a random S2Point that is often near the intersection of one of the
    // coodinates planes or coordinate axes with the unit sphere.  (It is possible
    // to represent very small perturbations near such points.)
    private static S2Point ChoosePoint()
    {
        var x = S2Testing.RandomPoint().ToArray();
        for (int i = 0; i < 3; ++i)
        {
            if (S2Testing.Random.OneIn(4))
            {  // Denormalized
                x[i] *= Math.Pow(2, -1022 - 53 * S2Testing.Random.RandDouble());
            }
            else if (S2Testing.Random.OneIn(3))
            {  // Zero when squared
                x[i] *= Math.Pow(2, -511 - 511 * S2Testing.Random.RandDouble());
            }
            else if (S2Testing.Random.OneIn(2))
            {  // Simply small
                x[i] *= Math.Pow(2, -100 * S2Testing.Random.RandDouble());
            }
        }
        var xx = new S2Point(x);
        if (xx.Norm2() >= MathUtils.Ldexp(1, -968))
        {
            return xx.Normalize();
        }
        return ChoosePoint();
    }

    // Perturbs the length of the given point slightly while ensuring that it still
    // satisfies the conditions of S2::IsUnitLength().
    private static S2Point PerturbLength(S2Point p)
    {
        S2Point q = p * S2Testing.Random.UniformDouble(1 - 2 * S2.DoubleEpsilon,
                                                     1 + 2 * S2.DoubleEpsilon);
        // S2::IsUnitLength() is not an exact test, since it tests the length using
        // ordinary double-precision arithmetic and allows for errors in its own
        // calculations.  This is fine for most purposes, but since here we are
        // explicitly trying to push the limits we instead use exact arithmetic to
        // test whether the point is within the error tolerance guaranteed by
        // Normalize().
        if (Math.Abs(q.ToExact().Norm2() - 1) <= 4 * S2.DoubleEpsilon) return q;
        return p;
    }

    [Fact]
    public void Test_S2_RobustCrossProdError()
    {
        // We repeatedly choose two points (usually linearly dependent, or almost so)
        // and measure the distance between the cross product returned by
        // RobustCrossProd and one returned by ExactCrossProd.
        PrecisionStats stats = new();
        for (int iter = 0; iter < 5000; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Point a, b;
            do
            {
                a = PerturbLength(ChoosePoint());
                S2Point dir = ChoosePoint();
                S1Angle r = S1Angle.FromRadians(S2.M_PI_2 * Math.Pow(2, -53 * S2Testing.Random.RandDouble()));
                // Occasionally perturb the point by a tiny distance.
                if (S2Testing.Random.OneIn(3)) r *= Math.Pow(2, -1022 * S2Testing.Random.RandDouble());
                b = PerturbLength(S2.InterpolateAtDistance(r, a, dir));
                if (S2Testing.Random.OneIn(2)) b = -b;
            } while (a == b);
            var prec = TestRobustCrossProdError(a, b);
            stats.Tally(prec);
        }
        // These stats are very heavily skewed towards degenerate cases.
        System.Diagnostics.Debug.WriteLine($"(ERROR) {stats}");
    }

    [Fact]
    public void Test_S2_AngleContainsVertex()
    {
        S2Point a = new(1, 0, 0), b = new(0, 1, 0);
        S2Point ref_b = S2.RefDir(b);

        // Degenerate angle ABA.
        Assert.False(S2.AngleContainsVertex(a, b, a));

        // An angle where A == RefDir(B).
        Assert.True(S2.AngleContainsVertex(ref_b, b, a));

        // An angle where C == RefDir(B).
        Assert.False(S2.AngleContainsVertex(a, b, ref_b));

        // Verify that when a set of polygons tile the region around the vertex,
        // exactly one of those polygons contains the vertex.
        var points = S2Testing.MakeRegularPoints(b, S1Angle.FromDegrees(10), 10);
        S2PointLoopSpan loop = new(points);
        int count = 0;
        for (int i = 0; i < loop.Count; ++i)
        {
            count += S2.AngleContainsVertex(loop[i + 1], b, loop[i]) ? 1 : 0;
        }
        Assert.Equal(count, 1);
    }

    // CrossingSign, VertexCrossing, SignedVertexCrossing, and EdgeOrVertexCrossing
    // are tested in s2edge_crosser_test.cc.

    [Fact]
    public void Test_S2_IntersectionError()
    {
        // We repeatedly construct two edges that cross near a random point "p", and
        // measure the distance from the actual intersection point "x" to the
        // exact intersection point and also to the edges.

        var tally_ = GetNewTally();
        S1Angle max_point_dist = S1Angle.Zero, max_edge_dist = S1Angle.Zero;
        for (int iter = 0; iter < 5000; ++iter)
        {
            // We construct two edges AB and CD that intersect near "p".  The angle
            // between AB and CD (expressed as a slope) is chosen randomly between
            // 1e-15 and 1e15 such that its logarithm is uniformly distributed.
            // Similarly, two edge lengths approximately between 1e-15 and 1 are
            // chosen.  The edge endpoints are chosen such that they are often very
            // close to the other edge (i.e., barely crossing).  Taken together this
            // ensures that we test both long and very short edges that intersect at
            // both large and very small angles.
            //
            // Sometimes the edges we generate will not actually cross, in which case
            // we simply try again.
            S2Testing.GetRandomFrame(out var p, out var d1, out var d2);
            double slope = 1e-15 * Math.Pow(1e30, S2Testing.Random.RandDouble());
            d2 = (d1 + slope * d2).Normalize();
            S2Point a, b, c, d;
            do
            {
                double ab_len = Math.Pow(1e-15, S2Testing.Random.RandDouble());
                double cd_len = Math.Pow(1e-15, S2Testing.Random.RandDouble());
                double a_fraction = Math.Pow(1e-5, S2Testing.Random.RandDouble());
                if (S2Testing.Random.OneIn(2)) a_fraction = 1 - a_fraction;
                double c_fraction = Math.Pow(1e-5, S2Testing.Random.RandDouble());
                if (S2Testing.Random.OneIn(2)) c_fraction = 1 - c_fraction;
                a = (p - a_fraction * ab_len * d1).Normalize();
                b = (p + (1 - a_fraction) * ab_len * d1).Normalize();
                c = (p - c_fraction * cd_len * d2).Normalize();
                d = (p + (1 - c_fraction) * cd_len * d2).Normalize();
            } while (S2.CrossingSign(a, b, c, d) <= 0);

            // Each constructed edge should be at most 1.5 * S2Constants.DoubleEpsilon away from the
            // original point P.
            Assert.True(S2.GetDistance(p, a, b) <=
                  S1Angle.FromRadians(1.5 * S2.DoubleEpsilon) + kGetDistanceAbsError);
            Assert.True(S2.GetDistance(p, c, d) <=
                      S1Angle.FromRadians(1.5 * S2.DoubleEpsilon) + kGetDistanceAbsError);

            // Verify that the expected intersection point is close to both edges and
            // also close to the original point P.  (It might not be very close to P
            // if the angle between the edges is very small.)
            S2Point expected = GetIntersectionExact(a, b, c, d);
            Assert.True(S2.GetDistance(expected, a, b) <=
                      S1Angle.FromRadians(3 * S2.DoubleEpsilon) + kGetDistanceAbsError);
            Assert.True(S2.GetDistance(expected, c, d) <=
                      S1Angle.FromRadians(3 * S2.DoubleEpsilon) + kGetDistanceAbsError);
            Assert.True(new S1Angle(expected, p) <=
                      S1Angle.FromRadians(3 * S2.DoubleEpsilon / slope) +
                      S2.kIntersectionErrorS1Angle);

            // Now we actually test the GetIntersection() method.
            S2Point actual = S2.GetIntersection(a, b, c, d, tally_);
            S1Angle dist_ab = S2.GetDistance(actual, a, b);
            S1Angle dist_cd = S2.GetDistance(actual, c, d);
            Assert.True(dist_ab <= S2.kIntersectionErrorS1Angle + kGetDistanceAbsError);
            Assert.True(dist_cd <= S2.kIntersectionErrorS1Angle + kGetDistanceAbsError);
            max_edge_dist = new[] { max_edge_dist, dist_ab, dist_cd }.Max();
            S1Angle point_dist = new(expected, actual);
            Assert.True(point_dist <= S2.kIntersectionErrorS1Angle);
            max_point_dist = new[] { max_point_dist, point_dist }.Max();
        }
        PrintIntersectionStats(tally_);
    }

    [Fact]
    public void Test_S2_GrazingIntersections()
    {
        // This test choose 5 points along a great circle (i.e., as collinear as
        // possible), and uses them to construct an edge AB and a triangle CDE such
        // that CD and CE both cross AB.  It then checks that the intersection
        // points returned by GetIntersection() have the correct relative ordering
        // along AB (to within kIntersectionError).
        var tally_ = GetNewTally();
        for (int iter = 0; iter < 1000; ++iter)
        {
            S2Testing.GetRandomFrame(out var x, out var y, out _);
            S2Point a, b, c, d, e, ab;
            do
            {
                a = ChooseSemicirclePoint(x, y);
                b = ChooseSemicirclePoint(x, y);
                c = ChooseSemicirclePoint(x, y);
                d = ChooseSemicirclePoint(x, y);
                e = ChooseSemicirclePoint(x, y);
                ab = (a - b).CrossProd(a + b);
            } while (ab.Norm() < 50 * S2.DoubleEpsilon ||
                     S2.CrossingSign(a, b, c, d) <= 0 ||
                     S2.CrossingSign(a, b, c, e) <= 0);
            S2Point xcd = S2.GetIntersection(a, b, c, d, tally_);
            S2Point xce = S2.GetIntersection(a, b, c, e, tally_);
            // Essentially this says that if CDE and CAB have the same orientation,
            // then CD and CE should intersect along AB in that order.
            ab = ab.Normalize();
            if (new S1Angle(xcd, xce) > 2 * S2.kIntersectionErrorS1Angle)
            {
                Assert.Equal(S2Pred.Sign(c, d, e) == S2Pred.Sign(c, a, b),
                          S2Pred.Sign(ab, xcd, xce) > 0);
            }
        }
        PrintIntersectionStats(tally_);
    }

    [Fact]
    public void Test_S2_ExactIntersectionUnderflow()
    {
        // Tests that a correct intersection is computed even when two edges are
        // exactly collinear and the normals of both edges underflow in double
        // precision when normalized (see the ToVector3 function for details).
        S2Point a0 = new(1, 0, 0), a1 = new(1, 2e-300, 0);
        S2Point b0 = new(1, 1e-300, 0), b1 = new(1, 3e-300, 0);
        Assert.Equal(new S2Point(1, 1e-300, 0), S2.GetIntersection(a0, a1, b0, b1, null));
    }

    [Fact]
    public void Test_S2_ExactIntersectionSign()
    {
        // Tests that a correct intersection is computed even when two edges are
        // exactly collinear and both edges have nearly antipodal endpoints.
        S2Point a0 = new(-1, -1.6065916409055676e-10, 0);
        S2Point a1 = new(1, 0, 0);
        S2Point b0 = new(1, -4.7617930898495072e-13, 0);
        S2Point b1 = new(-1, 1.2678623820887328e-09, 0);
        Assert.Equal(new S2Point(1, -4.7617930898495072e-13, 0),
                  S2.GetIntersection(a0, a1, b0, b1, null));
    }

    [Fact]
    public void Test_S2_GetIntersectionInvariants()
    {
        // Test that the result of GetIntersection does not change when the edges
        // are swapped and/or reversed.  The number of iterations is high because it
        // is difficult to generate test cases that show that CompareEdges() is
        // necessary and correct, for example.
        int kIters = 50000;
#if DEBUG
        kIters = 5000;
#endif
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Point a, b, c, d;
            do
            {
                // GetIntersectionStable() sorts the two edges by length, soruct
                // edges (a,b) and (c,d) that cross and have exactly the same length.
                // This can be done by swapping the "x" and "y" coordinates.
                // [Swapping other coordinate pairs doesn't work because it changes the
                // order of addition in Norm2() == (x**2 + y**2) + z**2.]
                a = c = S2Testing.RandomPoint();
                b = d = S2Testing.RandomPoint();
                c = new S2Point(c.Y, c.X, c.Z);
                d = new S2Point(d.Y, d.X, d.Z);
            } while (S2.CrossingSign(a, b, c, d) <= 0);
            Assert.Equal((a - b).Norm2(), (c - d).Norm2());

            // Now verify that GetIntersection returns exactly the same result when
            // the edges are swapped and/or reversed.
            S2Point result = S2.GetIntersection(a, b, c, d, null);
            if (S2Testing.Random.OneIn(2)) { var tmp = a; a = b; b = tmp; }
            if (S2Testing.Random.OneIn(2)) { var tmp = c; c = d; d = tmp; }
            if (S2Testing.Random.OneIn(2)) { var tmp = a; b = d; d = tmp; }
            Assert.Equal(result, S2.GetIntersection(a, b, c, d, null));
        }
    }

    // This returns the true intersection point of two line segments (a0,a1) and
    // (b0,b1), with a relative error of at most S2Constants.DoubleEpsilon in each coordinate
    // (i.e., one ulp, or twice the double precision rounding error).
    private static S2Point GetIntersectionExact(S2Point a0, S2Point a1, S2Point b0, S2Point b1)
    {
        S2Point x = GetIntersectionExact(a0, a1, b0, b1);
        if (x.DotProd((a0 + a1) + (b0 + b1)) < 0) x = -x;
        return x;
    }

    // Chooses a point in the XY plane that is separated from X by at least 1e-15
    // (to avoid choosing too many duplicate points) and by at most Pi/2 - 1e-3
    // (to avoid nearly-diametric edges, since the test below is not sophisticated
    // enough to test such edges).
    private static S2Point ChooseSemicirclePoint(S2Point x, S2Point y)
    {
        double sign = (2 * S2Testing.Random.Uniform(2)) - 1;
        return (x + sign * 1e3 * Math.Pow(1e-18, S2Testing.Random.RandDouble()) * y).Normalize();
    }

    // This method records statistics about the intersection methods used by
    // GetIntersection().
    private void PrintIntersectionStats(Dictionary<IntersectionMethod, int> tally_)
    {
        int total = 0;
        var totals = GetNewTally();
        foreach (var key in tally_.Keys)
        {
            total += tally_[key];
            totals[key] = total;
        }
        _logger.WriteLine($"{"Method":10} {"Successes":16} {"Attempts":16}  {"Rate":6}");
        foreach (var key in tally_.Keys)
        {
            var tal = tally_[key];
            if (tal == 0) continue;
            var name = GetIntersectionMethodName(key);
            var tot = totals[key];
            var suc = 100.0 * tal / total;
            var att = 100.0 * tot / total;
            var rat = 100.0 * tal / tot;
            _logger.WriteLine($"{name:10} {tal:9} {suc:5.1f}% {tot:9} {att:5.1f}%  {rat:5.1f}%");
        }
        foreach (var key in tally_.Keys)
            tally_[key] = 0;
    }
}
