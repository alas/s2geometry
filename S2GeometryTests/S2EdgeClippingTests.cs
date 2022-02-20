using System.Text;
using static S2Geometry.S2EdgeClipping;

namespace S2Geometry;

public class S2EdgeClippingTests
{
    private readonly ITestOutputHelper _logger;

    public S2EdgeClippingTests(ITestOutputHelper logger) { _logger = logger; }

    private void TestFaceClipping(S2Point a_raw, S2Point b_raw)
    {
        S2Point a = a_raw.Normalize();
        S2Point b = b_raw.Normalize();

        // First we test GetFaceSegments.
        FaceSegmentVector segments = new();
        GetFaceSegments(a, b, segments);
        int n = segments.Count;
        Assert.True(n >= 1);

        var msg = new StringBuilder($"\nA={a_raw}\nB={b_raw}\nN={S2.RobustCrossProd(a, b)}\nSegments:\n");
        int i1 = 0;
        foreach (var s in segments)
        {
            msg.AppendLine($"{i1++}: face={s.face}, a={s.a}, b={s.b}");
        }
        _logger.WriteLine(msg.ToString());

        R2Rect biunit = new(new R1Interval(-1, 1), new R1Interval(-1, 1));
        var kErrorRadians = kFaceClipErrorRadians;

        // The first and last vertices should approximately equal A and B.
        Assert.True(a.Angle(S2.FaceUVtoXYZ(segments[0].face, segments[0].a)) <=
                  kErrorRadians);
        Assert.True(b.Angle(S2.FaceUVtoXYZ(segments[n - 1].face, segments[n - 1].b)) <=
                  kErrorRadians);

        S2Point norm = S2.RobustCrossProd(a, b).Normalize();
        S2Point a_tangent = norm.CrossProd(a);
        S2Point b_tangent = b.CrossProd(norm);
        for (int i = 0; i < n; ++i)
        {
            // Vertices may not protrude outside the biunit square.
            Assert.True(biunit.Contains(segments[i].a));
            Assert.True(biunit.Contains(segments[i].b));
            if (i == 0) continue;

            // The two representations of each interior vertex (on adjacent faces)
            // must correspond to exactly the same S2Point.
            Assert.NotEqual(segments[i - 1].face, segments[i].face);
            Assert.Equal(S2.FaceUVtoXYZ(segments[i - 1].face, segments[i - 1].b),
                      S2.FaceUVtoXYZ(segments[i].face, segments[i].a));

            // Interior vertices should be in the plane containing A and B, and should
            // be contained in the wedge of angles between A and B (i.e., the dot
            // products with a_tangent and b_tangent should be non-negative).
            S2Point p = S2.FaceUVtoXYZ(segments[i].face, segments[i].a).Normalize();
            Assert.True(Math.Abs(p.DotProd(norm)) <= kErrorRadians);
            Assert.True(p.DotProd(a_tangent) >= -kErrorRadians);
            Assert.True(p.DotProd(b_tangent) >= -kErrorRadians);
        }

        // Now we test ClipToPaddedFace (sometimes with a padding of zero).  We do
        // this by defining an (x,y) coordinate system for the plane containing AB,
        // and converting points along the great circle AB to angles in the range
        // [-Pi, Pi].  We then accumulate the angle intervals spanned by each
        // clipped edge; the union over all 6 faces should approximately equal the
        // interval covered by the original edge.
        double padding = S2Testing.Random.OneIn(10) ? 0.0 : 1e-10 * Math.Pow(1e-5, S2Testing.Random.RandDouble());
        S2Point x_axis = a, y_axis = a_tangent;
        S1Interval expected_angles = new(0, a.Angle(b));
        S1Interval max_angles = expected_angles.Expanded(kErrorRadians);
        S1Interval actual_angles = new();
        for (int face = 0; face < 6; ++face)
        {
            if (ClipToPaddedFace(a, b, face, padding, out var a_uv, out var b_uv))
            {
                S2Point a_clip = S2.FaceUVtoXYZ(face, a_uv).Normalize();
                S2Point b_clip = S2.FaceUVtoXYZ(face, b_uv).Normalize();
                Assert.True(Math.Abs(a_clip.DotProd(norm)) <= kErrorRadians);
                Assert.True(Math.Abs(b_clip.DotProd(norm)) <= kErrorRadians);
                if (a_clip.Angle(a) > kErrorRadians)
                {
                    Assert2.DoubleEqual(1 + padding, Math.Max(Math.Abs(a_uv[0]), Math.Abs(a_uv[1])));
                }
                if (b_clip.Angle(b) > kErrorRadians)
                {
                    Assert2.DoubleEqual(1 + padding, Math.Max(Math.Abs(b_uv[0]), Math.Abs(b_uv[1])));
                }
                double a_angle = Math.Atan2(a_clip.DotProd(y_axis), a_clip.DotProd(x_axis));
                double b_angle = Math.Atan2(b_clip.DotProd(y_axis), b_clip.DotProd(x_axis));
                // Rounding errors may cause b_angle to be slightly less than a_angle.
                // We handle this by constructing the interval with FromPointPair(),
                // which is okay since the interval length is much less than Math.PI.
                S1Interval face_angles = S1Interval.FromPointPair(a_angle, b_angle);
                Assert.True(max_angles.Contains(face_angles));
                actual_angles = actual_angles.Union(face_angles);
            }
        }
        Assert.True(actual_angles.Expanded(kErrorRadians).Contains(expected_angles));
    }

    private void TestFaceClippingEdgePair(S2Point a, S2Point b)
    {
        TestFaceClipping(a, b);
        TestFaceClipping(b, a);
    }

    // This function is designed to choose line segment endpoints that are difficult
    // to handle correctly.  Given two adjacent cube vertices P and Q, it returns
    // either an edge midpoint, face midpoint, or corner vertex along the edge PQ
    // and then perturbs it slightly.  It also sometimes returns a random point from
    // anywhere on the sphere.
    private S2Point PerturbedCornerOrMidpoint(S2Point p, S2Point q)
    {
        S2Point a = (S2Testing.Random.Uniform(3) - 1) * p + (S2Testing.Random.Uniform(3) - 1) * q;
        if (S2Testing.Random.OneIn(10))
        {
            // This perturbation often has no effect except on coordinates that are
            // zero, in which case the perturbed value is so small that operations on
            // it often result in underflow.
            a += Math.Pow(1e-300, S2Testing.Random.RandDouble()) * S2Testing.RandomPoint();
        }
        else if (S2Testing.Random.OneIn(2))
        {
            // For coordinates near 1 (say > 0.5), this perturbation yields values
            // that are only a few representable values away from the initial value.
            a += 4 * S2.DoubleEpsilon * S2Testing.RandomPoint();
        }
        else
        {
            // A perturbation whose magnitude is in the range [1e-25, 1e-10].
            a += 1e-10 * Math.Pow(S2.DoubleError, S2Testing.Random.RandDouble()) * S2Testing.RandomPoint();
        }
        if (a.Norm2() < S2.DoubleMinNorm)
        {
            // If a.Norm2 is denormalized, Normalize() loses too much precision.
            return PerturbedCornerOrMidpoint(p, q);
        }
        return a;
    }

    [Fact]
    public void Test_S2_FaceClipping()
    {
        // Start with a few simple cases.
        // An edge that is entirely contained within one cube face:
        TestFaceClippingEdgePair(new S2Point(1, -0.5, -0.5), new S2Point(1, 0.5, 0.5));
        // An edge that crosses one cube edge:
        TestFaceClippingEdgePair(new S2Point(1, 0, 0), new S2Point(0, 1, 0));
        // An edge that crosses two opposite edges of face 0:
        TestFaceClippingEdgePair(new S2Point(0.75, 0, -1), new S2Point(0.75, 0, 1));
        // An edge that crosses two adjacent edges of face 2:
        TestFaceClippingEdgePair(new S2Point(1, 0, 0.75), new S2Point(0, 1, 0.75));
        // An edge that crosses three cube edges (four faces):
        TestFaceClippingEdgePair(new S2Point(1, 0.9, 0.95), new S2Point(-1, 0.95, 0.9));

        // Comprehensively test edges that are difficult to handle, especially those
        // that nearly follow one of the 12 cube edges.
        R2Rect biunit = new(new R1Interval(-1, 1), new R1Interval(-1, 1));
        int kIters = 1000;  // Test passes with 1e6 iterations
        for (int iter = 0; iter < kIters; ++iter)
        {
            _logger.WriteLine($"Iteration {iter}");
            // Choose two adjacent cube corners P and Q.
            int face = S2Testing.Random.Uniform(6);
            int i = S2Testing.Random.Uniform(4);
            int j = (i + 1) & 3;
            S2Point p = S2.FaceUVtoXYZ(face, biunit.GetVertex(i));
            S2Point q = S2.FaceUVtoXYZ(face, biunit.GetVertex(j));

            // Now choose two points that are nearly on the edge PQ, preferring points
            // that are near cube corners, face midpoints, or edge midpoints.
            S2Point a = PerturbedCornerOrMidpoint(p, q);
            S2Point b = PerturbedCornerOrMidpoint(p, q);
            TestFaceClipping(a, b);
        }
    }

    // Choose a random point in the rectangle defined by points A and B, sometimes
    // returning a point on the edge AB or the points A and B themselves.
    private static R2Point ChooseRectPoint(R2Point a, R2Point b)
    {
        if (S2Testing.Random.OneIn(5))
        {
            return S2Testing.Random.OneIn(2) ? a : b;
        }
        else if (S2Testing.Random.OneIn(3))
        {
            return a + S2Testing.Random.RandDouble() * (b - a);
        }
        else
        {
            // a[i] may be >, <, or == b[i], so we write it like this instead
            // of using UniformDouble.
            return new R2Point(a[0] + S2Testing.Random.RandDouble() * (b[0] - a[0]),
                           a[1] + S2Testing.Random.RandDouble() * (b[1] - a[1]));
        }
    }

    // Given a point X on the line AB (which is checked), return the fraction "t"
    // such that x = (1-t)*a + t*b.  Return 0 if A = B.
    private static double GetFraction(R2Point x, R2Point a, R2Point b)
    {
        // A bound for the error in edge clipping plus the error in the calculation
        // below (which is similar to IntersectsRect).
        double kError = kEdgeClipErrorUVDist + kIntersectsRectErrorUVDist;
        if (a == b) return 0.0;
        R2Point dir = R2Point.Normalize(b - a);
        Assert.True(Math.Abs((x - a).DotProd(dir.GetOrtho())) <= kError);
        return (x - a).DotProd(dir);
    }

    // Given a point P representing a possibly clipped endpoint A of an edge AB,
    // verify that "clip" contains P, and that if clipping occurred (i.e., P != A)
    // then P is on the boundary of "clip".
    private static void CheckPointOnBoundary(R2Point p, R2Point a, R2Rect clip)
    {
        Assert.True(clip.Contains(p));
        if (p != a)
        {
            Assert.False(clip.Contains(new R2Point(MathUtils.NextAfter(p[0], a[0]),
                                               MathUtils.NextAfter(p[1], a[1]))));
        }
    }

    // Given an edge AB and a rectangle "clip", verify that IntersectsRect(),
    // ClipEdge(), and ClipEdgeBound() produce consistent results.
    private static void TestClipEdge(R2Point a, R2Point b, R2Rect clip)
    {
        // A bound for the error in edge clipping plus the error in the
        // IntersectsRect calculation below.
        double kError = kEdgeClipErrorUVDist + kIntersectsRectErrorUVDist;
        if (!ClipEdge(a, b, clip, out var a_clipped, out var b_clipped))
        {
            Assert.False(IntersectsRect(a, b, clip.Expanded(-kError)));
        }
        else
        {
            Assert.True(IntersectsRect(a, b, clip.Expanded(kError)));
            // Check that the clipped points lie on the edge AB, and that the points
            // have the expected order along the segment AB.
            Assert.True(GetFraction(a_clipped, a, b) <= GetFraction(b_clipped, a, b));
            // Check that the clipped portion of AB is as large as possible.
            CheckPointOnBoundary(a_clipped, a, clip);
            CheckPointOnBoundary(b_clipped, b, clip);
        }
        // Choose a random initial bound to pass to ClipEdgeBound.
        R2Rect initial_clip = R2Rect.FromPointPair(ChooseRectPoint(a, b),
                                                    ChooseRectPoint(a, b));
        R2Rect bound = GetClippedEdgeBound(a, b, initial_clip);
        if (bound.IsEmpty()) return;  // Precondition of ClipEdgeBound not met
        R2Rect max_bound = bound.Intersection(clip);
        if (!ClipEdgeBound(a, b, clip, ref bound))
        {
            Assert.False(IntersectsRect(a, b, max_bound.Expanded(-kError)));
        }
        else
        {
            Assert.True(IntersectsRect(a, b, max_bound.Expanded(kError)));
            // Check that the bound is as large as possible.
            int ai = (a[0] > b[0]) ? 1 : 0, aj = (a[1] > b[1]) ? 1 : 0;
            CheckPointOnBoundary(bound.GetVertex(ai, aj), a, max_bound);
            CheckPointOnBoundary(bound.GetVertex(1 - ai, 1 - aj), b, max_bound);
        }
    }

    // Given an interval "clip", randomly choose either a value in the interval, a
    // value outside the interval, or one of the two interval endpoints, ensuring
    // that all cases have reasonable probability for any interval "clip".
    private static double ChooseEndpoint(R1Interval clip)
    {
        if (S2Testing.Random.OneIn(5))
        {
            return S2Testing.Random.OneIn(2) ? clip.Lo : clip.Hi;
        }
        else
        {
            return (S2Testing.Random.Uniform(3)) switch
            {
                0 => clip.Lo - S2Testing.Random.RandDouble(),
                1 => clip.Hi + S2Testing.Random.RandDouble(),
                _ => clip.Lo + S2Testing.Random.RandDouble() * clip.GetLength(),
            };
        }
    }

    // Given a rectangle "clip", choose a point that may lie in the rectangle
    // interior, along an extended edge, exactly at a vertex, or in one of the
    // eight regions exterior to "clip" that are separated by its extended edges.
    // Also sometimes return points that are exactly on one of the extended
    // diagonals of "clip".  All cases are reasonably likely to occur for any
    // given rectangle "clip".
    private static R2Point ChooseEndpoint(R2Rect clip)
    {
        if (S2Testing.Random.OneIn(10))
        {
            // Return a point on one of the two extended diagonals.
            int diag = S2Testing.Random.Uniform(2);
            double t = S2Testing.Random.UniformDouble(-1, 2);
            return (1 - t) * clip.GetVertex(diag) + t * clip.GetVertex(diag + 2);
        }
        else
        {
            return new R2Point(ChooseEndpoint(clip[0]), ChooseEndpoint(clip[1]));
        }
    }

    // Given a rectangle "clip", test the S2 edge clipping methods using
    // many edges that are randomly constructed to trigger special cases.
    private void TestEdgeClipping(R2Rect clip)
    {
        int kIters = 1000;  // Test passes with 1e6 iterations
        for (int iter = 0; iter < kIters; ++iter)
        {
            _logger.WriteLine($"Iteration {iter}");
            TestClipEdge(ChooseEndpoint(clip), ChooseEndpoint(clip), clip);
        }
    }

    [Fact]
    public void Test_S2_EdgeClipping()
    {
        // Test clipping against random rectangles.
        for (int i = 0; i < 5; ++i)
        {
            TestEdgeClipping(R2Rect.FromPointPair(
                new R2Point(S2Testing.Random.UniformDouble(-1, 1), S2Testing.Random.UniformDouble(-1, 1)),
                new R2Point(S2Testing.Random.UniformDouble(-1, 1), S2Testing.Random.UniformDouble(-1, 1))));
        }
        // Also clip against one-dimensional, singleton, and empty rectangles.
        TestEdgeClipping(new R2Rect(new R1Interval(-0.7, -0.7), new R1Interval(0.3, 0.35)));
        TestEdgeClipping(new R2Rect(new R1Interval(0.2, 0.5), new R1Interval(0.3, 0.3)));
        TestEdgeClipping(new R2Rect(new R1Interval(-0.7, 0.3), new R1Interval(0, 0)));
        TestEdgeClipping(R2Rect.FromPoint(new R2Point(0.3, 0.8)));
        TestEdgeClipping(R2Rect.Empty);
    }
}
