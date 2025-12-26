namespace S2Geometry;

using static S2BufferOperation;

// A callback that allows adding input to an S2BufferOperation.  This is
// convenient for testing the various input methods (AddPoint, AddShape, etc).
using InputCallback = Action<S2BufferOperation>;

public class S2BufferOperationTests(ITestOutputHelper logger)
{
    private readonly ITestOutputHelper _logger = logger;

    // Convenience function that calls the given lambda expression to add input to
    // an S2BufferOperation and returns the buffered result.
    private static S2LaxPolygonShape DoBuffer(
                    InputCallback input_callback, Options options)
    {
        S2LaxPolygonShape output = new();
        S2BufferOperation op = new(new LaxPolygonLayer(output), options);
        input_callback(op);
        Assert.True(op.Build(out _));// << error;
        /*if (S2_VLOG_IS_ON(1) && output.NumVertices < 1000)
        {
            _logger.WriteLine($"\nS2Polygon Error: {output}\n");
        }*/
        return output;
    }

    // Simpler version that accepts a buffer radius and error fraction.
    private static S2LaxPolygonShape DoBuffer(
                    InputCallback input_callback,
                    S1Angle buffer_radius, double error_fraction)
    {
        Options options = new()
        {
            BufferRadius = buffer_radius,
            ErrorFraction = error_fraction,
        };
        return DoBuffer(input_callback, options);
    }

    // Given a callback that adds empty geometry to the S2BufferOperation,
    // verifies that the result is empty after buffering.
    private static void TestBufferEmpty(InputCallback input)
    {
        // Test code paths that involve buffering by negative, zero, and positive
        // values, and also values where the result is usually empty or full.
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(-200), 0.1).IsEmpty());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(-1), 0.1).IsEmpty());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(0), 0.1).IsEmpty());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(1), 0.1).IsEmpty());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(200), 0.1).IsEmpty());
    }

    [Fact]
    internal void Test_S2BufferOperation_NoInput() =>
        TestBufferEmpty(op => { });

    [Fact]
    internal void Test_S2BufferOperation_EmptyPolyline() =>
        // Note that polylines with 1 vertex are defined to have no edges.
        TestBufferEmpty(op => op.AddPolyline(new List<S2Point>() { new(1, 0, 0) }));

    [Fact]
    internal void Test_S2BufferOperation_EmptyLoop() =>
        TestBufferEmpty(op => op.AddLoop([]));

    [Fact]
    internal void Test_S2BufferOperation_EmptyPointShape() =>
        TestBufferEmpty(op => op.AddShape(new S2PointVectorShape()));

    [Fact]
    internal void Test_S2BufferOperation_EmptyPolylineShape() =>
        TestBufferEmpty(op => { op.AddShape(MakeLaxPolylineOrDie("")); });

    [Fact]
    internal void Test_S2BufferOperation_EmptyPolygonShape() =>
        TestBufferEmpty(op => op.AddShape(MakeLaxPolygonOrDie("")));

    [Fact]
    internal void Test_S2BufferOperation_EmptyShapeIndex() =>
        TestBufferEmpty(op => op.AddShapeIndex(MakeIndexOrDie("# #")));

    [Fact]
    internal void Test_S2BufferOperation_Options()
    {
        // Provide test coverage for `options()`.
        Options options = new(S1Angle.FromRadians(1e-12));
        S2LaxPolygonShape output = new();
        S2BufferOperation op = new(new LaxPolygonLayer(output), options);
        Assert.Equal(options.BufferRadius, op.Options_.BufferRadius);
    }

    [Fact]
    internal void Test_S2BufferOperation_PoorlyNormalizedPoint()
    {
        // Verify that debugging assertions are not triggered when an input point is
        // not unit length (but within the limits guaranteed by S2Point.Normalize).
        //
        // The purpose of this test is not to check that the result is correct
        // (which is done elsewhere), simply that no assertions occur.
        DoBuffer(op =>
        {
            S2Point p = new(1 - 2 * S2.DoubleEpsilon, 0, 0);  // Maximum error allowed.
            Assert.True(p.IsUnitLength());
            op.AddPoint(p);
        }, S1Angle.FromDegrees(1), 0.01);
    }

    // Given a callback that adds the full polygon to the S2BufferOperation,
    // verifies that the result is full after buffering.
    private static void TestBufferFull(InputCallback input)
    {
        // Test code paths that involve buffering by negative, zero, and positive
        // values, and also values where the result is usually empty or full.
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(-200), 0.1).IsFull());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(-1), 0.1).IsFull());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(0), 0.1).IsFull());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(1), 0.1).IsFull());
        Assert.True(DoBuffer(input, S1Angle.FromDegrees(200), 0.1).IsFull());
    }

    [Fact]
    internal void Test_S2BufferOperation_FullPolygonShape() =>
        TestBufferFull(op => op.AddShape(MakeLaxPolygonOrDie("full")));

    [Fact]
    internal void Test_S2BufferOperation_FullShapeIndex() =>
        TestBufferFull(op => op.AddShapeIndex(MakeIndexOrDie("# # full")));

    [Fact]
    internal void Test_S2BufferOperation_PointsAndPolylinesAreRemoved()
    {
        // Test that points and polylines are removed with a negative buffer radius.
        var output = DoBuffer(op => op.AddShapeIndex(MakeIndexOrDie("0:0 # 2:2, 2:3#")), S1Angle.FromDegrees(-1), 0.1);
        Assert.True(output.IsEmpty());
    }

    [Fact]
    internal void Test_S2BufferOperation_BufferedPointsAreSymmetric()
    {
        // Test that points are buffered into regular polygons.  (This is not
        // guaranteed by the API but makes the output nicer to look at. :)
        var output = DoBuffer(op => op.AddPoint(new S2Point(1, 0, 0)), S1Angle.FromDegrees(5), 0.001234567);

        // We use the length of the last edge as our reference length.
        int n = output.NumVertices;
        S1Angle edge_len = new(output.LoopVertex(0, 0), output.LoopVertex(0, n - 1));
        for (int i = 1; i < n; ++i)
        {
            var testAngle = S1Angle.Abs(edge_len - new S1Angle(output.LoopVertex(0, i - 1),
                output.LoopVertex(0, i)));
            Assert.True(testAngle <= S1Angle.FromRadians(1e-14));
        }
    }

    [Fact]
    internal void Test_S2BufferOperation_SetCircleSegments()
    {
        // Test that when a point is buffered with a small radius the number of
        // edges matches options.circle_segments().  (This is not true for large
        // radii because large circles on the sphere curve less than 360 degrees.)
        // Using a tiny radius helps to catch rounding problems.
        Options options = new(S1Angle.FromRadians(1e-12));
        for (int circle_segments = 3; circle_segments <= 20; ++circle_segments)
        {
            options.CircleSegments = circle_segments;
            Assert.Equal(circle_segments, options.CircleSegments);
            var output = DoBuffer(op =>
            {
                op.AddPoint(new S2Point(1, 0, 0));
            }, options);
            Assert.Equal(output.NumVertices, circle_segments);
        }
    }

    [Fact]
    internal void Test_S2BufferOperation_SetSnapFunction()
    {
        // Verify that the snap function is passed through to S2Builder.
        // We use a buffer radius of zero to make the test simpler.
        Options options = new()
        {
            SnapFunction_ = new IntLatLngSnapFunction(0),
        };
        var output = DoBuffer(op => op.AddPoint(MakePointOrDie("0.1:-0.4")), options);
        Assert.Equal(output.NumVertices, 1);
        Assert.Equal(output.LoopVertex(0, 0), MakePointOrDie("0:0"));
    }

    [Fact]
    internal void Test_S2BufferOperation_NegativeBufferRadiusMultipleLayers()
    {
        // Verify that with a negative buffer radius, at most one polygon layer is
        // allowed.
        S2LaxPolygonShape output = new();
        S2BufferOperation op = new(
            new LaxPolygonLayer(output),
            new Options(S1Angle.FromRadians(-1)));
        op.AddLoop(new S2PointLoopSpan(ParsePointsOrDie("0:0, 0:1, 1:0")));
        op.AddShapeIndex(MakeIndexOrDie("# # 2:2, 2:3, 3:2"));
        Assert.False(op.Build(out S2Error error));
        Assert.Equal(error.Code, S2ErrorCode.FAILED_PRECONDITION);
    }

    // If buffer_radius > max_error, tests that "output" contains "input".
    // If buffer_radius < -max_error tests that "input" contains "output".
    // Otherwise does nothing.
    private static void TestContainment(MutableS2ShapeIndex input,
                    MutableS2ShapeIndex output,
                    S1Angle buffer_radius, S1Angle max_error)
    {
        S2BooleanOperation.Options options = new()
        {
            PolygonModel_ = S2BooleanOperation.PolygonModel.CLOSED,
            PolylineModel_ = S2BooleanOperation.PolylineModel.CLOSED
        };
        if (buffer_radius > max_error)
        {
            // For positive buffer radii, the output should contain the input.
            Assert.True(S2BooleanOperation.Contains(output, input, options));
        }
        else if (buffer_radius < -max_error)
        {
            // For negative buffer radii, the input should contain the output.
            Assert.True(S2BooleanOperation.Contains(input, output, options));
        }
    }

    // Tests that the minimum distance from the boundary of "output" to the
    // boundary of "input" is at least "min_dist" using exact predicates.
    private static void TestMinimumDistance(MutableS2ShapeIndex input,
                    MutableS2ShapeIndex output,
                    S1ChordAngle min_dist)
    {
        if (min_dist == S1ChordAngle.Zero) return;

        // We do one query to find the edges of "input" that might be too close to
        // "output", then for each such edge we do another query to find the edges
        // of "output" that might be too close to it.  Then we check the distance
        // between each edge pair (A, B) using exact predicates.

        // We make the distance limit big enough to find all edges whose true
        // distance might be less than "min_dist".
        S2ClosestEdgeQuery.Options query_options = new()
        {
            IncludeInteriors = false,
            MaxDistance = min_dist.PlusError(S2.GetUpdateMinDistanceMaxError(min_dist)),
        };

        S2ClosestEdgeQuery in_query = new(input, query_options);
        S2ClosestEdgeQuery.ShapeIndexTarget out_target = new(output)
        {
            IncludeInteriors = false,
        };
        S2ClosestEdgeQuery out_query = new(output, query_options);
        foreach (var in_result in in_query.FindClosestEdges(out_target))
        {
            var a = input.Shape(in_result.ShapeId).GetEdge(in_result.EdgeId);
            S2ClosestEdgeQuery.EdgeTarget in_target = new(a.V0, a.V1);
            foreach (var out_result in out_query.FindClosestEdges(in_target))
            {
                var b = output.Shape(out_result.ShapeId).GetEdge(out_result.EdgeId);
                Assert.True(S2Pred.CompareEdgePairDistance(
                    a.V0, a.V1, b.V0, b.V1, min_dist) >= 0);
            }
        }
    }

    // Tests that the Hausdorff distance from the boundary of "output" to the
    // boundary of "input" is at most (1 + error_fraction) * buffer_radius.  The
    // implementation approximates this by measuring the distance at a set of
    // points along the boundary of "output".
    private static void TestHausdorffDistance(MutableS2ShapeIndex input,
                    MutableS2ShapeIndex output,
                    S1ChordAngle max_dist)
    {
        S2ClosestEdgeQuery.Options query_options = new()
        {
            IncludeInteriors = false,
            MaxDistance = max_dist.PlusError(S2.GetUpdateMinDistanceMaxError(max_dist)),
        };

        S2ClosestEdgeQuery in_query = new(input, query_options);
        foreach (var out_shape in output)
        {
            for (int i = 0; i < out_shape.NumEdges(); ++i)
            {
                var e = out_shape.GetEdge(i);
                // Measure the distance at 5 points along the edge.
                for (double t = 0; t <= 1.0; t += 0.25)
                {
                    S2Point b = S2.Interpolate(e.V0, e.V1, t);
                    S2ClosestEdgeQuery.PointTarget out_target = new(b);
                    // We check the distance bound using exact predicates.
                    foreach (var in_result in in_query.FindClosestEdges(out_target))
                    {
                        var a = input.Shape(in_result.ShapeId).GetEdge(in_result.EdgeId);
                        Assert.True(S2Pred.CompareEdgeDistance(b, a.V0, a.V1, max_dist) <= 0);
                    }
                }
            }
        }
    }

    // Buffers the given input with the given buffer_radius and error_fraction and
    // verifies that the output is correct.
    private void TestBuffer(MutableS2ShapeIndex input, S1Angle buffer_radius,
                    double error_fraction)
    {
        // Ideally we would verify the correctness of buffering as follows.  Suppose
        // that B = Buffer(A, r) and let ~X denote the complement of region X.  Then
        // if r > 0, we would verify:
        //
        //   1a. Minimum distance between ~B and A >= r_min
        //   1b. Directed Hausdorff distance from B to A <= r_max
        //
        // Buffering A by r < 0 is equivalent to buffering ~A by |r|, so instead we
        // would verify the following (where r_min and r_max are based on |r|):
        //
        //   2a. Minimum distance between B and ~A >= r_min
        //   2b. Directed Hausdorff distance from ~B to ~A <= r_max
        //
        // Conditions 1a and 2a can be implemented as follows:
        //
        //   1a*: B.Contains(A) && minimum distance between @B and @A >= r_min
        //   2a*: A.Contains(B) && minimum distance between @B and @A >= r_min
        //
        // Note that if r_min <= 0 then there is nothing to be tested, since the
        // containment condition may not hold.  (Note that even when the specified
        // buffer radius is zero, edges can move slightly when crossing edges are
        // split during the snapping step.)  The correct approach would be to test
        // instead that the directed Hausdorff distance from A to ~B is at most
        // -r_min, but Hausdorff distance is not yet implemented.
        //
        // Similarly, conditions 1b and 2b need to be approximated because Hausdorff
        // distance is not yet implemented.  We do this by measuring the distance at
        // a set of points on the boundary of B:
        //
        //   1b*: Minimum distance from P to @A <= r_max for a set of points P on @B
        //   2b*: Minimum distance from P to @A <= r_max for a set of points P on @B
        //
        // This is not perfect (e.g., it won't detect cases where an entire boundary
        // loop of B is missing, such as returning a disc in the place of an
        // annulus) but it is sufficient to detect many types of errors.
        Options options = new()
        {
            BufferRadius = buffer_radius,
            ErrorFraction = error_fraction,
        };
        MutableS2ShapeIndex output =
        [
            DoBuffer(op => op.AddShapeIndex(input), options)
        ];

        _logger.WriteLine(@$"
radius = {buffer_radius.Radians:g17}, error_fraction = {error_fraction:g17}
input = {S2TextFormat.ToDebugString(input)}
output = {S2TextFormat.ToDebugString(output)}");

        // Check the 1a*/1b* condition above.
        S1Angle max_error = options.GetMaxError();
        TestContainment(input, output, buffer_radius, max_error);

        S1ChordAngle min_dist = new(S1Angle.Max(S1Angle.Zero, S1Angle.Abs(buffer_radius) - max_error));
        TestMinimumDistance(input, output, min_dist);

        // Check the 2a*/2b* condition (i.e., directed Hausdorff distance).
        S1ChordAngle max_dist = new(S1Angle.Abs(buffer_radius) + max_error);
        TestHausdorffDistance(input, output, max_dist);
    }

    // Convenience function that takes an S2ShapeIndex in s2textformat format.
    private void TestBuffer(string index_str, S1Angle buffer_radius,
                    double error_fraction) =>
        TestBuffer(MakeIndexOrDie(index_str), buffer_radius, error_fraction);

    // Convenience function that tests buffering using +/- the given radius.
    private void TestSignedBuffer(string index_str, S1Angle buffer_radius,
                    double error_fraction)
    {
        TestBuffer(index_str, buffer_radius, error_fraction);
        TestBuffer(index_str, -buffer_radius, error_fraction);
    }

    [Fact]
    internal void Test_S2BufferOperation_PointShell() =>
        TestSignedBuffer("# # 0:0", S1Angle.FromRadians(S2.M_PI_2), 0.01);

    [Fact]
    internal void Test_S2BufferOperation_SiblingPairShell() =>
        TestSignedBuffer("# # 0:0, 0:5", S1Angle.FromRadians(S2.M_PI_2), 0.01);

    [Fact]
    internal void Test_S2BufferOperation_SiblingPairHole() =>
        TestSignedBuffer("# # 0:0, 0:10, 7:7; 3:4, 3:6", S1Angle.FromDegrees(1), 0.01);

    [Fact]
    internal void Test_S2BufferOperation_Square()
    {
        TestSignedBuffer("# # -3:-3, -3:3, 3:3, 3:-3", S1Angle.FromDegrees(1), 0.01);
        TestSignedBuffer("# # -3:-3, -3:3, 3:3, 3:-3", S1Angle.FromDegrees(170), 1e-4);
    }

    [Fact]
    internal void Test_S2BufferOperation_HollowSquare() =>
        TestSignedBuffer("# # -3:-3, -3:3, 3:3, 3:-3; 2:2, -2:2, -2:-2, 2:-2",
            S1Angle.FromDegrees(1), 0.01);

    [Fact]
    internal void Test_S2BufferOperation_ZigZagLoop() =>
        TestSignedBuffer("# # 0:0, 0:7, 5:3, 5:10, 6:10, 6:1, 1:5, 1:0",
            S1Angle.FromDegrees(0.2), 0.01);

    [Fact]
    internal void Test_S2BufferOperation_Fractals()
    {
        foreach (double dimension in new[] { 1.02, 1.8 })
        {
            S2Testing.Fractal fractal = new();
            fractal.SetLevelForApproxMaxEdges(3 * 64);
            fractal.FractalDimension = dimension;
            var loop = fractal.MakeLoop(S2.GetFrame(new S2Point(1, 0, 0)),
                                            S1Angle.FromDegrees(10));
            MutableS2ShapeIndex input = [new S2Loop.Shape(loop)];
            TestBuffer(input, S1Angle.FromDegrees(0.4), 0.01);
        }
    }

    [Fact]
    internal void Test_S2BufferOperation_S2Curve()
    {
        // Tests buffering the S2 curve by an amount that yields the full polygon.
        const int kLevel = 2;  // Number of input edges == 6 * (4 ** kLevel)
        List<S2Point> points = [];
        for (S2CellId id = S2CellId.Begin(kLevel);
                id != S2CellId.End(kLevel); id = id.Next())
        {
            points.Add(id.ToPoint());
        }
        var pointsArr = points.ToArray();
        // Buffering by this amount or more is guaranteed to yield the full polygon.
        // (Note that the bound is not tight for S2CellIds at low levels.)
        S1Angle full_radius = S1Angle.FromRadians(0.5 * S2.kMaxDiag.GetValue(kLevel));
        Assert.True(DoBuffer(op =>
        {
            op.AddShape(new S2LaxClosedPolylineShape(pointsArr));
        }, full_radius, 0.1).IsFull());
    }

    // Tests buffering the given S2ShapeIndex with a variety of radii and error
    // fractions.  This method is intended to be used with relatively simple
    // shapes since calling it is quite expensive.
    private void TestRadiiAndErrorFractions(string index_str)
    {
        // Try the full range of radii with a representative error fraction.
        const double kFrac = 0.01;
        List<double> kTestRadiiRadians =
        [
            0,
            1e-300,
            1e-15,
            2e-15,
            3e-15,
            1e-5,
            0.01,
            0.1,
            1.0,
            (1 - kFrac) * S2.M_PI_2,
            S2.M_PI_2 - 1e-15,
            S2.M_PI_2,
            S2.M_PI_2 + 1e-15,
            (1 - kFrac) * S2.M_PI,
            S2.M_PI - 1e-6,
            S2.M_PI,
            1e300
        ];
        foreach (double radius in kTestRadiiRadians)
        {
            TestSignedBuffer(index_str, S1Angle.FromRadians(radius), kFrac);
        }

        // Now try the full range of error fractions with a few selected radii.
        List<double> kTestErrorFractions = [Options.kMinErrorFraction, 0.001, 0.01, 0.1, 1.0];
        foreach (double error_fraction in kTestErrorFractions)
        {
            TestBuffer(index_str, S1Angle.FromRadians(-1e-6), error_fraction);
            TestBuffer(index_str, S1Angle.FromRadians(1e-14), error_fraction);
            TestBuffer(index_str, S1Angle.FromRadians(1e-2), error_fraction);
            TestBuffer(index_str, S1Angle.FromRadians(S2.M_PI - 1e-3), error_fraction);
        }
    }

    [Fact]
    internal void Test_S2BufferOperation_RadiiAndErrorFractionCoverage()
    {
        // Test buffering simple shapes with a wide range of different buffer radii
        // and error fractions.

        // A single point.
        TestRadiiAndErrorFractions("1:1 # #");

        // A zig-zag polyline.
        TestRadiiAndErrorFractions("# 0:0, 0:30, 30:30, 30:60 #");

        // A triangular polygon with a triangular hole.  (The hole is clockwise.)
        TestRadiiAndErrorFractions("# # 0:0, 0:100, 70:50; 10:20, 50:50, 10:80");

        // A triangle with one very short and two very long edges.
        TestRadiiAndErrorFractions("# # 0:0, 0:179.99999999999, 1e-300:0");
    }

    internal class TestBufferPolyline
    {
        // Tests buffering a polyline with the given options.  This method is intended
        // only for testing Options.EndCapStyle and Options.PolylineSide; if these
        // options have their default values then TestBuffer() should be used
        // instead.  Similarly TestBuffer should be used to test negative buffer radii
        // and polylines with 0 or 1 vertices.
        internal TestBufferPolyline(string input_str,
                            Options options)
        {
            polyline_ = ParsePointsOrDie(input_str);
            buffer_radius_ = options.BufferRadius;
            max_error_ = options.GetMaxError();
            min_dist_ = new(S1Angle.Max(S1Angle.Zero, buffer_radius_ - max_error_));
            max_dist_ = new(buffer_radius_ + max_error_);
            round_ = options.EndCapStyle_ == EndCapStyle.ROUND;
            two_sided_ = options.PolylineSide_ == PolylineSide.BOTH;

            Assert.True(polyline_.Count >= 2);
            Assert.True(buffer_radius_ > S1Angle.Zero);

            MutableS2ShapeIndex input = [MakeLaxPolylineOrDie(input_str)];
            output_.Add(DoBuffer(
                op => op.AddShapeIndex(input), options));

            // Even with one-sided buffering and flat end caps the Hausdorff distance
            // criterion should still be true.  (This checks that the buffered result
            // is never further than (buffer_radius + max_error) from the input.)
            TestHausdorffDistance(input, output_, max_dist_);

            // However the minimum distance criterion is different; it only applies to
            // the portions of the boundary that are buffered using the given radius.
            // We check this approximately by walking along the polyline and checking
            // that (1) on portions of the polyline that should be buffered, the output
            // contains the offset point at distance (buffer_radius - max_error) and (2)
            // on portions of the polyline that should not be buffered, the output does
            // not contain the offset point at distance max_error.  The tricky part is
            // that both of these conditions have exceptions: (1) may not hold if the
            // test point is closest to the non-buffered side of the polyline (see the
            // last caveat in the documentation for Options.polyline_side), and (2)
            // may not hold if the test point is within (buffer_radius + max_error) of
            // the buffered side of any portion of the polyline.
            if (min_dist_ == S1ChordAngle.Zero) return;

            // Left-sided buffering is tested by reversing the polyline and then testing
            // whether it has been buffered correctly on the right.
            if (options.PolylineSide_ == PolylineSide.LEFT)
            {
                polyline_.Reverse();
            }

            int n = polyline_.Count;
            S2Point start0 = polyline_[0], start1 = polyline_[1];
            S2Point start_begin = GetEdgeAxis(start0, start1);
            S2Point start_mid = start0.CrossProd(start_begin);
            TestVertexArc(start0, start_begin, start_mid, round_ && two_sided_);
            TestVertexArc(start0, start_mid, -start_begin, round_);
            for (int i = 0; i < n - 2; ++i)
            {
                TestEdgeAndVertex(polyline_[i], polyline_[i + 1], polyline_[i + 2], true);
            }
            S2Point end0 = polyline_[n - 1], end1 = polyline_[n - 2];
            S2Point end_begin = GetEdgeAxis(end0, end1);
            S2Point end_mid = end0.CrossProd(end_begin);
            TestEdgeArc(end_begin, end1, end0, true);
            TestVertexArc(end0, end_begin, end_mid, round_);
            TestVertexArc(end0, end_mid, -end_begin, round_ && two_sided_);
            for (int i = n - 3; i >= 0; --i)
            {
                TestEdgeAndVertex(polyline_[i + 2], polyline_[i + 1], polyline_[i],
                                    two_sided_);
            }
            TestEdgeArc(start_begin, start1, start0, two_sided_);
        }

        private const double kArcLo = 0.001;
        private const double kArcHi = 0.999;
        private const int kArcSamples = 7;

        private static S2Point GetEdgeAxis(S2Point a, S2Point b)
        {
            return S2.RobustCrossProd(a, b).Normalize();
        }

        private bool PointBufferingUncertain(S2Point p, bool expect_contained)
        {
            // The only case where a point might be excluded from the buffered output is
            // if it is on the unbuffered side of the polyline.
            if (expect_contained && two_sided_) return false;

            int n = polyline_.Count;
            for (int i = 0; i < n - 1; ++i)
            {
                S2Point a = polyline_[i];
                S2Point b = polyline_[i + 1];
                if (!two_sided_)
                {
                    // Ignore points on the buffered side if expect_contained is true,
                    // and on the unbuffered side if expect_contained is false.
                    if ((S2Pred.Sign(a, b, p) < 0) == expect_contained) continue;
                }
                // TODO(ericv): Depending on how the erasing optimization is implemented,
                // it might be possible to add "&& expect_contained" to the test below.
                if (round_)
                {
                    if (S2.IsDistanceLess(p, a, b, max_dist_)) return true;
                }
                else
                {
                    if (S2.IsInteriorDistanceLess(p, a, b, max_dist_)) return true;
                    if (i > 0 && new S1ChordAngle(p, a) < max_dist_) return true;
                    if (i == n - 2 && new S1ChordAngle(p, b) < max_dist_) return true;
                }
            }
            return false;
        }

        private void TestPoint(S2Point p, S2Point dir, bool expect_contained)
        {
            S2Point x = S2.GetPointOnRay(
                p, dir, expect_contained ? buffer_radius_ - max_error_ : max_error_);
            if (!PointBufferingUncertain(x, expect_contained))
            {
                Assert.Equal(S2ContainsPointQueryFactory.MakeS2ContainsPointQuery(output_).Contains(x),
                            expect_contained);
            }
        }

        private void TestVertexArc(S2Point p, S2Point start, S2Point end,
                        bool expect_contained)
        {
            for (double t = kArcLo; t < 1; t += (kArcHi - kArcLo) / kArcSamples)
            {
                S2Point dir = S2.Interpolate(start, end, t);
                TestPoint(p, dir, expect_contained);
            }
        }

        private void TestEdgeArc(S2Point ba_axis, S2Point a, S2Point b,
                        bool expect_contained)
        {
            for (double t = kArcLo; t < 1; t += (kArcHi - kArcLo) / kArcSamples)
            {
                S2Point p = S2.Interpolate(a, b, t);
                TestPoint(p, ba_axis, expect_contained);
            }
        }

        private void TestEdgeAndVertex(S2Point a, S2Point b, S2Point c,
                            bool expect_contained)
        {
            S2Point ba_axis = GetEdgeAxis(b, a);
            S2Point cb_axis = GetEdgeAxis(c, b);
            TestEdgeArc(ba_axis, a, b, expect_contained);
            TestVertexArc(b, ba_axis, cb_axis, expect_contained);
        }

        private readonly List<S2Point> polyline_;
        private readonly MutableS2ShapeIndex output_ = [];
        private readonly S1Angle buffer_radius_, max_error_;
        private readonly S1ChordAngle min_dist_, max_dist_;
        private readonly bool round_, two_sided_;
    }

    [Fact]
    internal void Test_S2BufferOperation_ZigZagPolyline()
    {
        Options options = new(S1Angle.FromDegrees(1));
        foreach (var polyline_side in new[] { PolylineSide.LEFT, PolylineSide.RIGHT, PolylineSide.BOTH })
        {
            foreach (var end_cap_style in new[] { EndCapStyle.ROUND, EndCapStyle.FLAT })
            {
                _logger.WriteLine(
                    $"two_sided = {polyline_side == PolylineSide.BOTH}, round = {end_cap_style == EndCapStyle.ROUND}");
                options.PolylineSide_ = polyline_side;
                options.EndCapStyle_ = end_cap_style;
                _ = new TestBufferPolyline("0:0, 0:7, 5:3, 5:10", options);  // NOLINT
                _ = new TestBufferPolyline("10:0, 0:0, 5:1", options);       // NOLINT
            }
        }
    }
}
