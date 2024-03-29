﻿namespace S2Geometry;

using static S2WindingOperation;

public class S2WindingOperationTests(ITestOutputHelper logger)
{
    private readonly ITestOutputHelper _logger = logger;

    // Verifies that the S2WindingOperation with the given arguments produces the
    // given result.  In order to ensure that results are not affected by the
    // cyclic order of the loop vertices or the S2ShapeIndex loop ordering,
    // we compute the symmetric difference using S2BooleanOperation.  (We don't
    // use s2builderutil::IndexMatchingLayer because that class does not
    // distinguish empty from full polygons, and we don't need its ability to
    // match edge multiplicities here.)
    private void ExpectWindingResult(
        Options options, List<string> loop_strs,
        string ref_point_str, int ref_winding,
        WindingRule rule, string expected_str)
    {
        MutableS2ShapeIndex expected = [MakeLaxPolygonOrDie(expected_str)];
        MutableS2ShapeIndex actual = [];
        S2WindingOperation winding_op = new(
            new IndexedLaxPolygonLayer(actual), options);
        foreach (var loop_str in loop_strs)
        {
            winding_op.AddLoop(ParsePointsOrDie(loop_str));
        }
        Assert.True(winding_op.Build(MakePointOrDie(ref_point_str),
            ref_winding, rule, out _));
        _logger.WriteLine("(INFO) Actual: " + S2TextFormat.ToDebugString(actual));
        S2LaxPolygonShape difference = new();
        S2BooleanOperation diff_op = new(
            S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE,
            new LaxPolygonLayer(difference));
        Assert.True(diff_op.Build(actual, expected, out _));
        Assert.True(difference.IsEmpty());
    }

    // Like ExpectWindingResult(), but with two different expected results
    // depending on whether options.include_degeneracies() is false or true.
    private void ExpectDegenerateWindingResult(
                    Options options, List<string> loop_strs,
                    string ref_point_str, int ref_winding,
                    WindingRule rule,
                    string expected_str_false, string expected_str_true)
    {
        options.IncludeDegeneracies = false;
        ExpectWindingResult(options, loop_strs, ref_point_str, ref_winding, rule,
            expected_str_false);
        options.IncludeDegeneracies = true;
        ExpectWindingResult(options, loop_strs, ref_point_str, ref_winding, rule,
            expected_str_true);
    }

    [Fact]
    internal void Test_S2WindingOperation_Empty()
    {
        ExpectWindingResult(
            new(),
            [""], "5:5", 0, WindingRule.POSITIVE, "");
        ExpectWindingResult(
            new(),
            [""], "5:5", 1, WindingRule.POSITIVE, "full");
    }

    [Fact]
    internal void Test_S2WindingOperation_PointLoop()
    {
        ExpectDegenerateWindingResult(
            new(),
            ["2:2"], "5:5", 0, WindingRule.POSITIVE,
            "", "2:2");
    }

    [Fact]
    internal void Test_S2WindingOperation_S2WindingOperation_SiblingPairLoop()
    {
        ExpectDegenerateWindingResult(
            new(),
            ["2:2, 3:3"], "5:5", 0, WindingRule.POSITIVE,
            "", "2:2, 3:3");
    }

    [Fact]
    internal void Test_S2WindingOperation_Rectangle()
    {
        ExpectWindingResult(
            new(),
            ["0:0, 0:10, 10:10, 10:0"], "5:5", 1, WindingRule.POSITIVE,
            "0:0, 0:10, 10:10, 10:0");
        ExpectWindingResult(
            new(),
            ["0:0, 0:10, 10:10, 10:0"], "5:5", 1, WindingRule.NEGATIVE,
            "");
        ExpectWindingResult(
            new(),
            ["0:0, 0:10, 10:10, 10:0"], "5:5", 1, WindingRule.NON_ZERO,
            "0:0, 0:10, 10:10, 10:0");
        ExpectWindingResult(
            new(),
            ["0:0, 0:10, 10:10, 10:0"], "5:5", 1, WindingRule.ODD,
            "0:0, 0:10, 10:10, 10:0");
    }

    [Fact]
    internal void Test_S2WindingOperation_BowTie()
    {
        // Note that NEGATIVE, NON_ZERO, and ODD effectively reverse the orientation
        // of one of the two triangles that form the bow tie.
        ExpectWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(1))),
            ["5:-5, -5:5, 5:5, -5:-5"], "10:0", 0, WindingRule.POSITIVE,
            "0:0, -5:5, 5:5");
        ExpectWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(1))),
            ["5:-5, -5:5, 5:5, -5:-5"], "10:0", 0, WindingRule.NEGATIVE,
            "-5:-5, 0:0, 5:-5");
        ExpectWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(1))),
            ["5:-5, -5:5, 5:5, -5:-5"], "10:0", 0, WindingRule.NON_ZERO,
            "0:0, -5:5, 5:5; -5:-5, 0:0, 5:-5");
        ExpectWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(1))),
            ["5:-5, -5:5, 5:5, -5:-5"], "10:0", 0, WindingRule.ODD,
            "0:0, -5:5, 5:5; -5:-5, 0:0, 5:-5");
    }

    [Fact]
    internal void Test_S2WindingOperation_CollapsingShell()
    {
        ExpectDegenerateWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(5))),
            ["0:0, 0:3, 3:3"], "10:0", 0, WindingRule.POSITIVE,
            "", "0:0");
        ExpectDegenerateWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(5))),
            ["0:0, 0:3, 3:3"], "1:1", 1, WindingRule.POSITIVE,
            "", "0:0");
        ExpectWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(5))),
            ["0:0, 3:3, 0:3"], "10:0", 1, WindingRule.POSITIVE,
            "full");
        ExpectWindingResult(
            new(new IdentitySnapFunction(S1Angle.FromDegrees(5))),
            ["0:0, 3:3, 0:3"], "1:1", 0, WindingRule.POSITIVE,
            "full");
    }

    // Two triangles that touch along a common boundary.
    [Fact]
    internal void Test_S2WindingOperation_TouchingTriangles()
    {
        // The touch edges are considered to form a degenerate hole.  Such holes are
        // removed by WindingRule::POSITIVE since they are not needed for computing
        // N-way unions.  They are kept by WindingRule::ODD since they are needed in
        // order to compute N-way symmetric differences.
        ExpectWindingResult(
            new(),
            ["0:0, 0:8, 8:8", "0:0, 8:8, 8:0"], "1:1", 1, WindingRule.POSITIVE,
            "0:0, 0:8, 8:8, 8:0");
        ExpectDegenerateWindingResult(
            new(),
            ["0:0, 0:8, 8:8", "0:0, 8:8, 8:0"], "2:2", 1, WindingRule.ODD,
            "0:0, 0:8, 8:8, 8:0", "0:0, 0:8, 8:8; 0:0, 8:8, 8:0");
    }

    // Like the test above, but the triangles only touch after snapping.
    [Fact]
    internal void Test_S2WindingOperation_TouchingTrianglesAfterSnapping()
    {
        // The snap function below rounds coordinates to the nearest degree.
        ExpectWindingResult(
            new(new IntLatLngSnapFunction(0)),
            [
                "0.1:0.2, 0:7.8, 7.6:8.2",
                "0.3:0.2, 8.1:7.8, 7.6:0.4"
            ],
            "6:2", 1, WindingRule.POSITIVE,
            "0:0, 0:8, 8:8, 8:0");
        ExpectDegenerateWindingResult(
            new(new IntLatLngSnapFunction(0)),
            [
                "0.1:0.2, 0:7.8, 7.6:8.2",
                "0.3:0.2, 8.1:7.8, 7.6:0.4"
            ],
            "2:6", 1, WindingRule.ODD,
            "0:0, 0:8, 8:8, 8:0", "0:0, 0:8, 8:8; 0:0, 8:8, 8:0");
    }

    // This tests an N-way union of 5 overlapping squares forming a "staircase".
    [Fact]
    internal void Test_S2WindingOperation_UnionOfSquares()
    {
        ExpectWindingResult(
            new(new IntLatLngSnapFunction(1)),
            [
                "0:0, 0:4, 4:4, 4:0",
                "1:1, 1:5, 5:5, 5:1",
                "2:2, 2:6, 6:6, 6:2",
                "3:3, 3:7, 7:7, 7:3",
                "4:4, 4:8, 8:8, 8:4"
            ],
            "0.5:0.5", 1, WindingRule.POSITIVE,
            "7:4, 7:3, 6:3, 6:2, 5:2, 5:1, 4:1, 4:0, 0:0, 0:4, " +
            "1:4, 1:5, 2:5, 2:6, 3:6, 3:7, 4:7, 4:8, 8:8, 8:4");

        // This computes the region overlapped by at least two squares.
        ExpectWindingResult(
            new(new IntLatLngSnapFunction(1)),
            [
                "0:0, 0:4, 4:4, 4:0",
                "1:1, 1:5, 5:5, 5:1",
                "2:2, 2:6, 6:6, 6:2",
                "3:3, 3:7, 7:7, 7:3",
                "4:4, 4:8, 8:8, 8:4"
            ],
            "0.5:0.5", 0, WindingRule.POSITIVE,
            "6:4, 6:3, 5:3, 5:2, 4:2, 4:1, 1:1, 1:4, 2:4, 2:5, " +
            "3:5, 3:6, 4:6, 4:7, 7:7, 7:4");

        // This computes the region overlapped by at least three squares.
        ExpectWindingResult(
            new(new IntLatLngSnapFunction(1)),
            [
                "0:0, 0:4, 4:4, 4:0",
                "1:1, 1:5, 5:5, 5:1",
                "2:2, 2:6, 6:6, 6:2",
                "3:3, 3:7, 7:7, 7:3",
                "4:4, 4:8, 8:8, 8:4"
            ],
            "0.5:0.5", -1, WindingRule.POSITIVE,
            "5:4, 5:3, 4:3, 4:2, 2:2, 2:4, 3:4, 3:5, 4:5, 4:6, 6:6, 6:4");

        // This computes the region overlapped by at least four squares.
        ExpectWindingResult(
            new(new IntLatLngSnapFunction(1)),
            [
                "0:0, 0:4, 4:4, 4:0",
                "1:1, 1:5, 5:5, 5:1",
                "2:2, 2:6, 6:6, 6:2",
                "3:3, 3:7, 7:7, 7:3",
                "4:4, 4:8, 8:8, 8:4"
            ],
            "0.5:0.5", -2, WindingRule.POSITIVE,
            "3:3, 3:4, 4:4, 4:3; 4:4, 4:5, 5:5, 5:4");

        // WindingRule::ODD yields a pattern reminiscent of a checkerboard.
        ExpectWindingResult(
            new(new IntLatLngSnapFunction(1)),
            [
                "0:0, 0:4, 4:4, 4:0",
                "1:1, 1:5, 5:5, 5:1",
                "2:2, 2:6, 6:6, 6:2",
                "3:3, 3:7, 7:7, 7:3",
                "4:4, 4:8, 8:8, 8:4"
            ],
            "0.5:0.5", 1, WindingRule.ODD,
            "4:1, 4:0, 0:0, 0:4, 1:4, 1:1; " +
            "4:3, 4:2, 2:2, 2:4, 3:4, 3:3; " +
            "1:4, 1:5, 2:5, 2:4; " +
            "5:4, 5:3, 4:3, 4:4; " +
            "5:2, 5:1, 4:1, 4:2; " +
            "2:5, 2:6, 3:6, 3:5; " +
            "6:3, 6:2, 5:2, 5:3; " +
            "3:6, 3:7, 4:7, 4:6; " +
            "3:4, 3:5, 4:5, 4:4; " +
            "7:4, 7:3, 6:3, 6:4; " +
            "4:7, 4:8, 8:8, 8:4, 7:4, 7:7; " +
            "4:5, 4:6, 6:6, 6:4, 5:4, 5:5");
    }

    // This tests that WindingRule::ODD can be used to compute the symmtric
    // difference even for input geometries with degeneracies, e.g. one geometry
    // has a degenerate hole or degenerate shell that the other does not.
    [Fact]
    internal void Test_S2WindingOperation_SymmetricDifferenceDegeneracies()
    {
        ExpectDegenerateWindingResult(
            new(new IntLatLngSnapFunction(1)),
            [
                "0:0, 0:3, 3:3, 3:0",
                "1:1",
                "2:2",
                "4:4",   // Geometry 1
                "0:0, 0:3, 3:3, 3:0",
                "1:1",
                "4:4",
                "5:5"
            ],  // Geometry 2
            "10:10", 0, WindingRule.ODD,
            "", "2:2; 5:5");
    }

    [Fact]
    internal void Test_S2WindingOperationOptions_SetGetSnapFunction()
    {
        // Prevent these from being detected as dead code.
        S2WindingOperation.Options opts = new()
        {
            SnapFunction = new IdentitySnapFunction()
        };
        Assert.Equal(opts.SnapFunction.SnapRadius, S1Angle.Zero);
    }
}
