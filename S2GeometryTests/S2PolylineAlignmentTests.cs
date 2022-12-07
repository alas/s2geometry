namespace S2Geometry;

using WarpPath = List<(int start, int end)>;

public class S2PolylineAlignmentTests
{
    // PRIVATE API TESTS

    [Fact]
    internal void Test_S2PolylineAlignmentTest_CreatesWindowFromStrides() {
        //    0 1 2 3 4 5
        //  0 * * * . . .
        //  1 . * * * . .
        //  2 . . * * . .
        //  3 . . . * * *
        //  4 . . . . * *
        var strides = new [] { (0, 3), (1, 4), (2, 4), (3, 6), (4, 6) };
        var w = new S2PolylineAlignment.Window(strides);
        Assert.Equal(0, w.GetColumnStride(0).start);
        Assert.Equal(3, w.GetColumnStride(0).end);
        Assert.Equal(4, w.GetColumnStride(4).start);
        Assert.Equal(6, w.GetColumnStride(4).end);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_CreatesWindowFromWarpPath() {
        //   0 1 2 3 4 5
        // 0 * . . . . .
        // 1 * * . . . .
        // 2 . * . . . .
        // 3 . * * * . .
        // 4 . . . . * *
        var path = new WarpPath
            {
                (0, 0), (1, 0), (1, 1), (2, 1), (3, 1),
                (3, 2), (3, 3), (4, 4), (4, 5)
            };
        var w = new S2PolylineAlignment.Window(path);
        Assert.Equal(0, w.GetColumnStride(0).start);
        Assert.Equal(1, w.GetColumnStride(0).end);
        Assert.Equal(0, w.GetColumnStride(1).start);
        Assert.Equal(2, w.GetColumnStride(1).end);
        Assert.Equal(1, w.GetColumnStride(2).start);
        Assert.Equal(2, w.GetColumnStride(2).end);
        Assert.Equal(1, w.GetColumnStride(3).start);
        Assert.Equal(4, w.GetColumnStride(3).end);
        Assert.Equal(4, w.GetColumnStride(4).start);
        Assert.Equal(6, w.GetColumnStride(4).end);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_GeneratesWindowDebugString() {
        var strides = new[]{ (0, 4), (0, 4), (0, 4), (0, 4) };
        var w = new S2PolylineAlignment.Window(strides);
        string expected_output = @"
 * * * *
 * * * *
 * * * *
 * * * *
";
Assert.Equal(("\r\n" + w).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_UpsamplesWindowByFactorOfTwo() {
        //   0 1 2 3 4 5
        // 0 * * * . . .
        // 1 . * * * . .
        // 2 . . * * . .
        // 3 . . . * * *
        // 4 . . . . * *
        var strides = new []{
  (0, 3), (1, 4), (2, 4), (3, 6), (4, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_upscaled = w.Upsample(10, 12);
        string expected_output = @"
 * * * * * * . . . . . .
 * * * * * * . . . . . .
 . . * * * * * * . . . .
 . . * * * * * * . . . .
 . . . . * * * * . . . .
 . . . . * * * * . . . .
 . . . . . . * * * * * *
 . . . . . . * * * * * *
 . . . . . . . . * * * *
 . . . . . . . . * * * *
";
        Assert.Equal(("\r\n" + w_upscaled).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_UpsamplesWindowXAxisByFactorOfThree() {
        //   0 1 2 3 4 5
        // 0 * * * . . .
        // 1 . * * * . .
        // 2 . . * * . .  3 . . . * * *
        // 4 . . . . * *
        var strides = new[]{
  (0, 3), (1, 4), (2, 4), (3, 6), (4, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_upscaled = w.Upsample(5, 18);
        string expected_output = @"
 * * * * * * * * * . . . . . . . . .
 . . . * * * * * * * * * . . . . . .
 . . . . . . * * * * * * . . . . . .
 . . . . . . . . . * * * * * * * * *
 . . . . . . . . . . . . * * * * * *
";
Assert.Equal(("\r\n" + w_upscaled).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_UpsamplesWindowYAxisByFactorOfThree() {
        //   0 1 2 3 4 5
        // 0 * * * . . .
        // 1 . * * * . .
        // 2 . . * * . .
        // 3 . . . * * *
        // 4 . . . . * *
        var strides = new []{
  (0, 3), (1, 4), (2, 4), (3, 6), (4, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_upscaled = w.Upsample(15, 6);
        string expected_output = @"
 * * * . . .
 * * * . . .
 * * * . . .
 . * * * . .
 . * * * . .
 . * * * . .
 . . * * . .
 . . * * . .
 . . * * . .
 . . . * * *
 . . . * * *
 . . . * * *
 . . . . * *
 . . . . * *
 . . . . * *
";
Assert.Equal(("\r\n" + w_upscaled).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_UpsamplesWindowByNonInteger() {
        //   0 1 2 3 4 5
        // 0 * * * . . .
        // 1 . * * * . .
        // 2 . . * * . .
        // 3 . . . * * *
        // 4 . . . . * *
        var strides = new []{
  (0, 3), (1, 4), (2, 4), (3, 6), (4, 6)};
        var w = new S2PolylineAlignment.Window(strides);

        var w_upscaled = w.Upsample(19, 23);
        string expected_output = @"
 * * * * * * * * * * * * . . . . . . . . . . .
 * * * * * * * * * * * * . . . . . . . . . . .
 * * * * * * * * * * * * . . . . . . . . . . .
 * * * * * * * * * * * * . . . . . . . . . . .
 . . . . * * * * * * * * * * * . . . . . . . .
 . . . . * * * * * * * * * * * . . . . . . . .
 . . . . * * * * * * * * * * * . . . . . . . .
 . . . . * * * * * * * * * * * . . . . . . . .
 . . . . . . . . * * * * * * * . . . . . . . .
 . . . . . . . . * * * * * * * . . . . . . . .
 . . . . . . . . * * * * * * * . . . . . . . .
 . . . . . . . . . . . . * * * * * * * * * * *
 . . . . . . . . . . . . * * * * * * * * * * *
 . . . . . . . . . . . . * * * * * * * * * * *
 . . . . . . . . . . . . * * * * * * * * * * *
 . . . . . . . . . . . . . . . * * * * * * * *
 . . . . . . . . . . . . . . . * * * * * * * *
 . . . . . . . . . . . . . . . * * * * * * * *
 . . . . . . . . . . . . . . . * * * * * * * *
";
Assert.Equal(("\r\n" + w_upscaled).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_DilatesWindowByRadiusZero() {
        //   0 1 2 3 4 5
        // 0 * * * . . .
        // 1 . . * . . .
        // 2 . . * . . .
        // 3 . . * * . .
        // 4 . . . * * *
        var strides = new[]
        {(0, 3), (2, 3), (2, 3), (2, 4), (3, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_d = w.Dilate(0);
        string expected_output = @"
 * * * . . .
 . . * . . .
 . . * . . .
 . . * * . .
 . . . * * *
";
Assert.Equal(("\r\n" + w_d).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_DilatesWindowByRadiusOne() {
        //   0 1 2 3 4 5 (x's are the spots that we dilate into)
        // 0 * * * x . .
        // 1 x x * x . .
        // 2 . x * x x .
        // 3 . x * * x x
        // 4 . x x * * *
        var strides = new []
        {(0, 3), (2, 3), (2, 3), (2, 4), (3, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_d = w.Dilate(1);
        string expected_output = @"
 * * * * . .
 * * * * . .
 . * * * * .
 . * * * * *
 . * * * * *
";
Assert.Equal(("\r\n" + w_d).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_DilatesWindowByRadiusTwo() {
        //   0 1 2 3 4 5 (x's are the spots that we dilate into)
        // 0 * * * x x .
        // 1 x x * x x x
        // 2 x x * x x x
        // 3 x x * * x x
        // 4 x x x * * *
        var strides = new []
        {(0, 3), (2, 3), (2, 3), (2, 4), (3, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_d = w.Dilate(2);
        string expected_output = @"
 * * * * * .
 * * * * * *
 * * * * * *
 * * * * * *
 * * * * * *
";
Assert.Equal(("\r\n" + w_d).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_DilatesWindowByVeryLargeRadius() {
        var strides = new []{
            (0, 3), (2, 3), (2, 3), (2, 4), (3, 6)};
        var w = new S2PolylineAlignment.Window(strides);
        var w_d = w.Dilate(100);
        string expected_output = @"
 * * * * * *
 * * * * * *
 * * * * * *
 * * * * * *
 * * * * * *
";
Assert.Equal(("\r\n" + w_d).NoCR(), expected_output.NoCR());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_HalvesZeroLengthPolyline() {
        var line = MakePolylineOrDie("");
        var halved = S2PolylineAlignment.HalfResolution(line);
        var correct = MakePolylineOrDie("");
        Assert.Equal(halved.ToDebugString(), correct.ToDebugString());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_HalvesEvenLengthPolyline() {
        var line = MakePolylineOrDie("0:0, 0:1, 0:2, 1:2");
        var halved = S2PolylineAlignment.HalfResolution(line);
        var correct = MakePolylineOrDie("0:0, 0:2");
        Assert.Equal(halved.ToDebugString(), correct.ToDebugString());
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_HalvesOddLengthPolyline() {
        var line = MakePolylineOrDie("0:0, 0:1, 0:2, 1:2, 3:5");
        var halved = S2PolylineAlignment.HalfResolution(line);
        var correct = MakePolylineOrDie("0:0, 0:2, 3:5");
        Assert.Equal(halved.ToDebugString(), correct.ToDebugString());
    }

    // internal API TESTS

    [Fact]
    internal void Test_S2PolylineAlignmentDeathTest_ExactLengthZeroInputs() {
        var a = MakePolylineOrDie("");
        var b = MakePolylineOrDie("");
        var correct_path = new WarpPath{ };
        Assert.Throws<Exception>(() => VerifyPath(a, b, correct_path));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentDeathTest_ExactLengthZeroInputA() {
        var a = MakePolylineOrDie("");
        var b = MakePolylineOrDie("0:0, 1:1, 2:2");
        WarpPath correct_path = new();
        Assert.Throws<Exception>(() => VerifyPath(a, b, correct_path));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentDeathTest_ExactLengthZeroInputB() {
        var a = MakePolylineOrDie("0:0, 1:1, 2:2");
        var b = MakePolylineOrDie("");
        WarpPath correct_path = new();
        Assert.Throws<Exception>(() => VerifyPath(a, b, correct_path));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ExactLengthOneInputs() {
        var a = MakePolylineOrDie("1:1");
        var b = MakePolylineOrDie("2:2");
        var correct_path = new WarpPath{ (0, 0) };
        VerifyPath(a, b, correct_path);
        VerifyCost(a, b);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ExactLengthOneInputA() {
        var a = MakePolylineOrDie("0:0");
        var b = MakePolylineOrDie("0:0, 1:1, 2:2");
        var correct_path = new WarpPath{ (0, 0), (0, 1), (0, 2) };
        VerifyPath(a, b, correct_path);
        VerifyCost(a, b);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ExactLengthOneInputB() {
        var a = MakePolylineOrDie("0:0, 1:1, 2:2");
        var b = MakePolylineOrDie("0:0");
        var correct_path = new WarpPath{ (0, 0), (1, 0), (2, 0) };
        VerifyPath(a, b, correct_path);
        VerifyCost(a, b);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ExactHeaderFileExample() {
        var a = MakePolylineOrDie("1:0, 5:0, 6:0, 9:0");
        var b = MakePolylineOrDie("2:0, 7:0, 8:0");
        var correct_path = new WarpPath{ (0, 0), (1, 1), (2, 1), (3, 2) };
        VerifyPath(a, b, correct_path);
        VerifyCost(a, b);
    }

    // Take a small random selection of short correlated polylines and ensure that
    // the cost from the brute force solver equals the cost from the DP solvers.
    [Fact]
    internal void Test_S2PolylineAlignmentTest_FuzzedWithBruteForce() {
        int kNumPolylines = 10;
        int kNumVertices = 8;
        double kPerturbation = 1.5;
        var lines = GenPolylines(kNumPolylines, kNumVertices, kPerturbation);
        for (int i = 0; i < kNumPolylines; ++i) {
            for (int j = i + 1; j < kNumPolylines; ++j) {
                VerifyCost(lines[i], lines[j]);
            }
        }
    }

    // TESTS FOR TRAJECTORY CONSENSUS ALGORITHMS

    // Tests for GetMedoidPolyline
    [Fact]
    internal void Test_S2PolylineAlignmentDeathTest_MedoidPolylineNoPolylines() {
        S2Polyline[] polylines = Array.Empty<S2Polyline>();
        var default_opts = new S2PolylineAlignment.MedoidOptions();
        Assert.Throws<Exception>(() => S2PolylineAlignment.GetMedoidPolyline(polylines, default_opts));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_MedoidPolylineOnePolyline() {
        S2Polyline[] polylines = new S2Polyline[]{
        MakePolylineOrDie("5:0, 5:1, 5:2") };
        var default_opts = new S2PolylineAlignment.MedoidOptions();
        var medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, default_opts);
        Assert.Equal(0, medoid);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_MedoidPolylineTwoPolylines() {
        // Tie-breaking is contractually done by choosing the smallest tied index.
        // These inputs (really, any collection of two polylines) yield a tie.
        var polylines = new S2Polyline[]{
        MakePolylineOrDie("5:0, 5:1, 5:2"),
        MakePolylineOrDie("1:0, 1:1, 1:2")};

        var default_opts = new S2PolylineAlignment.MedoidOptions();
        var medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, default_opts);
        Assert.Equal(0, medoid);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_MedoidPolylineFewSmallPolylines() {
        var polylines = new S2Polyline[] {
            MakePolylineOrDie("5:0, 5:1, 5:2"),
            MakePolylineOrDie("3:0, 3:1, 3:2"),
            MakePolylineOrDie("1:0, 1:1, 1:2"),
        };

        var default_opts = new S2PolylineAlignment.MedoidOptions();
        var medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, default_opts);
        Assert.Equal(1, medoid);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_MedoidPolylineOverlappingPolylines() {
        // Given two identical polylines as input, break the tie with smallest index.
        var polylines = new S2Polyline[] {
        MakePolylineOrDie("1:0, 1:1, 1:2"),
        MakePolylineOrDie("1:0, 1:1, 1:2"),
        };
        var default_opts = new S2PolylineAlignment.MedoidOptions();
        var medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, default_opts);
        Assert.Equal(0, medoid);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_MedoidPolylineDifferentLengthPolylines() {
        var polylines = new S2Polyline[]
        {
            MakePolylineOrDie("5:0, 5:1, 5:2"),
            MakePolylineOrDie("3:0, 3:0.5, 3:1, 3:2"),
            MakePolylineOrDie("1:0, 1:0.5, 1:1, 1:1.5, 1:2")
        };

        var default_opts = new S2PolylineAlignment.MedoidOptions();
        var medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, default_opts);
        Assert.Equal(1, medoid);
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_MedoidPolylineFewLargePolylines() {
        // We pick num_vertices to be large so that the approx and exact vertex
        // alignment computations are likely to give different results.
        int num_polylines = 3;
        int num_vertices = 1024;
        double perturb = 0.9;
        var polylines = GenPolylines(num_polylines, num_vertices, perturb);

        // clang-format off
        double[] exact_costs = {
  S2PolylineAlignment.GetExactVertexAlignmentCost(polylines[0], polylines[1]) +
  S2PolylineAlignment.GetExactVertexAlignmentCost(polylines[0], polylines[2]),
  S2PolylineAlignment.GetExactVertexAlignmentCost(polylines[1], polylines[0]) +
  S2PolylineAlignment.GetExactVertexAlignmentCost(polylines[1], polylines[2]),
  S2PolylineAlignment.GetExactVertexAlignmentCost(polylines[2], polylines[0]) +
  S2PolylineAlignment.GetExactVertexAlignmentCost(polylines[2], polylines[1])
};
        double[] approx_costs = {
  S2PolylineAlignment.GetApproxVertexAlignment(polylines[0], polylines[1]).AlignmentCost +
  S2PolylineAlignment.GetApproxVertexAlignment(polylines[0], polylines[2]).AlignmentCost,
  S2PolylineAlignment.GetApproxVertexAlignment(polylines[1], polylines[0]).AlignmentCost +
  S2PolylineAlignment.GetApproxVertexAlignment(polylines[1], polylines[2]).AlignmentCost,
  S2PolylineAlignment.GetApproxVertexAlignment(polylines[2], polylines[0]).AlignmentCost +
  S2PolylineAlignment.GetApproxVertexAlignment(polylines[2], polylines[1]).AlignmentCost
};
        // clang-format on

        int exact_medoid_index = exact_costs.IndexOfMin();

        int approx_medoid_index = approx_costs.IndexOfMin();

        var options = new S2PolylineAlignment.MedoidOptions
        {
            Approx = false
        };
        var exact_medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, options);
        Assert.Equal(exact_medoid, exact_medoid_index);

        options.Approx = (true);
        var approx_medoid = S2PolylineAlignment.GetMedoidPolyline(polylines, options);
        Assert.Equal(approx_medoid, approx_medoid_index);
    }

    // Tests for GetConsensusPolyline

    [Fact]
    internal void Test_S2PolylineAlignmentDeathTest_ConsensusPolylineNoPolylines() {
        var polylines = Array.Empty<S2Polyline>();
        var default_opts = new S2PolylineAlignment.ConsensusOptions();
        Assert.Throws<Exception>(() => S2PolylineAlignment.GetConsensusPolyline(polylines, default_opts));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ConsensusPolylineOnePolyline() {
        var polylines = new[] { MakePolylineOrDie("3:0, 3:1, 3:2") };

        var default_opts = new S2PolylineAlignment.ConsensusOptions();
        var result = S2PolylineAlignment.GetConsensusPolyline(polylines, default_opts);
        var expected = MakePolylineOrDie("3:0, 3:1, 3:2");
        Assert.True(result.ApproxEquals(expected));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ConsensusPolylineTwoPolylines() {
        S2Polyline[] polylines = new S2Polyline[] 
        {
            MakePolylineOrDie("3:0, 3:1, 3:2"),
            MakePolylineOrDie("1:0, 1:1, 1:2"),
        };

        var default_opts = new S2PolylineAlignment.ConsensusOptions();
        var result = S2PolylineAlignment.GetConsensusPolyline(polylines, default_opts);
        var expected = MakePolylineOrDie("2:0, 2:1, 2:2");
        Assert.True(result.ApproxEquals(expected));
    }

    [Fact]
    internal void Test_S2PolylineAlignmentTest_ConsensusPolylineOverlappingPolylines() {
        var polylines = new S2Polyline[]
        {
            MakePolylineOrDie("1:0, 1:1, 1:2"),
            MakePolylineOrDie("1:0, 1:1, 1:2")
        };

        var default_opts = new S2PolylineAlignment.ConsensusOptions();
        var result = S2PolylineAlignment.GetConsensusPolyline(polylines, default_opts);
        var expected = MakePolylineOrDie("1:0, 1:1, 1:2");
        Assert.True(result.ApproxEquals(expected));
    }

    private static double[][] DistanceMatrix(S2Polyline a, S2Polyline b)
    {
        int a_n = a.NumVertices();
        int b_n = b.NumVertices();
        var table = new double[a_n][];
        for (int i = 0; i < a_n; ++i) {
            table[i] = new double[b_n];
            for (int j = 0; j < b_n; ++j) {
                table[i][j] = (a.Vertex(i) - b.Vertex(j)).Norm2();
            }
        }
        return table;
    }

    // Do some testing against random sequences with a brute-force solver.
    // Returns the optimal cost of alignment up until vertex i, j.
    private double GetBruteForceCost(double[][] table, int i, int j) {
        if (i == 0 && j == 0) {
            return table[0][0];
        } else if (i == 0) {
            return GetBruteForceCost(table, i, j - 1) + table[i][j];
        } else if (j == 0) {
            return GetBruteForceCost(table, i - 1, j) + table[i][j];
        } else {
            return new[]
                {
                    GetBruteForceCost(table, i - 1, j - 1),
                    GetBruteForceCost(table, i - 1, j),
                    GetBruteForceCost(table, i, j - 1)
                }.Min() + table[i][j];
        }
    }

    // Use Brute Force solver to verify exact Dynamic Programming solvers.
    private void VerifyCost(S2Polyline a, S2Polyline b) {
        int a_n = a.NumVertices();
        int b_n = b.NumVertices();
        double brute_cost =
            GetBruteForceCost(DistanceMatrix(a, b), a_n - 1, b_n - 1);
        double exact_cost = S2PolylineAlignment.GetExactVertexAlignmentCost(a, b);
        var exact_alignment = S2PolylineAlignment.GetExactVertexAlignment(a, b);
        Assert2.Near(brute_cost, exact_cost);
        Assert2.Near(brute_cost, exact_alignment.AlignmentCost);
    }

    // Check that the costs are the same between both exact computation methods, and
    // that the warp path matches the one given.
    private static void VerifyPath(S2Polyline a, S2Polyline b, WarpPath p) {
        double correct = 0;
        foreach (var (start, end) in p) {
            correct += (a.Vertex(start) - b.Vertex(end)).Norm2();
        }
        double exact_cost = S2PolylineAlignment.GetExactVertexAlignmentCost(a, b);
        var exact_alignment = S2PolylineAlignment.GetExactVertexAlignment(a, b);
        Assert2.Near(correct, exact_cost);
        Assert2.Near(correct, exact_alignment.AlignmentCost);
        Assert.Equal(exact_alignment.WarpPath_.Count, p.Count);
        for (int i = 0; i < exact_alignment.WarpPath_.Count; ++i) {
            Assert.Equal(exact_alignment.WarpPath_[i], p[i]);
        }
    }

    // Return vector of length `num_polylines` containing correlated random
    // polylines with `num_vertices` vertices each.
    //
    // First, we construct a regularly spaced base loop with `num_vertices`
    // vertices. Then, for each of `num_polylines` iterations, we construct an new
    // loop by uniformly perturbing each point in the base loop by an amount equal
    // to `perturbation` * edge_length in a spherical cap. If `perturbation` is less
    // than 0.5, then we can perturb each point in the second loop by up to 0.5 edge
    // lengths in any direction, which will leave that point with only one possible
    // closest vertex match in the base loop. On the other hand, if `perturbation`
    // is greater than 0.5, then each vertex in the additional loop will more than
    // one match (approximately 2*perturbation + 1) on average in the base loop. The
    // intent of this method is to provide a set of correlated testing lines for
    // benchmarks and fuzz tests.
    private static S2Polyline[] GenPolylines(int num_polylines, int num_vertices, double perturbation) {
        var kLoopRadius = S1Angle.FromRadians(0.01);
        var edge_length = S2.M_2_PI * kLoopRadius / num_vertices;
        var perturbation_radius = perturbation * edge_length;
        var center = S2Testing.RandomPoint();
        var loop = S2Testing.MakeRegularPoints(center, kLoopRadius, num_vertices);

        var polylines = new S2Polyline[num_polylines];

        for (int i = 0; i < num_polylines; ++i) {
            var pts = new List<S2Point>(num_vertices);
            for (int j = 0; j < num_vertices; ++j) {
                pts.Add(S2Testing.SamplePoint(new S2Cap(loop[j], perturbation_radius)));
            }
            polylines[i] = new S2Polyline(pts.ToArray());
        }
        return polylines;
    }
}


internal static class StringExtensions
{
    internal static string NoCR(this string s) => s.Replace("\r", "");
}
