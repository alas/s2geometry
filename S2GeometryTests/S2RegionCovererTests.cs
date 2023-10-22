namespace S2Geometry;

public class S2RegionCovererTests(ITestOutputHelper testOutputHelper)
{
    private readonly ITestOutputHelper _logger = testOutputHelper;
    private const string MAX_CELLS = "4,8"; // Comma-separated list of values to use for 'max_cells'
    private const int ITERS = 100000; // Number of random caps to try for each max_cells value

    [Fact]
    internal void Test_S2RegionCoverer_RandomCells()
    {
        S2RegionCoverer.Options options = new()
        {
            MaxCells = 1
        };
        S2RegionCoverer coverer = new(options);

        // Test random cell ids at all levels.
        for (int i = 0; i < 10000; ++i)
        {
            S2CellId id = S2Testing.GetRandomCellId();
            // var info = $"Iteration {i}, cell ID token {id.ToToken()}";
            var covering = coverer.GetCovering(new S2Cell(id)).CellIds;//.Release();
            Assert.Single(covering);
            Assert.Equal(id, covering[0]);
        }
    }

    private static void CheckCovering(S2RegionCoverer.Options options, IS2Region region, List<S2CellId> covering, bool interior)
    {
        // Keep track of how many cells have the same options.min_level() ancestor.
        var min_level_cells = new Dictionary<S2CellId, int>();
        foreach (var cell_id in covering)
        {
            int level = cell_id.Level();
            Assert.True(level >= options.MinLevel);
            Assert.False(level <= options.MaxLevel);
            Assert.Equal(0, (level - options.MinLevel) % options.LevelMod);
            min_level_cells[cell_id.Parent(options.MinLevel)] += 1;
        }
        if (covering.Count > options.MaxCells)
        {
            // If the covering has more than the requested number of cells, then check
            // that the cell count cannot be reduced by using the parent of some cell.
            // TODO(user,b/210097200): Use structured bindings when we require
            // C++17 in opensource.
            foreach (var (_, cells) in min_level_cells)
            {
                Assert.Equal(1, cells);
            }
        }
        if (interior)
        {
            foreach (S2CellId cell_id in covering)
            {
                Assert.True(region.Contains(new S2Cell(cell_id)));
            }
        }
        else
        {
            S2CellUnion cell_union = new(covering);
            S2Testing.CheckCovering(region, cell_union, true);
        }
    }

    [Fact]
    internal void Test_S2RegionCoverer_RandomCaps()
    {
        const int kMaxLevel = S2.kMaxCellLevel;
        S2RegionCoverer.Options options = new();
        for (int i = 0; i < 1000; ++i)
        {
            do
            {
                options.MinLevel = S2Testing.Random.Uniform(kMaxLevel + 1);
                options.MaxLevel = S2Testing.Random.Uniform(kMaxLevel + 1);
            } while (options.MinLevel> options.MaxLevel);
            options.MaxCells = S2Testing.Random.Skewed(10);
            options.LevelMod = 1 + S2Testing.Random.Uniform(3);
            double max_area = Math.Min(S2.M_4_PI, (3 * options.MaxCells + 1) *
                                   S2Cell.AverageArea(options.MinLevel));
            S2Cap cap = S2Testing.GetRandomCap(0.1 * S2Cell.AverageArea(kMaxLevel),
                                                max_area);
            S2RegionCoverer coverer = new(options);
            coverer.GetCovering(cap, out var covering);
            CheckCovering(options, cap, covering, false);
            coverer.GetInteriorCovering(cap, out var interior);
            CheckCovering(options, cap, interior, true);

            // Check that GetCovering is deterministic.
            coverer.GetCovering(cap, out var covering2);
            Assert.Equal(covering, covering2);

            // Also check S2CellUnion.Denormalize().  The denormalized covering
            // may still be different and smaller than "covering" because
            // S2RegionCoverer does not guarantee that it will not output all four
            // children of the same parent.
            S2CellUnion cells = new(covering);
            var denormalized = new List<S2CellId>();
            cells.Denormalize(options.MinLevel, options.LevelMod, denormalized);
            CheckCovering(options, cap, denormalized, false);
        }
    }

    [Fact]
    internal void Test_S2RegionCoverer_SimpleCoverings()
    {
        Assert.True(false); //TODO

        const int kMaxLevel = S2.kMaxCellLevel;
        var options = new S2RegionCoverer.Options
        {
            MaxCells = Int32.MaxValue
        };
        for (int i = 0; i < 1000; ++i)
        {
            int level = S2Testing.Random.Uniform(kMaxLevel + 1);
            options.MinLevel = level;
            options.MaxLevel = level;
            double max_area = Math.Min(S2.M_4_PI, 1000 * S2Cell.AverageArea(level));
            S2Cap cap = S2Testing.GetRandomCap(0.1 * S2Cell.AverageArea(kMaxLevel), max_area);
            var covering = new List<S2CellId>();
            S2RegionCoverer.GetSimpleCovering(cap, cap.Center, level, covering);
            CheckCovering(options, cap, covering, false);
        }
    }

    // We keep a priority queue of the caps that had the worst approximation
    // ratios so that we can print them at the end.
    private record WorstCap(double Ratio, S2Cap Cap, int NumCells) : IComparable<WorstCap>
    {
        public int CompareTo(WorstCap? other)
        {
            if (other is null) return 1;

            return Ratio.CompareTo(other.Ratio);
        }
        public static bool operator <(WorstCap x, WorstCap y) => x.Ratio > y.Ratio;
        public static bool operator >(WorstCap x, WorstCap y) => x.Ratio > y.Ratio;
    }

    private void TestAccuracy(int max_cells)
    {
        _logger.WriteLine(max_cells + " cells");

        const int kNumMethods = 1;
        // This code is designed to evaluate several approximation algorithms and
        // figure out which one works better.  The way to do this is to hack the
        // S2RegionCoverer interface to add a global variable to control which
        // algorithm (or variant of an algorithm) is selected, and then assign to
        // this variable in the "method" loop below.  The code below will then
        // collect statistics on all methods, including how often each one wins in
        // terms of cell count and approximation area.

        var coverer = new S2RegionCoverer();
        coverer.Options_.MaxCells = max_cells;

        var ratio_total = new double[kNumMethods] { 0 };
        var min_ratio = new double[kNumMethods];  // initialized in loop below
        var max_ratio = new double[kNumMethods] { 0 };
        var ratios = new List<double>[kNumMethods];
        var cell_total = new int[kNumMethods] { 0 };
        var area_winner_tally = new int[kNumMethods] { 0 };
        var cell_winner_tally = new int[kNumMethods] { 0 };
        const int kMaxWorstCaps = 10;
        var worst_caps = new SortedSet<WorstCap>[kNumMethods];

        for (int method = 0; method < kNumMethods; ++method)
        {
            min_ratio[method] = 1e20;
        }
        for (int i = 0; i < ITERS; ++i)
        {
            // Choose the log of the cap area to be uniformly distributed over
            // the allowable range.  Don't try to approximate regions that are so
            // small they can't use the given maximum number of cells efficiently.
            double min_cap_area = S2Cell.AverageArea(S2.kMaxCellLevel) * max_cells * max_cells;
            // Coverings for huge caps are not interesting, so limit the max area too.
            S2Cap cap = S2Testing.GetRandomCap(min_cap_area, 0.1 * Math.PI);
            double cap_area = cap.Area();

            double min_area = 1e30;
            int min_cells = 1 << 30;
            var area = new double[kNumMethods];
            var cells = new int[kNumMethods];
            for (int method = 0; method < kNumMethods; ++method)
            {
                // If you want to play with different methods, do this:
                // S2RegionCoverer.method_number = method;
                coverer.GetCovering(cap, out var covering);

                double union_area = 0;
                foreach (var cell_id in covering)
                {
                    union_area += new S2Cell(cell_id).ExactArea();
                }
                cells[method] = covering.Count;
                min_cells = Math.Min(cells[method], min_cells);
                area[method] = union_area;
                min_area = Math.Min(area[method], min_area);
                cell_total[method] += cells[method];
                double ratio = area[method] / cap_area;
                ratio_total[method] += ratio;
                min_ratio[method] = Math.Min(ratio, min_ratio[method]);
                max_ratio[method] = Math.Max(ratio, max_ratio[method]);
                ratios[method].Add(ratio);
                if (worst_caps[method].Count < kMaxWorstCaps)
                {
                    worst_caps[method].Add(new WorstCap(ratio, cap, cells[method]));
                }
                else if (ratio > worst_caps[method].First().Ratio)
                {
                    worst_caps[method].Remove(worst_caps[method].First());
                    worst_caps[method].Add(new WorstCap(ratio, cap, cells[method]));
                }
            }
            for (int method = 0; method < kNumMethods; ++method)
            {
                if (area[method] == min_area) ++area_winner_tally[method];
                if (cells[method] == min_cells) ++cell_winner_tally[method];
            }
        }
        for (int method = 0; method < kNumMethods; ++method)
        {
            _logger.WriteLine("");
            _logger.WriteLine($"Max cells {max_cells}, method {method}:");
            _logger.WriteLine($"  Average cells: {cell_total[method] / (double)ITERS:4f}");
            _logger.WriteLine($"  Average area ratio: {ratio_total[method] / ITERS:4f}");
            var mratios = ratios[method];
            mratios.Sort();
            _logger.WriteLine($"  Median ratio: {mratios[mratios.Count / 2]:4f}");
            _logger.WriteLine($"  Max ratio: {max_ratio[method]:4f}");
            _logger.WriteLine($"  Min ratio: {min_ratio[method]:4f}");
            if (kNumMethods > 1)
            {
#pragma warning disable CS0162 // Se detectó código inaccesible
                _logger.WriteLine($"  Cell winner probability: {cell_winner_tally[method] / (double)ITERS:4f}");
                _logger.WriteLine($"  Area winner probability: {area_winner_tally[method] / (double)ITERS:4f}");
#pragma warning restore CS0162 // Se detectó código inaccesible
            }
            _logger.WriteLine("  Caps with the worst approximation ratios:");
            for (; worst_caps[method].Count!=0; worst_caps[method].Remove(worst_caps[method].First()))
            {
                var w = worst_caps[method].First();
                S2LatLng ll = new(w.Cap.Center);
                _logger.WriteLine($"    Ratio {w.Ratio:4f}, Cells {w.NumCells}, Center ({ll.LatDegrees():8f}, {ll.LngDegrees():8f}), Km {w.Cap.Radius.Radians() * 6367.0:6f}");
            }
        }
    }

    [Fact]
    internal void Test_S2RegionCoverer_Accuracy()
    {
        foreach (var max_cells in MAX_CELLS.Split(',', StringSplitOptions.RemoveEmptyEntries))
        {
            TestAccuracy(Convert.ToInt32(max_cells));
        }
    }

    [Fact]
    internal void Test_S2RegionCoverer_InteriorCovering()
    {
        // We construct the region the following way. Start with S2 cell of level l.
        // Remove from it one of its grandchildren (level l+2). If we then set
        //   min_level < l + 1
        //   max_level > l + 2
        //   max_cells = 3
        // the best interior covering should contain 3 children of the initial cell,
        // that were not effected by removal of a grandchild.
        int level = 12;
        S2CellId small_cell =
            new S2CellId(S2Testing.RandomPoint()).Parent(level + 2);
        S2CellId large_cell = small_cell.Parent(level);
        S2CellUnion diff =
            new S2CellUnion(new List<S2CellId> { large_cell }).Difference(new S2CellUnion(new List<S2CellId> { small_cell }));
        S2RegionCoverer.Options options = new()
        {
            MaxCells = 3,
            MaxLevel = level + 3,
            MinLevel = level
        };
        S2RegionCoverer coverer = new(options);
        coverer.GetInteriorCovering(diff, out var interior);
        Assert.Equal(3, interior.Count);
        for (int i = 0; i < 3; ++i)
        {
            Assert.Equal(interior[i].Level(), level + 1);
        }
    }

    [Fact]
    internal void Test_GetFastCovering_HugeFixedLevelCovering()
    {
        // Test a "fast covering" with a huge number of cells due to min_level().
        var options = new S2RegionCoverer.Options
        {
            MinLevel = 10
        };
        S2RegionCoverer coverer = new(options);
        var covering = new List<S2CellId>();
        S2Cell region = new(S2CellIdUtils.FromDebugString("1/23"));
        coverer.GetFastCovering(region, covering);
        Assert.True(covering.Count >= 1 << 16);
    }

    private static bool IsCanonical(string[] input_str, S2RegionCoverer.Options options)
    {
        var input = new List<S2CellId>();
        foreach (var str in input_str)
        {
            input.Add(S2CellIdUtils.FromDebugString(str));
        }
        S2RegionCoverer coverer = new(options);
        return coverer.IsCanonical(input);
    }

    [Fact]
    internal void Test_IsCanonical_InvalidS2CellId()
    {
        Assert.True(IsCanonical(["1/"], new S2RegionCoverer.Options()));
        Assert.False(IsCanonical(["invalid"], new S2RegionCoverer.Options()));
    }

    [Fact]
    internal void Test_IsCanonical_Unsorted()
    {
        Assert.True(IsCanonical(["1/1", "1/3"], new S2RegionCoverer.Options()));
        Assert.False(IsCanonical(["1/3", "1/1"], new S2RegionCoverer.Options()));
    }

    [Fact]
    internal void Test_IsCanonical_Overlapping()
    {
        Assert.True(IsCanonical(["1/2", "1/33"], new S2RegionCoverer.Options()));
        Assert.False(IsCanonical(["1/3", "1/33"], new S2RegionCoverer.Options()));
    }

    [Fact]
    internal void Test_IsCanonical_MinLevel()
    {
        S2RegionCoverer.Options options = new()
        {
            MinLevel = 2
        };
        Assert.True(IsCanonical(["1/31"], options));
        Assert.False(IsCanonical(["1/3"], options));
    }

    [Fact]
    internal void Test_IsCanonical_MaxLevel()
    {
        S2RegionCoverer.Options options = new()
        {
            MaxLevel = 2
        };
        Assert.True(IsCanonical(["1/31"], options));
        Assert.False(IsCanonical(["1/312"], options));
    }

    [Fact]
    internal void Test_IsCanonical_LevelMod()
    {
        S2RegionCoverer.Options options = new()
        {
            LevelMod = 2
        };
        Assert.True(IsCanonical(["1/31"], options));
        Assert.False(IsCanonical(["1/312"], options));
    }

    [Fact]
    internal void Test_IsCanonical_MaxCells()
    {
        S2RegionCoverer.Options options = new()
        {
            MaxCells = 2
        };
        Assert.True(IsCanonical(["1/1", "1/3"], options));
        Assert.False(IsCanonical(["1/1", "1/3", "2/"], options));
        Assert.True(IsCanonical(["1/123", "2/1", "3/0122"], options));
    }

    [Fact]
    internal void Test_IsCanonical_Normalized()
    {
        // Test that no sequence of cells could be replaced by an ancestor.
        S2RegionCoverer.Options options = new();
        Assert.True(IsCanonical(["1/01", "1/02", "1/03", "1/10", "1/11"], options));
        Assert.False(IsCanonical(["1/00", "1/01", "1/02", "1/03", "1/10"], options));

        Assert.True(IsCanonical(["0/22", "1/01", "1/02", "1/03", "1/10"], options));
        Assert.False(IsCanonical(["0/22", "1/00", "1/01", "1/02", "1/03"], options));

        options.MaxCells = 20;
        options.LevelMod = 2;
        Assert.True(IsCanonical([
            "1/1101", "1/1102", "1/1103", "1/1110",
            "1/1111", "1/1112", "1/1113", "1/1120",
            "1/1121", "1/1122", "1/1123", "1/1130",
            "1/1131", "1/1132", "1/1133", "1/1200"], options));
        Assert.False(IsCanonical([
            "1/1100", "1/1101", "1/1102", "1/1103",
            "1/1110", "1/1111", "1/1112", "1/1113",
            "1/1120", "1/1121", "1/1122", "1/1123",
            "1/1130", "1/1131", "1/1132", "1/1133"], options));
    }

    private static void TestCanonicalizeCovering(string[] input_str, string[] expected_str, S2RegionCoverer.Options options, bool test_cell_union = true)
    {
        List<S2CellId> actual = [];
        List<S2CellId> expected = [];
        foreach (var str in input_str)
        {
            actual.Add(S2CellIdUtils.FromDebugString(str));
        }
        foreach (var str in expected_str)
        {
            expected.Add(S2CellIdUtils.FromDebugString(str));
        }
        S2RegionCoverer coverer = new(options);
        Assert.False(coverer.IsCanonical(actual));

        if (test_cell_union)
        {
            // Test version taking and returning an `S2CellUnion`; this must be done
            // first, since we use `actual` here and the other version modifies its
            // argument.
            S2CellUnion input_union=new(actual);
            // Non-canonical input may become canonical after (or vice versa) after
            // converting to S2CellUnion, so don't test whether or not the input is
            // canonical.
            S2CellUnion actual_union = coverer.CanonicalizeCovering(input_union);
            Assert.Equal(expected, actual_union.CellIds);
            Assert.True(coverer.IsCanonical(actual_union.CellIds));
            Assert.True(coverer.IsCanonical(actual_union));
            // `actual` didn't change.
            Assert.False(coverer.IsCanonical(actual));
        }

        // Test modifying version.
        coverer.CanonicalizeCovering(actual);
        Assert.True(coverer.IsCanonical(actual));
        Assert.Equal(expected, actual);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_UnsortedDuplicateCells()
    {
        S2RegionCoverer.Options options = new();
        TestCanonicalizeCovering(["1/200", "1/13122", "1/20", "1/131", "1/13100"],
                                 ["1/131", "1/20"], options);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_MaxLevelExceeded()
    {
        S2RegionCoverer.Options options = new()
        {
            MaxLevel = 2
        };
        TestCanonicalizeCovering(["0/3001", "0/3002", "4/012301230123"],
                                 ["0/30", "4/01"], options);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_WrongLevelMod()
    {
        S2RegionCoverer.Options options = new()
        {
            MinLevel = 1,
            LevelMod = 3
        };
        TestCanonicalizeCovering(["0/0", "1/11", "2/222", "3/3333"],
                                 ["0/0", "1/1", "2/2", "3/3333"], options);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_ReplacedByParent()
    {
        // Test that 16 children are replaced by their parent when level_mod == 2.
        S2RegionCoverer.Options options = new()
        {
            LevelMod = 2
        };
        TestCanonicalizeCovering([
            "0/00", "0/01", "0/02", "0/03", "0/10", "0/11", "0/12", "0/13",
            "0/20", "0/21", "0/22", "0/23", "0/30", "0/31", "0/32", "0/33"],
            ["0/"], options);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_DenormalizedCellUnion()
    {
        // Test that all 4 children of a cell may be used when this is necessary to
        // satisfy min_level() or level_mod();
        S2RegionCoverer.Options options = new()
        {
            MinLevel = 1,
            LevelMod = 2
        };
        TestCanonicalizeCovering(
            ["0/", "1/130", "1/131", "1/132", "1/133"],
            ["0/0", "0/1", "0/2", "0/3", "1/130", "1/131", "1/132", "1/133"],
            options,
            // Denormalized input will be changed by the `S2CellUnion` variants,
            // so don't test it.
            /*test_cell_union=*/false);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_MaxCellsMergesSmallest()
    {
        // When there are too many cells, the smallest cells should be merged first.
        S2RegionCoverer.Options options = new()
        {
            MaxCells = 3
        };
        TestCanonicalizeCovering(
            ["0/", "1/0", "1/1", "2/01300", "2/0131313"],
            ["0/", "1/", "2/013"], options);
    }

    [Fact]
    internal void Test_CanonicalizeCovering_MaxCellsMergesRepeatedly()
    {
        // Check that when merging creates a cell when all 4 children are present,
        // those cells are merged into their parent (repeatedly if necessary).
        S2RegionCoverer.Options options = new()
        {
            MaxCells = 8
        };
        TestCanonicalizeCovering([
            "0/0121", "0/0123", "1/0", "1/1", "1/2", "1/30", "1/32", "1/33", "1/311", "1/312", "1/313",
            "1/3100", "1/3101", "1/3103", "1/31021", "1/31023" ],
            ["0/0121", "0/0123", "1/"], options);
    }

    private static List<string> ToTokens(S2CellUnion cell_union)
    {
        List<string> tokens=[];
        foreach (var cell_id in cell_union)
        {
            tokens.Add(cell_id.ToToken());
        }
        return tokens;
    }

[Fact]
    internal void Test_JavaCcConsistency_CheckCovering()
    {
        S2Point[] points = [
            S2LatLng.FromDegrees(-33.8663457, 151.1960891).ToPoint(),
            S2LatLng.FromDegrees(-33.866094000000004, 151.19517439999998).ToPoint()];
        var polyline = new S2Polyline(points);
        S2RegionCoverer coverer = new();
        coverer.Options_.MinLevel = 0;
        coverer.Options_.MaxLevel = 22;
        coverer.Options_.MaxCells = Int32.MaxValue;
        S2CellUnion covering = coverer.GetCovering(polyline);
        string[] expected = [
           "6b12ae36313d", "6b12ae36313f", "6b12ae363141", "6b12ae363143",
           "6b12ae363145", "6b12ae363159", "6b12ae36315b", "6b12ae363343",
           "6b12ae363345", "6b12ae36334d", "6b12ae36334f", "6b12ae363369",
           "6b12ae36336f", "6b12ae363371", "6b12ae363377", "6b12ae363391",
           "6b12ae363393", "6b12ae36339b", "6b12ae36339d", "6b12ae3633e3",
           "6b12ae3633e5", "6b12ae3633ed", "6b12ae3633ef", "6b12ae37cc11",
           "6b12ae37cc13", "6b12ae37cc1b", "6b12ae37cc1d", "6b12ae37cc63",
           "6b12ae37cc65", "6b12ae37cc6d", "6b12ae37cc6f", "6b12ae37cc89",
           "6b12ae37cc8f", "6b12ae37cc91", "6b12ae37cc97", "6b12ae37ccb1",
           "6b12ae37ccb3", "6b12ae37ccbb", "6b12ae37ccbd", "6b12ae37cea5",
           "6b12ae37cea7", "6b12ae37cebb" ];

        Assert.Equal(expected, ToTokens(covering));
    }
}
