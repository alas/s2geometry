namespace S2Geometry;

public class S2RegionTermIndexerTests
{
    private const int iters = 400; // number of iterations for testing
    private readonly ITestOutputHelper _logger;

    internal S2RegionTermIndexerTests(ITestOutputHelper logger) { _logger = logger; }

    // We run one test case for each combination of space vs. time optimization,
    // and indexing regions vs. only points.

    [Fact]
    internal void Test_S2RegionTermIndexer_IndexRegionsQueryRegionsOptimizeTime()
    {
        var options = new S2RegionTermIndexer.Options();
        options.OptimizeForSpace = false;       // Optimize for time.
        options.MinLevel = (0);                    // Use face cells.
        options.MaxLevel = (16);
        options.MaxCells = (20);
        TestRandomCaps(options, QueryType.CAP);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_IndexRegionsQueryPointsOptimizeTime()
    {
        var options = new S2RegionTermIndexer.Options();
        options.OptimizeForSpace = (false);       // Optimize for time.
        options.MinLevel = (0);                    // Use face cells.
        options.MaxLevel = (16);
        options.MaxCells = (20);
        TestRandomCaps(options, QueryType.POINT);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_IndexRegionsQueryRegionsOptimizeTimeWithLevelMod()
    {
        var options = new S2RegionTermIndexer.Options();
        options.OptimizeForSpace = (false);       // Optimize for time.
        options.MinLevel = (6);                    // Constrain min/max levels.
        options.MaxLevel = (12);
        options.LevelMod = (3);
        TestRandomCaps(options, QueryType.CAP);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_IndexRegionsQueryRegionsOptimizeSpace()
    {
        var options = new S2RegionTermIndexer.Options();
        options.OptimizeForSpace = (true);        // Optimize for space.
        options.MinLevel = (4);
        options.MaxLevel = (S2.kMaxCellLevel);  // Use leaf cells.
        options.MaxCells = (8);
        TestRandomCaps(options, QueryType.CAP);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_IndexPointsQueryRegionsOptimizeTime()
    {
        var options = new S2RegionTermIndexer.Options();
        options.OptimizeForSpace = (false);       // Optimize for time.
        options.MinLevel = (0);                    // Use face cells.
        options.MaxLevel = (S2.kMaxCellLevel);
        options.LevelMod = (2);
        options.MaxCells = (20);
        options.IndexContainsPointsOnly = (true);
        TestRandomCaps(options, QueryType.CAP);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_IndexPointsQueryRegionsOptimizeSpace()
    {
        var options = new S2RegionTermIndexer.Options
        {
            OptimizeForSpace = true,        // Optimize for space.
            IndexContainsPointsOnly = true,
        };
        // Use default parameter values.
        TestRandomCaps(options, QueryType.CAP);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_MarkerCharacter()
    {
        S2RegionTermIndexer.Options options = new();
        options.MinLevel = 20;
        options.MaxLevel = 20;

        S2RegionTermIndexer indexer = new(options);
        S2Point point = S2LatLng.FromDegrees(10, 20).ToPoint();
        Assert.Equal(indexer.Options_.MarkerCharacter, '$');
        Assert.Equal(indexer.GetQueryTerms(point, ""),
            new List<string>(){ "11282087039", "$11282087039"});

        indexer.Options_.MarkerCharacter = ':';
        Assert.Equal(indexer.Options_.MarkerCharacter, ':');
        Assert.Equal(indexer.GetQueryTerms(point, ""),
            new List<string>(){ "11282087039", ":11282087039"});
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_MaxLevelSetLoosely()
    {
        // Test that correct terms are generated even when (max_level - min_level)
        // is not a multiple of level_mod.
        var options = new S2RegionTermIndexer.Options();
        options.MinLevel = (1);
        options.LevelMod = (2);
        options.MaxLevel = (19);
        var indexer1 = new S2RegionTermIndexer(options);
        options.MaxLevel = (20);
        var indexer2 = new S2RegionTermIndexer(options);

        S2Point point = S2Testing.RandomPoint();
        Assert.Equal(indexer1.GetIndexTerms(point, ""),
                  indexer2.GetIndexTerms(point, ""));
        Assert.Equal(indexer1.GetQueryTerms(point, ""),
                  indexer2.GetQueryTerms(point, ""));

        S2Cap cap = S2Testing.GetRandomCap(0.0, 1.0);  // Area range.
        Assert.Equal(indexer1.GetIndexTerms(cap, ""),
                  indexer2.GetIndexTerms(cap, ""));
        Assert.Equal(indexer1.GetQueryTerms(cap, ""),
                  indexer2.GetQueryTerms(cap, ""));
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_MoveConstructor()
    {
        var x = new S2RegionTermIndexer(new S2RegionTermIndexer.Options());
        x.Options_.MaxCells = (12345);
        var y = x;
        Assert.Equal(12345, y.Options_.MaxCells);
    }

    [Fact]
    internal void Test_S2RegionTermIndexer_MoveAssignmentOperator()
    {
        var x = new S2RegionTermIndexer(new S2RegionTermIndexer.Options());
        x.Options_.MaxCells = (12345);
        var y = new S2RegionTermIndexer(new S2RegionTermIndexer.Options());
        y.Options_.MaxCells = (0);
        y = x;
        Assert.Equal(12345, y.Options_.MaxCells);
    }

    private void TestRandomCaps(S2RegionTermIndexer.Options options, QueryType query_type)
    {
        // This function creates an index consisting either of points (if
        // options.index_contains_points_only() is true) or S2Caps of random size.
        // It then executes queries consisting of points (if query_type == POINT)
        // or S2Caps of random size (if query_type == CAP).
        var indexer = new S2RegionTermIndexer(options);
        var coverer = new S2RegionCoverer(options);
        var caps = new List<S2Cap>();
        var coverings = new List<S2CellUnion>();
        var index = new Dictionary<string, List<int>>();
        int index_terms = 0, query_terms = 0;
        for (int i = 0; i < iters; ++i)
        {
            // Choose the region to be indexed: either a single point or a cap
            // of random size (up to a full sphere).
            S2Cap cap;
            List<string> terms;
            if (options.IndexContainsPointsOnly)
            {
                cap = S2Cap.FromPoint(S2Testing.RandomPoint());
                terms = indexer.GetIndexTerms(cap.Center, "");
            }
            else
            {
                cap = S2Testing.GetRandomCap(
                    0.3 * S2Cell.AverageArea(options.MaxLevel),
                    4.0 * S2Cell.AverageArea(options.MinLevel));
                terms = indexer.GetIndexTerms(cap, "");
            }
            caps.Add(cap);
            coverings.Add(coverer.GetCovering(cap));
            foreach (var term in terms)
            {
                if (!index.ContainsKey(term))
                    index.Add(term, new List<int>());

                index[term].Add(i);
            }
            index_terms += terms.Count;
        }
        for (int i = 0; i < iters; ++i)
        {
            // Choose the region to be queried: either a random point or a cap of
            // random size.
            S2Cap cap;
            List<string> terms;
            if (query_type == QueryType.CAP)
            {
                cap = S2Cap.FromPoint(S2Testing.RandomPoint());
                terms = indexer.GetQueryTerms(cap.Center, "");
            }
            else
            {
                cap = S2Testing.GetRandomCap(
                    0.3 * S2Cell.AverageArea(options.MaxLevel),
                    4.0 * S2Cell.AverageArea(options.MinLevel));
                terms = indexer.GetQueryTerms(cap, "");
            }
            // Compute the expected results of the S2Cell query by brute force.
            S2CellUnion covering = coverer.GetCovering(cap);
            var expected = new List<int>();
            var actual = new List<int>();
            for (int j = 0; j < caps.Count; ++j)
            {
                if (covering.Intersects(coverings[j]))
                {
                    expected.Add(j);
                }
            }
            foreach (var term in terms)
            {
                actual.AddRange(index[term]);
            }
            Assert.Equal(expected, actual);
            query_terms += terms.Count;
        }
        _logger.WriteLine($"Index terms/doc: {((double)index_terms) / iters:2f},  Query terms/doc: {((double)query_terms) / iters:2f}");
    }

    private enum QueryType { POINT, CAP };
}
