namespace S2Geometry;

using LabelledCell = S2CellIndex.LabelledCell;
using Target = S2MinDistanceTarget;

public class S2ClosestCellQueryTests
{
    // The approximate radius of S2Cap from which query cells are chosen.
    private static readonly S1Angle kTestCapRadius = S2Testing.KmToAngle(10);
    private const int kNumIndexes = 20;
    private const int kNumCells = 100;
    private const int kNumQueries = 100;
    private static ITestOutputHelper _logger;

    public S2ClosestCellQueryTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    public void Test_S2ClosestCellQuery_NoCells() {
        S2CellIndex index = new();
        index.Build();
        S2ClosestCellQuery query = new(index);
        S2ClosestCellQuery.PointTarget target = new(new S2Point(1, 0, 0));
        var result = query.FindClosestCell(target);
        Assert.Equal(S1ChordAngle.Infinity, result.Distance);
        Assert.Equal(S2CellId.None, result.CellId);
        Assert.Equal(-1, result.Label);
        Assert.True(result.IsEmpty());
        Assert.Equal(S1ChordAngle.Infinity, query.GetDistance(target));
    }

    [Fact]
    public void Test_S2ClosestCellQuery_OptionsNotModified() {
        // Tests that FindClosestCell(), GetDistance(), and IsDistanceLess() do not
        // modify query.Options_, even though all of these methods have their own
        // specific options requirements.
        S2ClosestCellQuery.Options options = new()
        {
            MaxResults = 3,
            MaxDistance = S1ChordAngle.FromDegrees(3),
            MaxError = S1ChordAngle.FromDegrees(0.001)
        };
        S2CellIndex index = new();
        index.Add(new S2CellId(MakePointOrDie("1:1")), 1);
        index.Add(new S2CellId(MakePointOrDie("1:2")), 2);
        index.Add(new S2CellId(MakePointOrDie("1:3")), 3);
        index.Build();
        S2ClosestCellQuery query = new(index, options);
        S2ClosestCellQuery.PointTarget target = new(MakePointOrDie("2:2"));
        Assert.Equal(2, query.FindClosestCell(target).Label);
        Assert2.Near(1.0, query.GetDistance(target).Degrees(), 1e-7);
        Assert.True(query.IsDistanceLess(target, S1ChordAngle.FromDegrees(1.5)));

        // Verify that none of the options above were modified.
        Assert.Equal(options.MaxResults, query.Options_.MaxResults);
        Assert.Equal(options.MaxDistance, query.Options_.MaxDistance);
        Assert.Equal(options.MaxError, query.Options_.MaxError);
    }

    [Fact]
    public void Test_S2ClosestCellQuery_DistanceEqualToLimit() {
        // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
        // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
        // the distance to the target exactly equals the chosen limit.
        S2CellId id0 = new(MakePointOrDie("23:12")), id1 = new(MakePointOrDie("47:11"));
        S2CellIndex index = new();
        index.Add(id0, 0);
        index.Build();
        S2ClosestCellQuery query = new(index);

        // Start with two identical cells and a zero distance.
        S2ClosestCellQuery.CellTarget target0 = new(new S2Cell(id0));
        S1ChordAngle dist0 = S1ChordAngle.Zero;
        Assert.False(query.IsDistanceLess(target0, dist0));
        Assert.True(query.IsDistanceLessOrEqual(target0, dist0));
        Assert.True(query.IsConservativeDistanceLessOrEqual(target0, dist0));

        // Now try two cells separated by a non-zero distance.
        S2ClosestCellQuery.CellTarget target1 = new(new S2Cell(id1));
        S1ChordAngle dist1 = new S2Cell(id0).Distance(new S2Cell(id1));
        Assert.False(query.IsDistanceLess(target1, dist1));
        Assert.True(query.IsDistanceLessOrEqual(target1, dist1));
        Assert.True(query.IsConservativeDistanceLessOrEqual(target1, dist1));
    }

    [Fact]
    public void Test_S2ClosestCellQuery_TargetPointInsideIndexedCell() {
        // Tests a target point in the interior of an indexed cell.
        S2CellId cell_id = MakeCellIdOrDie("4/012");
        S2CellIndex index = new();
        index.Add(cell_id, 1);
        index.Build();
        S2ClosestCellQuery query = new(index);
        S2ClosestCellQuery.PointTarget target = new(cell_id.ToPoint());
        var result = query.FindClosestCell(target);
        Assert.Equal(S1ChordAngle.Zero, result.Distance);
        Assert.Equal(cell_id, result.CellId);
        Assert.Equal(1, result.Label);
    }

    [Fact]
    public void Test_S2ClosestCellQuery_EmptyTargetOptimized() {
        // Ensure that the optimized algorithm handles empty targets when a distance
        // limit is specified.
        S2CellIndex index = new();
        for (int i = 0; i < 1000; ++i) {
            index.Add(S2Testing.GetRandomCellId(), i);
        }
        index.Build();
        S2ClosestCellQuery query = new(index);
        query.Options_.MaxDistance = new S1ChordAngle(S1Angle.FromRadians(1e-5));
        MutableS2ShapeIndex target_index=new();
        S2ClosestCellQuery.ShapeIndexTarget target = new(target_index);
        Assert.Empty(query.FindClosestCells(target));
    }

    [Fact]
    public void Test_S2ClosestCellQuery_EmptyCellUnionTarget() {
        // Verifies that distances are measured correctly to empty S2CellUnion
        // targets.
        S2ClosestCellQuery.CellUnionTarget target = new(new S2CellUnion());

        S2CellIndex empty_index = new();
        empty_index.Build();
        S2ClosestCellQuery empty_query = new(empty_index);
        Assert.Equal(S1ChordAngle.Infinity, empty_query.GetDistance(target));

        S2CellIndex one_cell_index = new();
        one_cell_index.Add(MakeCellIdOrDie("1/123123"), 1);
        one_cell_index.Build();
        S2ClosestCellQuery one_cell_query = new(one_cell_index);
        Assert.Equal(S1ChordAngle.Infinity, one_cell_query.GetDistance(target));
    }

    [Fact]
    public void Test_S2ClosestCellQuery_PointCloudCells() {
        TestWithIndexFactory(new PointCloudCellIndexFactory(),
                             kNumIndexes, kNumCells, kNumQueries);
    }

    [Fact]
    public void Test_S2ClosestCellQuery_CapsCells() {
        TestWithIndexFactory(new CapsCellIndexFactory(16 /*max_cells_per_cap*/, 0.1 /*density*/),
            kNumIndexes, kNumCells, kNumQueries);
    }

    [Fact]
    public void Test_S2ClosestCellQuery_ConservativeCellDistanceIsUsed() {
        // Don't use google.FlagSaver, so it works in opensource without gflags.
        int saved_seed = S2Testing.Random.RandomSeed;
        // These specific test cases happen to fail if max_error() is not properly
        // taken into account when measuring distances to S2ShapeIndex cells.
        foreach (int seed in new[] { 32, 109, 253, 342, 948, 1535, 1884, 1887, 2133 }) {
            S2Testing.Random.RandomSeed = seed;
            TestWithIndexFactory(new PointCloudCellIndexFactory(), 5, 100, 10);
        }
        S2Testing.Random.RandomSeed = saved_seed;
    }

    // Use "query" to find the closest cell(s) to the given target, and extract
    // the query results into the given vector.  Also verify that the results
    // satisfy the search criteria.
    private static void GetClosestCells(Target target, S2ClosestCellQuery query, List<(S1ChordAngle, LabelledCell)> results)
    {
        var query_results = query.FindClosestCells(target);
        Assert.True(query_results.Count <= query.Options_.MaxResults);
        var region = query.Options_.Region;
        if (region == null && query.Options_.MaxDistance== S1ChordAngle.Infinity)
        {
            // We can predict exactly how many cells should be returned.
            Assert.Equal(Math.Min(query.Options_.MaxResults,
                               query.Index().NumCells()),
                      query_results.Count);
        }
        foreach (var result in query_results)
        {
            // Check that the cell satisfies the region() condition.
            if (region != null) Assert.True(region.MayIntersect(new S2Cell(result.CellId)));

            // Check that it satisfies the max_distance() condition.
            Assert.True(result.Distance < query.Options_.MaxDistance);
            results.Add((result.Distance, new LabelledCell(result.CellId, result.Label)));
        }
    }

    private static void TestFindClosestCells(Target target, S2ClosestCellQuery query)
    {
        List<(S1ChordAngle, LabelledCell)> expected = new(), actual = new();
        query.Options_.UseBruteForce = (true);
        GetClosestCells(target, query, expected);
        query.Options_.UseBruteForce = (false);
        GetClosestCells(target, query, actual);
        Assert.True(S2TestingCheckDistance<LabelledCell, S1ChordAngle>
            .CheckDistanceResults(expected, actual,
                query.Options_.MaxResults,
                query.Options_.MaxDistance,
                query.Options_.MaxError,
                _logger.WriteLine));

        if (!expected.Any()) return;

        // Note that when options.max_error() > 0, expected[0].distance() may not be
        // the minimum distance.  It is never larger by more than max_error(), but
        // the actual value also depends on max_results().
        //
        // Here we verify that GetDistance() and IsDistanceLess() return results
        // that are consistent with the max_error() setting.
        S1ChordAngle max_error = query.Options_.MaxError;
        S1ChordAngle min_distance = expected[0].Item1;
        Assert.True(query.GetDistance(target) <= min_distance + max_error);

        // Test IsDistanceLess().
        Assert.False(query.IsDistanceLess(target, min_distance - max_error));
        Assert.True(query.IsConservativeDistanceLessOrEqual(target, min_distance));
    }

    // The running time of this test is proportional to
    //    (num_indexes + num_queries) * num_cells.
    // (Note that every query is checked using the brute force algorithm.)
    private static void TestWithIndexFactory(ICellIndexFactory factory,
        int num_indexes, int num_cells, int num_queries)
    {
        // Build a set of S2CellIndexes containing the desired geometry.
        List<S2Cap> index_caps = new();
        List<S2CellIndex> indexes = new();
        for (int i = 0; i < num_indexes; ++i)
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
            index_caps.Add(new S2Cap(S2Testing.RandomPoint(), kTestCapRadius));
            var index = new S2CellIndex();
            factory.AddCells(index_caps.Last(), num_cells, index);
            index.Build();
            indexes.Add(index);
        }
        for (int i = 0; i < num_queries; ++i)
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
            int i_index = S2Testing.Random.Uniform(num_indexes);
            S2Cap index_cap = index_caps[i_index];

            // Choose query points from an area approximately 4x larger than the
            // geometry being tested.
            S1Angle query_radius = 2 * index_cap.RadiusAngle();
            S2Cap query_cap = new(index_cap.Center, query_radius);
            S2ClosestCellQuery query = new(indexes[i_index]);

            // Occasionally we don't set any limit on the number of result cells.
            // (This may return all cells if we also don't set a distance limit.)
            if (!S2Testing.Random.OneIn(10))
            {
                query.Options_.MaxResults = (1 + S2Testing.Random.Uniform(10));
            }
            // We set a distance limit 2/3 of the time.
            if (!S2Testing.Random.OneIn(3))
            {
                query.Options_.MaxDistance = new S1ChordAngle(
                    S2Testing.Random.RandDouble() * query_radius);
            }
            if (S2Testing.Random.OneIn(2))
            {
                // Choose a maximum error whose logarithm is uniformly distributed over
                // a reasonable range, except that it is sometimes zero.
                query.Options_.MaxError = new S1ChordAngle(S1Angle.FromRadians(
                    Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius.Radians));
            }
            S2LatLngRect filter_rect = S2LatLngRect.FromCenterSize(
                new S2LatLng(S2Testing.SamplePoint(query_cap)),
                new S2LatLng(S2Testing.Random.RandDouble() * kTestCapRadius,
                         S2Testing.Random.RandDouble() * kTestCapRadius));
            if (S2Testing.Random.OneIn(5))
            {
                query.Options_.Region = (filter_rect);
            }
            int target_type = S2Testing.Random.Uniform(5);
            if (target_type == 0)
            {
                // Find the cells closest to a given point.
                S2Point point = S2Testing.SamplePoint(query_cap);
                S2ClosestCellQuery.PointTarget target = new(point);
                TestFindClosestCells(target, query);
            }
            else if (target_type == 1)
            {
                // Find the cells closest to a given edge.
                S2Point a = S2Testing.SamplePoint(query_cap);
                S2Point b = S2Testing.SamplePoint(
                    new S2Cap(a, Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius));
                S2ClosestCellQuery.EdgeTarget target = new(a, b);
                TestFindClosestCells( target, query);
            }
            else if (target_type == 2)
            {
                // Find the cells closest to a given cell.
                int min_level = S2.kMaxDiag.GetLevelForMaxValue(query_radius.Radians);
                int level = min_level + S2Testing.Random.Uniform(
                    S2.kMaxCellLevel - min_level + 1);
                S2Point a = S2Testing.SamplePoint(query_cap);
                S2Cell cell = new(new S2CellId(a).Parent(level));
                S2ClosestCellQuery.CellTarget target = new(cell);
                TestFindClosestCells(target, query);
            }
            else if (target_type == 3)
            {
                // Find the cells closest to an S2Cap covering.
                S2Cap cap = new(S2Testing.SamplePoint(query_cap),
                          0.1 * Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius);
                S2RegionCoverer coverer = new();
                coverer.Options_.MaxCells=(16);
                S2ClosestCellQuery.CellUnionTarget target = new(coverer.GetCovering(cap));
                TestFindClosestCells(target, query);
            }
            else
            {
                Assert.Equal(4, target_type);
                MutableS2ShapeIndex target_index = new();
                new FractalLoopShapeIndexFactory().AddEdges(index_cap, 100, target_index);
                S2ClosestCellQuery.ShapeIndexTarget target = new(target_index)
                {
                    IncludeInteriors = S2Testing.Random.OneIn(2)
                };
                TestFindClosestCells(target, query);
            }
        }
    }

    // An abstract class that adds cells to an S2CellIndex for benchmarking.
    private interface ICellIndexFactory
    {
        // Requests that approximately "num_cells" cells located within the given
        // S2Cap bound should be added to "index".
        void AddCells(S2Cap index_cap, int num_cells, S2CellIndex index);
    }

    // Generates a cloud of points that approximately fills the given S2Cap, and
    // adds a leaf S2CellId for each one.
    private struct PointCloudCellIndexFactory : ICellIndexFactory
    {
        public void AddCells(S2Cap index_cap, int num_cells, S2CellIndex index)
        {
            for (int i = 0; i < num_cells; ++i)
            {
                index.Add(new S2CellId(S2Testing.SamplePoint(index_cap)), i);
            }
        }
    }

    // Generates a collection of S2Caps that are approximately within the given
    // "index_cap", generates a covering with "max_cells_per_cap" for each one,
    // and adds the coverings to the index.  The radius of each cap is chosen
    // randomly such that the total area of the coverings is approximately
    // "cap_density" times the area of "index_cap".  In other words, a random
    // point inside "index_cap" is likely to intersect about "cap_density"
    // coverings (within a factor of 2 or so).
    private struct CapsCellIndexFactory : ICellIndexFactory
    {
        public CapsCellIndexFactory(int max_cells_per_cap, double cap_density)
        {
            max_cells_per_cap_ = max_cells_per_cap;
            cap_density_ = cap_density;
        }

        public void AddCells(S2Cap index_cap, int num_cells, S2CellIndex index)
        {
            // All of this math is fairly approximate, since the coverings don't have
            // exactly the given number of cells, etc.
            int num_caps = (num_cells - 1) / max_cells_per_cap_ + 1;
            double max_area = index_cap.Area() * cap_density_ / num_caps;
            for (int i = 0; i < num_caps; ++i)
            {
                // The coverings are bigger than the caps, so we compensate for this by
                // choosing the cap area randomly up to the limit value.
                var cap = S2Cap.FromCenterArea(S2Testing.SamplePoint(index_cap),
                                                 S2Testing.Random.RandDouble() * max_area);
                S2RegionCoverer coverer = new();
                coverer.Options_.MaxCells = (max_cells_per_cap_);
                index.Add(coverer.GetCovering(cap), i);
            }
        }

        private readonly int max_cells_per_cap_;
        private readonly double cap_density_;
    }
}
