namespace S2Geometry;

using Target = S2DistanceTarget<S2MaxDistance>;

public class S2FurthestEdgeQueryTests
{
    private const int kNumIndexes = 50;
    private const int kNumEdges = 100;
    private const int kNumQueries = 200;

    // The approximate radius of S2Cap from which query edges are chosen.
    private static readonly S1Angle kTestCapRadius = S2Testing.KmToAngle(10);
    private static ITestOutputHelper _logger;

    public S2FurthestEdgeQueryTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_NoEdges()
    {
        var index = new MutableS2ShapeIndex();
        var query = new S2FurthestEdgeQuery(index);
        var target = new S2FurthestEdgeQuery.PointTarget(new S2Point(1, 0, 0));
        var edge = query.FindFurthestEdge(target);
        Assert.Equal(S1ChordAngle.Negative, edge.Distance);
        Assert.Equal(-1, edge.EdgeId);
        Assert.Equal(-1, edge.ShapeId);
        Assert.False(edge.IsInterior());
        Assert.True(edge.IsEmpty());
        Assert.Equal(S1ChordAngle.Negative, query.GetDistance(target));
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_OptionsNotModified()
    {
        // Tests that FindFurthestEdge(), GetDistance(), and IsDistanceGreater() do
        // not modify query.options(), even though all of these methods have their
        // own specific options requirements.
        S2FurthestEdgeQuery.Options options = new()
        {
            MaxResults = 3,
            MinDistance = S1ChordAngle.FromDegrees(1),
            MaxError = S1ChordAngle.FromDegrees(0.001)
        };
        var index = MakeIndexOrDie("0:1 | 0:2 | 0:3 # #");
        S2FurthestEdgeQuery query = new(index, options);
        S2FurthestEdgeQuery.PointTarget target = new(MakePointOrDie("0:4"));
        Assert.Equal(0, query.FindFurthestEdge(target).EdgeId);
        Assert2.Near(3.0, query.GetDistance(target).Degrees(), S2.DoubleError);
        Assert.True(query.IsDistanceGreater(target, S1ChordAngle.FromDegrees(1.5)));

        // Verify that none of the options above were modified.
        Assert.Equal(options.MaxResults, query.Options_.MaxResults);
        Assert.Equal(options.MinDistance, query.Options_.MinDistance);
        Assert.Equal(options.MaxError, query.Options_.MaxError);
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_DistanceEqualToLimit()
    {
        // Tests the behavior of IsDistanceGreater, IsDistanceGreaterOrEqual, and
        // IsConservativeDistanceGreaterOrEqual (and the corresponding Options) when
        // the distance to the target exactly equals the chosen limit.
        var p0 = MakePointOrDie("23:12");
        S2Point p1 = MakePointOrDie("47:11");
        S2Point[] index_points = [p0];
        MutableS2ShapeIndex index = [new S2PointVectorShape(index_points)];
        S2FurthestEdgeQuery query = new(index);

        // Start with antipodal points and a maximum (180 degrees) distance.
        S2FurthestEdgeQuery.PointTarget target0 = new(-p0);
        S1ChordAngle dist_max = S1ChordAngle.Straight;
        Assert.False(query.IsDistanceGreater(target0, dist_max));
        Assert.True(query.IsDistanceGreaterOrEqual(target0, dist_max));
        Assert.True(query.IsConservativeDistanceGreaterOrEqual(target0, dist_max));

        // Now try two points separated by a non-maximal distance.
        S2FurthestEdgeQuery.PointTarget target1 = new(-p1);
        S1ChordAngle dist1 = GetMaxDistanceToEdge(p0, -p1, -p1);
        Assert.False(query.IsDistanceGreater(target1, dist1));
        Assert.True(query.IsDistanceGreaterOrEqual(target1, dist1));
        Assert.True(query.IsConservativeDistanceGreaterOrEqual(target1, dist1));
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_TrueDistanceGreaterThanS1ChordAngleDistance()
    {
        // Tests that IsConservativeDistanceGreaterOrEqual returns points where the
        // true distance is slightly greater than the one computed by S1ChordAngle.
        //
        // The points below had the worst error from among 1x10^6 random pairs.
        S2Point p0 = new(0.72362949088190598, -0.39019820403414807, -0.56930283812266336);
        S2Point p1 = new(0.54383822931548842, 0.758981734255934404, 0.35803171284238039);

        // The S1ChordAngle distance is ~3 ulps greater than the true distance.
        S1ChordAngle dist1 = GetMaxDistanceToEdge(p0, p1, p1);
        var limit = dist1.Successor().Successor().Successor();
        Assert.True(S2Pred.CompareDistance(p0, p1, limit) > 0);

        // Verify that IsConservativeDistanceGreaterOrEqual() still returns "p1".
        S2Point[] index_points = [p0];
        MutableS2ShapeIndex index = [new S2PointVectorShape(index_points)];
        S2FurthestEdgeQuery query = new(index);
        S2FurthestEdgeQuery.PointTarget target1 = new(p1);
        Assert.False(query.IsDistanceGreater(target1, limit));
        Assert.False(query.IsDistanceGreaterOrEqual(target1, limit));
        Assert.True(query.IsConservativeDistanceGreaterOrEqual(target1, limit));
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_AntipodalPointInsideIndexedPolygon()
    {
        // Tests a target point antipodal to the interior of an indexed polygon.
        // (The index also includes a polyline loop with no interior.)
        var index = MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
        S2FurthestEdgeQuery.Options options = new()
        {
            // First check that with include_interiors set to true, the distance is 180.
            IncludeInteriors = true,
            MinDistance = new(S1Angle.FromDegrees(178))
        };
        S2FurthestEdgeQuery query = new(index, options);
        S2FurthestEdgeQuery.PointTarget target = new(-MakePointOrDie("2:12"));
        var results = query.FindFurthestEdges(target);
        Assert.True(results.Count > 0);
        Assert.Equal(S1ChordAngle.Straight, results[0].Distance);
        // Should find the polygon shape (id = 1).
        Assert.Equal(1, results[0].ShapeId);
        // Should find the interior, so no specific edge id.
        Assert.Equal(-1, results[0].EdgeId);
        Assert.True(results[0].IsInterior());
        Assert.False(results[0].IsEmpty());

        // Next check that with include_interiors set to false, the distance is less
        // than 180 for the same target and index.
        query.Options_.IncludeInteriors = false;
        results = query.FindFurthestEdges(target);
        Assert.True(results.Count > 0);
        Assert.True(results[0].Distance <= S1ChordAngle.Straight);
        Assert.Equal(1, results[0].ShapeId);
        // Found a specific edge, so id should be positive.
        Assert.Equal(3, results[0].EdgeId);
        Assert.False(results[0].IsInterior());
        Assert.False(results[0].IsEmpty());
        S2Shape.Edge e0 = query.GetEdge(results[0]);
        Assert.True(S2.ApproxEquals(e0.V0, S2LatLng.FromDegrees(5, 10).ToPoint()));
            //<< S2LatLng(e0.v0);
        Assert.True(S2.ApproxEquals(e0.V1, S2LatLng.FromDegrees(0, 10).ToPoint()));
            //<< S2LatLng(e0.v1);
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_AntipodalPointOutsideIndexedPolygon()
    {
        // Tests a target point antipodal to the interior of a polyline loop with no
        // interior.  The index also includes a polygon almost antipodal to the
        // target, but with all edges closer than the min_distance threshold.
        var index = MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
        S2FurthestEdgeQuery.Options options = new()
        {
            IncludeInteriors = true,
            MinDistance = new(S1Angle.FromDegrees(179))
        };
        S2FurthestEdgeQuery query = new(index, options);
        S2FurthestEdgeQuery.PointTarget target = new(-MakePointOrDie("2:2"));
        var results = query.FindFurthestEdges(target);
        Assert.Empty(results);
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_TargetPolygonContainingIndexedPoints()
    {
        // Two points are contained within a polyline loop (no interior) and two
        // points are contained within a polygon.
        var index = MakeIndexOrDie("2:2 | 4:4 | 1:11 | 3:12 # #");
        S2FurthestEdgeQuery query = new(index);
        query.Options_.UseBruteForce = false;
        var target_index = MakeIndexOrDie(
            "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
        S2FurthestEdgeQuery.ShapeIndexTarget target = new(target_index)
        {
            IncludeInteriors = true,
            UseBruteForce = true
        };
        var results1 = query.FindFurthestEdges(target);
        // All points should be returned since we did not specify max_results.
        Assert.Equal(4, results1.Count);
        Assert.NotEqual(S1ChordAngle.Zero, results1[0].Distance);
        Assert.Equal(0, results1[0].ShapeId);
        Assert.Equal(0, results1[0].EdgeId);  // 2:2 (to 5:15)
        Assert.NotEqual(S1ChordAngle.Zero, results1[1].Distance);
        Assert.Equal(0, results1[1].ShapeId);
        Assert.Equal(3, results1[1].EdgeId);  // 3:12 (to 0:0)
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_AntipodalPolygonContainingIndexedPoints()
    {
        // Two antipodal points are contained within a polyline loop (no interior)
        // and two antipodal points are contained within a polygon.
        var antipodal_points = (
            from p in ParsePointsOrDie("2:2, 3:3, 1:11, 3:13")
            select -p).ToArray();
        MutableS2ShapeIndex index = [new S2PointVectorShape(antipodal_points)];

        S2FurthestEdgeQuery query = new(index);
        query.Options_.MinDistance = new(S1Angle.FromDegrees(179));
        var target_index = MakeIndexOrDie(
            "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
        S2FurthestEdgeQuery.ShapeIndexTarget target = new(target_index)
        {
            IncludeInteriors = true
        };
        var results = query.FindFurthestEdges(target);
        Assert.Equal(2, results.Count);
        Assert.Equal(S1ChordAngle.Straight, results[0].Distance);
        Assert.Equal(0, results[0].ShapeId);
        Assert.Equal(2, results[0].EdgeId);  // 1:11
        Assert.Equal(S1ChordAngle.Straight, results[1].Distance);
        Assert.Equal(0, results[1].ShapeId);
        Assert.Equal(3, results[1].EdgeId);  // 3:13
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_EmptyPolygonTarget()
    {
        // Verifies that distances are measured correctly to empty polygon targets.
        var empty_polygon_index = MakeIndexOrDie("# # empty");
        var point_index = MakeIndexOrDie("1:1 # #");
        var full_polygon_index = MakeIndexOrDie("# # full");
        S2FurthestEdgeQuery.ShapeIndexTarget target = new(empty_polygon_index)
        {
            IncludeInteriors = true
        };

        S2FurthestEdgeQuery empty_query = new(empty_polygon_index);
        empty_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Negative, empty_query.GetDistance(target));

        S2FurthestEdgeQuery point_query = new(point_index);
        point_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Negative, point_query.GetDistance(target));

        S2FurthestEdgeQuery full_query = new(full_polygon_index);
        full_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Negative, full_query.GetDistance(target));
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_FullLaxPolygonTarget()
    {
        // Verifies that distances are measured correctly to full LaxPolygon targets.
        var empty_polygon_index = MakeIndexOrDie("# # empty");
        var point_index = MakeIndexOrDie("1:1 # #");
        var full_polygon_index = MakeIndexOrDie("# # full");
        S2FurthestEdgeQuery.ShapeIndexTarget target = new(full_polygon_index)
        {
            IncludeInteriors = true
        };

        S2FurthestEdgeQuery empty_query = new(empty_polygon_index);
        empty_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Negative, empty_query.GetDistance(target));

        S2FurthestEdgeQuery point_query = new(point_index);
        point_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Straight, point_query.GetDistance(target));

        S2FurthestEdgeQuery full_query = new(full_polygon_index);
        full_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Straight, full_query.GetDistance(target));
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_FullS2PolygonTarget()
    {
        // Verifies that distances are measured correctly to full S2Polygon targets
        // (which use a different representation of "full" than LaxPolygon does).
        var empty_polygon_index = MakeIndexOrDie("# # empty");
        var point_index = MakeIndexOrDie("1:1 # #");
        var full_polygon_index = MakeIndexOrDie("# #");
        full_polygon_index.Add(new S2Polygon.OwningShape(
            MakePolygonOrDie("full")));

        S2FurthestEdgeQuery.ShapeIndexTarget target = new(full_polygon_index)
        {
            IncludeInteriors = true
        };

        S2FurthestEdgeQuery empty_query = new(empty_polygon_index);
        empty_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Negative, empty_query.GetDistance(target));

        S2FurthestEdgeQuery point_query = new(point_index);
        point_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Straight, point_query.GetDistance(target));

        S2FurthestEdgeQuery full_query = new(full_polygon_index);
        full_query.Options_.IncludeInteriors = true;
        Assert.Equal(S1ChordAngle.Straight, full_query.GetDistance(target));
    }

    //////////////////////////////////////////////////////////////////////////////
    //  General query testing by comparing with brute force method.
    //////////////////////////////////////////////////////////////////////////////

    [Fact]
    internal void Test_S2FurthestEdgeQuery_CircleEdges()
    {
        TestWithIndexFactory(new RegularLoopShapeIndexFactory(),
                             kNumIndexes, kNumEdges, kNumQueries);
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_FractalEdges()
    {
        TestWithIndexFactory(new FractalLoopShapeIndexFactory(),
                             kNumIndexes, kNumEdges, kNumQueries);
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_PointCloudEdges()
    {
        TestWithIndexFactory(new PointCloudShapeIndexFactory(),
                             kNumIndexes, kNumEdges, kNumQueries);
    }

    [Fact]
    internal void Test_S2FurthestEdgeQuery_OptionsS1AngleSetters()
    {
        // Verify that the S1Angle and S1ChordAngle versions do the same thing.
        // This is mainly to prevent the (so far unused) S1Angle versions from
        // being detected as dead code.
        S2FurthestEdgeQuery.Options angle_options = new()
        {
            MinDistance = new(S1Angle.FromDegrees(1))
        };
        S2FurthestEdgeQuery.Options chord_angle_options = new()
        {
            MinDistance = S1ChordAngle.FromDegrees(1)
        };
        Assert.Equal(chord_angle_options.MinDistance, angle_options.MinDistance);

        angle_options.InclusiveMinDistance = new(S1Angle.FromDegrees(1));
        chord_angle_options.InclusiveMinDistance = S1ChordAngle.FromDegrees(1);
        Assert.Equal(chord_angle_options.MinDistance, angle_options.MinDistance);

        angle_options.ConservativeMinDistance = new(S1Angle.FromDegrees(1));
        chord_angle_options.ConservativeMinDistance = S1ChordAngle.FromDegrees(1);
        Assert.Equal(chord_angle_options.MinDistance, angle_options.MinDistance);
    }

    // In furthest edge queries, the following distance computation is used when
    // updating max distances.
    private static S1ChordAngle GetMaxDistanceToEdge(S2Point x, S2Point y0, S2Point y1)
    {
        S1ChordAngle dist = S1ChordAngle.Negative;
        S2.UpdateMaxDistance(x, y0, y1, ref dist);
        return dist;
    }


    // Converts to the format required by CheckDistanceResults() in S2Testing.h
    // TODO(user): When S2ClosestEdgeQuery.Result is made into a class, some
    // of the following code may become redundant with that in
    // s2closest_edge_query.cc.
    private static List<(S2MaxDistance, Edge)> ConvertResults(List<S2FurthestEdgeQuery.Result> edges)
    {
        return edges.Select(edge => (new S2MaxDistance(edge.Distance),
                        new Edge(edge.ShapeId, edge.EdgeId))).ToList();
    }

    // Use "query" to find the furthest edge(s) to the given target.  Verify that
    // the results satisfy the search criteria.
    private static void GetFurthestEdges(Target target, S2FurthestEdgeQuery query, List<S2FurthestEdgeQuery.Result> edges)
    {
        query.FindFurthestEdges(target, edges);
        Assert.True(edges.Count <= query.Options_.MaxResults);
        if (query.Options_.MinDistance == S1ChordAngle.Negative)
        {
            int min_expected = Math.Min(query.Options_.MaxResults,
                                   query.Index().GetCountEdges());
            if (!query.Options_.IncludeInteriors)
            {
                // We can predict exactly how many edges should be returned.
                Assert.Equal(min_expected, edges.Count);
            }
            else
            {
                // All edges should be returned, and possibly some shape interiors.
                Assert.True(min_expected <= edges.Count);
            }
        }
        foreach (var edge in edges)
        {
            // Check that the edge satisfies the min_distance() condition.
            Assert.True(edge.Distance >= query.Options_.MinDistance);
        }
    }

    private static void TestFindFurthestEdges(Target target, S2FurthestEdgeQuery query)
    {
        var expected = new List<S2FurthestEdgeQuery.Result>();
        var actual = new List<S2FurthestEdgeQuery.Result>();

        query.Options_.UseBruteForce = true;
        GetFurthestEdges(target, query, expected);
        query.Options_.UseBruteForce = false;
        GetFurthestEdges(target, query, actual);

        S1ChordAngle min_distance = query.Options_.MinDistance;
        S1ChordAngle max_error = query.Options_.MaxError;
        Assert.True(S2TestingCheckDistance<Edge, S2MaxDistance>.CheckDistanceResults(
            ConvertResults(expected),
            ConvertResults(actual),
            query.Options_.MaxResults,
            new S2MaxDistance(min_distance),
            max_error,
            _logger.WriteLine));

        if (expected.Count==0)
        {
            return;
        }

        // Note that when options.max_error() > 0, expected[0].distance may not be
        // the maximum distance.  It is never smaller by more than max_error(), but
        // the actual value also depends on max_results().
        //
        // Here we verify that GetDistance() and IsDistanceGreater() return results
        // that are consistent with the max_error() setting.
        S1ChordAngle expected_distance = expected[0].Distance;
        Assert.True(query.GetDistance(target) >= expected_distance - max_error);

        // Test IsDistanceGreater().
        Assert.False(query.IsDistanceGreater(
            target, expected_distance + max_error));
        Assert.True(query.IsDistanceGreater(
            target, expected[0].Distance.Predecessor()));
    }

    // The running time of this test is proportional to
    //    (num_indexes + num_queries) * num_edges.
    // (Note that every query is checked using the brute force algorithm.)
    private static void TestWithIndexFactory(IShapeIndexFactory factory,
                             int num_indexes, int num_edges,
                             int num_queries)
    {
        // Build a set of MutableS2ShapeIndexes containing the desired geometry.
        List<S2Cap> index_caps = [];
        List<MutableS2ShapeIndex> indexes = [];
        for (int i = 0; i < num_indexes; ++i)
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
            index_caps.Add(new S2Cap(S2Testing.RandomPoint(), kTestCapRadius));
            indexes.Add([]);
            factory.AddEdges(index_caps.Last(), num_edges, indexes.Last());
        }

        for (int i = 0; i < num_queries; ++i)
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
            int i_index = S2Testing.Random.Uniform(num_indexes);
            var index_cap = index_caps[i_index];

            // Choose query points from an area approximately 4x larger than the
            // geometry being tested.
            S1Angle query_radius = 2 * index_cap.RadiusAngle();
            S2FurthestEdgeQuery query = new(indexes[i_index]);

            // Exercise the opposite-hemisphere code 1/5 of the time.
            int antipodal = S2Testing.Random.OneIn(5) ? -1 : 1;
            S2Cap query_cap = new(antipodal * index_cap.Center, query_radius);

            // Occasionally we don't set any limit on the number of result edges.
            // (This may return all edges if we also don't set a distance limit.)
            if (!S2Testing.Random.OneIn(5))
            {
                query.Options_.MaxResults = 1 + S2Testing.Random.Uniform(10);
            }
            // We set a distance limit 2/3 of the time.
            if (!S2Testing.Random.OneIn(3))
            {
                query.Options_.MinDistance = new(S2Testing.Random.RandDouble() * query_radius);
            }
            if (S2Testing.Random.OneIn(2))
            {
                // Choose a maximum error whose logarithm is uniformly distributed over
                // a reasonable range, except that it is sometimes zero.
                query.Options_.MaxError = new(S1Angle.FromRadians(
                    Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius.Radians));
            }
            query.Options_.IncludeInteriors = S2Testing.Random.OneIn(2);
            int target_type = S2Testing.Random.Uniform(4);
            if (target_type == 0)
            {
                // Find the edges furthest from a given point.
                S2Point point = S2Testing.SamplePoint(query_cap);
                S2FurthestEdgeQuery.PointTarget target = new(point);
                TestFindFurthestEdges(target, query);
            }
            else if (target_type == 1)
            {
                // Find the edges furthest from a given edge.
                S2Point a = S2Testing.SamplePoint(query_cap);
                S2Point b = S2Testing.SamplePoint(
                    new S2Cap(a, Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius));
                var target = new S2FurthestEdgeQuery.EdgeTarget(a, b);
                TestFindFurthestEdges(target, query);
            }
            else if (target_type == 2)
            {
                // Find the edges furthest from a given cell.
                int min_level = S2.kMaxDiag.GetLevelForMaxValue(query_radius.Radians);
                int level = min_level + S2Testing.Random.Uniform(
                    S2.kMaxCellLevel - min_level + 1);
                S2Point a = S2Testing.SamplePoint(query_cap);
                S2Cell cell = new(new S2CellId(a).Parent(level));
                S2FurthestEdgeQuery.CellTarget target = new(cell);
                TestFindFurthestEdges(target, query);
            }
            else
            {
                Assert.Equal(3, target_type);
                // Use another one of the pre-built indexes as the target.
                int j_index = S2Testing.Random.Uniform(num_indexes);
                S2FurthestEdgeQuery.ShapeIndexTarget target = new(indexes[j_index])
                {
                    IncludeInteriors = S2Testing.Random.OneIn(2)
                };
                TestFindFurthestEdges(target, query);
            }
        }
    }
}
