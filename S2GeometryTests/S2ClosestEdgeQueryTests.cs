using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace S2Geometry
{
    using S2ShapeUtil;
    using Distance = S1ChordAngle; // S2MinDistance
    using Result = S2ClosestEdgeQueryBase<S1ChordAngle>.Result;
    using Target = S2DistanceTarget<S1ChordAngle>; // S2MinDistanceTarget

    using TestingResult = ValueTuple<S1ChordAngle, S2ShapeUtil.Edge>;
    using TestingDistance = S2TestingCheckDistance<S2ShapeUtil.Edge, S1ChordAngle>;

    public class S2ClosestEdgeQueryTests
    {
        // The approximate radius of S2Cap from which query edges are chosen.
        private static readonly S1Angle kTestCapRadius = S2Testing.KmToAngle(10);
        // An approximate bound on the distance measurement error for "reasonable"
        // distances (say, less than Pi/2) due to using S1ChordAngle.
        private static readonly double kTestChordAngleError = S2Constants.DoubleError;
        private const int kNumIndexes = 50;
        private const int kNumEdges = 100;
        private const int kNumQueries = 200;

        [Fact]
        public void Test_S2ClosestEdgeQuery_NoEdges()
        {
            MutableS2ShapeIndex index = new();
            S2ClosestEdgeQuery query = new(index);
            S2ClosestEdgeQuery.PointTarget target = new(new S2Point(1, 0, 0));
            var edge = query.FindClosestEdge(target);
            Assert.Equal(S1ChordAngle.Infinity, edge.Distance);
            Assert.Equal(-1, edge.ShapeId);
            Assert.Equal(-1, edge.EdgeId);
            Assert.False(edge.IsInterior);
            Assert.True(edge.IsEmpty);
            Assert.Equal(S1ChordAngle.Infinity, query.GetDistance(target));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_OptionsNotModified()
        {
            // Tests that FindClosestEdge(), GetDistance(), and IsDistanceLess() do not
            // modify query.Options_, even though all of these methods have their own
            // specific options requirements.
            S2ClosestEdgeQuery.Options options = new();
            options.MaxResults = (3);
            options.MaxDistance = (S1ChordAngle.FromDegrees(3));
            options.MaxError = (S1ChordAngle.FromDegrees(0.001));
            var index = S2TextFormat.MakeIndexOrDie("1:1 | 1:2 | 1:3 # #");
            S2ClosestEdgeQuery query = new(index, options);
            S2ClosestEdgeQuery.PointTarget target = new(S2TextFormat.MakePointOrDie("2:2"));
            Assert.Equal(1, query.FindClosestEdge(target).EdgeId);
            Assert2.Near(1.0, query.GetDistance(target).Degrees, S2Constants.DoubleError);
            Assert.True(query.IsDistanceLess(target, S1ChordAngle.FromDegrees(1.5)));

            // Verify that none of the options above were modified.
            Assert.Equal(options.MaxResults, query.Options_.MaxResults);
            Assert.Equal(options.MaxDistance, query.Options_.MaxDistance);
            Assert.Equal(options.MaxError, query.Options_.MaxError);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_DistanceEqualToLimit()
        {
            // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
            // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
            // the distance to the target exactly equals the chosen limit.
            S2Point p0 = (S2TextFormat.MakePointOrDie("23:12")), p1 = (S2TextFormat.MakePointOrDie("47:11"));
            var index_points = new[] { p0 };
            MutableS2ShapeIndex index = new();
            index.Add(new S2PointVectorShape(index_points));
            S2ClosestEdgeQuery query = new(index);

            // Start with two identical points and a zero distance.
            S2ClosestEdgeQuery.PointTarget target0 = new(p0);
            S1ChordAngle dist0 = S1ChordAngle.Zero;
            Assert.False(query.IsDistanceLess(target0, dist0));
            Assert.True(query.IsDistanceLessOrEqual(target0, dist0));
            Assert.True(query.IsConservativeDistanceLessOrEqual(target0, dist0));

            // Now try two points separated by a non-zero distance.
            S2ClosestEdgeQuery.PointTarget target1 = new(p1);
            S1ChordAngle dist1 = new(p0, p1);
            Assert.False(query.IsDistanceLess(target1, dist1));
            Assert.True(query.IsDistanceLessOrEqual(target1, dist1));
            Assert.True(query.IsConservativeDistanceLessOrEqual(target1, dist1));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_TrueDistanceLessThanS1ChordAngleDistance()
        {
            // Tests that IsConservativeDistanceLessOrEqual returns points where the
            // true distance is slightly less than the one computed by S1ChordAngle.
            //
            // The points below had the worst error from among 100,000 random pairs.
            S2Point p0 = new(0.78516762584829192, -0.50200400690845970, -0.36263449417782678);
            S2Point p1 = new(0.78563011732429433, -0.50187655940493503, -0.36180828883938054);

            // The S1ChordAngle distance is ~4 ulps greater than the true distance.
            S1ChordAngle dist1 = new(p0, p1);
            var limit = dist1.Predecessor.Predecessor.Predecessor.Predecessor;
            Assert.True(S2Pred.CompareDistance(p0, p1, limit) < 0);

            // Verify that IsConservativeDistanceLessOrEqual() still returns "p1".
            S2Point[] index_points = new[] { p0 };
            MutableS2ShapeIndex index = new();
            index.Add(new S2PointVectorShape(index_points));
            S2ClosestEdgeQuery query = new(index);
            S2ClosestEdgeQuery.PointTarget target1 = new(p1);
            Assert.False(query.IsDistanceLess(target1, limit));
            Assert.False(query.IsDistanceLessOrEqual(target1, limit));
            Assert.True(query.IsConservativeDistanceLessOrEqual(target1, limit));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_TestReuseOfQuery()
        {
            // Tests that between queries, the internal mechanism for de-duplicating
            // results is re-set.  See b/71646017.
            var index = S2TextFormat.MakeIndexOrDie("2:2 # #");
            S2ClosestEdgeQuery query = new(index);
            query.Options_.MaxError = new(S1Angle.FromDegrees(1));
            var target_index = S2TextFormat.MakeIndexOrDie("## 0:0, 0:5, 5:5, 5:0");
            S2ClosestEdgeQuery.ShapeIndexTarget target = new(target_index);
            var results1 = query.FindClosestEdges(target);
            var results2 = query.FindClosestEdges(target);
            Assert.Equal(results1.Count, results2.Count);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_TargetPointInsideIndexedPolygon()
        {
            // Tests a target point in the interior of an indexed polygon.
            // (The index also includes a polyline loop with no interior.)
            var index = S2TextFormat.MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
            S2ClosestEdgeQuery.Options options = new();
            options.IncludeInteriors = (true);
            options.MaxDistance = new(S1Angle.FromDegrees(1));
            S2ClosestEdgeQuery query = new(index, options);
            S2ClosestEdgeQuery.PointTarget target = new(S2TextFormat.MakePointOrDie("2:12"));
            var results = query.FindClosestEdges(target);
            Assert.Single(results);
            Assert.Equal(S1ChordAngle.Zero, results[0].Distance);
            Assert.Equal(1, results[0].ShapeId);
            Assert.Equal(-1, results[0].EdgeId);
            Assert.True(results[0].IsInterior);
            Assert.False(results[0].IsEmpty);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_TargetPointOutsideIndexedPolygon()
        {
            // Tests a target point in the interior of a polyline loop with no
            // interior.  (The index also includes a nearby polygon.)
            var index = S2TextFormat.MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
            S2ClosestEdgeQuery.Options options = new();
            options.IncludeInteriors = (true);
            options.MaxDistance = new(S1Angle.FromDegrees(1));
            S2ClosestEdgeQuery query = new(index, options);
            S2ClosestEdgeQuery.PointTarget target = new(S2TextFormat.MakePointOrDie("2:2"));
            var results = query.FindClosestEdges(target);
            Assert.Empty(results);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_TargetPolygonContainingIndexedPoints()
        {
            // Two points are contained within a polyline loop (no interior) and two
            // points are contained within a polygon.
            var index = S2TextFormat.MakeIndexOrDie("2:2 | 3:3 | 1:11 | 3:13 # #");
            S2ClosestEdgeQuery query = new(index);
            query.Options_.MaxDistance = new(S1Angle.FromDegrees(1));
            var target_index = S2TextFormat.MakeIndexOrDie(
                "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
            S2ClosestEdgeQuery.ShapeIndexTarget target = new(target_index);
            target.IncludeInteriors = (true);
            var results = query.FindClosestEdges(target);
            Assert.Equal(2, results.Count);
            Assert.Equal(S1ChordAngle.Zero, results[0].Distance);
            Assert.Equal(0, results[0].ShapeId);
            Assert.Equal(2, results[0].EdgeId);  // 1:11
            Assert.Equal(S1ChordAngle.Zero, results[1].Distance);
            Assert.Equal(0, results[1].ShapeId);
            Assert.Equal(3, results[1].EdgeId);  // 3:13
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_EmptyTargetOptimized()
        {
            // Ensure that the optimized algorithm handles empty targets when a distance
            // limit is specified.
            MutableS2ShapeIndex index = new();
            index.Add(new S2Polygon.OwningShape(new S2Polygon(
                S2Loop.MakeRegularLoop(new S2Point(1, 0, 0), S1Angle.FromRadians(0.1), 1000))));
            S2ClosestEdgeQuery query = new(index);
            query.Options_.MaxDistance = new(S1Angle.FromRadians(1e-5));
            MutableS2ShapeIndex target_index = new();
            S2ClosestEdgeQuery.ShapeIndexTarget target = new(target_index);
            Assert.Empty(query.FindClosestEdges(target));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_EmptyPolygonTarget()
        {
            // Verifies that distances are measured correctly to empty polygon targets.
            var empty_polygon_index = S2TextFormat.MakeIndexOrDie("# # empty");
            var point_index = S2TextFormat.MakeIndexOrDie("1:1 # #");
            var full_polygon_index = S2TextFormat.MakeIndexOrDie("# # full");
            S2ClosestEdgeQuery.ShapeIndexTarget target = new(empty_polygon_index);
            target.IncludeInteriors = (true);

            S2ClosestEdgeQuery empty_query = new(empty_polygon_index);
            empty_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Infinity, empty_query.GetDistance(target));

            S2ClosestEdgeQuery point_query = new(point_index);
            point_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Infinity, point_query.GetDistance(target));

            S2ClosestEdgeQuery full_query = new(full_polygon_index);
            full_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Infinity, full_query.GetDistance(target));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_FullLaxPolygonTarget()
        {
            // Verifies that distances are measured correctly to full LaxPolygon targets.
            var empty_polygon_index = S2TextFormat.MakeIndexOrDie("# # empty");
            var point_index = S2TextFormat.MakeIndexOrDie("1:1 # #");
            var full_polygon_index = S2TextFormat.MakeIndexOrDie("# # full");
            S2ClosestEdgeQuery.ShapeIndexTarget target = new(full_polygon_index);
            target.IncludeInteriors = (true);

            S2ClosestEdgeQuery empty_query = new(empty_polygon_index);
            empty_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Infinity, empty_query.GetDistance(target));

            S2ClosestEdgeQuery point_query = new(point_index);
            point_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Zero, point_query.GetDistance(target));

            S2ClosestEdgeQuery full_query = new(full_polygon_index);
            full_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Zero, full_query.GetDistance(target));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_FullS2PolygonTarget()
        {
            // Verifies that distances are measured correctly to full S2Polygon targets
            // (which use a different representation of "full" than LaxPolygon does).
            var empty_polygon_index = S2TextFormat.MakeIndexOrDie("# # empty");
            var point_index = S2TextFormat.MakeIndexOrDie("1:1 # #");
            var full_polygon_index = S2TextFormat.MakeIndexOrDie("# #");
            full_polygon_index.Add(new S2Polygon.OwningShape(
                S2TextFormat.MakePolygonOrDie("full")));

            S2ClosestEdgeQuery.ShapeIndexTarget target = new(full_polygon_index);
            target.IncludeInteriors = (true);

            S2ClosestEdgeQuery empty_query = new(empty_polygon_index);
            empty_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Infinity, empty_query.GetDistance(target));

            S2ClosestEdgeQuery point_query = new(point_index);
            point_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Zero, point_query.GetDistance(target));

            S2ClosestEdgeQuery full_query = new(full_polygon_index);
            full_query.Options_.IncludeInteriors = (true);
            Assert.Equal(S1ChordAngle.Zero, full_query.GetDistance(target));
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_IsConservativeDistanceLessOrEqual()
        {
            // Test
            int num_tested = 0;
            int num_conservative_needed = 0;
            for (int iter = 0; iter < 1000; ++iter)
            {
                S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
                S2Point x = S2Testing.RandomPoint();
                S2Point dir = S2Testing.RandomPoint();
                S1Angle r = S1Angle.FromRadians(Math.PI * Math.Pow(1e-30, S2Testing.Random.RandDouble()));
                S2Point y = S2EdgeDistances.InterpolateAtDistance(r, x, dir);
                S1ChordAngle limit = new(r);
                if (S2Pred.CompareDistance(x, y, limit) <= 0)
                {
                    MutableS2ShapeIndex index = new();
                    index.Add(new S2PointVectorShape(new S2Point[] { x }));
                    S2ClosestEdgeQuery query = new(index);
                    S2ClosestEdgeQuery.PointTarget target = new(y);
                    Assert.True(query.IsConservativeDistanceLessOrEqual(target, limit));
                    ++num_tested;
                    if (!query.IsDistanceLess(target, limit)) ++num_conservative_needed;
                }
            }
            // Verify that in most test cases, the distance between the target points
            // was close to the desired value.  Also verify that at least in some test
            // cases, the conservative distance test was actually necessary.
            Assert.True(num_tested >= 300);
            Assert.True(num_tested <= 700);
            Assert.True(num_conservative_needed >= 25);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_CircleEdges()
        {
            TestWithIndexFactory(new RegularLoopShapeIndexFactory(),
                                 kNumIndexes, kNumEdges, kNumQueries);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_FractalEdges()
        {
            TestWithIndexFactory(new FractalLoopShapeIndexFactory(),
                                 kNumIndexes, kNumEdges, kNumQueries);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_PointCloudEdges()
        {
            TestWithIndexFactory(new PointCloudShapeIndexFactory(),
                                 kNumIndexes, kNumEdges, kNumQueries);
        }

        [Fact]
        public void Test_S2ClosestEdgeQuery_ConservativeCellDistanceIsUsed()
        {
            // Don't use google.FlagSaver, so it works in opensource without gflags.
            int saved_seed = S2Testing.Random.RandomSeed;
            // These specific test cases happen to fail if max_error() is not properly
            // taken into account when measuring distances to S2ShapeIndex cells.
            foreach (int seed in new[] { 42, 681, 894, 1018, 1750, 1759, 2401 })
            {
                S2Testing.Random.RandomSeed = seed;
                TestWithIndexFactory(new FractalLoopShapeIndexFactory(),
                                     5, 100, 10);
            }
            S2Testing.Random.RandomSeed = saved_seed;
        }

        // Converts to the format required by CheckDistanceResults() in S2Testing.h.
        private static List<TestingResult> ConvertResults(List<Result> results)
        {
            List<TestingResult> testing_results = new();
            foreach (var result in results)
            {
                testing_results.Add(new(result.Distance, new(result.ShapeId, result.EdgeId)));
            }
            return testing_results;
        }

        // Use "query" to find the closest edge(s) to the given target.  Verify that
        // the results satisfy the search criteria.
        private static void GetClosestEdges(Target target,
                            S2ClosestEdgeQuery query,
                            List<Result> edges)
        {
            query.FindClosestEdges(target, edges);
            Assert.True(edges.Count <= query.Options_.MaxResults);
            if (query.Options_.MaxDistance==
                Distance.Infinity)
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
                // Check that the edge satisfies the max_distance() condition.
                Assert.True(edge.Distance < query.Options_.MaxDistance);
            }
        }

        private static Result TestFindClosestEdges(Target target, S2ClosestEdgeQuery query)
        {
            List<Result> expected = new(), actual = new();
            query.Options_.UseBruteForce = (true);
            GetClosestEdges(target, query, expected);
            query.Options_.UseBruteForce = (false);
            GetClosestEdges(target, query, actual);
            Assert.True(TestingDistance.CheckDistanceResults(ConvertResults(expected),
                    ConvertResults(actual),
                    query.Options_.MaxResults,
                    query.Options_.MaxDistance,
                    query.Options_.MaxError));

            if (!expected.Any()) return new Result();

            // Note that when options.max_error() > 0, expected[0].distance() may not
            // be the minimum distance.  It is never larger by more than max_error(),
            // but the actual value also depends on max_results().
            //
            // Here we verify that GetDistance() and IsDistanceLess() return results
            // that are consistent with the max_error() setting.
            var max_error = query.Options_.MaxError;
            var min_distance = expected[0].Distance;
            Assert.True(query.GetDistance(target) <= min_distance + max_error);

            // Test IsDistanceLess().
            Assert.False(query.IsDistanceLess(target, min_distance - max_error));
            Assert.True(query.IsConservativeDistanceLessOrEqual(target, min_distance));

            // Return the closest edge result so that we can also test Project.
            return expected[0];
        }

        // The running time of this test is proportional to
        //    (num_indexes + num_queries) * num_edges.
        // (Note that every query is checked using the brute force algorithm.)
        private static void TestWithIndexFactory(IShapeIndexFactory factory,
                                 int num_indexes, int num_edges,
                                 int num_queries)
        {
            // Build a set of MutableS2ShapeIndexes containing the desired geometry.
            List<S2Cap> index_caps = new();
            List<MutableS2ShapeIndex> indexes = new();
            for (int i = 0; i < num_indexes; ++i)
            {
                S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
                index_caps.Add(new S2Cap(S2Testing.RandomPoint(), kTestCapRadius));
                indexes.Add(new MutableS2ShapeIndex());
                factory.AddEdges(index_caps.Last(), num_edges, indexes.Last());
            }
            for (int i = 0; i < num_queries; ++i)
            {
                S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
                int i_index = S2Testing.Random.Uniform(num_indexes);
                var index_cap = index_caps[i_index];

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                var query_radius = 2 * index_cap.RadiusAngle;
                S2Cap query_cap = new(index_cap.Center, query_radius);
                S2ClosestEdgeQuery query = new(indexes[i_index]);

                // Occasionally we don't set any limit on the number of result edges.
                // (This may return all edges if we also don't set a distance limit.)
                if (!S2Testing.Random.OneIn(5))
                {
                    query.Options_.MaxResults = (1 + S2Testing.Random.Uniform(10));
                }
                // We set a distance limit 2/3 of the time.
                if (!S2Testing.Random.OneIn(3))
                {
                    query.Options_.MaxDistance = new(S2Testing.Random.RandDouble() * query_radius);
                }
                if (S2Testing.Random.OneIn(2))
                {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.Options_.MaxError = new(S1Angle.FromRadians(
                        Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius.Radians));
                }
                query.Options_.IncludeInteriors = (S2Testing.Random.OneIn(2));
                int target_type = S2Testing.Random.Uniform(4);
                if (target_type == 0)
                {
                    // Find the edges closest to a given point.
                    S2Point point = S2Testing.SamplePoint(query_cap);
                    S2ClosestEdgeQuery.PointTarget target = new(point);
                    var closest = TestFindClosestEdges(target, query);
                    if (!closest.Distance.IsInfinity)
                    {
                        // Also test the Project method.
                        Assert2.Near(
                            closest.Distance.ToAngle().Radians,
                            new S1Angle(point, query.Project(point, closest)).Radians,
                            kTestChordAngleError);
                    }
                }
                else if (target_type == 1)
                {
                    // Find the edges closest to a given edge.
                    S2Point a = S2Testing.SamplePoint(query_cap);
                    S2Point b = S2Testing.SamplePoint(
                        new S2Cap(a, Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius));
                    S2ClosestEdgeQuery.EdgeTarget target = new(a, b);
                    TestFindClosestEdges(target, query);
                }
                else if (target_type == 2)
                {
                    // Find the edges closest to a given cell.
                    int min_level = S2Metrics.kMaxDiag.GetLevelForMaxValue(query_radius.Radians);
                    int level = min_level + S2Testing.Random.Uniform(
                        S2Constants.kMaxCellLevel - min_level + 1);
                    S2Point a = S2Testing.SamplePoint(query_cap);
                    S2Cell cell = new(new S2CellId(a).Parent(level));
                    S2ClosestEdgeQuery.CellTarget target = new(cell);
                    TestFindClosestEdges(target, query);
                }
                else
                {
                    Assert.Equal(3, target_type);
                    // Use another one of the pre-built indexes as the target.
                    int j_index = S2Testing.Random.Uniform(num_indexes);
                    S2ClosestEdgeQuery.ShapeIndexTarget target = new(indexes[j_index]);
                    target.IncludeInteriors = (S2Testing.Random.OneIn(2));
                    TestFindClosestEdges(target, query);
                }
            }
        }
    }
}
