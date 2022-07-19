using System.Collections.Generic;
using System.Linq;
using System;
using Xunit;
using TestIndex = S2Geometry.S2PointIndex<int>;
using TestQuery = S2Geometry.S2ClosestPointQuery<int>;
// The result format required by CheckDistanceResults() in S2Testing.h.
using Target = S2Geometry.S2DistanceTarget<S2Geometry.S1ChordAngle>;

namespace S2Geometry
{
    public class S2ClosestPointQueryTests
    {
        // The approximate radius of S2Cap from which query points are chosen.
        private static readonly S1Angle kTestCapRadius = S2Testing.KmToAngle(10);
        private const int kNumIndexes = 10;
        private const int kNumPoints = 1000;
        private const int kNumQueries = 50;

        [Fact]
        public void Test_S2ClosestPointQuery_NoPoints() {
            TestIndex index=new();
            TestQuery query = new(index);
            TestQuery.PointTarget target = new(new S2Point(1, 0, 0));
            var results = query.FindClosestPoints(target);
            Assert.Empty(results);
        }

        [Fact]
        public void Test_S2ClosestPointQuery_ManyDuplicatePoints() {
            const int kNumPoints = 10000;
            S2Point kTestPoint=new(1, 0, 0);
            TestIndex index=new();
            for (int i = 0; i < kNumPoints; ++i) {
                index.Add(kTestPoint, i);
            }
            TestQuery query=new(index);
            TestQuery.PointTarget target=new(kTestPoint);
            var results = query.FindClosestPoints(target);
            Assert.Equal(kNumPoints, results.Count);
        }

        [Fact]
        public void Test_S2ClosestPointQuery_EmptyTargetOptimized() {
            // Ensure that the optimized algorithm handles empty targets when a distance
            // limit is specified.
            TestIndex index=new();
            for (int i = 0; i < 1000; ++i) {
                index.Add(S2Testing.RandomPoint(), i);
            }
            TestQuery query=new(index);
            query.Options_.MaxDistance = new S1ChordAngle(S1Angle.FromRadians(1e-5));
            MutableS2ShapeIndex target_index=new();
            TestQuery.ShapeIndexTarget target=new(target_index);
            Assert.Empty(query.FindClosestPoints(target));
        }

        [Fact]
        public void Test_S2ClosestPointQueryTest_CirclePoints() {
            TestWithIndexFactory(new CirclePointIndexFactory(),
                                 kNumIndexes, kNumPoints, kNumQueries);
        }

        [Fact]
        public void Test_S2ClosestPointQueryTest_FractalPoints() {
            TestWithIndexFactory(new FractalPointIndexFactory(),
                                 kNumIndexes, kNumPoints, kNumQueries);
        }

        [Fact]
        public void Test_S2ClosestPointQueryTest_GridPoints() {
            TestWithIndexFactory(new GridPointIndexFactory(),
                                 kNumIndexes, kNumPoints, kNumQueries);
        }

        [Fact]
        public void Test_S2ClosestPointQueryTest_ConservativeCellDistanceIsUsed() {
            int saved_seed = S2Testing.Random.RandomSeed;
            // These specific test cases happen to fail if max_error() is not properly
            // taken into account when measuring distances to S2PointIndex cells.  They
            // all involve S2ShapeIndexTarget, which takes advantage of max_error() to
            // optimize its distance calculation.
            foreach (int seed in new[] { 16, 586, 589, 822, 1959, 2298, 3155, 3490, 3723, 4953 }) {
                S2Testing.Random.RandomSeed = seed;
                TestWithIndexFactory(new FractalPointIndexFactory(), 5, 100, 10);
            }
            S2Testing.Random.RandomSeed = saved_seed;
        }

        // An abstract class that adds points to an S2PointIndex for benchmarking.
        private interface IPointIndexFactory {

            // Requests that approximately "num_points" points located within the given
            // S2Cap bound should be added to "index".
            void AddPoints(S2Cap index_cap, int num_points, TestIndex index);
        }

        // Generates points that are regularly spaced along a circle.  Points along a
        // circle are nearly the worst case for distance calculations, since many
        // points are nearly equidistant from any query point that is not immediately
        // adjacent to the circle.
        private class CirclePointIndexFactory : IPointIndexFactory {
            public void AddPoints(S2Cap index_cap, int num_points, TestIndex index) {
                var points = S2Testing.MakeRegularPoints(
                    index_cap.Center, index_cap.RadiusAngle(), num_points);
                for (int i = 0; i < points.Length; ++i) {
                    index.Add(points[i], i);
                }
            }
        };

        // Generates the vertices of a fractal whose convex hull approximately
        // matches the given cap.
        private class FractalPointIndexFactory : IPointIndexFactory {
            public void AddPoints(S2Cap index_cap, int num_points, TestIndex index) {
                S2Testing.Fractal fractal = new();
                fractal.SetLevelForApproxMaxEdges(num_points);
                fractal.FractalDimension = (1.5);
                var loop = (
                    fractal.MakeLoop(S2Testing.GetRandomFrameAt(index_cap.Center),
                                     index_cap.RadiusAngle()));
                for (int i = 0; i < loop.NumVertices; ++i) {
                    index.Add(loop.Vertex(i), i);
                }
            }
        }

        // Generates vertices on a square grid that includes the entire given cap.
        private class GridPointIndexFactory : IPointIndexFactory {
            public void AddPoints(S2Cap index_cap, int num_points, TestIndex index) {
                int sqrt_num_points = (int)Math.Ceiling(Math.Sqrt(num_points));
                S2PointS2Point frame = S2Testing.GetRandomFrameAt(index_cap.Center);
                double radius = index_cap.RadiusAngle().Radians;
                double spacing = 2 * radius / sqrt_num_points;
                for (int i = 0; i < sqrt_num_points; ++i) {
                    for (int j = 0; j < sqrt_num_points; ++j) {
                        S2Point point = new(Math.Tan((i + 0.5) * spacing - radius),
                                      Math.Tan((j + 0.5) * spacing - radius), 1.0);
                        index.Add(S2.FromFrame(frame, point.Normalize()),
                                   i * sqrt_num_points + j);
                    }
                }
            }
        }

        // Use "query" to find the closest point(s) to the given target, and extract
        // the query results into the given vector.  Also verify that the results
        // satisfy the search criteria.
        private static void GetClosestPoints(Target target, TestQuery query, List<(S1ChordAngle, int)> results) {
            var query_results = query.FindClosestPoints(target);
            Assert.True(query_results.Count <= query.Options_.MaxResults);
            var region = query.Options_.Region;
            if (region == null && query.Options_.MaxDistance.Equals(S1ChordAngle.Infinity)) {
                // We can predict exactly how many points should be returned.
                Assert.Equal(Math.Min(query.Options_.MaxResults,
                                   query.Index().NumPoints()),
                          query_results.Count);
            }
            foreach (var result in query_results) {
                // Check that the point satisfies the region() condition.
                if (region != null) Assert.True(region.Contains(result.Point));

                // Check that it satisfies the max_distance() condition.
                Assert.True(result.Distance.IsLessThan(query.Options_.MaxDistance));
                results.Add((result.Distance.ToS1ChordAngle(), result.Data));
            }
        }

        private static void TestFindClosestPoints(Target target, TestQuery query) {
            List<(S1ChordAngle, int)> expected = new(), actual = new();
            query.Options_.UseBruteForce = (true);
            GetClosestPoints(target, query, expected);
            query.Options_.UseBruteForce = (false);
            GetClosestPoints(target, query, actual);
            Assert.True(S2TestingCheckDistance<int, S1ChordAngle>.CheckDistanceResults(expected, actual,
                                             query.Options_.MaxResults,
                                             query.Options_.MaxDistance.ToS1ChordAngle(),
                                             query.Options_.MaxError));

            if (!expected.Any()) return;

            // Note that when options.max_error() > 0, expected[0].distance may not be
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

        // (Note that every query is checked using the brute force algorithm.)
        private static void TestWithIndexFactory(IPointIndexFactory factory,
            int num_indexes, int num_points, int num_queries) {
            // Build a set of S2PointIndexes containing the desired geometry.
            List<S2Cap> index_caps = new();
            List<TestIndex> indexes = new();
            for (int i = 0; i < num_indexes; ++i) {
                S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
                index_caps.Add(new S2Cap(S2Testing.RandomPoint(), kTestCapRadius));
                indexes.Add(new TestIndex());
                factory.AddPoints(index_caps.Last(), num_points, indexes.Last());
            }
            for (int i = 0; i < num_queries; ++i) {
                S2Testing.Random.Reset(S2Testing.Random.RandomSeed + i);
                int i_index = S2Testing.Random.Uniform(num_indexes);
                S2Cap index_cap = index_caps[i_index];

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                S1Angle query_radius = 2 * index_cap.RadiusAngle();
                S2Cap query_cap = new(index_cap.Center, query_radius);
                TestQuery query = new(indexes[i_index]);

                // Occasionally we don't set any limit on the number of result points.
                // (This may return all points if we also don't set a distance limit.)
                if (!S2Testing.Random.OneIn(5)) {
                    query.Options_.MaxResults = (1 + S2Testing.Random.Uniform(10));
                }
                // We set a distance limit 2/3 of the time.
                if (!S2Testing.Random.OneIn(3)) {
                    query.Options_.MaxDistance = new S1ChordAngle(
                        S2Testing.Random.RandDouble() * query_radius);
                }
                if (S2Testing.Random.OneIn(2)) {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.Options_.MaxError = new S1ChordAngle(S1Angle.FromRadians(
                        Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius.Radians));
                }
                S2LatLngRect filter_rect = S2LatLngRect.FromCenterSize(
                    new S2LatLng(S2Testing.SamplePoint(query_cap)),
                    new S2LatLng(S2Testing.Random.RandDouble() * kTestCapRadius,
                             S2Testing.Random.RandDouble() * kTestCapRadius));
                if (S2Testing.Random.OneIn(5)) {
                    query.Options_.Region = (filter_rect);
                }
                int target_type = S2Testing.Random.Uniform(4);
                if (target_type == 0) {
                    // Find the points closest to a given point.
                    S2Point point = S2Testing.SamplePoint(query_cap);
                    TestQuery.PointTarget target = new(point);
                    TestFindClosestPoints(target, query);
                } else if (target_type == 1) {
                    // Find the points closest to a given edge.
                    S2Point a = S2Testing.SamplePoint(query_cap);
                    S2Point b = S2Testing.SamplePoint(
                        new S2Cap(a, Math.Pow(1e-4, S2Testing.Random.RandDouble()) * query_radius));
                    TestQuery.EdgeTarget target = new(a, b);
                    TestFindClosestPoints(target, query);
                } else if (target_type == 2) {
                    // Find the points closest to a given cell.
                    int min_level = S2.kMaxDiag.GetLevelForMaxValue(query_radius.Radians);
                    int level = min_level + S2Testing.Random.Uniform(
                        S2.kMaxCellLevel - min_level + 1);
                    S2Point a = S2Testing.SamplePoint(query_cap);
                    S2Cell cell = new(new S2CellId(a).Parent(level));
                    TestQuery.CellTarget target = new(cell);
                    TestFindClosestPoints(target, query);
                } else {
                    Assert.Equal(3, target_type);
                    MutableS2ShapeIndex target_index=new();
                    new FractalLoopShapeIndexFactory().AddEdges(index_cap, 100, target_index);
                    TestQuery.ShapeIndexTarget target = new(target_index);
                    target.IncludeInteriors = (S2Testing.Random.OneIn(2));
                    TestFindClosestPoints(target, query);
                }
            }
        }
    }
}
