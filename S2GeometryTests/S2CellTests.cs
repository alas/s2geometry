using System;
using System.Collections.Generic;
using Xunit;
using Xunit.Abstractions;

namespace S2Geometry
{
    public class S2CellTests
    {
        private readonly ITestOutputHelper _logger;

        public S2CellTests(ITestOutputHelper logger) { _logger = logger; }

        [Fact]
        public void Test_S2Cell_TestFaces()
        {
            var edge_counts = new Dictionary<S2Point, int>();
            var vertex_counts = new Dictionary<S2Point, int>();
            for (int face = 0; face < 6; ++face)
            {
                S2CellId id = S2CellId.FromFace(face);
                S2Cell cell = new S2Cell(id);
                Assert.Equal(id, cell.Id);
                Assert.Equal(face, cell.Face);
                Assert.Equal(0, cell.Level);
                // Top-level faces have alternating orientations to get RHS coordinates.
                Assert.Equal(face & S2Constants.kSwapMask, cell.Orientation);
                Assert.False(cell.IsLeaf);
                for (int k = 0; k < 4; ++k)
                {
                    var key = cell.GetEdgeRaw(k);
                    if (edge_counts.ContainsKey(key))
                        edge_counts[key] += 1;
                    else
                        edge_counts[key] = 1;
                    var key2 = cell.GetVertexRaw(k);
                    if (vertex_counts.ContainsKey(key2))
                        vertex_counts[key2] += 1;
                    else
                        vertex_counts[key2] = 1;
                    Assert2.Near(0.0, cell.GetVertexRaw(k).DotProd(key));
                    Assert2.Near(0.0, cell.GetVertexRaw(k + 1).DotProd(key));
                    Assert2.Near(1.0, cell.GetVertexRaw(k).CrossProd(cell.GetVertexRaw(k + 1))
                        .Normalized.DotProd(cell.GetEdge(k)));
                }
            }
            // Check that edges have multiplicity 2 and vertices have multiplicity 3.
            foreach (var p in edge_counts)
            {
                Assert.Equal(2, p.Value);
            }
            foreach (var p in vertex_counts)
            {
                Assert.Equal(3, p.Value);
            }
        }

        private class LevelStats
        {
            public double Count = 0;
            public double Min_area = 0, Max_area = 0, Avg_area = 0;
            public double Min_width = 0, Max_width = 0, Avg_width = 0;
            public double Min_edge = 0, Max_edge = 0, Avg_edge = 0, Max_edge_aspect = 0;
            public double Min_diag = 0, Max_diag = 0, Avg_diag = 0, Max_diag_aspect = 0;
            public double Min_angle_span = 0, Max_angle_span = 0, Avg_angle_span = 0;
            public double Min_approx_ratio = 0, Max_approx_ratio = 0;
        
            public static LevelStats GetDefault()
            {
                return new LevelStats {
                    Min_area = 100,
                    Min_width = 100,
                    Min_edge = 100,
                    Min_diag = 100,
                    Min_angle_span = 100,
                    Min_approx_ratio = 100,
                };
            }
        }
        private static readonly LevelStats[] level_stats = new LevelStats[S2Constants.kMaxCellLevel + 1];

        private static void GatherStats(S2Cell cell)
        {
            var s = level_stats[cell.Level];
            double exact_area = cell.ExactArea();
            double approx_area = cell.ApproxArea();
            double min_edge = 100, max_edge = 0, avg_edge = 0;
            double min_diag = 100, max_diag = 0;
            double min_width = 100, max_width = 0;
            double min_angle_span = 100, max_angle_span = 0;
            for (int i = 0; i < 4; ++i)
            {
                double edge = cell.GetVertexRaw(i).Angle(cell.GetVertexRaw(i + 1));
                min_edge = Math.Min(edge, min_edge);
                max_edge = Math.Max(edge, max_edge);
                avg_edge += 0.25 * edge;
                S2Point mid = cell.GetVertexRaw(i) + cell.GetVertexRaw(i + 1);
                double width = S2Constants.M_PI_2 - mid.Angle(cell.GetEdgeRaw(i + 2));
                min_width = Math.Min(width, min_width);
                max_width = Math.Max(width, max_width);
                if (i < 2)
                {
                    double diag = cell.GetVertexRaw(i).Angle(cell.GetVertexRaw(i + 2));
                    min_diag = Math.Min(diag, min_diag);
                    max_diag = Math.Max(diag, max_diag);
                    double angle_span = cell.GetEdgeRaw(i).Angle(-cell.GetEdgeRaw(i + 2));
                    min_angle_span = Math.Min(angle_span, min_angle_span);
                    max_angle_span = Math.Max(angle_span, max_angle_span);
                }
            }
            s.Count += 1;
            s.Min_area = Math.Min(exact_area, s.Min_area);
            s.Max_area = Math.Max(exact_area, s.Max_area);
            s.Avg_area += exact_area;
            s.Min_width = Math.Min(min_width, s.Min_width);
            s.Max_width = Math.Max(max_width, s.Max_width);
            s.Avg_width += 0.5 * (min_width + max_width);
            s.Min_edge = Math.Min(min_edge, s.Min_edge);
            s.Max_edge = Math.Max(max_edge, s.Max_edge);
            s.Avg_edge += avg_edge;
            s.Max_edge_aspect = Math.Max(max_edge / min_edge, s.Max_edge_aspect);
            s.Min_diag = Math.Min(min_diag, s.Min_diag);
            s.Max_diag = Math.Max(max_diag, s.Max_diag);
            s.Avg_diag += 0.5 * (min_diag + max_diag);
            s.Max_diag_aspect = Math.Max(max_diag / min_diag, s.Max_diag_aspect);
            s.Min_angle_span = Math.Min(min_angle_span, s.Min_angle_span);
            s.Max_angle_span = Math.Max(max_angle_span, s.Max_angle_span);
            s.Avg_angle_span += 0.5 * (min_angle_span + max_angle_span);
            double approx_ratio = approx_area / exact_area;
            s.Min_approx_ratio = Math.Min(approx_ratio, s.Min_approx_ratio);
            s.Max_approx_ratio = Math.Max(approx_ratio, s.Max_approx_ratio);
        }

        static void TestSubdivide(S2Cell cell)
        {
            GatherStats(cell);
            if (cell.IsLeaf) return;

            var children = new S2Cell[4];
            Assert.True(cell.Subdivide(children));
            S2CellId child_id = cell.Id.ChildBegin();
            double exact_area = 0;
            double approx_area = 0;
            double average_area = 0;
            for (int i = 0; i < 4; ++i, child_id = child_id.Next)
            {
                exact_area += children[i].ExactArea();
                approx_area += children[i].ApproxArea();
                average_area += children[i].AverageArea();

                // Check that the child geometry is consistent with its cell ID.
                Assert.Equal(child_id, children[i].Id);
                Assert.True(S2PointUtil.ApproxEquals(children[i].GetCenter(), child_id.ToPoint()));
                S2Cell direct = new S2Cell(child_id);
                Assert.Equal(direct.Face, children[i].Face);
                Assert.Equal(direct.Level, children[i].Level);
                Assert.Equal(direct.Orientation, children[i].Orientation);
                Assert.Equal(direct.GetCenterRaw(), children[i].GetCenterRaw());
                for (int k = 0; k < 4; ++k)
                {
                    Assert.Equal(direct.GetVertexRaw(k), children[i].GetVertexRaw(k));
                    Assert.Equal(direct.GetEdgeRaw(k), children[i].GetEdgeRaw(k));
                }

                // Test Contains() and MayIntersect().
                Assert.True(cell.Contains(children[i]));
                Assert.True(cell.MayIntersect(children[i]));
                Assert.False(children[i].Contains(cell));
                Assert.True(cell.Contains(children[i].GetCenterRaw()));
                for (int j = 0; j < 4; ++j)
                {
                    Assert.True(cell.Contains(children[i].GetVertexRaw(j)));
                    if (j != i)
                    {
                        Assert.False(children[i].Contains(children[j].GetCenterRaw()));
                        Assert.False(children[i].MayIntersect(children[j]));
                    }
                }

                // Test GetCapBound and GetRectBound.
                S2Cap parent_cap = cell.GetCapBound();
                S2LatLngRect parent_rect = cell.GetRectBound();
                if (cell.Contains(new S2Point(0, 0, 1)) || cell.Contains(new S2Point(0, 0, -1)))
                {
                    Assert.True(parent_rect.Lng.IsFull);
                }
                S2Cap child_cap = children[i].GetCapBound();
                S2LatLngRect child_rect = children[i].GetRectBound();
                Assert.True(child_cap.Contains(children[i].GetCenter()));
                Assert.True(child_rect.Contains(children[i].GetCenterRaw()));
                Assert.True(parent_cap.Contains(children[i].GetCenter()));
                Assert.True(parent_rect.Contains(children[i].GetCenterRaw()));
                for (int j = 0; j < 4; ++j)
                {
                    Assert.True(child_cap.Contains(children[i].GetVertex(j)));
                    Assert.True(child_rect.Contains(children[i].GetVertex(j)));
                    Assert.True(child_rect.Contains(children[i].GetVertexRaw(j)));
                    Assert.True(parent_cap.Contains(children[i].GetVertex(j)));
                    Assert.True(parent_rect.Contains(children[i].GetVertex(j)));
                    Assert.True(parent_rect.Contains(children[i].GetVertexRaw(j)));
                    if (j != i)
                    {
                        // The bounding caps and rectangles should be tight enough so that
                        // they exclude at least two vertices of each adjacent cell.
                        int cap_count = 0;
                        int rect_count = 0;
                        for (int k = 0; k < 4; ++k)
                        {
                            if (child_cap.Contains(children[j].GetVertex(k)))
                                ++cap_count;
                            if (child_rect.Contains(children[j].GetVertexRaw(k)))
                                ++rect_count;
                        }
                        Assert.True(cap_count <= 2);
                        if (child_rect.LatLo.Radians > -S2Constants.M_PI_2 &&
                            child_rect.LatHi.Radians < S2Constants.M_PI_2)
                        {
                            // Bounding rectangles may be too large at the poles because the
                            // pole itself has an arbitrary fixed longitude.
                            Assert.True(rect_count <= 2);
                        }
                    }
                }

                // Check all children for the first few levels, and then sample randomly.
                // We also always subdivide the cells containing a few chosen points so
                // that we have a better chance of sampling the minimum and maximum metric
                // values.  kMaxSizeUV is the absolute value of the u- and v-coordinate
                // where the cell size at a given level is maximal.
                double kMaxSizeUV = 0.3964182625366691;
                R2Point[] special_uv = {
     new R2Point(S2Constants.DoubleEpsilon, S2Constants.DoubleEpsilon),  // Face center
     new R2Point(S2Constants.DoubleEpsilon, 1),            // Edge midpoint
     new R2Point(1, 1),                      // Face corner
     new R2Point(kMaxSizeUV, kMaxSizeUV),    // Largest cell area
     new R2Point(S2Constants.DoubleEpsilon, kMaxSizeUV),   // Longest edge/diagonal
    };
                bool force_subdivide = false;
                foreach (R2Point uv in special_uv)
                {
                    if (children[i].BoundUV.Contains(uv))
                        force_subdivide = true;
                }

                var debugFlag =
#if s2debug
                true;
#else
                false;
#endif

                if (force_subdivide ||
        cell.Level < (debugFlag ? 5 : 6) ||
        S2Testing.Random.OneIn(debugFlag ? 5 : 4))
                {
                    TestSubdivide(children[i]);
                }
            }

            // Check sum of child areas equals parent area.
            //
            // For ExactArea(), the best relative error we can expect is about 1e-6
            // because the precision of the unit vector coordinates is only about 1e-15
            // and the edge length of a leaf cell is about 1e-9.
            //
            // For ApproxArea(), the areas are accurate to within a few percent.
            //
            // For AverageArea(), the areas themselves are not very accurate, but
            // the average area of a parent is exactly 4 times the area of a child.

            Assert.True(Math.Abs(Math.Log(exact_area / cell.ExactArea())) <= Math.Abs(Math.Log((1 + 1e-6))));
            Assert.True(Math.Abs(Math.Log((approx_area / cell.ApproxArea()))) <= Math.Abs(Math.Log((1.03))));
            Assert.True(Math.Abs(Math.Log((average_area / cell.AverageArea()))) <= Math.Abs(Math.Log((1 + 1e-15))));
        }

        private void CheckMinMaxAvg(
            string label, int level, double count, double abs_error,
            double min_value, double max_value, double avg_value,
            S2Metrics.Metric min_metric,
            S2Metrics.Metric max_metric,
            S2Metrics.Metric avg_metric)
        {

            // All metrics are minimums, maximums, or averages of differential
            // quantities, and therefore will not be exact for cells at any finite
            // level.  The differential minimum is always a lower bound, and the maximum
            // is always an upper bound, but these minimums and maximums may not be
            // achieved for two different reasons.  First, the cells at each level are
            // sampled and we may miss the most extreme examples.  Second, the actual
            // metric for a cell is obtained by integrating the differential quantity,
            // which is notant across the cell.  Therefore cells at low levels
            // (bigger cells) have smaller variations.
            //
            // The "tolerance" below is an attempt to model both of these effects.
            // At low levels, error is dominated by the variation of differential
            // quantities across the cells, while at high levels error is dominated by
            // the effects of random sampling.
            double tolerance = (max_metric.GetValue(level) - min_metric.GetValue(level)) /
                               Math.Sqrt(Math.Min(count, 0.5 * (1 << level)));
            if (tolerance == 0) tolerance = abs_error;

            double min_error = min_value - min_metric.GetValue(level);
            double max_error = max_metric.GetValue(level) - max_value;
            double avg_error = Math.Abs(avg_metric.GetValue(level) - avg_value);
            _logger.WriteLine("%-10s (%6.0f samples, tolerance %8.3g) - min %9.4g (%9.3g : %9.3g) " +
         "max %9.4g (%9.3g : %9.3g), avg %9.4g (%9.3g : %9.3g)",
         label, count, tolerance,
         min_value, min_error / min_value, min_error / tolerance,
         max_value, max_error / max_value, max_error / tolerance,
         avg_value, avg_error / avg_value, avg_error / tolerance);

            Assert.True(min_metric.GetValue(level) <= min_value + abs_error);
            Assert.True(min_metric.GetValue(level) >= min_value - tolerance);
            Assert.True(max_metric.GetValue(level) <= max_value + tolerance);
            Assert.True(max_metric.GetValue(level) >= max_value - abs_error);
            Assert2.Near(avg_metric.GetValue(level), avg_value, 10 * tolerance);
        }

        [Fact]
        public void Test_S2Cell_TestSubdivide()
        {
            // Only test a sample of faces to reduce the runtime.
            TestSubdivide(S2Cell.FromFace(0));
            TestSubdivide(S2Cell.FromFace(3));
            TestSubdivide(S2Cell.FromFace(5));

            #region Print Header
            
            // This table is useful in evaluating the quality of the various S2
            // projections.
            //
            // The maximum edge *ratio* is the ratio of the longest edge of any cell to
            // the shortest edge of any cell at the same level (and similarly for the
            // maximum diagonal ratio).
            //
            // The maximum edge *aspect* is the maximum ratio of the longest edge of a
            // cell to the shortest edge of that same cell (and similarly for the
            // maximum diagonal aspect).
            _logger.WriteLine(
@"Ratio:  (Max value for any cell) / (Min value for any cell)
Aspect: (Max value / min value) for any cell
                   Edge          Diag       Approx Area/    Avg Area/
         Area     Length        Length       Exact Area    Exact Area
Level   Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max
--------------------------------------------------------------------"); 

            #endregion

            for (int i = 0; i <= S2Constants.kMaxCellLevel; ++i)
            {
                var s = level_stats[i];
                if (s.Count > 0)
                {
                    s.Avg_area /= s.Count;
                    s.Avg_width /= s.Count;
                    s.Avg_edge /= s.Count;
                    s.Avg_diag /= s.Count;
                    s.Avg_angle_span /= s.Count;
                }
                _logger.WriteLine($"{i:d5}  {s.Max_area / s.Min_area:f6.3} {s.Max_edge / s.Min_edge:f6.3} {s.Max_edge_aspect:f6.3} {s.Max_diag / s.Min_diag:f6.3} {s.Max_diag_aspect:f6.3} {s.Min_approx_ratio:f6.3} {s.Max_approx_ratio:f6.3} {S2Cell.AverageArea(i) / s.Max_area:f6.3} {S2Cell.AverageArea(i) / s.Min_area:f6.3}");
            }

            // Now check the validity of the S2 length and area metrics.
            for (int i = 0; i <= S2Constants.kMaxCellLevel; ++i)
            {
                var s = level_stats[i];
                if (s.Count == 0) continue;

                _logger.WriteLine($"Level {i:2d} - metric value (error/actual : error/tolerance)");

                // The various length calculations are only accurate to 1e-15 or so,
                // so we need to allow for this amount of discrepancy with the theoretical
                // minimums and maximums.  The area calculation is accurate to about 1e-15
                // times the cell width.
                CheckMinMaxAvg("area", i, s.Count, 1e-15 * s.Min_width,
                               s.Min_area, s.Max_area, s.Avg_area,
                               S2Metrics.kMinArea, S2Metrics.kMaxArea, S2Metrics.kAvgArea);
                CheckMinMaxAvg("width", i, s.Count, 1e-15,
                               s.Min_width, s.Max_width, s.Avg_width,
                               S2Metrics.kMinWidth, S2Metrics.kMaxWidth, S2Metrics.kAvgWidth);
                CheckMinMaxAvg("edge", i, s.Count, 1e-15,
                               s.Min_edge, s.Max_edge, s.Avg_edge,
                               S2Metrics.kMinEdge, S2Metrics.kMaxEdge, S2Metrics.kAvgEdge);
                CheckMinMaxAvg("diagonal", i, s.Count, 1e-15,
                               s.Min_diag, s.Max_diag, s.Avg_diag,
                               S2Metrics.kMinDiag, S2Metrics.kMaxDiag, S2Metrics.kAvgDiag);
                CheckMinMaxAvg("angle span", i, s.Count, 1e-15,
                               s.Min_angle_span, s.Max_angle_span, s.Avg_angle_span,
                               S2Metrics.kMinAngleSpan, S2Metrics.kMaxAngleSpan, S2Metrics.kAvgAngleSpan);

                // The aspect ratio calculations are ratios of lengths and are therefore
                // less accurate at higher subdivision levels.
                Assert.True(s.Max_edge_aspect <= S2Metrics.kMaxEdgeAspect + 1e-15 * (1 << i));
                Assert.True(s.Max_diag_aspect <= S2Metrics.kMaxDiagAspect + 1e-15 * (1 << i));
            }
        }

        [Fact]
        public void Test_S2Cell_CellVsLoopRectBound()
        {
            // This test verifies that the S2Cell and S2Loop bounds contain each other
            // to within their maximum errors.
            //
            // The S2Cell and S2Loop calculations for the latitude of a vertex can differ
            // by up to 2 * S2Constants.DoubleEpsilon, therefore the S2Cell bound should never exceed
            // the S2Loop bound by more than this (the reverse is not true, because the
            // S2Loop code sometimes thinks that the maximum occurs along an edge).
            // Similarly, the longitude bounds can differ by up to 4 * S2Constants.DoubleEpsilon since
            // the S2Cell bound has an error of 2 * S2Constants.DoubleEpsilon and then expands by this
            // amount, while the S2Loop bound does no expansion at all.

            for (int iter = 0; iter < 1000; ++iter)
            {
                S2Cell cell = new S2Cell(S2Testing.GetRandomCellId());
                S2Loop loop = new S2Loop(cell);
                S2LatLngRect cell_bound = cell.GetRectBound();
                S2LatLngRect loop_bound = loop.GetRectBound();
                Assert.True(loop_bound.Expanded(kCellError).Contains(cell_bound));
                Assert.True(cell_bound.Expanded(kLoopError).Contains(loop_bound));
            }
        }
        // Possible additional S2Cell error compared to S2Loop error:
        private static readonly S2LatLng kCellError = S2LatLng.FromRadians(2 * S2Constants.DoubleEpsilon, 4 * S2Constants.DoubleEpsilon);
        // Possible additional S2Loop error compared to S2Cell error:
        private static readonly S2LatLng kLoopError = S2LatLngRectBounder.MaxErrorForTests();

        [Fact]
        public void Test_S2Cell_RectBoundIsLargeEnough()
        {
            // Construct many points that are nearly on an S2Cell edge, and verify that
            // whenever the cell contains a point P then its bound contains S2LatLng(P).

            for (int iter = 0; iter < 1000; /* advanced in loop below */)
            {
                S2Cell cell = new S2Cell(S2Testing.GetRandomCellId());
                int i = S2Testing.Random.Uniform(4);
                S2Point v1 = cell.GetVertex(i);
                S2Point v2 = S2Testing.SamplePoint(
                    new S2Cap(cell.GetVertex(i + 1), S1Angle.FromRadians(1e-15)));
                S2Point p = S2EdgeDistances.Interpolate(S2Testing.Random.RandDouble(), v1, v2);
                if (new S2Loop(cell).Contains(p))
                {
                    Assert.True(cell.GetRectBound().Contains(new S2LatLng(p)));
                    ++iter;
                }
            }
        }

        [Fact]
        public void Test_S2Cell_ConsistentWithS2CellIdFromPoint()
        {
            // Construct many points that are nearly on an S2Cell edge, and verify that
            // S2Cell(S2CellId(p)).Contains(p) is always true.

            for (int iter = 0; iter < 1000; ++iter)
            {
                S2Cell cell = new S2Cell(S2Testing.GetRandomCellId());
                int i = S2Testing.Random.Uniform(4);
                S2Point v1 = cell.GetVertex(i);
                S2Point v2 = S2Testing.SamplePoint(
                    new S2Cap(cell.GetVertex(i + 1), S1Angle.FromRadians(1e-15)));
                S2Point p = S2EdgeDistances.Interpolate(S2Testing.Random.RandDouble(), v1, v2);
                Assert.True(new S2Cell(new S2CellId(p)).Contains(p));
            }
        }

        [Fact]
        public void Test_S2CellId_AmbiguousContainsPoint()
        {
            // This tests a case where S2CellId returns the "wrong" cell for a point
            // that is very close to the cell edge. (ConsistentWithS2CellIdFromPoint
            // generates more examples like this.)
            //
            // The S2Point below should have x = 0, but conversion from latlng to
            // (x,y,z) gives x = 6.1e-17.  When xyz is converted to uv, this gives u =
            // -6.1e-17.  However when converting to st, which is centered at 0.5 rather
            // than 0, the low precision bits of u are lost and we wind up with s = 0.5.
            // S2CellId(S2Point) then chooses an arbitrary neighboring cell.
            //
            // This tests that S2Cell.Contains() expands the cell bounds sufficiently
            // so that the returned cell is still considered to contain "p".
            S2Point p = S2LatLng.FromDegrees(-2, 90).ToPoint();
            S2CellId cell_id = new S2CellId(p).Parent(1);
            S2Cell cell = new S2Cell(cell_id);
            Assert.True(cell.Contains(p));
        }

        private static S1ChordAngle GetDistanceToPointBruteForce(S2Cell cell, S2Point target)
        {
            S1ChordAngle min_distance = S1ChordAngle.Infinity;
            for (int i = 0; i < 4; ++i)
            {
                S2EdgeDistances.UpdateMinDistance(target, cell.GetVertex(i),
                                              cell.GetVertex(i + 1), ref min_distance);
            }
            return min_distance;
        }

        private static S1ChordAngle GetMaxDistanceToPointBruteForce(S2Cell cell, S2Point target)
        {
            if (cell.Contains(-target))
            {
                return S1ChordAngle.Straight;
            }
            S1ChordAngle max_distance = S1ChordAngle.Negative;
            for (int i = 0; i < 4; ++i)
            {
                S2EdgeDistances.UpdateMaxDistance(target, cell.GetVertex(i),
                                          cell.GetVertex(i + 1), ref max_distance);
            }
            return max_distance;
        }

        [Fact]
        public void Test_S2Cell_GetDistanceToPoint()
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int iter = 0; iter < 1000; ++iter)
            {
                _logger.WriteLine($"Iteration {iter}");
                S2Cell cell = new S2Cell(S2Testing.GetRandomCellId());
                S2Point target = S2Testing.RandomPoint();
                S1Angle expected_to_boundary =
                    GetDistanceToPointBruteForce(cell, target).ToAngle();
                S1Angle expected_to_interior =
                    cell.Contains(target) ? S1Angle.Zero : expected_to_boundary;
                S1Angle expected_max =
                    GetMaxDistanceToPointBruteForce(cell, target).ToAngle();
                S1Angle actual_to_boundary = cell.GetBoundaryDistance(target).ToAngle();
                S1Angle actual_to_interior = cell.GetDistance(target).ToAngle();
                S1Angle actual_max = cell.GetMaxDistance(target).ToAngle();
                // The error has a peak near Pi/2 for edge distance, and another peak near
                // Pi for vertex distance.
                Assert2.Near(expected_to_boundary.Radians,
                            actual_to_boundary.Radians, 1e-12);
                Assert2.Near(expected_to_interior.Radians,
                            actual_to_interior.Radians, 1e-12);
                Assert2.Near(expected_max.Radians,
                            actual_max.Radians, 1e-12);
                if (expected_to_boundary.Radians <= Math.PI / 3)
                {
                    Assert2.Near(expected_to_boundary.Radians,
                                actual_to_boundary.Radians, S2Constants.DoubleError);
                    Assert2.Near(expected_to_interior.Radians,
                                actual_to_interior.Radians, S2Constants.DoubleError);
                }
                if (expected_max.Radians <= Math.PI / 3)
                {
                    Assert2.Near(expected_max.Radians, actual_max.Radians, S2Constants.DoubleError);
                }
            }
        }

        private static void ChooseEdgeNearCell(S2Cell cell, out S2Point a, out S2Point b)
        {
            S2Cap cap = cell.GetCapBound();
            if (S2Testing.Random.OneIn(5))
            {
                // Choose a point anywhere on the sphere.
                a = S2Testing.RandomPoint();
            }
            else
            {
                // Choose a point inside or somewhere near the cell.
                a = S2Testing.SamplePoint(new S2Cap(cap.Center, 1.5 * cap.RadiusAngle));
            }
            // Now choose a maximum edge length ranging from very short to very long
            // relative to the cell size, and choose the other endpoint.
            double max_length = Math.Min(100 * Math.Pow(1e-4, S2Testing.Random.RandDouble()) *
                                    cap.Radius.Radians, S2Constants.M_PI_2);
            b = S2Testing.SamplePoint(new S2Cap(a, S1Angle.FromRadians(max_length)));

            if (S2Testing.Random.OneIn(20))
            {
                // Occasionally replace edge with antipodal edge.
                a = -a;
                b = -b;
            }
        }

        private static S1ChordAngle GetDistanceToEdgeBruteForce(S2Cell cell, S2Point a, S2Point b)
        {
            if (cell.Contains(a) || cell.Contains(b))
            {
                return S1ChordAngle.Zero;
            }

            S1ChordAngle min_dist = S1ChordAngle.Infinity;
            for (int i = 0; i < 4; ++i)
            {
                S2Point v0 = cell.GetVertex(i);
                S2Point v1 = cell.GetVertex(i + 1);
                // If the edge crosses through the cell, max distance is 0.
                if (S2EdgeCrossings.CrossingSign(a, b, v0, v1) >= 0)
                {
                    return S1ChordAngle.Zero;
                }
                S2EdgeDistances.UpdateMinDistance(a, v0, v1, ref min_dist);
                S2EdgeDistances.UpdateMinDistance(b, v0, v1, ref min_dist);
                S2EdgeDistances.UpdateMinDistance(v0, a, b, ref min_dist);
            }
            return min_dist;
        }

        private static S1ChordAngle GetMaxDistanceToEdgeBruteForce(S2Cell cell, S2Point a, S2Point b)
        {
            // If any antipodal endpoint is within the cell, the max distance is Pi.
            if (cell.Contains(-a) || cell.Contains(-b))
            {
                return S1ChordAngle.Straight;
            }

            S1ChordAngle max_dist = S1ChordAngle.Negative;
            for (int i = 0; i < 4; ++i)
            {
                S2Point v0 = cell.GetVertex(i);
                S2Point v1 = cell.GetVertex(i + 1);
                // If the antipodal edge crosses through the cell, max distance is Pi.
                if (S2EdgeCrossings.CrossingSign(-a, -b, v0, v1) >= 0)
                {
                    return S1ChordAngle.Straight;
                }
                S2EdgeDistances.UpdateMaxDistance(a, v0, v1, ref max_dist);
                S2EdgeDistances.UpdateMaxDistance(b, v0, v1, ref max_dist);
                S2EdgeDistances.UpdateMaxDistance(v0, a, b, ref max_dist);
            }
            return max_dist;
        }

        [Fact]
        public void Test_S2Cell_GetDistanceToEdge()
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int iter = 0; iter < 1000; ++iter)
            {
                _logger.WriteLine($"Iteration {iter}");
                S2Cell cell = new S2Cell(S2Testing.GetRandomCellId());
                ChooseEdgeNearCell(cell, out var a, out var b);
                S1Angle expected_min = GetDistanceToEdgeBruteForce(cell, a, b).ToAngle();
                S1Angle expected_max =
                    GetMaxDistanceToEdgeBruteForce(cell, a, b).ToAngle();
                S1Angle actual_min = cell.GetDistance(a, b).ToAngle();
                S1Angle actual_max = cell.GetMaxDistance(a, b).ToAngle();
                // The error has a peak near Pi/2 for edge distance, and another peak near
                // Pi for vertex distance.
                if (expected_min.Radians > S2Constants.M_PI_2)
                {
                    // Max error for S1ChordAngle as it approaches Pi is about 2e-8.
                    Assert2.Near(expected_min.Radians, actual_min.Radians, 2e-8);
                }
                else if (expected_min.Radians <= Math.PI / 3)
                {
                    Assert2.Near(expected_min.Radians, actual_min.Radians, S2Constants.DoubleError);
                }
                else
                {
                    Assert2.Near(expected_min.Radians, actual_min.Radians, 1e-12);
                }

                Assert2.Near(expected_max.Radians, actual_max.Radians, 1e-12);
                if (expected_max.Radians <= Math.PI / 3)
                {
                    Assert2.Near(expected_max.Radians, actual_max.Radians, S2Constants.DoubleError);
                }
            }
        }

        [Fact]
        public void Test_S2Cell_GetMaxDistanceToEdge()
        {
            // Test an edge for which its antipode crosses the cell. Validates both the
            // standard and brute force implementations for this case.
            S2Cell cell = S2Cell.FromFacePosLevel(0, 0, 20);
            S2Point a = -S2EdgeDistances.Interpolate(2.0, cell.GetCenter(), cell.GetVertex(0));
            S2Point b = -S2EdgeDistances.Interpolate(2.0, cell.GetCenter(), cell.GetVertex(2));

            S1ChordAngle actual = cell.GetMaxDistance(a, b);
            S1ChordAngle expected = GetMaxDistanceToEdgeBruteForce(cell, a, b);

            Assert2.Near(expected.Radians, S1ChordAngle.Straight.Radians, S2Constants.DoubleError);
            Assert2.Near(actual.Radians, S1ChordAngle.Straight.Radians, S2Constants.DoubleError);
        }

        [Fact]
        public void Test_S2Cell_GetMaxDistanceToCellAntipodal()
        {
            S2Point p = S2TextFormat.MakePointOrDie("0:0");
            S2Cell cell = new S2Cell(p);
            S2Cell antipodal_cell = new S2Cell(-p);
            S1ChordAngle dist = cell.GetMaxDistance(antipodal_cell);
            Assert.Equal(S1ChordAngle.Straight, dist);
        }

        [Fact]
        public void Test_S2Cell_GetMaxDistanceToCell()
        {
            for (int i = 0; i < 1000; i++)
            {
                S2Cell cell = new S2Cell(S2Testing.GetRandomCellId());
                S2Cell test_cell = new S2Cell(S2Testing.GetRandomCellId());
                S2CellId antipodal_leaf_id = new S2CellId(-test_cell.GetCenter());
                S2Cell antipodal_test_cell = new S2Cell(antipodal_leaf_id.Parent(test_cell.Level));

                S1ChordAngle dist_from_min = S1ChordAngle.Straight -
                    cell.GetDistance(antipodal_test_cell);
                S1ChordAngle dist_from_max = cell.GetMaxDistance(test_cell);
                Assert2.Near(dist_from_min.Radians, dist_from_max.Radians, 1e-8);
            }
        }

        [Fact]
        public void Test_S2Cell_EncodeDecode()
        {
            S2Cell orig_cell = new S2Cell(S2LatLng.FromDegrees(40.7406264, -74.0029963));
            Encoder encoder = new Encoder();
            orig_cell.Encode(encoder);
            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            var (_, decoded_cell) = S2Cell.DecodeStatic(decoder);

            Assert.Equal(orig_cell, decoded_cell);
            Assert.Equal(orig_cell.Face, decoded_cell.Face);
            Assert.Equal(orig_cell.Level, decoded_cell.Level);
            Assert.Equal(orig_cell.Orientation, decoded_cell.Orientation);
            Assert.Equal(orig_cell.Id, decoded_cell.Id);
            Assert.Equal(orig_cell.BoundUV, decoded_cell.BoundUV);
        }
    }
}
