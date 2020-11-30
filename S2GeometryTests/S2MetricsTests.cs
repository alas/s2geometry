using System;
using Xunit;
using Metric = S2Geometry.S2Metrics.Metric;

namespace S2Geometry
{
    public class S2MetricsTests
    {
        [Fact]
        public void Test_S2_Metrics()
        {
            var angle_span = new MetricBundle(S2Metrics.kMinAngleSpan, S2Metrics.kMaxAngleSpan, S2Metrics.kAvgAngleSpan);
            var widthBundle = new MetricBundle(S2Metrics.kMinWidth, S2Metrics.kMaxWidth, S2Metrics.kAvgWidth);
            var edge = new MetricBundle(S2Metrics.kMinEdge, S2Metrics.kMaxEdge, S2Metrics.kAvgEdge);
            var diag = new MetricBundle(S2Metrics.kMinDiag, S2Metrics.kMaxDiag, S2Metrics.kAvgDiag);
            var areaBundle = new MetricBundle(S2Metrics.kMinArea, S2Metrics.kMaxArea, S2Metrics.kAvgArea);

            // First, check that min <= avg <= max for each metric.
            CheckMinMaxAvg(angle_span);
            CheckMinMaxAvg(widthBundle);
            CheckMinMaxAvg(edge);
            CheckMinMaxAvg(diag);
            CheckMinMaxAvg(areaBundle);

            // Check that the maximum aspect ratio of an individual cell is consistent
            // with the global minimums and maximums.
            Assert.True(S2Metrics.kMaxEdgeAspect >= 1);
            Assert.True(S2Metrics.kMaxEdgeAspect <=
                S2Metrics.kMaxEdge.Deriv/ S2Metrics.kMinEdge.Deriv);
            Assert.True(S2Metrics.kMaxDiagAspect >= 1);
            Assert.True(S2Metrics.kMaxDiagAspect <=
                S2Metrics.kMaxDiag.Deriv/ S2Metrics.kMinDiag.Deriv);

            // Check various conditions that are provable mathematically.
            CheckLessOrEqual(widthBundle, angle_span);
            CheckLessOrEqual(widthBundle, edge);
            CheckLessOrEqual(edge, diag);

            Assert.True(S2Metrics.kMinArea.Deriv>=
                S2Metrics.kMinWidth.Deriv* S2Metrics.kMinEdge.Deriv- S2Constants.DoubleError);
            Assert.True(S2Metrics.kMaxArea.Deriv<=
                S2Metrics.kMaxWidth.Deriv* S2Metrics.kMaxEdge.Deriv+ S2Constants.DoubleError);

            // GetLevelForMaxValue() and friends have built-in assertions, we just need
            // to call these functions to test them.
            //
            // We don't actually check that the metrics are correct here, e.g. that
            // GetMinWidth(10) is a lower bound on the width of cells at level 10.
            // It is easier to check these properties in s2cell_test, since
            // S2Cell has methods to compute the cell vertices, etc.

            for (int level = -2; level <= S2Constants.kMaxCellLevel + 3; ++level)
            {
                var width = S2Metrics.kMinWidth.Deriv* Math.Pow(2, -level);
                if (level >= S2Constants.kMaxCellLevel + 3) width = 0;

                // Check boundary cases (exactly equal to a threshold value).
                int expected_level = Math.Max(0, Math.Min(S2Constants.kMaxCellLevel, level));
                Assert.Equal(S2Metrics.kMinWidth.GetLevelForMaxValue(width), expected_level);
                Assert.Equal(S2Metrics.kMinWidth.GetLevelForMinValue(width), expected_level);
                Assert.Equal(S2Metrics.kMinWidth.GetClosestLevel(width), expected_level);

                // Also check non-boundary cases.
                Assert.Equal(S2Metrics.kMinWidth.GetLevelForMaxValue(1.2 * width), expected_level);
                Assert.Equal(S2Metrics.kMinWidth.GetLevelForMinValue(0.8 * width), expected_level);
                Assert.Equal(S2Metrics.kMinWidth.GetClosestLevel(1.2 * width), expected_level);
                Assert.Equal(S2Metrics.kMinWidth.GetClosestLevel(0.8 * width), expected_level);

                // Same thing for area.
                double area = S2Metrics.kMinArea.Deriv* Math.Pow(4, -level);
                if (level <= -3) area = 0;
                Assert.Equal(S2Metrics.kMinArea.GetLevelForMaxValue(area), expected_level);
                Assert.Equal(S2Metrics.kMinArea.GetLevelForMinValue(area), expected_level);
                Assert.Equal(S2Metrics.kMinArea.GetClosestLevel(area), expected_level);
                Assert.Equal(S2Metrics.kMinArea.GetLevelForMaxValue(1.2 * area), expected_level);
                Assert.Equal(S2Metrics.kMinArea.GetLevelForMinValue(0.8 * area), expected_level);
                Assert.Equal(S2Metrics.kMinArea.GetClosestLevel(1.2 * area), expected_level);
                Assert.Equal(S2Metrics.kMinArea.GetClosestLevel(0.8 * area), expected_level);
            }
        }

        private static void CheckMinMaxAvg(MetricBundle bundle)
        {
            Assert.True(bundle.min_.Deriv<= bundle.avg_.Deriv);
            Assert.True(bundle.avg_.Deriv<= bundle.max_.Deriv);
        }

        private static void CheckLessOrEqual(MetricBundle a, MetricBundle b)
        {
            Assert.True(a.min_.Deriv<= b.min_.Deriv);
            Assert.True(a.max_.Deriv<= b.max_.Deriv);
            Assert.True(a.avg_.Deriv<= b.avg_.Deriv);
        }

        // Note: obviously, I could have defined a bundle of metrics like this in the
        // S2 class itself rather than just for testing.  However, it's not clear that
        // this is useful other than for testing purposes, and I find
        // S2::kMinWidth.GetLevelForMinValue(width) to be slightly more readable than
        // than S2::kWidth.min().GetLevelForMinValue(width).  Also, there is no
        // fundamental reason that we need to analyze the minimum, maximum, and average
        // values of every metric; it would be perfectly reasonable to just define
        // one of these.
        private class MetricBundle
        {
            public MetricBundle(Metric min, Metric max, Metric avg)
            { min_ = min; max_ = max; avg_ = avg; }

            public Metric min_;
            public Metric max_;
            public Metric avg_;
        }
    }
}
