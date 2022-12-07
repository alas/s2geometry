namespace S2Geometry;

using Metric = S2.Metric;

public class S2MetricsTests
{
    [Fact]
    internal void Test_S2_Metrics()
    {
        var angle_span = new MetricBundle(S2.kMinAngleSpan, S2.kMaxAngleSpan, S2.kAvgAngleSpan);
        var widthBundle = new MetricBundle(S2.kMinWidth, S2.kMaxWidth, S2.kAvgWidth);
        var edge = new MetricBundle(S2.kMinEdge, S2.kMaxEdge, S2.kAvgEdge);
        var diag = new MetricBundle(S2.kMinDiag, S2.kMaxDiag, S2.kAvgDiag);
        var areaBundle = new MetricBundle(S2.kMinArea, S2.kMaxArea, S2.kAvgArea);

        // First, check that min <= avg <= max for each metric.
        CheckMinMaxAvg(angle_span);
        CheckMinMaxAvg(widthBundle);
        CheckMinMaxAvg(edge);
        CheckMinMaxAvg(diag);
        CheckMinMaxAvg(areaBundle);

        // Check that the maximum aspect ratio of an individual cell is consistent
        // with the global minimums and maximums.
        Assert.True(S2.kMaxEdgeAspect >= 1);
        Assert.True(S2.kMaxEdgeAspect <=
            S2.kMaxEdge.Deriv/ S2.kMinEdge.Deriv);
        Assert.True(S2.kMaxDiagAspect >= 1);
        Assert.True(S2.kMaxDiagAspect <=
            S2.kMaxDiag.Deriv/ S2.kMinDiag.Deriv);

        // Check various conditions that are provable mathematically.
        CheckLessOrEqual(widthBundle, angle_span);
        CheckLessOrEqual(widthBundle, edge);
        CheckLessOrEqual(edge, diag);

        Assert.True(S2.kMinArea.Deriv>=
            S2.kMinWidth.Deriv* S2.kMinEdge.Deriv- S2.DoubleError);
        Assert.True(S2.kMaxArea.Deriv<=
            S2.kMaxWidth.Deriv* S2.kMaxEdge.Deriv+ S2.DoubleError);

        // The minimum level for which the minimum or maximum width of a cell is at
        // most 0 is kMaxCellLevel, because no cell at any level has width less than
        // or equal to zero.
        Assert.Equal(S2.kMinWidth.GetLevelForMaxValue(0), S2.kMaxCellLevel);
        Assert.Equal(S2.kMaxWidth.GetLevelForMaxValue(0), S2.kMaxCellLevel);

        // The maximum level for which the minimum or maximum width of a cell is at
        // least 4 is 0, because no cell at any level has width greater than 4.
        Assert.Equal(S2.kMinWidth.GetLevelForMinValue(4), 0);
        Assert.Equal(S2.kMaxWidth.GetLevelForMinValue(4), 0);

        // GetLevelForMaxValue() and friends have built-in assertions, we just need
        // to call these functions to test them.
        //
        // We don't actually check that the metrics are correct here, e.g. that
        // GetLevelForMaxValue(1) is a lower bound on the level of cells with width 1.
        // It is easier to check these properties in s2cell_test, since S2Cell has
        // methods to compute the cell vertices, etc.

        for (int level = -2; level <= S2.kMaxCellLevel + 3; ++level)
        {
            var width = S2.kMinWidth.Deriv* Math.Pow(2, -level);
            if (level >= S2.kMaxCellLevel + 3) width = 0;

            // Check boundary cases (exactly equal to a threshold value).
            int expected_level = Math.Max(0, Math.Min(S2.kMaxCellLevel, level));
            Assert.Equal(S2.kMinWidth.GetLevelForMaxValue(width), expected_level);
            Assert.Equal(S2.kMinWidth.GetLevelForMinValue(width), expected_level);
            Assert.Equal(S2.kMinWidth.GetClosestLevel(width), expected_level);

            // Also check non-boundary cases.
            Assert.Equal(S2.kMinWidth.GetLevelForMaxValue(1.2 * width), expected_level);
            Assert.Equal(S2.kMinWidth.GetLevelForMinValue(0.8 * width), expected_level);
            Assert.Equal(S2.kMinWidth.GetClosestLevel(1.2 * width), expected_level);
            Assert.Equal(S2.kMinWidth.GetClosestLevel(0.8 * width), expected_level);

            // Same thing for area.
            double area = S2.kMinArea.Deriv* Math.Pow(4, -level);
            if (level <= -3) area = 0;
            Assert.Equal(S2.kMinArea.GetLevelForMaxValue(area), expected_level);
            Assert.Equal(S2.kMinArea.GetLevelForMinValue(area), expected_level);
            Assert.Equal(S2.kMinArea.GetClosestLevel(area), expected_level);
            Assert.Equal(S2.kMinArea.GetLevelForMaxValue(1.2 * area), expected_level);
            Assert.Equal(S2.kMinArea.GetLevelForMinValue(0.8 * area), expected_level);
            Assert.Equal(S2.kMinArea.GetClosestLevel(1.2 * area), expected_level);
            Assert.Equal(S2.kMinArea.GetClosestLevel(0.8 * area), expected_level);
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
    // S2.kMinWidth.GetLevelForMinValue(width) to be slightly more readable than
    // than S2.kWidth.min().GetLevelForMinValue(width).  Also, there is no
    // fundamental reason that we need to analyze the minimum, maximum, and average
    // values of every metric; it would be perfectly reasonable to just define
    // one of these.
    private class MetricBundle
    {
        internal MetricBundle(Metric min, Metric max, Metric avg)
        { min_ = min; max_ = max; avg_ = avg; }

        internal Metric min_;
        internal Metric max_;
        internal Metric avg_;
    }
}
