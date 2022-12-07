namespace S2Geometry;

public class S2PaddedCellTests
{
    [Fact]
    internal void Test_S2PaddedCell_S2CellMethods()
    {
        // Test the S2PaddedCell methods that have approximate S2Cell equivalents.
        int kIters = 1000;
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2CellId id = S2Testing.GetRandomCellId();
            double padding = Math.Pow(1e-15, S2Testing.Random.RandDouble());
            S2Cell cell = new(id);
            S2PaddedCell pcell = new(id, padding);
            CompareS2CellToPadded(cell, pcell, padding);
            if (!id.IsLeaf())
            {
                var children = new S2Cell[4];
                Assert.True(cell.Subdivide(children));
                for (int pos = 0; pos < 4; ++pos)
                {
                    pcell.GetChildIJ(pos, out var i, out var j);
                    CompareS2CellToPadded(children[pos], new S2PaddedCell(pcell, i, j), padding);
                }
            }
        }
    }

    [Fact]
    internal void Test_S2PaddedCell_GetEntryExitVertices()
    {
        int kIters = 1000;
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2CellId id = S2Testing.GetRandomCellId();
            // Check that entry/exit vertices do not depend on padding.
            Assert.Equal(new S2PaddedCell(id, 0).GetEntryVertex(),
                      new S2PaddedCell(id, 0.5).GetEntryVertex());
            Assert.Equal(new S2PaddedCell(id, 0).GetExitVertex(),
                      new S2PaddedCell(id, 0.5).GetExitVertex());

            // Check that the exit vertex of one cell is the same as the entry vertex
            // of the immediately following cell.  (This also tests wrapping from the
            // end to the start of the S2CellId curve with high probability.)
            Assert.Equal(new S2PaddedCell(id, 0).GetExitVertex(),
                      new S2PaddedCell(id.NextWrap(), 0).GetEntryVertex());

            // Check that the entry vertex of a cell is the same as the entry vertex
            // of its first child, and similarly for the exit vertex.
            if (!id.IsLeaf())
            {
                Assert.Equal(new S2PaddedCell(id, 0).GetEntryVertex(),
                          new S2PaddedCell(id.Child(0), 0).GetEntryVertex());
                Assert.Equal(new S2PaddedCell(id, 0).GetExitVertex(),
                          new S2PaddedCell(id.Child(3), 0).GetExitVertex());
            }
        }
    }

    [Fact]
    internal void Test_S2PaddedCell_ShrinkToFit()
    {
        const int kIters = 1000;
        for (int iter = 0; iter < kIters; ++iter)
        {
            // Start with the desired result and work backwards.
            S2CellId result = S2Testing.GetRandomCellId();
            R2Rect result_uv = result.BoundUV();
            R2Point size_uv = result_uv.GetSize();

            // Find the biggest rectangle that fits in "result" after padding.
            // (These calculations ignore numerical errors.)
            double max_padding = 0.5 * Math.Min(size_uv[0], size_uv[1]);
            double padding = max_padding * S2Testing.Random.RandDouble();
            R2Rect max_rect = result_uv.Expanded(-padding);

            // Start with a random subset of the maximum rectangle.
            R2Point a = new(SampleInterval(max_rect[0]), SampleInterval(max_rect[1]));
            R2Point b = new(SampleInterval(max_rect[0]), SampleInterval(max_rect[1]));
            if (!result.IsLeaf())
            {
                // If the result is not a leaf cell, we must ensure that no child of
                // "result" also satisfies the conditions of ShrinkToFit().  We do this
                // by ensuring that "rect" intersects at least two children of "result"
                // (after padding).
                int axis = S2Testing.Random.Uniform(2);
                double center = result.CenterUV()[axis];

                // Find the range of coordinates that are shared between child cells
                // along that axis.
                R1Interval shared = new(center - padding, center + padding);
                double mid = SampleInterval(shared.Intersection(max_rect[axis]));
                a = a.SetAxis(axis, SampleInterval(new R1Interval(max_rect[axis].Lo, mid)));
                b = b.SetAxis(axis, SampleInterval(new R1Interval(mid, max_rect[axis].Hi)));
            }
            R2Rect rect = R2Rect.FromPointPair(a, b);

            // Choose an arbitrary ancestor as the S2PaddedCell.
            S2CellId initial_id = result.Parent(S2Testing.Random.Uniform(result.Level() + 1));
            Assert.Equal(result, new S2PaddedCell(initial_id, padding).ShrinkToFit(rect));
        }
    }

    private static void CompareS2CellToPadded(S2Cell cell, S2PaddedCell pcell, double padding)
    {
        Assert.Equal(cell.Id, pcell.Id);
        Assert.Equal(cell.Level, pcell.Level);
        Assert.Equal(padding, pcell.Padding);
        Assert.Equal(cell.BoundUV.Expanded(padding), pcell.Bound);
        var center_uv = cell.Id.CenterUV();
        Assert.Equal(R2Rect.FromPoint(center_uv).Expanded(padding), pcell.Middle);
        Assert.Equal(cell.Center(), pcell.GetCenter());
    }

    private static double SampleInterval(R1Interval x)
        => S2Testing.Random.UniformDouble(x.Lo, x.Hi);
}
