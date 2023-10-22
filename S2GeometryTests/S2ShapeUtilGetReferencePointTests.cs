namespace S2Geometry;

public class S2ShapeUtilGetReferencePointTests
{
    [Fact]
    internal void Test_GetReferencePoint_EmptyPolygon()
    {
        S2LaxPolygonShape shape = new(new S2Polygon());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_GetReferencePoint_FullPolygon()
    {
        S2LaxPolygonShape shape = new(new S2Polygon(MakeLoopOrDie("full")));
        Assert.True(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_GetReferencePoint_DegenerateLoops()
    {
        List<List<S2Point>> loops = [
            ParsePointsOrDie("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
            ParsePointsOrDie("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
            ParsePointsOrDie("5:5, 6:6")
        ];
        S2LaxPolygonShape shape = new(loops);
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_GetReferencePoint_InvertedLoops()
    {
        List<List<S2Point>> loops = [
            ParsePointsOrDie("1:2, 1:1, 2:2"),
            ParsePointsOrDie("3:4, 3:3, 4:4")
        ];
        var shape = new S2LaxPolygonShape(loops);
        Assert.True(shape.ContainsBruteForce(S2.Origin));
    }

    [Fact]
    internal void Test_GetReferencePoint_PartiallyDegenerateLoops()
    {
        for (var iter = 0; iter < 100; ++iter)
        {
            // First we construct a long convoluted edge chain that follows the
            // S2CellId Hilbert curve.  At some random point along the curve, we
            // insert a small triangular loop.
            List<List<S2Point>> loops = new(1) { new() };
            var loop = loops[0];
            int num_vertices = 100;
            var start = S2Testing.GetRandomCellId(S2.kMaxCellLevel - 1);
            var end = start.AdvanceWrap(num_vertices);
            var loop_cellid = start.AdvanceWrap(
                S2Testing.Random.Uniform(num_vertices - 2) + 1);
            var triangle = new List<S2Point>();
            for (var cellid = start; cellid != end; cellid = cellid.NextWrap())
            {
                if (cellid == loop_cellid)
                {
                    // Insert a small triangular loop.  We save the loop so that we can
                    // test whether it contains the origin later.
                    triangle.Add(cellid.Child(0).ToPoint());
                    triangle.Add(cellid.Child(1).ToPoint());
                    triangle.Add(cellid.Child(2).ToPoint());
                    loop.AddRange(triangle);
                    loop.Add(cellid.Child(0).ToPoint());
                }
                else
                {
                    loop.Add(cellid.ToPoint());
                }
            }
            // Now we retrace our steps, except that we skip the three edges that form
            // the triangular loop above.
            for (S2CellId cellid = end; cellid != start; cellid = cellid.PrevWrap())
            {
                if (cellid == loop_cellid)
                {
                    loop.Add(cellid.Child(0).ToPoint());
                }
                else
                {
                    loop.Add(cellid.ToPoint());
                }
            }
            S2LaxPolygonShape shape = new(loops);
            S2Loop triangle_loop = new(triangle);
            var rp = shape.GetReferencePoint();
            Assert.Equal(triangle_loop.Contains(rp.Point), rp.Contained);
        }
    }
}
