namespace S2Geometry;

using S2MaxDistanceTargets = S2DistanceTarget<S2MaxDistance>;

public class S2MaxDistanceTargetsTests
{
    [Fact]
    public void Test_CellTarget_GetCapBound()
    {
        for (int i = 0; i < 100; i++)
        {
            S2Cell cell = new(S2Testing.GetRandomCellId());
            S2MaxDistanceCellTarget target = new(cell);
            var cap = target.GetCapBound();

            for (int j = 0; j < 100; j++)
            {
                S2Point p_test = S2Testing.RandomPoint();
                // Check points outside of cap to be away from S2MaxDistance.Zero().
                if (!(cap.Contains(p_test)))
                {
                    S1ChordAngle dist = cell.MaxDistance(p_test);
                    Assert.True(S2MaxDistance.Zero < new S2MaxDistance(dist));
                }
            }
        }
    }

    [Fact]
    public void Test_IndexTarget_GetCapBound()
    {
        MutableS2ShapeIndex index = new();

        S2Polygon polygon = new(new S2Cell(S2Testing.GetRandomCellId()));
        index.Add(new S2Polygon.Shape(polygon));

        S2Point p = S2Testing.RandomPoint();
        S2Point[] pts = { p };
        index.Add(new S2PointVectorShape(pts));

        S2MaxDistanceShapeIndexTarget target = new(index);
        var cap = target.GetCapBound();

        for (int j = 0; j < 100; j++)
        {
            var p_test = S2Testing.RandomPoint();
            // Check points outside of cap to be away from S2MaxDistance.Zero().
            if (!(cap.Contains(p_test)))
            {
                var cur_dist = S2MaxDistance.Infinity;
                Assert.True(target.UpdateMinDistance(p_test, ref cur_dist));
                Assert.True(S2MaxDistance.Zero < cur_dist);
            }
        }
    }

    [Fact]
    public void Test_PointTarget_UpdateMaxDistance()
    {
        S2MaxDistancePointTarget target = new(MakePointOrDie("0:0"));
        S2MaxDistance dist0 = new(S1ChordAngle.FromDegrees(0));
        S2MaxDistance dist10 = new(S1ChordAngle.FromDegrees(10));

        // Update max distance target to point.
        var p = MakePointOrDie("1:0");
        Assert.True(target.UpdateMinDistance(p, ref dist0));
        Assert2.Near(1.0, dist0.ToS1ChordAngle().Degrees(), S2.DoubleError);
        Assert.False(target.UpdateMinDistance(p, ref dist10));

        // Reset dist0 which was updated.
        dist0 = new(S1ChordAngle.FromDegrees(0));
        // Test for edges.
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist0));
        Assert2.Near(1.0, dist0.ToS1ChordAngle().Degrees(), S2.DoubleError);
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist10));

        // Reset dist0 which was updated.
        dist0 = new(S1ChordAngle.FromDegrees(0));
        // Test for cell.
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist0));
        // Leaf cell will be tiny compared to 10 degrees - expect no update.
        Assert.False(target.UpdateMinDistance(cell, ref dist10));
    }

    [Fact]
    public void Test_PointTarget_UpdateMaxDistanceToEdgeWhenEqual()
    {
        // Verifies that UpdateMinDistance only returns true when the new distance
        // is greater than the old distance (not greater than or equal to).
        S2MaxDistancePointTarget target = new(MakePointOrDie("1:0"));
        var dist = S2MaxDistance.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_PointTarget_UpdateMaxDistanceToCellWhenEqual()
    {
        S2MaxDistancePointTarget target = new(MakePointOrDie("1:0"));
        var dist = S2MaxDistance.Infinity;
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_EdgeTarget_UpdateMaxDistance()
    {
        var target_edge = ParsePointsOrDie("0:-1, 0:1");
        S2MaxDistanceEdgeTarget target = new(target_edge[0], target_edge[1]);
        S2MaxDistance dist0 = new(S1ChordAngle.FromDegrees(0));
        S2MaxDistance dist10 = new(S1ChordAngle.FromDegrees(10));

        // Update max distance target to point.
        S2Point p = MakePointOrDie("0:2");
        Assert.True(target.UpdateMinDistance(p, ref dist0));
        Assert2.Near(3.0, dist0.ToS1ChordAngle().Degrees(), S2.DoubleError);
        Assert.False(target.UpdateMinDistance(p, ref dist10));

        // Reset dist0 which was updated.
        dist0 = new(S1ChordAngle.FromDegrees(0));
        // Test for edges.
        var test_edge = ParsePointsOrDie("0:2, 0:3");
        Assert.True(target.UpdateMinDistance(test_edge[0], test_edge[1], ref dist0));
        Assert2.Near(4.0, dist0.ToS1ChordAngle().Degrees(), S2.DoubleError);
        Assert.False(target.UpdateMinDistance(test_edge[0], test_edge[1], ref dist10));

        // Reset dist0 which was updated.
        dist0 = new(S1ChordAngle.FromDegrees(0));
        // Test for cell.
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist0));
        // Leaf cell will be tiny compared to 10 degrees - expect no update.
        Assert.False(target.UpdateMinDistance(cell, ref dist10));
    }

    [Fact]
    public void Test_EdgeTarget_UpdateMaxDistanceToEdgeWhenEqual()
    {
        S2MaxDistanceEdgeTarget target = new(MakePointOrDie("1:0"), MakePointOrDie("1:1"));
        var dist = S2MaxDistance.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_EdgeTarget_UpdateMaxDistanceToEdgeAntipodal()
    {
        S2MaxDistanceEdgeTarget target = new(MakePointOrDie("0:89"), MakePointOrDie("0:91"));
        var dist = S2MaxDistance.Infinity;
        var edge = ParsePointsOrDie("1:-90, -1:-90");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.Equal(dist.ToS1ChordAngle(), S1ChordAngle.Straight);
    }

    [Fact]
    public void Test_EdgeTarget_UpdateMaxDistanceToCellWhenEqual()
    {
        S2MaxDistanceEdgeTarget target = new(MakePointOrDie("1:0"), MakePointOrDie("1:1"));
        var dist = S2MaxDistance.Infinity;
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_CellTarget_UpdateMaxDistance()
    {
        S2MaxDistanceCellTarget target = new(new S2Cell(new S2CellId(MakePointOrDie("0:1"))));
        S2MaxDistance dist0 = new(S1ChordAngle.FromDegrees(0));
        S2MaxDistance dist10 = new(S1ChordAngle.FromDegrees(10));

        // Update max distance target to point.
        var p = MakePointOrDie("0:0");
        Assert.True(target.UpdateMinDistance(p, ref dist0));
        Assert.False(target.UpdateMinDistance(p, ref dist10));

        // Reset dist0 which was updated.
        dist0 = new(S1ChordAngle.FromDegrees(0));
        // Test for edges.
        var test_edge = ParsePointsOrDie("0:2, 0:3");
        Assert.True(target.UpdateMinDistance(test_edge[0], test_edge[1], ref dist0));
        Assert.False(target.UpdateMinDistance(test_edge[0], test_edge[1], ref dist10));

        // Reset dist0 which was updated.
        dist0 = new(S1ChordAngle.FromDegrees(0));
        // Test for cell.
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist0));
        // Leaf cell extent will be tiny compared to 10 degrees - expect no update.
        Assert.False(target.UpdateMinDistance(cell, ref dist10));
    }

    [Fact]
    public void Test_CellTarget_UpdateMaxDistanceToEdgeWhenEqual()
    {
        S2MaxDistanceCellTarget target = new(new S2Cell(new S2CellId(MakePointOrDie("0:1"))));
        var dist = S2MaxDistance.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_CellTarget_UpdateMaxDistanceToCellWhenEqual()
    {
        S2MaxDistanceCellTarget target = new(new S2Cell(new S2CellId(MakePointOrDie("0:1"))));
        var dist = S2MaxDistance.Infinity;
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_ShapeIndexTarget_UpdateMaxDistanceToEdgeWhenEqual()
    {
        var target_index = MakeIndexOrDie("1:0 # #");
        S2MaxDistanceShapeIndexTarget target = new(target_index);
        var dist = S2MaxDistance.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_ShapeIndexTarget_UpdateMaxDistanceToCellWhenEqual()
    {
        var target_index = MakeIndexOrDie("1:0 # #");
        S2MaxDistanceShapeIndexTarget target = new(target_index);
        var dist = S2MaxDistance.Infinity;
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_CellTarget_UpdateMaxDistanceToCellAntipodal()
    {
        var p = MakePointOrDie("0:0");
        S2MaxDistanceCellTarget target = new(new S2Cell(p));
        var dist = S2MaxDistance.Infinity;
        S2Cell cell = new(-p);
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.Equal(dist.ToS1ChordAngle(), S1ChordAngle.Straight);
        // Expect a second update to do nothing.
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_PointTarget_VisitContainingShapes()
    {
        // Only shapes 2 and 4 should contain the target point.
        var index = MakeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
        S2Point p = MakePointOrDie("1:1");
        // Test against antipodal point.
        S2MaxDistancePointTarget target = new(-p);
        Assert.Equal((new int[] { 2 }), GetContainingShapes(target, index, 1));
        Assert.Equal((new int[] { 2, 4 }), GetContainingShapes(target, index, 5));
    }

    [Fact]
    public void Test_EdgeTarget_VisitContainingShapes()
    {
        // Only shapes 2 and 4 should contain the target edge.
        var index = MakeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
        // Test against antipodal edge.
        var edge = ParsePointsOrDie("1:2, 2:1");
        S2MaxDistanceEdgeTarget target = new(-edge[0], -edge[1]);
        Assert.Equal((new int[] { 2 }), GetContainingShapes(target, index, 1));
        Assert.Equal((new int[] { 2, 4 }), GetContainingShapes(target, index, 5));
    }

    [Fact]
    public void Test_CellTarget_VisitContainingShapes()
    {
        var index = MakeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
        // Only shapes 2 and 4 should contain a very small cell near
        // the antipode of 1:1.
        S2CellId cellid1 = new(-MakePointOrDie("1:1"));
        S2MaxDistanceCellTarget target1 = new(new S2Cell(cellid1));
        Assert.Equal((new int[] { 2 }), GetContainingShapes(target1, index, 1));
        Assert.Equal((new int[] { 2, 4 }), GetContainingShapes(target1, index, 5));

        // For a larger antipodal cell that properly contains one or more index
        // cells, all shapes that intersect the first such cell in S2CellId order are
        // returned.  In the test below, this happens to again be the 1st and 3rd
        // polygons (whose shape_ids are 2 and 4).
        var cellid2 = cellid1.Parent(5);
        S2MaxDistanceCellTarget target2 = new(new S2Cell(cellid2));
        Assert.Equal((new int[] { 2, 4 }), GetContainingShapes(target2, index, 5));
    }

    [Fact]
    public void Test_ShapeIndexTarget_VisitContainingShapes()
    {
        // Create an index containing a repeated grouping of one point, one
        // polyline, and one polygon.
        var index = MakeIndexOrDie("1:1 | 4:4 | 7:7 | 10:10 # " +
            "1:1, 1:2 | 4:4, 4:5 | 7:7, 7:8 | 10:10, 10:11 # " +
            "0:0, 0:3, 3:0 | 3:3, 3:6, 6:3 | 6:6, 6:9, 9:6 | 9:9, 9:12, 12:9");

        // Construct a target consisting of one point, one polyline, and one polygon
        // with two loops where only the second loop is contained by a polygon in
        // the index above.
        MutableS2ShapeIndex target_index = new();

        var pts = Reflect(ParsePointsOrDie("1:1")).ToArray();
        target_index.Add(new S2PointVectorShape(pts));

        S2Polyline line = new(Reflect(ParsePointsOrDie("4:5, 5:4")).ToArray());
        target_index.Add(new S2Polyline.Shape(line));

        var loop1 = Reflect(ParsePointsOrDie("20:20, 20:21, 21:20"));
        var loop2 = Reflect(ParsePointsOrDie("10:10, 10:11, 11:10"));
        List<List<S2Point>> loops = new(){ loop1, loop2 };
        target_index.Add(new S2LaxPolygonShape(loops));

        S2MaxDistanceShapeIndexTarget target = new(target_index);
        // These are the shape_ids of the 1st, 2nd, and 4th polygons of "index"
        // (noting that the 4 points are represented by one S2PointVectorShape).
        Assert.Equal((new int[] { 5, 6, 8 }), GetContainingShapes(target, index, 5));
    }

    [Fact]
    public void Test_ShapeIndexTarget_VisitContainingShapesEmptyAndFull()
    {
        // Verify that VisitContainingShapes never returns empty polygons and always
        // returns full polygons (i.e., those containing the entire sphere).

        // Creating an index containing one empty and one full polygon.
        var index = MakeIndexOrDie("# # empty | full");

        // Check only the full polygon is returned for a point target.
        var point_index = MakeIndexOrDie("1:1 # #");
        S2MaxDistanceShapeIndexTarget point_target = new(point_index);
        Assert.Equal((new int[] { 1 }), GetContainingShapes(point_target, index, 5));

        // Check only the full polygon is returned for a full polygon target.
        var full_polygon_index = MakeIndexOrDie("# # full");
        S2MaxDistanceShapeIndexTarget full_target = new(full_polygon_index);
        Assert.Equal((new int[] { 1 }), GetContainingShapes(full_target, index, 5));

        // Check that nothing is returned for an empty polygon target.  (An empty
        // polygon has no connected components and does not intersect anything, so
        // according to the API of GetContainingShapes nothing should be returned.)
        var empty_polygon_index = MakeIndexOrDie("# # empty");
        S2MaxDistanceShapeIndexTarget empty_target = new(empty_polygon_index);
        Assert.Equal((Array.Empty<int>()), GetContainingShapes(empty_target, index, 5));
    }

    // Negates S2 points to reflect them through the sphere.
    private static List<S2Point> Reflect(List<S2Point> pts)
    {
        List<S2Point> negative_pts = new();
        foreach (var p in pts)
        {
            negative_pts.Add(-p);
        }
        return negative_pts;
    }

    private static int[] GetContainingShapes(S2MaxDistanceTargets target, S2ShapeIndex index, int max_shapes)
    {
        SortedSet<Int32> shape_ids = new();
        target.VisitContainingShapes(
            index, (S2Shape containing_shape, S2Point target_point) =>
            {
                shape_ids.Add(containing_shape.Id);
                return shape_ids.Count < max_shapes;
            });
        return shape_ids.ToArray();
    }
}

public static class S2MaxDistanceExtensions
{
    public static S1ChordAngle ToS1ChordAngle(this S2MaxDistance me) => me.Distance;
}
