namespace S2Geometry;

using S2MinDistanceTarget = S2DistanceTarget<S1ChordAngle>;

public class S2MinDistanceTargetsTests
{
    [Fact]
    public void Test_PointTarget_UpdateMinDistanceToEdgeWhenEqual() {
        // Verifies that UpdateMinDistance only returns true when the new distance
        // is less than the old distance (not less than or equal to).
        var target = new S2MinDistancePointTarget(MakePointOrDie("1:0"));
        S1ChordAngle dist = S1ChordAngle.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_PointTarget_UpdateMinDistanceToCellWhenEqual() {
        // Verifies that UpdateMinDistance only returns true when the new distance
        // is less than the old distance (not less than or equal to).
        var target = new S2MinDistancePointTarget(MakePointOrDie("1:0"));
        var dist = S1ChordAngle.Infinity;
        S2Cell cell = new(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_EdgeTarget_UpdateMinDistanceToEdgeWhenEqual() {
        var target = new S2MinDistanceEdgeTarget(
            MakePointOrDie("1:0"), MakePointOrDie("1:1"));
        var dist = S1ChordAngle.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_EdgeTarget_UpdateMinDistanceToCellWhenEqual() {
        var target = new S2MinDistanceEdgeTarget(
            MakePointOrDie("1:0"), MakePointOrDie("1:1"));
        var dist = S1ChordAngle.Infinity;
        var cell = new S2Cell(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_CellTarget_UpdateMinDistanceToEdgeWhenEqual() {
        var target = new S2MinDistanceCellTarget(new S2Cell(
            new S2CellId(MakePointOrDie("0:1"))));
        var dist = S1ChordAngle.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_CellTarget_UpdateMinDistanceToCellWhenEqual() {
        var target = new S2MinDistanceCellTarget(new S2Cell(
            new S2CellId(MakePointOrDie("0:1"))));
        var dist = S1ChordAngle.Infinity;
        var cell = new S2Cell(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_CellUnionTarget_UpdateMinDistanceToEdgeWhenEqual() {
        var target = new S2MinDistanceCellUnionTarget(new S2CellUnion(
            new List<S2CellId>{new S2CellId(MakePointOrDie("0:1"))}));
        var dist = S1ChordAngle.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_CellUnionTarget_UpdateMinDistanceToCellWhenEqual()
    {
        var target = new S2MinDistanceCellUnionTarget(new S2CellUnion(
            new List<S2CellId>{new S2CellId(MakePointOrDie("0:1"))}));
        var dist = S1ChordAngle.Infinity;
        var cell = new S2Cell(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_ShapeIndexTarget_UpdateMinDistanceToEdgeWhenEqual() {
        var target_index = MakeIndexOrDie("1:0 # #");
        var target = new S2MinDistanceShapeIndexTarget(target_index);
        var dist = S1ChordAngle.Infinity;
        var edge = ParsePointsOrDie("0:-1, 0:1");
        Assert.True(target.UpdateMinDistance(edge[0], edge[1], ref dist));
        Assert.False(target.UpdateMinDistance(edge[0], edge[1], ref dist));
    }

    [Fact]
    public void Test_ShapeIndexTarget_UpdateMinDistanceToCellWhenEqual() {
        var target_index = MakeIndexOrDie("1:0 # #");
        var target = new S2MinDistanceShapeIndexTarget(target_index);
        var dist = S1ChordAngle.Infinity;
        var cell = new S2Cell(new S2CellId(MakePointOrDie("0:0")));
        Assert.True(target.UpdateMinDistance(cell, ref dist));
        Assert.False(target.UpdateMinDistance(cell, ref dist));
    }

    [Fact]
    public void Test_PointTarget_VisitContainingShapes() {
        // Only shapes 2 and 4 should contain the target point.
        var index = MakeIndexOrDie(
            "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
        var target = new S2MinDistancePointTarget(MakePointOrDie("1:1"));
        Assert.True(IsSubsetOfSize(GetContainingShapes(target, index, 1),
                                   new int[]{ 2, 4 }, 1));
        Assert.Equal((new int[]{ 2, 4}), GetContainingShapes(target, index, 5));
    }

    [Fact]
    public void Test_EdgeTarget_VisitContainingShapes() {
        // Only shapes 2 and 4 should contain the target point.
        var index = MakeIndexOrDie(
            "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
        var target = new S2MinDistanceEdgeTarget(MakePointOrDie("1:2"), MakePointOrDie("2:1"));
        Assert.True(IsSubsetOfSize(GetContainingShapes(target, index, 1),
                                   new int[]{ 2, 4}, 1));
        Assert.Equal((new int[]{ 2, 4}), GetContainingShapes(target, index, 5));
    }

    [Fact]
    public void Test_CellTarget_VisitContainingShapes() {
        var index = MakeIndexOrDie(
            "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
        // Only shapes 2 and 4 should contain a very small cell near 1:1.
        var cellid1 = new S2CellId(MakePointOrDie("1:1"));
        var target1 = new S2MinDistanceCellTarget(new S2Cell(cellid1));
        Assert.True(IsSubsetOfSize(GetContainingShapes(target1, index, 1),
                                   new int[]{ 2, 4}, 1));
        Assert.Equal(new int[]{ 2, 4}, GetContainingShapes(target1, index, 5));

        // For a larger cell that properly contains one or more index cells, all
        // shapes that intersect the first such cell in S2CellId order are returned.
        // In the test below, this happens to again be the 1st and 3rd polygons
        // (whose shape_ids are 2 and 4).
        var cellid2 = cellid1.Parent(5);
        var target2 = new S2MinDistanceCellTarget(new S2Cell(cellid2));
        Assert.Equal(new int[]{ 2, 4}, GetContainingShapes(target2, index, 5));
    }

    [Fact]
    public void Test_CellUnionTarget_VisitContainingShapes() {
        var index = MakeIndexOrDie(
            "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
        // Shapes 2 and 4 contain the leaf cell near 1:1, while shape 3 contains the
        // leaf cell near 7:7.
        var cellid1 = new S2CellId(MakePointOrDie("1:1"));
        var cellid2 = new S2CellId(MakePointOrDie("7:7"));
        var target1 = new S2MinDistanceCellUnionTarget(new S2CellUnion(
            new List<S2CellId>{ cellid1, cellid2}));
        Assert.True(IsSubsetOfSize(GetContainingShapes(target1, index, 1),
                                   new int[]{ 2, 3, 4}, 1));
        Assert.Equal((new int[]{ 2, 3, 4}), GetContainingShapes(target1, index, 5));
    }

    [Fact]
    public void Test_ShapeIndexTarget_VisitContainingShapes() {
        // Create an index containing a repeated grouping of one point, one
        // polyline, and one polygon.
        var index = MakeIndexOrDie("1:1 | 4:4 | 7:7 | 10:10 # " +
            "1:1, 1:2 | 4:4, 4:5 | 7:7, 7:8 | 10:10, 10:11 # " +
            "0:0, 0:3, 3:0 | 3:3, 3:6, 6:3 | 6:6, 6:9, 9:6 | 9:9, 9:12, 12:9");

        // Construct a target consisting of one point, one polyline, and one polygon
        // with two loops where only the second loop is contained by a polygon in
        // the index above.
        var target_index = MakeIndexOrDie(
            "1:1 # 4:5, 5:4 # 20:20, 20:21, 21:20; 10:10, 10:11, 11:10");

        var target = new S2MinDistanceShapeIndexTarget(target_index);
        // These are the shape_ids of the 1st, 2nd, and 4th polygons of "index"
        // (noting that the 4 points are represented by one S2PointVectorShape).
        Assert.Equal((new int[]{ 5, 6, 8}), GetContainingShapes(target, index, 5));
    }

    [Fact]
    public void Test_ShapeIndexTarget_VisitContainingShapesEmptyAndFull() {
        // Verify that VisitContainingShapes never returns empty polygons and always
        // returns full polygons (i.e., those containing the entire sphere).

        // Creating an index containing one empty and one full polygon.
        var index = MakeIndexOrDie("# # empty | full");

        // Check only the full polygon is returned for a point target.
        var point_index = MakeIndexOrDie("1:1 # #");
        var point_target = new S2MinDistanceShapeIndexTarget(point_index);
        Assert.Equal((new int[]{ 1}), GetContainingShapes(point_target, index, 5));

        // Check only the full polygon is returned for a full polygon target.
        var full_polygon_index = MakeIndexOrDie("# # full");
        var full_target = new S2MinDistanceShapeIndexTarget(full_polygon_index);
        Assert.Equal((new int[]{ 1}), GetContainingShapes(full_target, index, 5));

        // Check that nothing is returned for an empty polygon target.  (An empty
        // polygon has no connected components and does not intersect anything, so
        // according to the API of GetContainingShapes nothing should be returned.)
        var empty_polygon_index = MakeIndexOrDie("# # empty");
        var empty_target = new S2MinDistanceShapeIndexTarget(empty_polygon_index);
        Assert.Equal((new int[]{ }), GetContainingShapes(empty_target, index, 5));
    }

    private int[] GetContainingShapes(S2MinDistanceTarget target, S2ShapeIndex index, int max_shapes)
    {
        var shape_ids = new SortedSet<Int32>();
        target.VisitContainingShapes(
            index, (S2Shape containing_shape, S2Point target_point) => {
                shape_ids.Add(containing_shape.Id);
                return shape_ids.Count < max_shapes;
            });
        return shape_ids.ToArray();
    }

    // Given two sorted vectors "x" and "y", returns true if x is a subset of y
    // and x.size() == x_size.
    private bool IsSubsetOfSize(int[] x, int[] y, int x_size) {
        if (x.Length != x_size) return false;
        return y.IsSubsetOf(x);
    }
}
