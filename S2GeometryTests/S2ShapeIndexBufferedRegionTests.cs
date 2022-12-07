namespace S2Geometry;

public class S2ShapeIndexBufferedRegionTests
{
    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_EmptyIndex()
    {
        // Test buffering an empty S2ShapeIndex.
        var index = new MutableS2ShapeIndex();
        var radius = new S1ChordAngle(S1Angle.FromDegrees(2));
        var region = new S2ShapeIndexBufferedRegion(index, radius);
        var coverer = new S2RegionCoverer();
        S2CellUnion covering = coverer.GetCovering(region);
        Assert.True(covering.IsEmpty());
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_InitEmptyIndex()
    {
        // As above, but with Init().  This is mainly to prevent Init() from being
        // detected as dead code.
        MutableS2ShapeIndex index=new();
        S1ChordAngle radius=new(S1Angle.FromDegrees(2));
        S2ShapeIndexBufferedRegion region=new(index, radius);
        S2RegionCoverer coverer=new();
        S2CellUnion covering = coverer.GetCovering(region);
        Assert.True(covering.IsEmpty());
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_FullPolygon()
    {
        // Test buffering an S2ShapeIndex that contains a full polygon.
        var index = MakeIndexOrDie("# # full");
        var radius = new S1ChordAngle(S1Angle.FromDegrees(2));
        var region = new S2ShapeIndexBufferedRegion(index, radius);
        var coverer = new S2RegionCoverer();
        S2CellUnion covering = coverer.GetCovering(region);
        Assert.Equal(6, covering.Size());
        foreach (S2CellId id in covering)
        {
            Assert.True(id.IsFace());
        }
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_FullAfterBuffering()
    {
        // Test a region that becomes the full polygon after buffering.
        var index = MakeIndexOrDie("0:0 | 0:90 | 0:180 | 0:-90 | 90:0 | -90:0 # #");
        var radius = new S1ChordAngle(S1Angle.FromDegrees(60));
        var region = new S2ShapeIndexBufferedRegion(index, radius);
        var coverer = new S2RegionCoverer();
        coverer.Options_.MaxCells = 1000;
        S2CellUnion covering = coverer.GetCovering(region);
        Assert.Equal(6, covering.Size());
        foreach (S2CellId id in covering)
        {
            Assert.True(id.IsFace());
        }
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_PointZeroRadius()
    {
        // Test that buffering a point using a zero radius produces a non-empty
        // covering.  (This requires using "less than or equal to" distance tests.)
        var index = MakeIndexOrDie("34:25 # #");
        var region = new S2ShapeIndexBufferedRegion(index, S1ChordAngle.Zero);
        var coverer = new S2RegionCoverer();
        S2CellUnion covering = coverer.GetCovering(region);
        Assert.Equal(1, covering.Size());
        foreach (S2CellId id in covering)
        {
            Assert.True(id.IsLeaf());
        }
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_BufferedPointVsCap()
    {
        // Compute an S2Cell covering of a buffered S2Point, then make sure that the
        // covering is equivalent to the corresponding S2Cap.
        var index = MakeIndexOrDie("3:5 # #");
        S2Point point = MakePointOrDie("3:5");
        var radius = new S1ChordAngle(S1Angle.FromDegrees(2));
        var region = new S2ShapeIndexBufferedRegion(index, radius);
        S2RegionCoverer coverer = new();
        coverer.Options_.MaxCells = 50;
        S2CellUnion covering = coverer.GetCovering(region);
        S2Cap equivalent_cap = new(point, radius);
        S2Testing.CheckCovering(equivalent_cap, covering, true);
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_PointSet()
    {
        // Test buffering a set of points.
        var coverer = new S2RegionCoverer();
        coverer.Options_.MaxCells = 100;
        TestBufferIndex("10:20 | 10:23 | 10:26 # #", S1Angle.FromDegrees(5), coverer);
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_Polyline()
    {
        // Test buffering a polyline.
        S2RegionCoverer coverer = new();
        coverer.Options_.MaxCells = (100);
        TestBufferIndex("# 10:5, 20:30, -10:60, -60:100 #",
                        S1Angle.FromDegrees(2), coverer);
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_PolygonWithHole()
    {
        // Test buffering a polygon with a hole.
        S2RegionCoverer coverer = new();
        coverer.Options_.MaxCells = (100);
        TestBufferIndex("# # 10:10, 10:100, 70:0; 11:11, 69:0, 11:99",
                        S1Angle.FromDegrees(2), coverer);
    }

    [Fact]
    internal void Test_S2ShapeIndexBufferedRegion_HugeBufferRadius()
    {
        // Test buffering a set of points.
        S2RegionCoverer coverer = new();
        coverer.Options_.MaxCells = 100;
        TestBufferIndex("10:20 # #", S1Angle.FromDegrees(200), coverer);
    }

    // Verifies that an arbitrary S2ShapeIndex is buffered correctly, by first
    // converting the covering to an S2Polygon and then checking that (a) the
    // S2Polygon contains the original geometry and (b) the distance between the
    // original geometry and the boundary of the S2Polygon is at least "radius".
    //
    // The "radius" parameter is an S1Angle for convenience.
    // TODO(ericv): Add Degrees, Radians, etc, methods to S1ChordAngle?
    private static void TestBufferIndex(string index_str, S1Angle radius_angle, S2RegionCoverer coverer)
    {
        var index = MakeIndexOrDie(index_str);
        S1ChordAngle radius = new(radius_angle);
        var region = new S2ShapeIndexBufferedRegion(index, radius);
        S2CellUnion covering = coverer.GetCovering(region);
        // Compute an S2Polygon representing the union of the cells in the covering.
        S2Polygon covering_polygon = new();
        covering_polygon.InitToCellUnionBorder(covering);
        MutableS2ShapeIndex covering_index = new();
        covering_index.Add(new S2Polygon.Shape(covering_polygon));

        // (a) Check that the covering contains the original index.
        Assert.True(S2BooleanOperation.Contains(covering_index, index));

        // (b) Check that the distance between the boundary of the covering and the
        // the original indexed geometry is at least "radius".
        S2ClosestEdgeQuery query = new(covering_index);
        query.Options_.IncludeInteriors = (false);
        var target = new S2ClosestEdgeQuery.ShapeIndexTarget(index);
        Assert.False(query.IsDistanceLess(target, radius));
    }
}
