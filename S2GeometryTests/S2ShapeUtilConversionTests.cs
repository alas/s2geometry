namespace S2Geometry;

using static S2ShapeUtil;
using Loop = List<S2Point>;

public class S2ShapeUtilConversionTests
{
    // ShapeToS2Points tests

    [Fact]
    internal void Test_S2ShapeConversionUtilTest_PointVectorShapeToPoints()
    {
        var points = ParsePointsOrDie("11:11, 10:0, 5:5");
        S2PointVectorShape point_vector = new(points.ToArray());
        var extract = ShapeToS2Points(point_vector);
        // TODO(user,b/205813109): Use gmock ASSERT_THAT.
        Assert.Equal(extract.Count, 3);
        for (int i = 0; i < extract.Count; i++)
        {
            Assert.Equal(extract[i], points[i]);
        }
    }

    // ShapeToS2Polyline tests

    [Fact]
    internal void Test_SS2ShapeConversionUtilTest_LineToS2Polyline()
    {
        var points = ParsePointsOrDie("11:11, 10:0, 5:5");
        S2LaxPolylineShape lax_polyline = new(points.ToArray());
        var polyline = ShapeToS2Polyline(lax_polyline);
        Assert.Equal(polyline.NumVertices(), 3);
        for (int i = 0; i < polyline.NumVertices(); i++)
        {
            Assert.Equal(polyline.Vertex(i), points[i]);
        }
    }

    [Fact]
    internal void Test_S2ShapeConversionUtilTest_ClosedLineToS2Polyline()
    {
        var points = ParsePointsOrDie("0:0, 0:10, 10:10, 0:0");
        S2LaxPolylineShape lax_polyline = new(points.ToArray());
        var polyline = ShapeToS2Polyline(lax_polyline);
        Assert.Equal(polyline.NumVertices(), 4);
        for (int i = 0; i < polyline.NumVertices(); i++)
        {
            Assert.Equal(polyline.Vertex(i), points[i]);
        }
    }

    // ShapeToS2Polygon tests

    // Creates a (lax) polygon shape from the provided loops, and ensures that the
    // S2Polygon produced by ShapeToS2Polygon represents the same polygon.
    private static void VerifyShapeToS2Polygon(List<Loop> loops,
        int expected_num_loops, int expected_num_vertices)
    {
        S2LaxPolygonShape lax_polygon = new(loops);
        var polygon = ShapeToS2Polygon(lax_polygon);

        Assert.Equal(polygon.NumLoops(), expected_num_loops);
        Assert.Equal(polygon.NumVertices, expected_num_vertices);
        for (int i = 0; i < polygon.NumLoops(); i++)
        {
            var loop = polygon.Loop(i);
            for (int j = 0; j < loop.NumVertices; j++)
            {
                Assert.Equal(loop.OrientedVertex(j), loops[i][j]);
            }
        }
    }

    [Fact]
    internal void Test_S2ShapeConversionUtilTest_PolygonWithHoleToS2Polygon()
    {
        // a polygon with one shell and one hole
        Loop shell = new(ParsePointsOrDie("0:0, 0:10, 10:10, 10:0"));
        Loop hole = new(ParsePointsOrDie("4:4, 6:4, 6:6, 4:6"));
        List<Loop> loops = new() { shell, hole };

        VerifyShapeToS2Polygon(loops, 2, 8);
    }

    [Fact]
    internal void Test_S2ShapeConversionUtilTest_MultiPolygonToS2Polygon()
    {
        // a polygon with multiple shells
        Loop shell1 = new(ParsePointsOrDie("0:0, 0:2, 2:2, 2:0"));
        Loop shell2 = new(ParsePointsOrDie("0:4, 0:6, 3:6"));
        List<Loop> loops = new() { shell1, shell2 };

        VerifyShapeToS2Polygon(loops, 2, 7);
    }

    [Fact]
    internal void Test_S2ShapeConversionUtilTest_TwoHolesToS2Polygon()
    {
        // a polygon shell with two holes
        Loop shell = new(ParsePointsOrDie("0:0, 0:10, 10:10, 10:0"));
        Loop hole1 = new(ParsePointsOrDie("1:1, 3:3, 1:3"));
        Loop hole2 = new(ParsePointsOrDie("2:6, 4:7, 2:8"));
        List<Loop> loops = new() { shell, hole1, hole2 };

        VerifyShapeToS2Polygon(loops, 3, 10);
    }

    [Fact]
    internal void Test_S2ShapeConversionUtilTest_FullPolygonToS2Polygon()
    {
        S2Point kFullVertex = new(0, 0, -1);
        var q1 = new S2Point[] { kFullVertex };
        S2Loop kFull = new(q1);
        var q2 = kFull.CloneVertices();
        Loop q3 = new(q2);
        // verify that a full polygon is converted correctly
        Loop loop1 = new(S2Loop.kFull.CloneVertices());
        List<Loop> loops = new() { loop1 };

        var lax_polygon = MakeLaxPolygonOrDie("full");
        var polygon = ShapeToS2Polygon(lax_polygon);
        Assert.Equal(polygon.NumLoops(), 1);
        Assert.Equal(polygon.NumVertices, 1);
        Assert.True(polygon.IsFull());
        for (int i = 0; i < polygon.NumLoops(); i++)
        {
            var loop = polygon.Loop(i);
            for (int j = 0; j < loop.NumVertices; j++)
            {
                Assert.Equal(loop.OrientedVertex(j), loops[i][j]);
            }
        }
    }
}
