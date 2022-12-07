namespace S2Geometry;

using DirectedResult = S2HausdorffDistanceQuery.DirectedResult;
using Result = S2HausdorffDistanceQuery.Result;
using Options = S2HausdorffDistanceQuery.Options;

public class S2HausdorffDistanceQueryTests
{
    // Test the constructors and accessors of Result and DirectedResult.
    [Fact]
    internal void Test_S2HausdorffDistanceQueryTest_ResultConstructorsAndAccessorsWork()
    {
        S2Point point1 = S2LatLng.FromDegrees(3, 4).ToPoint();
        S2Point point2 = S2LatLng.FromDegrees(5, 6).ToPoint();
        S1ChordAngle distance1 = S1ChordAngle.FromDegrees(5);
        S1ChordAngle distance2 = S1ChordAngle.FromDegrees(5);

        DirectedResult directed_result1 = new(distance1, point1);
        DirectedResult directed_result2 = new(distance2, point2);
        Result result12 = new(directed_result1, directed_result2);

        Assert.Equal(directed_result1.TargetPoint, point1);
        Assert.Equal(directed_result1.Distance, distance1);
        Assert.Equal(directed_result2.TargetPoint, point2);
        Assert.Equal(directed_result2.Distance, distance2);

        Assert.Equal(result12.TargetToSource.TargetPoint, point1);
        Assert.Equal(result12.SourceToTarget.TargetPoint, point2);
        Assert.Equal(result12.GetDistance(), directed_result2.Distance);
    }

    // Test the constructors and accessors of the Options.
    [Fact]
    internal void Test_S2HausdorffDistanceQueryTest_OptionsConstructorsAndAccessorsWork()
    {
        Options default_options = new();
        Options options = new()
        {
            IncludeInteriors = !default_options.IncludeInteriors
        };

        Assert.True(default_options.IncludeInteriors);
        Assert.False(options.IncludeInteriors);
    }

    // Test the constructors and accessors of the Options.
    [Fact]
    internal void Test_S2HausdorffDistanceQueryTest_QueryOptionsAccessorsWorks()
    {
        S2HausdorffDistanceQuery query = new();
        var default_include_interiors = query.Options_.IncludeInteriors;

        query.Options_.IncludeInteriors = !default_include_interiors;
        var modified_include_interiors = query.Options_.IncludeInteriors;

        Assert.True(default_include_interiors);
        Assert.False(modified_include_interiors);
    }

    // Test involving 2 simple polyline shape indexes..
    [Fact]
    internal void Test_S2HausdorffDistanceQueryTest_SimplePolylineQueriesSucceed()
    {
        var a0 = ParsePointsOrDie("0:0, 0:1, 0:1.5");
        var a1 = ParsePointsOrDie("0:2, 0:1.5, -10:1");
        var b0 = ParsePointsOrDie("1:0, 1:1, 3:2");

        // Setup the shape indexes.
        MutableS2ShapeIndex empty_index = new();

        // Shape index a consists or 2 polylines, a0 and a1.
        MutableS2ShapeIndex a = new();
        a.Add(new S2LaxPolylineShape(a0));
        a.Add(new S2LaxPolylineShape(a1));

        // Shape index b consists or 1 polylines: b0.
        MutableS2ShapeIndex b = new();
        b.Add(new S2LaxPolylineShape(b0));

        // Calculate expected distances.
        // Directed a to b HD is achieved at the vertex 2 of a1 and vertex 1 of b0.
        S1ChordAngle expected_a_to_b = new(a1[2], b0[1]);
        // Directed b to a HD is achieved at the vertex 2 of b0 and vertex 0 of a1.
        S1ChordAngle expected_b_to_a = new(b0[2], a1[0]);

        Options options = new();
        S2HausdorffDistanceQuery query = new() { Options_ = options };

        var directed_empty_to_a = query.GetDirectedResult(empty_index, a);
        var directed_a_to_empty = query.GetDirectedResult(a, empty_index);
        var directed_a_to_empty_distance = query.GetDirectedDistance(a, empty_index);

        var directed_a_to_b = query.GetDirectedResult(a, b);
        var directed_b_to_a = query.GetDirectedResult(b, a);
        var directed_a_to_b_distance = query.GetDirectedDistance(a, b);

        var a_to_b = query.GetResult(a, b);
        var b_to_a = query.GetResult(b, a);
        var b_to_a_distance = query.GetDistance(b, a);
        var bb = query.GetResult(b, b);

        Assert.Null(directed_empty_to_a);
        Assert.Null(directed_a_to_empty);
        Assert.True(directed_a_to_empty_distance.IsInfinity());

        Assert.NotNull(directed_a_to_b);
        Assert.NotNull(directed_b_to_a);

        Assert.Equal(directed_a_to_b.Distance.Degrees(),
                     expected_a_to_b.Degrees());
        Assert.Equal(directed_a_to_b_distance.Degrees(),
                     expected_a_to_b.Degrees());
        Assert.Equal(directed_b_to_a.Distance.Degrees(),
                     expected_b_to_a.Degrees());

        Assert.NotNull(a_to_b);
        Assert.NotNull(b_to_a);
        Assert.NotNull(bb);

        Assert.Equal(a_to_b.GetDistance().Degrees(), b_to_a.GetDistance().Degrees());
        Assert.Equal(bb.GetDistance().Degrees(), 0);
        Assert.Equal(
            a_to_b.GetDistance().Degrees(),
            Math.Max(a_to_b.GetDistance().Degrees(), b_to_a.GetDistance().Degrees()));
        Assert.Equal(b_to_a_distance.Degrees(), b_to_a.GetDistance().Degrees());
    }

    // Test involving a polyline shape (dimension == 1) and a point shape (dimension
    // == 0).
    [Fact]
    internal void Test_S2HausdorffDistanceQueryTest_PointVectorShapeQueriesSucceed()
    {
        // Points for the polyline shape.
        var a_points = ParsePointsOrDie("2:0, 0:1, 1:2, 0:3, 0:4");
        // Points for the polint vector shape.
        var b_points = ParsePointsOrDie("-1:2, -0.5:0.5, -0.5:3.5");

        MutableS2ShapeIndex a = new();
        a.Add(new S2LaxPolylineShape(a_points));

        MutableS2ShapeIndex b = new();
        b.Add(new S2PointVectorShape(b_points.ToArray()));

        Options options = new();
        S2HausdorffDistanceQuery query = new() { Options_ = options };

        // Directed Hausdorff distance from a to b is achieved at the vertex 0 of a
        // and vertex 1 of b.
        S1ChordAngle expected_a_to_b = new(a_points[0], b_points[1]);

        // Directed Hausdorff distance from b to a is achieved at the vertex 0 of b
        // and vertex 3 of a.
        S1ChordAngle expected_b_to_a = new(b_points[0], a_points[3]);

        // Unddrected Hausdorff distance between a and b is the maximum of the two
        // directed Hausdorff distances.
        var expected_a_b = S1ChordAngle.Max(expected_a_to_b, expected_b_to_a);

        var directed_a_to_b = query.GetDirectedResult(a, b);
        var directed_b_to_a = query.GetDirectedResult(b, a);
        S1ChordAngle undirected_a_b = query.GetDistance(a, b);

        Assert.NotNull(directed_a_to_b);
        Assert.NotNull(directed_b_to_a);
        Assert.False(undirected_a_b.IsInfinity());

        Assert.Equal(directed_a_to_b!.Distance.Degrees(),
                     expected_a_to_b.Degrees());
        Assert.Equal(directed_a_to_b.TargetPoint, a_points[0]);

        Assert.Equal(directed_b_to_a.Distance.Degrees(),
                     expected_b_to_a.Degrees());
        Assert.Equal(directed_b_to_a.TargetPoint, b_points[0]);

        Assert.Equal(undirected_a_b.Degrees(), expected_a_b.Degrees());
    }

    // Test involving partially overlapping polygons.
    [Fact]
    internal void Test_S2HausdorffDistanceQueryTest_OverlappingPolygons()
    {
        // The first polygon is a triangle. It's first two vertices are inside the
        // quadrangle b (defined below), and the last vertex is outside of b.
        MutableS2ShapeIndex a = new();
        a.Add(MakeLaxPolygonOrDie("1:1, 1:2, 3.5:1.5"));

        // The other polygon is a quadraangle.
        MutableS2ShapeIndex b = new();
        b.Add(MakeLaxPolygonOrDie("0:0, 0:3, 3:3, 3:0"));

        // The first query does not include the interiors.
        Options options = new()
        {
            IncludeInteriors = false
        };
        S2HausdorffDistanceQuery query_1 = new() { Options_ = options };
        var a_to_b_1 = query_1.GetDirectedResult(a, b);

        // The second query has include_interiors set to true.
        options.IncludeInteriors = true;
        S2HausdorffDistanceQuery query_2 = new() { Options_ = options };
        var a_to_b_2 = query_2.GetDirectedResult(a, b);

        // Error tolerance to account for the difference between the northern edge of
        // the quadrangle, which is a geodesic line, and the parallel lat=3 connecting
        // the verices of that edge.
        const double kEpsilon = 3.0e-3;

        // The directed Hausdorff distance from the first query is achieved on the
        // vertex of the triangle that is inside the quadrangle, and is approximately
        // 1 degree away from the nearest edge of the quadrangle.
        S2Point expected_target_point_1 = S2LatLng.FromDegrees(1, 2).ToPoint();

        // The directed Hausdorff distance from the second query is achieved on the
        // last vertex of the triangle that is outside the quadrangle, and is about
        // 0.5 degrees away from the nearest edge of the quadrangle.
        S2Point expected_target_point_2 = S2LatLng.FromDegrees(3.5, 1.5).ToPoint();

        Assert.NotNull(a_to_b_1);
        Assert2.Near(a_to_b_1.Distance.Degrees(), 1, kEpsilon);
        Assert.Equal(a_to_b_1.TargetPoint, expected_target_point_1);

        Assert.NotNull(a_to_b_2);
        Assert2.Near(a_to_b_2.Distance.Degrees(), 0.5, kEpsilon);
        Assert.Equal(a_to_b_2.TargetPoint, expected_target_point_2);
    }
}
