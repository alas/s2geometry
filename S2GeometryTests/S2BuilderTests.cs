namespace S2Geometry;

using System.Text;
using S2BuilderUtil;
using static S2Builder;
using static S2Builder.GraphOptions;
using static S2TextFormat;

public class S2BuilderTests(ITestOutputHelper logger)
{
    private const int iteration_multiplier = 1; // Iteration multiplier for randomized tests
    private const S2ErrorCode INPUT_EDGE_ID_MISMATCH = S2ErrorCode.USER_DEFINED_START;
    private readonly ITestOutputHelper _logger = logger;

    [Fact]
    internal void Test_S2Builder_AddShape()
    {
        S2Builder builder = new(new Options());
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        var input = MakePolygonOrDie("0:0, 0:5, 5:5, 5:0; 1:1, 1:4, 4:4, 4:1");
        builder.AddShape(input.Index.Shape(0));
        Assert.True(builder.Build(out _));
        ExpectPolygonsEqual(input, output);
    }

    [Fact]
    internal void Test_S2Builder_SimpleVertexMerging()
    {
        // When IdentitySnapFunction is used (i.e., no special requirements on
        // vertex locations), check that vertices closer together than the snap
        // radius are merged together.

        S1Angle snap_radius = S1Angle.FromDegrees(0.5);
        S2Builder builder = new(new Options(new IdentitySnapFunction(snap_radius)));
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        var input = MakePolygonOrDie(
            "0:0, 0.2:0.2, 0.1:0.2, 0.1:0.9, 0:1, 0.1:1.1, 0.9:1, 1:1, 1:0.9");
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        var expected = MakePolygonOrDie("0:0, 0:1, 1:0.9");
        ExpectPolygonsApproxEqual(expected, output, snap_radius);
    }

    [Fact]
    internal void Test_S2Builder_SimpleS2CellIdSnapping()
    {
        // When S2CellIdSnapFunction is used, check that all output vertices are the
        // centers of S2CellIds at the specified level level.

        int level = S2CellIdSnapFunction.LevelForMaxSnapRadius(S1Angle.FromDegrees(1));
        S2CellIdSnapFunction snap_function = new(level);
        S2Builder builder = new(new Options(snap_function));
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        var input = MakePolygonOrDie(
            "2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3");
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        Assert.Equal(1, output.NumLoops());
        var loop = output.Loop(0);
        for (int i = 0; i < loop.NumVertices; ++i)
        {
            Assert.Equal(new S2CellId(loop.Vertex(i)).Parent(level).ToPoint(),
                      loop.Vertex(i));
        }
        ExpectPolygonsApproxEqual(input, output, snap_function.SnapRadius);
    }

    [Fact]
    internal void Test_S2Builder_SimpleIntLatLngSnapping()
    {
        S2Builder builder = new(new Options(new IntLatLngSnapFunction(0)));  // E0 coords
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        var input = MakePolygonOrDie(
            "2.01:2.09, 3.24:4.49, 1.78:6.25, 3.51:5.49, 6.11:6.11, " +
            "5.22:3.88, 5.55:2.49, 4.49:2.51");
        var expected = MakePolygonOrDie(
            "2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3");
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        Assert.Equal(1, output.NumLoops());
        ExpectPolygonsEqual(expected, output);
    }

    [Fact]
    internal void Test_S2Builder_VerticesMoveLessThanSnapRadius()
    {
        // Check that chains of closely spaced vertices do not collapse into a
        // single vertex.

        S1Angle snap_radius = S1Angle.FromDegrees(1);
        S2Builder builder = new(new Options(new IdentitySnapFunction(snap_radius)));
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        // The spacing between input vertices is about 2*pi*20/1000 = 0.125 degrees.
        // The output vertices are spaced between 1 and 2 degrees apart; the average
        // spacing is about 1.33 degrees.
        S2Polygon input = new(
            S2Loop.MakeRegularLoop(new S2Point(1, 0, 0), S1Angle.FromDegrees(20), 1000));
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        Assert.Equal(1, output.NumLoops());
        Assert.True(output.Loop(0).NumVertices >= 90);
        Assert.True(output.Loop(0).NumVertices <= 100);
        Assert.True(output.BoundaryNear(input, snap_radius));
    }

    [Fact]
    internal void Test_S2Builder_MinEdgeVertexSeparation()
    {
        // Check that edges are separted from non-incident vertices by at least
        // min_edge_vertex_separation().  This requires adding new vertices (not
        // present in the input) in some cases.

        // The input is a skinny right triangle with two legs of length 10 and 1,
        // and whose diagonal is subdivided into 10 short edges.  Using a snap
        // radius of 0.5, about half of the long leg is snapped onto the diagonal
        // (which causes that part of the polygon to be removed).  But the real
        // problem is that the remaining part of the long leg gets too close to the
        // remaining vertices on the diagonal, i.e. it would violate the minimum
        // edge-vertex separation guarantee.  S2Builder handles this by creating at
        // least one vertex along the original long leg, to keep the snapped edge
        // far enough away from the diagonal.
        var input = MakePolygonOrDie(
            "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 6:.4, 7:.3, 8:.2, 9:.1, 10:0");
        var expected = MakePolygonOrDie(
            "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 4.00021862252687:0");
        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(0.5)));
        S2Builder builder = new(options);
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        ExpectPolygonsApproxEqual(expected, output, S1Angle.FromRadians(S2.DoubleError));
    }

    [Fact]
    internal void Test_S2Builder_MaxEdgeDeviation()
    {
        // Test that when long edges are snapped to nearby sites, appropriate extra
        // vertices are added so that the snapped edge chain stays within
        // options.max_edge_deviation() of the original edge.
        //
        // We do this by constructing an edge AB with two nearly antipodal vertices,
        // then adding another vertex C slightly perturbed from A.  Frequently AB
        // will snap to C producing the edge chain ACB, however edge CB diverges
        // very far from the original edge AB.  This requires S2Builder to add a new
        // vertex M near the middle of AB so that ACMB stays close everywhere to AB.
        S2Builder.Options options = new()
        {
            SplitCrossingEdges = true,
            Idempotent = false
        };
        S2Builder builder=new(options);

        // Even though we are using the default snap radius of zero, the edge snap
        // radius is S2::kIntersectionError because split_crossing_edges() is true.
        Assert.Equal(builder.Options_.EdgeSnapRadius(), S2.kIntersectionErrorS1Angle);
        S1Angle max_deviation = builder.Options_.MaxEdgeDeviation();
        const int kIters = 50 * iteration_multiplier;

        // Test cases are constructed randomly not all tests are effective (i.e., AB
        // might not snap to the perturbed vertex C).  Here we keep track of the
        // number of effective tests.
        int num_effective = 0;
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2LaxPolylineShape output=new();
            builder.StartLayer(new LaxPolylineLayer(output));
            S2Point a = S2Testing.RandomPoint();

            // B is slightly perturbed from -A, and C is slightly perturbed from A.
            // The perturbation amount is not critical except that it should not be
            // too much larger than edge_snap_radius() (otherwise AB will rarely snap
            // to C) and it should also not be so small that A and C are likely to be
            // the same point.
            S2Point b = (-a + 5e-16 * S2Testing.RandomPoint()).Normalize();
            S2Point c = (a + 5e-16 * S2Testing.RandomPoint()).Normalize();
            builder.AddEdge(a, b);
            builder.ForceVertex(c);
            Assert.True(builder.Build(out _));// << error;
            int n = output.NumVertices();
            Assert.Equal(a, output.Vertex(0));
            Assert.Equal(b, output.Vertex(n - 1));
            for (int i = 0; i + 1 < n; ++i)
            {
                Assert.True(S2.IsEdgeBNearEdgeA(
                    a, b, output.Vertex(i), output.Vertex(i + 1), max_deviation));
                _logger.WriteLine(string.Format($"Iteration {0}: ({1}, {2}), {3}",
                    iter,
                    S2TextFormat.ToDebugString(a),
                    S2TextFormat.ToDebugString(S2.Interpolate(a, b, 0.5)),
                    S2TextFormat.ToDebugString(output)));
            }
            if (n > 2) ++num_effective;
        }
        // We require at least 20% of the test cases to be successful.
        Assert.True(num_effective >= 10 * iteration_multiplier);
    }

    [Fact]
    internal void Test_S2Builder_IdempotencySnapsInadequatelySeparatedVertices()
    {
        // This test checks that when vertices are closer together than
        // min_vertex_separation() then they are snapped together even when
        // options.idempotent() is true.
        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(1.0)));
        S2Builder builder = new(options);
        S2Polyline output = new();
        builder.StartLayer(new S2PolylineLayer(output));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:0.9, 0:2"));
        Assert.True(builder.Build(out _));
        string expected = "0:0, 0:2";
        Assert.Equal(expected, output.ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_IdempotencySnapsIdenticalVerticesWithZeroSnapRadius()
    {
        // This test checks that even when the snap radius is zero, identical
        // vertices are snapped together.
        S2Builder builder = new(new Options());
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        builder.AddPolyline(MakePolylineOrDie("0:1, 1:0"));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:1"));
        builder.AddEdge(MakePointOrDie("0:1"), MakePointOrDie("0:1"));
        builder.AddPolyline(MakePolylineOrDie("1:0, 0:0"));
        Assert.True(builder.Build(out _));
        string expected = "0:0, 0:1, 1:0";
        Assert.Equal(expected, output.ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_IdempotencySnapsIdenticalVerticesWithZeroSnapRadiusEdgeSplitting()
    {
        // This test checks that identical vertices are snapped together even when
        // the snap radius is zero and options.split_crossing_edges() is true.
        Options options = new()
        {
            SplitCrossingEdges = true
        };
        S2Builder builder = new(options);
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        builder.AddPolyline(MakePolylineOrDie("0:1, 1:0"));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:1"));
        builder.AddEdge(MakePointOrDie("0:1"), MakePointOrDie("0:1"));
        builder.AddPolyline(MakePolylineOrDie("1:0, 0:0"));
        Assert.True(builder.Build(out _));
        string expected = "0:0, 0:1, 1:0";
        Assert.Equal(expected, output.ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_IdempotencySnapsUnsnappedVertices()
    {
        // When idempotency is requested, no snapping is done unless S2Builder finds
        // at least one vertex or edge that could not be the output of a previous
        // snapping operation.  This test checks that S2Builder detects vertices
        // that are not at a valid location returned by the given snap function.

        // In this example we snap two vertices to integer lat/lng coordinates.  The
        // two vertices are far enough apart (more than min_vertex_separation) so
        // that they might be the result of a previous snapping operation, but one
        // of the two vertices does not have integer lat/lng coordinates.  We use
        // internal knowledge of how snap sites are chosen (namely, that candidates
        // are considered in S2CellId order) to construct two different cases, one
        // where the snapped vertex is processed first and one where the unsnapped
        // vertex is processed first.  This exercises two different code paths.
        IntLatLngSnapFunction snap_function = new(0);
        Assert.True(snap_function.SnapRadius>= S1Angle.FromDegrees(0.7));
        Assert.True(snap_function.MinVertexSeparation() <= S1Angle.FromDegrees(0.35));
        S2Builder builder = new(new Options(snap_function));

        // In this example, the snapped vertex (0, 0) is processed first and is
        // selected as a Voronoi site (i.e., output vertex).  The second vertex is
        // closer than min_(), therefore it is snapped to the first vertex
        // and the polyline becomes degenerate.
        S2Point a = S2LatLng.FromDegrees(0, 0).ToPoint();
        S2Point b = S2LatLng.FromDegrees(0.01, 0.6).ToPoint();
        Assert.True(new S2CellId(a) < new S2CellId(b));
        S2Polyline input1 = new(new[] { a, b }), output1 = new();
        builder.StartLayer(new S2PolylineLayer(output1));
        builder.AddPolyline(input1);
        Assert.True(builder.Build(out _));
        Assert.Equal("0:0, 0:1", output1.ToDebugString());

        // In this example the unsnapped vertex is processed first and is snapped to
        // (0, 0).  The second vertex is further than snap_radius() away, so it is
        // also snapped (which does nothing) and is left at (0, 1).
        S2Point c = S2LatLng.FromDegrees(0.01, 0.4).ToPoint();
        S2Point d = S2LatLng.FromDegrees(0, 1).ToPoint();
        Assert.True(new S2CellId(c) < new S2CellId(d));
        S2Polyline input2 = new(new[] { c, d }), output2 = new();
        builder.StartLayer(new S2PolylineLayer(output2));
        builder.AddPolyline(input2);
        Assert.True(builder.Build(out _));
        Assert.Equal("0:0, 0:1", output2.ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_IdempotencySnapsEdgesWithTinySnapRadius()
    {
        // When idempotency is requested, no snapping is done unless S2Builder finds
        // at least one vertex or edge that could not be the output of a previous
        // snapping operation.  This test checks that S2Builder detects edges that
        // are too close to vertices even when the snap radius is very small
        // (e.g., S2EdgeCrossings.kIntersectionError).
        //
        // Previously S2Builder used a conservative approximation to decide whether
        // edges were too close to vertices; unfortunately this meant that when the
        // snap radius was very small then no snapping would be done at all, because
        // even an edge/vertex distance of zero was considered far enough apart.
        //
        // This tests that the current code (which uses exact predicates) handles
        // this situation correctly (i.e., that an edge separated from a
        // non-incident vertex by a distance of zero cannot be the output of a
        // previous snapping operation).
        Options options = new(new IdentitySnapFunction(S2.kIntersectionErrorS1Angle));
        S2PolylineVectorLayer.Options layer_options = new()
        {
            DuplicateEdges_ = DuplicateEdges.MERGE
        };
        S2Builder builder = new(options);
        List<S2Polyline> output = [];
        builder.StartLayer(
            new S2PolylineVectorLayer(output, layer_options));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:10"));
        builder.AddPolyline(MakePolylineOrDie("0:5, 0:7"));
        Assert.True(builder.Build(out _));
        Assert.Single(output);
        Assert.Equal("0:0, 0:5, 0:7, 0:10", output[0].ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_IdempotencyDoesNotSnapAdequatelySeparatedEdges()
    {
        // When idempotency is requested, no snapping is done unless S2Builder finds
        // at least one vertex or edge that could not be the output of a previous
        // snapping operation.  This test checks that when an edge is further away
        // than min_edge_vertex_separation() then no snapping is done.
        Options options = new(new IntLatLngSnapFunction(0))
        {
            Idempotent = true  // Test fails if this is "false".
        };
        S2Builder builder = new(options);
        S2Polygon output1 = new(), output2 = new();
        builder.StartLayer(new S2PolygonLayer(output1));
        builder.AddPolygon(MakePolygonOrDie("1.49:0, 0:2, 0.49:3"));
        Assert.True(builder.Build(out _));
        string expected = "1:0, 0:2, 0:3";
        Assert.Equal(expected, output1.ToDebugString());
        builder.StartLayer(new S2PolygonLayer(output2));
        builder.AddPolygon(output1);
        Assert.True(builder.Build(out _));
        Assert.Equal(expected, output2.ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_NearbyVerticesSnappedWithZeroSnapRadiusEdgeSplitting()
    {
        // Verify that even when the split_crossing_edges() option is used with a snap
        // radius of zero, edges are snapped to nearby vertices (those within a
        // distance of S2::kIntersectionError).  This is necessary for correctness
        // whenever new intersection vertices are created (since such vertices are
        // only within S2::kIntersectionError of the true intersection point), and it
        // is necessary for consistency even when they are not.
        S2Builder.Options options = new()
        {
            SplitCrossingEdges = true
        };
        S2Builder builder=new(options);
        S2PolylineVectorLayer.Options layer_options = new()
        {
            PolylineType_ = Graph.PolylineType.WALK
        };
        List<S2Polyline> output=[];
        builder.StartLayer(
            new S2PolylineVectorLayer(output, layer_options));
        builder.AddPolyline(MakePolylineOrDie("0:180, 0:3"));
        // The second value below produces an S2Point that is distinct from 0:180
        // and yet is so close that 0:180 is the nearest representable value when
        // it is converted back to an S2LatLng.
        builder.AddPolyline(MakePolylineOrDie("90:180, 0:179.9999999999999"));
        builder.AddPolyline(MakePolylineOrDie("10:10, 1e-15:10"));
        Assert.True(builder.Build(out S2Error error));
        Assert.Equal(3, output.Count);
        // The first two points below are not duplicates (see above).
        Assert.Equal("0:180, 0:180, 1e-15:10, 0:3", S2TextFormat.ToDebugString(output[0]));
        Assert.Equal("90:180, 0:180", S2TextFormat.ToDebugString(output[1]));
        Assert.Equal("10:10, 1e-15:10", S2TextFormat.ToDebugString(output[2]));
    }

    [Fact]
    internal void Test_S2Builder_NearbyIntersectionSnappedWithZeroSnapRadius()
    {
        // A simpler version of the test above that uses AddIntersection() rather
        // than split_crossing_edges().
        S2Builder.Options options = new()
        {
            IntersectionTolerance = S2.kIntersectionErrorS1Angle
        };
        S2Builder builder=new(options);
        S2Polyline output=new();
        builder.StartLayer(new S2PolylineLayer(output));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:10"));
        builder.AddIntersection(MakePointOrDie("1e-16:5"));
        Assert.True(builder.Build(out _));
        Assert.Equal("0:0, 1e-16:5, 0:10", S2TextFormat.ToDebugString(output));
    }

    [Fact]
    internal void Test_S2Builder_TopologyPreservedWithZeroSnapRadiusEdgeSplitting()
    {
        // Verify that even when the split_crossing_edges() option is used with a snap
        // radius of zero, snapped edges do not cross vertices (i.e., the input
        // topology is preserved up the creation of degeneracies).  In this situation
        // the snap radius for vertices is zero (i.e. vertices are never merged) but
        // the snap radius for edges is S2::kIntersectionError (i.e., edges are
        // snapped to vertices up to this distance away).
        //
        // For this test, we create an edge AB that is 47 degrees long.  We then add
        // two vertices X,Y that are 1 degree away from the endpoints of A,B and
        // offset from AB by a distance of 0.99 * S2::kIntersectionError.  This will
        // cause AB to be snapped to AXYB, where the segment XY is 45 degrees long.
        // Because all edges are geodesics, the distance from the midpoint M of XY to
        // the original edge AB is approximately 1.07 * S2::kIntersectionError (see
        // the comments in S2Builder::Options::max_edge_deviation()).  Let vertex C be
        // offset from the midpoint of AB by 1.03 * S2::kIntersectionError.  Then C is
        // too far away from AB to be snapped to it, but it is closer to AB than M.
        // Therefore if nothing is done about it, the input topology would change.
        //
        // This test checks that the algorithm detects this situation and adds another
        // vertex Z near the projection of C onto AB.  This causes AB to be snapped to
        // AXZYB so that the snapped edge passes on the correct side of C.
        //
        // Vertex C is represented in the form of an edge CD perpendicular to AB so
        // that we don't need to use more than one output layer.  We also need to turn
        // off the idempotent() option to ensure that AB is snapped to X and Y, since
        // otherwise the algorithm correctly determines that there is no need for
        // snapping (since there are no crossing edges in this test case).
        S2Builder.Options options = new()
        {
            SplitCrossingEdges = true,
            Idempotent = false
        };
        S2Builder builder=new(options);
        S2PolylineVectorLayer.Options layer_options = new()
        {
            PolylineType_ = Graph.PolylineType.WALK
        };
        List<S2Polyline> output=[];
        builder.StartLayer(
            new S2PolylineVectorLayer(output, layer_options));
        var kEdgeSnapRadDegrees = S2.kIntersectionErrorS1Angle.GetDegrees();
        var a = S2LatLng.FromDegrees(0, -1).ToPoint();
        var b = S2LatLng.FromDegrees(0, 46).ToPoint();
        var x = S2LatLng.FromDegrees(0.99 * kEdgeSnapRadDegrees, 0).ToPoint();
        var y = S2LatLng.FromDegrees(0.99 * kEdgeSnapRadDegrees, 45).ToPoint();
        var c = S2LatLng.FromDegrees(1.03 * kEdgeSnapRadDegrees, 22.5).ToPoint();
        var d = S2LatLng.FromDegrees(10, 22.5).ToPoint();
        builder.AddEdge(a, b);
        builder.ForceVertex(x);
        builder.ForceVertex(y);
        builder.AddEdge(c, d);
        Assert.True(builder.Build(out S2Error error));
        Assert.Equal(2, output.Count);
        Assert.Equal(  // The snapped edge AXZYB described above.
            "0:-1, 5.03799861543821e-14:0, 0:22.5, 5.03799861543821e-14:45, 0:46",
            S2TextFormat.ToDebugString(output[0]));
        Assert.Equal(  // The input edge CD.
            "5.24155411505188e-14:22.5, 10:22.5",
            S2TextFormat.ToDebugString(output[1]));
        Assert.True(S2.CrossingSign(output[0].Vertex(1), output[0].Vertex(2),
                                   output[1].Vertex(0), output[1].Vertex(1)) < 0);
    }

    [Fact]
    internal void Test_S2Builder_TopologyPreservedWithForcedVertices()
    {
        // Verify that the input topology is preserved around vertices added using
        // ForceVertex (even though there are no minimum separation guarantees for
        // such vertices).
        //
        // This test is the same as the one above except for the following:
        //  - split_crossing_edges() is false
        //  - we use a snap raidus of S2::kIntersectionError rather than zero
        //  - vertex C is added using ForceVertex().
        S2Builder.Options options = new(new IdentitySnapFunction(S2.kIntersectionErrorS1Angle))
        {
            Idempotent = false
        };
        S2Builder builder = new(options);
        S2PolylineVectorLayer.Options layer_options = new()
        {
            PolylineType_ = Graph.PolylineType.WALK
        };
        List<S2Polyline> output=[];
        builder.StartLayer(
            new S2PolylineVectorLayer(output, layer_options));
        var kEdgeSnapRadDegrees = S2.kIntersectionErrorS1Angle.GetDegrees();
        var a = S2LatLng.FromDegrees(0, -1).ToPoint();
        var b = S2LatLng.FromDegrees(0, 46).ToPoint();
        var x = S2LatLng.FromDegrees(0.99 * kEdgeSnapRadDegrees, 0).ToPoint();
        var y = S2LatLng.FromDegrees(0.99 * kEdgeSnapRadDegrees, 45).ToPoint();
        var c = S2LatLng.FromDegrees(1.03 * kEdgeSnapRadDegrees, 22.5).ToPoint();
        var d = S2LatLng.FromDegrees(10, 22.5).ToPoint();
        builder.AddEdge(a, b);
        builder.ForceVertex(x);
        builder.ForceVertex(y);
        builder.ForceVertex(c);
        builder.AddEdge(c, d);
        Assert.True(builder.Build(out S2Error error));
        Assert.Equal(2, output.Count);
        Assert.Equal(  // The snapped edge AXZYB described above.
            "0:-1, 5.03799861543821e-14:0, 0:22.5, 5.03799861543821e-14:45, 0:46",
            S2TextFormat.ToDebugString(output[0]));
        Assert.Equal(  // The input edge CD.
            "5.24155411505188e-14:22.5, 10:22.5",
            S2TextFormat.ToDebugString(output[1]));
        Assert.True(S2.CrossingSign(output[0].Vertex(1), output[0].Vertex(2),
                                   output[1].Vertex(0), output[1].Vertex(1)) < 0);
    }

    [Fact]
    internal void Test_S2Builder_kMaxSnapRadiusCanSnapAtLevel0()
    {
        // Verify that kMaxSnapRadius will allow snapping at S2CellId level 0.
        Assert.True(S2CellIdSnapFunction.MinSnapRadiusForLevel(0) <=
                  SnapFunction.kMaxSnapRadius);
    }

    [Fact]
    internal void Test_S2Builder_S2CellIdSnappingAtAllLevels()
    {
        S2Polygon input = MakePolygonOrDie(
            "0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            S2CellIdSnapFunction snap_function = new(level);
            S2Builder builder = new(new Options(snap_function));
            S2Polygon output = new();
            builder.StartLayer(new S2PolygonLayer(output));
            builder.AddPolygon(input);
            Assert.True(builder.Build(out _));
            Assert.True(output.IsValid());
            // The ApproxContains calls below are not guaranteed to succeed in general
            // because ApproxContains works by snapping both polygons together using
            // the given tolerance and then checking for containment.  Since
            // ApproxContains snaps to an arbitrary subset of the input vertices
            // rather than to S2CellId centers at the current level, this means that
            // corresponding vertices in "input" and "output" can snap to different
            // sites, which causes the containment test to fail.  Nevertheless, by
            // using a larger tolerance of 2 * snap_radius, all calls in this test
            // succeed (and would be likely to succeed in other similar tests).
            // (To guarantee correctness we would need to use S2CellIdSnapFunction
            // within the ApproxContains implementation.)
            S1Angle tolerance = new[]{ 2 * snap_function.SnapRadius,
                        SnapFunction.kMaxSnapRadius }.Min();
            Assert.True(output.ApproxContains(input, tolerance));
            Assert.True(input.ApproxContains(output, tolerance));
        }
    }

    [Fact]
    internal void Test_S2Builder_SnappingDoesNotRotateVertices()
    {
        // This is already tested extensively elsewhere.
        var input = MakePolygonOrDie(
            "49.9305505:-124.8345463, 49.9307448:-124.8299657, " +
            "49.9332101:-124.8301996, 49.9331224:-124.8341368; " +
            "49.9311087:-124.8327042, 49.9318176:-124.8312621, " +
            "49.9318866:-124.8334451");
        Options options = new(new S2CellIdSnapFunction());
        S2Builder builder = new(options);
        S2Polygon output1 = new(), output2 = new();
        builder.StartLayer(new S2PolygonLayer(output1));
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        // This checks that the vertices are in the same cyclic order, and that
        // vertices have not moved by more than "snap_radius".
        ExpectPolygonsApproxEqual(input, output1,
                            options.SnapFunction.SnapRadius);

        // Check that snapping twice doesn't rotate the vertices.  This also
        // verifies that S2Builder can be used again after Build() is called.
        builder.StartLayer(new S2PolygonLayer(output2));
        builder.AddPolygon(output1);
        Assert.True(builder.Build(out _));
        ExpectPolygonsEqual(output1, output2);
    }

    [Fact]
    internal void Test_S2Builder_SelfIntersectingPolyline()
    {
        // Check that when two edges of a polyline cross, the intersection point is
        // added to both edges.

        Options options = new();
        IntLatLngSnapFunction snap_function = new(1);  // Snap to E1 coordinates
        options.SnapFunction = snap_function;
        options.SplitCrossingEdges = true;
        S2Builder builder = new(options);
        S2Polyline output = new();
        builder.StartLayer(new S2PolylineLayer(output));
        var input = MakePolylineOrDie("3:1, 1:3, 1:1, 3:3");
        var expected = MakePolylineOrDie("3:1, 2:2, 1:3, 1:1, 2:2, 3:3");
        builder.AddPolyline(input);
        Assert.True(builder.Build(out _));
        ExpectPolylinesEqual(expected, output);
    }

    [Fact]
    internal void Test_S2Builder_SelfIntersectingPolygon()
    {
        // Check that when two edge of a polygon cross, the intersection point is
        // added to both edges, and that the resulting (undirected) edges can be
        // assembled into a valid polygon.

        IntLatLngSnapFunction snap_function = new(1);  // Snap to E1 coordinates
        Options options = new()
        {
            SnapFunction = snap_function,
            SplitCrossingEdges = true
        };
        S2Builder builder = new(options);
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(
            output, new S2PolygonLayer.Options(EdgeType.UNDIRECTED)));
        var input = MakePolylineOrDie("3:1, 1:3, 1:1, 3:3, 3:1");
        var expected = MakePolygonOrDie("1:1, 1:3, 2:2; 3:3, 3:1, 2:2");
        builder.AddPolyline(input);
        Assert.True(builder.Build(out _));
        ExpectPolygonsEqual(expected, output);
    }

    [Fact]
    internal void Test_S2Builder_TieBreakingIsConsistent()
    {
        // Check that when an edge passes between two equally distant vertices, that
        // the choice of which one to snap to does not depend on the edge direction.

        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(2)))
        {
            Idempotent = false
        };
        S2Builder builder = new(options);
        builder.ForceVertex(S2LatLng.FromDegrees(1, 0).ToPoint());
        builder.ForceVertex(S2LatLng.FromDegrees(-1, 0).ToPoint());
        S2Polyline output1 = new(), output2 = new();
        builder.StartLayer(new S2PolylineLayer(output1));
        builder.AddPolyline(MakePolylineOrDie("0:-5, 0:5"));
        builder.StartLayer(new S2PolylineLayer(output2));
        builder.AddPolyline(MakePolylineOrDie("0:5, 0:-5"));
        Assert.True(builder.Build(out _));
        Assert.Equal(3, output1.NumVertices());
        Assert.Equal(3, output2.NumVertices());
        for (int i = 0; i < 3; ++i)
        {
            Assert.Equal(output1.Vertex(i), output2.Vertex(2 - i));
        }
    }

    [Fact]
    internal void Test_S2Builder_GraphPersistence()
    {
        // Ensure that the Graph objects passed to S2Builder.Layer.Build() methods
        // remain valid until all layers have been built.
        List<Graph> graphs = [];
        List<GraphClone> clones = [];
        S2Builder builder = new(new Options());
        for (int i = 0; i < 20; ++i)
        {
            builder.StartLayer(new GraphPersistenceLayer(
                new GraphOptions(), graphs, clones));
            for (int n = S2Testing.Random.Uniform(10); n > 0; --n)
            {
                builder.AddEdge(S2Testing.RandomPoint(), S2Testing.RandomPoint());
            }
        }
        Assert.True(builder.Build(out _));
    }

    [Fact]
    internal void Test_S2Builder_SimplifyOneEdge()
    {
        // Simplify a perturbed edge chain into a single edge.

        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(1)))
        {
            SimplifyEdgeChains = true
        };
        TestPolylineLayersBothEdgeTypes(["0:0, 1:0.5, 2:-0.5, 3:0.5, 4:-0.5, 5:0"],
                                        ["0:0, 5:0"],
                                        new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyNearlyAntipodal()
    {
        // Verify that nothing goes wrong when attempting to simplify a nearly
        // antipodal edge.

        S2Builder.Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(1)))
        {
            SimplifyEdgeChains = true
        };
        TestPolylineLayersBothEdgeTypes(["0:180, 0:1e-09, 32:32"],
                              ["0:180, 0:1e-09, 32:32"],
                              new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyTwoLayers()
    {
        // Construct two layers, each containing a polyline that could be simplified
        // to a single edge on its own.  However the two polylines actually cross,
        // so make sure that the output still contains the intersection vertex.

        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(0.5)))
        {
            SplitCrossingEdges = true,
            SimplifyEdgeChains = true
        };
        TestPolylineLayersBothEdgeTypes(
            ["-2:-1, -1:0, 1:0, 2:1", "1:-2, 0:-1, 0:1, -1:2"],
            ["-2:-1, 0:0, 2:1", "1:-2, 0:0, -1:2"],
            new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyOneLoop()
    {
        // Simplify a regular loop with 1000 vertices and a radius of 20 degrees.
        // Turning on edge chain simplification yields a dramatically smaller number
        // of vertices than snapping alone (10 vertices vs 95 vertices using a snap
        // radius of 1 degree).  This is because snapping alone yields vertices that
        // stay within 1 degree of the input *vertices*, while simplifying edge
        // chains yields edges that stay within 1 degree of the input *edges*.

        for (int i = 0; i < 2; ++i)
        {
            EdgeType edge_type = (EdgeType)i;
            S1Angle snap_radius = S1Angle.FromDegrees(1);
            Options options = new(new IdentitySnapFunction(snap_radius))
            {
                SimplifyEdgeChains = true
            };
            S2Builder builder = new(options);
            S2Polygon output = new();
            builder.StartLayer(new S2PolygonLayer(
                output, new S2PolygonLayer.Options(edge_type)));
            // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
            S2Polygon input = new(
                S2Loop.MakeRegularLoop(new S2Point(1, 0, 0), S1Angle.FromDegrees(20), 1000));
            builder.AddPolygon(input);
            Assert.True(builder.Build(out _));
            Assert.Equal(1, output.NumLoops());
            Assert.True(output.Loop(0).NumVertices >= 10);
            Assert.True(output.Loop(0).NumVertices <= 12);
            Assert.True(output.BoundaryNear(input, snap_radius));
        }
    }

    [Fact]
    internal void Test_S2Builder_SimplifyOppositeDirections()
    {
        // We build two layers with two polylines that follow the same circular arc
        // in opposite directions, and verify that they are snapped identically.
        // (The snap radius is adjusted so that the arc is simplified into a long
        // edge and a short edge, and therefore we would get a different result if
        // the two layers followed the edge chain in different directions.)

        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(0.5)))
        {
            SimplifyEdgeChains = true
        };
        TestPolylineLayersBothEdgeTypes(
            ["-4:0.83, -3:0.46, -2:0.2, -1:0.05, 0:0, 1:0.5, 2:0.2, 3:0.46, 4:0.83",
        "4:.83, 3:.46, 2:.2, 1:.05, 0:0, -1:.5, -2:.2, -3:.46, -4:.83"],
            ["-4:0.83, -2:0.2, 4:0.83", "4:0.83, -2:0.2, -4:0.83"],
            new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyKeepsEdgeVertexSeparation()
    {
        // We build two layers each containing a polyline, such that the polyline in
        // the first layer could be simplified to a straight line except that then
        // it would approach the second polyline too closely.

        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(1.0)))
        {
            SimplifyEdgeChains = true
        };
        TestPolylineLayersBothEdgeTypes(
            ["0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"],
            ["0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"],
            new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyBacktrackingEdgeChain()
    {
        // Test simplifying an edge chain that backtracks on itself.  (This should
        // prevent simplification, since edge chains are approximated parametrically
        // rather than geometrically.)

        Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(0.5)))
        {
            SimplifyEdgeChains = true
        };
        TestPolylineLayersBothEdgeTypes(
            ["0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 4:0, 3:0, "+
        "2:0, 3:0, 4:0, 5:0, 6:0, 7:0"],
            ["0:0, 2:0, 5:0, 2:0, 5:0, 7:0"],
            new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyAvoidsBacktrackingVertices()
    {
        // As noted above, edge chains must proceed monotonically away from the
        // source vertex in order to be simplified.  However in rare cases, adding a
        // new vertex to the chain (e.g. extending ABC with a vertex D) can require
        // us to avoid a nearby vertex (E) that is closer than the previous endpoint
        // of the chain (C).  A previous version of the algorithm did not handle
        // this case correctly.
        //
        // The first polyline (ABC) below cannot be simplified to AC that edge would
        // pass too close to the first vertex of the second polyline DE.  Vertex D
        // is not avoided while processing the edge AB because AD > AB; instead it
        // should be processed when edge BC is added to the simplified chain.
        S2Builder.Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(1.0)))
        {
            SimplifyEdgeChains = true
        };
        Assert.True(S2.GetDistance(MakePointOrDie("0:1.05"),
                                 MakePointOrDie("0:0"), MakePointOrDie("1:2")) <
                 options.SnapFunction.MinEdgeVertexSeparation());
        TestPolylineLayersBothEdgeTypes(
          ["0:0, 1:0.1, 1:2", "0:1.05, -10:1.05"],
  ["0:0, 1:0.1, 1:2", "0:1.05, -10:1.05"],
  new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyLimitsEdgeDeviation()
    {
        // Make sure that simplification does not create long edges such that the
        // midpoint of the edge might be further than max_edge_deviation() from an
        // input edge.  In the example below, vertices are snapped to integer
        // lat/lng coordinates, and the snap radius is approximately 0.707 degrees.
        // Snapping moves the input vertices perpendicular to the input edge by just
        // slightly less than the snap radius (0.693 degrees).  Now the midpoint of
        // the snapped edge is about 0.98 degrees from the input edge, which causes
        // an extra site to be added at the midpoint of the original edge.
        //
        // When simplify_edge_chains() is enabled, then usually an extra site like
        // this would be simplified away (because the simplified edge would still be
        // within snap_radius() of all the input vertices) except that there is an
        // explicit check in S2Builder that prevents this.  (If the check is removed
        // then this test fails.)

        Options options = new(new IntLatLngSnapFunction(0))
        {
            SimplifyEdgeChains = true
        };  // E0 coordinates
        TestPolylineLayersBothEdgeTypes(
            ["-30.49:-29.51, 29.51:30.49"], ["-30:-30, -1:1, 30:30"],
            new S2PolylineLayer.Options(), options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyPreservesTopology()
    {
        // Crate several nested concentric loops, and verify that the loops are
        // still nested after simplification.

        int kNumLoops = 20;
        int kNumVerticesPerLoop = 1000;
        S1Angle kBaseRadius = S1Angle.FromDegrees(5);
        S1Angle kSnapRadius = S1Angle.FromDegrees(0.1);
        Options options = new(new IdentitySnapFunction(kSnapRadius))
        {
            SimplifyEdgeChains = true
        };
        S2Builder builder = new(options);
        List<S2Polygon> input = [], output = [];
        for (int j = 0; j < kNumLoops; ++j)
        {
            // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
            S1Angle radius = kBaseRadius + 0.7 * j * j / kNumLoops * kSnapRadius;
            input.Add(new S2Polygon(S2Loop.MakeRegularLoop(
                new S2Point(1, 0, 0), radius, kNumVerticesPerLoop)));
            output.Add(new S2Polygon());
            builder.StartLayer(new S2PolygonLayer(output.Last()));
            builder.AddPolygon(input.Last());
        }
        Assert.True(builder.Build(out _));
        for (int j = 0; j < kNumLoops; ++j)
        {
            Assert.True(output[j].BoundaryNear(input[j], kSnapRadius));
            if (j > 0) Assert.True(output[j].Contains(output[j - 1]));
        }
    }

    [Fact]
    internal void Test_S2Builder_SimplifyRemovesSiblingPairs()
    {
        Options options = new(new IntLatLngSnapFunction(0));  // E0 coords
        S2PolylineVectorLayer.Options layer_options = new()
        {
            SiblingPairs = SiblingPairs.DISCARD
        };

        // Check that there is no sibling pair without simplification.
        TestPolylineVector(
  ["0:0, 0:10", "0:10, 0.6:5, 0:0"],
  ["0:0, 0:10, 1:5, 0:0"], layer_options, options);

        // Now check that (1) simplification produces a sibling pair,
        // and (2) the sibling pair is removed (since we requested it).
        options.SimplifyEdgeChains = true;
        TestPolylineVector(
["0:0, 0:10", "0:10, 0.6:5, 0:0"],
[], layer_options, options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyMergesDuplicateEdges()
    {
        Options options = new(new IntLatLngSnapFunction(0));  // E0 coords
        S2PolylineVectorLayer.Options layer_options = new()
        {
            DuplicateEdges_ = DuplicateEdges.MERGE
        };

        // Check that there are no duplicate edges without simplification.
        TestPolylineVector(
  ["0:0, 0:10", "0:0, 0.6:5, 0:10"],
  ["0:0, 0:10", "0:0, 1:5, 0:10"], layer_options, options);

        // Now check that (1) simplification produces a duplicate edge pair,
        // and (2) the duplicate pair is merged (since we requested it).
        options.SimplifyEdgeChains = true;
        TestPolylineVector(
  ["0:0, 0:10", "0:0, 0.6:5, 0:10"],
  ["0:0, 0:10"], layer_options, options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyKeepsForcedVertices()
    {
        Options options = new(new IdentitySnapFunction(S1Angle.FromRadians(S2.DoubleError)))
        {
            SimplifyEdgeChains = true
        };
        S2Builder builder = new(options);
        S2Polyline output = new();
        builder.StartLayer(new S2PolylineLayer(output));
        builder.AddPolyline(MakePolylineOrDie("0:0, 0:1, 0:2, 0:3"));
        builder.ForceVertex(MakePointOrDie("0:1"));
        Assert.True(builder.Build(out _));
        Assert.Equal("0:0, 0:1, 0:3", output.ToDebugString());
    }

    [Fact]
    internal void Test_S2Builder_InputEdgeIdAssignment()
    {
        // Check that input edge ids are assigned in order.
        TestInputEdgeIds(["0:0, 0:1, 0:2"], [
  ("0:0, 0:1", new[]{ 0 }), ("0:1, 0:2", new[]{ 1 })],
  new GraphOptions(), new Options());
    }

    [Fact]
    internal void Test_S2Builder_UndirectedSiblingsDontHaveInputEdgeIds()
    {
        // Check that the siblings of undirected edges do not have InputEdgeIds.
        GraphOptions graph_options = new()
        {
            EdgeType_ = EdgeType.UNDIRECTED
        };
        TestInputEdgeIds(["0:0, 0:1, 0:2"], [
  ("0:0, 0:1", new[]{0}), ("0:1, 0:2", new[]{1}),
  ("0:1, 0:0", Array.Empty<int>()), ("0:2, 0:1", Array.Empty<int>())],
               graph_options, new Options());
    }

    [Fact]
    internal void Test_S2Builder_CreatedSiblingsDontHaveInputEdgeIds()
    {
        // Check that edges created by SiblingPairs.CREATE do not have
        // InputEdgeIds.
        GraphOptions graph_options = new()
        {
            SiblingPairs_ = SiblingPairs.CREATE
        };
        TestInputEdgeIds(["0:0, 0:1, 0:2"], [
  ("0:0, 0:1", new[]{0}), ("0:1, 0:2", new[]{1})],
  new GraphOptions(), new Options());
    }

    [Fact]
    internal void Test_S2Builder_EdgeMergingDirected()
    {
        // Tests that input edge ids are merged when directed edges are merged.
        GraphOptions graph_options = new()
        {
            DuplicateEdges_ = DuplicateEdges.MERGE
        };
        TestInputEdgeIds(["0:0, 0:1", "0:0, 0:1"], [("0:0, 0:1", new[] { 0, 1 })],
               graph_options, new Options());
    }

    [Fact]
    internal void Test_S2Builder_EdgeMergingUndirected()
    {
        // Tests that input edge ids are merged when undirected edges are merged.
        GraphOptions graph_options = new()
        {
            DuplicateEdges_ = DuplicateEdges.MERGE,
            SiblingPairs_ = SiblingPairs.KEEP
        };
        TestInputEdgeIds(["0:0, 0:1, 0:2", "0:0, 0:1", "0:2, 0:1"], [
  ("0:0, 0:1", new[]{0, 2}), ("0:1, 0:2", new[]{1}), ("0:2, 0:1", new[]{3})
], graph_options, new Options());
    }

    [Fact]
    internal void Test_S2Builder_SimplifyDegenerateEdgeMergingEasy()
    {
        // Check that when an input edge is snapped to a chain that includes
        // degenerate edges, and the edge chain is simplified, that the InputEdgeIds
        // attached to those degenerate edges are transferred to the simplified
        // edge.  For example (using integers for vertices), an edge chain 1.2,
        // 2.2, 2.3 that is simplified to 1.3 should get the InputEdgeIds
        // associated with all three original edges.  (This ensures that the labels
        // attached to those edges are also transferred.)
        //
        // This also tests that degenerate edges at the start and end of the
        // simplified chain are *not* merged.  (It's up to the output layer to
        // decide what to do with these edges.  The only reason we merge degenerate
        // edges in the interior of the interior of the simplified edge is because
        // those edges are being removed from the graph.)
        GraphOptions graph_options = new()
        {
            DegenerateEdges_ = DegenerateEdges.KEEP
        };
        Options options = new(new IntLatLngSnapFunction(0))
        {
            SimplifyEdgeChains = true
        };
        TestInputEdgeIds(["0:0, 0:0.1, 0:1.1, 0:1, 0:0.9, 0:2, 0:2.1"], [
  ("0:0, 0:0", new[]{0}), ("0:0, 0:2", new[]{1, 2, 3, 4}), ("0:2, 0:2", new[]{5})
], graph_options, options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyDegenerateEdgeMergingHard()
    {
        // This is a harder version of the test above.  Now there are several edge
        // chains that overlap each other in both directions, and several degenerate
        // edges at that middle vertex.  This tests that if exactly one edge chain
        // contains a degenerate edge in input edge order (e.g., the input order was
        // AB, BB, BC), then the degenerate edge is assigned to that edge chain.
        // Otherwise the edge is assigned to an arbitrary chain.
        GraphOptions graph_options = new();  // Default options keep everything.
        Options options = new(new IntLatLngSnapFunction(0))
        {
            SimplifyEdgeChains = true
        };
        var input = new List<string>{
"0:1, 0:1.1", "0:0, 0:1, 0:2",  // Degenerate edge defined before chain
"0:0, 0:0.9, 0:1, 0:1.1, 0:2",  // Degenerate edge defined in chain
"0:2, 0:1, 0:0.9, 0:0",         // Defined in chain, chain reversed
"0:2, 0:1, 0:0", "0:1.1, 0:1", "0:1, 0:1.1",  // Defined after chain
};
        var expected = new EdgeInputEdgeIds{
("0:0, 0:2", new[]{0, 1, 2}), ("0:0, 0:2", new[]{3, 4, 5, 6}),
("0:2, 0:0", new[]{7, 8, 9}), ("0:2, 0:0", new[]{10, 11, 12, 13})
};
        TestInputEdgeIds(input, expected, graph_options, options);

        // Now try the same test with undirected edges.  This results in four more
        // simplified edges that are not labelled with any input edge ids.
        expected.AddRange(new EdgeInputEdgeIds{
  ("0:0, 0:2", Array.Empty<int>()), ("0:0, 0:2", Array.Empty<int>()),
  ("0:2, 0:0", Array.Empty<int>()), ("0:2, 0:0", Array.Empty<int>())});
        graph_options.EdgeType_ = EdgeType.UNDIRECTED;
        TestInputEdgeIds(input, expected, graph_options, options);
    }

    [Fact]
    internal void Test_S2Builder_SimplifyDegenerateEdgeMergingMultipleLayers()
    {
        // Check that degenerate edges are assigned to an edge in the correct layer
        // when multiple edge chains in different layers are simplified in the same
        // way (i.e., yielding a set of identical or reversed edges in different
        // layers).
        GraphOptions graph_options = new();  // Default options keep everything.
        Options options = new(new IntLatLngSnapFunction(0))
        {
            SimplifyEdgeChains = true
        };

        // Note below that the edge chains in different layers have different vertex
        // locations, different number of interior vertices, different degenerate
        // edges, etc, and yet they can all be simplified together.
        var input = new List<List<string>> { new() {
  "0.1:5, 0:5.2", "0.1:0, 0:9.9",   // Defined before chain
  "0:10.1, 0:0.1", "0:3.1, 0:2.9",  // Defined after chain
}, new() {
  "0.1:3, 0:3.2", "-0.1:0, 0:4.1, 0:9.9",  // Defined before chain
  "0.1:9.9, 0:7, 0.1:6.9, 0.1:0.2",        // Defined inside chain
}, new() {
  "0.2:0.3, 0.1:6, 0:5.9, 0.1:10.2",       // Defined inside chain
  "0.1:0.1, 0:9.8", "0.1:2, 0:2.1",        // Defined after chain
}
};
        var expected = new List<EdgeInputEdgeIds> {
  new() {
  ("0:0, 0:10", new[]{0, 1}), ("0:10, 0:0", new[]{2, 3})
}, new() {
  ("0:0, 0:10", new[]{4, 5, 6}), ("0:10, 0:0", new[]{7, 8, 9})
}, new() {
  ("0:0, 0:10", new[]{10, 11, 12}), ("0:0, 0:10", new[]{13, 14})
}
};
        S2Builder builder = new(options);
        for (int i = 0; i < input.Count; ++i)
        {
            builder.StartLayer(new InputEdgeIdCheckingLayer(expected[i], graph_options));
            foreach (var input_str in input[i])
            {
                builder.AddPolyline(MakePolylineOrDie(input_str));
            }
        }
        Assert.True(builder.Build(out _));
    }

    [Fact]
    internal void Test_S2Builder_HighPrecisionPredicates()
    {
        // To produce correct output in this example, the algorithm needs fall back
        // to high precision predicates when the output of the normal predicates is
        // uncertain.
        var vertices = new S2Point[] {
new(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
new(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508),
new(-0.10531192039116592, -0.80522217309701472, 0.58354661457028933),
};
        S2Polyline input = new(vertices);
        S1Angle snap_radius = S2.kIntersectionMergeRadiusS1Angle;
        Options options = new(new IdentitySnapFunction(snap_radius))
        {
            Idempotent = false
        };
        S2Builder builder = new(options);
        S2Polyline output = new();
        builder.StartLayer(new S2PolylineLayer(output));
        builder.ForceVertex(new S2Point(
            -0.10531192039134191, -0.80522217309705857, 0.58354661457019719));
        builder.AddPolyline(input);
        Assert.True(builder.Build(out _));
    }

    [Fact]
    internal void Test_S2Builder_HighPrecisionStressTest()
    {
        // This testructs many small, random inputs such that the output is
        // likely to be inconsistent unless high-precision predicates are used.

        S1Angle snap_radius = S2.kIntersectionMergeRadiusS1Angle;
        // Some S2Builder calculations use an upper bound that takes into account
        // S1ChordAngle errors.  We sometimes try perturbing points by very close to
        // that distance in an attempt to expose errors.
        S1ChordAngle ca = new(snap_radius);
        S1Angle snap_radius_with_error = ca.PlusError(
            ca.S1AngleConstructorMaxError +
            S2.GetUpdateMinDistanceMaxError(ca)).ToAngle();

        int non_degenerate = 0;
        int kIters = 8000 * iteration_multiplier;
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.

            // We construct a nearly degenerate triangle where one of the edges is
            // sometimes very short.  Then we add a forced vertex somewhere near the
            // shortest edge.  Then after snapping, we check that (1) the edges still
            // form a loop, and (2) if the loop is non-degenerate, then it has the
            // same orientation as the original triangle.
            //
            // v1 is located randomly.  (v0,v1) is the longest of the three edges.
            // v2 is located along (v0,v1) but is perturbed by up to 2 * snap_radius.
            S2Point v1 = ChoosePoint(), v0_dir = ChoosePoint();
            double d0 = Math.Pow(1e-16, S2Testing.Random.RandDouble());
            S2Point v0 = S2.GetPointOnLine(v1, v0_dir, S1Angle.FromRadians(d0));
            double d2 = 0.5 * d0 * Math.Pow(1e-16, Math.Pow(S2Testing.Random.RandDouble(), 2));
            S2Point v2 = S2.GetPointOnLine(v1, v0_dir, S1Angle.FromRadians(d2));
            v2 = S2Testing.SamplePoint(new S2Cap(v2, 2 * snap_radius));
            // Vary the edge directions by randomly swapping v0 and v2.
            if (S2Testing.Random.OneIn(2)) { var tmp = v0; v0 = v2; v2 = tmp; }

            // The forced vertex (v3) is either located near the (v1, v2) edge.
            // We perturb it either in a random direction from v1 or v2, or
            // perpendicular to (v1, v2) starting from an interior edge point.
            S1Angle d3 = S2Testing.Random.OneIn(2) ? snap_radius : snap_radius_with_error;
            if (S2Testing.Random.OneIn(3)) d3 = 1.5 * S2Testing.Random.RandDouble() * d3;
            S2Point v3;
            if (S2Testing.Random.OneIn(5))
            {
                v3 = S2Testing.Random.OneIn(2) ? v1 : v2;
                v3 = S2.GetPointOnLine(v3, ChoosePoint(), d3);
            }
            else
            {
                v3 = S2.Interpolate(v1, v2, Math.Pow(1e-16, S2Testing.Random.RandDouble()));
                v3 = S2.GetPointOnLine(v3, v1.CrossProd(v2).Normalize(), d3);
            }
            Options options = new(new IdentitySnapFunction(snap_radius))
            {
                Idempotent = false
            };
            S2Builder builder = new(options);
            S2Polygon output = new()
            {
                S2DebugOverride = S2Debug.DISABLE
            };
            builder.StartLayer(new S2PolygonLayer(output));
            builder.ForceVertex(v3);
            builder.AddEdge(v0, v1);
            builder.AddEdge(v1, v2);
            builder.AddEdge(v2, v0);
            Assert.False(builder.Build(out var error));
            if (error.IsOk() && !output.IsEmpty())
            {
                Assert.Equal(1, output.NumLoops());
                if (output.NumLoops() == 1)
                {
                    Assert.True(output.IsValid());
                    Assert.Equal(S2Pred.Sign(v0, v1, v2) > 0, output.Loop(0).IsNormalized());
                    // "d0=" << d0 << ", d2=" << d2 << ", d3=" << d3;
                    ++non_degenerate;
                }
            }
        }
        Assert.True(non_degenerate >= kIters / 10);
    }

    [Fact]
    internal void Test_S2Builder_SelfIntersectionStressTest()
    {
        int kIters = 50 * iteration_multiplier;
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            var dt = DateTime.UtcNow;

            // The minimum radius is about 36cm on the Earth's surface.  The
            // performance is reduced for radii much smaller than this because
            // S2ShapeIndex only indexes regions down to about 1cm across.
            S2Cap cap = S2Testing.GetRandomCap(1e-14, 1e-2);

            Options options = new()
            {
                SplitCrossingEdges = true
            };
            if (S2Testing.Random.OneIn(2))
            {
                S1Angle radius = cap.RadiusAngle();
                int min_exp = IntLatLngSnapFunction.ExponentForMaxSnapRadius(radius);
                int exponent = Math.Min(IntLatLngSnapFunction.kMaxExponent,
                                   min_exp + S2Testing.Random.Uniform(5));
                options.SnapFunction = new IntLatLngSnapFunction(exponent);
            }
            S2Builder builder = new(options);

            // Note that the number of intersections (and the running time) is
            // quadratic in the number of vertices.  With 200 input vertices, the
            // output consists of about 2300 loops and 9000 vertices.
            S2Polygon output = new();
            builder.StartLayer(new S2PolygonLayer(
                output, new S2PolygonLayer.Options(EdgeType.UNDIRECTED)));
#if DEBUG
            var vertices = new S2Point[50];
#else
            var vertices = new S2Point[200];
#endif
            for (var i = 0; i < vertices.Length; i++)
            {
                vertices[i] = S2Testing.SamplePoint(cap);
            }
            vertices[^1] = vertices.First();
            S2Polyline input = new(vertices);
            builder.AddPolyline(input);
            Assert.True(builder.Build(out _));
            Assert.False(output.FindValidationError(out _));
            if (iter == -1)
            {
                _logger.WriteLine($"S2Polyline: {input.ToDebugString()}");
                _logger.WriteLine($"S2Polyline: {output.ToDebugString()}");
            }
            if (iter < 50)
            {
                var dif = DateTime.UtcNow.Subtract(dt).TotalMilliseconds;
                _logger.WriteLine($"iter={iter:d4}: ms={dif:d5}, radius={cap.RadiusAngle().Radians:g8.3}, loops={output.NumLoops()}, vertices={output.NumVertices}");
            }
        }
    }

    [Fact]
    internal void Test_S2Builder_FractalStressTest()
    {
#if DEBUG
        int kIters = 100 * iteration_multiplier;
        int levelForApproxMaxEdges = 800;
#else
        int kIters = 1000 * iteration_multiplier;
        int levelForApproxMaxEdges = 12800;
#endif
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
            S2Testing.Fractal fractal = new();
            fractal.SetLevelForApproxMaxEdges(levelForApproxMaxEdges);
            fractal.SetLevelForApproxMinEdges(12);
            fractal.FractalDimension = 1.5 + 0.5 * S2Testing.Random.RandDouble();
            S2Polygon input = new(fractal.MakeLoop(S2Testing.GetRandomFrame(),
                                             S1Angle.FromDegrees(20)));
            Options options = new();
            if (S2Testing.Random.OneIn(3))
            {
                int exponent = S2Testing.Random.Uniform(11);
                options.SnapFunction = new IntLatLngSnapFunction(exponent);
            }
            else if (S2Testing.Random.OneIn(2))
            {
                int level = S2Testing.Random.Uniform(20);
                options.SnapFunction = new S2CellIdSnapFunction(level);
            }
            else
            {
                options.SnapFunction = new IdentitySnapFunction(
                    S1Angle.FromDegrees(10 * Math.Pow(1e-4, S2Testing.Random.RandDouble())));
            }
            S2Builder builder = new(options);
            S2Polygon output = new();
            builder.StartLayer(new S2PolygonLayer(output));
            builder.AddPolygon(input);
            Assert.True(builder.Build(out _));
            Assert.False(output.FindValidationError(out _));
            if (iter == -1)
            {
                _logger.WriteLine($"S2Polygon: {input.ToDebugString()}");
                _logger.WriteLine($"S2Polygon: {output.ToDebugString()}");
            }
            if (iter < 50)
            {
                _logger.WriteLine($"iter={iter:d4}: in_vertices={input.NumVertices}, out_vertices={output.NumVertices}");
            }
        }
    }

    [Fact]
    internal void Test_S2Builder_AdjacentCoverageIntervalsSpanMoreThan90Degrees()
    {
        // The test for whether one Voronoi site excludes another along a given
        // input edge boils down to a test of whether two angle intervals "a" and
        // "b" overlap.  Let "ra" and "rb" be the semi-widths of the two intervals,
        // and let "d" be the angle between their centers.  Then "a" contains "b" if
        // (rb + d <= ra), and "b" contains "a" if (rb - d >= ra).  However the
        // actual code uses the sines of the angles, e.g.  sin(rb + d) <= sin(ra).
        // This works fine most of the time, but the first condition (rb + d <= ra)
        // also needs to check that rb + d < 90 degrees.  This test verifies that
        // case.

        // The following 3 tests have d < 90, d = 90, and d > 90 degrees, but in all
        // 3 cases rb + d > 90 degrees.
        TestSnappingWithForcedVertices("0:0, 0:80", S1Angle.FromDegrees(60),
                                       "0:0, 0:70", "0:0, 0:70");
        TestSnappingWithForcedVertices("0:0, 0:80", S1Angle.FromDegrees(60),
                                       "0:0, 0:90", "0:0, 0:90");
        TestSnappingWithForcedVertices("0:0, 0:80", S1Angle.FromDegrees(60),
                                       "0:0, 0:110", "0:0, 0:110");

        // This test has d = 180 degrees, i.e. the two sites project to points that
        // are 180 degrees apart along the input edge.  The snapped edge doesn't
        // stay within max_edge_deviation() of the input edge, so an extra site is
        // added and it is snapped again (yielding two edges).  The case we are
        // testing here is the first call to SnapEdge() before adding the site.
        TestSnappingWithForcedVertices("0:10, 0:170", S1Angle.FromDegrees(50),
                                       "47:0, 49:180", "47:0, 0:90, 49:180");

        // This test has d = 220 degrees, i.e. when the input edge is snapped it
        // goes the "wrong way" around the sphere.  Again, the snapped edge is too
        // far from the input edge so an extra site is added and it is resnapped.
        TestSnappingWithForcedVertices("0:10, 0:170", S1Angle.FromDegrees(70),
                                       "0:-20, 0:-160", "0:-20, 0:90, 0:-160");

        // Without using forced vertices, the maximum angle between the coverage
        // interval centers is d = 300 degrees.  This would use an edge 180 degrees
        // long, and then place two sites 60 degrees past either endpoint.  With
        // forced vertices we can increase the snap radius to 70 degrees and get an
        // angle of up to d = 320 degrees, but the sites are only 40 degrees apart
        // (which is why it requires forced vertices).  The test below is an
        // approximation of this situation with d = 319.6 degrees.
        TestSnappingWithForcedVertices("0:0.1, 0:179.9", S1Angle.FromDegrees(70),
                                       "0:-69.8, 0:-110.2",
                                       "0:-69.8, 0:90, 0:-110.2");

        static void TestSnappingWithForcedVertices(string input_str,
            S1Angle snap_radius, string vertices_str, string expected_str)
        {
            S2Builder builder = new(new Options(new IdentitySnapFunction(snap_radius)));
            var vertices = ParsePointsOrDie(vertices_str);
            foreach (var vertex in vertices)
            {
                builder.ForceVertex(vertex);
            }
            S2Polyline output = new();
            builder.StartLayer(new S2PolylineLayer(output));
            builder.AddPolyline(MakePolylineOrDie(input_str));
            Assert.True(builder.Build(out _));
            Assert.Equal(expected_str, output.ToDebugString());
        }
    }

    // The following test requires internal debugging checks to be skipped.
#if true///NDEBUG

    [Fact]
    internal void Test_S2Builder_NaNVertices()
    {
        // Test that S2Builder operations don't crash when some vertices are NaN.
        List<List<S2Point>> loops = [
new(){ new(double.NaN, double.NaN, double.NaN), new(double.NaN, double.NaN, double.NaN), new(double.NaN, double.NaN, double.NaN) }
];
        S2LaxPolygonShape output = new();
        S2Builder builder = new(
            new Options(new IdentitySnapFunction(S1Angle.FromRadians(1e-15)))
        );
        builder.StartLayer(new LaxPolygonLayer(output));
        builder.AddShape(new S2LaxPolygonShape(loops));

        // The operation fails with S2Error::BUILDER_SNAP_RADIUS_TOO_SMALL because
        // the distance between two NaN vertices happens to be computed as Pi.  We
        // aren't concerned about the actual error (or even whether there is an
        // error), we simply don't want to crash in this case.
        Assert.False(builder.Build(out S2Error error));
        Assert.Equal(output.NumLoops, 0);
    }
#endif

    [Fact]
    internal void Test_S2Builder_OldS2PolygonBuilderBug()
    {
        // This is a polygon that caused the obsolete S2PolygonBuilder class to
        // generate an invalid output polygon (duplicate edges).
        S2Polygon input = MakePolygonOrDie(
            "32.2983095:72.3416582, 32.2986281:72.3423059, " +
            "32.2985238:72.3423743, 32.2987176:72.3427807, " +
            "32.2988174:72.3427056, 32.2991269:72.3433480, " +
            "32.2991881:72.3433077, 32.2990668:72.3430462, " +
            "32.2991745:72.3429778, 32.2995078:72.3436725, " +
            "32.2996075:72.3436269, 32.2985465:72.3413832, " +
            "32.2984558:72.3414530, 32.2988015:72.3421839, " +
            "32.2991552:72.3429416, 32.2990498:72.3430073, " +
            "32.2983764:72.3416059");
        Assert.True(input.IsValid());

        S1Angle snap_radius = S2Testing.MetersToAngle(20 / 0.866);
        S2Builder builder = new(new Options(new IdentitySnapFunction(snap_radius)));
        S2Polygon output = new();
        builder.StartLayer(new S2PolygonLayer(output));
        builder.AddPolygon(input);
        Assert.True(builder.Build(out _));
        Assert.True(output.IsValid());
        S2Polygon expected = MakePolygonOrDie(
            "32.2991552:72.3429416, 32.2991881:72.3433077, 32.2996075:72.3436269; " +
            "32.2988015:72.3421839, 32.2985465:72.3413832, 32.2983764:72.3416059, " +
            "32.2985238:72.3423743, 32.2987176:72.3427807");
        ExpectPolygonsEqual(expected, output);
    }

    [Fact]
    internal void Test_S2Builder_SeparationSitesRegressionBug()
    {
        // This test uses a snap radius of zero with edge splitting, and has some very
        // closely spaced vertices (much smaller than S2::kIntersectionError).  It
        // used to fail (on the PPC architecture only) because there was a missing
        // test in S2Builder::MaybeAddExtraSites().
        S2Builder.Options options = new()
        {
            SplitCrossingEdges = true
        };
        S2Builder builder = new(options);
        S2PolylineVectorLayer.Options layer_options = new()
        {
            PolylineType_ = Graph.PolylineType.WALK
        };
        List<S2Polyline> output = [];
        builder.StartLayer(
            new S2PolylineVectorLayer(output, layer_options));
        List<List<S2Point>> input_polylines = [
            new(){new(0.99482894039096326, 0.087057485575229562, 0.05231035811301657),
                  new(0.19008255728509718, 0.016634125542513145, 0.98162718344766398)},
            new(){new(0.99802098666373784, 0.052325259429907504, 0.034873735164620751),
                  new(0.99585181570926085, 0.087146997393412709, 0.026164135641767797),
                  new(0.99939172130835197, 6.9770704216017258e-20, 0.034873878194564757),
                  new(0.99939172130835197, 1.7442676054004314e-202, 0.034873878194564757),
                  new(0.99939172130835197, 2.4185105853059967e-57, 0.034873878194564757),
                  new(0.99939091697091686, 0, 0.034896920724182809),
                  new(0.99543519482327569, 0.088840224357046416, 0.034873879097925588)},
            new(){new(-0.86549861898490243, 0.49969586065415578, 0.034873878194564757),
                  new(0.99939172130835197, 1.542605867912342e-181, 0.034873878194564757),
                  new(0.99939172130835197, 1.5426058679123417e-281, 0.034873878194564757),
                  new(0.99939172130835197, 1.5426058504696658e-231, 0.034873878194564757),
                  new(0.19080899537654492, 3.3302452117433465e-113, 0.98162718344766398)},
            new(){new(0.99802098660295513, 0.052325259426720727, 0.034873736908888363),
                  new(0.99558688908226523, 0.08712381366290145, 0.034873878194564757),
                  new(0.99939172130835197, 1.0221039496805218e-23, 0.034873878194564757),
                  new(0.99939172127682907, 3.4885352106908273e-20, 0.034873879097925602),
                  new(0.99391473614090387, 0.10448593114531293, 0.03487387954694085)}
        ];
        foreach (var polyline in input_polylines) {
            for (int i = 0; i + 1 < polyline.Count; ++i)
            {
                builder.AddEdge(polyline[i], polyline[i + 1]);
            }
        }
        Assert.True(builder.Build(out S2Error error));
    }

    [Fact]
    internal void Test_S2Builder_VoronoiSiteExclusionBug1()
    {
        // This reproduces a former bug where the input edge could sometimes snap
        // incorrectly to an extra vertex.  This could only happen when the sum of
        // the edge length and the snap radius is more than 180 degrees.
        S2Builder.Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(64.83)));
        S2Builder builder = new(options);
        S2Polyline output = new();
        builder.StartLayer(new S2PolylineLayer(output));
        // The first vertex of this edge snaps to the first forced vertex below.
        // The edge should not snap to the second forced vertex.
        builder.AddPolyline(MakePolylineOrDie("29.40:173.03, -18.02:-5.83"));
        builder.ForceVertex(MakePointOrDie("25.84:131.46"));
        builder.ForceVertex(MakePointOrDie("-29.23:-166.58"));
        Assert.True(builder.Build(out _));
        var expected = "25.84:131.46, -18.02:-5.83";
        Assert.Equal(expected, S2TextFormat.ToDebugString(output));
    }

    [Fact]
    internal void Test_S2Builder_VoronoiSiteExclusionBug2()
    {
        // This reproduces a former bug where the input edge could sometimes snap
        // incorrectly to an extra vertex.  This could only happen when the sum of
        // the edge length and the snap radius is more than 180 degrees.
        S2Builder.Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(67.75)));
        S2Builder builder = new(options);
        S2Polyline output=new();
        builder.StartLayer(new S2PolylineLayer(output));
        builder.AddPolyline(MakePolylineOrDie("47.06:-175.17, -47.59:10.57"));
        builder.ForceVertex(MakePointOrDie("36.36:47.63"));
        builder.ForceVertex(MakePointOrDie("-28.34:-72.46"));
        Assert.True(builder.Build(out _));
        // Snapping to the given vertices would cause the snapped edge to deviate
        // too far from the input edge, so S2Builder adds an extra site.  Given the
        // new site, snapping to the
        var expected = "47.06:-175.17, -34.4968065428191:69.7125289482374";
        Assert.Equal(expected, S2TextFormat.ToDebugString(output));
    }

    [Fact]
    internal void Test_S2Builder_HausdorffDistanceBug()
    {
        // This reproduces a former bug where a very long edge (close to 180
        // degrees) could sometimes snap incorrectly due to a bug in
        // S2::IsEdgeBNearEdgeA().
        S2Builder.Options options = new(new IdentitySnapFunction(S1Angle.FromDegrees(70)));
        S2Builder builder= new(options);
        S2Polygon output=new();
        builder.StartLayer(new S2PolygonLayer(output));
        builder.AddShape(MakeLaxPolygonOrDie(
            "35:17; -40:88, 68:-161, 48:-156, -45:-10, -40:88"));
        Assert.True(builder.Build(out _));
        Assert.Equal(output.NumLoops(), 1);
    }

    [Fact]
    internal void Test_S2Builder_IncorrectSeparationSiteBug()
    {
        // This reproduces a bug that used to attempt to create a separation site
        // where none was necessary.
        S2Builder.Options options = new()
        {
            Idempotent = false,
            SplitCrossingEdges = true
        };
        S2Builder builder= new (options);
        S2Polyline output = new();
        builder.StartLayer(new S2PolylineLayer(output));
        builder.AddEdge(new S2Point(-0.50094438964076704, -0.86547947317509455, 0),
                        new S2Point(1, 1.7786363250284876e-322, 4.7729929394856611e-65));
        builder.ForceVertex(new S2Point(1, 0, -4.7729929394856611e-65));
        builder.ForceVertex(
            new S2Point(1, 2.2603503297237029e-320, 4.7729929394856619e-65));
        Assert.True(builder.Build(out _));
    }

    [Fact]
    internal void Test_S2Builder_PushPopLabel()
    {
        // TODO(b/232074544): Test more thoroughly.
        S2Builder builder=new();
        builder.PushLabel(1);
        builder.PopLabel();
    }

    private static void ExpectPolygonsEqual(S2Polygon expected, S2Polygon actual)
    {
        Assert.True(expected == actual);
    }

    private static void ExpectPolygonsApproxEqual(S2Polygon expected, S2Polygon actual, S1Angle tolerance)
    {
        Assert.True(expected.BoundaryApproxEquals(actual, tolerance));
        // expected.ToDebugString()
        // actual.ToDebugString()
        // tolerance.Degrees
    }

    private static void ExpectPolylinesEqual(S2Polyline expected, S2Polyline actual)
    {
        Assert.True(expected == actual);
    }

    // Verifies that two graphs have the same vertices and edges.
    private static void ExpectGraphsEqual(Graph expected, Graph actual)
    {
        Assert.Equal(expected.Vertices, actual.Vertices);
        Assert.Equal(expected.Edges, actual.Edges);
        Assert.Equal(expected.InputEdgeIdSetIds, actual.InputEdgeIdSetIds);
    }

    private static void TestPolylineLayers(List<string> input_strs, List<string> expected_strs,
        S2PolylineLayer.Options layer_options, Options? builder_options = null)
    {
        S2Builder builder = new(builder_options ?? new Options());
        List<S2Polyline> output = [];
        foreach (var input_str in input_strs)
        {
            output.Add(new S2Polyline());
            builder.StartLayer(new S2PolylineLayer(output.Last(), layer_options));
            builder.AddPolyline(MakePolylineOrDie(input_str));
        }
        Assert.True(builder.Build(out _));
        List<string> output_strs = [];
        foreach (var polyline in output)
        {
            output_strs.Add(polyline.ToDebugString());
        }
        Assert.Equal(string.Join("; ", expected_strs),
                  string.Join("; ", output_strs));
    }

    private static void TestPolylineVector(List<string> input_strs, List<string> expected_strs,
        S2PolylineVectorLayer.Options layer_options, Options? builder_options = null)
    {
        S2Builder builder = new(builder_options ?? new Options());
        List<S2Polyline> output = [];
        builder.StartLayer(
            new S2PolylineVectorLayer(output, layer_options));
        foreach (var input_str in input_strs)
        {
            builder.AddPolyline(MakePolylineOrDie(input_str));
        }
        Assert.True(builder.Build(out _));
        List<string> output_strs = [];
        foreach (var polyline in output)
        {
            output_strs.Add(polyline.ToDebugString());
        }
        Assert.Equal(string.Join("; ", expected_strs),
                  string.Join("; ", output_strs));
    }

    private static void TestPolylineLayersBothEdgeTypes(List<string> input_strs,
        List<string> expected_strs, S2PolylineLayer.Options layer_options, Options? builder_options = null)
    {
        builder_options ??= new Options();
        layer_options.EdgeType = EdgeType.DIRECTED;
        TestPolylineLayers(input_strs, expected_strs, layer_options, builder_options);
        layer_options.EdgeType = EdgeType.UNDIRECTED;
        TestPolylineLayers(input_strs, expected_strs, layer_options, builder_options);
    }

    private static void TestInputEdgeIds(List<string> input_strs,
        EdgeInputEdgeIds expected, GraphOptions graph_options, Options options)
    {
        S2Builder builder = new(options);
        builder.StartLayer(new InputEdgeIdCheckingLayer(expected, graph_options));
        foreach (var input_str in input_strs)
        {
            builder.AddPolyline(MakePolylineOrDie(input_str));
        }
        Assert.True(builder.Build(out _));
    }

    // Chooses a random S2Point that is often near the intersection of one of the
    // coodinates planes or coordinate axes with the unit sphere.  (It is possible
    // to represent very small perturbations near such points.)
    private static S2Point ChoosePoint()
    {
        S2Point x = S2Testing.RandomPoint();
        for (int i = 0; i < 3; ++i)
        {
            if (S2Testing.Random.OneIn(3))
            {
                x.SetAxis(i, x[i] * Math.Pow(1e-50, S2Testing.Random.RandDouble()));
            }
        }
        return x.Normalize();
    }

    private class InputEdgeIdCheckingLayer : Layer
    {
        internal InputEdgeIdCheckingLayer(EdgeInputEdgeIds expected, GraphOptions graph_options)
        { expected_ = expected; graph_options_ = graph_options; }

        public override GraphOptions GraphOptions_()
        { return graph_options_; }
        private readonly GraphOptions graph_options_;

        public override void Build(Graph g, out S2Error error)
        {
            EdgeInputEdgeIds actual = [];
            List<S2Point> vertices = [];
            for (var e = 0; e < g.NumEdges; ++e)
            {
                vertices.Clear();
                vertices.Add(g.Vertex(g.GetEdge(e).ShapeId));
                vertices.Add(g.Vertex(g.GetEdge(e).EdgeId));
                string edge = new S2Point[]{
                g.Vertex(g.GetEdge(e).ShapeId),
                    g.Vertex(g.GetEdge(e).EdgeId)}.ToDebugString();
                var ids = g.InputEdgeIds(e).ToArray();
                actual.Add((edge, ids));
            }
            // This comparison doesn't consider multiplicity, but that's fine.
            StringBuilder missing = new(), extra = new();
            foreach (var p in expected_)
            {
                if (actual.Any(t => t == p)) continue;
                missing.Append(ToDebugString(p));
            }
            foreach (var p in actual)
            {
                if (expected_.Any(t => t == p)) continue;
                extra.Append(ToDebugString(p));
            }
            if (missing.Length > 0 || extra.Length > 0)
            {
                error = new(INPUT_EDGE_ID_MISMATCH, $"Missing:\n{missing}Extra:\n{extra}\n");
            }
            else
            {
                error = S2Error.OK;
            }
        }

        private static string ToDebugString((string, int[]) p)
        {
            var sb = new StringBuilder($"  ({p.Item1})={{");
            if (p.Item2.Length!=0)
            {
                foreach (int id in p.Item2)
                {
                    sb.Append($"{id}, ");
                }
                sb.Remove(sb.Length - 2, 2);
            }
            sb.Append("}\n");
            return sb.ToString();
        }

        private readonly EdgeInputEdgeIds expected_;
    }

    // This layer makes both a shallow and a deep copy of the Graph object passed
    // to its Build() method and appends them to two vectors.  Furthermore, it
    // verifies that the shallow and deep copies of any graphs previously appended
    // to those vectors are still identical.
    private class GraphPersistenceLayer(GraphOptions graph_options, List<Graph> graphs, List<GraphClone> clones) : Layer
    {
        public override GraphOptions GraphOptions_() => graph_options_;
        private readonly GraphOptions graph_options_ = graph_options;

        public override void Build(Graph g, out S2Error error)
        {
            error = S2Error.OK;
            // Verify that all graphs built so far are unchanged.
            for (int i = 0; i < graphs_.Count; ++i)
            {
                ExpectGraphsEqual(clones_[i].Graph(), graphs_[i]);
            }
            graphs_.Add(g);
            clones_.Add(new GraphClone(g));
        }

        private readonly List<Graph> graphs_ = graphs;       // Shallow copies.
        private readonly List<GraphClone> clones_ = clones;  // Deep copies.
    }
}
