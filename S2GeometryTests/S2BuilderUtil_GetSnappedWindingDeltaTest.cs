global using S2EdgeCrosser = S2Geometry.S2EdgeCrosserBase;
global using S2CopyingEdgeCrosser = S2Geometry.S2EdgeCrosserBase;

namespace S2Geometry;

using DegenerateEdges = GraphOptions.DegenerateEdges;
using DuplicateEdges = GraphOptions.DuplicateEdges;
using SiblingPairs = GraphOptions.SiblingPairs;

using InputEdgeFilter = Func<InputEdgeId, bool>;

// Used to build a histogram of winding numbers.
using WindingTally = Dictionary<int, int>;

public class S2BuilderUtil_GetSnappedWindingDeltaTest
{
    // This S2Builder layer simply calls s2builderutil::GetSnappedWindingDelta()
    // with the given "ref_input_edge_id" and compares the result to
    // "expected_winding_delta".
    public class WindingNumberComparingLayer : S2Builder.Layer 
    {
        public WindingNumberComparingLayer(InputEdgeId ref_input_edge_id,
                                       S2Builder builder,
                                       int expected_winding_delta)
        {
                ref_input_edge_id_ = ref_input_edge_id;
                builder_ = builder;
                expected_winding_delta_ = expected_winding_delta;
        }

        public override GraphOptions GraphOptions_()
        {
            return new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                DuplicateEdges.MERGE, SiblingPairs.KEEP);
        }

        public override void Build(Graph g, out S2Error error)
        {
            // Find the reference vertex before and after snapping.
            var ref_in = builder_.InputEdge(ref_input_edge_id_).V0;
            var ref_v = SnappedWindingDelta.FindFirstVertexId(ref_input_edge_id_, g);
            Assert.True(ref_v >= 0);
            int winding_delta = SnappedWindingDelta.GetSnappedWindingDelta(
                ref_in, ref_v, null, builder_, g, out error);
            Assert.True(error.IsOk());
            Assert.Equal(winding_delta, expected_winding_delta_);
        }

        private readonly InputEdgeId ref_input_edge_id_;
        private readonly S2Builder builder_;
        private readonly int expected_winding_delta_;
    }


    // Given a set of loops, a set of forced vertices, and a snap radius in
    // degrees, verifies that the change in winding number computed by
    // s2builderutil::GetSnappedWindingDelta() for the degenerate edge
    // "ref_input_edge_id" is "expected_winding_delta".
    public static void ExpectWindingDelta(
                    string loops_str, string forced_vertices_str,
                    double snap_radius_degrees, InputEdgeId ref_input_edge_id,
                    int expected_winding_delta)
{
    S2Builder builder= new(new Options(new IdentitySnapFunction(
        S1Angle.FromDegrees(snap_radius_degrees))));
    builder.StartLayer(new WindingNumberComparingLayer(
        ref_input_edge_id, builder, expected_winding_delta));
    foreach (var v in S2TextFormat.ParsePointsOrDie(forced_vertices_str)) {
    builder.ForceVertex(v);
}
builder.AddShape(S2TextFormat.MakeLaxPolygonOrDie(loops_str));
var ref_edge = builder.InputEdge(ref_input_edge_id);
        Assert.True(ref_edge.V0 == ref_edge.V1); // "Reference edge not degenerate";
        Assert.True(builder.Build(out _)); // << error;
    }

    // Since the GetSnappedWindingDelta() algorithm doesn't depend on absolute
    // vertex positions, we use "0:0" as the snapped vertex position in the tests
    // below (using ForceVertex() to ensure this).  Most tests use a snap radius
    // of 10 degrees since this makes it convenient to construct suitable tests.
    //
    // FYI the S2::RefDir() direction for "0:0" (used to determine whether a loop
    // contains one of its vertices) is just slightly north of due west,
    // i.e. approximately 179.7 degrees CCW from the positive longitude axis.

    // No edges except the degenerate edges that defines the reference vertex.
    [Fact]
    public void Test_GetSnappedWindingDelta_NoOtherEdges()
    {
        ExpectWindingDelta("0:0", "0:0", 10.0, 0, 0);
}

    // Degenerate input loops.
    [Fact]
    public void Test_GetSnappedWindingDelta_DegenerateInputLoops() {
        ExpectWindingDelta("0:0; 1:1; 2:2", "0:0", 10.0, 0, 0);
}

    // Duplicate degenerate input loops.
    [Fact]
    public void Test_GetSnappedWindingDelta_DuplicateDegenerateInputLoops() {
        ExpectWindingDelta("0:0; 0:0; 1:1; 1:1", "0:0", 10.0, 0, 0);
}

    // A shell around the reference vertex that collapses to a single point.
    [Fact]
    public void Test_GetSnappedWindingDelta_CollapsingShell() {
        ExpectWindingDelta("0:0; 1:1, 1:-2, -2:1", "0:0", 10.0, 0, -1);
}

    // A hole around the reference vertex that collapses to a single point.
    [Fact]
    public void Test_GetSnappedWindingDelta_CollapsingHole() {
        ExpectWindingDelta("0:0; 1:1, -2:1, 1:-2", "0:0", 10.0, 0, +1);
}

    // A single "shell" that winds around the reference vertex twice.
    [Fact]
    public void Test_GetSnappedWindingDelta_CollapsingDoubleShell() {
        ExpectWindingDelta("0:0; 1:1, 1:-2, -2:1, 2:2, 2:-3, -3:2",
                       "0:0", 10.0, 0, -2);
}

    // A loop that enters the Voronoi region of the snapped reference vertex and
    // then leaves again, where the reference vertex is not contained by the loop
    // and does not move during snapping.
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalLoopRefVertexStaysOutside() {
        ExpectWindingDelta("0:0; 20:0, 0:0, 0:20", "0:0", 10.0, 0, 0);
}

    // Like the above, except that the reference vertex is contained by the loop.
    // (See S2::RefDir comments above.)
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalLoopRefVertexStaysInside() {
        ExpectWindingDelta("0:0; 0:-20, 0:0, 20:0", "0:0", 10.0, 0, 0);
}

    // The reference vertex moves from outside to inside an external loop during
    // snapping.
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalLoopRefVertexMovesInside() {
        ExpectWindingDelta("1:1; 0:-20, 1:-1, 20:0", "0:0", 10.0, 0, +1);
}

    // A single loop edge crosses the Voronoi region of the reference vertex and the
    // reference vertex stays outside the loop during snapping.
    [Fact]
    public void Test_GetSnappedWindingDelta_CrossingEdgeRefVertexStaysOutside() {
        ExpectWindingDelta("-1:-1; 20:-20, -20:20, 20:20", "0:0", 10.0, 0, 0);
}

    // A single loop edge crosses the Voronoi region of the reference vertex and the
    // reference vertex moves outside the loop during snapping.
    [Fact]
    public void Test_GetSnappedWindingDelta_CrossingEdgeRefVertexMovesOutside() {
        ExpectWindingDelta("1:1; 20:-20, -20:20, 20:20", "0:0", 10.0, 0, -1);
}

    // An external loop that winds CW around the reference vertex twice, where the
    // reference vertex moves during snapping, and where the reference vertex is
    // outside the loop after snapping (so that its winding number only increases
    // by 1).
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalLoopDoubleHoleToSingleHole() {
        ExpectWindingDelta("4:4; 0:20, 3:3, 6:3, 2:7, 2:2, 2:20", "0:0", 10.0, 0, +1);
}

    // An external loop that winds CW around the reference vertex twice, where the
    // reference vertex moves during snapping, and where the reference vertex is
    // inside the loop after snapping (so that its winding number increases by 3).
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalLoopDoubleHoleToSingleShell() {
        ExpectWindingDelta("4:4; 0:-20, 6:2, 2:6, 2:2, 6:2, 2:6, 2:2, 20:0",
                       "0:0", 10.0, 0, +3);
}

    // This and the following tests vertify that the partial loops formed by the
    // local input and output edges are closed consistently with each other (such
    // that the hypothetical connecting edges can deform from one to the other
    // without passing through the reference vertex).
    //
    // An external loop where the input edges that enter/exit the Voronoi region
    // cross, but the snapped edges do not.  (This can happen when the input edges
    // snap to multiple edges and the crossing occurs outside the Voronoi region
    // of the reference vertex.)  In this particular test, the entering/exiting
    // edges snap to the same adjacent Voronoi site so that the snapped edges form
    // a loop with one external vertex.
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalEdgesCrossSnapToSameVertex() {
        ExpectWindingDelta("1:1; -5:30, 7:-3, -7:-3, 5:30",
                       "0:0, 0:15", 10.0, 0, -1);
}

    // This test is similar except that the entering/exiting edges snap to two
    // different external Voronoi sites.  Again, the input edges cross but the
    // snapped edges do not.
    [Fact]
    public void Test_GetSnappedWindingDelta_ExternalEdgesCrossSnapToDifferentVertices() {
        ExpectWindingDelta("1:1; -5:40, 7:-3, -7:-3, 5:40",
                       "0:0, 6:10, -6:10", 10.0, 0, -1);
}

    // Test cases where the winding numbers of the reference points Za and Zb in
    // the algorithm description change due to snapping.  (The points Za and Zb
    // are the centers of the great circles defined by the first and last input
    // edges.)  For their winding numbers to change, the input loop needs to cross
    // these points as it deforms during snapping.
    //
    // In all the test below we use perturbations of 70:-180, 5:0 as the first
    // input edge (which yields points close to 0:90 as Za) and perturbations of
    // 0:5, 0:110 as the last input edge (which yields points close to 90:0 as Zb).
    // The first/last vertex can be adjusted slightly to control which side of each
    // edge Za/Zb is on.
    [Fact]
    public void Test_GetSnappedWindingDelta_ReferencePointWindingNumbersChange() {
        // Winding number of Za ~= 0.01:90 changes.
        ExpectWindingDelta("1:1; 70:-179.99, 5:0, 0:5, -0.01:110",
                       "0:0, 1:90", 10.0, 0, 0);

        // Winding number of Zb ~= 89.99:90 changes.
        ExpectWindingDelta("1:1; 70:-179.99, 5:0, 0:5, -0.01:110",
                       "0:0, 89:90", 10.0, 0, 0);

        // Winding numbers of Za and Zb both change.
        ExpectWindingDelta("1:1; 70:-179.99, 5:0, 0:5, -0.01:110",
                       "0:0, 1:90, 89:90", 10.0, 0, 0);

        // Winding number of Za ~= -0.01:90 changes in the opposite direction.
        ExpectWindingDelta("1:1; 70:179.99, 5:0, 0:5, 0:110",
                       "0:0, -1:20, 1:90", 10.0, 0, 0);
}

    // This test demonstrates that a connecting vertex may be necessary in order
    // to ensure that the loops L and L' used to compute the change in winding
    // number for reference points Za and Zb do not pass through those points
    // during snapping.  (This can only happen when the edges A0A1 or B0B1 snap
    // to an edge chain longer than 180 degrees, i.e. where the shortest edge
    // between their new endpoints goes the wrong way around the sphere.)
    [Fact]
    public void Test_GetSnappedWindingDelta_ReferenceLoopsTopologicallyConsistent() {
        // A0A1 follows the equator, Za is at the north pole.  A0A1 snaps to the
        // polyline (0:148, 0:74, 44:-39, -31:-48), where the last two vertices are
        // A0' and R' respectively.  (The perpendicular bisector of A0' and R' just
        // barely intersects A0A1 and therefore it snaps to both vertices.)  A
        // connecting vertex is needed between A0 = 0:148 and A0' = 44:-39 to ensure
        // that this edge stays within the snap radius of A0A1.
        ExpectWindingDelta("-45:24; 0:148, 0:0, -31:-48, 44:-39, -59:0",
                       "-31:-48, 44:-39", 60.0, 0, -1);

        // This tests the symmetric case where a connecting vertex is needed between
        // B1 and B1' to ensure that the edge stays within the snap radius of B0B1.
        ExpectWindingDelta("-45:24;  -59:0, 44:-39, -31:-48, 0:0, 0:148",
                       "-31:-48, 44:-39", 60.0, 0, 1);
}

    // A complex example with multiple loops that combines many of the situations
    // tested individually above.
    [Fact]
    public void Test_GetSnappedWindingDelta_ComplexExample() {
        ExpectWindingDelta("1:1; "+  
                       "70:179.99, 5:0, 0:5, 0:110; "+  
                       "70:179.99, 0:0, 0:3, 3:0, 0:-1, 0:110; "+  
                       "10:-10, -10:10, 10:10; "+  
                       "2:2, 1:-2, -1:2, 2:2, 1:-2, -1:2 ",
                       "0:0, -1:90, 1:90, 45:-5", 10.0, 0, -5);
}

    // This test demonstrates the necessity of the algorithm step that reverses
    // the sign of Za, Zb if necessary to point away from the Voronoi site R'.
    // Further examples can be generated by running RandomLoops test below for
    // enough iterations.
    [Fact]
    public void Test_GetSnappedWindingDelta_EnsureZaZbNotInVoronoiRegion() {
        ExpectWindingDelta(
        "30:42, 30:42; -27:52, 66:131, 30:-93", "", 67.0, 0, -1);
}

    // This test demonstrates the necessity of closing the "chain_diff" loop used
    // by the algorithm.  Further examples can be found by running RandomLoops.
    [Fact]
    public void Test_GetSnappedWindingDelta_EnsureChainDiffLoopIsClosed() {
        ExpectWindingDelta(
        "8:26, 8:26; -36:70, -64:-35, -41:48", "", 66, 0, 0);
}

    // This test previously failed due to a bug in GetVoronoiSiteExclusion()
    // involving long edges (near 180 degrees) and large snap radii.
    [Fact]
    public void Test_GetSnappedWindingDelta_VoronoiExclusionBug() {
        ExpectWindingDelta(
        "24.97:102.02, 24.97:102.02; "+  
        "25.84:131.46, -29.23:-166.58, 29.40:173.03, -18.02:-5.83",
        "", 64.83, 0, -1);
}

// This S2Builder::Layer checks that the change in winding number due to
// snapping computed by s2builderutil::GetSnappedWindingDelta() is correct for
// the given configuration of input edges.
//
// "ref_input_edge_id" should be a degenerate edge SS that specifies the
// reference vertex R whose change in winding number is verified.
//
// "isolated_input_edge_id" should be a degenerate edge II that is not
// expected to snap together with any other edges.  (This ensures that the
// winding number of its vertex I does not change due to snapping.)  I should
// be chosen to be as far away as possible from other vertices and edges used
// for testing purposes.  If more than one edge happens to snap to this
// vertex, S2Error::FAILED_PRECONDITION is returned.
public class WindingNumberCheckingLayer : S2Builder.Layer
{
public  WindingNumberCheckingLayer(InputEdgeId ref_input_edge_id,
                                      InputEdgeId isolated_input_edge_id,
                                      S2Builder builder,
                                      WindingTally winding_tally)
       {
            ref_input_edge_id_ = ref_input_edge_id;
            isolated_input_edge_id_ = isolated_input_edge_id;
            builder_ = builder; winding_tally_ = winding_tally;
}

public override GraphOptions GraphOptions_() 
        {
    // Some of the graph options are chosen randomly.
    return new GraphOptions(
        EdgeType.DIRECTED, DegenerateEdges.KEEP,
      S2Testing.Random.OneIn(2) ? DuplicateEdges.KEEP : DuplicateEdges.MERGE,
      S2Testing.Random.OneIn(2) ? SiblingPairs.KEEP : SiblingPairs.CREATE);
  }

public override void Build(Graph g, out S2Error error)
{
    // First we locate the vertices R, R', I, I'.
    S2Point ref_in = builder_.InputEdge(ref_input_edge_id_).V0;
    VertexId ref_v = SnappedWindingDelta.FindFirstVertexId(ref_input_edge_id_, g);
    S2Point ref_out = g.Vertex(ref_v);

    S2Point iso_in = builder_.InputEdge(isolated_input_edge_id_).V0;
    VertexId iso_v = SnappedWindingDelta.FindFirstVertexId(isolated_input_edge_id_, g);
    S2Point iso_out = g.Vertex(iso_v);

    // If more than one edge snapped to the isolated vertex I, we skip this test
    // since we have no way to  independently verify correctness of the results.
    Graph.VertexOutMap out_map=new(g);
    if (out_map.Degree(iso_v) != 1 ||
        g.InputEdgeIds(out_map.EdgeIds(iso_v).First()).Count != 1)
    {
        error = new S2Error(S2ErrorCode.FAILED_PRECONDITION, "Isolated vertex not isolated");
        return;
    }

    // Next we compute the winding number of R relative to I by counting signed
    // edge crossings of the input edges, and the winding number of R' related
    // to I' by counting signed edge crossings of the output edges.  (In order
    // to support DuplicateEdges::MERGE and SiblingEdges::CREATE, we also need
    // to take into account the number of input edges that snapped to each
    // output edge.)
    int winding_in = 0;
    S2CopyingEdgeCrosser crosser=new(iso_in, ref_in);
    for (int e = 0; e < builder_.NumInputEdges(); ++e)
    {
        var edge = builder_.InputEdge(e);
        winding_in += crosser.SignedEdgeOrVertexCrossing(edge.V0, edge.V1);
    }
    int winding_out = 0;
    crosser.Init(iso_out, ref_out);
    for (EdgeId e = 0; e < g.NumEdges; ++e)
    {
        var edge = g.GetEdge(e);
        winding_out += g.InputEdgeIds(e).Count *
                       crosser.SignedEdgeOrVertexCrossing(g.Vertex(edge.ShapeId),
                                                          g.Vertex(edge.EdgeId));
    }

    // Finally, check that s2builderutil::GetSnappedWindingDelta() computes the
    // difference between the two (using only local snapping information).
    int winding_delta = SnappedWindingDelta.GetSnappedWindingDelta(
        ref_in, ref_v, null, builder_, g, out error);
    Assert.True(error.IsOk());
    Assert.Equal(winding_delta, winding_out - winding_in);
    winding_tally_[winding_delta] += 1;
}

private readonly InputEdgeId ref_input_edge_id_;
private readonly InputEdgeId isolated_input_edge_id_;
private readonly S2Builder builder_;
private readonly WindingTally winding_tally_;
}

    [Fact]
    public void TestGetSnappedWindingDelta_RandomLoops() {
        // Count the number of tests for each winding number result and also the
        // number of tests where the isolated vertex was not isolated, to verify
        // that the test is working as intended.
        const int numIters = 1000;  // Passes with 10,000,000 iterations.
        int num_not_isolated = 0;
        WindingTally winding_tally = new();
        for (int iter = 0; iter < numIters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // For reproducability.
            System.Diagnostics.Debug.WriteLine("Iteration " + iter);

            // Choose a random snap radius up to the allowable maximum.
            S1Angle snap_radius = S2Testing.Random.RandDouble() *
                    SnapFunction.kMaxSnapRadius;
            S2Builder builder = new(new Options(new IdentitySnapFunction(snap_radius)));
            builder.StartLayer(new WindingNumberCheckingLayer(
                0 /*ref_input_edge_id*/, 1 /*isolated_input_edge_id*/,
                builder, winding_tally));

            // Choose a random reference vertex, and an isolated vertex that is as far
            // away as possible.  (The small amount of perturbation reduces the number
            // of calls to S2::ExpensiveSign() and is not necessary for correctness.)
            S2Point ref_ = S2Testing.RandomPoint();
            S2Point isolated = (-ref_ + 1e-12 * S2.Ortho(ref_)).Normalize();
            builder.AddEdge(ref_, ref_);            // Reference vertex edge.
            builder.AddEdge(isolated, isolated);  // Isolated vertex edge.

            // Now repeatedly build and add loops.  Loops consist of 1 or more random
            // vertices where approximately 2/3 are expected to be within snap_radius
            // of the reference vertex.  Some vertices are duplicates of previous
            // vertices.  Vertices are kept at least snap_radius away from the
            // isolated vertex to reduce the chance that edges will snap to it.
            // (This can still happen with long edges, or because the reference
            // vertex snapped to a new location far away from its original location.)
            List<S2Point> vertices_used = new(), loop = new();
            for (int num_loops = S2Testing.Random.Uniform(5) + 1; --num_loops >= 0;)
            {
                for (int num_vertices = S2Testing.Random.Uniform(9) + 1; --num_vertices >= 0;)
                {
                    if (vertices_used.Any() && S2Testing.Random.OneIn(4))
                    {
                        loop.Add(vertices_used[S2Testing.Random.Uniform(vertices_used.Count)]);
                    }
                    else if (S2Testing.Random.OneIn(3))
                    {
                        loop.Add(S2Testing.SamplePoint(
                            new S2Cap(ref_, S1Angle.FromRadians(S2.M_PI) - snap_radius)));
                    }
                    else
                    {
                        loop.Add(S2Testing.SamplePoint(new S2Cap(ref_, snap_radius)));
                    }
                }
                builder.AddShape(new S2LaxLoopShape(loop.ToArray()));
                loop.Clear();
            }

            if (!builder.Build(out _)) ++num_not_isolated;
        }
        // We expect that at most 20% of tests will result in an isolated vertex.
        Assert.True(num_not_isolated <= 0.2 * numIters);
        System.Diagnostics.Debug.WriteLine("Histogram of winding number deltas tested:");
        foreach (var entry in winding_tally) {
        System.Diagnostics.Debug.WriteLine($"{entry.Key} : {entry.Value}");
    }
}


}
