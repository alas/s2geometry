namespace S2Geometry;

using static S2Geometry.S2ShapeUtil;

public class S2CrossingEdgeQueryTests(ITestOutputHelper logger)
{
    private readonly ITestOutputHelper _logger = logger;

    // Test edges that lie in the plane of one of the S2 cube edges.  Such edges
    // may lie on the boundary between two cube faces, or pass through a cube
    // vertex, or follow a 45 diagonal across a cube face toward its center.
    //
    // This test is sufficient to demonstrate that padding the cell boundaries is
    // necessary for correctness.  (It fails if MutableS2ShapeIndex.kCellPadding
    // is set to zero.)
    [Fact]
    internal void Test_GetCrossingCandidates_PerturbedCubeEdges()
    {
        List<(S2Point, S2Point)> edges = [];
        for (int iter = 0; iter < 10; ++iter)
        {
            int face = S2Testing.Random.Uniform(6);
            double scale = Math.Pow(S2.DoubleError, S2Testing.Random.RandDouble());
            R2Point uv = new(2 * S2Testing.Random.Uniform(2) - 1, 2 * S2Testing.Random.Uniform(2) - 1);  // vertex
            S2Point a0 = S2.FaceUVtoXYZ(face, scale * uv);
            S2Point b0 = a0 - 2 * S2.GetNorm(face);
            // TODO(ericv): This test is currently slow because *every* crossing test
            // needs to invoke S2Pred.ExpensiveSign().
            GetPerturbedSubEdges(a0, b0, 30, edges);
            TestAllCrossings(edges);
        }
    }

    // Test edges that lie in the plane of one of the S2 cube face axes.  These
    // edges are special because one coordinate is zero, and they lie on the
    // boundaries between the immediate child cells of the cube face.
    [Fact]
    internal void Test_GetCrossingCandidates_PerturbedCubeFaceAxes()
    {
        List<(S2Point, S2Point)> edges = [];
        for (int iter = 0; iter < 5; ++iter)
        {
            int face = S2Testing.Random.Uniform(6);
            double scale = Math.Pow(S2.DoubleError, S2Testing.Random.RandDouble());
            S2Point axis = S2.GetUVWAxis(face, S2Testing.Random.Uniform(2));
            S2Point a0 = scale * axis + S2.GetNorm(face);
            S2Point b0 = scale * axis - S2.GetNorm(face);
            GetPerturbedSubEdges(a0, b0, 30, edges);
            TestAllCrossings(edges);
        }
    }

    [Fact]
    internal void Test_GetCrossingCandidates_CapEdgesNearCubeVertex()
    {
        // Test a random collection of edges near the S2 cube vertex where the
        // Hilbert curve starts and ends.
        List<(S2Point, S2Point)> edges = [];
        GetCapEdges(new S2Cap(new S2Point(-1, -1, 1).Normalize(), S1Angle.FromRadians(1e-3)),
                    S1Angle.FromRadians(1e-4), 1000, edges);
        TestAllCrossings(edges);
    }

    [Fact]
    internal void Test_GetCrossingCandidates_DegenerateEdgeOnCellVertexIsItsOwnCandidate()
    {
        for (int i = 0; i < 100; ++i)
        {
            List<(S2Point, S2Point)> edges = [];
            S2Cell cell = new(S2Testing.GetRandomCellId());
            edges.Add((cell.Vertex(0), cell.Vertex(0)));
            TestAllCrossings(edges);
        }
    }

    [Fact]
    internal void Test_GetCrossingCandidates_CollinearEdgesOnCellBoundaries()
    {
        int kNumEdgeIntervals = 8;  // 9*8/2 = 36 edges
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            S2Cell cell = new(S2Testing.GetRandomCellId(level));
            int j2 = S2Testing.Random.Uniform(4);
            S2Point p1 = cell.VertexRaw(j2);
            S2Point p2 = cell.VertexRaw(j2 + 1);
            S2Point delta = (p2 - p1) / kNumEdgeIntervals;
            List<(S2Point, S2Point)> edges = [];
            for (int i = 0; i <= kNumEdgeIntervals; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    edges.Add(((p1 + i * delta).Normalize(),
                                                   (p1 + j * delta).Normalize()));
                }
            }
            TestAllCrossings(edges);
        }
    }

    [Fact]
    internal void Test_GetCrossings_PolylineCrossings()
    {
        MutableS2ShapeIndex index =
        [
            // Three zig-zag lines near the equator.
            new S2Polyline.OwningShape(
                MakePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")),
            new S2Polyline.OwningShape(
                MakePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")),
            new S2Polyline.OwningShape(
                MakePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")),
        ];
        TestPolylineCrossings(index, MakePointOrDie("1:0"), MakePointOrDie("1:4"));
        TestPolylineCrossings(index, MakePointOrDie("5:5"), MakePointOrDie("6:6"));
    }

    [Fact]
    internal void Test_GetCrossings_ShapeIdsAreCorrect()
    {
        // This tests that when some index cells contain only one shape, the
        // intersecting edges are returned with the correct shape id.
        MutableS2ShapeIndex index =
        [
            new S2Polyline.OwningShape(
                new S2Polyline(S2Testing.MakeRegularPoints(
                    MakePointOrDie("0:0"), S1Angle.FromDegrees(5), 100))),
            new S2Polyline.OwningShape(
                new S2Polyline(S2Testing.MakeRegularPoints(
                    MakePointOrDie("0:20"), S1Angle.FromDegrees(5), 100))),
        ];
        TestPolylineCrossings(index, MakePointOrDie("1:-10"), MakePointOrDie("1:30"));
    }

    // Verifies that when VisitCells() is called with a specified root cell and a
    // query edge that barely intersects that cell, that at least one cell is
    // visited.  (At one point this was not always true, because when the query edge
    // is clipped to the index cell boundary without using any padding then the
    // result is sometimes empty, i.e., the query edge appears not to intersect the
    // specifed root cell.  The code now uses an appropriate amount of padding,
    // i.e. S2EdgeClipping.kFaceClipErrorUVCoord.)
    [Fact]
    internal void Test_VisitCells_QueryEdgeOnFaceBoundary()
    {
        int kIters = 100;
        for (int iter = 0; iter < kIters; ++iter)
        {
            _logger.WriteLine("Iteration " + iter);

            // Choose an edge AB such that B is nearly on the edge between two S2 cube
            // faces, and such that the result of clipping AB to the face that nominally
            // contains B (according to S2.GetFace) is empty when no padding is used.
            int a_face, b_face;
            S2Point a, b;
            R2Point a_uv, b_uv;
            do
            {
                a_face = S2Testing.Random.Uniform(6);
                a = S2.FaceUVtoXYZ(a_face, S2Testing.Random.UniformDouble(-1, 1),
                                    S2Testing.Random.UniformDouble(-1, 1)).Normalize();
                b_face = S2.GetUVWFace(a_face, 0, 1);  // Towards positive u-axis
                var uTmp = 1 - S2Testing.Random.Uniform(2) * 0.5 * S2.DoubleEpsilon;
                b = S2.FaceUVtoXYZ(b_face, uTmp, S2Testing.Random.UniformDouble(-1, 1)).Normalize();
            } while (S2.GetFace(b) != b_face ||
                     S2EdgeClipping.ClipToFace(a, b, b_face, out a_uv, out b_uv));

            // Verify that the clipping result is non-empty when a padding of
            // S2EdgeClipping.kFaceClipErrorUVCoord is used instead.
            Assert.True(S2EdgeClipping.ClipToPaddedFace(a, b, b_face, S2EdgeClipping.kFaceClipErrorUVCoord,
                                             out a_uv, out b_uv));

            // Create an S2ShapeIndex containing a single edge BC, where C is on the
            // same S2 cube face as B (which is different than the face containing A).
            S2Point c = S2.FaceUVtoXYZ(b_face, S2Testing.Random.UniformDouble(-1, 1),
                                        S2Testing.Random.UniformDouble(-1, 1)).Normalize();
            MutableS2ShapeIndex index =
            [
                new S2Polyline.OwningShape(
                    new S2Polyline([b, c])),
            ];

            // Check that the intersection between AB and BC is detected when the face
            // containing BC is specified as a root cell.  (Note that VisitCells()
            // returns false only if the CellVisitor returns false, and otherwise
            // returns true.)
            S2CrossingEdgeQuery query = new(index);
            S2PaddedCell root = new(S2CellId.FromFace(b_face), 0);
            Assert.False(query.VisitCells(a, b, root, x => false));
        }
    }

    private static S2Point PerturbAtDistance(S1Angle distance, S2Point a0, S2Point b0)
    {
        S2Point x = S2.GetPointOnLine(a0, b0, distance);
        if (S2Testing.Random.OneIn(2))
        {
            for (int i = 0; i < 3; ++i)
            {
                x = x.SetAxis(i, MathUtils.NextAfter(x[i], S2Testing.Random.OneIn(2) ? 1 : -1));
            }
            x = x.Normalize();
        }
        return x;
    }

    // Generate sub-edges of some given edge (a0,b0).  The length of the sub-edges
    // is distributed exponentially over a large range, and the endpoints may be
    // slightly perturbed to one side of (a0,b0) or the other.
    private static void GetPerturbedSubEdges(S2Point a0, S2Point b0, int count, List<(S2Point, S2Point)> edges)
    {
        edges.Clear();
        a0 = a0.Normalize();
        b0 = b0.Normalize();
        S1Angle length0 = new(a0, b0);
        for (int i = 0; i < count; ++i)
        {
            S1Angle length = length0 * Math.Pow(S2.DoubleError, S2Testing.Random.RandDouble());
            S1Angle offset = (length0 - length) * S2Testing.Random.RandDouble();
            edges.Add((PerturbAtDistance(offset, a0, b0),
                PerturbAtDistance(offset + length, a0, b0)));
        }
    }

    // Generate edges whose center is randomly chosen from the given S2Cap, and
    // whose length is randomly chosen up to "max_length".
    private static void GetCapEdges(S2Cap center_cap, S1Angle max_length, int count, List<(S2Point, S2Point)> edges)
    {
        edges.Clear();
        for (int i = 0; i < count; ++i)
        {
            S2Point center = S2Testing.SamplePoint(center_cap);
            S2Cap edge_cap = new(center, 0.5 * max_length);
            S2Point p1 = S2Testing.SamplePoint(edge_cap);
            // Compute p1 reflected through "center", and normalize for good measure.
            S2Point p2 = (2 * p1.DotProd(center) * center - p1).Normalize();
            edges.Add(new(p1, p2));
        }
    }

    // Project ShapeEdges to ShapeEdgeIds.  Useful because
    // ShapeEdge does not have operator==, but ShapeEdgeId does.
    private static List<Edge> GetShapeEdgeIds(List<ShapeEdge> shape_edges)
    {
        List<Edge> shape_edge_ids = [];
        foreach (var shape_edge in shape_edges)
        {
            shape_edge_ids.Add(shape_edge.Id);
        }
        return shape_edge_ids;
    }

    private void TestAllCrossings(List<(S2Point, S2Point)> edges)
    {
        var shape = new S2EdgeVectorShape();  // raw pointer since "shape" used below
        foreach (var edge in edges)
        {
            shape.Add(edge.Item1, edge.Item2);
        }
        // Force more subdivision than usual to make the test more challenging.
        MutableS2ShapeIndex.Options options = new()
        {
            MaxEdgesPerCell = 1
        };
        MutableS2ShapeIndex index = new(options);
        int shape_id = index.Add(shape);
        Assert.Equal(0, shape_id);
        // To check that candidates are being filtered reasonably, we count the
        // total number of candidates that the total number of edge pairs that
        // either intersect or are very close to intersecting.
        int num_candidates = 0, num_nearby_pairs = 0;
        int i2 = 0;
        foreach (var edge in edges)
        {
            _logger.WriteLine("Iteration " + i2++);
            S2Point a = edge.Item1;
            S2Point b = edge.Item2;
            S2CrossingEdgeQuery query = new(index);
            var candidates = query.GetCandidates(a, b, shape);

            // Verify that the second version of GetCandidates returns the same result.
            var edge_candidates = query.GetCandidates(a, b);
            Assert.Equal(candidates, edge_candidates);
            Assert.True(candidates.Count!=0);

            // Now check the actual candidates.
            Assert.True(candidates.IsSorted());
            Assert.Equal(0, candidates.Last().ShapeId);  // Implies all shape_ids are 0.
            Assert.True(candidates.First().EdgeId >= 0);
            Assert.True(candidates.Last().EdgeId < shape.NumEdges());
            num_candidates += candidates.Count;
            var missing_candidates = string.Empty;
            List<Edge> expected_crossings = [], expected_interior_crossings = [];
            for (int i = 0; i < shape.NumEdges(); ++i)
            {
                var edge2 = shape.GetEdge(i);
                S2Point c = edge2.V0;
                S2Point d = edge2.V1;
                int sign = S2.CrossingSign(a, b, c, d);
                if (sign >= 0)
                {
                    expected_crossings.Add(new(0, i));
                    if (sign > 0)
                    {
                        expected_interior_crossings.Add(new(0, i));
                    }
                    ++num_nearby_pairs;
                    if (candidates.BinarySearch(new(0, i)) < 0)
                    {
                        missing_candidates += " " + i;
                    }
                }
                else
                {
                    double kMaxDist = S2.kMaxDiag.GetValue(S2.kMaxCellLevel);
                    if (S2.GetDistance(a, c, d).Radians < kMaxDist ||
                        S2.GetDistance(b, c, d).Radians < kMaxDist ||
                        S2.GetDistance(c, a, b).Radians < kMaxDist ||
                        S2.GetDistance(d, a, b).Radians < kMaxDist)
                    {
                        ++num_nearby_pairs;
                    }
                }
            }
            Assert.True(missing_candidates.Length==0);

            // Test that GetCrossings() returns only the actual crossing edges.
            var actual_crossings =
                query.GetCrossingEdges(a, b, shape, CrossingType.ALL);
            Assert.Equal(expected_crossings, GetShapeEdgeIds(actual_crossings));

            // Verify that the second version of GetCrossings returns the same result.
            var actual_edge_crossings =
                query.GetCrossingEdges(a, b, CrossingType.ALL);
            Assert.Equal(expected_crossings, GetShapeEdgeIds(actual_edge_crossings));

            // Verify that CrossingType.INTERIOR returns only the interior crossings.
            var actual_interior_crossings =
                query.GetCrossingEdges(a, b, shape, CrossingType.INTERIOR);
            Assert.Equal(expected_interior_crossings,
                      GetShapeEdgeIds(actual_interior_crossings));
        }
        // There is nothing magical about this particular ratio; this check exists
        // to catch changes that dramatically increase the number of candidates.
        Assert.True(num_candidates <= 3 * num_nearby_pairs);
    }

    // This is the example from the header file, with a few extras.
    private static void TestPolylineCrossings(S2ShapeIndex index, S2Point a0, S2Point a1)
    {
        S2CrossingEdgeQuery query = new(index);
        var edges =
            query.GetCrossingEdges(a0, a1, CrossingType.ALL);
        if (edges.Count==0) return;

        foreach (var edge in edges)
        {
            Assert.True(S2.CrossingSign(a0, a1, edge.V0, edge.V1) >= 0);
        }
        // Also test that no edges are missing.
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = (S2Polyline.Shape)index.Shape(i)!;
            var polyline = shape.Polyline;
            for (int e = 0; e < polyline.NumVertices() - 1; ++e)
            {
                if (S2.CrossingSign(a0, a1, polyline.Vertex(e), polyline.Vertex(e + 1)) >= 0)
                {
                    Assert.Equal(1, edges.Count(edge => edge.Id == new Edge(i, e)));
                }
            }
        }
    }
}
