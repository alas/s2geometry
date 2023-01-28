namespace S2Geometry;

public class S2LoopTests
{
    private static readonly S2LatLng kRectError = S2LatLngRectBounder.MaxErrorForTests();

    private readonly ITestOutputHelper _logger;

    // The set of all loops declared below.
    private readonly List<S2Loop> all_loops = new();

    // Some standard loops to use in the tests (see descriptions below).
    private readonly S2Loop empty_;
    private readonly S2Loop full_;
    private readonly S2Loop north_hemi_;
    private readonly S2Loop north_hemi3_;
    private readonly S2Loop south_hemi_;
    private readonly S2Loop west_hemi_;
    private readonly S2Loop east_hemi_;
    private readonly S2Loop near_hemi_;
    private readonly S2Loop far_hemi_;
    private readonly S2Loop candy_cane_;
    private readonly S2Loop small_ne_cw_;
    private readonly S2Loop arctic_80_;
    private readonly S2Loop antarctic_80_;
    private readonly S2Loop line_triangle_;
    private readonly S2Loop skinny_chevron_;
    private readonly S2Loop loop_a_;
    private readonly S2Loop loop_b_;
    private readonly S2Loop a_intersect_b_;
    private readonly S2Loop a_union_b_;
    private readonly S2Loop a_minus_b_;
    private readonly S2Loop b_minus_a_;
    private readonly S2Loop loop_c_;
    private readonly S2Loop loop_d_;
    private readonly S2Loop loop_e_;
    private readonly S2Loop loop_f_;
    private readonly S2Loop loop_g_;
    private readonly S2Loop loop_h_;
    private readonly S2Loop loop_i_;
    private readonly S2Loop snapped_loop_a_;

    public S2LoopTests(ITestOutputHelper logger)
    {
        _logger = logger;
        // The empty loop.
        empty_ = AddLoop(S2Loop.kEmpty);

        // The full loop.
        full_ = AddLoop(S2Loop.kFull);

        // The northern hemisphere, defined using two pairs of antipodal points.
        north_hemi_ = AddLoop("0:-180, 0:-90, 0:0, 0:90");

        // The northern hemisphere, defined using three points 120 degrees apart.
        north_hemi3_ = AddLoop("0:-180, 0:-60, 0:60");

        // The southern hemisphere, defined using two pairs of antipodal points.
        south_hemi_ = AddLoop("0:90, 0:0, 0:-90, 0:-180");

        // The western hemisphere, defined using two pairs of antipodal points.
        west_hemi_ = AddLoop("0:-180, -90:0, 0:0, 90:0");

        // The eastern hemisphere, defined using two pairs of antipodal points.
        east_hemi_ = AddLoop("90:0, 0:0, -90:0, 0:-180");

        // The "near" hemisphere, defined using two pairs of antipodal points.
        near_hemi_ = AddLoop("0:-90, -90:0, 0:90, 90:0");

        // The "far" hemisphere, defined using two pairs of antipodal points.
        far_hemi_ = AddLoop("90:0, 0:90, -90:0, 0:-90");

        // A spiral stripe that slightly over-wraps the equator.
        candy_cane_ = AddLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70");

        // A small clockwise loop in the northern & eastern hemisperes.
        small_ne_cw_ = AddLoop("35:20, 45:20, 40:25");

        // Loop around the north pole at 80 degrees.
        arctic_80_ = AddLoop("80:-150, 80:-30, 80:90");

        // Loop around the south pole at 80 degrees.
        antarctic_80_ = AddLoop("-80:120, -80:0, -80:-120");

        // A completely degenerate triangle along the equator that Sign()
        // considers to be CCW.
        line_triangle_ = AddLoop("0:1, 0:2, 0:3");

        // A nearly-degenerate CCW chevron near the equator with very long sides
        // (about 80 degrees).  Its area is less than 1e-640, which is too small
        // to represent in double precision.
        skinny_chevron_ = AddLoop("0:0, -1e-320:80, 0:1e-320, 1e-320:80");

        // A diamond-shaped loop around the point 0:180.
        loop_a_ = AddLoop("0:178, -1:180, 0:-179, 1:-180");

        // Another diamond-shaped loop around the point 0:180.
        loop_b_ = AddLoop("0:179, -1:180, 0:-178, 1:-180");

        // The intersection of A and B.
        a_intersect_b_ = AddLoop("0:179, -1:180, 0:-179, 1:-180");

        // The union of A and B.
        a_union_b_ = AddLoop("0:178, -1:180, 0:-178, 1:-180");

        // A minus B (concave).
        a_minus_b_ = AddLoop("0:178, -1:180, 0:179, 1:-180");

        // B minus A (concave).
        b_minus_a_ = AddLoop("0:-179, -1:180, 0:-178, 1:-180");

        // A shape gotten from A by adding a triangle to one edge, and
        // subtracting a triangle from the opposite edge.
        loop_c_ = AddLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180");

        // A shape gotten from A by adding a triangle to one edge, and
        // adding another triangle to the opposite edge.
        loop_d_ = AddLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180");

        //   3------------2
        //   |            |               ^
        //   |  7-8  b-c  |               |
        //   |  | |  | |  |      Latitude |
        //   0--6-9--a-d--1               |
        //   |  | |       |               |
        //   |  f-e       |               +----------.
        //   |            |                 Longitude
        //   4------------5
        //
        // Important: It is not okay to skip over collinear vertices when
        // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
        // uses symbolic perturbations to ensure that no three vertices are
        // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
        // collinear).  In other words, it is unpredictable (modulo knowing the
        // details of the symbolic perturbations) whether 0123 contains 06123,
        // for example.
        //
        // Loop E:  0,6,9,a,d,1,2,3
        // Loop F:  0,4,5,1,d,a,9,6
        // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
        // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
        // Loop I:  7,6,f,e,9,8
        loop_e_ = AddLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30");
        loop_f_ = AddLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34");
        loop_g_ = AddLoop("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, " +
                        "10:41, 0:41, 0:44, 30:44, 30:30");
        loop_h_ = AddLoop("0:30, 0:34, -10:34, -10:36, 0:36, 0:39, " +
                        "10:39, 10:41, 0:41, 0:44, 30:44, 30:30");
        loop_i_ = AddLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36");

        // Like loop_a, but the vertices are at leaf cell centers.
        var snapped_loop_a_vertices = new S2Point[]{
                new S2CellId(MakePointOrDie("0:178")).ToPoint(),
                new S2CellId(MakePointOrDie("-1:180")).ToPoint(),
                new S2CellId(MakePointOrDie("0:-179")).ToPoint(),
                new S2CellId(MakePointOrDie("1:-180")).ToPoint()};
        snapped_loop_a_ = AddLoop(new S2Loop(snapped_loop_a_vertices));
    }

    [Fact]
    internal void Test_S2LoopTestBase_GetRectBound()
    {
        Assert.True(empty_.GetRectBound().IsEmpty());
        Assert.True(full_.GetRectBound().IsFull());
        Assert.True(candy_cane_.GetRectBound().Lng.IsFull());
        Assert.True(candy_cane_.GetRectBound().LatLo().GetDegrees() < -20);
        Assert.True(candy_cane_.GetRectBound().LatHi().GetDegrees() > 10);
        Assert.True(small_ne_cw_.GetRectBound().IsFull());
        Assert.True(arctic_80_.GetRectBound().ApproxEquals(
            new S2LatLngRect(S2LatLng.FromDegrees(80, -180),
                         S2LatLng.FromDegrees(90, 180)), kRectError));
        Assert.True(antarctic_80_.GetRectBound().ApproxEquals(
            new S2LatLngRect(S2LatLng.FromDegrees(-90, -180),
                         S2LatLng.FromDegrees(-80, 180)), kRectError));

        // Create a loop that contains the complement of the "arctic_80" loop.
        var arctic_80_inv = (S2Loop)arctic_80_.CustomClone();
        arctic_80_inv.Invert();
        // The highest latitude of each edge is attained at its midpoint.
        S2Point mid = 0.5 * (arctic_80_inv.Vertex(0) + arctic_80_inv.Vertex(1));
        Assert2.Near(arctic_80_inv.GetRectBound().LatHi().Radians,
                    new S2LatLng(mid).LatRadians, kRectError.LatRadians);

        Assert.True(south_hemi_.GetRectBound().Lng.IsFull());
        Assert.True(south_hemi_.GetRectBound().Lat.ApproxEquals(
            new R1Interval(-S2.M_PI_2, 0), kRectError.LatRadians));
    }

    [Fact]
    internal void Test_S2LoopTestBase_AreaConsistentWithCurvature()
    {
        // Check that the area computed using GetArea() is consistent with the
        // curvature of the loop computed using GetTurnAngle().  According to
        // the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
        // minus its curvature.
        foreach (var loop in all_loops)
        {
            var area = loop.Area();
            var gauss_area = S2.M_2_PI - loop.Curvature();
            // The error bound is sufficient for current tests but not guaranteed.
            _logger.WriteLine("loop: " + loop.ToDebugString());
            Assert.True(Math.Abs(area - gauss_area) <= 1e-14);
        }
    }

    [Fact]
    internal void Test_S2LoopTestBase_GetAreaConsistentWithSign()
    {
        // Test that GetArea() returns an area near 0 for degenerate loops that
        // contain almost no points, and an area near 4*Pi for degenerate loops that
        // contain almost all points.

        const int kMaxVertices = 6;
        for (int i = 0; i < 50; ++i)
        {
            int num_vertices = 3 + S2Testing.Random.Uniform(kMaxVertices - 3 + 1);
            // Repeatedly choose N vertices that are exactly on the equator until we
            // find some that form a valid loop.
            S2Loop loop;
            //loop.set_s2debug_override(S2Debug.DISABLE);
            do
            {
                var vertices = new List<S2Point>();
                for (int i2 = 0; i2 < num_vertices; ++i2)
                {
                    // We limit longitude to the range [0, 90] to ensure that the loop is
                    // degenerate (as opposed to following the entire equator).
                    vertices.Add(S2LatLng.FromRadians(0,
                        S2Testing.Random.RandDouble() * S2.M_PI_2).ToPoint());
                }
                loop = new S2Loop(vertices, S2Debug.DISABLE);
            } while (!loop.IsValid());
            bool ccw = loop.IsNormalized();
            _logger.WriteLine("loop: " + loop.ToDebugString());
            Assert2.Near(ccw ? 0 : S2.M_4_PI, loop.Area(), S2.DoubleError);
            Assert.Equal(!ccw, loop.Contains(new S2Point(0, 0, 1)));
        }
    }

    [Fact]
    internal void Test_S2LoopTestBase_GetAreaAccuracy()
    {
        // TODO(b/200091211): Test that GetArea() has an accuracy significantly better
        // than 1e-15 on loops whose area is small.
    }

    [Fact]
    internal void Test_S2LoopTestBase_GetAreaAndCentroid()
    {
        Assert.Equal(0.0, empty_.Area());
        Assert.Equal(S2.M_4_PI, full_.Area());
        Assert.Equal(S2Point.Empty, empty_.Centroid());
        Assert.Equal(S2Point.Empty, full_.Centroid());

        Assert2.Near(north_hemi_.Area(), S2.M_2_PI);
        Assert2.Near(east_hemi_.Area(), S2.M_2_PI, 1e-15);

        // Construct spherical caps of random height, and approximate their boundary
        // with closely spaces vertices.  Then check that the area and centroid are
        // correct.

        for (int i = 0; i < 50; ++i)
        {
            // Choose a coordinate frame for the spherical cap.
            S2Testing.GetRandomFrame(out var x, out var y, out var z);

            // Given two points at latitude phi and whose longitudes differ by dtheta,
            // the geodesic between the two points has a maximum latitude of
            // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
            // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
            //
            // We want to position the vertices close enough together so that their
            // maximum distance from the boundary of the spherical cap is kMaxDist.
            // Thus we want Math.Abs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
            const double kMaxDist = 1e-6;
            double height = 2 * S2Testing.Random.RandDouble();
            double phi = Math.Asin(1 - height);
            double max_dtheta = 2 * Math.Acos(Math.Tan(Math.Abs(phi)) / Math.Tan(Math.Abs(phi) + kMaxDist));
            max_dtheta = Math.Min(Math.PI, max_dtheta);  // At least 3 vertices.

            var vertices = new List<S2Point>();
            for (double theta = 0; theta < S2.M_2_PI;
                 theta += S2Testing.Random.RandDouble() * max_dtheta)
            {
                vertices.Add(Math.Cos(theta) * Math.Cos(phi) * x +
                                   Math.Sin(theta) * Math.Cos(phi) * y +
                                   Math.Sin(phi) * z);
            }
            S2Loop loop = new(vertices);
            double area = loop.Area();
            S2Point centroid = loop.Centroid();
            double expected_area = S2.M_2_PI * height;
            Assert.True(Math.Abs(area - expected_area) <= S2.M_2_PI * kMaxDist);
            S2Point expected_centroid = expected_area * (1 - 0.5 * height) * z;
            Assert.True((centroid - expected_centroid).Norm() <= 2 * kMaxDist);
        }
    }

    [Fact]
    internal void Test_S2LoopTestBase_GetCurvature()
    {
        Assert.Equal(S2.M_2_PI, empty_.Curvature());
        Assert.Equal(-S2.M_2_PI, full_.Curvature());

        Assert2.Near(0, north_hemi3_.Curvature(), 1e-15);
        CheckCurvatureInvariants(north_hemi3_);

        Assert2.Near(0, west_hemi_.Curvature(), 1e-15);
        CheckCurvatureInvariants(west_hemi_);

        // We don't have an easy way to estimate the curvature of this loop, but
        // we can still check that the expected invariants hold.
        CheckCurvatureInvariants(candy_cane_);

        Assert2.Near(S2.M_2_PI, line_triangle_.Curvature());
        CheckCurvatureInvariants(line_triangle_);

        Assert2.Near(S2.M_2_PI, skinny_chevron_.Curvature());
        CheckCurvatureInvariants(skinny_chevron_);

        // Build a narrow spiral loop starting at the north pole.  This is designed
        // to test that the error in GetCurvature is linear in the number of
        // vertices even when the partial sum of the curvatures gets very large.
        // The spiral consists of two "arms" defining opposite sides of the loop.
        int kArmPoints = 10000;    // Number of vertices in each "arm"
        double kArmRadius = 0.01;  // Radius of spiral.
        var vertices = new S2Point[2 * kArmPoints];
        vertices[kArmPoints] = new S2Point(0, 0, 1);
        for (int i = 0; i < kArmPoints; ++i)
        {
            double angle = (S2.M_2_PI / 3) * i;
            double x = Math.Cos(angle);
            double y = Math.Sin(angle);
            double r1 = i * kArmRadius / kArmPoints;
            double r2 = (i + 1.5) * kArmRadius / kArmPoints;
            vertices[kArmPoints - i - 1] = new S2Point(r1 * x, r1 * y, 1).Normalize();
            vertices[kArmPoints + i] = new S2Point(r2 * x, r2 * y, 1).Normalize();
        }
        // This is a pathological loop that contains many long parallel edges, and
        // takes tens of seconds to validate in debug mode.
        S2Loop spiral = new(vertices, S2Debug.DISABLE);

        // Check that GetCurvature() is consistent with GetArea() to within the
        // error bound of the former.  We actually use a tiny fraction of the
        // worst-case error bound, since the worst case only happens when all the
        // roundoff errors happen in the same direction and this test is not
        // designed to achieve that.  The error in GetArea() can be ignored for the
        // purposes of this test since it is generally much smaller.
        Assert2.Near(S2.M_2_PI - spiral.Area(), spiral.Curvature(),
                    0.01 * spiral.CurvatureMaxError());
    }

    [Fact]
    internal void Test_S2LoopTestBase_NormalizedCompatibleWithContains()
    {
        CheckNormalizeAndContains(line_triangle_);
        CheckNormalizeAndContains(skinny_chevron_);
    }

    [Fact]
    internal void Test_S2LoopTestBase_Contains()
    {
        // Check the full and empty loops have the correct containment relationship
        // with the special "vertex" that defines them.
        Assert.False(empty_.Contains(S2Loop.kEmpty[0]));
        Assert.True(full_.Contains(S2Loop.kFull[0]));

        Assert.True(candy_cane_.Contains(S2LatLng.FromDegrees(5, 71).ToPoint()));

        // Create copies of these loops so that we can change the vertex order.
        var north_copy = (S2Loop)north_hemi_.CustomClone();
        var south_copy = (S2Loop)south_hemi_.CustomClone();
        var west_copy = (S2Loop)west_hemi_.CustomClone();
        var east_copy = (S2Loop)east_hemi_.CustomClone();
        for (int i = 0; i < 4; ++i)
        {
            Assert.True(north_copy.Contains(new S2Point(0, 0, 1)));
            Assert.False(north_copy.Contains(new S2Point(0, 0, -1)));
            Assert.False(south_copy.Contains(new S2Point(0, 0, 1)));
            Assert.True(south_copy.Contains(new S2Point(0, 0, -1)));
            Assert.False(west_copy.Contains(new S2Point(0, 1, 0)));
            Assert.True(west_copy.Contains(new S2Point(0, -1, 0)));
            Assert.True(east_copy.Contains(new S2Point(0, 1, 0)));
            Assert.False(east_copy.Contains(new S2Point(0, -1, 0)));
            Rotate(ref north_copy);
            Rotate(ref south_copy);
            Rotate(ref east_copy);
            Rotate(ref west_copy);
        }

        // This code checks each cell vertex is contained by exactly one of
        // the adjacent cells.
        for (int level = 0; level < 3; ++level)
        {
            var loops = new List<S2Loop>();
            var loop_vertices = new List<S2Point>();
            var points = new List<S2Point>();
            for (S2CellId id = S2CellId.Begin(level);
                 id != S2CellId.End(level); id = id.Next())
            {
                S2Cell cell = new(id);
                points.Add(cell.Center());
                for (int k = 0; k < 4; ++k)
                {
                    loop_vertices.Add(cell.Vertex(k));
                    points.Add(cell.Vertex(k));
                }
                loops.Add(new S2Loop(loop_vertices));
                loop_vertices.Clear();
            }
            foreach (var point in points)
            {
                int count = 0;
                foreach (var loop in loops)
                {
                    if (loop.Contains(point)) ++count;
                }
                Assert.Equal(1, count);
            }
        }
    }

    [Fact]
    internal void Test_S2Loop_ContainsMatchesCrossingSign()
    {
        // This test demonstrates a former incompatibility between CrossingSign()
        // and Contains(S2Point).  Itructs an S2Cell-based loop L and
        // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
        // previously, Contains() returned false for both endpoints of E.
        //
        // The reason for the bug was that the loop bound was sometimes too tight.
        // The Contains() code for a0 bailed out early because a0 was found not to
        // be inside the bound of L.

        // Start with a cell that ends up producing the problem.
        S2CellId cell_id = new S2CellId(new S2Point(1, 1, 1)).Parent(21);

        var children = new S2Cell[4];
        new S2Cell(cell_id).Subdivide(children);

        var points = new S2Point[4];
        for (int i = 0; i < 4; ++i)
        {
            // Note extra normalization. GetCenter() is already normalized.
            // The test results will no longer be inconsistent if the extra
            // Normalize() is removed.
            points[i] = children[i].Center().Normalize();
        }

        S2Loop loop = new(points);

        // Get a vertex from a grandchild cell.
        // +---------------+---------------+
        // |               |               |
        // |    points[3]  |   points[2]   |
        // |       v       |       v       |
        // |       +-------+------ +       |
        // |       |       |       |       |
        // |       |       |       |       |
        // |       |       |       |       |
        // +-------+-------+-------+-------+
        // |       |       |       |       |
        // |       |    <----------------------- grandchild_cell
        // |       |       |       |       |
        // |       +-------+------ +       |
        // |       ^       |       ^       | <-- cell
        // | points[0]/a0  |     points[1] |
        // |               |               |
        // +---------------+---------------+
        S2Cell grandchild_cell = new(cell_id.Child(0).Child(2));
        S2Point a0 = grandchild_cell.Vertex(0);

        // If this doesn't hold, the rest of the test is pointless.
        // This test depends on rounding errors that should make
        // a0 slightly different from points[0]
        // setprecision(20)
        Assert.NotEqual(points[0], a0);

        // The edge from a0 to the origin crosses one boundary.
        Assert.Equal(-1, S2.CrossingSign(a0, S2.Origin, loop.Vertex(0), loop.Vertex(1)));
        Assert.Equal(1, S2.CrossingSign(a0, S2.Origin, loop.Vertex(1), loop.Vertex(2)));
        Assert.Equal(-1, S2.CrossingSign(a0, S2.Origin, loop.Vertex(2), loop.Vertex(3)));
        Assert.Equal(-1, S2.CrossingSign(a0, S2.Origin, loop.Vertex(3), loop.Vertex(4)));

        // Contains should return false for the origin, and true for a0.
        Assert.False(loop.Contains(S2.Origin));
        Assert.True(loop.Contains(a0));

        // Since a0 is inside the loop, it should be inside the bound.
        S2LatLngRect bound = loop.GetRectBound();
        Assert.True(bound.Contains(a0));
    }

    [Fact]
    internal void Test_S2LoopTestBase_LoopRelations()
    {
        // Check full and empty relationships with normal loops and each other.
        TestRelation(full_, full_, RelationFlags.CONTAINS | RelationFlags.CONTAINED | RelationFlags.COVERS, true);
        TestRelation(full_, north_hemi_, RelationFlags.CONTAINS | RelationFlags.COVERS, false);
        TestRelation(full_, empty_, RelationFlags.CONTAINS | RelationFlags.DISJOINT | RelationFlags.COVERS, false);
        TestRelation(north_hemi_, full_, RelationFlags.CONTAINED | RelationFlags.COVERS, false);
        TestRelation(north_hemi_, empty_, RelationFlags.CONTAINS | RelationFlags.DISJOINT, false);
        TestRelation(empty_, full_, RelationFlags.CONTAINED | RelationFlags.DISJOINT | RelationFlags.COVERS, false);
        TestRelation(empty_, north_hemi_, RelationFlags.CONTAINED | RelationFlags.DISJOINT, false);
        TestRelation(empty_, empty_, RelationFlags.CONTAINS | RelationFlags.CONTAINED | RelationFlags.DISJOINT, false);

        TestRelation(north_hemi_, north_hemi_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(north_hemi_, south_hemi_, RelationFlags.DISJOINT | RelationFlags.COVERS, true);
        TestRelation(north_hemi_, east_hemi_, 0, false);
        TestRelation(north_hemi_, arctic_80_, RelationFlags.CONTAINS, false);
        TestRelation(north_hemi_, antarctic_80_, RelationFlags.DISJOINT, false);
        TestRelation(north_hemi_, candy_cane_, 0, false);

        // We can't compare north_hemi3 vs. north_hemi or south_hemi because the
        // result depends on the "simulation of simplicity" implementation details.
        TestRelation(north_hemi3_, north_hemi3_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(north_hemi3_, east_hemi_, 0, false);
        TestRelation(north_hemi3_, arctic_80_, RelationFlags.CONTAINS, false);
        TestRelation(north_hemi3_, antarctic_80_, RelationFlags.DISJOINT, false);
        TestRelation(north_hemi3_, candy_cane_, 0, false);

        TestRelation(south_hemi_, north_hemi_, RelationFlags.DISJOINT | RelationFlags.COVERS, true);
        TestRelation(south_hemi_, south_hemi_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(south_hemi_, far_hemi_, 0, false);
        TestRelation(south_hemi_, arctic_80_, RelationFlags.DISJOINT, false);
        TestRelation(south_hemi_, antarctic_80_, RelationFlags.CONTAINS, false);
        TestRelation(south_hemi_, candy_cane_, 0, false);

        TestRelation(candy_cane_, north_hemi_, 0, false);
        TestRelation(candy_cane_, south_hemi_, 0, false);
        TestRelation(candy_cane_, arctic_80_, RelationFlags.DISJOINT, false);
        TestRelation(candy_cane_, antarctic_80_, RelationFlags.DISJOINT, false);
        TestRelation(candy_cane_, candy_cane_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);

        TestRelation(near_hemi_, west_hemi_, 0, false);

        TestRelation(small_ne_cw_, south_hemi_, RelationFlags.CONTAINS, false);
        TestRelation(small_ne_cw_, west_hemi_, RelationFlags.CONTAINS, false);

        TestRelation(small_ne_cw_, north_hemi_, RelationFlags.COVERS, false);
        TestRelation(small_ne_cw_, east_hemi_, RelationFlags.COVERS, false);

        TestRelation(loop_a_, loop_a_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(loop_a_, loop_b_, 0, false);
        TestRelation(loop_a_, a_intersect_b_, RelationFlags.CONTAINS, true);
        TestRelation(loop_a_, a_union_b_, RelationFlags.CONTAINED, true);
        TestRelation(loop_a_, a_minus_b_, RelationFlags.CONTAINS, true);
        TestRelation(loop_a_, b_minus_a_, RelationFlags.DISJOINT, true);

        TestRelation(loop_b_, loop_a_, 0, false);
        TestRelation(loop_b_, loop_b_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(loop_b_, a_intersect_b_, RelationFlags.CONTAINS, true);
        TestRelation(loop_b_, a_union_b_, RelationFlags.CONTAINED, true);
        TestRelation(loop_b_, a_minus_b_, RelationFlags.DISJOINT, true);
        TestRelation(loop_b_, b_minus_a_, RelationFlags.CONTAINS, true);

        TestRelation(a_intersect_b_, loop_a_, RelationFlags.CONTAINED, true);
        TestRelation(a_intersect_b_, loop_b_, RelationFlags.CONTAINED, true);
        TestRelation(a_intersect_b_, a_intersect_b_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(a_intersect_b_, a_union_b_, RelationFlags.CONTAINED, false);
        TestRelation(a_intersect_b_, a_minus_b_, RelationFlags.DISJOINT, true);
        TestRelation(a_intersect_b_, b_minus_a_, RelationFlags.DISJOINT, true);

        TestRelation(a_union_b_, loop_a_, RelationFlags.CONTAINS, true);
        TestRelation(a_union_b_, loop_b_, RelationFlags.CONTAINS, true);
        TestRelation(a_union_b_, a_intersect_b_, RelationFlags.CONTAINS, false);
        TestRelation(a_union_b_, a_union_b_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(a_union_b_, a_minus_b_, RelationFlags.CONTAINS, true);
        TestRelation(a_union_b_, b_minus_a_, RelationFlags.CONTAINS, true);

        TestRelation(a_minus_b_, loop_a_, RelationFlags.CONTAINED, true);
        TestRelation(a_minus_b_, loop_b_, RelationFlags.DISJOINT, true);
        TestRelation(a_minus_b_, a_intersect_b_, RelationFlags.DISJOINT, true);
        TestRelation(a_minus_b_, a_union_b_, RelationFlags.CONTAINED, true);
        TestRelation(a_minus_b_, a_minus_b_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
        TestRelation(a_minus_b_, b_minus_a_, RelationFlags.DISJOINT, false);

        TestRelation(b_minus_a_, loop_a_, RelationFlags.DISJOINT, true);
        TestRelation(b_minus_a_, loop_b_, RelationFlags.CONTAINED, true);
        TestRelation(b_minus_a_, a_intersect_b_, RelationFlags.DISJOINT, true);
        TestRelation(b_minus_a_, a_union_b_, RelationFlags.CONTAINED, true);
        TestRelation(b_minus_a_, a_minus_b_, RelationFlags.DISJOINT, false);
        TestRelation(b_minus_a_, b_minus_a_, RelationFlags.CONTAINS | RelationFlags.CONTAINED, true);
    }

    // Make sure the relations are correct if the loop crossing happens on
    // two ends of a shared boundary segment.
    [Fact]
    internal void Test_S2LoopTestBase_LoopRelationsWhenSameExceptPiecesStickingOutAndIn()
    {
        TestRelation(loop_a_, loop_c_, 0, true);
        TestRelation(loop_c_, loop_a_, 0, true);
        TestRelation(loop_a_, loop_d_, RelationFlags.CONTAINED, true);
        TestRelation(loop_d_, loop_a_, RelationFlags.CONTAINS, true);
        TestRelation(loop_e_, loop_f_, RelationFlags.DISJOINT, true);
        TestRelation(loop_e_, loop_g_, RelationFlags.CONTAINS, true);
        TestRelation(loop_e_, loop_h_, 0, true);
        TestRelation(loop_e_, loop_i_, 0, false);
        TestRelation(loop_f_, loop_g_, RelationFlags.DISJOINT, true);
        TestRelation(loop_f_, loop_h_, 0, true);
        TestRelation(loop_f_, loop_i_, 0, false);
        TestRelation(loop_g_, loop_h_, RelationFlags.CONTAINED, true);
        TestRelation(loop_h_, loop_g_, RelationFlags.CONTAINS, true);
        TestRelation(loop_g_, loop_i_, RelationFlags.DISJOINT, true);
        TestRelation(loop_h_, loop_i_, RelationFlags.CONTAINS, true);
    }

    [Fact]
    internal void Test_S2Loop_LoopRelations2()
    {
        // Construct polygons consisting of a sequence of adjacent cell ids
        // at some fixed level.  Comparing two polygons at the same level
        // ensures that there are no T-vertices.
        for (int iter = 0; iter < 1000; ++iter)
        {
            S2CellId begin = new(S2Testing.Random.Rand64() | 1);
            if (!begin.IsValid()) continue;
            begin = begin.Parent(S2Testing.Random.Uniform(S2.kMaxCellLevel));
            S2CellId a_begin = begin.Advance(S2Testing.Random.Skewed(6));
            S2CellId a_end = a_begin.Advance(S2Testing.Random.Skewed(6) + 1);
            S2CellId b_begin = begin.Advance(S2Testing.Random.Skewed(6));
            S2CellId b_end = b_begin.Advance(S2Testing.Random.Skewed(6) + 1);
            if (!a_end.IsValid() || !b_end.IsValid()) continue;

            var a = MakeCellLoop(a_begin, a_end);
            var b = MakeCellLoop(b_begin, b_end);
            if (a is not null && b is not null)
            {
                bool contained = (a_begin <= b_begin && b_end <= a_end);
                bool intersects = (a_begin < b_end && b_begin < a_end);
                Assert.Equal(a.Contains(b), contained);
                Assert.Equal(a.Intersects(b), intersects);
            }
            else
            {
                Assert.True(false); //MakeCellLoop failed to create a loop.
            }
        }
    }

    [Fact]
    internal void Test_S2Loop_BoundsForLoopContainment()
    {
        // To reliably test whether one loop contains another, the bounds of the
        // outer loop are expanded slightly.  This testructs examples where
        // this expansion is necessary and verifies that it is sufficient.

        for (int iter = 0; iter < 1000; ++iter)
        {
            // We construct a triangle ABC such that A,B,C are nearly colinear, B is
            // the point of maximum latitude, and the edge AC passes very slightly
            // below B (i.e., ABC is CCW).
            var b = (S2Testing.RandomPoint() + new S2Point(0, 0, 1)).Normalize();
            var v = b.CrossProd(new S2Point(0, 0, 1)).Normalize();
            var a = S2.Interpolate(-v, b, S2Testing.Random.RandDouble());
            var c = S2.Interpolate(b, v, S2Testing.Random.RandDouble());
            if (S2Pred.Sign(a, b, c) < 0)
            {
                --iter; continue;
            }
            // Now construct another point D directly below B, and create two loops
            // ABCD and ACD.
            S2Point d = new S2Point(b.X, b.Y, 0).Normalize();
            var vertices = new S2Point[] { c, d, a, b };  // Reordered for convenience
            S2Loop outer = new(vertices.Take(4));
            S2Loop inner = new(vertices.Take(3));
            // Now because the bounds calculation is less accurate when the maximum is
            // attained along an edge (rather than at a vertex), sometimes the inner
            // loop will have a *larger* bounding box than the outer loop.  We look
            // only for those cases.
            if (outer.GetRectBound().Contains(inner.GetRectBound()))
            {
                --iter; continue;
            }
            Assert.True(outer.Contains(inner));
        }
    }

    [Fact]
    internal void Test_S2Loop_BoundaryNear()
    {
        S1Angle degree = S1Angle.FromDegrees(1);

        TestNear("0:0, 0:10, 5:5",
                 "0:0.1, -0.1:9.9, 5:5.2",
                 0.5 * degree, true);
        TestNear("0:0, 0:3, 0:7, 0:10, 3:7, 5:5",
                 "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1",
                 S1Angle.FromRadians(1e-3), true);

        // All vertices close to some edge, but not equivalent.
        TestNear("0:0, 0:2, 2:2, 2:0",
                 "0:0, 1.9999:1, 0:2, 2:2, 2:0",
                 0.5 * degree, false);

        // Two triangles that backtrack a bit on different edges.  A simple
        // greedy matching algorithm would fail on this example.
        string t1 = "0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, " +
               "2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2";
        string t2 = "0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, " +
               "0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4";
        TestNear(t1, t2, 1.5 * degree, true);
        TestNear(t1, t2, 0.5 * degree, false);
    }

    [Fact]
    internal void Test_S2Loop_EncodeDecode()
    {
        var l = MakeLoopOrDie("30:20, 40:20, 39:43, 33:35");
        l.Depth = 3;
        TestEncodeDecode(l);

        S2Loop empty = S2Loop.kEmpty;
        TestEncodeDecode(empty);
        S2Loop full = S2Loop.kFull;
        TestEncodeDecode(full);

        S2Loop uninitialized = new(Array.Empty<S2Point>());
        TestEncodeDecode(uninitialized);
    }

    [Fact]
    internal void Test_S2Loop_Moveable()
    {
        // We'll need a couple identical copies of a reference loop to compare.
        var loop_factory = () => {
            var loop = MakeLoopOrDie("30:20, 40:20, 39:43, 33:35");
            loop!.Depth = 3;
            return loop;
        };

        // Check for move-constructability.
        {
            var loop = loop_factory();
            S2Loop a = loop;
            CheckIdentical(a, loop_factory());
        }

        // Check for move-assignability.
        {
            var loop = loop_factory();
            S2Loop a;
            a = loop;
            CheckIdentical(a, loop_factory());
        }
    }

    [Fact]
    internal void Test_S2Loop_EmptyFullLossyConversions()
    {
        // Verify that the empty and full loops can be encoded lossily.
        S2Loop empty = S2Loop.kEmpty;
        TestEmptyFullConversions(empty);

        S2Loop full = S2Loop.kFull;
        TestEmptyFullConversions(full);
    }

    [Fact]
    internal void Test_S2Loop_EncodeDecodeWithinScope()
    {
        S2Loop l = MakeLoopOrDie("30:20, 40:20, 39:43, 33:35");
        l.Depth = 3;
        Encoder encoder = new();
        l.Encode(encoder);
        var decoder1 = encoder.GetDecoder();

        // Initialize the loop using DecodeWithinScope and check that it is the
        // same as the original loop.
        var (success1, loop1) = S2Loop.Decode(decoder1);
        Assert.True(success1);
        Assert.True(l.BoundaryEquals(loop1));
        Assert.Equal(l.Depth, loop1.Depth);
        Assert.Equal(l.GetRectBound(), loop1.GetRectBound());

        // Initialize the same loop using Init with a vector of vertices, and
        // check that it doesn't deallocate the original memory.
        S2Point[] vertices = { loop1.Vertex(0), loop1.Vertex(2), loop1.Vertex(3) };
        loop1 = new S2Loop(vertices);
        var decoder2 = encoder.GetDecoder();
        var (success2, loop2) = S2Loop.Decode(decoder2);
        Assert.True(success2);
        Assert.True(l.BoundaryEquals(loop2!));
        Assert.Equal(l.Vertex(1), loop2.Vertex(1));
        Assert.NotEqual(loop1.Vertex(1), loop2.Vertex(1));

        // Initialize loop2 using Decode with a decoder on different data.
        // Check that the original memory is not deallocated or overwritten.
        var l2 = MakeLoopOrDie("30:40, 40:75, 39:43, 80:35");
        l2.Depth = 2;
        Encoder encoder2 = new();
        l2.Encode(encoder2);
        var decoder3 = encoder2.GetDecoder();
        var (success3, loop3) = S2Loop.Decode(decoder3);
        Assert.True(success3);
        var decoder4 = encoder.GetDecoder();
        var (success4, loop4) = S2Loop.Decode(decoder4);
        Assert.True(success4);
        Assert.True(l.BoundaryEquals(loop4));
        Assert.Equal(l.Vertex(1), loop4.Vertex(1));
        Assert.NotEqual(loop4.Vertex(1), loop3.Vertex(1));
    }

    [Fact]
    internal void Test_S2LoopTestBase_FourVertexCompressedLoopRequires36Bytes()
    {
        Encoder encoder = new();
        TestEncodeCompressed(snapped_loop_a_, S2.kMaxCellLevel, encoder);

        // 1 byte for num_vertices
        // 1 byte for origin_inside and boolean indicating we did not
        //   encode the bound
        // 1 byte for depth
        // Vertices:
        // 1 byte for faces
        // 8 bytes for each vertex.
        // 1 byte indicating that there is no unsnapped vertex.
        Assert.Equal(37, encoder.Length());
    }

    [Fact]
    internal void Test_S2LoopTestBase_CompressedEncodedLoopDecodesApproxEqual()
    {
        var loop = (S2Loop)snapped_loop_a_.CustomClone();
        loop.Depth = 3;

        Encoder encoder = new();
        TestEncodeCompressed(loop, S2.kMaxCellLevel, encoder);
        TestDecodeCompressed(encoder, S2.kMaxCellLevel, out var decoded_loop);
        CheckIdentical(loop, decoded_loop);
    }

    // This test checks that S2Loops created directly from S2Cells behave
    // identically to S2Loops created from the vertices of those cells; this
    // previously was not the case, because S2Cells calculate their bounding
    // rectangles slightly differently, and S2Loops created from them just copied
    // the S2Cell bounds.
    [Fact]
    internal void Test_S2Loop_S2CellConstructorAndContains()
    {
        S2Cell cell = new(new S2CellId(S2LatLng.FromE6(40565459, -74645276)));
        S2Loop cell_as_loop = new(cell);

        var vertices = new List<S2Point>();
        for (int i = 0; i < cell_as_loop.NumVertices; ++i)
        {
            vertices.Add(cell_as_loop.Vertex(i));
        }
        S2Loop loop_copy = new(vertices);
        Assert.True(loop_copy.Contains(cell_as_loop));
        Assert.True(cell_as_loop.Contains(loop_copy));

        // Demonstrates the reason for this test; the cell bounds are more
        // conservative than the resulting loop bounds.
        Assert.False(loop_copy.GetRectBound().Contains(cell.GetRectBound()));
    }

    [Fact]
    internal void Test_S2Loop_IsValidDetectsInvalidLoops()
    {
        // Not enough vertices.  Note that all single-vertex loops are valid; they
        // are interpreted as being either empty or full.
        CheckLoopIsInvalid("", "at least 3 vertices");
        CheckLoopIsInvalid("20:20, 21:21", "at least 3 vertices");

        // There is a degenerate edge
        CheckLoopIsInvalid("20:20, 20:20, 20:21", "degenerate");
        CheckLoopIsInvalid("20:20, 20:21, 20:20", "degenerate");

        // There is a duplicate vertex
        CheckLoopIsInvalid("20:20, 21:21, 21:20, 20:20, 20:21", "duplicate vertex");

        // Some edges cross
        CheckLoopIsInvalid("20:20, 21:21, 21:20.5, 21:20, 20:21", "crosses");

        // Adjacent antipodal vertices
        CheckLoopIsInvalid(new[] { new S2Point(1, 0, 0), new S2Point(-1, 0, 0), new S2Point(0, 0, 1) }, "antipodal");
    }

#if GTEST_HAS_DEATH_TEST
    [Fact]
    internal void Test_S2LoopDeathTest_IsValidDetectsInvalidLoops()
    {
        // Points with non-unit length (triggers S2_DCHECK failure in debug)
        EXPECT_DEBUG_DEATH(CheckLoopIsInvalid(new[] { new S2Point(2, 0, 0), new S2Point(0, 1, 0), new S2Point(0, 0, 1) }, "unit length"), "IsUnitLength");
    }
#endif

    [Fact]
    internal void Test_S2LoopTestBase_DistanceMethods()
    {
        // S2ClosestEdgeQuery is already tested, so just do a bit of sanity checking.

        // The empty and full loops don't have boundaries.
        TestDistanceMethods(empty_, new S2Point(0, 1, 0), new S2Point());
        TestDistanceMethods(full_, new S2Point(0, 1, 0), new S2Point());

        // A CCW square around the S2LatLng point (0,0).  Note that because lines of
        // latitude are curved on the sphere, it is not straightforward to project
        // points onto any edge except along the equator.  (The equator is the only
        // line of latitude that is also a geodesic.)
        var square = MakeLoopOrDie("-1:-1, -1:1, 1:1, 1:-1");
        Assert.True(square.IsNormalized());

        // A vertex.
        TestDistanceMethods(square, S2LatLng.FromDegrees(1, -1).ToPoint(), new S2Point());
        // A point on one of the edges.
        TestDistanceMethods(square, S2LatLng.FromDegrees(0.5, 1).ToPoint(), new S2Point());
        // A point inside the square.
        TestDistanceMethods(square, S2LatLng.FromDegrees(0, 0.5).ToPoint(), S2LatLng.FromDegrees(0, 1).ToPoint());
        // A point outside the square that projects onto an edge.
        TestDistanceMethods(square, S2LatLng.FromDegrees(0, -2).ToPoint(), S2LatLng.FromDegrees(0, -1).ToPoint());
        // A point outside the square that projects onto a vertex.
        TestDistanceMethods(square, S2LatLng.FromDegrees(3, 4).ToPoint(), S2LatLng.FromDegrees(1, 1).ToPoint());
    }

    [Fact]
    internal void Test_S2LoopTestBase_MakeRegularLoop()
    {
        S2Point center = S2LatLng.FromDegrees(80, 135).ToPoint();
        S1Angle radius = S1Angle.FromDegrees(20);
        var loop = S2Loop.MakeRegularLoop(center, radius, 4);

        Assert.Equal(4, loop.NumVertices);
        S2Point p0 = loop.Vertex(0);
        S2Point p1 = loop.Vertex(1);
        S2Point p2 = loop.Vertex(2);
        S2Point p3 = loop.Vertex(3);
        // Make sure that the radius is correct.
        Assert2.Near(20.0, new S2LatLng(center).GetDistance(new S2LatLng(p0)).GetDegrees());
        Assert2.Near(20.0, new S2LatLng(center).GetDistance(new S2LatLng(p1)).GetDegrees());
        Assert2.Near(20.0, new S2LatLng(center).GetDistance(new S2LatLng(p2)).GetDegrees());
        Assert2.Near(20.0, new S2LatLng(center).GetDistance(new S2LatLng(p3)).GetDegrees());
        // Make sure that all angles of the polygon are the same.
        Assert2.Near(S2.M_PI_2, (p1 - p0).Angle(p3 - p0));
        Assert2.Near(S2.M_PI_2, (p2 - p1).Angle(p0 - p1));
        Assert2.Near(S2.M_PI_2, (p3 - p2).Angle(p1 - p2));
        Assert2.Near(S2.M_PI_2, (p0 - p3).Angle(p2 - p3));
        // Make sure that all edges of the polygon have the same length.
        Assert2.Near(27.990890717782829, new S2LatLng(p0).GetDistance(new S2LatLng(p1)).GetDegrees());
        Assert2.Near(27.990890717782829, new S2LatLng(p1).GetDistance(new S2LatLng(p2)).GetDegrees());
        Assert2.Near(27.990890717782829, new S2LatLng(p2).GetDistance(new S2LatLng(p3)).GetDegrees());
        Assert2.Near(27.990890717782829, new S2LatLng(p3).GetDistance(new S2LatLng(p0)).GetDegrees());

        // Check actual coordinates. This may change if we switch the algorithm
        // intentionally.
        Assert2.Near(62.162880741097204, new S2LatLng(p0).Lat().GetDegrees());
        Assert2.Near(103.11051028343407, new S2LatLng(p0).Lng().GetDegrees());
        Assert2.Near(61.955157772928345, new S2LatLng(p1).Lat().GetDegrees());
        Assert2.Near(165.25681963683536, new S2LatLng(p1).Lng().GetDegrees());
        Assert2.Near(75.139812547718478, new S2LatLng(p2).Lat().GetDegrees());
        Assert2.Near(-119.13042521187423, new S2LatLng(p2).Lng().GetDegrees());
        Assert2.Near(75.524190079054392, new S2LatLng(p3).Lat().GetDegrees());
        Assert2.Near(26.392175948257943, new S2LatLng(p3).Lng().GetDegrees());
    }

    [Fact]
    internal void Test_S2LoopShape_Basic()
    {
        var loop = MakeLoopOrDie("0:0, 0:1, 1:0");
        S2Loop.Shape shape = new(loop);
        Assert.Equal(loop, shape.Loop);
        Assert.Equal(3, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.Equal(0, shape.GetChain(0).Start);
        Assert.Equal(3, shape.GetChain(0).Length);
        var edge2 = shape.GetEdge(2);
        Assert.Equal("1:0", edge2.V0.ToDebugString());
        Assert.Equal("0:0", edge2.V1.ToDebugString());
        Assert.Equal(2, shape.Dimension());
        Assert.False(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2LoopShape_EmptyLoop()
    {
        S2Loop loop = S2Loop.kEmpty;
        var shape = new S2Loop.Shape(loop);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(0, shape.NumChains());
        Assert.True(shape.IsEmpty());
        Assert.False(shape.IsFull());
        Assert.False(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2LoopShape_FullLoop()
    {
        S2Loop loop = S2Loop.kFull;
        S2Loop.Shape shape = new(loop);
        Assert.Equal(0, shape.NumEdges());
        Assert.Equal(1, shape.NumChains());
        Assert.False(shape.IsEmpty());
        Assert.True(shape.IsFull());
        Assert.True(shape.GetReferencePoint().Contained);
    }

    [Fact]
    internal void Test_S2LoopOwningShape_Ownership()
    {
        // Debug mode builds will catch any memory leak below.
        var loop = S2Loop.kEmpty;
        var shape = new S2Loop.Shape(loop);
        _logger.WriteLine("sghape.Id: " + shape.Id);
    }

    private S2Loop AddLoop(string str)
    {
        return AddLoop(MakeLoopOrDie(str));
    }

    private S2Loop AddLoop(S2Loop loop)
    {
        all_loops.Add(loop);
        return loop;
    }

    // Wrapper function that encodes "loop" into "encoder" using the private
    // EncodeCompressed() method.
    private static void TestEncodeCompressed(S2Loop loop, int level, Encoder encoder)
    {
        var points = new S2PointCompression.S2XYZFaceSiTi[loop.NumVertices];
        loop.GetXYZFaceSiTiVertices(points, 0);
        loop.EncodeCompressed(encoder, points, 0, level);
    }

    // Wrapper function that decodes the contents of "encoder" into "loop" using
    // the private DecodeCompressed() method.
    private static void TestDecodeCompressed(Encoder encoder, int level, out S2Loop loop)
    {
        var decoder = encoder.GetDecoder();
        var (success, loop_) = S2Loop.DecodeCompressed(decoder, level);
        Assert.True(success);
        loop = loop_!;
    }

    private static void Rotate(ref S2Loop loop)
    {
        var vertices = new List<S2Point>();
        for (int i = 1; i < loop.NumVertices; ++i)
        {
            vertices.Add(loop.Vertex(i));
        }
        vertices.Add(loop.Vertex(0));
        loop = new S2Loop(vertices);
    }

    // Check that the curvature is *identical* when the vertex order is
    // rotated, and that the sign is inverted when the vertices are reversed.
    private static void CheckCurvatureInvariants(S2Loop loop)
    {
        var expected = loop.Curvature();
        var loop_copy = (S2Loop)loop.CustomClone();
        for (int i = 0; i < loop.NumVertices; ++i)
        {
            loop_copy.Invert();
            Assert.Equal(-expected, loop_copy.Curvature());
            loop_copy.Invert();
            Rotate(ref loop_copy);
            Assert.Equal(expected, loop_copy.Curvature());
        }
    }

    // Checks that if a loop is normalized, it doesn't contain a
    // point outside of it, and vice versa.
    private static void CheckNormalizeAndContains(S2Loop loop)
    {
        S2Point p = MakePointOrDie("40:40");

        var flip = (S2Loop)loop.CustomClone();
        flip.Invert();
        Assert.True(loop.IsNormalized() ^ loop.Contains(p));
        Assert.True(flip.IsNormalized() ^ flip.Contains(p));

        Assert.True(loop.IsNormalized() ^ flip.IsNormalized());

        flip.Normalize();
        Assert.False(flip.Contains(p));
    }

    // Given a pair of loops where A contains B, check various identities.
    private static void TestOneNestedPair(S2Loop a, S2Loop b)
    {
        Assert.True(a.Contains(b));
        Assert.Equal(a.BoundaryEquals(b), b.Contains(a));
        Assert.Equal(!b.IsEmpty(), a.Intersects(b));
        Assert.Equal(!b.IsEmpty(), b.Intersects(a));
    }

    // Given a pair of disjoint loops A and B, check various identities.
    private static void TestOneDisjointPair(S2Loop a, S2Loop b)
    {
        Assert.False(a.Intersects(b));
        Assert.False(b.Intersects(a));
        Assert.Equal(b.IsEmpty(), a.Contains(b));
        Assert.Equal(a.IsEmpty(), b.Contains(a));
    }

    // Given loops A and B whose union covers the sphere, check various identities.
    private static void TestOneCoveringPair(S2Loop a, S2Loop b)
    {
        Assert.Equal(a.IsFull(), a.Contains(b));
        Assert.Equal(b.IsFull(), b.Contains(a));
        var a1 = (S2Loop)a.CustomClone();
        a1.Invert();
        bool complementary = a1.BoundaryEquals(b);
        Assert.Equal(!complementary, a.Intersects(b));
        Assert.Equal(!complementary, b.Intersects(a));
    }

    // Given loops A and B such that both A and its complement intersect both B
    // and its complement, check various identities.
    private static void TestOneOverlappingPair(S2Loop a, S2Loop b)
    {
        Assert.False(a.Contains(b));
        Assert.False(b.Contains(a));
        Assert.True(a.Intersects(b));
        Assert.True(b.Intersects(a));
    }

    // Given a pair of loops where A contains B, test various identities
    // involving A, B, and their complements.
    private static void TestNestedPair(S2Loop a, S2Loop b)
    {
        var a1 = (S2Loop)a.CustomClone();
        var b1 = (S2Loop)b.CustomClone();
        a1.Invert();
        b1.Invert();
        TestOneNestedPair(a, b);
        TestOneNestedPair(b1, a1);
        TestOneDisjointPair(a1, b);
        TestOneCoveringPair(a, b1);
    }

    // Given a pair of disjoint loops A and B, test various identities
    // involving A, B, and their complements.
    private static void TestDisjointPair(S2Loop a, S2Loop b)
    {
        var a1 = (S2Loop)a.CustomClone();
        a1.Invert();
        TestNestedPair(a1, b);
    }

    // Given loops A and B whose union covers the sphere, test various identities
    // involving A, B, and their complements.
    private static void TestCoveringPair(S2Loop a, S2Loop b)
    {
        var b1 = (S2Loop)b.CustomClone();
        b1.Invert();
        TestNestedPair(a, b1);
    }

    // Given loops A and B such that both A and its complement intersect both B
    // and its complement, test various identities involving these four loops.
    private static void TestOverlappingPair(S2Loop a, S2Loop b)
    {
        var a1 = (S2Loop)a.CustomClone();
        var b1 = (S2Loop)b.CustomClone();
        a1.Invert();
        b1.Invert();
        TestOneOverlappingPair(a, b);
        TestOneOverlappingPair(a1, b1);
        TestOneOverlappingPair(a1, b);
        TestOneOverlappingPair(a, b1);
    }

    private void TestRelation(S2Loop a, S2Loop b, RelationFlags flags, bool shared_edge)
    {
        TestRelationWithDesc(a, b, flags, shared_edge, $"args {a}, {b}");
    }

    // Verify the relationship between two loops A and B.  "flags" is the set of
    // RelationFlags that apply.  "shared_edge" means that the loops share at
    // least one edge (possibly reversed).
    private void TestRelationWithDesc(S2Loop a, S2Loop b, RelationFlags flags, bool shared_edge, string test_description)
    {
        _logger.WriteLine(test_description);
        if ((flags & RelationFlags.CONTAINS) != 0)
        {
            TestNestedPair(a, b);
        }
        if ((flags & RelationFlags.CONTAINED) != 0)
        {
            TestNestedPair(b, a);
        }
        if ((flags & RelationFlags.COVERS) != 0)
        {
            TestCoveringPair(a, b);
        }
        if ((flags & RelationFlags.DISJOINT) != 0)
        {
            TestDisjointPair(a, b);
        }
        else if ((flags & (RelationFlags.CONTAINS | RelationFlags.CONTAINED | RelationFlags.COVERS)) == 0)
        {
            TestOverlappingPair(a, b);
        }
        if (!shared_edge && ((flags & (RelationFlags.CONTAINS | RelationFlags.CONTAINED | RelationFlags.DISJOINT)) != 0))
        {
            Assert.Equal(a.Contains(b), a.ContainsNested(b));
        }
        // A contains the boundary of B if either A contains B, or the two loops
        // contain each other's boundaries and there are no shared edges (since at
        // least one such edge must be reversed, and therefore is not considered to
        // be contained according to the rules of CompareBoundary).
        int comparison = 0;
        if (((flags & RelationFlags.CONTAINS) != 0) || (((flags & RelationFlags.COVERS) != 0) && !shared_edge))
        {
            comparison = 1;
        }
        // Similarly, A excludes the boundary of B if either A and B are disjoint,
        // or B contains A and there are no shared edges (since A is considered to
        // contain such edges according to the rules of CompareBoundary).
        if (((flags & RelationFlags.DISJOINT) != 0) || (((flags & RelationFlags.CONTAINED) != 0) && !shared_edge))
        {
            comparison = -1;
        }
        // CompareBoundary requires that neither loop is empty.
        if (!a.IsEmpty() && !b.IsEmpty())
        {
            Assert.Equal(comparison, a.CompareBoundary(b));
        }
    }

    private static S2Loop MakeCellLoop(S2CellId begin, S2CellId end)
    {
        // Construct a CCW polygon whose boundary is the union of the cell ids
        // in the range [begin, end).  We add the edges one by one, removing
        // any edges that are already present in the opposite direction.

        var edges = new Dictionary<S2Point, List<S2Point>>();
        for (S2CellId id = begin; id != end; id = id.Next())
        {
            S2Cell cell = new(id);
            for (int k = 0; k < 4; ++k)
            {
                S2Point a = cell.Vertex(k);
                S2Point b = cell.Vertex(k + 1);
                if (!edges[b].Remove(a))
                {
                    edges[a].Add(b);
                }
                else if (!edges[b].Any())
                {
                    edges.Remove(b);
                }
            }
        }

        // The remaining edges form a single loop.  We simply follow it starting
        // at an arbitrary vertex and build up a list of vertices.

        var vertices = new List<S2Point>();
        S2Point p = edges.First().Value.First();
        while (edges.Any())
        {
            Assert.Single(edges[p]);
            S2Point next = edges[p].First();
            vertices.Add(p);
            edges.Remove(p);
            p = next;
        }

        return new S2Loop(vertices);
    }

    private static void TestNear(string a_str, string b_str, S1Angle max_error, bool expected)
    {
        var a = MakeLoopOrDie(a_str);
        var b = MakeLoopOrDie(b_str);
        Assert.Equal(a.BoundaryNear(b, max_error), expected);
        Assert.Equal(b.BoundaryNear(a, max_error), expected);
    }

    private static void CheckIdentical(S2Loop loop, S2Loop loop2)
    {
        Assert.Equal(loop.Depth, loop2.Depth);
        Assert.Equal(loop.NumVertices, loop2.NumVertices);
        for (int i = 0; i < loop.NumVertices; ++i)
        {
            Assert.Equal(loop.Vertex(i), loop2.Vertex(i));
        }
        Assert.Equal(loop.IsEmpty(), loop2.IsEmpty());
        Assert.Equal(loop.IsFull(), loop2.IsFull());
        Assert.Equal(loop.Depth, loop2.Depth);
        Assert.Equal(loop.IsNormalized(), loop2.IsNormalized());
        Assert.Equal(loop.Contains(S2.Origin), loop2.Contains(S2.Origin));
        Assert.Equal(loop.GetRectBound(), loop2.GetRectBound());
    }

    private static void TestEncodeDecode(S2Loop loop)
    {
        Encoder encoder = new();
        loop.Encode(encoder);
        var decoder = encoder.GetDecoder();
        var (success, loop2) = S2Loop.Decode(decoder);
        Assert.True(success);
        CheckIdentical(loop, loop2);
    }

    private static void TestEmptyFullSnapped(S2Loop loop, int level)
    {
        Assert.True(loop.IsEmptyOrFull());
        S2CellId cellid = new S2CellId(loop.Vertex(0)).Parent(level);
        S2Point[] vertices = { cellid.ToPoint() };
        S2Loop loop2 = new(vertices);
        Assert.True(loop.BoundaryEquals(loop2));
        Assert.True(loop.BoundaryApproxEquals(loop2));
        Assert.True(loop.BoundaryNear(loop2));
    }

    // Test converting the empty/full loops to S2LatLng representations.  (We
    // don't bother testing E5/E6/E7 because that test is less demanding.)
    private static void TestEmptyFullLatLng(S2Loop loop)
    {
        Assert.True(loop.IsEmptyOrFull());
        S2Point[] vertices = { new S2LatLng(loop.Vertex(0)).ToPoint() };
        S2Loop loop2 = new(vertices);
        Assert.True(loop.BoundaryEquals(loop2));
        Assert.True(loop.BoundaryApproxEquals(loop2));
        Assert.True(loop.BoundaryNear(loop2));
    }

    private static void TestEmptyFullConversions(S2Loop loop)
    {
        TestEmptyFullSnapped(loop, S2.kMaxCellLevel);
        TestEmptyFullSnapped(loop, 1);  // Worst case for approximation
        TestEmptyFullSnapped(loop, 0);
        TestEmptyFullLatLng(loop);
    }

    // Construct a loop using S2TextFormat.MakeLoop(str) and check that it
    // produces a validation error that includes "snippet".
    private static void CheckLoopIsInvalid(string str, string snippet)
    {
        var loop = MakeLoopOrDie(str, S2Debug.DISABLE);
        Assert.True(loop.FindValidationError(out var error));
        Assert.NotEqual(-1, error.Text.IndexOf(snippet));
    }

    private static void CheckLoopIsInvalid(S2Point[] points, string snippet)
    {
        S2Loop l = new(points, S2Debug.DISABLE);
        Assert.True(l.FindValidationError(out var error));
        Assert.Contains(snippet, error.Text);
    }

    // Helper function for testing the distance methods.  "boundary_x" is the
    // expected result of projecting "x" onto the loop boundary.  For convenience
    // it can be set to S2Point() to indicate that (boundary_x == x).
    private static void TestDistanceMethods(S2Loop loop, S2Point x, S2Point boundary_x)
    {
        // This error is not guaranteed by the implementation but is okay for tests.
        S1Angle kMaxError = S1Angle.FromRadians(1e-15);

        if (boundary_x == new S2Point()) boundary_x = x;
        Assert.True(new S1Angle(boundary_x, loop.ProjectToBoundary(x)) <= kMaxError);

        if (loop.IsEmptyOrFull())
        {
            Assert.Equal(S1Angle.Infinity, loop.DistanceToBoundary(x));
        }
        else
        {
            // Assert.True.Near only works with doubles.
            Assert2.Near(new S1Angle(x, boundary_x).GetDegrees(), loop.DistanceToBoundary(x).GetDegrees(), kMaxError.GetDegrees());
        }
        if (loop.Contains(x))
        {
            Assert.Equal(S1Angle.Zero, loop.Distance(x));
            Assert.Equal(x, loop.Project(x));
        }
        else
        {
            Assert.Equal(loop.DistanceToBoundary(x), loop.Distance(x));
            Assert.Equal(loop.ProjectToBoundary(x), loop.Project(x));
        }
    }

    private enum RelationFlags
    {
        CONTAINS = 0x01,  // A contains B
        CONTAINED = 0x02,  // B contains A
        DISJOINT = 0x04,  // A and B are disjoint (intersection is empty)
        COVERS = 0x08,  // (A union B) covers the entire sphere
    };
}
