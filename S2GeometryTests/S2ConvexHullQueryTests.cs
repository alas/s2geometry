namespace S2Geometry;

public class S2ConvexHullQueryTests
{
    private readonly ITestOutputHelper _logger;

    public S2ConvexHullQueryTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    public void Test_S2ConvexHullQuery_NoPoints()
    {
        S2ConvexHullQuery query = new();
        var result = (query.GetConvexHull());
        Assert.True(result.IsEmpty());
    }

    [Fact]
    public void Test_S2ConvexHullQuery_OnePoint()
    {
        S2ConvexHullQuery query = new();
        var p = new S2Point(0, 0, 1);
        query.AddPoint(p);
        var result = query.GetConvexHull();
        Assert.Equal(3, result.NumVertices);
        Assert.True(result.IsNormalized());
        Assert.True(LoopHasVertex(result, p));
        // Add some duplicate points and check that the result is the same.
        query.AddPoint(p);
        query.AddPoint(p);
        var result2 = (query.GetConvexHull());
        Assert.True(result2 == result);
    }

    [Fact]
    public void Test_S2ConvexHullQuery_TwoPoints()
    {
        S2ConvexHullQuery query = new();
        S2Point p = new(0, 0, 1);
        S2Point q = new(0, 1, 0);
        query.AddPoint(p);
        query.AddPoint(q);
        var result = (query.GetConvexHull());
        Assert.Equal(3, result.NumVertices);
        Assert.True(result.IsNormalized());
        Assert.True(LoopHasVertex(result, p));
        Assert.True(LoopHasVertex(result, q));
        // Add some duplicate points and check that the result is the same.
        query.AddPoint(q);
        query.AddPoint(p);
        query.AddPoint(p);
        var result2 = query.GetConvexHull();
        Assert.True(result2 == result);
    }

    [Fact]
    public void Test_S2ConvexHullQuery_TwoAntipodalPoints()
    {
        S2ConvexHullQuery query=new();
        query.AddPoint(new S2Point(0, 0, 1));
        query.AddPoint(new S2Point(0, 0, -1));
        var result = query.GetConvexHull();
        Assert.True(result.IsFull());
    }

    [Fact]
    public void Test_S2ConvexHullQuery_EmptyLoop()
    {
        S2ConvexHullQuery query = new();
        S2Loop empty = (S2Loop.kEmpty);
        query.AddLoop(empty);
        var result = (query.GetConvexHull());
        Assert.True(result.IsEmpty());
    }

    [Fact]
    public void Test_S2ConvexHullQuery_FullLoop()
    {
        S2ConvexHullQuery query = new();
        S2Loop full = (S2Loop.kFull);
        query.AddLoop(full);
        var result = (query.GetConvexHull());
        Assert.True(result.IsFull());
    }

    [Fact]
    public void Test_S2ConvexHullQuery_EmptyPolygon()
    {
        S2ConvexHullQuery query = new();
        List<S2Loop> loops = new();
        S2Polygon empty = new(loops);
        query.AddPolygon(empty);
        var result = (query.GetConvexHull());
        Assert.True(result.IsEmpty());
    }

    [Fact]
    public void Test_S2ConvexHullQuery_NonConvexPoints()
    {
        // Generate a point set such that the only convex region containing them is
        // the entire sphere.  In other words, you can generate any point on the
        // sphere by repeatedly linearly interpolating between the points.  (The
        // four points of a tetrahedron would also work, but this is easier.)
        S2ConvexHullQuery query = new();
        for (int face = 0; face < 6; ++face)
        {
            query.AddPoint(S2CellId.FromFace(face).ToPoint());
        }
        var result = (query.GetConvexHull());
        Assert.True(result.IsFull());
    }

    [Fact]
    public void Test_S2ConvexHullQuery_SimplePolyline()
    {
        // A polyline is handling identically to a point set, so there is no need
        // for special testing other than code coverage.
        var polyline = (MakePolylineOrDie(
  "0:1, 0:9, 1:6, 2:6, 3:10, 4:10, 5:5, 4:0, 3:0, 2:5, 1:5"));
        S2ConvexHullQuery query = new();
        query.AddPolyline(polyline);
        var result = (query.GetConvexHull());
        var expected_result = (
            MakeLoopOrDie("0:1, 0:9, 3:10, 4:10, 5:5, 4:0, 3:0"));
        Assert.True(result.BoundaryEquals(expected_result));
    }

    [Fact]
    public void Test_S2ConvexHullQuery_LoopsAroundNorthPole()
    {
        // Test loops of various sizes around the north pole.
        TestNorthPoleLoop(S1Angle.FromDegrees(1), 3);
        TestNorthPoleLoop(S1Angle.FromDegrees(89), 3);

        // The following two loops should yield the full loop.
        TestNorthPoleLoop(S1Angle.FromDegrees(91), 3);
        TestNorthPoleLoop(S1Angle.FromDegrees(179), 3);

        TestNorthPoleLoop(S1Angle.FromDegrees(10), 100);
        TestNorthPoleLoop(S1Angle.FromDegrees(89), 1000);
    }

    [Fact]
    public void Test_S2ConvexHullQuery_PointsInsideHull()
    {
        // Repeatedly build the convex hull of a set of points, then add more points
        // inside that loop and build the convex hull again.  The result should
        // always be the same.
        int kIters = 1000;
        for (int iter = 0; iter < kIters; ++iter)
        {
            S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.

            // Choose points from within a cap of random size, up to but not including
            // an entire hemisphere.
            S2Cap cap = S2Testing.GetRandomCap(S2.DoubleError, 1.999 * Math.PI);
            S2ConvexHullQuery query = new();
            int num_points1 = S2Testing.Random.Uniform(100) + 3;
            for (int i = 0; i < num_points1; ++i)
            {
                query.AddPoint(S2Testing.SamplePoint(cap));
            }
            var hull = query.GetConvexHull();

            // When the convex hull is nearly a hemisphere, the algorithm sometimes
            // returns a full cap instead.  This is because it first computes a
            // bounding rectangle for all the input points/edges and then converts it
            // to a bounding cap, which sometimes yields a non-convex cap (radius
            // larger than 90 degrees).  This should not be a problem in practice
            // (since most convex hulls are not hemispheres), but in order make this
            // test pass reliably it means that we need to reject convex hulls whose
            // bounding cap (when computed from a bounding rectangle) is not convex.
            //
            // TODO(b/203702905): This test can still fail (about 1 iteration in
            // 500,000) because the S2LatLngRect::GetCapBound implementation does not
            // guarantee that A.Contains(B) implies
            // A.GetCapBound().Contains(B.GetCapBound()).
            if (hull.GetCapBound().Height() >= 1) continue;

            // Otherwise, add more points inside the convex hull.
            int num_points2 = 1000;
            for (int i = 0; i < num_points2; ++i)
            {
                S2Point p = S2Testing.SamplePoint(cap);
                if (hull.Contains(p))
                {
                    query.AddPoint(p);
                }
            }
            // Finally, build a new convex hull and check that it hasn't changed.
            var hull2 = (query.GetConvexHull());
            _logger.WriteLine("Iteration: " + iter);
            Assert.True(hull2.BoundaryEquals(hull));
        }
    }

    [Fact]
    public void Test_S2ConvexHullQuery_CapBoundExpandedToHemisphere()
    {
        // The following 3 points yield an S2Cap bound that is slightly smaller than
        // a hemisphere.  Here we test that the cap is expanded using a conservative
        // error bound to yield a hemisphere, which causes the convex hull algorithm
        // to return the full sphere.
        S2ConvexHullQuery query=new();
        query.AddPoint(MakePointOrDie("0:0"));
        query.AddPoint(MakePointOrDie("0:45"));
        query.AddPoint(MakePointOrDie("0:-135"));
        var result = query.GetConvexHull();
        Assert.True(result.IsFull());
    }

    private static bool LoopHasVertex(S2Loop loop, S2Point p)
    {
        for (int i = 0; i < loop.NumVertices; ++i)
        {
            if (loop.Vertex(i) == p) return true;
        }
        return false;
    }

    private static void TestNorthPoleLoop(S1Angle radius, int num_vertices)
    {
        // If the radius is very close to 90, then it's hard to predict whether the
        // result will be the full loop or not.
        Assert.True(Math.Abs(radius.Radians - S2.M_PI_2) >= S2.DoubleError);

        S2ConvexHullQuery query = new();
        var loop = S2Loop.MakeRegularLoop(new S2Point(0, 0, 1), radius, num_vertices);
        query.AddLoop(loop);
        var result = query.GetConvexHull();
        if (radius > S1Angle.FromRadians(S2.M_PI_2))
        {
            Assert.True(result.IsFull());
        }
        else
        {
            Assert.True(result.BoundaryEquals(loop));
        }
    }
}
