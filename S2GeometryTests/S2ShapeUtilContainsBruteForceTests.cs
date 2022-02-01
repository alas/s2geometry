namespace S2Geometry;

public class S2ShapeUtilContainsBruteForceTests
{
    [Fact]
    public void Test_ContainsBruteForce_NoInterior()
    {
        // Defines a polyline that almost entirely encloses the point 0:0.
        var polyline = MakeLaxPolylineOrDie("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
        Assert.False(polyline.ContainsBruteForce(MakePointOrDie("0:0")));
    }

    [Fact]
    public void Test_ContainsBruteForce_ContainsReferencePoint()
    {
        // Checks that ContainsBruteForce agrees with GetReferencePoint.
        var polygon = MakeLaxPolygonOrDie("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
        var rp = polygon.GetReferencePoint();
        Assert.Equal(rp.Contained, polygon.ContainsBruteForce(rp.Point));
    }

    [Fact]
    public void Test_ContainsBruteForce_ConsistentWithS2Loop()
    {
        // Checks that ContainsBruteForce agrees with S2Loop.Contains().
        var loop = S2Loop.MakeRegularLoop(MakePointOrDie("89:-179"), S1Angle.FromDegrees(10), 100);
        S2Loop.Shape shape = new(loop);
        for (int i = 0; i < loop.NumVertices; ++i)
        {
            Assert.Equal(loop.Contains(loop.Vertex(i)),
                      shape.ContainsBruteForce(loop.Vertex(i)));
        }
    }
}
