namespace S2Geometry;

public class S2CentroidTests
{
    [Fact]
    internal void Test_PlanarCentroid_SemiEquator()
    {
        // Test the centroid of polyline ABC that follows the equator and consists
        // of two 90 degree edges (i.e., C = -A).  The centroid should point toward
        // B and have a norm of 1/3.  This is not a thorough test of
        // PlanarCentroid; it is only intended to prevent it from being detected as
        // dead code.
        S2Point a=new(0, -1, 0), b=new(1, 0, 0), c=new(0, 1, 0);
        S2Point centroid = S2Centroid.PlanarCentroid(a, b, c);
        Assert.True(S2.ApproxEquals(b, centroid.Normalize()));
        Assert.Equal(1 / 3.0, centroid.Norm());
    }

    [Fact]
    internal void Test_TriangleTrueCentroid_SmallTriangles() {
        // Test TrueCentroid() with very small triangles.  This test assumes that
        // the triangle is small enough so that it is nearly planar.
        for (int iter = 0; iter < 100; ++iter) {
            S2Testing.GetRandomFrame(out var p, out var x, out var y);
            double d = 1e-4 * Math.Pow(1e-4, S2Testing.Random.RandDouble());
            S2Point p0 = (p - d * x).Normalize();
            S2Point p1 = (p + d * x).Normalize();
            S2Point p2 = (p + 3 * d * y).Normalize();
            S2Point centroid = S2Centroid.TrueCentroid(p0, p1, p2).Normalize();

            // The centroid of a planar triangle is at the intersection of its
            // medians, which is two-thirds of the way along each median.
            S2Point expected_centroid = (p + d * y).Normalize();
            Assert.True(centroid.Angle(expected_centroid) <= 2e-8);
        }
    }

    [Fact]
    internal void Test_EdgeTrueCentroid_SemiEquator() {
        // Test the centroid of polyline ABC that follows the equator and consists
        // of two 90 degree edges (i.e., C = -A).  The centroid (multiplied by
        // length) should point toward B and have a norm of 2.0.  (The centroid
        // itself has a norm of 2/Pi, and the total edge length is Pi.)
        S2Point a = new(0, -1, 0), b = new(1, 0, 0), c = new(0, 1, 0);
        S2Point centroid = S2Centroid.TrueCentroid(a, b) + S2Centroid.TrueCentroid(b, c);
        Assert.True(S2.ApproxEquals(b, centroid.Normalize()));
        Assert2.DoubleEqual(2.0, centroid.Norm());
    }

    [Fact]
    internal void Test_EdgeTrueCentroid_GreatCircles() {
        // Construct random great circles and divide them randomly into segments.
        // Then make sure that the centroid is approximately at the center of the
        // sphere.  Note that because of the way the centroid is computed, it does
        // not matter how we split the great circle into segments.
        //
        // Note that this is a direct test of the properties that the centroid
        // should have, rather than a test that it matches a particular formula.

        for (int iter = 0; iter < 100; ++iter) {
            // Choose a coordinate frame for the great circle.
            S2Point centroid = S2Point.Empty;
            S2Testing.GetRandomFrame(out var x, out var y, out _);

            S2Point v0 = x;
            for (double theta = 0; theta < S2.M_2_PI;
                 theta += Math.Pow(S2Testing.Random.RandDouble(), 10)) {
                S2Point v1 = Math.Cos(theta) * x + Math.Sin(theta) * y;
                centroid += S2Centroid.TrueCentroid(v0, v1);
                v0 = v1;
            }
            // Close the circle.
            centroid += S2Centroid.TrueCentroid(v0, x);
            Assert.True(centroid.Norm() <= 2e-14);
        }
    }
}
