// There are several notions of the "centroid" of a triangle.  First, there
// is the planar centroid, which is simply the centroid of the ordinary
// (non-spherical) triangle defined by the three vertices.  Second, there is
// the surface centroid, which is defined as the intersection of the three
// medians of the spherical triangle.  It is possible to show that this
// point is simply the planar centroid projected to the surface of the
// sphere.  Finally, there is the true centroid (mass centroid), which is
// defined as the surface integral over the spherical triangle of (x,y,z)
// divided by the triangle area.  This is the point that the triangle would
// rotate around if it was spinning in empty space.
//
// The best centroid for most purposes is the true centroid.  Unlike the
// planar and surface centroids, the true centroid behaves linearly as
// regions are added or subtracted.  That is, if you split a triangle into
// pieces and compute the average of their centroids (weighted by triangle
// area), the result equals the centroid of the original triangle.  This is
// not true of the other centroids.
//
// Also note that the surface centroid may be nowhere near the intuitive
// "center" of a spherical triangle.  For example, consider the triangle
// with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
// The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
// within a distance of 2*eps of the vertex B.  Note that the median from A
// (the segment connecting A to the midpoint of BC) passes through S, since
// this is the shortest path connecting the two endpoints.  On the other
// hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
// the surface is a much more reasonable interpretation of the "center" of
// this triangle.

namespace S2Geometry;

public static class S2Centroid
{
    // Returns the centroid of the planar triangle ABC.  This can be normalized to
    // unit length to obtain the "surface centroid" of the corresponding spherical
    // triangle, i.e. the intersection of the three medians.  However, note that
    // for large spherical triangles the surface centroid may be nowhere near the
    // intuitive "center" (see example above).
    public static S2Point PlanarCentroid(S2Point a, S2Point b, S2Point c)
    {
        return 1.0 / 3 * (a + b + c);
    }

    // Returns the true centroid of the spherical triangle ABC multiplied by the
    // signed area of spherical triangle ABC.  The reasons for multiplying by the
    // signed area are (1) this is the quantity that needs to be summed to compute
    // the centroid of a union or difference of triangles, and (2) it's actually
    // easier to calculate this way.  All points must have unit length.
    //
    // Note that the result of this function is defined to be S2Point.Empty if
    // the triangle is degenerate (and that this is intended behavior).
    public static S2Point TrueCentroid(S2Point a, S2Point b, S2Point c)
    {
        MyDebug.Assert(a.IsUnitLength());
        MyDebug.Assert(b.IsUnitLength());
        MyDebug.Assert(c.IsUnitLength());

        // I couldn't find any references for computing the true centroid of a
        // spherical triangle...  I have a truly marvellous demonstration of this
        // formula which this margin is too narrow to contain :)

        // Use Angle() in order to get accurate results for small triangles.
        var angle_a = b.Angle(c);
        var angle_b = c.Angle(a);
        var angle_c = a.Angle(b);
        var ra = (angle_a == 0) ? 1 : (angle_a / Math.Sin(angle_a));
        var rb = (angle_b == 0) ? 1 : (angle_b / Math.Sin(angle_b));
        var rc = (angle_c == 0) ? 1 : (angle_c / Math.Sin(angle_c));

        // Now compute a point M such that:
        //
        //  [Ax Ay Az] [Mx]                       [ra]
        //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
        //  [Cx Cy Cz] [Mz]                       [rc]
        //
        // To improve the numerical stability we subtract the first row (A) from the
        // other two rows; this reduces the cancellation error when A, B, and C are
        // very close together.  Then we solve it using Cramer's rule.
        //
        // The result is the true centroid of the triangle multiplied by the
        // triangle's area.
        //
        // TODO(b/205027737): This code still isn't as numerically stable as it could
        // be.  The biggest potential improvement is to compute B-A and C-A more
        // accurately so that (B-A)x(C-A) is always inside triangle ABC.
        var x = new S2Point(a.X, b.X - a.X, c.X - a.X);
        var y = new S2Point(a.Y, b.Y - a.Y, c.Y - a.Y);
        var z = new S2Point(a.Z, b.Z - a.Z, c.Z - a.Z);
        var r = new S2Point(ra, rb - ra, rc - ra);
        return 0.5 * new S2Point(
            y.CrossProd(z).DotProd(r),
            z.CrossProd(x).DotProd(r),
            x.CrossProd(y).DotProd(r));
    }

    // Returns the true centroid of the spherical geodesic edge AB multiplied by
    // the length of the edge AB.  As with triangles, the true centroid of a
    // collection of edges may be computed simply by summing the result of this
    // method for each edge.
    //
    // Note that the planar centroid of a geodesic edge simply 0.5 * (a + b),
    // while the surface centroid is (a + b).Normalize().  However neither of
    // these values is appropriate for computing the centroid of a collection of
    // edges (such as a polyline).
    //
    // Also note that the result of this function is defined to be S2Point.Empty
    // if the edge is degenerate (and that this is intended behavior).
    public static S2Point TrueCentroid(S2Point a, S2Point b)
    {
        // The centroid (multiplied by length) is a vector toward the midpoint
        // of the edge, whose length is twice the sine of half the angle between
        // the two vertices.  Defining theta to be this angle, we have:
        var vdiff = a - b;  // Length == 2*sin(theta)
        var vsum = a + b;   // Length == 2*cos(theta)
        var sin2 = vdiff.Norm2();
        var cos2 = vsum.Norm2();
        if (cos2 == 0) return S2Point.Empty;  // Ignore antipodal edges.
        return Math.Sqrt(sin2 / cos2) * vsum;  // Length == 2*sin(theta)
    }
}
