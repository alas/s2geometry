// Defines various measures for polylines on the sphere.  These are low-level
// methods that work directly with arrays of S2Points.  They are used to
// implement the methods in s2shapeindex_measures.h, s2shape_measures.h, and
// s2polyline.h.
//
// See s2loop_measures.h, s2edge_distances.h, and s2measures.h for additional
// low-level methods.

namespace S2Geometry;

public static class S2PolylineMeasures
{
    // Returns the length of the polyline.  Returns zero for polylines with fewer
    // than two vertices.
    public static S1Angle GetLength(IList<S2Point> polyline)
    {
        var length = S1Angle.Zero;
        for (int i = 1; i < polyline.Count; ++i)
        {
            length += new S1Angle(polyline[i - 1], polyline[i]);
        }
        return length;
    }

    // Returns the true centroid of the polyline multiplied by the length of the
    // polyline (see s2centroids.h for details on centroids).  The result is not
    // unit length, so you may want to normalize it.
    //
    // Scaling by the polyline length makes it easy to compute the centroid of
    // several polylines (by simply adding up their centroids).
    //
    // CAVEAT: Returns S2Point() for degenerate polylines (e.g., AA).  [Note that
    // this answer is correct; the result of this function is a line integral over
    // the polyline, whose value is always zero if the polyline is degenerate.]
    public static S2Point GetCentroid(IList<S2Point> polyline)
    {
        var centroid = S2Point.Empty;
        for (int i = 1; i < polyline.Count; ++i)
        {
            centroid += S2Centroid.TrueCentroid(polyline[i - 1], polyline[i]);
        }
        return centroid;
    }
}
