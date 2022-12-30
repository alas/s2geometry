// Defines various angle and area measures for S2ShapeIndex objects.  In
// general, these methods return the sum of the corresponding measure for the
// S2Shapes in the index.

namespace S2Geometry;

public static class S2ShapeIndexMeasures
{
    // Returns the maximum dimension of any shape in the index.  Returns -1 if the
    // index does not contain any shapes.
    //
    // Note that the dimension does *not* depend on whether the shapes in the
    // index contain any points; for example, the dimension of an empty point set
    // is 0, and the dimension of an empty polygon is 2.
    public static int GetDimension(S2ShapeIndex index)
    {
        int dim = -1;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null) dim = Math.Max(dim, shape.Dimension());
        }
        return dim;
    }

    // Returns the number of points (objects of dimension zero) in the index.
    // Note that polyline and polygon vertices are *not* included in this count.
    public static int GetNumPoints(S2ShapeIndex index)
    {
        int count = 0;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null && shape.Dimension() == 0)
            {
                count += shape.NumEdges();
            }
        }
        return count;
    }

    // Returns the total length of all polylines in the index.  Returns zero if no
    // polylines are present.
    //
    // All edges are modeled as spherical geodesics.  The result can be converted
    // to a distance on the Earth's surface (with a worst-case error of 0.562%
    // near the equator) using the functions in s2earth.h.
    public static S1Angle GetLength(S2ShapeIndex index)
    {
        var length = S1Angle.Zero;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null) length += S2.GetLength(shape);
        }
        return length;
    }

    // Returns the total perimeter of all polygons in the index (including both
    // "shells" and "holes").  Returns zero if no polygons are present.
    //
    // All edges are modeled as spherical geodesics.  The result can be converted
    // to a distance on the Earth's surface (with a worst-case error of 0.562%
    // near the equator) using the functions in s2earth.h.
    public static S1Angle GetPerimeter(S2ShapeIndex index)
    {
        var perimeter = S1Angle.Zero;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null) perimeter += S2.GetPerimeter(shape);
        }
        return perimeter;
    }

    // Returns the total area of all polygons in the index.  Returns zero if no
    // polygons are present.  This method has good relative accuracy for both very
    // large and very small regions.  Note that the result may exceed 4*Pi if the
    // index contains overlapping polygons.
    //
    // All edges are modeled as spherical geodesics.  The result can be converted
    // to an area on the Earth's surface (with a worst-case error of 0.900% near
    // the poles) using the functions in s2earth.h.
    public static double GetArea(S2ShapeIndex index)
    {
        double area = 0;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null) area += S2.GetArea(shape);
        }
        return area;
    }

    // Like GetArea(), except that this method is faster and has more error.  The
    // additional error is at most 2.22e-15 steradians per vertex, which works out
    // to about 0.09 square meters per vertex on the Earth's surface.  For
    // example, a loop with 100 vertices has a maximum error of about 9 square
    // meters.  (The actual error is typically much smaller than this.)
    public static double GetApproxArea(S2ShapeIndex index)
    {
        double area = 0;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null) area += S2.GetApproxArea(shape);
        }
        return area;
    }

    // Returns the centroid of all shapes whose dimension is maximal within the
    // index, multiplied by the measure of those shapes.  For example, if the
    // index contains points and polylines, then the result is defined as the
    // centroid of the polylines multiplied by the total length of those
    // polylines.  The points would be ignored when computing the centroid.
    //
    // The measure of a given shape is defined as follows:
    //
    //  - For dimension 0 shapes, the measure is shape.num_edges().
    //  - For dimension 1 shapes, the measure is GetLength(shape).
    //  - For dimension 2 shapes, the measure is GetArea(shape).
    //
    // Note that the centroid is not unit length, so you may need to call
    // Normalize() before passing it to other S2 functions.  Note that this
    // function returns (0, 0, 0) if the index contains no geometry.
    //
    // The centroid is scaled by the total measure of the shapes for two reasons:
    // (1) it is cheaper to compute this way, and (2) this makes it easier to
    // compute the centroid of a collection of shapes (since the individual
    // centroids can simply be summed).
    public static S2Point GetCentroid(S2ShapeIndex index)
    {
        int dim = GetDimension(index);
        var centroid = S2Point.Empty;
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var shape = index.Shape(i);
            if (shape is not null && shape.Dimension() == dim)
            {
                centroid += S2.GetCentroid(shape);
            }
        }
        return centroid;
    }
}
