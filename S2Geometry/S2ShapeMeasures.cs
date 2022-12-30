// Defines various angle and area measures for S2Shape objects.  Unlike the
// built-in S2Polygon and S2Polyline methods, these methods allow the
// underlying data to be represented arbitrarily.

namespace S2Geometry;

public static partial class S2
{
    // For shapes of dimension 1, returns the sum of all polyline lengths on the
    // unit sphere.  Otherwise returns zero.  (See GetPerimeter for shapes of
    // dimension 2.)
    //
    // All edges are modeled as spherical geodesics.  The result can be converted
    // to a distance on the Earth's surface (with a worst-case error of 0.562%
    // near the equator) using the functions in s2earth.h.
    public static S1Angle GetLength(S2Shape shape)
    {
        if (shape.Dimension() != 1) return S1Angle.Zero;
        var length = S1Angle.Zero;
        int num_chains = shape.NumChains();
        for (int chain_id = 0; chain_id < num_chains; ++chain_id)
        {
            GetChainVertices(shape, chain_id, out var vertices);
            length += S2PolylineMeasures.GetLength(vertices);
        }
        return length;
    }

    // For shapes of dimension 2, returns the sum of all loop perimeters on the
    // unit sphere.  Otherwise returns zero.  (See GetLength for shapes of
    // dimension 1.)
    //
    // All edges are modeled as spherical geodesics.  The result can be converted
    // to a distance on the Earth's surface (with a worst-case error of 0.562%
    // near the equator) using the functions in s2earth.h.
    public static S1Angle GetPerimeter(S2Shape shape)
    {
        if (shape.Dimension() != 2) return S1Angle.Zero;
        var perimeter = S1Angle.Zero;
        int num_chains = shape.NumChains();
        for (int chain_id = 0; chain_id < num_chains; ++chain_id)
        {
            GetChainVertices(shape, chain_id, out var vertices);
            perimeter += S2.GetPerimeter(vertices);
        }
        return perimeter;
    }

    // For shapes of dimension 2, returns the area of the shape on the unit
    // sphere.  The result is between 0 and 4*Pi steradians.  Otherwise returns
    // zero.  This method has good relative accuracy for both very large and very
    // small regions.
    //
    // All edges are modeled as spherical geodesics.  The result can be converted
    // to an area on the Earth's surface (with a worst-case error of 0.900% near
    // the poles) using the functions in s2earth.h.
    public static double GetArea(S2Shape shape)
    {
        if (shape.Dimension() != 2) return 0.0;

        // Since S2Shape uses the convention that the interior of the shape is to
        // the left of all edges, in theory we could compute the area of the polygon
        // by simply adding up all the loop areas modulo 4*Pi.  The problem with
        // this approach is that polygons holes typically have areas near 4*Pi,
        // which can create large cancellation errors when computing the area of
        // small polygons with holes.  For example, a shell with an area of 4 square
        // meters (1e-13 steradians) surrounding a hole with an area of 3 square
        // meters (7.5e-14 sterians) would lose almost all of its accuracy if the
        // area of the hole was computed as 12.566370614359098.
        //
        // So instead we use S2.GetSignedArea() to ensure that all loops have areas
        // in the range [-2*Pi, 2*Pi].
        //
        // TODO(ericv): Rarely, this function returns the area of the complementary
        // region (4*Pi - area).  This can only happen when the true area is very
        // close to zero or 4*Pi and the polygon has multiple loops.  To make this
        // function completely robust requires checking whether the signed area sum is
        // ambiguous, and if so, determining the loop nesting structure.  This allows
        // the sum to be evaluated in a way that is guaranteed to have the correct
        // sign.
        double area = 0;
        double max_error = 0;
        int num_chains = shape.NumChains();
        for (int chain_id = 0; chain_id < num_chains; ++chain_id)
        {
            GetChainVertices(shape, chain_id, out var vertices);
            area += S2.GetSignedArea(vertices.ToList());
#if DEBUG
            max_error += S2.GetCurvatureMaxError(new(vertices));
#endif
        }
        // Note that S2.GetSignedArea() guarantees that the full loop (containing
        // all points on the sphere) has a very small negative area.
        Debug.Assert(Math.Abs(area) <= S2.M_4_PI + max_error);
        if (area < 0.0) area += S2.M_4_PI;
        return area;
    }

    // Like GetArea(), except that this method is faster and has more error.  The
    // additional error is at most 2.22e-15 steradians per vertex, which works out
    // to about 0.09 square meters per vertex on the Earth's surface.  For
    // example, a loop with 100 vertices has a maximum error of about 9 square
    // meters.  (The actual error is typically much smaller than this.)
    public static double GetApproxArea(S2Shape shape)
    {
        if (shape.Dimension() != 2) return 0.0;

        double area = 0;
        int num_chains = shape.NumChains();
        for (int chain_id = 0; chain_id < num_chains; ++chain_id)
        {
            GetChainVertices(shape, chain_id, out var vertices);
            area += S2.GetApproxArea(vertices.ToList()); //todo
        }
        // Special case to ensure that full polygons are handled correctly.
        if (area <= S2.M_4_PI) return area;
        return area % S2.M_4_PI;
    }

    // Returns the centroid of the shape multiplied by the measure of the shape,
    // which is defined as follows:
    //
    //  - For dimension 0 shapes, the measure is shape.num_edges().
    //  - For dimension 1 shapes, the measure is GetLength(shape).
    //  - For dimension 2 shapes, the measure is GetArea(shape).
    //
    // Note that the result is not unit length, so you may need to call
    // Normalize() before passing it to other S2 functions.
    //
    // The result is scaled by the measure defined above for two reasons: (1) it
    // is cheaper to compute this way, and (2) this makes it easier to compute the
    // centroid of a collection of shapes.  (This requires simply summing the
    // centroids of all shapes in the collection whose dimension is maximal.)
    public static S2Point GetCentroid(S2Shape shape)
    {
        var centroid = S2Point.Empty;
        S2Point[] vertices;
        int dimension = shape.Dimension();
        int num_chains = shape.NumChains();
        for (int chain_id = 0; chain_id < num_chains; ++chain_id)
        {
            switch (dimension)
            {
                case 0:
                    centroid += shape.GetEdge(chain_id).V0;
                    break;
                case 1:
                    GetChainVertices(shape, chain_id, out vertices);
                    centroid += S2PolylineMeasures.GetCentroid(vertices);
                    break;
                default:
                    GetChainVertices(shape, chain_id, out vertices);
                    centroid += S2.GetCentroid(vertices.ToList());
                    break;
            }
        }
        return centroid;
    }

    // Overwrites "vertices" with the vertices of the given edge chain of "shape".
    // If dimension == 1, the chain will have (chain.length + 1) vertices, and
    // otherwise it will have (chain.length) vertices.
    //
    // This is a low-level helper method used in the implementations of some of
    // the methods above.
    public static void GetChainVertices(S2Shape shape, int chain_id, out S2Point[] vertices)
    {
        S2Shape.Chain chain = shape.GetChain(chain_id);
        int num_vertices = chain.Length + (shape.Dimension() == 1 ? 1 : 0);
        int e = 0;
        var verts = new List<S2Point>();
        if ((num_vertices & 1) != 0)
        {
            verts.Add(shape.ChainEdge(chain_id, e++).V0);
        }
        for (; e < num_vertices; e += 2)
        {
            var edge = shape.ChainEdge(chain_id, e);
            verts.Add(edge.V0);
            verts.Add(edge.V1);
        }

        vertices = verts.ToArray();
    }
}
