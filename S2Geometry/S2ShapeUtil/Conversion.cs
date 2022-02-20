// Helper functions for converting S2Shapes to typed shapes: S2Point,
// S2Polyline and S2Polygon.

namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // Converts a 0-dimensional S2Shape into a list of S2Points.
    // This method does allow an empty shape (i.e. a shape with no vertices).
    public static List<S2Point> ShapeToS2Points(S2Shape multipoint)
    {
        System.Diagnostics.Debug.Assert(multipoint.Dimension() == 0);
        List<S2Point> points = new();
        points.Capacity = multipoint.NumEdges();
        for (int i = 0; i < multipoint.NumEdges(); ++i)
        {
            points.Add(multipoint.GetEdge(i).V0);
        }
        return points;
    }

    // Converts a 1-dimensional S2Shape into an S2Polyline. Converts the first
    // chain in the shape to a vector of S2Point vertices and uses that to construct
    // the S2Polyline. Note that the input 1-dimensional S2Shape must contain at
    // most 1 chain, and that this method does not accept an empty shape (i.e. a
    // shape with no vertices).
    public static S2Polyline ShapeToS2Polyline(S2Shape line)
    {
        System.Diagnostics.Debug.Assert(line.Dimension() == 1);
        System.Diagnostics.Debug.Assert(line.NumChains() == 1);
        S2Point[] vertices;
        S2.GetChainVertices(line, 0, out vertices);
        return new S2Polyline(vertices);
    }

    // Converts a 2-dimensional S2Shape into an S2Polygon. Each chain is converted
    // to an S2Loop and the vector of loops is used to construct the S2Polygon.
    public static S2Polygon ShapeToS2Polygon(S2Shape poly)
    {
        if (poly.IsFull())
        {
            return new S2Polygon(S2Loop.kFull);
        }
        System.Diagnostics.Debug.Assert(poly.Dimension() == 2);
        List<S2Loop> loops = new();
        for (int i = 0; i < poly.NumChains(); ++i)
        {
            S2.GetChainVertices(poly, i, out var vertices);
            loops.Add(new S2Loop(vertices));
        }
        var output_poly = new S2Polygon(loops, initOriented: true);

        return output_poly;
    }
}
