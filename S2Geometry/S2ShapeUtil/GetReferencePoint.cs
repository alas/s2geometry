namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // This is a helper function for implementing S2Shape.GetReferencePoint().
    //
    // Given a shape consisting of closed polygonal loops, the interior of the
    // shape is defined as the region to the left of all edges (which must be
    // oriented consistently).  This function then chooses an arbitrary point and
    // returns true if that point is contained by the shape.
    //
    // Unlike S2Loop and S2Polygon, this method allows duplicate vertices and
    // edges, which requires some extra care with definitions.  The rule that we
    // apply is that an edge and its reverse edge "cancel" each other: the result
    // is the same as if that edge pair were not present.  Therefore shapes that
    // consist only of degenerate loop(s) are either empty or full; by convention,
    // the shape is considered full if and only if it contains an empty loop (see
    // S2LaxPolygonShape for details).
    //
    // Determining whether a loop on the sphere contains a point is harder than
    // the corresponding problem in 2D plane geometry.  It cannot be implemented
    // just by counting edge crossings because there is no such thing as a "point
    // at infinity" that is guaranteed to be outside the loop.
    public static S2Shape.ReferencePoint GetReferencePoint(this S2Shape shape)
    {
        MyDebug.Assert(shape.Dimension() == 2);
        if (shape.NumEdges() == 0)
        {
            // A shape with no edges is defined to be full if and only if it
            // contains at least one chain.
            return S2Shape.ReferencePoint.FromContained(shape.NumChains() > 0);
        }
        // Define a "matched" edge as one that can be paired with a corresponding
        // reversed edge.  Define a vertex as "balanced" if all of its edges are
        // matched. In order to determine containment, we must find an unbalanced
        // vertex.  Often every vertex is unbalanced, so we start by trying an
        // arbitrary vertex.
        var edge = shape.GetEdge(0);
        if (GetReferencePointAtVertex(shape, edge.V0, out var result))
        {
            return result;
        }
        // That didn't work, so now we do some extra work to find an unbalanced
        // vertex (if any).  Essentially we gather a list of edges and a list of
        // reversed edges, and then sort them.  The first edge that appears in one
        // list but not the other is guaranteed to be unmatched.
        int n = shape.NumEdges();
        var edges = new S2Shape.Edge[n];
        var rev_edges = new S2Shape.Edge[n];
        for (int i = 0; i < n; ++i)
        {
            var edge2 = shape.GetEdge(i);
            edges[i] = edge2;
            rev_edges[i] = new S2Shape.Edge(edge2.V1, edge2.V0);
        }
        Array.Sort(edges);
        Array.Sort(rev_edges);
        for (int i = 0; i < n; ++i)
        {
            if (edges[i] < rev_edges[i])
            {  // edges[i] is unmatched
                MyDebug.Assert(GetReferencePointAtVertex(shape, edges[i].V0, out result));
                return result;
            }
            if (rev_edges[i] < edges[i])
            {  // rev_edges[i] is unmatched
                MyDebug.Assert(GetReferencePointAtVertex(shape, rev_edges[i].V0, out result));
                return result;
            }
        }
        // All vertices are balanced, so this polygon is either empty or full except
        // for degeneracies.  By convention it is defined to be full if it contains
        // any chain with no edges.
        for (int i = 0; i < shape.NumChains(); ++i)
        {
            if (shape.GetChain(i).Length == 0) return S2Shape.ReferencePoint.FromContained(true);
        }
        return S2Shape.ReferencePoint.FromContained(false);
    }

    // This is a helper function for GetReferencePoint() below.
    //
    // If the given vertex "vtest" is unbalanced (see definition below), sets
    // "result" to a ReferencePoint indicating whther "vtest" is contained and
    // returns true.  Otherwise returns false.
    private static bool GetReferencePointAtVertex(S2Shape shape, S2Point vtest, out S2Shape.ReferencePoint result)
    {
        // Let P be an unbalanced vertex.  Vertex P is defined to be inside the
        // region if the region contains a particular direction vector starting from
        // P, namely the direction S2::RefDir(P).  This can be calculated using
        // S2ContainsVertexQuery.
        var contains_query = new S2ContainsVertexQuery(vtest);
        int n = shape.NumEdges();
        for (int e = 0; e < n; ++e)
        {
            var edge = shape.GetEdge(e);
            if (edge.V0 == vtest) contains_query.AddEdge(edge.V1, 1);
            if (edge.V1 == vtest) contains_query.AddEdge(edge.V0, -1);
        }
        int contains_sign = contains_query.ContainsSign();
        if (contains_sign == 0)
        {
            result = new S2Shape.ReferencePoint(S2Point.Empty, false);
            return false;  // There are no unmatched edges incident to this vertex.
        }
        result = new S2Shape.ReferencePoint(vtest, contains_sign > 0);
        return true;
    }
}
