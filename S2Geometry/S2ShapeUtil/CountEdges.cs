namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // Returns the total number of edges in all indexed shapes.  This method takes
    // time linear in the number of shapes.
    public static int GetCountEdges(this S2ShapeIndex index)
    {
        return GetCountEdgesUpTo(index, int.MaxValue);
    }

    // Like CountEdges(), but stops once "max_edges" edges have been found (in
    // which case the current running total is returned).
    public static int GetCountEdgesUpTo(this S2ShapeIndex index, int max_edges)
    {
        int num_edges = 0;
        foreach (var shape in index)
        {
            if (shape == null) continue;
            num_edges += shape.NumEdges();
            if (num_edges >= max_edges) break;
        }
        return num_edges;
    }
}
