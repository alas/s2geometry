namespace S2Geometry;

using Options = S2ContainsPointQueryOptions;

public static class S2ContainsPointQueryFactory
{
    // Returns an S2ContainsPointQuery for the given S2ShapeIndex.  Note that
    // it is efficient to return S2ContainsPointQuery objects by value.
    public static S2ContainsPointQuery<T> MakeS2ContainsPointQuery<T>(this T index, Options? options = null) where T : S2ShapeIndex
    {
        return new S2ContainsPointQuery<T>(index,
            options ?? new Options());
    }
}

// S2ContainsPointQuery determines whether one or more shapes in an
// S2ShapeIndex contain a given S2Point.  The S2ShapeIndex may contain any
// number of points, polylines, and/or polygons (possibly overlapping).
// Shape boundaries may be modeled as OPEN, SEMI_OPEN, or CLOSED (this affects
// whether or not shapes are considered to contain their vertices).
//
// Example usage:
//   var query = index.MakeS2ContainsPointQuery(S2VertexModel.CLOSED);
//   return query.Contains(point);
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
//
// However, note that if you need to do a large number of point containment
// tests, it is more efficient to re-use the S2ContainsPointQuery object
// rather than constructing a new one each time.
public class S2ContainsPointQuery<TIndex> where TIndex : S2ShapeIndex
{
    // Rather than calling this constructor, which requires specifying the
    // IndexType template argument explicitly, the preferred idiom is to call
    // MakeS2ContainsPointQuery() instead.  For example:
    //
    //   return index.MakeS2ContainsPointQuery.Contains(p);
    public S2ContainsPointQuery(TIndex index, Options? options = null)
    {
        Index = index;
        Options_ = options ?? new Options();
    }

    // Convenience constructor that accepts the S2VertexModel directly.
    public S2ContainsPointQuery(TIndex index, S2VertexModel vertex_model)
        : this(index, new Options(vertex_model)) { }

    public TIndex Index { get; }
    public Options Options_ { get; }

    // Returns true if any shape in the given index() contains the point "p"
    // under the vertex model specified (OPEN, SEMI_OPEN, or CLOSED).
    public bool Contains(S2Point p)
    {
        var (pos, found) = Index.LocatePoint(p);
        if (!found) return false;

        var icell = Index.GetIndexCell(pos);
        var cell = icell.Value.Item2;
        int num_clipped = cell.NumClipped();
        for (int s = 0; s < num_clipped; ++s)
        {
            if (ShapeContains(icell.Value.Item1, cell.Clipped(s), p)) return true;
        }
        return false;
    }

    // Returns true if the given shape contains the point "p" under the vertex
    // model specified (OPEN, SEMI_OPEN, or CLOSED).
    //
    // REQUIRES: "shape" belongs to index().
    public bool ShapeContains(S2Shape shape, S2Point p)
    {
        var (pos, found) = Index.LocatePoint(p);
        if (!found) return false;

        var icell = Index.GetIndexCell(pos);
        var cell = icell.Value.Item2;
        var clipped = cell.FindClipped(shape.Id);
        if (clipped == null) return false;
        return ShapeContains(icell.Value.Item1, clipped, p);
    }

    // Visits all shapes in the given index() that contain the given point "p",
    // terminating early if the given ShapeVisitor function returns false (in
    // which case VisitContainingShapes returns false as well).  Each shape is
    // visited at most once.
    //
    // Note that the API allows non-const access to the visited shapes.
    //
    // Also see S2ShapeIndexRegion::VisitIntersectingShapes() which allows
    // visiting all shapes in an S2ShapeIndex that intersect or contain a given
    // target S2CellId.
    //
    // ENSURES: shape != nullptr
    public delegate bool ShapeVisitor(S2Shape shape);

    public bool VisitContainingShapes(S2Point p, ShapeVisitor visitor)
    {
        // This function returns "false" only if the algorithm terminates early
        // because the "visitor" function returned false.
        var (pos, found) = Index.LocatePoint(p);
        if (!found) return true;

        var icell = Index.GetIndexCell(pos);
        var cell = icell.Value.Item2;
        int num_clipped = cell.NumClipped();
        for (int s = 0; s < num_clipped; ++s)
        {
            var clipped = cell.Clipped(s);
            if (ShapeContains(icell.Value.Item1, clipped, p) &&
                !visitor(Index.Shape(clipped.ShapeId)))
            {
                return false;
            }
        }
        return true;
    }

    // Convenience function that returns all the shapes that contain the given
    // point "p".
    public List<S2Shape> GetContainingShapes(S2Point p)
    {
        var results = new List<S2Shape>();
        VisitContainingShapes(p, (S2Shape shape) =>
        {
            results.Add(shape);
            return true;
        });
        return results;
    }

    // Visits all edges in the given index() that are incident to the point "p"
    // (i.e., "p" is one of the edge endpoints), terminating early if the given
    // EdgeVisitor function returns false (in which case VisitIncidentEdges
    // returns false as well).  Each edge is visited at most once.
    public delegate bool EdgeVisitor(S2ShapeUtil.ShapeEdge shapeEdge);

    public bool VisitIncidentEdges(S2Point p, EdgeVisitor visitor)
    {
        // This function returns "false" only if the algorithm terminates early
        // because the "visitor" function returned false.
        var (pos, found) = Index.LocatePoint(p);
        if (!found) return true;

        var icell = Index.GetIndexCell(pos);
        var cell = icell.Value.Item2;
        int num_clipped = cell.NumClipped();
        for (int s = 0; s < num_clipped; ++s)
        {
            var clipped = cell.Clipped(s);
            int num_edges = clipped.NumEdges;
            if (num_edges == 0) continue;
            var shape = Index.Shape(clipped.ShapeId);
            for (int i = 0; i < num_edges; ++i)
            {
                int edge_id = clipped.Edge(i);
                var edge = shape.GetEdge(edge_id);
                if ((edge.V0 == p || edge.V1 == p) &&
                    !visitor(new(shape.Id, edge_id, edge)))
                {
                    return false;
                }
            }
        }
        return true;
    }

    /////////////////////////// Low-Level Methods ////////////////////////////
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use.

    // Low-level helper method that returns true if the given S2ClippedShape in
    // the given cell contains the point "p".
    public bool ShapeContains(S2CellId cellId, S2ClippedShape clipped, S2Point p)
    {
        bool inside = clipped.ContainsCenter;
        int num_edges = clipped.NumEdges;
        if (num_edges > 0)
        {
            var shape = Index.Shape(clipped.ShapeId);
            if (shape.Dimension() < 2)
            {
                // Points and polylines can be ignored unless the vertex model is CLOSED.
                if (Options_.VertexModel != S2VertexModel.CLOSED) return false;

                // Otherwise, the point is contained if and only if it matches a vertex.
                for (int i = 0; i < num_edges; ++i)
                {
                    var edge = shape.GetEdge(clipped.Edge(i));
                    if (edge.V0 == p || edge.V1 == p) return true;
                }
                return false;
            }
            // Test containment by drawing a line segment from the cell center to the
            // given point and counting edge crossings.
            var crosser = new S2CopyingEdgeCrosser(cellId.ToPoint(), p);
            for (int i = 0; i < num_edges; ++i)
            {
                var edge = shape.GetEdge(clipped.Edge(i));
                int sign = crosser.CrossingSign(edge.V0, edge.V1);
                if (sign < 0) continue;
                if (sign == 0)
                {
                    // For the OPEN and CLOSED models, check whether "p" is a vertex.
                    if (Options_.VertexModel != S2VertexModel.SEMI_OPEN &&
                        (edge.V0 == p || edge.V1 == p))
                    {
                        return (Options_.VertexModel == S2VertexModel.CLOSED);
                    }
                    sign = S2.VertexCrossing(crosser.A, crosser.B, edge.V0, edge.V1) ? 1 : 0;
                }
                inside = ((inside ? 1 : 0) ^ sign) != 0;
            }
        }
        return inside;
    }
}

// This class defines the options supported by S2ContainsPointQuery.
public class S2ContainsPointQueryOptions
{
    public S2ContainsPointQueryOptions() { }

    // Convenience constructor that sets the vertex_model() option.
    public S2ContainsPointQueryOptions(S2VertexModel vertex_model)
    {
        VertexModel = vertex_model;
    }

    // Controls whether shapes are considered to contain their vertices (see
    // definitions above).  By default the SEMI_OPEN model is used.
    //
    // DEFAULT: S2VertexModel.SEMI_OPEN
    public S2VertexModel VertexModel { get; set; } = S2VertexModel.SEMI_OPEN;
}

// Defines whether shapes are considered to contain their vertices.  Note that
// these definitions differ from the ones used by S2BooleanOperation.
//
//  - In the OPEN model, no shapes contain their vertices (not even points).
//    Therefore Contains(S2Point) returns true if and only if the point is
//    in the interior of some polygon.
//
//  - In the SEMI_OPEN model, polygon point containment is defined such that
//    if several polygons tile the region around a vertex, then exactly one of
//    those polygons contains that vertex.  Points and polylines still do not
//    contain any vertices.
//
//  - In the CLOSED model, all shapes contain their vertices (including points
//    and polylines).
//
// Note that points other than vertices are never contained by polylines.
// If you need this behavior, use S2ClosestEdgeQuery::IsDistanceLess()
// with a suitable distance threshold instead.
public enum S2VertexModel : byte { OPEN, SEMI_OPEN, CLOSED }
