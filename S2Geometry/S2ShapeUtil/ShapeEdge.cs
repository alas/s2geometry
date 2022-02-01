namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // A class representing a ShapeEdgeId together with the two endpoints of that
    // edge.  It should be passed by reference.
    public record class ShapeEdge(Edge Id, S2Point V0, S2Point V1)
    {
        public ShapeEdge(S2Shape shape, Int32 edge_id)
            : this(shape.Id, edge_id, shape.GetEdge(edge_id)) { }
        public ShapeEdge(Int32 shape_id, Int32 edge_id, S2Shape.Edge edge)
            : this(new(shape_id, edge_id), edge.V0, edge.V1) { }
    }
}
