using System;

namespace S2Geometry.S2ShapeUtil
{
    // A class representing a ShapeEdgeId together with the two endpoints of that
    // edge.  It should be passed by reference.
    public record ShapeEdge(Edge Id, S2Point V0, S2Point V1)
    {
        public ShapeEdge(S2Shape shape, Int32 edge_id)
            : this(shape.Id, edge_id, shape.GetEdge(edge_id)) { }
        public ShapeEdge(Int32 shape_id, Int32 edge_id, S2Shape.Edge edge)
            : this(new(shape_id, edge_id), edge.V0, edge.V1) { }
    }

    public readonly struct Edge : IEquatable<Edge>, IComparable<Edge>
    {
        public readonly int ShapeId;
        public readonly int EdgeId;

        public Edge(int shapeId, int edgeId) { ShapeId = shapeId; EdgeId = edgeId; }
        public void Deconstruct(out int shapeId, out int edgeId) => (shapeId, edgeId) = (ShapeId, EdgeId);

        #region IEquatable

        public override bool Equals(object obj) => obj is Edge edge && Equals(edge);
        public bool Equals(Edge edge) => ShapeId == edge.ShapeId && EdgeId == edge.EdgeId;
        public override int GetHashCode() => HashCode.Combine(ShapeId, EdgeId);

        public static bool operator ==(Edge x, Edge y) => Equals(x, y);
        public static bool operator !=(Edge x, Edge y) => !Equals(x, y);

        #endregion

        #region IComparable

        public int CompareTo(Edge other)
        {
            if (ShapeId.CompareTo(other.ShapeId) != 0)
                return ShapeId.CompareTo(other.ShapeId);

            return EdgeId.CompareTo(other.EdgeId);
        }

        #endregion
    }

    public readonly struct EdgePair
    {
        public readonly Edge Item1 { get; init; }
        public readonly Edge Item2 { get; init; }
        public EdgePair(Edge item1, Edge item2) { Item1 = item1; Item2 = item2; }
        public override string ToString() => $"({Item1},{Item2})";
    }
}
