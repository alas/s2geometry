namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // A class representing a ShapeEdgeId together with the two endpoints of that
    // edge.  It should be passed by reference.
    public readonly record struct ShapeEdgeId(Int32 ShapeId, Int32 EdgeId)
        : IComparable<ShapeEdgeId>
    {
        public override string ToString() => $"{ShapeId}:{EdgeId}";

        public int CompareTo(ShapeEdgeId other)
        {
            if (ShapeId < other.ShapeId) return -1;
            if (ShapeId > other.ShapeId) return 1;
            if (EdgeId < other.EdgeId) return -1;
            if (EdgeId > other.EdgeId) return 1;
            return 0;
        }

        public static bool operator <(ShapeEdgeId x, ShapeEdgeId y) => x.CompareTo(y) < 0;
        public static bool operator >(ShapeEdgeId x, ShapeEdgeId y) => x.CompareTo(y) > 0;
        public static bool operator <=(ShapeEdgeId x, ShapeEdgeId y) => x.CompareTo(y) <= 0;
        public static bool operator >=(ShapeEdgeId x, ShapeEdgeId y) => x.CompareTo(y) >= 0;
    }


    public record class ShapeEdgeIdHash
    {
        public long Hash { get; private set; }

        public ShapeEdgeIdHash(Int32 ShapeId, Int32 EdgeId)
        {
            Hash = (long)(((ulong)ShapeId << 32) | (uint)EdgeId);
        }
    }
}
