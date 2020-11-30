using System;

namespace S2Geometry
{
    // SourceId identifies an edge from one of the two input S2ShapeIndexes.
    // It consists of a region id (0 or 1), a shape id within that region's
    // S2ShapeIndex, and an edge id within that shape.
    //
    // Note(Alas): extracted from S2BooleanOperation
    public readonly struct SourceId : IEquatable<SourceId>, IComparable<SourceId>
    {
        #region Fields, Constants

        public readonly int RegionId;
        public readonly int ShapeId;
        public readonly int EdgeId;

        public static readonly SourceId Zero = new SourceId(0, 0, -1);

        #endregion

        #region Constructors

        public SourceId(int region_id, Int32 shape_id, Int32 edge_id)
        {
            RegionId = region_id; ShapeId = shape_id; EdgeId = edge_id;
        }
        public SourceId(Int32 special_edge_id)
        {
            RegionId = 0; ShapeId = 0; EdgeId = special_edge_id;
        }

        #endregion

        #region IEquatable

        public override int GetHashCode() => HashCode.Combine(RegionId, ShapeId, EdgeId);
        public bool Equals(SourceId other) =>
            RegionId == other.RegionId &&
            ShapeId == other.ShapeId &&
            EdgeId == other.EdgeId;
        public override bool Equals(object other) => other is SourceId sid && Equals(sid);
        public static bool operator ==(SourceId x, SourceId y) => Equals(x, y);
        public static bool operator !=(SourceId x, SourceId y) => !Equals(x, y);

        #endregion

        #region IComparable

        public int CompareTo(SourceId other)
        {
            if (RegionId.CompareTo(other.RegionId) != 0)
                return RegionId.CompareTo(other.RegionId);

            if (ShapeId.CompareTo(other.ShapeId) != 0)
                return ShapeId.CompareTo(other.ShapeId);

            return EdgeId.CompareTo(other.EdgeId);
        }
        public static bool operator <(SourceId x, SourceId y) => x.CompareTo(y) < 0;
        public static bool operator >(SourceId x, SourceId y) => x.CompareTo(y) > 0;

        #endregion
    }
}
