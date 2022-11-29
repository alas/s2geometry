namespace S2Geometry;

/// <summary>
/// SourceId identifies an edge from one of the two input S2ShapeIndexes.
/// It consists of a region id (0 or 1), a shape id within that region's
/// S2ShapeIndex, and an edge id within that shape.
///
/// Note(Alas): extracted from S2BooleanOperation
/// </summary>
public readonly record struct SourceId(int RegionId, int ShapeId, int EdgeId) : IComparable<SourceId>
{
    #region Fields, Constants

    public static readonly SourceId Zero = new(0, 0, -1);

    #endregion

    #region Constructors

    public SourceId(Int32 special_edge_id)
        : this(0, 0, special_edge_id)
    {
    }

    #endregion

    #region IComparable

    public int CompareTo(SourceId other)
    {
        if (RegionId < other.RegionId) return -1;
        if (RegionId > other.RegionId) return 1;
        if (ShapeId < other.ShapeId) return -1;
        if (ShapeId > other.ShapeId) return 1;
        if (EdgeId < other.EdgeId) return -1;
        if (EdgeId > other.EdgeId) return 1;
        return 0;
    }

    public static bool operator <(SourceId x, SourceId y) => x.CompareTo(y) < 0;
    public static bool operator >(SourceId x, SourceId y) => x.CompareTo(y) > 0;
    public static bool operator <=(SourceId x, SourceId y) => x.CompareTo(y) <= 0;
    public static bool operator >=(SourceId x, SourceId y) => x.CompareTo(y) >= 0;

    #endregion
}
