// S2PointVectorShape is an S2Shape representing a set of S2Points. Each point
// is represented as a degenerate edge with the same starting and ending
// vertices.
//
// This class is useful for adding a collection of points to an S2ShapeIndex.

namespace S2Geometry;

public class S2PointVectorShape(S2Point[] points) : S2Shape, IInitEncoder<S2PointVectorShape>, IEquatable<S2PointVectorShape>
{
    #region Fields, Constants

    public const TypeTag kTypeTag = TypeTag.S2PointVectorShape;

    private readonly S2Point[] Points = points;

    #endregion

    #region Constructors

    // Constructs an empty point vector.
    public S2PointVectorShape() : this([]) { }

    #endregion

    #region S2PointVectorShape

    public int NumPoints => Points.Length;
    public S2Point Point(int i) => Points[i];

    #endregion

    #region IInitEncoder

    // Decodes an S2PointVectorShape, returning true on success.  (The method
    // name is chosen for compatibility with EncodedS2PointVectorShape below.)
    public static (bool, S2PointVectorShape?) Init(Decoder decoder)
    {
        var (success, points) = EncodedS2PointVector.Init(decoder);
        if (!success) return (false, null);

        S2PointVectorShape shape = new(points!.Decode());
        return (true, shape);
    }

    // Appends an encoded representation of the S2PointVectorShape to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT) =>
        EncodedS2PointVector.EncodeS2PointVector(Points, hint, encoder);

    #endregion

    #region S2Shape interface.

    // Returns the number of points.
    public sealed override int NumEdges() => NumPoints;

    // Returns a point represented as a degenerate edge.
    public sealed override Edge GetEdge(int e) => new(Points[e], Points[e]);
    public sealed override int Dimension() => 0;
    public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public sealed override int NumChains() => NumPoints;
    public sealed override Chain GetChain(int i) => new(i, 1);
    public sealed override Edge ChainEdge(int i, int j)
    { MyDebug.Assert(j == 0); return new Edge(Points[i], Points[i]); }
    public sealed override ChainPosition GetChainPosition(int e) => new(e, 0);

    // Define as enum so we don't have to declare storage.
    // TODO(user, b/210097200): Use static constexpr when C++17 is allowed
    // in opensource.
    public override TypeTag GetTypeTag() => kTypeTag;

    #endregion

    #region IEquatable

    public override bool Equals(object? other) => other is S2PointVectorShape shape && Equals(shape);

    public bool Equals(S2PointVectorShape? other) =>
        other is not null && other.Id == Id && Enumerable.SequenceEqual(other.Points, Points);

    public override int GetHashCode() =>
        HashCode.Combine(Id.GetHashCode(), Points.GetHashCode());

    #endregion
}

// Exactly like S2PointVectorShape, except that the points are kept in an
// encoded form and are decoded only as they are accessed.  This allows for
// very fast initialization and no additional memory use beyond the encoded
// data.  The encoded data is not owned by this class; typically it points
// into a large contiguous buffer that contains other encoded data as well.
public class EncodedS2PointVectorShape(EncodedS2PointVector points) : S2Shape, IInitEncoder<EncodedS2PointVectorShape>
{
    #region Fields, Constants

    public EncodedS2PointVector Points { private get; init; } = points;

    #endregion

    #region IInitEncoder

    // Initializes an EncodedS2PointVectorShape.
    //
    // REQUIRES: The Decoder data buffer must outlive this object.
    public static (bool, EncodedS2PointVectorShape?) Init(Decoder decoder)
    {
        var (success, points_) = EncodedS2PointVector.Init(decoder);
        if (!success) return (false, null);

        EncodedS2PointVectorShape shape = new(points_!);
        return (true, shape);
    }

    // Appends an encoded representation of the S2LaxPolygonShape to "encoder".
    // The coding hint is ignored, and whatever method was originally used to
    // encode the shape is preserved.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        Points.Encode(encoder);
    }

    #endregion

    #region EncodedS2PointVectorShape

    public int NumPoints => Points.Count();
    public S2Point Point(int i) => Points[i];

    #endregion

    #region S2Shape interface:

    public sealed override int NumEdges() => NumPoints;

    public sealed override Edge GetEdge(int e) => new(Points[e], Points[e]);
    public sealed override int Dimension() => 0;
    public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public sealed override int NumChains() => NumPoints;
    public sealed override Chain GetChain(int i) => new(i, 1);
    public sealed override Edge ChainEdge(int i, int j)
    { MyDebug.Assert(j == 0); return new Edge(Points[i], Points[i]); }
    public sealed override ChainPosition GetChainPosition(int e) => new(e, 0);
    public override TypeTag GetTypeTag() => S2PointVectorShape.kTypeTag;

    #endregion
}
