namespace S2Geometry;

// S2PointVectorShape is an S2Shape representing a set of S2Points. Each point
// is reprsented as a degenerate edge with the same starting and ending
// vertices.
//
// This class is useful for adding a collection of points to an S2ShapeIndex.
public class S2PointVectorShape : S2Shape, IInitEncoder<S2PointVectorShape>
{
    public const TypeTag kTypeTag = TypeTag.S2PointVectorShape;

    // Constructs an empty point vector.
    public S2PointVectorShape() { }

    // Constructs an S2PointVectorShape from a vector of points.
    public S2PointVectorShape(S2Point[] points)
    {
        points_ = points;
    }

    public int NumPoints => points_.Length;
    public S2Point Point(int i) { return points_[i]; }

    // Appends an encoded representation of the S2PointVectorShape to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        EncodedS2PointVector.EncodeS2PointVector(points_, hint, encoder);
    }

    // Decodes an S2PointVectorShape, returning true on success.  (The method
    // name is chosen for compatibility with EncodedS2PointVectorShape below.)
    public static (bool, S2PointVectorShape?) Init(Decoder decoder)
    {
        EncodedS2PointVector points = new();
        if (!points.Init(decoder)) return (false, null);

        S2PointVectorShape shape = new() { points_ = points.Decode() };
        return (true, shape);
    }

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return NumPoints;
    }

    public sealed override Edge GetEdge(int e) => new(points_[e], points_[e]);
    public sealed override int Dimension() { return 0; }
    public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public sealed override int NumChains() { return NumPoints; }
    public sealed override Chain GetChain(int i) { return new Chain(i, 1); }
    public sealed override Edge ChainEdge(int i, int j)
    { Assert.True(j == 0); return new Edge(points_[i], points_[i]); }
    public sealed override ChainPosition GetChainPosition(int e) => new(e, 0);
    public override TypeTag GetTypeTag() => kTypeTag;

    private S2Point[] points_;
}

// Exactly like S2PointVectorShape, except that the points are kept in an
// encoded form and are decoded only as they are accessed.  This allows for
// very fast initialization and no additional memory use beyond the encoded
// data.  The encoded data is not owned by this class; typically it points
// into a large contiguous buffer that contains other encoded data as well.
public class EncodedS2PointVectorShape : S2Shape, IEncoder, IInitEncoder<EncodedS2PointVectorShape>
{
    // Constructs an uninitialized object; requires Init() to be called.
    public EncodedS2PointVectorShape() => points_ = new();

    // Initializes an EncodedS2PointVectorShape.
    //
    // REQUIRES: The Decoder data buffer must outlive this object.
    public static (bool, EncodedS2PointVectorShape?) Init(Decoder decoder)
    {
        EncodedS2PointVectorShape shape = new();
        if (!shape.points_.Init(decoder)) return (false, null);

        return (true, shape);
    }

    // Appends an encoded representation of the S2LaxPolygonShape to "encoder".
    // The coding hint is ignored, and whatever method was originally used to
    // encode the shape is preserved.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint)
    {
        points_.Encode(encoder);
    }
    public int NumPoints => points_.Count();
    public S2Point Point(int i) { return points_[i]; }

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return NumPoints;
    }

    public sealed override Edge GetEdge(int e) => new(points_[e], points_[e]);
    public sealed override int Dimension() { return 0; }
    public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public sealed override int NumChains() { return NumPoints; }
    public sealed override Chain GetChain(int i) { return new Chain(i, 1); }
    public sealed override Edge ChainEdge(int i, int j)
    { Assert.True(j == 0); return new Edge(points_[i], points_[i]); }
    public sealed override ChainPosition GetChainPosition(int e) => new(e, 0);
    public override TypeTag GetTypeTag() => S2PointVectorShape.kTypeTag;

    private readonly EncodedS2PointVector points_;
}
