// S2LaxPolylineShape represents a polyline.  It is similar to
// S2Polyline::Shape except that adjacent vertices are allowed to be identical
// or antipodal, and the representation is slightly more compact.
//
// Polylines may have any number of vertices, but note that polylines with
// fewer than 2 vertices do not define any edges.  (To create a polyline
// consisting of a single degenerate edge, either repeat the same vertex twice
// or use S2LaxClosedPolylineShape defined in s2_lax_loop_shape.h.)

namespace S2Geometry;

public class S2LaxPolylineShape : S2Shape, IInitEncoder<S2LaxPolylineShape>
{
    public const TypeTag kTypeTag = TypeTag.S2LaxPolylineShape;

    private S2Point[]? vertices_;

    #region Constructors

    // Constructs an empty polyline.
    public S2LaxPolylineShape() { }

    // Constructs an S2LaxPolylineShape with the given vertices, by copying
    // its data.
    public S2LaxPolylineShape(IEnumerable<S2Point> vertices) => vertices_ = vertices.ToArray();

    // Constructs an S2LaxPolylineShape from the given S2Polyline.
    public S2LaxPolylineShape(S2Polyline polyline) : this(polyline.Vertices) { }

    public S2LaxPolylineShape(S2LaxPolylineShape other)
    {
        vertices_ = other.vertices_;
        other.vertices_ = null;
    }

    #endregion

    public int NumVertices() => vertices_?.Length ?? 0;
    public S2Point Vertex(int i) => vertices_[i];
    public S2Point[]? Vertices() => (S2Point[]?)(vertices_?.Clone());

    // Appends an encoded representation of the S2LaxPolylineShape to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT) =>
        EncodedS2PointVector.EncodeS2PointVector(vertices_, hint, encoder);

    // Decodes an S2LaxPolylineShape, returning true on success.  (The method
    // name is chosen for compatibility with EncodedS2LaxPolylineShape below.)
    public static (bool, S2LaxPolylineShape?) Init(Decoder decoder)
    {
        EncodedS2PointVector vertices = new();
        if (!vertices.Init(decoder)) return (false, null);

        var mumVertices = vertices.Count();
        var verticesArr = new S2Point[vertices.Count()];
        for (int i = 0; i < mumVertices; ++i)
        {
            verticesArr[i] = vertices[i];
        }
        
        return (true, new(verticesArr));
    }

    // S2Shape interface:
    public sealed override int NumEdges() => Math.Max(0, NumVertices() - 1);

    public sealed override Edge GetEdge(int e)
    {
        Debug.Assert(e < NumEdges());
        return new Edge(Vertex(e), Vertex(e + 1));
    }
    public sealed override int Dimension() => 1;
    public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public sealed override int NumChains() => Math.Min(1, NumEdges());
    public sealed override Chain GetChain(int i) => new Chain(0, NumEdges());
    public sealed override Edge ChainEdge(int i, int j)
    {
        Debug.Assert(i == 0);
        Debug.Assert(j < NumEdges());
        return new Edge(Vertex(j), Vertex(j + 1));
    }
    public sealed override ChainPosition GetChainPosition(int e) => new ChainPosition(0, e);

    // Define as enum so we don't have to declare storage.
    // TODO(user, b/210097200): Use static constexpr when C++17 is allowed
    // in opensource.
    public override TypeTag GetTypeTag() => kTypeTag;

}

// Exactly like S2LaxPolylineShape, except that the vertices are kept in an
// encoded form and are decoded only as they are accessed.  This allows for
// very fast initialization and no additional memory use beyond the encoded
// data.  The encoded data is not owned by this class; typically it points
// into a large contiguous buffer that contains other encoded data as well.
public class EncodedS2LaxPolylineShape : S2Shape, IInitEncoder<EncodedS2LaxPolylineShape>
{
    private readonly EncodedS2PointVector vertices_;

    // Constructs an uninitialized object; requires Init() to be called.
    public EncodedS2LaxPolylineShape() { vertices_ = new EncodedS2PointVector(); }

    // Initializes an EncodedS2LaxPolylineShape.
    //
    // REQUIRES: The Decoder data buffer must outlive this object.
    public static (bool, EncodedS2LaxPolylineShape?) Init(Decoder decoder)
    {
        EncodedS2LaxPolylineShape shape = new();
        if (!shape.vertices_.Init(decoder)) return (false, null);

        return (true, shape);
    }

    // Appends an encoded representation of the S2LaxPolylineShape to "encoder".
    // The coding hint is ignored, and whatever method was originally used to
    // encode the shape is preserved.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    //
    // The encoding must be identical to S2LaxPolylineShape::Encode().
    public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        vertices_.Encode(encoder);
    }

    public int NumVertices => vertices_.Count();
    public S2Point Vertex(int i) { return vertices_[i]; }

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return Math.Max(0, NumVertices - 1);
    }

    public sealed override Edge GetEdge(int e)
    {
        Debug.Assert(e < NumEdges());
        return new Edge(Vertex(e), Vertex(e + 1));
    }
    public sealed override int Dimension() { return 1; }
    public sealed override ReferencePoint GetReferencePoint() { return ReferencePoint.FromContained(false); }
    public sealed override int NumChains()
    {
        return Math.Min(1, NumEdges());
    }
    public sealed override Chain GetChain(int i)
    {
        return new Chain(0, NumEdges());
    }
    public sealed override Edge ChainEdge(int i, int j)
    {
        Debug.Assert(i == 0);
        Debug.Assert(j < NumEdges());
        return new Edge(Vertex(j), Vertex(j + 1));
    }
    public sealed override ChainPosition GetChainPosition(int e)
    {
        return new ChainPosition(0, e);
    }

    // Define as enum so we don't have to declare storage.
    // TODO(user, b/210097200): Use static constexpr when C++17 is allowed
    // in opensource.
    public override TypeTag GetTypeTag() => S2LaxPolylineShape.kTypeTag;
}
