using System;
namespace S2Geometry
{
    // S2LaxPolylineShape represents a polyline.  It is similar to
    // S2Polyline.Shape except that duplicate vertices are allowed, and the
    // representation is slightly more compact.
    //
    // Polylines may have any number of vertices, but note that polylines with
    // fewer than 2 vertices do not define any edges.  (To create a polyline
    // consisting of a single degenerate edge, either repeat the same vertex twice
    // or use S2LaxClosedPolylineShape defined in s2_lax_loop_shape.h.)
    public class S2LaxPolylineShape : S2Shape, IEncodeInit
    {
        public int NumVertices { get; private set; }

        private S2Point[] vertices_;

        #region Constructors

        // Constructs an empty polyline.
        public S2LaxPolylineShape() { NumVertices = 0; }

        // Constructs an S2LaxPolylineShape with the given vertices.
        public S2LaxPolylineShape(S2Point[] vertices)
        {
            Init(vertices);
        }

        // Constructs an S2LaxPolylineShape from the given S2Polyline, by copying
        // its data.
        public S2LaxPolylineShape(S2Polyline polyline)
        {
            Init(polyline);
        } 

        #endregion

        // Initializes an S2LaxPolylineShape with the given vertices.
        public void Init(S2Point[] vertices)
        {
            NumVertices = vertices.Length;
            vertices_ = (S2Point[])vertices.Clone();
        }

        // Initializes an S2LaxPolylineShape from the given S2Polyline, by copying
        // its data.
        public void Init(S2Polyline polyline)
        {
            NumVertices = polyline.NumVertices;
            vertices_ = polyline.VerticesClone;
        }

        public S2Point Vertex(int i) { return vertices_[i]; }
        public S2Point[] Vertices() { return (S2Point[])vertices_.Clone(); }

        // Appends an encoded representation of the S2LaxPolylineShape to "encoder".
        //
        // REQUIRES: "encoder" uses the defaultructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public void Encode(Encoder encoder)
        {
            Encode(encoder, CodingHint.COMPACT);
        }

        public void Encode(Encoder encoder, CodingHint hint)
        {
            EncodedS2PointVector.EncodeS2PointVector(vertices_, hint, encoder);
        }

        // Decodes an S2LaxPolylineShape, returning true on success.  (The method
        // name is chosen for compatibility with EncodedS2LaxPolylineShape below.)
        public bool Init(Decoder decoder)
        {
            var vertices = new EncodedS2PointVector();
            if (!vertices.Init(decoder)) return false;
            NumVertices = vertices.Count();
            vertices_ = new S2Point[vertices.Count()];
            for (int i = 0; i < NumVertices; ++i)
            {
                vertices_[i] = vertices[i];
            }
            return true;
        }

        // S2Shape interface:
        public sealed override int NumEdges => Math.Max(0, NumVertices - 1);
        public sealed override Edge GetEdge(int e)
        {
            Assert.True(e < NumEdges);
            return new Edge(Vertex(e), Vertex(e + 1));
        }
        public sealed override int Dimension() { return 1; }
        public sealed override ReferencePoint GetReferencePoint() { return ReferencePoint.FromContained(false); }
        public sealed override int NumChains()
        {
            return Math.Min(1, NumEdges);
        }
        public sealed override Chain GetChain(int i)
        {
            return new Chain(0, NumEdges);
        }
        public sealed override Edge ChainEdge(int i, int j)
        {
            Assert.True(i == 0);
            Assert.True(j < NumEdges);
            return new Edge(Vertex(j), Vertex(j + 1));
        }
        public sealed override ChainPosition GetChainPosition(int e)
        {
            return new S2Shape.ChainPosition(0, e);
        }
        public override S2Shape.TypeTag GetTypeTag() { return S2Shape.TypeTag.S2LaxPolylineShape; }
    }

    // Exactly like S2LaxPolylineShape, except that the vertices are kept in an
    // encoded form and are decoded only as they are accessed.  This allows for
    // very fast initialization and no additional memory use beyond the encoded
    // data.  The encoded data is not owned by this class; typically it points
    // into a large contiguous buffer that contains other encoded data as well.
    public class EncodedS2LaxPolylineShape : S2Shape
    {
        // Constructs an uninitialized object; requires Init() to be called.
        public EncodedS2LaxPolylineShape() { vertices_ = new EncodedS2PointVector(); }

        // Initializes an EncodedS2LaxPolylineShape.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder)
        {
            return vertices_.Init(decoder);
        }

        public int NumVertices => vertices_.Count();
        public S2Point Vertex(int i) { return vertices_[i]; }

        // S2Shape interface:
        public sealed override int NumEdges => Math.Max(0, NumVertices - 1);
        public sealed override Edge GetEdge(int e)
        {
            Assert.True(e < NumEdges);
            return new Edge(Vertex(e), Vertex(e + 1));
        }
        public sealed override int Dimension() { return 1; }
        public sealed override ReferencePoint GetReferencePoint() { return ReferencePoint.FromContained(false); }
        public sealed override int NumChains()
        {
            return Math.Min(1, NumEdges);
        }
        public sealed override Chain GetChain(int i)
        {
            return new Chain(0, NumEdges);
        }
        public sealed override Edge ChainEdge(int i, int j)
        {
            Assert.True(i == 0);
            Assert.True(j < NumEdges);
            return new Edge(Vertex(j), Vertex(j + 1));
        }
        public sealed override ChainPosition GetChainPosition(int e)
        {
            return new S2Shape.ChainPosition(0, e);
        }

        private readonly EncodedS2PointVector vertices_;
    }
}
