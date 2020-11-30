using System;
namespace S2Geometry
{
    // S2PointVectorShape is an S2Shape representing a set of S2Points. Each point
    // is reprsented as a degenerate edge with the same starting and ending
    // vertices.
    //
    // This class is useful for adding a collection of points to an S2ShapeIndex.
    public class S2PointVectorShape : S2Shape, IEncodeInit
    {
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
        // REQUIRES: "encoder" uses the defaultructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public void Encode(Encoder encoder)
        {
            Encode(encoder, CodingHint.COMPACT);
        }
        public void Encode(Encoder encoder, CodingHint hint)
        {
            EncodedS2PointVector.EncodeS2PointVector(points_, hint, encoder);
        }

        // Decodes an S2PointVectorShape, returning true on success.  (The method
        // name is chosen for compatibility with EncodedS2PointVectorShape below.)
        public bool Init(Decoder decoder)
        {
            var points = new EncodedS2PointVector();
            if (!points.Init(decoder)) return false;
            points_ = points.Decode();
            return true;
        }

        // S2Shape interface:
        public sealed override int NumEdges => NumPoints;
        public sealed override Edge GetEdge(int e) => new Edge(points_[e], points_[e]);
        public sealed override int Dimension() { return 0; }
        public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
        public sealed override int NumChains() { return NumPoints; }
        public sealed override Chain GetChain(int i) { return new Chain(i, 1); }
        public sealed override Edge ChainEdge(int i, int j)
        { Assert.True(j == 0); return new Edge(points_[i], points_[i]); }
        public sealed override ChainPosition GetChainPosition(int e) => new ChainPosition(e, 0);
        public override TypeTag GetTypeTag() => TypeTag.S2PointVectorShape;

        private S2Point[] points_;
    }

    // Exactly like S2PointVectorShape, except that the points are kept in an
    // encoded form and are decoded only as they are accessed.  This allows for
    // very fast initialization and no additional memory use beyond the encoded
    // data.  The encoded data is not owned by this class; typically it points
    // into a large contiguous buffer that contains other encoded data as well.
    public class EncodedS2PointVectorShape : S2Shape
    {
        // Constructs an uninitialized object; requires Init() to be called.
        public EncodedS2PointVectorShape() { points_ = new EncodedS2PointVector(); }

        // Initializes an EncodedS2PointVectorShape.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder) { return points_.Init(decoder); }

        public int NumPoints => points_.Count();
        public S2Point Point(int i) { return points_[i]; }

        // S2Shape interface:
        public sealed override int NumEdges => NumPoints;
        public sealed override Edge GetEdge(int e) => new Edge(points_[e], points_[e]);
        public sealed override int Dimension() { return 0; }
        public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
        public sealed override int NumChains() { return NumPoints; }
        public sealed override Chain GetChain(int i) { return new Chain(i, 1); }
        public sealed override Edge ChainEdge(int i, int j)
        { Assert.True(j == 0); return new Edge(points_[i], points_[i]); }
        public sealed override ChainPosition GetChainPosition(int e) => new ChainPosition(e, 0);

        private readonly EncodedS2PointVector points_;
    }
}
