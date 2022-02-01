// An S2Shape that simply some other shape.  This is useful for adding
// an existing S2Shape to a new S2ShapeIndex without needing to copy
// its underlying data.
//
// Also see s2shapeutil::WrappedShapeFactory in s2shapeutil_coding.h, which
// is useful for testing S2ShapeIndex coding.

namespace S2Geometry;

public class S2WrappedShape : S2Shape
{
    private readonly S2Shape shape_;

    public S2WrappedShape(S2Shape shape) => shape_ = shape;

    // S2Shape interface:
    public sealed override int NumEdges() => shape_.NumEdges();
    public sealed override Edge GetEdge(int e) => shape_.GetEdge(e);
    public sealed override int Dimension() => shape_.Dimension();
    public sealed override ReferencePoint GetReferencePoint() => shape_.GetReferencePoint();
    public sealed override int NumChains() => shape_.NumChains();
    public sealed override Chain GetChain(int i) => shape_.GetChain(i);
    public sealed override Edge ChainEdge(int i, int j) => shape_.ChainEdge(i, j);
    public sealed override ChainPosition GetChainPosition(int e) => shape_.GetChainPosition(e);
}
