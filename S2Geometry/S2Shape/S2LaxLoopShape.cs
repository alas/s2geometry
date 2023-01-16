// S2LaxLoopShape represents a closed loop of edges surrounding an interior
// region.  It is similar to S2Loop.Shape except that this class allows
// duplicate vertices and edges.  Loops may have any number of vertices,
// including 0, 1, or 2.  (A one-vertex loop defines a degenerate edge
// consisting of a single point.)
//
// Note that S2LaxLoopShape is faster to initialize and more compact than
// S2Loop.Shape, but does not support the same operations as S2Loop.

namespace S2Geometry;

public class S2LaxLoopShape : S2Shape, IEquatable<S2LaxLoopShape>
{
    #region Fields, Constants

    private readonly S2Point[] vertices_;
    public int NumVertices { get; private set; }

    #endregion

    #region Constructors

    // Constructs an empty loop.
    public S2LaxLoopShape() { NumVertices = 0; vertices_ = Array.Empty<S2Point>(); }

    // Constructs an S2LaxLoopShape with the given vertices.
    //
    // Initializes an S2LaxLoopShape with the given vertices.
    public S2LaxLoopShape(S2Point[] vertices)
    {
        NumVertices = vertices.Length;
        vertices_ = vertices;
    }

    // Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
    //
    // Initializes an S2LaxLoopShape from the given S2Loop, by copying its data.
    //
    // REQUIRES: !loop.IsFull
    //           [Use S2LaxPolygonShape if you need to represent a full loop.]
    public S2LaxLoopShape(S2Loop loop) 
    {
        MyDebug.Assert(!loop.IsFull()); // Full loops not supported; use S2LaxPolygonShape
        if (loop.IsEmpty())
        {
            NumVertices = 0;
            vertices_ = Array.Empty<S2Point>();
        }
        else
        {
            NumVertices = loop.NumVertices;
            vertices_ = loop.CloneVertices();
        }
    }

    #endregion

    #region S2LaxLoopShape

    public S2Point Vertex(int i) => vertices_[i];

    #endregion

    #region S2Shape

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return NumVertices;
    }

    public sealed override Edge GetEdge(int e0)
    {
        MyDebug.Assert(e0 < NumEdges());
        int e1 = e0 + 1;
        if (e1 == NumVertices) e1 = 0;
        return new Edge(Vertex(e0), Vertex(e1));
    }
    // Not final; overridden by S2LaxClosedPolylineShape.
    public override int Dimension() { return 2; }
    // Not final; overridden by S2LaxClosedPolylineShape.
    public override ReferencePoint GetReferencePoint() => S2ShapeUtil.GetReferencePoint(this);
    public sealed override int NumChains() => Math.Min(1, NumVertices);
    public sealed override Chain GetChain(int i) => new(0, NumVertices);
    public sealed override Edge ChainEdge(int i, int j)
    {
        MyDebug.Assert(i == 0);
        MyDebug.Assert(j < NumEdges());
        int k = (j + 1 == NumVertices) ? 0 : j + 1;
        return new Edge(Vertex(j), Vertex(k));
    }
    public sealed override ChainPosition GetChainPosition(int e) => new(0, e);

    #endregion

    #region IEquatable

    public override bool Equals(object? other) => other is S2LaxLoopShape shape && Equals(shape);

    public bool Equals(S2LaxLoopShape? other) =>
        other is not null && other.Id == Id && other.NumVertices == NumVertices
        && Enumerable.SequenceEqual(other.vertices_, vertices_);

    public override int GetHashCode() => vertices_.GetHashCode();

    #endregion
}

// S2LaxClosedPolylineShape is like S2LaxPolylineShape except that the last
// vertex is implicitly joined to the first.  It is also like S2LaxLoopShape
// except that it does not have an interior (which makes it more efficient to
// index).
public class S2LaxClosedPolylineShape : S2LaxLoopShape
{
    #region Constructors

    public S2LaxClosedPolylineShape() : base() { }

    public S2LaxClosedPolylineShape(S2Point[] vertices) : base(vertices) { }

    public S2LaxClosedPolylineShape(S2Loop loop) : base(loop) { }

    #endregion

    // See S2LaxLoopShape for constructors.
    public sealed override int Dimension() => 1;
    public sealed override ReferencePoint GetReferencePoint() =>
        ReferencePoint.FromContained(false);
}

// S2VertexIdLaxLoopShape is just like S2LaxLoopShape, except that vertices are
// specified as indices into a vertex array.  This representation can be more
// compact when many loops are arranged in a mesh structure.
public class S2VertexIdLaxLoopShape : S2Shape, IEquatable<S2VertexIdLaxLoopShape>
{
    #region Fields, Constants

    private readonly Int32[] vertex_ids_;
    private readonly S2Point[] vertex_array_;

    #endregion

    #region Constructors

    // Constructs an empty loop.
    public S2VertexIdLaxLoopShape()
    {
        NumVertices = 0;
        vertex_ids_ = Array.Empty<Int32>();
        vertex_array_ = Array.Empty<S2Point>();
    }

    // Constructs the shape from the given vertex array and indices.
    // "vertex_ids" is a vector of indices into "vertex_array".
    //
    // ENSURES:  loop.vertex(i) == (*vertex_array)[vertex_ids[i]]
    // REQUIRES: "vertex_array" persists for the lifetime of this object.
    //
    // Initializes the shape from the given vertex array and indices.
    // "vertex_ids" is a vector of indices into "vertex_array".
    public S2VertexIdLaxLoopShape(Int32[] vertex_ids, S2Point[]? vertex_array)
    {
        NumVertices = vertex_ids.Length;
        vertex_ids_ = (Int32[])vertex_ids.Clone();
        vertex_array_ = vertex_array ?? Array.Empty<S2Point>();
    }

    #endregion

    #region S2VertexIdLaxLoopShape

    // Returns the number of vertices in the loop.
    public int NumVertices { get; private set; }
    public Int32 VertexId(int i) => vertex_ids_[i];
    public S2Point Vertex(int i) => vertex_array_[VertexId(i)];

    #endregion

    #region S2Shape

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return NumVertices;
    }

    public sealed override Edge GetEdge(int e0)
    {
        MyDebug.Assert(e0 < NumEdges());
        int e1 = e0 + 1;
        if (e1 == NumVertices) e1 = 0;
        return new Edge(Vertex(e0), Vertex(e1));
    }
    public sealed override int Dimension() => 2;
    public sealed override ReferencePoint GetReferencePoint()
    {
        // GetReferencePoint interprets a loop with no vertices as "full".
        if (NumVertices == 0) return ReferencePoint.FromContained(false);
        return S2ShapeUtil.GetReferencePoint(this);
    }
    public sealed override int NumChains() => Math.Min(1, NumVertices);
    public sealed override Chain GetChain(int i) => new(0, NumVertices);
    public sealed override Edge ChainEdge(int i, int j)
    {
        MyDebug.Assert(i == 0);
        MyDebug.Assert(j < NumEdges());
        int k = (j + 1 == NumVertices) ? 0 : j + 1;
        return new Edge(Vertex(j), Vertex(k));
    }
    public sealed override ChainPosition GetChainPosition(int e) => new(0, e);

    #endregion

    #region IEquatable

    public override bool Equals(object? other) => other is S2VertexIdLaxLoopShape shape && Equals(shape);

    public bool Equals(S2VertexIdLaxLoopShape? other) =>
        other is not null && other.Id == Id
        && Enumerable.SequenceEqual(other.vertex_ids_, vertex_ids_)
        && Enumerable.SequenceEqual(other.vertex_array_, vertex_array_);

    public override int GetHashCode() =>
        HashCode.Combine(Id.GetHashCode(), vertex_ids_.GetHashCode(), vertex_array_.GetHashCode());

    #endregion
}

