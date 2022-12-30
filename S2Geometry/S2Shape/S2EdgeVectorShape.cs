// S2EdgeVectorShape is an S2Shape representing an arbitrary set of edges.  It
// is mainly used for testing, but it can also be useful if you have, say, a
// collection of polylines and don't care about memory efficiency (since this
// class would store most of the vertices twice).
//
// Note that if you already have data stored in an S2Loop, S2Polyline, or
// S2Polygon, then you would be better off using the "Shape" class defined
// within those classes (e.g., S2Loop.Shape).  Similarly, if the vertex data
// is stored in your own data structures, you can easily write your own
// subclass of S2Shape that points to the existing vertex data rather than
// copying it.

namespace S2Geometry;

public class S2EdgeVectorShape : S2Shape
{
    #region Fields, Constants

    private readonly List<(S2Point, S2Point)> edges_ = new();

    #endregion

    #region Constructors

    // Constructs an empty edge vector.
    public S2EdgeVectorShape() { }

    // Constructs an S2EdgeVectorShape from a vector of edges.
    public S2EdgeVectorShape(List<(S2Point, S2Point)> edges) => edges_ = edges;

    // Creates an S2EdgeVectorShape containing a single edge.
    public S2EdgeVectorShape(S2Point a, S2Point b) => edges_.Add((a, b));

    #endregion

    #region S2EdgeVectorShape

    // Adds an edge to the vector.
    //
    // IMPORTANT: This method should only be called *before* adding the
    // S2EdgeVectorShape to an S2ShapeIndex.  S2Shapes can only be modified by
    // removing them from the index, making changes, and adding them back again.
    public void Add(S2Point a, S2Point b) => edges_.Add((a, b));

    #endregion

    #region S2Shape

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return edges_.Count;
    }

    public sealed override Edge GetEdge(int e) => new(edges_[e].Item1, edges_[e].Item2);
    public sealed override int Dimension() => 1;
    public sealed override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
    public sealed override int NumChains() => edges_.Count;
    public sealed override Chain GetChain(int i) => new(i, 1);
    public sealed override Edge ChainEdge(int i, int j)
    {
        Debug.Assert(j == 0);
        return new Edge(edges_[i].Item1, edges_[i].Item2);
    }
    public sealed override ChainPosition GetChainPosition(int e) => new(e, 0);

    #endregion
}
