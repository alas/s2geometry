// A 32-bit tag that can be used to identify the type of an encoded S2Shape.
// All encodable types have a non-zero type tag.  The tag associated with a
// given shape type can be accessed as Shape.kTypeTag, while the tag
// associated with a given object can be accessed as shape.type_tag().
//
// Type tags in the range 0..8191 are reserved for use by the S2 library.
//
// The purpose of S2Shape is to represent polygonal geometry in a flexible
// way.  It is organized as a collection of edges that optionally defines an
// interior.  All geometry represented by an S2Shape must have the same
// dimension, which means that an S2Shape can represent either a set of
// points, a set of polylines, or a set of polygons.
//
// S2Shape is defined as an abstract base class in order to give clients
// control over the underlying data representation.  Sometimes an S2Shape does
// not have any data of its own, but instead "wraps" some other class.  There
// are various useful subtypes defined in *_shape.h, and some S2 classes also
// have a nested "Shape" class (e.g., S2Polygon.Shape).  It is easy for
// clients to implement their own subtypes, since the interface is minimal.
//
// S2Shape operations are typically defined on S2ShapeIndex objects rather
// than individual shapes.  An S2ShapeIndex is simply a collection of
// S2Shapes, possibly of different dimensions (e.g. 10 points and 3 polygons),
// organized into a data structure for efficient edge access.
//
// The edges of an S2Shape are identified by a contiguous range of "edge ids"
// starting at 0.  The edges are further subdivided into "chains", where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// For example, an S2Shape representing two polylines AB and CDE would have
// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
// Similarly, an S2Shape representing 5 points would have 5 chains consisting
// of one edge each.
//
// S2Shape has methods that allow edges to be accessed either using the global
// numbering (edge id) or within a particular chain.  The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see S2BooleanOperation).

namespace S2Geometry;

public abstract class S2Shape
{
    #region Fields, Constants

    // Assigned by S2ShapeIndex when the shape is added.
    // A unique id assigned to this shape by S2ShapeIndex.  Shape ids are
    // assigned sequentially starting from 0 in the order shapes are added.
    public int Id { get; private set; } = -1;
    public void SetId(int id) { Id = id; }

    #endregion

    #region S2Shape

    // Returns the number of edges in this shape.  Edges have ids ranging from 0
    // to num_edges() - 1.
    public abstract int NumEdges();

    // Returns the endpoints of the given edge id.
    //
    // REQUIRES: 0 <= id < num_edges()
    public abstract Edge GetEdge(int edge_id);

    // Returns the dimension of the geometry represented by this shape.
    //
    //  0 - Point geometry.  Each point is represented as a degenerate edge.
    //
    //  1 - Polyline geometry.  Polyline edges may be degenerate.  A shape may
    //      represent any number of polylines.  Polylines edges may intersect.
    //
    //  2 - Polygon geometry.  Edges should be oriented such that the polygon
    //      interior is always on the left.  In theory the edges may be returned
    //      in any order, but typically the edges are organized as a collection
    //      of edge chains where each chain represents one polygon loop.
    //      Polygons may have degeneracies (e.g., degenerate edges or sibling
    //      pairs consisting of an edge and its corresponding reversed edge).
    //      A polygon loop may also be full (containing all points on the
    //      sphere); by convention this is represented as a chain with no edges.
    //      (See S2LaxPolygonShape for details.)
    //
    // Note that this method allows degenerate geometry of different dimensions
    // to be distinguished, e.g. it allows a point to be distinguished from a
    // polyline or polygon that has been simplified to a single point.
    public abstract int Dimension();

    // Returns true if the shape contains no points.  (Note that the full
    // polygon is represented as a chain with zero edges.)
    public bool IsEmpty() => NumEdges() == 0 && (Dimension() < 2 || NumChains() == 0);

    // Returns true if the shape contains all points on the sphere.
    public bool IsFull() => NumEdges() == 0 && Dimension() == 2 && NumChains() > 0;

    // Returns an arbitrary point P along with a boolean indicating whether P is
    // contained by the shape.  (The boolean value must be false for shapes that
    // do not have an interior.)
    //
    // This ReferencePoint may then be used to compute the containment of other
    // points by counting edge crossings.
    public abstract ReferencePoint GetReferencePoint();

    // Returns the number of contiguous edge chains in the shape.  For example,
    // a shape whose edges are [AB, BC, CD, AE, EF] would consist of two chains
    // (AB,BC,CD and AE,EF).  Every chain is assigned a "chain id" numbered
    // sequentially starting from zero.
    //
    // Note that it is always acceptable to implement this method by returning
    // num_edges() (i.e. every chain consists of a single edge), but this may
    // reduce the efficiency of some algorithms.
    public abstract int NumChains();

    // Returns the range of edge ids corresponding to the given edge chain.  The
    // edge chains must form contiguous, non-overlapping ranges that cover the
    // entire range of edge ids.  This is spelled out more formally below:
    //
    // REQUIRES: 0 <= i < num_chains()
    // REQUIRES: chain(i).length >= 0, for all i
    // REQUIRES: chain(0).start == 0
    // REQUIRES: chain(i).start + chain(i).length == chain(i+1).start,
    //           for i < num_chains() - 1
    // REQUIRES: chain(i).start + chain(i).length == num_edges(),
    //           for i == num_chains() - 1
    public abstract Chain GetChain(int chain_id);

    // Returns the edge at offset "offset" within edge chain "chain_id".
    // Equivalent to "shape.edge(shape.chain(chain_id).start + offset)"
    // but may be more efficient.
    public abstract Edge ChainEdge(int chain_id, int offset);

    // Finds the chain containing the given edge, and returns the position of
    // that edge as a (chain_id, offset) pair.
    //
    // REQUIRES: shape.chain(pos.chain_id).start + pos.offset == edge_id
    // REQUIRES: shape.chain(pos.chain_id + 1).start > edge_id
    //
    // where     pos == shape.chain_position(edge_id).
    public abstract ChainPosition GetChainPosition(int edge_id);

    // Returns an integer that can be used to identify the type of an encoded
    // S2Shape (see TypeTag below).
    public virtual TypeTag GetTypeTag() => TypeTag.kNoTypeTag;

    #endregion

    // Appends an encoded representation of the S2Shape to "encoder".  Note that
    // the encoding does *not* include the type_tag(), so the tag will need to
    // be encoded separately if the shape type will be unknown at decoding time
    // (see s2shapeutil::EncodeTaggedShapes() and related functions).
    //
    // The encoded representation should satisfy the following:
    //
    //  - It should include a version number, so that the encoding may be changed
    //    or improved in the future.
    //
    //  - If "hint" is CodingHint::FAST, the encoding should be optimized for
    //    decoding performance.  If "hint" is CodingHint::COMPACT, the encoding
    //    should be optimized for space.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public virtual void Encode(Encoder encoder, CodingHint hint)
    {
        throw new NotImplementedException("(FATAL) Encoding not implemented for this S2Shape type");
    }

    public enum TypeTag : System.UInt32
    {
        // Indicates that a given S2Shape type cannot be encoded.
        kNoTypeTag = 0,
        S2Polygon = 1,
        S2Polyline = 2,
        S2PointVectorShape = 3,
        S2LaxPolylineShape = 4,
        S2LaxPolygonShape = 5,

        // The following constant should be updated whenever new types are added.
        kNextAvailableTypeTag = 6,

        // The minimum allowable tag for user-defined S2Shape types.
        kMinUserTypeTag = 8192,
    }

    // An edge, consisting of two vertices "v0" and "v1".  Zero-length edges are
    // allowed, and can be used to represent points.
    public readonly record struct Edge(S2Point V0, S2Point V1) : IComparable<Edge>
    {
        #region IComparable

        public int CompareTo(Edge other)
        {
            var c = V0.CompareTo(other.V0);
            if (c != 0)
                return c;

            return V1.CompareTo(other.V1);
        }

        public static bool operator <(Edge x, Edge y) => x.CompareTo(y) < 0;
        public static bool operator >(Edge x, Edge y) => x.CompareTo(y) > 0;
        public static bool operator <=(Edge x, Edge y) => x.CompareTo(y) <= 0;
        public static bool operator >=(Edge x, Edge y) => x.CompareTo(y) >= 0;

        #endregion
    }

    // A range of edge ids corresponding to a chain of zero or more connected
    // edges, specified as a (start, length) pair.  The chain is defined to
    // consist of edge ids {start, start + 1, ..., start + length - 1}.
    public readonly record struct Chain(Int32 Start, Int32 Length);

    // The position of an edge within a given edge chain, specified as a
    // (chain_id, offset) pair.  Chains are numbered sequentially starting from
    // zero, and offsets are measured from the start of each chain.
    public readonly record struct ChainPosition(Int32 ChainId, Int32 Offset);

    // A ReferencePoint consists of a point P and a boolean indicating whether P
    // is contained by a particular shape.
    public record ReferencePoint(S2Point Point, bool Contained)
    {
        // Returns a ReferencePoint with the given "contained" value and a default
        // "point".  It should be used when all points or no points are contained.
        public static ReferencePoint FromContained(bool _contained) => new(S2.Origin, _contained);
    }
}

// Controls whether to optimize for speed or size when encoding shapes.  (Note
// that encoding is always lossless, and that compact encodings are currently
// only possible when points have been snapped to S2CellId centers.)
public enum CodingHint : byte { FAST, COMPACT }
