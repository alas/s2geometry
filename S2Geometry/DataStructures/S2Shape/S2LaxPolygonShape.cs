using System;
using System.Collections.Generic;
using System.Linq;
using S2Geometry.S2ShapeUtil;

namespace S2Geometry
{
    // S2LaxPolygonShape represents a region defined by a collection of zero or
    // more closed loops.  The interior is the region to the left of all loops.
    // This is similar to S2Polygon.Shape except that this class supports
    // polygons with degeneracies.  Degeneracies are of two types: degenerate
    // edges (from a vertex to itself) and sibling edge pairs (consisting of two
    // oppositely oriented edges).  Degeneracies can represent either "shells" or
    // "holes" depending on the loop they are contained by.  For example, a
    // degenerate edge or sibling pair contained by a "shell" would be interpreted
    // as a degenerate hole.  Such edges form part of the boundary of the polygon.
    //
    // Loops with fewer than three vertices are interpreted as follows:
    //  - A loop with two vertices defines two edges (in opposite directions).
    //  - A loop with one vertex defines a single degenerate edge.
    //  - A loop with no vertices is interpreted as the "full loop" containing
    //    all points on the sphere.  If this loop is present, then all other loops
    //    must form degeneracies (i.e., degenerate edges or sibling pairs).  For
    //    example, two loops {} and {X} would be interpreted as the full polygon
    //    with a degenerate single-point hole at X.
    //
    // S2LaxPolygonShape does not have any error checking, and it is perfectly
    // fine to create S2LaxPolygonShape objects that do not meet the requirements
    // below (e.g., in order to analyze or fix those problems).  However,
    // S2LaxPolygonShapes must satisfy some additional conditions in order to
    // perform certain operations:
    //
    //  - In order to be valid for point containment tests, the polygon must
    //    satisfy the "interior is on the left" rule.  This means that there must
    //    not be any crossing edges, and if there are duplicate edges then all but
    //    at most one of thm must belong to a sibling pair (i.e., the number of
    //    edges in opposite directions must differ by at most one).
    //
    //  - To be valid for boolean operations (S2BooleanOperation), degenerate
    //    edges and sibling pairs cannot coincide with any other edges.  For
    //    example, the following situations are not allowed:
    //
    //      {AA, AA}      // degenerate edge coincides with another edge
    //      {AA, AB}      // degenerate edge coincides with another edge
    //      {AB, BA, AB}  // sibling pair coincides with another edge
    //
    // Note that S2LaxPolygonShape is must faster to initialize and is more
    // compact than S2Polygon, but unlike S2Polygon it does not have any built-in
    // operations.  Instead you should use S2ShapeIndex operations
    // (S2BooleanOperation, S2ClosestEdgeQuery, etc).
    public class S2LaxPolygonShape : S2Shape, IEncodeInit
    {
        #region Fields, Constants

        // Returns the number of loops.
        public int NumLoops { get; private set; }

        private S2Point[] vertices_;

        // If num_loops_ <= 1, this union stores the number of vertices.
        // Otherwise it points to an array of size (num_loops + 1) where element "i"
        // is the total number of vertices in loops 0..i-1.
        private Int32 num_vertices_;

        private UInt32[] cumulative_vertices_;

#pragma warning disable IDE1006 // Estilos de nombres
        public const byte kCurrentEncodingVersionNumber = 1;
#pragma warning restore IDE1006 // Estilos de nombres

        #endregion

        #region Constructors

        // Constructs an empty polygon.
        public S2LaxPolygonShape()
        {
            NumLoops = 0;
            num_vertices_ = 0;
        }

        // Constructs an S2LaxPolygonShape from the given vertex loops.
        public S2LaxPolygonShape(List<List<S2Point>> loops)
        {
            Init(loops);
        }

        // Constructs an S2LaxPolygonShape from an S2Polygon, by copying its data.
        // Full and empty S2Polygons are supported.
        public S2LaxPolygonShape(S2Polygon polygon)
        {
            Init(polygon);
        }

        #endregion

        #region S2LaxPolygonShape

        // Initializes an S2LaxPolygonShape from an S2Polygon, by copying its data.
        // Full and empty S2Polygons are supported.
        public void Init(S2Polygon polygon)
        {
            var spans = new List<List<S2Point>>();
            for (int i = 0; i < polygon.NumLoops(); ++i)
            {
                var loop = polygon.Loop(i);
                if (loop.IsFull)
                {
                    spans.Add(new List<S2Point>());
                }
                else
                {
                    spans.Add(loop.CloneVertices().ToList());
                }
            }
            Init(spans);

            // S2Polygon and S2LaxPolygonShape holes are oriented oppositely, so we need
            // to reverse the orientation of any loops representing holes.
            for (int i = 0; i < polygon.NumLoops(); ++i)
            {
                if (polygon.Loop(i).IsHole)
                {
                    Array.Reverse(vertices_, (int)cumulative_vertices_[i], NumLoopVertices(i));
                }
            }
        }

        // Decodes an S2LaxPolygonShape, returning true on success.  (The method
        // name is chosen for compatibility with EncodedS2LaxPolygonShape below.)
        public bool Init(Decoder decoder)
        {
            if (decoder.Avail() < 1) return false;
            var version = decoder.Get8();
            if (version != kCurrentEncodingVersionNumber) return false;

            if (!decoder.TryGetVarUInt32(out var num_loops)) return false;
            NumLoops = (int)num_loops;
            var vertices = new EncodedS2PointVector();
            if (!vertices.Init(decoder)) return false;

            if (NumLoops == 0)
            {
                num_vertices_ = 0;
                vertices_ = null;
            }
            else
            {
                vertices_ = new S2Point[vertices.Count()];
                for (int i = 0; i < vertices.Count(); ++i)
                {
                    vertices_[i] = vertices[i];
                }
                if (NumLoops == 1)
                {
                    num_vertices_ = vertices.Count();
                }
                else
                {
                    var cumulative_vertices = new EncodedUintVector_UInt32(); ;
                    if (!cumulative_vertices.Init(decoder)) return false;
                    cumulative_vertices_ = new UInt32[cumulative_vertices.Count];
                    for (int i = 0; i < cumulative_vertices.Count; ++i)
                    {
                        cumulative_vertices_[i] = cumulative_vertices[i];
                    }
                }
            }
            return true;
        }

        public void Init(List<List<S2Point>> loops)
        {
            NumLoops = loops.Count;
            if (NumLoops == 0)
            {
                num_vertices_ = 0;
                vertices_ = null;
            }
            else if (NumLoops == 1)
            {
                num_vertices_ = loops[0].Count;
                vertices_ = loops[0].ToArray();
            }
            else
            {
                cumulative_vertices_ = new UInt32[NumLoops + 1];
                Int32 num_vertices = 0;
                for (int i = 0; i < NumLoops; ++i)
                {
                    cumulative_vertices_[i] = (uint)num_vertices;
                    num_vertices += loops[i].Count;
                }
                cumulative_vertices_[NumLoops] = (uint)num_vertices;
                vertices_ = loops.SelectMany(t => t).ToArray();
            }
        }

        // Returns the total number of vertices in all loops.
        public int NumVertices
        {
            get
            {
                if (NumLoops <= 1)
                {
                    return num_vertices_;
                }
                else
                {
                    return (int)cumulative_vertices_[NumLoops];
                }
            }
        }

        // Returns the number of vertices in the given loop.
        public int NumLoopVertices(int i)
        {
            Assert.True(i < NumLoops);
            if (NumLoops == 1)
            {
                return num_vertices_;
            }
            else
            {
                return (int)(cumulative_vertices_[i + 1] - cumulative_vertices_[i]);
            }
        }

        // Appends an encoded representation of the S2LaxPolygonShape to "encoder".
        //
        // REQUIRES: "encoder" uses the defaultructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public void Encode(Encoder encoder)
        {
            Encode(encoder, CodingHint.COMPACT);
        }
        public void Encode(Encoder encoder, CodingHint hint)
        {
            encoder.Ensure(1 + Encoder.kVarintMax32);
            encoder.Put8(kCurrentEncodingVersionNumber);
            encoder.PutVarInt32(NumLoops);
            EncodedS2PointVector.EncodeS2PointVector(vertices_, hint, encoder);
            if (NumLoops > 1)
            {
                EncodedUintVector.EncodeUintVector(cumulative_vertices_, encoder);
            }
        }

#if s2debug
        // Returns the vertex from loop "i" at index "j".
        // REQUIRES: 0 <= i < num_loops()
        // REQUIRES: 0 <= j < num_loop_vertices(i)
        public IEnumerable<S2Point> LoopVertex(int i, int j)
        {
            Assert.True(i < NumLoops);
            Assert.True(j < NumLoopVertices(i));
            if (NumLoops == 1)
            {
                return vertices_.Skip(j);
            }
            else
            {
                return vertices_.Skip((int)cumulative_vertices_[i] + j);
            }
        }
#endif

        #endregion

        #region S2Shape

        // S2Shape interface:
        public sealed override int NumEdges => NumVertices;
        public sealed override Edge GetEdge(int e0)
        {
            Assert.True(e0 < NumEdges);
            int e1 = e0 + 1;
            if (NumLoops== 1)
            {
                if (e1 == num_vertices_) { e1 = 0; }
            }
            else
            {
                // Find the index of the first vertex of the loop following this one.
                int kMaxLinearSearchLoops = 12;  // From benchmarks.
                var next = 1;
                if (NumLoops<= kMaxLinearSearchLoops)
                {
                    while (cumulative_vertices_[next] <= e0) ++next;
                }
                else
                {
                    next = cumulative_vertices_.GetLowerBound((uint)e1, next, next + NumLoops);
                }
                // Wrap around to the first vertex of the loop if necessary.
                if (e1 == cumulative_vertices_[next]) { e1 = (int)cumulative_vertices_[next - 1]; }
            }
            return new Edge(vertices_[e0], vertices_[e1]);
        }
        public sealed override int Dimension() => 2;
        public sealed override ReferencePoint GetReferencePoint()
            => S2ShapeX.GetReferencePoint(this);
        public sealed override int NumChains() => NumLoops;
        public sealed override Chain GetChain(int i)
        {
            Assert.True(i < NumLoops);
            if (NumLoops== 1)
            {
                return new Chain(0, num_vertices_);
            }
            else
            {
                int start = (int)cumulative_vertices_[i];
                return new Chain(start, (int)(cumulative_vertices_[i + 1] - start));
            }
        }
        public sealed override Edge ChainEdge(int i, int j)
        {
            Assert.True(i < NumLoops);
            Assert.True(j < NumLoopVertices(i));
            int n = NumLoopVertices(i);
            int k = (j + 1 == n) ? 0 : j + 1;
            if (NumLoops== 1)
            {
                return new Edge(vertices_[j], vertices_[k]);
            }
            else
            {
                int base_ = (int)cumulative_vertices_[i];
                return new Edge(vertices_[base_ + j], vertices_[base_ + k]);
            }
        }
        public sealed override ChainPosition GetChainPosition(int e)
        {
            Assert.True(e < NumEdges);
            int kMaxLinearSearchLoops = 12;  // From benchmarks.
            if (NumLoops== 1)
            {
                return new S2Shape.ChainPosition(0, e);
            }
            else
            {
                // Find the index of the first vertex of the loop following this one.
                int next = 1;
                if (NumLoops<= kMaxLinearSearchLoops)
                {
                    while (cumulative_vertices_[next] <= e) ++next;
                }
                else
                {
                    next = cumulative_vertices_.GetLowerBound((uint)(e + 1), next, next + NumLoops);
                }
                return new S2Shape.ChainPosition(next - 1, e - (int)cumulative_vertices_[next - 1]);
            }
        }
        public override TypeTag GetTypeTag() => TypeTag.S2LaxPolygonShape;

        #endregion
    }

    // Exactly like S2LaxPolygonShape, except that the vertices are kept in an
    // encoded form and are decoded only as they are accessed.  This allows for
    // very fast initialization and no additional memory use beyond the encoded
    // data.  The encoded data is not owned by this class; typically it points
    // into a large contiguous buffer that contains other encoded data as well.
    public class EncodedS2LaxPolygonShape : S2Shape
    {
        // Constructs an uninitialized object; requires Init() to be called.
        public EncodedS2LaxPolygonShape() { }

        // Initializes an EncodedS2LaxPolygonShape.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder)
        {
            if (decoder.Avail() < 1) return false;
            var version = decoder.Get8();
            if (version != S2LaxPolygonShape.kCurrentEncodingVersionNumber) return false;

            if (!decoder.TryGetVarUInt32(out var num_loops)) return false;
            NumLoops = (int)num_loops;

            if (!vertices_.Init(decoder)) return false;

            if (NumLoops > 1)
            {
                if (!cumulative_vertices_.Init(decoder)) return false;
            }
            return true;
        }

        public int NumLoops { get; private set; }
        public int NumVertices {
            get
            {
                if (NumLoops<= 1)
                {
                    return vertices_.Count();
                }
                else
                {
                    return (int)cumulative_vertices_[NumLoops];
                }
            }
        }
        public int NumLoopVertices(int i)
        {
            Assert.True(i < NumLoops);
            if (NumLoops== 1)
            {
                return vertices_.Count();
            }
            else
            {
                return (int)(cumulative_vertices_[i + 1] - cumulative_vertices_[i]);
            }
        }
        public S2Point LoopVertex(int i, int j)
        {
            Assert.True(i < NumLoops);
            Assert.True(j < NumLoopVertices(i));
            if (NumLoops== 1)
            {
                return vertices_[j];
            }
            else
            {
                return vertices_[(int)cumulative_vertices_[i] + j];
            }
        }

        // S2Shape interface:
        public sealed override int NumEdges => NumVertices;
        public sealed override Edge GetEdge(int e)
        {
            Assert.True(e < NumEdges);
            int e1 = e + 1;
            if (NumLoops== 1)
            {
                if (e1 == vertices_.Count()) { e1 = 0; }
            }
            else
            {
                // Find the index of the first vertex of the loop following this one.
                int kMaxLinearSearchLoops = 12;  // From benchmarks.
                int next = 1;
                if (NumLoops<= kMaxLinearSearchLoops)
                {
                    while (cumulative_vertices_[next] <= e) ++next;
                }
                else
                {
                    next = cumulative_vertices_.LowerBound((uint)e1);
                }
                // Wrap around to the first vertex of the loop if necessary.
                if (e1 == (int)cumulative_vertices_[next])
                {
                    e1 = (int)cumulative_vertices_[next - 1];
                }
            }
            return new Edge(vertices_[e], vertices_[e1]);
        }
        public sealed override int Dimension() { return 2; }
        public sealed override ReferencePoint GetReferencePoint()
        {
            return S2ShapeX.GetReferencePoint(this);
        }
        public sealed override int NumChains() { return NumLoops; }
        public sealed override Chain GetChain(int i)
        {
            Assert.True(i < NumLoops);
            if (NumLoops== 1)
            {
                return new Chain(0, vertices_.Count());
            }
            else
            {
                int start = (int)cumulative_vertices_[i];
                return new Chain(start, (int)(cumulative_vertices_[i + 1] - start));
            }
        }
        public sealed override Edge ChainEdge(int i, int j)
        {
            Assert.True(i < NumLoops);
            Assert.True(j < NumLoopVertices(i));
            int n = NumLoopVertices(i);
            int k = (j + 1 == n) ? 0 : j + 1;
            if (NumLoops== 1)
            {
                return new Edge(vertices_[j], vertices_[k]);
            }
            else
            {
                int base_ = (int)cumulative_vertices_[i];
                return new Edge(vertices_[base_ + j], vertices_[base_ + k]);
            }
        }
        public sealed override ChainPosition GetChainPosition(int e)
        {
            Assert.True(e < NumEdges);
            int kMaxLinearSearchLoops = 12;  // From benchmarks.
            if (NumLoops== 1)
            {
                return new S2Shape.ChainPosition(0, e);
            }
            else
            {
                // Find the index of the first vertex of the loop following this one.
                int next = 1;
                if (NumLoops<= kMaxLinearSearchLoops)
                {
                    while (cumulative_vertices_[next] <= e) ++next;
                }
                else
                {
                    next = cumulative_vertices_.LowerBound((uint)(e + 1));
                }
                return new S2Shape.ChainPosition(next - 1, e - (int)cumulative_vertices_[next - 1]);
            }
        }

        private readonly EncodedS2PointVector vertices_ = new EncodedS2PointVector();
        private readonly EncodedUintVector_UInt32 cumulative_vertices_ = new EncodedUintVector_UInt32();
    }
}
