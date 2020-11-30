using System;
using S2Geometry.S2ShapeUtil;

namespace S2Geometry
{
    // S2LaxLoopShape represents a closed loop of edges surrounding an interior
    // region.  It is similar to S2Loop.Shape except that this class allows
    // duplicate vertices and edges.  Loops may have any number of vertices,
    // including 0, 1, or 2.  (A one-vertex loop defines a degenerate edge
    // consisting of a single point.)
    //
    // Note that S2LaxLoopShape is faster to initialize and more compact than
    // S2Loop.Shape, but does not support the same operations as S2Loop.
    public class S2LaxLoopShape : S2Shape
    {
        #region Fields, Constants
        
        private S2Point[] vertices_;
        public int NumVertices { get; private set; } 

        #endregion

        #region Constructors

        // Constructs an empty loop.
        public S2LaxLoopShape() { NumVertices = 0; }

        // Constructs an S2LaxLoopShape with the given vertices.
        public S2LaxLoopShape(S2Point[] vertices) { Init(vertices); }

        // Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
        public S2LaxLoopShape(S2Loop loop) { Init(loop); }

        #endregion

        #region S2LaxLoopShape
        
        // Initializes an S2LaxLoopShape with the given vertices.
        public void Init(S2Point[] vertices)
        {
            NumVertices = vertices.Length;
            vertices_ = vertices;
        }

        // Initializes an S2LaxLoopShape from the given S2Loop, by copying its data.
        //
        // REQUIRES: !loop.IsFull
        //           [Use S2LaxPolygonShape if you need to represent a full loop.]
        public void Init(S2Loop loop)
        {
            Assert.True(!loop.IsFull); // Full loops not supported; use S2LaxPolygonShape
            if (loop.IsEmpty)
            {
                NumVertices = 0;
                vertices_ = null;
            }
            else
            {
                NumVertices = loop.NumVertices;
                vertices_ = loop.CloneVertices();
            }
        }

        public S2Point Vertex(int i) => vertices_[i]; 

        #endregion

        #region S2Shape

        // S2Shape interface:
        public sealed override int NumEdges => NumVertices;
        public sealed override Edge GetEdge(int e0)
        {
            Assert.True(e0 < NumEdges);
            int e1 = e0 + 1;
            if (e1 == NumVertices) e1 = 0;
            return new Edge(vertices_[e0], vertices_[e1]);
        }
        // Not final; overridden by S2LaxClosedPolylineShape.
        public override int Dimension() { return 2; }
        // Not final; overridden by S2LaxClosedPolylineShape.
        public override ReferencePoint GetReferencePoint()
        {
            return S2ShapeX.GetReferencePoint(this);
        }
        public sealed override int NumChains() => Math.Min(1, NumVertices);
        public sealed override Chain GetChain(int i) => new Chain(0, NumVertices);
        public sealed override Edge ChainEdge(int i, int j)
        {
            Assert.True(i == 0);
            Assert.True(j < NumEdges);
            int k = (j + 1 == NumVertices) ? 0 : j + 1;
            return new Edge(vertices_[j], vertices_[k]);
        }
        public sealed override ChainPosition GetChainPosition(int e) => new ChainPosition(0, e); 

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
        public sealed override ReferencePoint GetReferencePoint()
            => ReferencePoint.FromContained(false);
    }

    // S2VertexIdLaxLoopShape is just like S2LaxLoopShape, except that vertices are
    // specified as indices into a vertex array.  This representation can be more
    // compact when many loops are arranged in a mesh structure.
    public class S2VertexIdLaxLoopShape : S2Shape
    {
        #region Fields, Constants

        private Int32[] vertex_ids_;
        private S2Point[] vertex_array_; 

        #endregion

        #region Constructors

        // Constructs an empty loop.
        public S2VertexIdLaxLoopShape() { NumVertices = 0; }

        // Constructs the shape from the given vertex array and indices.
        // "vertex_ids" is a vector of indices into "vertex_array".
        //
        // ENSURES:  loop.vertex(i) == (*vertex_array)[vertex_ids[i]]
        // REQUIRES: "vertex_array" persists for the lifetime of this object.
        public S2VertexIdLaxLoopShape(Int32[] vertex_ids, S2Point[] vertex_array)
        {
            Init(vertex_ids, vertex_array);
        }

        #endregion

        #region S2VertexIdLaxLoopShape

        // Initializes the shape from the given vertex array and indices.
        // "vertex_ids" is a vector of indices into "vertex_array".
        public void Init(Int32[] vertex_ids, S2Point[] vertex_array)
        {
            NumVertices = vertex_ids.Length;
            vertex_ids_ = (Int32[])vertex_ids.Clone();
            vertex_array_ = vertex_array;
        }

        // Returns the number of vertices in the loop.
        public int NumVertices { get; private set; }
        public Int32 VertexId(int i) => vertex_ids_[i];
        public S2Point Vertex(int i) => vertex_array_[VertexId(i)]; 

        #endregion

        #region S2Shape

        // S2Shape interface:
        public sealed override int NumEdges => NumVertices;
        public sealed override Edge GetEdge(int e0)
        {
            Assert.True(e0 < NumEdges);
            int e1 = e0 + 1;
            if (e1 == NumVertices) e1 = 0;
            return new Edge(Vertex(e0), Vertex(e1));
        }
        public sealed override int Dimension() => 2;
        public sealed override ReferencePoint GetReferencePoint()
        {
            // GetReferencePoint interprets a loop with no vertices as "full".
            if (NumVertices == 0) return ReferencePoint.FromContained(false);
            return S2ShapeX.GetReferencePoint(this);
        }
        public sealed override int NumChains() => Math.Min(1, NumVertices);
        public sealed override Chain GetChain(int i) => new Chain(0, NumVertices);
        public sealed override Edge ChainEdge(int i, int j)
        {
            Assert.True(i == 0);
            Assert.True(j < NumEdges);
            int k = (j + 1 == NumVertices) ? 0 : j + 1;
            return new Edge(Vertex(j), Vertex(k));
        }
        public sealed override ChainPosition GetChainPosition(int e) => new ChainPosition(0, e);

        #endregion
    }
}

