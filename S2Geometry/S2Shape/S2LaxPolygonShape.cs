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
//    at most one of them must belong to a sibling pair (i.e., the number of
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
// Note that S2LaxPolygonShape is much faster to initialize and is more
// compact than S2Polygon, but unlike S2Polygon it does not have any built-in
// operations.  Instead you should use S2ShapeIndex operations
// (S2BooleanOperation, S2ClosestEdgeQuery, etc).

namespace S2Geometry;

using Loop = List<S2Point>;

public class S2LaxPolygonShape : S2Shape, IInitEncoder<S2LaxPolygonShape>
{
    #region Fields, Constants

    public const TypeTag kTypeTag = TypeTag.S2LaxPolygonShape;

    // Returns the number of loops.
    //
    // Note that the parent class has a 4-byte S2Shape::id_ field so there is no
    // wasted space in the following layout.
    public int NumLoops { get; private set; }

    // The loop that contained the edge returned by the previous call to the
    // edge() method.  This is used as a hint to speed up edge location when
    // there are many loops.  Benchmarks indicate that the improved locality
    // this provides can speed up chain position lookup by 1.7-4.7x.
    private int prev_loop_ = 0;

    public int NumVertices { get; private set; }
    private S2Point[]? vertices_;

    // When num_loops_ > 1, stores an array of size (num_loops_ + 1) where
    // element "i" represents the total number of vertices in loops 0..i-1.
    private UInt32[]? loop_starts_;

    public const byte kCurrentEncodingVersionNumber = 1;

    #endregion

    #region Constructors

    // Constructs an empty polygon.
    public S2LaxPolygonShape()
    {
        NumLoops = 0;
        NumVertices = 0;
    }

    // Constructs an S2LaxPolygonShape from the given vertex loops.
    public S2LaxPolygonShape(List<Loop> loops)
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

    // Initializes an S2LaxPolygonShape from an S2Polygon by copying its data.
    // Full and empty S2Polygons are supported.
    public void Init(S2Polygon polygon)
    {
        var spans = new List<List<S2Point>>();
        for (int i = 0; i < polygon.NumLoops(); ++i)
        {
            var loop = polygon.Loop(i);
            if (loop.IsFull())
            {
                spans.Add([]);
            }
            else
            {
                spans.Add([.. loop.CloneVertices()]);
            }
        }
        Init(spans);

        // S2Polygon and S2LaxPolygonShape holes are oriented oppositely, so we need
        // to reverse the orientation of any loops representing holes.
        for (int i = 0; i < polygon.NumLoops(); ++i)
        {
            if (polygon.Loop(i).IsHole())
            {
                Array.Reverse(vertices_, (int)loop_starts_[i], NumLoopVertices(i));
            }
        }
    }

    public void Init(List<Loop> loops)
    {
        NumLoops = loops.Count;
        if (NumLoops == 0)
        {
            NumVertices = 0;
        }
        else if (NumLoops == 1)
        {
            NumVertices = loops[0].Count;
            // TODO(ericv): Use std::allocator to obtain uninitialized memory instead.
            // This would avoid default-constructing all the elements before we
            // overwrite them, and it would also save 8 bytes of memory allocation
            // since "new T[]" stores its own copy of the array size.
            //
            // Note that even absl::make_unique_for_overwrite<> and c++20's
            // std::make_unique_for_overwrite<T[]> default-construct all elements when
            // T is a class type.
            vertices_ = [.. loops[0]];
        }
        else
        {
            loop_starts_ = new UInt32[NumLoops + 1];
            NumVertices = 0;
            for (int i = 0; i < NumLoops; ++i)
            {
                loop_starts_[i] = (uint)NumVertices;
                NumVertices += loops[i].Count;
            }
            loop_starts_[NumLoops] = (uint)NumVertices;
            vertices_ = loops.SelectMany(t => t).ToArray();  // TODO(see above)
        }
    }

    // Returns the number of vertices in the given loop.
    public int NumLoopVertices(int i)
    {
        MyDebug.Assert(i < NumLoops);
        if (NumLoops == 1)
        {
            return NumVertices;
        }
        else
        {
            return (int)(loop_starts_![i + 1] - loop_starts_[i]);
        }
    }

#if s2debug

    // Returns the vertex from loop "i" at index "j".
    // REQUIRES: 0 <= i < num_loops()
    // REQUIRES: 0 <= j < num_loop_vertices(i)
    public S2Point LoopVertex(int i, int j)
    {
        MyDebug.Assert(i < NumLoops);
        MyDebug.Assert(j < NumLoopVertices(i));
        if (i == 0)
        {
            return vertices_.Skip(j).FirstOrDefault();
        }
        else
        {
            return vertices_.Skip((int)loop_starts_![i] + j).FirstOrDefault();
        }
    }

    // Returns the vertex from loop "i" at index "j".
    // REQUIRES: 0 <= i < num_loops()
    // REQUIRES: 0 <= j < num_loop_vertices(i)
    public IEnumerable<S2Point> LoopVertices(int i, int j)
    {
        MyDebug.Assert(i < NumLoops);
        MyDebug.Assert(j < NumLoopVertices(i));
        if (NumLoops == 1)
        {
            return vertices_.Skip(j);
        }
        else
        {
            return vertices_.Skip((int)loop_starts_[i] + j);
        }
    }

#endif

    #endregion

    #region IInitEncoder

    // Decodes an S2LaxPolygonShape, returning true on success.  (The method
    // name is chosen for compatibility with EncodedS2LaxPolygonShape below.)
    public static (bool, S2LaxPolygonShape?) Init(Decoder decoder)
    {
        if (decoder.Avail() < 1) return (false, null);
        var version = decoder.Get8();
        if (version != kCurrentEncodingVersionNumber) return (false, null);

        if (!decoder.TryGetVarUInt32(out var numLoops)) return (false, null);
        var (success, vertices) = EncodedS2PointVector.Init(decoder);
        if (!success) return (false, null);

        var numVertices = 0;
        S2Point[]? vertices_ = null;
        UInt32[]? loop_starts_ = null;
        if (numLoops != 0)
        {
            numVertices = vertices!.Count();
            vertices_ = new S2Point[numVertices];  // TODO(see above)
            for (int i = 0; i < numVertices; ++i)
            {
                vertices_[i] = vertices[i];
            }

            if (numLoops > 1)
            {
                var (success2, loop_starts) = EncodedUIntVector<UInt32>.Init(decoder);
                if (!success2) return (false, null);
                loop_starts_ = new UInt32[loop_starts!.Count];
                for (int i = 0; i < loop_starts.Count; ++i)
                {
                    loop_starts_[i] = loop_starts[i];
                }
            }
        }

        S2LaxPolygonShape shape = new()
        {
            NumLoops = (int)numLoops,
            NumVertices = numVertices,
            vertices_ = vertices_,
            loop_starts_ = loop_starts_,
        };
        return (true, shape);
    }

    // Appends an encoded representation of the S2LaxPolygonShape to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(1 + Encoder.kVarintMax32);
        encoder.Put8(kCurrentEncodingVersionNumber);
        encoder.PutVarInt32(NumLoops);
        EncodedS2PointVector.EncodeS2PointVector(vertices_, hint, encoder);
        if (NumLoops > 1)
        {
            EncodedUIntVector<uint>.EncodeUIntVector(loop_starts_!, encoder);
        }
    } 

    #endregion

    #region S2Shape

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return NumVertices;
    }

    public sealed override Edge GetEdge(int e)
    {
        // Method names are fully specified to enable inlining.
        ChainPosition pos = GetChainPositionInternal(e);
        return ChainEdgeInternal(pos.ChainId, pos.Offset);
    }
    public sealed override int Dimension() => 2;
    public sealed override ReferencePoint GetReferencePoint() =>
        S2ShapeUtil.GetReferencePoint(this);
    public sealed override int NumChains() => NumLoops;
    public sealed override Chain GetChain(int i)
    {
        MyDebug.Assert(i < NumLoops);
        if (NumLoops == 1)
        {
            return new Chain(0, NumVertices);
        }
        else
        {
            int start = (int)loop_starts_[i];
            return new Chain(start, (int)(loop_starts_[i + 1] - start));
        }
    }
    public sealed override Edge ChainEdge(int i, int j) => ChainEdgeInternal(i, j);
    private Edge ChainEdgeInternal(int i, int j)
    {
        MyDebug.Assert(i < NumLoops);
        MyDebug.Assert(j < NumLoopVertices(i));
        int n = NumLoopVertices(i);
        int k = (j + 1 == n) ? 0 : j + 1;
        if (NumLoops == 1)
        {
            return new Edge(vertices_[j], vertices_[k]);
        }
        else
        {
            int start = (int)loop_starts_[i];
            return new Edge(vertices_[start + j], vertices_[start + k]);
        }
    }
    public sealed override ChainPosition GetChainPosition(int e) => GetChainPositionInternal(e);
    private ChainPosition GetChainPositionInternal(int e)
    {
        MyDebug.Assert(e < NumEdges());
        if (NumLoops == 1)
        {
            return new ChainPosition(0, e);
        }
        // Test if this edge belongs to the loop returned by the previous call.
        var start = 0 + prev_loop_; //.load(std::memory_order_relaxed);
        if (e >= loop_starts_[start] && e < loop_starts_[start + 1])
        {
            // This edge belongs to the same loop as the previous call.
        }
        else
        {
            if (e == loop_starts_[start + 1])
            {
                // This is the edge immediately following the previous loop.
                do { ++start; } while (e == loop_starts_[start + 1]);
            }
            else
            {
                start = 0;
                const int kMaxLinearSearchLoops = 12;  // From benchmarks.
                if (NumLoops <= kMaxLinearSearchLoops)
                {
                    while (loop_starts_[start + 1] <= e) ++start;
                }
                else
                {
                    start = loop_starts_.GetUpperBound((uint)e, start + 1, NumLoops) - 1;
                }
            }
            prev_loop_ = start;
        }
        return new ChainPosition(start, (int)((uint)e - loop_starts_[start]));
    }

    // Define as enum so we don't have to declare storage.
    // TODO(user, b/210097200): Use static constexpr when C++17 is allowed
    // in opensource.
    public override TypeTag GetTypeTag() => kTypeTag;

    #endregion
}

// Exactly like S2LaxPolygonShape, except that the vertices are kept in an
// encoded form and are decoded only as they are accessed.  This allows for
// very fast initialization and no additional memory use beyond the encoded
// data.  The encoded data is not owned by this class; typically it points
// into a large contiguous buffer that contains other encoded data as well.
public class EncodedS2LaxPolygonShape : S2Shape, IInitEncoder<EncodedS2LaxPolygonShape>
{
    // The loop that contained the edge returned by the previous call to the
    // edge() method.  This is used as a hint to speed up edge location when
    // there are many loops.
    private int prev_loop_ = 0;

    public required EncodedS2PointVector Vertices { private get; init; }
    public EncodedUIntVector<UInt32>? LoopStarts { private get; init; }

    // Constructs an uninitialized object; requires Init() to be called.
    public EncodedS2LaxPolygonShape() { }

    // Initializes an EncodedS2LaxPolygonShape.
    //
    // REQUIRES: The Decoder data buffer must outlive this object.
    public static (bool, EncodedS2LaxPolygonShape?) Init(Decoder decoder)
    {
        if (decoder.Avail() < 1) return (false, null);
        var version = decoder.Get8();
        if (version != S2LaxPolygonShape.kCurrentEncodingVersionNumber) return (false, null);

        if (!decoder.TryGetVarUInt32(out var num_loops)) return (false, null);

        var (success_e, v) = EncodedS2PointVector.Init(decoder);
        if (!success_e) return (false, null);
        EncodedUIntVector<UInt32>? loop_starts_ = null;
        if (num_loops > 1)
        {
            var (success, ls) = EncodedUIntVector<UInt32>.Init(decoder);
            if (!success) return (false, null);
            loop_starts_ = ls;
        }
        EncodedS2LaxPolygonShape shape = new()
        {
            NumLoops = (int)num_loops,
            LoopStarts = loop_starts_,
            Vertices = v!,
        };
        return (true, shape);
    }

    // Appends an encoded representation of the S2LaxPolygonShape to "encoder".
    // The coding hint is ignored, and whatever method was originally used to
    // encode the shape is preserved.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    //
    // The encoding must be identical to S2LaxPolygonShape::Encode().
    public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(1 + Encoder.kVarintMax32);
        encoder.Put8(S2LaxPolygonShape.kCurrentEncodingVersionNumber);
        encoder.PutVarInt32(NumLoops);
        Vertices.Encode(encoder);
        if (NumLoops > 1)
        {
            LoopStarts.Encode(encoder);
        }
    }

    public int NumLoops { get; private set; }
    public int NumVertices
    {
        get
        {
            if (NumLoops <= 1)
            {
                return Vertices.Count();
            }
            else
            {
                return (int)LoopStarts[NumLoops];
            }
        }
    }
    public int NumLoopVertices(int i)
    {
        MyDebug.Assert(i < NumLoops);
        if (NumLoops == 1)
        {
            return Vertices.Count();
        }
        else
        {
            return (int)(LoopStarts[i + 1] - LoopStarts[i]);
        }
    }
    public S2Point LoopVertex(int i, int j)
    {
        MyDebug.Assert(i < NumLoops);
        MyDebug.Assert(j < NumLoopVertices(i));
        if (NumLoops == 1)
        {
            return Vertices[j];
        }
        else
        {
            return Vertices[(int)LoopStarts[i] + j];
        }
    }

    // S2Shape interface:
    public sealed override int NumEdges()
    {
        return NumVertices;
    }

    public sealed override Edge GetEdge(int e)
    {
        MyDebug.Assert(e < NumEdges());
        int e1 = e + 1;
        if (NumLoops == 1)
        {
            if (e1 == Vertices.Count()) { e1 = 0; }

            return new Edge(Vertices[e], Vertices[e1]);
        }
        else
        {
            // Method names are fully specified to enable inlining.
            ChainPosition pos = GetChainPositionInternal(e);
            return ChainEdgeInternal(pos.ChainId, pos.Offset);
        }
    }
    public sealed override int Dimension() { return 2; }
    public sealed override ReferencePoint GetReferencePoint()
    {
        return S2ShapeUtil.GetReferencePoint(this);
    }
    public sealed override int NumChains() { return NumLoops; }
    public sealed override Chain GetChain(int i)
    {
        MyDebug.Assert(i < NumLoops);
        if (NumLoops == 1)
        {
            return new Chain(0, Vertices.Count());
        }
        else
        {
            int start = (int)LoopStarts[i];
            return new Chain(start, (int)(LoopStarts[i + 1] - start));
        }
    }
    public sealed override Edge ChainEdge(int i, int j) => ChainEdgeInternal(i, j);
    private Edge ChainEdgeInternal(int i, int j)
    {
        MyDebug.Assert(i < NumLoops);
        MyDebug.Assert(j < NumLoopVertices(i));
        int n = NumLoopVertices(i);
        int k = (j + 1 == n) ? 0 : j + 1;
        if (NumLoops == 1)
        {
            return new Edge(Vertices[j], Vertices[k]);
        }
        else
        {
            int start = (int)LoopStarts[i];
            return new Edge(Vertices[start + j], Vertices[start + k]);
        }
    }
    public sealed override ChainPosition GetChainPosition(int e) => GetChainPositionInternal(e);
    private ChainPosition GetChainPositionInternal(int e)
    {
        MyDebug.Assert(e < NumEdges());
        if (NumLoops == 1)
        {
            return new ChainPosition(0, e);
        }
        const int kMaxLinearSearchLoops = 12;  // From benchmarks.
        int i = prev_loop_;//.load(std::memory_order_relaxed);
        if (i == 0 && e < LoopStarts[1])
        {
            return new ChainPosition(0, e);  // Optimization for first loop.
        }

        if (e >= LoopStarts[i] && e < LoopStarts[i + 1])
        {
            // This edge belongs to the same loop as the previous call.
        }
        else
        {
            if (e == LoopStarts[i + 1])
            {
                // This is the edge immediately following the previous loop.
                do { ++i; } while (e == LoopStarts[i + 1]);
            }
            else if (NumLoops <= kMaxLinearSearchLoops)
            {
                for (i = 0; LoopStarts[i + 1] <= e; ++i) { }
            }
            else
            {
                i = LoopStarts.LowerBound((uint)(e + 1)) - 1;
            }
            prev_loop_ = i; //.store(i, std::memory_order_relaxed);
        }
        return new ChainPosition(i, e - (int)LoopStarts[i]);
    }

    public override TypeTag GetTypeTag() => S2LaxPolygonShape.kTypeTag;
}
