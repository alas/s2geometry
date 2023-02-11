// An S2Loop represents a simple spherical polygon.  It consists of a single
// chain of vertices where the first vertex is implicitly connected to the
// last. All loops are defined to have a CCW orientation, i.e. the interior of
// the loop is on the left side of the edges.  This implies that a clockwise
// loop enclosing a small area is interpreted to be a CCW loop enclosing a
// very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not).  Non-adjacent edges are not allowed to intersect, and furthermore edges
// of length 180 degrees are not allowed (i.e., adjacent vertices cannot be
// antipodal).  Loops must have at least 3 vertices (except for the empty and
// full loops discussed below).  Although these restrictions are not enforced
// in optimized code, you may get unexpected results if they are violated.
//
// There are two special loops: the "empty loop" contains no points, while the
// "full loop" contains all points.  These loops do not have any edges, but to
// preserve the invariant that every loop can be represented as a vertex
// chain, they are defined as having exactly one vertex each (see kEmpty and
// kFull).
//
// Point containment of loops is defined such that if the sphere is subdivided
// into faces (loops), every point is contained by exactly one face.  This
// implies that loops do not necessarily contain their vertices.
//
// Note: The reason that duplicate vertices and intersecting edges are not
// allowed is that they make it harder to define and implement loop
// relationships, e.g. whether one loop contains another.  If your data does
// not satisfy these restrictions, you can use S2Builder to normalize it.

namespace S2Geometry;

using System;
using System.Runtime.InteropServices;

public sealed record class S2Loop : IS2Region<S2Loop>, IComparable<S2Loop>, IDecoder<S2Loop>
{
    #region Fields, Constants

    // The depth of a loop is defined as its nesting level within its containing
    // polygon.  "Outer shell" loops have depth 0, holes within those loops have
    // depth 1, shells within those holes have depth 2, etc.  This field is only
    // used by the S2Polygon implementation.
    //
    // The nesting depth, if this field belongs to an S2Polygon.  We define it
    // here to optimize field packing.
    public int Depth { get; set; }

    // Returns true if this loop contains S2.Origin.
    public bool ContainsOrigin { get; private set; } = false;

    public int NumVertices { get; private set; } = 0;

    private int _firstLogicalVertex;

    // We store the vertices in an array rather than a vector because we don't
    // need any STL methods, and computing the number of vertices using size()
    // would be relatively expensive (due to division by SizeHelper.SizeOf(typeof(S2Point)) == 24).
    public S2Point[] Vertices { get; init; }

    // In general we build the index the first time it is needed, but we make an
    // exception for Contains(S2Point) because this method has a simple brute
    // force implementation that is also relatively cheap.  For this one method
    // we keep track of the number of calls made and only build the index once
    // enough calls have been made that we think an index would be worthwhile.
    private Int32 _unindexedContainsCalls_;

    // "bound_" is a conservative bound on all points contained by this loop:
    // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
    private S2LatLngRect _bound;

    // Since "bound_" is not exact, it is possible that a loop A contains
    // another loop B whose bounds are slightly larger.  "subregion_bound_"
    // has been expanded sufficiently to account for this error, i.e.
    // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
    private S2LatLngRect _subregionBound;

    // Spatial index for this loop.
    private MutableS2ShapeIndex _index = new();

    //    "Build the S2ShapeIndex only when it is first needed.  This can save "
    //    "significant amounts of memory and time when geometry is constructed but "
    //    "never queried, for example when loops are passed directly to S2Polygon, "
    //    "or when geometry is being converted from one format to another.");

    // The maximum number of vertices we'll allow when decoding a loop.
    // The default value of 50 million is about 30x bigger than the number of
    private const int s2polygon_decode_max_num_vertices = 50000000;

    // Any single-vertex loop is interpreted as being either the empty loop or the
    // full loop, depending on whether the vertex is in the northern or southern
    // hemisphere respectively.

    // The single vertex in the "empty loop" vertex chain.
    private static readonly S2Point kEmptyVertex = new(0, 0, 1);

    // The single vertex in the "full loop" vertex chain.
    private static readonly S2Point kFullVertex = new(0, 0, -1);

    // A special vertex chain of length 1 that creates an empty loop (i.e., a
    // loop with no edges that contains no points).  Example usage:
    //
    //    S2Loop empty(S2Loop.kEmpty());
    //
    // The loop may be safely encoded lossily (e.g. by snapping it to an S2Cell
    // center) as long as its position does not move by 90 degrees or more.
    public static S2Loop kEmpty => new(new S2Point[] { kEmptyVertex });

    // A special vertex chain of length 1 that creates a full loop (i.e., a loop
    // with no edges that contains all points).  See kEmpty() for details.
    public static S2Loop kFull => new(new S2Point[] { kFullVertex });

    // Aux key for Dictionaries
    public static S2Loop NullLoop() => new(new S2Point[] { S2Point.Empty }, S2Debug.DISABLE);

    // Allows overriding the automatic validity checks controlled by the
    // --s2debug flag.  If this flag is true, then loops are automatically
    // checked for validity as they are initialized.  The main reason to disable
    // this flag is if you intend to call IsValid() explicitly, like this:
    //
    //   S2Loop loop;
    //   loop.set_s2debug_override(S2Debug::DISABLE);
    //   loop.Init(...);
    //   if (!loop.IsValid()) { ... }
    //
    // Without the call to set_s2debug_override(), invalid data would cause a
    // fatal error in Init() whenever the --s2debug flag is enabled.
    //
    // This setting is preserved across calls to Init() and Decode().
    public S2Debug S2DebugOverride { get; set; } = S2Debug.ALLOW;

    #endregion

    #region Constructors

    // NOTE(Alas): Merged the 2 constructors and Init
    //
    // ---
    //
    // Convenience constructor that calls Init() with the given vertices.
    //
    // ---
    //
    // Convenience constructor to disable the automatic validity checking
    // controlled by the --s2debug flag.  Example:
    //
    //   S2Loop* loop = new S2Loop(vertices, S2Debug::DISABLE);
    //
    // This is equivalent to:
    //
    //   S2Loop* loop = new S2Loop;
    //   loop->set_s2debug_override(S2Debug::DISABLE);
    //   loop->Init(vertices);
    //
    // The main reason to use this constructor is if you intend to call
    // IsValid() explicitly.  See set_s2debug_override() for details.
    //
    // ---
    //
    // Initialize a loop with given vertices.  The last vertex is implicitly
    // connected to the first.  All points should be unit length.  Loops must
    // have at least 3 vertices (except for the empty and full loops, see
    // kEmpty and kFull).  This method may be called multiple times.
    public S2Loop(IEnumerable<S2Point> vertices, S2Debug override_ = S2Debug.ALLOW)
    {
        S2DebugOverride = override_;
        ClearIndex();
        NumVertices = vertices.Count();
        Vertices = vertices.ToArray();
        InitOriginAndBound();
        InitFirstLogicalVertex();
    }

    // for Clone and Decode
    private S2Loop() { Vertices = Array.Empty<S2Point>(); }

    // Construct a loop corresponding to the given cell.
    //
    // Note that the loop and cell *do not* contain exactly the same set of
    // points, because S2Loop and S2Cell have slightly different definitions of
    // point containment.  For example, an S2Cell vertex is contained by all
    // four neighboring S2Cells, but it is contained by exactly one of four
    // S2Loops constructed from those cells.  As another example, the S2Cell
    // coverings of "cell" and "S2Loop(cell)" will be different, because the
    // loop contains points on its boundary that actually belong to other cells
    // (i.e., the covering will include a layer of neighboring cells).
    public S2Loop(S2Cell cell)
    {
        Depth = 0;
        NumVertices = 4;
        Vertices = new S2Point[NumVertices];
        S2DebugOverride = S2Debug.ALLOW;
        _unindexedContainsCalls_ = 0;

        for (int i = 0; i < 4; ++i)
        {
            Vertices[i] = cell.Vertex(i);
        }
        // We recompute the bounding rectangle ourselves, since S2Cell uses a
        // different method and we need all the bounds to be consistent.
        InitOriginAndBound();
        InitFirstLogicalVertex();
    }

    #endregion

    #region S2Loop

    public S2Point this[int index]
    {
        get
        {
            if (index < 0 || index >= Vertices.Length)
                throw new ArgumentOutOfRangeException(nameof(index));

            return Vertices[index];
        }
    }

    // Returns true if this is a valid loop.  Note that validity is checked
    // automatically during initialization when --s2debug is enabled (true by
    // default in debug binaries).
    public bool IsValid() => !FindValidationError(out _);

    // Returns true if this is *not* a valid loop and sets "error"
    // appropriately.  Otherwise returns false and leaves "error" unchanged.
    //
    // REQUIRES: error != null
    public bool FindValidationError(out S2Error error) => FindValidationErrorNoIndex(out error) || S2ShapeUtil.EdgePairs.FindSelfIntersection(_index, out error);

    // For convenience, we make two entire copies of the vertex list available:
    // vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == NumVertices.
    //
    // REQUIRES: 0 <= i < 2 * NumVertices
    public S2Point Vertex(int i)
    {
        MyDebug.Assert(i >= 0);
        MyDebug.Assert(i < (2 * NumVertices));
        int j = i - NumVertices;
        return Vertices[j < 0 ? i : j];
    }

    public S2Point[] CloneVertices()
    {
        var v = new S2Point[NumVertices];
        Array.Copy(Vertices, v, NumVertices);
        return v;
    }

    // Like vertex(), but this method returns vertices in reverse order if the
    // loop represents a polygon hole.  For example, arguments 0, 1, 2 are
    // mapped to vertices n-1, n-2, n-3, where n == NumVertices.  This
    // ensures that the interior of the polygon is always to the left of the
    // vertex chain.
    //
    // REQUIRES: 0 <= i < 2 * NumVertices
    public S2Point OrientedVertex(int i)
    {
        MyDebug.Assert(i >= 0);
        MyDebug.Assert(i < 2 * NumVertices);
        int j = i - NumVertices;
        if (j < 0) j = i;
        if (IsHole()) j = NumVertices - 1 - j;
        return Vertices[j];
    }

    // Returns an S2PointLoopSpan containing the loop vertices, for use with the
    // functions defined in s2loop_measures.h.
    public S2PointLoopSpan VerticesSpan() => new(Vertices);

    // Returns true if this is the special empty loop that contains no points.
    public bool IsEmpty() => IsEmptyOrFull() && !ContainsOrigin;

    // Returns true if this is the special full loop that contains all points.
    public bool IsFull() => IsEmptyOrFull() && ContainsOrigin;

    // Returns true if this loop is either empty or full.
    public bool IsEmptyOrFull() => NumVertices == 1;

    // Returns true if this loop represents a hole in its containing polygon.
    public bool IsHole() => (Depth & 1) != 0;

    // The sign of a loop is -1 if the loop represents a hole in its containing
    // polygon, and +1 otherwise.
    public int Sign() => IsHole() ? -1 : 1;

    // Returns true if the loop area is at most 2*Pi.  Degenerate loops are
    // handled consistently with S2Pred.Sign(), i.e., if a loop can be
    // expressed as the union of degenerate or nearly-degenerate CCW triangles,
    // then it will always be considered normalized.
    public bool IsNormalized()
    {
        // Optimization: if the longitude span is less than 180 degrees, then the
        // loop covers less than half the sphere and is therefore normalized.
        if (_bound.Lng.GetLength() < Math.PI) return true;

        return S2.IsNormalized(Vertices.ToList());
    }

    // Invert the loop if necessary so that the area enclosed by the loop is at
    // most 2*Pi.
    public void Normalize()
    {
        if (!IsNormalized()) Invert();
        MyDebug.Assert(IsNormalized());
    }

    // Reverse the order of the loop vertices, effectively complementing the
    // region represented by the loop.  For example, the loop ABCD (with edges
    // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
    // Notice that the last edge is the same in both cases except that its
    // direction has been reversed.
    public void Invert()
    {
        ClearIndex();
        if (IsEmptyOrFull())
        {
            Vertices[0] = IsFull() ? kEmptyVertex : kFullVertex;
        }
        else
        {
            Array.Reverse(Vertices, 0, NumVertices);
        }
        // origin_inside_ must be set correctly before building the S2ShapeIndex.
        ContainsOrigin ^= true;
        if (_bound.Lat.Lo > -S2.M_PI_2 && _bound.Lat.Hi < S2.M_PI_2)
        {
            // The complement of this loop contains both poles.
            _subregionBound = _bound = S2LatLngRect.Full;
        }
        else
        {
            InitBound();
        }
        InitIndex();
        InitFirstLogicalVertex();
    }

    /// <summary>
    /// Calculates firstLogicalVertex, the vertex in this loop that comes first in
    /// a total ordering of all vertices (by way of S2Point's compareTo function).
    /// </summary>
    private void InitFirstLogicalVertex()
    {
        var first = 0;
        for (var i = 1; i < NumVertices; ++i)
        {
            if (Vertex(i).CompareTo(Vertex(first)) < 0)
            {
                first = i;
            }
        }
        _firstLogicalVertex = first;
    }

    // Returns the area of the loop interior, i.e. the region on the left side of
    // the loop.  The return value is between 0 and 4*Pi.  (Note that the return
    // value is not affected by whether this loop is a "hole" or a "shell".)
    public double Area()
    {
        // S2Loop has its own convention for empty and full loops.
        if (IsEmptyOrFull())
        {
            return ContainsOrigin ? S2.M_4_PI : 0;
        }
        return S2.GetArea(new S2PointLoopSpan(Vertices));
    }

    // Returns the true centroid of the loop multiplied by the area of the loop
    // (see s2centroids.h for details on centroids).  The result is not unit
    // length, so you may want to normalize it.  Also note that in general, the
    // centroid may not be contained by the loop.
    //
    // We prescale by the loop area for two reasons: (1) it is cheaper to
    // compute this way, and (2) it makes it easier to compute the centroid of
    // more complicated shapes (by splitting them into disjoint regions and
    // adding their centroids).
    //
    // Note that the return value is not affected by whether this loop is a
    // "hole" or a "shell".
    public S2Point Centroid()
    {
        // Empty and full loops are handled correctly.
        return S2.GetCentroid(Vertices.ToList());
    }

    // Returns the geodesic curvature of the loop, defined as the sum of the turn
    // angles at each vertex (see S2.TurnAngle).  The result is positive if the
    // loop is counter-clockwise, negative if the loop is clockwise, and zero if
    // the loop is a great circle.  The geodesic curvature is equal to 2*Pi minus
    // the area of the loop.
    //
    // Degenerate and nearly-degenerate loops are handled consistently with
    // S2Pred.Sign().  So for example, if a loop has zero area (i.e., it is a
    // very small CCW loop) then its geodesic curvature will always be positive.
    public double Curvature()
    {
        // S2Loop has its own convention for empty and full loops.  For such loops,
        // we return the limit value as the area approaches 0 or 4*Pi respectively.
        if (IsEmptyOrFull())
        {
            return ContainsOrigin ? (-S2.M_2_PI) : S2.M_2_PI;
        }
        return S2.GetCurvature(Vertices.ToList());
    }

    // Returns the maximum error in GetCurvature().  The return value is not
    // constant; it depends on the loop.
    public double CurvatureMaxError()
    {
        return S2.GetCurvatureMaxError(Vertices.ToList());
    }

    // Returns the distance from the given point to the loop interior.  If the
    // loop is empty, return S1Angle.Infinity().  "x" should be unit length.
    public S1Angle Distance(S2Point x)
    {
        // Note that S2Loop.Contains(S2Point) is slightly more efficient than the
        // generic version used by S2ClosestEdgeQuery.
        if (Contains(x)) return S1Angle.Zero;
        return DistanceToBoundary(x);
    }

    // Returns the distance from the given point to the loop boundary.  If the
    // loop is empty or full, return S1Angle.Infinity() (since the loop has no
    // boundary).  "x" should be unit length.
    public S1Angle DistanceToBoundary(S2Point x)
    {
        var options = new S2ClosestEdgeQuery.Options
        {
            IncludeInteriors = false
        };
        var t = new S2ClosestEdgeQuery.PointTarget(x);
        return new S2ClosestEdgeQuery(_index, options).GetDistance(t).ToAngle();
    }

    // If the given point is contained by the loop, return it.  Otherwise return
    // the closest point on the loop boundary.  If the loop is empty, return the
    // input argument.  Note that the result may or may not be contained by the
    // loop.  "x" should be unit length.
    public S2Point Project(S2Point x)
    {
        if (Contains(x)) return x;
        return ProjectToBoundary(x);
    }

    // Returns the closest point on the loop boundary to the given point.  If the
    // loop is empty or full, return the input argument (since the loop has no
    // boundary).  "x" should be unit length.
    public S2Point ProjectToBoundary(S2Point x)
    {
        S2ClosestEdgeQuery.Options options = new();
        options.IncludeInteriors = false;
        var q = new S2ClosestEdgeQuery(_index, options);
        var target = new S2ClosestEdgeQuery.PointTarget(x);
        var edge = q.FindClosestEdge(target);
        return q.Project(x, edge);
    }

    // Returns true if the region contained by this loop is a superset of the
    // region contained by the given other loop.
    public bool Contains(S2Loop b)
    {
        // For this loop A to contains the given loop B, all of the following must
        // be true:
        //
        //  (1) There are no edge crossings between A and B except at vertices.
        //
        //  (2) At every vertex that is shared between A and B, the local edge
        //      ordering implies that A contains B.
        //
        //  (3) If there are no shared vertices, then A must contain a vertex of B
        //      and B must not contain a vertex of A.  (An arbitrary vertex may be
        //      chosen in each case.)
        //
        // The second part of (3) is necessary to detect the case of two loops whose
        // union is the entire sphere, i.e. two loops that contains each other's
        // boundaries but not each other's interiors.
        if (!_subregionBound.Contains(b._bound)) return false;

        // Special cases to handle either loop being empty or full.
        if (IsEmptyOrFull() || b.IsEmptyOrFull())
        {
            return IsFull() || b.IsEmpty();
        }

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        var relation = new ContainsRelation();
        if (HasCrossingRelation(this, b, relation)) return false;

        // There are no crossings, and if there are any shared vertices then A
        // contains B locally at each shared vertex.
        if (relation.FoundSharedVertex) return true;

        // Since there are no edge intersections or shared vertices, we just need to
        // test condition (3) above.  We can skip this test if we discovered that A
        // contains at least one point of B while checking for edge crossings.
        if (!Contains(b.Vertex(0))) return false;

        // We still need to check whether (A union B) is the entire sphere.
        // Normally this check is very cheap due to the bounding box precondition.
        if ((b._subregionBound.Contains(_bound) ||
             b._bound.Union(_bound).IsFull()) && b.Contains(Vertex(0)))
        {
            return false;
        }
        return true;
    }

    // Returns true if the region contained by this loop intersects the region
    // contained by the given other loop.
    public bool Intersects(S2Loop b)
    {
        // a.Intersects(b) if and only if !a.Complement().Contains(b).
        // This code is similar to Contains(), but is optimized for the case
        // where both loops enclose less than half of the sphere.
        if (!_bound.Intersects(b._bound)) return false;

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        var relation = new IntersectsRelation();
        if (HasCrossingRelation(this, b, relation)) return true;
        if (relation.FoundSharedVertex) return false;

        // Since there are no edge intersections or shared vertices, the loops
        // intersect only if A contains B, B contains A, or the two loops contain
        // each other's boundaries.  These checks are usually cheap because of the
        // bounding box preconditions.  Note that neither loop is empty (because of
        // the bounding box check above), so it is safe to access vertex(0).

        // Check whether A contains B, or A and B contain each other's boundaries.
        // (Note that A contains all the vertices of B in either case.)
        if (_subregionBound.Contains(b._bound) ||
            _bound.Union(b._bound).IsFull())
        {
            if (Contains(b.Vertex(0))) return true;
        }
        // Check whether B contains A.
        if (b._subregionBound.Contains(_bound))
        {
            if (b.Contains(Vertex(0))) return true;
        }
        return false;
    }

    // Returns true if two loops have the same boundary.  This is true if and
    // only if the loops have the same vertices in the same cyclic order (i.e.,
    // the vertices may be cyclically rotated).  The empty and full loops are
    // considered to have different boundaries.
    public bool BoundaryEquals(S2Loop b)
    {
        if (NumVertices != b.NumVertices) return false;

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if (IsEmptyOrFull()) return IsEmpty() == b.IsEmpty();

        for (int offset = 0; offset < NumVertices; ++offset)
        {
            if (Vertex(offset) == b.Vertex(0))
            {
                // There is at most one starting offset since loop vertices are unique.
                for (int i = 0; i < NumVertices; ++i)
                {
                    if (Vertex(i + offset) != b.Vertex(i)) return false;
                }
                return true;
            }
        }
        return false;
    }

    // Returns true if two loops have the same boundary except for vertex
    // perturbations.  More precisely, the vertices in the two loops must be in
    // the same cyclic order, and corresponding vertex pairs must be separated
    // by no more than "max_error".
    public bool BoundaryApproxEquals(S2Loop b)
    {
        return BoundaryApproxEquals(b, S1Angle.FromRadians(S2.DoubleError));
    }

    public bool BoundaryApproxEquals(S2Loop b, S1Angle max_error)
    {
        if (NumVertices != b.NumVertices) return false;

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if (IsEmptyOrFull()) return IsEmpty() == b.IsEmpty();

        for (int offset = 0; offset < NumVertices; ++offset)
        {
            if (S2.ApproxEquals(Vertex(offset), b.Vertex(0), max_error))
            {
                bool success = true;
                for (int i = 0; i < NumVertices; ++i)
                {
                    if (!S2.ApproxEquals(Vertex(i + offset), b.Vertex(i), max_error))
                    {
                        success = false;
                        break;
                    }
                }
                if (success) return true;
                // Otherwise continue looping.  There may be more than one candidate
                // starting offset since vertices are only matched approximately.
            }
        }
        return false;
    }

    // Returns true if the two loop boundaries are within "max_error" of each
    // other along their entire lengths.  The two loops may have different
    // numbers of vertices.  More precisely, this method returns true if the two
    // loops have parameterizations a:[0,1] . S^2, b:[0,1] . S^2 such that
    // distance(a(t), b(t)) <= max_error for all t.  You can think of this as
    // testing whether it is possible to drive two cars all the way around the
    // two loops such that no car ever goes backward and the cars are always
    // within "max_error" of each other.
    public bool BoundaryNear(S2Loop b)
    {
        return BoundaryNear(b, S1Angle.FromRadians(S2.DoubleError));
    }
    public bool BoundaryNear(S2Loop b, S1Angle max_error)
    {
        // Special case to handle empty or full loops.
        if (IsEmptyOrFull() || b.IsEmptyOrFull())
        {
            return (IsEmpty() && b.IsEmpty()) || (IsFull() && b.IsFull());
        }

        for (int a_offset = 0; a_offset < NumVertices; ++a_offset)
        {
            if (MatchBoundaries(this, b, a_offset, max_error)) return true;
        }
        return false;
    }

    // This method computes the oriented surface integral of some quantity f(x)
    // over the loop interior, given a function f_tri(A,B,C) that returns the
    // corresponding integral over the spherical triangle ABC.  Here "oriented
    // surface integral" means:
    //
    // (1) f_tri(A,B,C) must be the integral of f if ABC is counterclockwise,
    //     and the integral of -f if ABC is clockwise.
    //
    // (2) The result of this function is *either* the integral of f over the
    //     loop interior, or the integral of (-f) over the loop exterior.
    //
    // Note that there are at least two common situations where it easy to work
    // around property (2) above:
    //
    //  - If the integral of f over the entire sphere is zero, then it doesn't
    //    matter which case is returned because they are always equal.
    //
    //  - If f is non-negative, then it is easy to detect when the integral over
    //    the loop exterior has been returned, and the integral over the loop
    //    interior can be obtained by adding the integral of f over the entire
    //    unit sphere (a constant) to the result.
    //
    // Also requires that the default constructor for T must initialize the
    // value to zero.  (This is true for built-in types such as "double".)
    public double GetSurfaceIntegral(Func<S2Point, S2Point, S2Point, double> f_tri)
    {
        return S2.GetSurfaceIntegral(Vertices.ToList(), f_tri);
    }
    public S2Point GetSurfaceIntegral(Func<S2Point, S2Point, S2Point, S2Point> f_tri)
    {
        return S2.GetSurfaceIntegral(Vertices.ToList(), f_tri);
    }

    // Constructs a regular polygon with the given number of vertices, all
    // located on a circle of the specified radius around "center".  The radius
    // is the actual distance from "center" to each vertex.
    public static S2Loop MakeRegularLoop(S2Point center, S1Angle radius, int num_vertices)
    {
        return MakeRegularLoop(S2.GetFrame(center), radius, num_vertices);
    }

    // Like the function above, but this version constructs a loop centered
    // around the z-axis of the given coordinate frame, with the first vertex in
    // the direction of the positive x-axis.  (This allows the loop to be
    // rotated for testing purposes.)
    public static S2Loop MakeRegularLoop(S2PointS2Point frame, S1Angle radius, int num_vertices)
    {
        // We construct the loop in the given frame coordinates, with the center at
        // (0, 0, 1).  For a loop of radius "r", the loop vertices have the form
        // (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r).  The distance on the
        // sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
        double z = Math.Cos(radius.Radians);
        double r = Math.Sin(radius.Radians);
        double radian_step = S2.M_2_PI / num_vertices;
        var vertices = new List<S2Point>();
        for (int i = 0; i < num_vertices; ++i)
        {
            double angle = i * radian_step;
            var p = new S2Point(r * Math.Cos(angle), r * Math.Sin(angle), z);
            vertices.Add(S2.FromFrame(frame, p).Normalize());
        }
        return new S2Loop(vertices);
    }

    // Returns the total number of bytes used by the loop.
    public int SpaceUsed()
    {
        int size = Marshal.SizeOf(this);
        size += NumVertices * SizeHelper.SizeOf(typeof(S2Point));
        // index_ itself is already included in sizeof(*this).
        size += _index.SpaceUsed() - Marshal.SizeOf(_index);
        return size;
    }

    ////////////////////////////////////////////////////////////////////////
    // Methods intended primarily for use by the S2Polygon implementation:

    // Given two loops of a polygon, return true if A contains B.  This version
    // of Contains() is cheap because it does not test for edge intersections.
    // The loops must meet all the S2Polygon requirements; for example this
    // implies that their boundaries may not cross or have any shared edges
    // (although they may have shared vertices).
    public bool ContainsNested(S2Loop b)
    {
        if (!_subregionBound.Contains(b._bound)) return false;

        // Special cases to handle either loop being empty or full.  Also bail out
        // when B has no vertices to avoid heap overflow on the vertex(1) call
        // below.  (This method is called during polygon initialization before the
        // client has an opportunity to call IsValid.)
        if (IsEmptyOrFull() || b.NumVertices < 2)
        {
            return IsFull() || b.IsEmpty();
        }

        // We are given that A and B do not share any edges, and that either one
        // loop contains the other or they do not intersect.
        int m = FindVertex(b.Vertex(1));
        if (m < 0)
        {
            // Since b.vertex(1) is not shared, we can check whether A contains it.
            return Contains(b.Vertex(1));
        }
        // Check whether the edge order around b.vertex(1) is compatible with
        // A containing B.
        return S2WedgeRelations.WedgeContains(Vertex(m - 1), Vertex(m), Vertex(m + 1), b.Vertex(0), b.Vertex(2));
    }

    // Returns +1 if A contains the boundary of B, -1 if A excludes the boundary
    // of B, and 0 if the boundaries of A and B cross.  Shared edges are handled
    // as follows: If XY is a shared edge, define Reversed(XY) to be true if XY
    // appears in opposite directions in A and B.  Then A contains XY if and
    // only if Reversed(XY) == B.is_hole().  (Intuitively, this checks whether
    // A contains a vanishingly small region extending from the boundary of B
    // toward the interior of the polygon to which loop B belongs.)
    //
    // This method is used for testing containment and intersection of
    // multi-loop polygons.  Note that this method is not symmetric, since the
    // result depends on the direction of loop A but not on the direction of
    // loop B (in the absence of shared edges).
    //
    // REQUIRES: neither loop is empty.
    // REQUIRES: if b.IsFull, then !b.is_hole().
    public int CompareBoundary(S2Loop b)
    {
        MyDebug.Assert(!IsEmpty() && !b.IsEmpty());
        MyDebug.Assert(!b.IsFull() || !b.IsHole());

        // The bounds must intersect for containment or crossing.
        if (!_bound.Intersects(b._bound)) return -1;

        // Full loops are handled as though the loop surrounded the entire sphere.
        if (IsFull()) return 1;
        if (b.IsFull()) return -1;

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        var relation = new CompareBoundaryRelation(b.IsHole());
        if (HasCrossingRelation(this, b, relation)) return 0;
        if (relation.FoundSharedVertex)
        {
            return relation.ContainsEdge ? 1 : -1;
        }

        // There are no edge intersections or shared vertices, so we can check
        // whether A contains an arbitrary vertex of B.
        return Contains(b.Vertex(0)) ? 1 : -1;
    }

    // Given two loops whose boundaries do not cross (see CompareBoundary),
    // return true if A contains the boundary of B.  If "reverse_b" is true, the
    // boundary of B is reversed first (which only affects the result when there
    // are shared edges).  This method is cheaper than CompareBoundary() because
    // it does not test for edge intersections.
    //
    // REQUIRES: neither loop is empty.
    // REQUIRES: if b.IsFull, then reverse_b == false.
    public bool ContainsNonCrossingBoundary(S2Loop b, bool reverse_b)
    {
        MyDebug.Assert(!IsEmpty() && !b.IsEmpty());
        MyDebug.Assert(!b.IsFull() || !reverse_b);

        // The bounds must intersect for containment.
        if (!_bound.Intersects(b._bound)) return false;

        // Full loops are handled as though the loop surrounded the entire sphere.
        if (IsFull()) return true;
        if (b.IsFull()) return false;

        int m = FindVertex(b.Vertex(0));
        if (m < 0)
        {
            // Since vertex b0 is not shared, we can check whether A contains it.
            return Contains(b.Vertex(0));
        }
        // Otherwise check whether the edge (b0, b1) is contained by A.
        return WedgeContainsSemiwedge(Vertex(m - 1), Vertex(m), Vertex(m + 1),
                                      b.Vertex(1), reverse_b);
    }

    private void InitOriginAndBound()
    {
        if (NumVertices < 3)
        {
            // Check for the special empty and full loops (which have one vertex).
            if (!IsEmptyOrFull())
            {
                ContainsOrigin = false;
                _bound = new S2LatLngRect(new R1Interval(0, 0), new S1Interval(0, 0));
                return;  // Bail out without trying to access non-existent vertices.
            }
            // If the vertex is in the southern hemisphere then the loop is full,
            // otherwise it is empty.
            ContainsOrigin = Vertex(0).Z < 0;
        }
        else
        {
            // The brute force point containment algorithm works by counting edge
            // crossings starting at a fixed reference point (chosen as S2::Origin()
            // for historical reasons).  Loop initialization would be more efficient
            // if we used a loop vertex such as vertex(0) as the reference point
            // instead, however making this change would be a lot of work because
            // origin_inside_ is currently part of the Encode() format.
            //
            // In any case, we initialize origin_inside_ by first guessing that it is
            // outside, and then seeing whether we get the correct containment result
            // for vertex 1.  If the result is incorrect, the origin must be inside
            // the loop instead.  Note that the S2Loop is not necessarily valid and so
            // we need to check the requirements of S2::AngleContainsVertex() first.
            bool v1_inside = Vertex(0) != Vertex(1) && Vertex(2) != Vertex(1) &&
                             S2.AngleContainsVertex(Vertex(0), Vertex(1), Vertex(2));
            ContainsOrigin = false;  // Initialize before calling Contains().

            // Note that Contains(S2Point) only does a bounds check once InitIndex()
            // has been called, so it doesn't matter that bound_ is undefined here.
            if (v1_inside != Contains(Vertex(1)))
            {
                ContainsOrigin = true;
            }
        }
        // We *must* call InitBound() before InitIndex(), because InitBound() calls
        // Contains(S2Point), and Contains(S2Point) does a bounds check whenever the
        // index is not fresh (i.e., the loop has been added to the index but the
        // index has not been updated yet).
        //
        // TODO(ericv): When fewer S2Loop methods depend on internal bounds checks,
        // consider computing the bound on demand as well.
        InitBound();
        InitIndex();
    }
    private void InitBound()
    {
        // Check for the special empty and full loops.
        if (IsEmptyOrFull())
        {
            if (IsEmpty())
            {
                _subregionBound = _bound = S2LatLngRect.Empty;
            }
            else
            {
                _subregionBound = _bound = S2LatLngRect.Full;
            }
            return;
        }

        // The bounding rectangle of a loop is not necessarily the same as the
        // bounding rectangle of its vertices.  First, the maximal latitude may be
        // attained along the interior of an edge.  Second, the loop may wrap
        // entirely around the sphere (e.g. a loop that defines two revolutions of a
        // candy-cane stripe).  Third, the loop may include one or both poles.
        // Note that a small clockwise loop near the equator contains both poles.

        S2LatLngRectBounder bounder = new();
        for (int i = 0; i < NumVertices; ++i)
        {
            bounder.AddPoint(Vertex(i));
        }
        var b = bounder.GetBound();
        if (Contains(new S2Point(0, 0, 1)))
        {
            b = new S2LatLngRect(new R1Interval(b.Lat.Lo, S2.M_PI_2), S1Interval.Full);
        }
        // If a loop contains the south pole, then either it wraps entirely
        // around the sphere (full longitude range), or it also contains the
        // north pole in which case b.lng().IsFull due to the test above.
        // Either way, we only need to do the south pole containment test if
        // b.lng().IsFull.
        if (b.Lng.IsFull() && Contains(new S2Point(0, 0, -1)))
        {
            b = new S2LatLngRect(new R1Interval(-S2.M_PI_2, b.Lat.Hi), b.Lng);
        }
        _bound = b;
        _subregionBound = S2LatLngRectBounder.ExpandForSubregions(_bound);
    }

    private void InitIndex()
    {
        _index.Add(new Shape(this));
#if s2loop_not_lazy_indexing
        index_.ForceBuild();
#endif
#if s2debug
        if (S2DebugOverride == S2Debug.ALLOW)
        {
            // Note that s2debug is false in optimized builds (by default).
            MyDebug.Assert(IsValid());
        }
#endif
    }

    // A version of Contains(S2Point) that does not use the S2ShapeIndex.
    // Used by the S2Polygon implementation.
    public bool BruteForceContains(S2Point p)
    {
        // Empty and full loops don't need a special case, but invalid loops with
        // zero vertices do, so we might as well handle them all at once.
        if (NumVertices < 3) return ContainsOrigin;

        var crosser = new S2EdgeCrosser(S2.Origin, p, Vertex(0));
        bool inside = ContainsOrigin;
        for (int i = 1; i <= NumVertices; ++i)
        {
            inside ^= crosser.EdgeOrVertexCrossing(Vertex(i));
        }
        return inside;
    }

    // Like FindValidationError(), but skips any checks that would require
    // building the S2ShapeIndex (i.e., self-intersection tests).  This is used
    // by the S2Polygon implementation, which uses its own index to check for
    // loop self-intersections.
    public bool FindValidationErrorNoIndex(out S2Error error)
    {
        // subregion_bound_ must be at least as large as bound_.  (This is an
        // internal consistency check rather than a test of client data.)
        MyDebug.Assert(_subregionBound.Contains(_bound));

        // All vertices must be unit length.  (Unfortunately this check happens too
        // late in debug mode, because S2Loop construction calls S2Pred.Sign which
        // expects vertices to be unit length.  But it is still a useful check in
        // optimized builds.)
        for (int i = 0; i < NumVertices; ++i)
        {
            if (!Vertex(i).IsUnitLength())
            {
                error = new(S2ErrorCode.NOT_UNIT_LENGTH, $"Vertex {i} is not unit length");
                return true;
            }
        }
        // Loops must have at least 3 vertices (except for the empty and full loops).
        if (NumVertices < 3)
        {
            if (IsEmptyOrFull())
            {
                error = S2Error.OK;
                return false;  // Skip remaining tests.
            }
            error = new(S2ErrorCode.LOOP_NOT_ENOUGH_VERTICES, "Non-empty, non-full loops must have at least 3 vertices");
            return true;
        }
        // Loops are not allowed to have any duplicate vertices or edge crossings.
        // We split this check into two parts.  First we check that no edge is
        // degenerate (identical endpoints).  Then we check that there are no
        // intersections between non-adjacent edges (including at vertices).  The
        // second part needs the S2ShapeIndex, so it does not fall within the scope
        // of this method.
        for (int i = 0; i < NumVertices; ++i)
        {
            if (Vertex(i) == Vertex(i + 1))
            {
                error = new(S2ErrorCode.DUPLICATE_VERTICES, $"Edge {i} is degenerate (duplicate vertex)");
                return true;
            }
            if (Vertex(i) == -Vertex(i + 1))
            {
                error = new(S2ErrorCode.ANTIPODAL_VERTICES, $"Vertices {i} and {(i + 1) % NumVertices} are antipodal");
                return true;
            }
        }
        error = S2Error.OK;
        return false;
    }

    // Converts the loop vertices to the S2XYZFaceSiTi format and store the result
    // in the given array, which must be large enough to store all the vertices.
    public void GetXYZFaceSiTiVertices(S2PointCompression.S2XYZFaceSiTi[] vertices, int offset)
    {
        for (int i = 0; i < NumVertices; ++i)
        {
            var xyz = Vertex(i);
            var cellLevel = S2.XYZtoFaceSiTi(xyz, out int face, out uint si, out uint ti);
            vertices[offset + i] = new S2PointCompression.S2XYZFaceSiTi(xyz, face, si, ti, cellLevel);
        }
    }

    // Encode the loop's vertices using S2EncodePointsCompressed.  Uses
    // approximately 8 bytes for the first vertex, going down to less than 4 bytes
    // per vertex on Google's geographic repository, plus 24 bytes per vertex that
    // does not correspond to the center of a cell at level 'snap_level'. The loop
    // vertices must first be converted to the S2XYZFaceSiTi format with
    // GetXYZFaceSiTiVertices.
    //
    // REQUIRES: the loop is initialized and valid.
    public void EncodeCompressed(Encoder encoder, S2PointCompression.S2XYZFaceSiTi[] vertices, int offset, int snap_level)
    {
        // Ensure enough for the data we write before S2EncodePointsCompressed.
        // S2EncodePointsCompressed ensures its space.
        encoder.Ensure(Encoder.kVarintMax32);
        encoder.PutVarUInt32((uint)NumVertices);

        S2PointCompression.S2EncodePointsCompressed(vertices.Skip(offset).Take(NumVertices).ToArray(), snap_level, encoder);

        var properties = GetCompressedEncodingProperties();

        // Ensure enough only for what we write.  Let the bound ensure its own
        // space.
        encoder.Ensure(2 * Encoder.kVarintMax32);
        encoder.PutVarUInt32((uint)properties);
        encoder.PutVarUInt32((uint)Depth);
        if (BitsUtils.IsBitSet((uint)properties, (int)CompressedLoopProperty.kBoundEncoded))
        {
            (_bound as IEncoder).Encode(encoder);
        }
        MyDebug.Assert(encoder.Avail() >= 0);
    }

    // Decode a loop encoded with EncodeCompressed. The parameters must be the
    // same as the one used when EncodeCompressed was called.
    public static (bool success, S2Loop? output) DecodeCompressed(Decoder decoder, int snap_level)
    {
        // TryGetVarUInt32 takes a UInt32*, but NumVertices is signed.
        // Decode to a temporary variable to avoid reinterpret_cast.
        if (!decoder.TryGetVarUInt32(out var unsigned_num_vertices))
        {
            return (false, null);
        }
        if (unsigned_num_vertices == 0 ||
            unsigned_num_vertices > s2polygon_decode_max_num_vertices)
        {
            return (false, null);
        }

        var vertices = new S2Point[unsigned_num_vertices];

        if (!S2PointCompression.S2DecodePointsCompressed(decoder, snap_level, vertices, 0))
        {
            return (false, null);
        }
        if (!decoder.TryGetVarUInt32(out var properties_UInt32))
        {
            return (false, null);
        }
        var containsOrigin = BitsUtils.IsBitSet(properties_UInt32, (int)CompressedLoopProperty.kOriginInside);

        if (!decoder.TryGetVarUInt32(out var unsigned_depth))
        {
            return (false, null);
        }
        var depth = (int)unsigned_depth;

        var initBoundFlag = false;
        S2LatLngRect bound = S2LatLngRect.Empty;
        S2LatLngRect subregionBound = S2LatLngRect.Empty;
        if (BitsUtils.IsBitSet(properties_UInt32, (int)CompressedLoopProperty.kBoundEncoded))
        {
            var (success_bound, bound2) = S2LatLngRect.Decode(decoder);
            if (!success_bound)
            {
                return (false, null);
            }
            bound = bound2;
            subregionBound = S2LatLngRectBounder.ExpandForSubregions(bound);
        }
        else
        {
            initBoundFlag = true;
        }

        var output = new S2Loop
        {
            NumVertices = (int)unsigned_num_vertices,
            Vertices = vertices,
            ContainsOrigin = containsOrigin,
            Depth = depth,
            _bound = bound,
            _subregionBound = subregionBound,
        };
        if (initBoundFlag)
        {
            output.InitBound();
        }

        output.InitIndex();
        return (true, output);
    }

    // Returns a bitset of properties used by EncodeCompressed
    // to efficiently encode boolean values.  Properties are
    // origin_inside and whether the bound was encoded.
    private int GetCompressedEncodingProperties()
    {
        var properties = 0;
        if (ContainsOrigin)
        {
            properties = BitsUtils.SetBit(properties, (int)CompressedLoopProperty.kOriginInside);
        }

        // Write whether there is a bound so we can change the threshold later.
        // Recomputing the bound multiplies the decode time taken per vertex
        // by a factor of about 3.5.  Without recomputing the bound, decode
        // takes approximately 125 ns / vertex.  A loop with 63 vertices
        // encoded without the bound will take ~30us to decode, which is
        // acceptable.  At ~3.5 bytes / vertex without the bound, adding
        // the bound will increase the size by <15%, which is also acceptable.
        const int kMinVerticesForBound = 64;
        if (NumVertices >= kMinVerticesForBound)
        {
            properties = BitsUtils.SetBit(properties, (int)CompressedLoopProperty.kBoundEncoded);
        }
        return properties & 3;
    }

    // Given an iterator that is already positioned at the S2ShapeIndexCell
    // containing "p", returns Contains(p).
    private bool Contains(MutableS2ShapeIndex.Enumerator it, S2Point p)
    {
        // Test containment by drawing a line segment from the cell center to the
        // given point and counting edge crossings.
        S2ClippedShape a_clipped = it.Cell.Clipped(0);
        bool inside = a_clipped.ContainsCenter;
        int a_num_edges = a_clipped.NumEdges;
        if (a_num_edges > 0)
        {
            S2Point center = it.Center;
            var crosser = new S2EdgeCrosser(center, p);
            int ai_prev = -2;
            for (int i = 0; i < a_num_edges; ++i)
            {
                int ai = a_clipped.Edges[i];
                if (ai != ai_prev + 1) crosser.RestartAt(Vertex(ai));
                ai_prev = ai;
                inside ^= crosser.EdgeOrVertexCrossing(Vertex(ai + 1));
            }
        }
        return inside;
    }

    // Returns true if the loop boundary intersects "target".  It may also
    // return true when the loop boundary does not intersect "target" but
    // some edge comes within the worst-case error tolerance.
    //
    // REQUIRES: it.id().contains(target.id())
    // [This condition is true whenever it.Locate(target) returns INDEXED.]
    private bool BoundaryApproxIntersects(MutableS2ShapeIndex.Enumerator it, S2Cell target)
    {
        MyDebug.Assert(it.Id.Contains(target.Id));
        S2ClippedShape a_clipped = it.Cell.Clipped(0);
        int a_num_edges = a_clipped.NumEdges;

        // If there are no edges, there is no intersection.
        if (a_num_edges == 0) return false;

        // We can save some work if "target" is the index cell itself.
        if (it.Id == target.Id) return true;

        // Otherwise check whether any of the edges intersect "target".
        var kMaxError = S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kIntersectsRectErrorUVDist;
        R2Rect bound = target.BoundUV.Expanded(kMaxError);
        for (int i = 0; i < a_num_edges; ++i)
        {
            int ai = a_clipped.Edges[i];
            var clip = S2EdgeClipping.ClipToPaddedFace(Vertex(ai), Vertex(ai + 1), target.Face, kMaxError, out R2Point v0, out R2Point v1);
            if (clip && S2EdgeClipping.IntersectsRect(v0, v1, bound))
            {
                return true;
            }
        }
        return false;
    }

    // Returns an index "first" and a direction "dir" such that the vertex
    // sequence (first, first + dir, ..., first + (n - 1) * dir) does not change
    // when the loop vertex order is rotated or reversed.  This allows the loop
    // vertices to be traversed in a canonical order.
    public S2.LoopOrder GetCanonicalLoopOrder()
    {
        return S2.GetCanonicalLoopOrder(Vertices.ToList());
    }

    // Returns the index of a vertex at point "p", or -1 if not found.
    // The return value is in the range 1..NumVertices if found.
    private int FindVertex(S2Point p)
    {
        if (NumVertices < 10)
        {
            // Exhaustive search.  Return value must be in the range [1..N].
            for (int i = 1; i <= NumVertices; ++i)
            {
                if (Vertex(i) == p) return i;
            }
            return -1;
        }
        MutableS2ShapeIndex.Enumerator it = new(_index);
        if (!it.Locate(p)) return -1;

        S2ClippedShape a_clipped = it.Cell.Clipped(0);
        for (int i = a_clipped.NumEdges - 1; i >= 0; --i)
        {
            int ai = a_clipped.Edges[i];
            // Return value must be in the range [1..N].
            if (Vertex(ai) == p) return (ai == 0) ? NumVertices : ai;
            if (Vertex(ai + 1) == p) return ai + 1;
        }
        return -1;
    }

    // This method checks all edges of loop A for intersection against all edges
    // of loop B.  If there is any shared vertex, the wedges centered at this
    // vertex are sent to "relation".
    //
    // If the two loop boundaries cross, this method is guaranteed to return
    // true.  It also returns true in certain cases if the loop relationship is
    // equivalent to crossing.  For example, if the relation is Contains() and a
    // point P is found such that B contains P but A does not contain P, this
    // method will return true to indicate that the result is the same as though
    // a pair of crossing edges were found (since Contains() returns false in
    // both cases).
    //
    // See Contains(), Intersects() and CompareBoundary() for the three uses of
    // this function.
    private static bool HasCrossingRelation(S2Loop a, S2Loop b, LoopRelation relation)
    {
        // We look for S2CellId ranges where the indexes of A and B overlap, and
        // then test those edges for crossings.
        var ai = new RangeEnumerator(a._index);
        var bi = new RangeEnumerator(b._index);
        var ab = new LoopCrosser(a, b, relation, false);  // Tests edges of A against B
        var ba = new LoopCrosser(b, a, relation, true);   // Tests edges of B against A
        while (!ai.Done() || !bi.Done())
        {
            if (ai.RangeMax < bi.RangeMin)
            {
                // The A and B cells don't overlap, and A precedes B.
                ai.SeekTo(bi);
            }
            else if (bi.RangeMax < ai.RangeMin)
            {
                // The A and B cells don't overlap, and B precedes A.
                bi.SeekTo(ai);
            }
            else
            {
                // One cell contains the other.  Determine which cell is larger.
                var ab_relation = (long)(ai.Id.LowestOnBit() - bi.Id.LowestOnBit());
                if (ab_relation > 0)
                {
                    // A's index cell is larger.
                    if (ab.HasCrossingRelation(ai, bi)) return true;
                }
                else if (ab_relation < 0)
                {
                    // B's index cell is larger.
                    if (ba.HasCrossingRelation(bi, ai)) return true;
                }
                else
                {
                    // The A and B cells are the same.  Since the two cells have the same
                    // center point P, check whether P satisfies the crossing targets.
                    if (ai.ContainsCenter() == (ab.ACrossingTarget != 0) &&
                        bi.ContainsCenter() == (ab.BCrossingTarget != 0))
                    {
                        return true;
                    }
                    // Otherwise test all the edge crossings directly.
                    if (ai.NumEdges() > 0 && bi.NumEdges() > 0 &&
                        ab.CellCrossesCell(ai.Clipped(), bi.Clipped()))
                    {
                        return true;
                    }
                    ai.MoveNext();
                    bi.MoveNext();
                }
            }
        }
        return false;
    }

    // When the loop is modified (Invert(), or Init() called again) then the
    // indexing structures need to be cleared since they become invalid.
    private void ClearIndex()
    {
        _unindexedContainsCalls_ = 0;
        _index.Clear();
    }

    // Returns true if the wedge (a0, ab1, a2) contains the "semiwedge" defined as
    // any non-empty open set of rays immediately CCW from the edge (ab1, b2).  If
    // "reverse_b" is true, then substitute "clockwise" for "CCW"; this simulates
    // what would happen if the direction of loop B was reversed.
    private static bool WedgeContainsSemiwedge(S2Point a0, S2Point ab1, S2Point a2, S2Point b2, bool reverse_b)
    {
        if (b2 == a0 || b2 == a2)
        {
            // We have a shared or reversed edge.
            return (b2 == a0) == reverse_b;
        }
        else
        {
            return S2Pred.OrderedCCW(a0, a2, b2, ab1);
        }
    }

    private static bool MatchBoundaries(S2Loop a, S2Loop b, int a_offset, S1Angle max_error)
    {
        // The state consists of a pair (i,j).  A state transition consists of
        // incrementing either "i" or "j".  "i" can be incremented only if
        // a(i+1+a_offset) is near the edge from b(j) to b(j+1), and a similar rule
        // applies to "j".  The function returns true iff we can proceed all the way
        // around both loops in this way.
        //
        // Note that when "i" and "j" can both be incremented, sometimes only one
        // choice leads to a solution.  We handle this using a stack and
        // backtracking.  We also keep track of which states have already been
        // explored to avoid duplicating work.

        var pending = new List<(int i, int j)>();
        var done = new SortedSet<(int i, int j)>();
        pending.Add((0, 0));
        while (pending.Any())
        {
            var (i, j) = pending.Last();
            pending.RemoveAt(pending.Count - 1);
            if (i == a.NumVertices && j == b.NumVertices)
            {
                return true;
            }
            done.Add((i, j));

            // If (i == na && offset == na-1) where na == a.NumVertices, then
            // then (i+1+offset) overflows the [0, 2*na-1] range allowed by vertex().
            // So we reduce the range if necessary.
            int io = i + a_offset;
            if (io >= a.NumVertices) io -= a.NumVertices;

            if (i < a.NumVertices && !done.Contains((i + 1, j)) &&
                S2.GetDistance(a.Vertex(io + 1), b.Vertex(j),
                                        b.Vertex(j + 1)) <= max_error)
            {
                pending.Add((i + 1, j));
            }
            if (j < b.NumVertices && !done.Contains((i, j + 1)) &&
                S2.GetDistance(b.Vertex(j + 1), a.Vertex(io),
                                        a.Vertex(io + 1)) <= max_error)
            {
                pending.Add((i, j + 1));
            }
        }
        return false;
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    // GetRectBound() returns essentially tight results, while GetCapBound()
    // might have a lot of extra padding.  Both bounds are conservative in that
    // if the loop contains a point P, then the bound contains P also.
    public S2Cap GetCapBound()
    {
        return _bound.GetCapBound();
    }
    public S2LatLngRect GetRectBound()
    {
        return _bound;
    }
    public bool Contains(S2Cell target)
    {
        MutableS2ShapeIndex.Enumerator it = new(_index);
        S2CellRelation relation = it.Locate(target.Id);

        // If "target" is disjoint from all index cells, it is not contained.
        // Similarly, if "target" is subdivided into one or more index cells then it
        // is not contained, since index cells are subdivided only if they (nearly)
        // intersect a sufficient number of edges.  (But note that if "target" itself
        // is an index cell then it may be contained, since it could be a cell with
        // no edges in the loop interior.)
        if (relation != S2CellRelation.INDEXED) return false;

        // Otherwise check if any edges intersect "target".
        if (BoundaryApproxIntersects(it, target)) return false;

        // Otherwise check if the loop contains the center of "target".
        return Contains(it, target.Center());
    }
    public bool MayIntersect(S2Cell target)
    {
        MutableS2ShapeIndex.Enumerator it = new(_index);
        S2CellRelation relation = it.Locate(target.Id);

        // If "target" does not overlap any index cell, there is no intersection.
        if (relation == S2CellRelation.DISJOINT) return false;

        // If "target" is subdivided into one or more index cells, there is an
        // intersection to within the S2ShapeIndex error bound (see Contains).
        if (relation == S2CellRelation.SUBDIVIDED) return true;

        // If "target" is an index cell, there is an intersection because index cells
        // are created only if they have at least one edge or they are entirely
        // contained by the loop.
        if (it.Id == target.Id) return true;

        // Otherwise check if any edges intersect "target".
        if (BoundaryApproxIntersects(it, target)) return true;

        // Otherwise check if the loop contains the center of "target".
        return Contains(it, target.Center());
    }
    // The point 'p' does not need to be normalized.
    public bool Contains(S2Point p)
    {
        // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
        // worthwhile only when it might allow us to delay building the index.
        if (!_index.IsFresh() && !_bound.Contains(p)) return false;

        // For small loops it is faster to just check all the crossings.  We also
        // use this method during loop initialization because InitOriginAndBound()
        // calls Contains() before InitIndex().  Otherwise, we keep track of the
        // number of calls to Contains() and only build the index when enough calls
        // have been made so that we think it is worth the effort.  Note that the
        // code below is structured so that if many calls are made in parallel only
        // one thread builds the index, while the rest continue using brute force
        // until the index is actually available.
        //
        // The constants below were tuned using the benchmarks.  It turns out that
        // building the index costs roughly 50x as much as Contains().  (The ratio
        // increases slowly from 46x with 64 points to 61x with 256k points.)  The
        // textbook approach to this problem would be to wait until the cumulative
        // time we would have saved with an index approximately equals the cost of
        // building the index, and then build it.  (This gives the optimal
        // competitive ratio of 2; look up "competitive algorithms" for details.)
        // We set the limit somewhat lower than this (20 rather than 50) because
        // building the index may be forced anyway by other API calls, and so we
        // want to err on the side of building it too early.

        const int kMaxBruteForceVertices = 32;
        const int kMaxUnindexedContainsCalls = 20;  // See notes above.
        if (_index.NumShapeIds() == 0 ||  // InitIndex() not called yet
            NumVertices <= kMaxBruteForceVertices ||
            (!_index.IsFresh() &&
             ++_unindexedContainsCalls_ != kMaxUnindexedContainsCalls))
        {
            return BruteForceContains(p);
        }
        // Otherwise we look up the S2ShapeIndex cell containing this point.  Note
        // the index is built automatically the first time an iterator is created.
        MutableS2ShapeIndex.Enumerator it = new(_index);
        if (!it.Locate(p)) return false;

        return Contains(it, p);
    }

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2Loop to "encoder".
    //
    // Generally clients should not use S2Loop.Encode().  Instead they should
    // encode an S2Polygon, which unlike this method supports (lossless)
    // compression.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(NumVertices * SizeHelper.SizeOf(typeof(S2Point)) + 20);  // sufficient

        encoder.Put8(S2.kCurrentLosslessEncodingVersionNumber);
        encoder.Put32(NumVertices);
        encoder.PutPoints(Vertices);
        encoder.Put8((byte)(ContainsOrigin ? 1 : 0));
        encoder.Put32(Depth);
        MyDebug.Assert(encoder.Avail() >= 0);

        _bound.Encode(encoder, hint);
    }

    // Decodes a loop encoded with Encode() or the private method
    // EncodeCompressed() (used by the S2Polygon encoder).  Returns true on
    // success.
    //
    // This method may be called with loops that have already been initialized.
    public static (bool, S2Loop?) Decode(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte)) return (false, null);
        byte version = decoder.Get8();
        if (version != S2.kCurrentLosslessEncodingVersionNumber)
            return (false, null);

        // Perform all checks before modifying vertex state. Empty loops are
        // explicitly allowed here: a newly created loop has zero vertices
        // and such loops encode and decode properly.
        if (decoder.Avail() < sizeof(UInt32)) return (false, null);
        var num_vertices = (int)decoder.Get32();
        if (num_vertices > s2polygon_decode_max_num_vertices)
        {
            return (false, null);
        }

        var pvSize = SizeHelper.SizeOf(typeof(S2Point));

        if (decoder.Avail() < (num_vertices * pvSize + sizeof(byte) + sizeof(UInt32)))
        {
            return (false, null);
        }

        var newVertices = new S2Point[num_vertices];
        decoder.GetPoints(newVertices, 0, num_vertices);
        var originInside = decoder.Get8() != 0;
        var depth_ = (int)decoder.Get32();
        var (success_bound, bound) = S2LatLngRect.Decode(decoder);
        if (!success_bound) return (false, null);

        var subregionBound = S2LatLngRectBounder.ExpandForSubregions(bound);
        var loop = new S2Loop()
        {
            NumVertices = num_vertices,
            Vertices = newVertices,
            ContainsOrigin = originInside,
            Depth = depth_,
            _bound = bound,
            _subregionBound = subregionBound,
        };
        loop.InitFirstLogicalVertex();

        // An initialized loop will have some non-zero count of vertices. A default
        // (uninitialized) has zero vertices. This code supports encoding and
        // decoding of uninitialized loops, but we only want to call InitIndex for
        // initialized loops. Otherwise we defer InitIndex until the call to Init().
        if (num_vertices > 0)
        {
            loop.InitIndex();
        }
        return (true, loop);
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone()
        => new S2Loop
        {
            Depth = Depth,
            NumVertices = NumVertices,
            ContainsOrigin = ContainsOrigin,
            _unindexedContainsCalls_ = 0,
            _bound = _bound,
            _subregionBound = _subregionBound,
            Vertices = (S2Point[])Vertices.Clone(),
            _index = _index,
            _firstLogicalVertex = _firstLogicalVertex,
        };

    #endregion

    #region IEquatable

    // Returns true if two loops have the same vertices in the same linear order
    // (i.e., cyclic rotations are not allowed).

    public bool Equals(S2Loop? b) => b is not null && Vertices.SequenceEqual(b.Vertices);

    public override int GetHashCode() => LinqUtils.GetSequenceHashCode(Vertices);

    #endregion

    #region IComparable

    public int CompareTo(S2Loop? other)
    {
        if (other is null) return 1;

        if (NumVertices != other.NumVertices)
        {
            return NumVertices - other.NumVertices;
        }
        // Compare the two loops' vertices, starting with each loop's
        // firstLogicalVertex. This allows us to always catch cases where logically
        // identical loops have different vertex orderings (e.g. ABCD and BCDA).
        var maxVertices = NumVertices;
        var iThis = _firstLogicalVertex;
        var iOther = other._firstLogicalVertex;
        for (var i = 0; i < maxVertices; ++i, ++iThis, ++iOther)
        {
            var compare = Vertex(iThis).CompareTo(other.Vertex(iOther));
            if (compare != 0)
            {
                return compare;
            }
        }
        return 0;
    }

    public static bool operator <(S2Loop x, S2Loop y) => x.CompareTo(y) < 0;
    public static bool operator >(S2Loop x, S2Loop y) => x.CompareTo(y) > 0;
    public static bool operator <=(S2Loop x, S2Loop y) => x.CompareTo(y) <= 0;
    public static bool operator >=(S2Loop x, S2Loop y) => x.CompareTo(y) >= 0;

    #endregion

    #region Object

    public override string ToString()
    {
        var sb = new System.Text.StringBuilder("S2Loop, ");

        sb.Append(Vertices.Length).Append(" points. [");

        foreach (var v in Vertices)
        {
            sb.Append(v.ToString()).Append(' ');
        }
        sb.Append(']');

        return sb.ToString();
    }

    #endregion

    //    "The upper limit on the number of loops that are allowed by the "
    //    "S2Polygon.Decode method.");

    // Boolean properties for compressed loops.
    // See GetCompressedEncodingProperties.
    private enum CompressedLoopProperty
    {
        kOriginInside,
        kBoundEncoded,
        kNumProperties
    }

    // Wrapper class for indexing a loop (see S2ShapeIndex).  Once this object
    // is inserted into an S2ShapeIndex it is owned by that index, and will be
    // automatically deleted when no longer needed by the index.  Note that this
    // class does not take ownership of the loop itself (see OwningShape below).
    // You can also subtype this class to store additional data (see S2Shape for
    // details).
    public class Shape : S2Shape
    {
        public Shape(S2Loop loop) { Loop = loop; }

        public S2Loop Loop { get; private set; }

        // S2Shape interface:
        public override int NumEdges()
        {
            return Loop.IsEmptyOrFull() ? 0 : Loop.NumVertices;
        }

        public override Edge GetEdge(int e)
        {
            return new Edge(Loop.Vertex(e), Loop.Vertex(e + 1));
        }

        public override int Dimension()
        {
            return 2;
        }

        public override ReferencePoint GetReferencePoint()
        {
            return new ReferencePoint(S2.Origin, Loop.ContainsOrigin);
        }

        public override int NumChains()
        {
            return Loop.IsEmpty() ? 0 : 1;
        }

        public override Chain GetChain(int i)
        {
            MyDebug.Assert(i == 0);
            return new Chain(0, NumEdges());
        }

        public override Edge ChainEdge(int i, int j)
        {
            MyDebug.Assert(i == 0);
            return new Edge(Loop.Vertex(j), Loop.Vertex(j + 1));
        }

        public override ChainPosition GetChainPosition(int e)
        {
            return new ChainPosition(0, e);
        }
    }

    // LoopRelation is an abstract class that defines a relationship between two
    // loops (Contains, Intersects, or CompareBoundary).
    internal abstract class LoopRelation
    {
        public LoopRelation() { }

        // Optionally, a_target() and b_target() can specify an early-exit condition
        // for the loop relation.  If any point P is found such that
        //
        //   A.Contains(P) == a_crossing_target() &&
        //   B.Contains(P) == b_crossing_target()
        //
        // then the loop relation is assumed to be the same as if a pair of crossing
        // edges were found.  For example, the Contains() relation has
        //
        //   a_crossing_target() == 0
        //   b_crossing_target() == 1
        //
        // because if A.Contains(P) == 0 (false) and B.Contains(P) == 1 (true) for
        // any point P, then it is equivalent to finding an edge crossing (i.e.,
        // since Contains() returns false in both cases).
        //
        // Loop relations that do not have an early-exit condition of this form
        // should return -1 for both crossing targets.
        public abstract int ACrossingTarget();
        public abstract int BCrossingTarget();

        // Given a vertex "ab1" that is shared between the two loops, return true if
        // the two associated wedges (a0, ab1, a2) and (b0, ab1, b2) are equivalent
        // to an edge crossing.  The loop relation is also allowed to maintain its
        // own internal state, and can return true if it observes any sequence of
        // wedges that are equivalent to an edge crossing.
        public abstract bool WedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2);
    }

    // RangeIterator is a wrapper over MutableS2ShapeIndex::Iterator with extra
    // methods that are useful for merging the contents of two or more
    // S2ShapeIndexes.
    public class RangeEnumerator : IReversableEnumerator<S2ShapeIndexCell>
    {
        private readonly S2ShapeIndex _index;
        private readonly MutableS2ShapeIndex.Enumerator _it;
        // The min and max leaf cell ids covered by the current cell.  If done() is
        // true, these methods return a value larger than any valid cell id.
        public S2CellId RangeMin { get; private set; }
        public S2CellId RangeMax { get; private set; }

        // Construct a new RangeIterator positioned at the first cell of the index.
        public RangeEnumerator(MutableS2ShapeIndex index)
        {
            _index = index;
            _it = new(index, S2ShapeIndex.InitialPosition.BEGIN);
            Refresh();
        }

        // The current S2CellId and cell contents.
        public S2CellId Id => _it.Id;
        public S2ShapeIndexCell Cell => _it.Cell;

        // Various other convenience methods for the current cell.
        public S2ClippedShape Clipped() { return Cell.Clipped(0); }
        public int NumEdges() { return Clipped().NumEdges; }
        public bool ContainsCenter() { return Clipped().ContainsCenter; }

        public S2ShapeIndexCell Current => _it.Current;
        object System.Collections.IEnumerator.Current => _it.Current;

        public bool MoveNext() { var has = _it.MoveNext(); Refresh(); return has; }
        public bool MovePrevious() { var has = _it.MovePrevious(); Refresh(); return has; }
        public void Reset() => _it.Reset();
        public bool Done() => _it.Done();
        public void SetPosition(int position) => _it.SetPosition(position);
        public void Dispose() { GC.SuppressFinalize(this); }

        // Position the iterator at the first cell that overlaps or follows
        // "target", i.e. such that RangeMax >= target.RangeMin.
        public void SeekTo(RangeEnumerator target)
        {
            _it.SetPosition(_index.SeekCell(target.RangeMin).pos);
            // If the current cell does not overlap "target", it is possible that the
            // previous cell is the one we are looking for.  This can only happen when
            // the previous cell contains "target" but has a smaller S2CellId.
            if (_it.Done() || _it.Id.RangeMin() > target.RangeMax)
            {
                if (_it.MovePrevious() && _it.Id.RangeMax() < target.Id) _it.MoveNext();
            }
            Refresh();
        }

        // Position the iterator at the first cell that follows "target", i.e. the
        // first cell such that RangeMin > target.RangeMax.
        public void SeekBeyond(RangeEnumerator target)
        {
            _it.SetPosition(_index.SeekCell(target.RangeMax.Next()).pos);
            if (!_it.Done() && _it.Id.RangeMin() <= target.RangeMax)
            {
                _it.MoveNext();
            }
            Refresh();
        }

        // Updates internal state after the iterator has been repositioned.
        private void Refresh()
        {
            RangeMin = Id.RangeMin();
            RangeMax = Id.RangeMax();
        }
    }

    // LoopCrosser is a helper class for determining whether two loops cross.
    // It is instantiated twice for each pair of loops to be tested, once for the
    // pair (A,B) and once for the pair (B,A), in order to be able to process
    // edges in either loop nesting order.
    internal class LoopCrosser
    {
        #region Fields, Constants

        // Return the crossing targets for the loop relation, taking into account
        // whether the loops have been swapped.
        public readonly int ACrossingTarget;
        public readonly int BCrossingTarget;

        private readonly S2Loop a_;
        private readonly S2Loop b_;
        private readonly LoopRelation relation_;
        private readonly bool swapped_;

        // State maintained by StartEdge() and EdgeCrossesCell().
        private readonly S2EdgeCrosser crosser_ = new();
        private int aj_, bj_prev_;

        // Temporary data declared here to avoid repeated memory allocations.
        private readonly S2CrossingEdgeQuery b_query_;
        private readonly List<S2ShapeIndexCell> b_cells_ = new();

        #endregion

        #region Constructors

        // If "swapped" is true, the loops A and B have been swapped.  This affects
        // how arguments are passed to the given loop relation, since for example
        // A.Contains(B) is not the same as B.Contains(A).
        public LoopCrosser(S2Loop a, S2Loop b, LoopRelation relation, bool swapped)
        {
            a_ = a;
            b_ = b;
            relation_ = relation;
            swapped_ = swapped;
            ACrossingTarget = relation.ACrossingTarget();
            BCrossingTarget = relation.BCrossingTarget();
            b_query_ = new S2CrossingEdgeQuery(b._index);
            if (swapped)
            {
                (BCrossingTarget, ACrossingTarget) = (ACrossingTarget, BCrossingTarget);
            }
        }

        #endregion

        #region LoopCrosser

        // Given two iterators positioned such that ai.id().Contains(bi.id()),
        // return true if there is a crossing relationship anywhere within ai.id().
        // Specifically, this method returns true if there is an edge crossing, a
        // wedge crossing, or a point P that matches both "crossing targets".
        // Advances both iterators past ai.id().
        public bool HasCrossingRelation(RangeEnumerator ai, RangeEnumerator bi)
        {
            MyDebug.Assert(ai.Id.Contains(bi.Id));
            if (ai.NumEdges() == 0)
            {
                if (ai.ContainsCenter() == (ACrossingTarget != 0))
                {
                    // All points within ai.id() satisfy the crossing target for A, so it's
                    // worth iterating through the cells of B to see whether any cell
                    // centers also satisfy the crossing target for B.
                    do
                    {
                        if (bi.ContainsCenter() == (BCrossingTarget != 0)) return true;
                        bi.MoveNext();
                    } while (bi.Id <= ai.RangeMax);
                }
                else
                {
                    // The crossing target for A is not satisfied, so we skip over the cells
                    // of B using binary search.
                    bi.SeekBeyond(ai);
                }
            }
            else
            {
                // The current cell of A has at least one edge, so check for crossings.
                if (HasCrossing(ai, bi)) return true;
            }
            ai.MoveNext();
            return false;
        }

        // Given two index cells, return true if there are any edge crossings or
        // wedge crossings within those cells.
        public bool CellCrossesCell(S2ClippedShape a_clipped, S2ClippedShape b_clipped)
        {
            // Test all edges of "a_clipped" against all edges of "b_clipped".
            int a_num_edges = a_clipped.NumEdges;
            for (int i = 0; i < a_num_edges; ++i)
            {
                StartEdge(a_clipped.Edges[i]);
                if (EdgeCrossesCell(b_clipped)) return true;
            }
            return false;
        }

        // Given two iterators positioned such that ai.id().Contains(bi.id()),
        // return true if there is an edge crossing or wedge crosssing anywhere
        // within ai.id().  Advances "bi" (only) past ai.id().
        private bool HasCrossing(RangeEnumerator ai, RangeEnumerator bi)
        {
            MyDebug.Assert(ai.Id.Contains(bi.Id));
            // If ai.id() intersects many edges of B, then it is faster to use
            // S2CrossingEdgeQuery to narrow down the candidates.  But if it intersects
            // only a few edges, it is faster to check all the crossings directly.
            // We handle this by advancing "bi" and keeping track of how many edges we
            // would need to test.

            const int kEdgeQueryMinEdges = 20;  // Tuned using benchmarks.
            int total_edges = 0;
            b_cells_.Clear();
            do
            {
                if (bi.NumEdges() > 0)
                {
                    total_edges += bi.NumEdges();
                    if (total_edges >= kEdgeQueryMinEdges)
                    {
                        // There are too many edges to test them directly, so use
                        // S2CrossingEdgeQuery.
                        if (CellCrossesAnySubcell(ai.Clipped(), ai.Id)) return true;
                        bi.SeekBeyond(ai);
                        return false;
                    }
                    b_cells_.Add(bi.Cell);
                }
                bi.MoveNext();
            } while (bi.Id <= ai.RangeMax);

            // Test all the edge crossings directly.
            foreach (var b_cell in b_cells_)
            {
                if (CellCrossesCell(ai.Clipped(), b_cell.Clipped(0)))
                {
                    return true;
                }
            }
            return false;
        }

        // Given an index cell of A, return true if there are any edge or wedge
        // crossings with any index cell of B contained within "b_id".
        private bool CellCrossesAnySubcell(S2ClippedShape a_clipped, S2CellId b_id)
        {
            // Test all edges of "a_clipped" against all edges of B.  The relevant B
            // edges are guaranteed to be children of "b_id", which lets us find the
            // correct index cells more efficiently.
            var b_root = new S2PaddedCell(b_id, 0);
            int a_num_edges = a_clipped.NumEdges;
            for (int i = 0; i < a_num_edges; ++i)
            {
                int aj = a_clipped.Edges[i];
                // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
                // of B that might contain crossing edges.
                b_query_.GetCells(a_.Vertex(aj), a_.Vertex(aj + 1), b_root, b_cells_);
                if (!b_cells_.Any()) continue;
                StartEdge(aj);
                foreach (var b_cell in b_cells_)
                {
                    if (EdgeCrossesCell(b_cell.Clipped(0))) return true;
                }
            }
            return false;
        }

        // Prepare to check the given edge of loop A for crossings.
        private void StartEdge(int aj)
        {
            // Start testing the given edge of A for crossings.
            crosser_.Init(a_.Vertex(aj), a_.Vertex(aj + 1));
            aj_ = aj;
            bj_prev_ = -2;
        }

        // Check the current edge of loop A for crossings with all edges of the
        // given index cell of loop B.
        private bool EdgeCrossesCell(S2ClippedShape b_clipped)
        {
            // Test the current edge of A against all edges of "b_clipped".
            int b_num_edges = b_clipped.NumEdges;
            for (int j = 0; j < b_num_edges; ++j)
            {
                int bj = b_clipped.Edges[j];
                if (bj != bj_prev_ + 1) crosser_.RestartAt(b_.Vertex(bj));
                bj_prev_ = bj;
                int crossing = crosser_.CrossingSign(b_.Vertex(bj + 1));
                if (crossing < 0) continue;
                if (crossing > 0) return true;
                // We only need to check each shared vertex once, so we only
                // consider the case where a_vertex(aj_+1) == b_.vertex(bj+1).
                if (a_.Vertex(aj_ + 1) == b_.Vertex(bj + 1))
                {
                    if (swapped_ ?
                        relation_.WedgesCross(
                            b_.Vertex(bj), b_.Vertex(bj + 1), b_.Vertex(bj + 2),
                            a_.Vertex(aj_), a_.Vertex(aj_ + 2)) :
                        relation_.WedgesCross(
                            a_.Vertex(aj_), a_.Vertex(aj_ + 1), a_.Vertex(aj_ + 2),
                            b_.Vertex(bj), b_.Vertex(bj + 2)))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        #endregion
    }

    // Loop relation for Contains().
    internal class ContainsRelation : LoopRelation
    {
        public bool FoundSharedVertex { get; private set; }

        public ContainsRelation()
        {
            FoundSharedVertex = false;
        }

        // If A.Contains(P) == false && B.Contains(P) == true, it is equivalent to
        // having an edge crossing (i.e., Contains returns false).
        public override int ACrossingTarget() { return 0; }
        public override int BCrossingTarget() { return 1; }

        public override bool WedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2)
        {
            FoundSharedVertex = true;
            return !S2WedgeRelations.WedgeContains(a0, ab1, a2, b0, b2);
        }
    }

    // Loop relation for Intersects().
    internal class IntersectsRelation : LoopRelation
    {
        public IntersectsRelation() { FoundSharedVertex = false; }
        public bool FoundSharedVertex { get; private set; }

        // If A.Contains(P) == true && B.Contains(P) == true, it is equivalent to
        // having an edge crossing (i.e., Intersects returns true).
        public override int ACrossingTarget() { return 1; }
        public override int BCrossingTarget() { return 1; }

        public override bool WedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2)
        {
            FoundSharedVertex = true;
            return S2WedgeRelations.WedgeIntersects(a0, ab1, a2, b0, b2);
        }
    }

    // Loop relation for CompareBoundary().
    internal class CompareBoundaryRelation : LoopRelation
    {
        public bool FoundSharedVertex { get; private set; }
        public bool ContainsEdge { get; private set; }
        private readonly bool reverse_b_;      // True if loop B should be reversed.
        private bool excludes_edge_;        // True if any edge of B is excluded by A.

        public CompareBoundaryRelation(bool reverse_b)
        {
            reverse_b_ = reverse_b;
            FoundSharedVertex = false;
            ContainsEdge = false;
            excludes_edge_ = false;
        }

        // The CompareBoundary relation does not have a useful early-exit condition,
        // so we return -1 for both crossing targets.
        //
        // Aside: A possible early exit condition could be based on the following.
        //   If A contains a point of both B and ~B, then A intersects Boundary(B).
        //   If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
        //   So if the intersections of {A, ~A} with {B, ~B} are all non-empty,
        //   the return value is 0, i.e., Boundary(A) intersects Boundary(B).
        // Unfortunately it isn't worth detecting this situation because by the
        // time we have seen a point in all four intersection regions, we are also
        // guaranteed to have seen at least one pair of crossing edges.
        public override int ACrossingTarget() { return -1; }
        public override int BCrossingTarget() { return -1; }

        public override bool WedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2)
        {
            // Because we don't care about the interior of B, only its boundary, it is
            // sufficient to check whether A contains the semiwedge (ab1, b2).
            FoundSharedVertex = true;
            if (WedgeContainsSemiwedge(a0, ab1, a2, b2, reverse_b_))
            {
                ContainsEdge = true;
            }
            else
            {
                excludes_edge_ = true;
            }
            return ContainsEdge & excludes_edge_;
        }
    }
}
