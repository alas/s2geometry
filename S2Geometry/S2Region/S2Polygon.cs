// An S2Polygon is an S2Region object that represents a polygon.  A polygon is
// defined by zero or more loops; recall that the interior of a loop is
// defined to be its left-hand side (see S2Loop).  There are two different
// conventions for creating an S2Polygon:
//
//   - InitNested() expects the input loops to be nested hierarchically.  The
//     polygon interior then consists of the set of points contained by an odd
//     number of loops.  So for example, a circular region with a hole in it
//     would be defined as two CCW loops, with one loop containing the other.
//     The loops can be provided in any order.
//
//     When the orientation of the input loops is unknown, the nesting
//     requirement is typically met by calling S2Loop.Normalize() on each
//     loop (which inverts the loop if necessary so that it encloses at most
//     half the sphere).  But in fact any set of loops can be used as long as
//     (1) there is no pair of loops that cross, and (2) there is no pair of
//     loops whose union is the entire sphere.
//
//   - InitOriented() expects the input loops to be oriented such that the
//     polygon interior is on the left-hand side of every loop.  So for
//     example, a circular region with a hole in it would be defined using a
//     CCW outer loop and a CW inner loop.  The loop orientations must all be
//     consistent; for example, it is not valid to have one CCW loop nested
//     inside another CCW loop, because the region between the two loops is on
//     the left-hand side of one loop and the right-hand side of the other.
//
// Most clients will not call these methods directly; instead they should use
// S2Builder, which has better support for dealing with imperfect data.
//
// When the polygon is initialized, the given loops are automatically
// converted into a canonical form consisting of "shells" and "holes".  Shells
// and holes are both oriented CCW, and are nested hierarchically.  The loops
// are reordered to correspond to a preorder traversal of the nesting
// hierarchy; InitOriented may also invert some loops. The set of input S2Loop
// pointers is always preserved; the caller can use this to determine how the
// loops were reordered if desired.
//
// Polygons may represent any region of the sphere with a polygonal boundary,
// including the entire sphere (known as the "full" polygon).  The full
// polygon consists of a single full loop (see S2Loop), whereas the empty
// polygon has no loops at all.
//
// Polygons have the following restrictions:
//
//  - Loops may not cross, i.e. the boundary of a loop may not intersect
//    both the interior and exterior of any other loop.
//
//  - Loops may not share edges, i.e. if a loop contains an edge AB, then
//    no other loop may contain AB or BA.
//
//  - Loops may share vertices, however no vertex may appear twice in a
//    single loop (see S2Loop).
//
//  - No loop may be empty.  The full loop may appear only in the full polygon.

namespace S2Geometry;

using System.Runtime.InteropServices;

// A map from each loop to its immediate children with respect to nesting.
// This map is built during initialization of multi-loop polygons to
// determine which are shells and which are holes, and then discarded.
using LoopMap = Dictionary<S2Loop, List<S2Loop>>;
using S2BuilderUtil;

public sealed record class S2Polygon : IS2Region<S2Polygon>, IDecoder<S2Polygon>, IDisposable
{
    #region Fields, Constants

    // Returns the built-in S2ShapeIndex associated with every S2Polygon.  This
    // can be used in conjunction with the various S2ShapeIndex query classes
    // (S2ClosestEdgeQuery, S2BooleanOperation, etc) to do things beyond what is
    // possible with S2Polygon built-in convenience methods.
    //
    // For example, to measure the distance from one S2Polygon to another, you
    // can write:
    //   S2ClosestEdgeQuery query(out polygon1.index());
    //   S2ClosestEdgeQuery.ShapeIndexTarget target(out polygon2.index());
    //   S1ChordAngle distance = query.GetDistance(out target);
    //
    // The index contains a single S2Polygon.Shape object.
    public MutableS2ShapeIndex Index { get; } = [];

    // Total number of vertices in all loops.
    // TODO(ericv): Change to num_edges() for consistency with the S2Shape API.
    public int NumVertices { get; private set; }

    // Allows overriding the automatic validity checking controlled by the
    // --s2debug flag.
    //
    // Allows overriding the automatic validity checks controlled by
    // --s2debug (which is true by default in non-optimized builds).
    // When this flag is enabled, a fatal error is generated whenever
    // an invalid polygon is constructed.  The main reason to disable
    // this flag is if you intend to call IsValid() explicitly, like this:
    //
    //   S2Polygon polygon;
    //   polygon.set_s2debug_override(S2Debug::DISABLE);
    //   polygon.Init(...);
    //   if (!polygon.IsValid()) { ... }
    //
    // This setting is preserved across calls to Init() and Decode().
    public S2Debug S2DebugOverride { get; set; }

    private readonly List<S2Loop> loops_ = [];

    // True if InitOriented() was called and the given loops had inconsistent
    // orientations (i.e., it is not possible to construct a polygon such that
    // the interior is on the left-hand side of all loops).  We need to remember
    // this error so that it can be returned later by FindValidationError(),
    // since it is not possible to detect this error once the polygon has been
    // initialized.  This field is not preserved by Encode/Decode.
    private bool error_inconsistent_loop_orientations_;

    // In general we build the index the first time it is needed, but we make an
    // exception for Contains(S2Point) because this method has a simple brute
    // force implementation that is also relatively cheap.  For this one method
    // we keep track of the number of calls made and only build the index once
    // enough calls have been made that we think an index would be worthwhile.
    private Label unindexed_contains_calls_;

    // "bound_" is a conservative bound on all points contained by this polygon:
    // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
    private S2LatLngRect bound_;

    // Since "bound_" is not exact, it is possible that a polygon A contains
    // another polygon B whose bounds are slightly larger.  "subregion_bound_"
    // has been expanded sufficiently to account for this error, i.e.
    // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
    private S2LatLngRect subregion_bound_;

    // The maximum number of loops we'll allow when decoding a polygon.
    // The default value of 10 million is 200x bigger than the number of
    private const int s2polygon_decode_max_num_loops = 10000000; // The upper limit on the number of loops that are allowed by the S2Polygon.Decode method.

    // When adding a new encoding, be aware that old binaries will not
    // be able to decode it.
    private const byte kCurrentUncompressedEncodingVersionNumber = 1;
    private const byte kCurrentCompressedEncodingVersionNumber = 4;

    #endregion

    #region Constructors

    // The default constructor creates an empty polygon.  It can be made
    // non-empty by calling Init(), Decode(), etc.
    public S2Polygon()
    {
        S2DebugOverride = S2Debug.ALLOW;
        error_inconsistent_loop_orientations_ = false;
        unindexed_contains_calls_ = 0;
    }

    // Convenience constructor that calls InitNested() with the given loops.
    public S2Polygon(List<S2Loop> loops, bool initOriented = false, S2Debug s2debug_override = S2Debug.ALLOW)
    {
        S2DebugOverride = s2debug_override;
        loops_ = loops;
        if (initOriented) InitOriented();
        else InitNested();
    }

    // Convenience constructor that creates a polygon with a single loop
    // corresponding to the given cell.
    public S2Polygon(S2Cell cell)
    {
        S2DebugOverride = S2Debug.ALLOW;
        Init(new S2Loop(cell));
    }

    // Convenience constructor that calls Init(S2Loop*).  Note that this method
    // automatically converts the special empty loop (see S2Loop) into an empty
    // polygon, unlike the vector-of-loops constructor which does not allow
    // empty loops at all.
    public S2Polygon(S2Loop loop, S2Debug s2debug_override = S2Debug.ALLOW)
    {
        S2DebugOverride = s2debug_override;
        Init(loop);
    }

    private S2Polygon(List<S2Loop> loops, S2LatLngRect bound, S2LatLngRect subregion_bound)
    {
        loops_ = loops;
        bound_ = bound;
        subregion_bound_ = subregion_bound;
        InitIndex();
    }

    #endregion

    #region S2Polygon

    // Create a polygon from a set of hierarchically nested loops.  The polygon
    // interior consists of the points contained by an odd number of loops.
    // (Recall that a loop contains the set of points on its left-hand side.)
    //
    // This method figures out the loop nesting hierarchy and assigns every
    // loop a depth.  Shells have even depths, and holes have odd depths.  Note
    // that the loops are reordered so the hierarchy can be traversed more
    // easily (see GetParent, GetLastDescendant(), and S2Loop.depth).
    //
    // This method may be called more than once, in which case any existing
    // loops are deleted before being replaced by the input loops.
    private void InitNested()
    {
        if (NumLoops() == 1)
        {
            InitOneLoop();
            return;
        }
        var loop_map = new LoopMap();
        for (var i = 0; i < NumLoops(); ++i)
        {
            InsertLoop(Loop(i), S2Loop.NullLoop(), loop_map);
        }
        loops_.Clear();
        InitLoops(loop_map);

        // Compute num_vertices_, bound_, subregion_bound_.
        InitLoopProperties();
    }

    // Like InitNested(), but expects loops to be oriented such that the polygon
    // interior is on the left-hand side of all loops.  This implies that shells
    // and holes should have opposite orientations in the input to this method.
    // (During initialization, loops representing holes will automatically be
    // inverted.)
    public void InitOriented()
    {
        // Here is the algorithm:
        //
        // 1. Remember which of the given loops contain S2.Origin.
        //
        // 2. Invert loops as necessary to ensure that they are nestable (i.e., no
        //    loop contains the complement of any other loop).  This may result in a
        //    set of loops corresponding to the complement of the given polygon, but
        //    we will fix that problem later.
        //
        //    We make the loops nestable by first normalizing all the loops (i.e.,
        //    inverting any loops whose curvature is negative).  This handles
        //    all loops except those whose curvature is very close to zero
        //    (within the maximum error tolerance).  Any such loops are inverted if
        //    and only if they contain S2.Origin.  (In theory this step is only
        //    necessary if there are at least two such loops.)  The resulting set of
        //    loops is guaranteed to be nestable.
        //
        // 3. Build the polygon.  This yields either the desired polygon or its
        //    complement.
        //
        // 4. If there is at least one loop, we find a loop L that is adjacent to
        //    S2.Origin (where "adjacent" means that there exists a path
        //    connecting S2.Origin to some vertex of L such that the path does
        //    not cross any loop).  There may be a single such adjacent loop, or
        //    there may be several (in which case they should all have the same
        //    contains_origin() value).  We choose L to be the loop containing the
        //    origin whose depth is greatest, or loop(0) (a top-level shell) if no
        //    such loop exists.
        //
        // 5. If (L originally contained origin) != (polygon contains origin), we
        //    invert the polygon.  This is done by inverting a top-level shell whose
        //    curvature is minimal and then fixing the nesting hierarchy.  Note
        //    that because we normalized all the loops initially, this step is only
        //    necessary if the polygon requires at least one non-normalized loop to
        //    represent it.

        SortedSet<S2Loop> contained_origin = [];
        foreach (var loop in loops_)
        {
            if (loop.ContainsOrigin)
            {
                contained_origin.Add(loop);
            }
            var angle = loop.Curvature();
            if (Math.Abs(angle) > loop.CurvatureMaxError())
            {
                // Normalize the loop.
                if (angle < 0) loop.Invert();
            }
            else
            {
                // Ensure that the loop does not contain the origin.
                if (loop.ContainsOrigin) loop.Invert();
            }
        }
        InitNested();
        if (NumLoops() > 0)
        {
            var origin_loop = Loop(0);
            var polygon_contains_origin = false;
            for (var i = 0; i < NumLoops(); ++i)
            {
                if (Loop(i).ContainsOrigin)
                {
                    polygon_contains_origin ^= true;
                    origin_loop = Loop(i);
                }
            }
            if (contained_origin.Contains(origin_loop) != polygon_contains_origin)
            {
                Invert();
            }
        }
        // Verify that the original loops had consistent shell/hole orientations.
        // Each original loop L should have been inverted if and only if it now
        // represents a hole.
        foreach (var loop in loops_)
        {
            if (contained_origin.Contains(loop) != loop.ContainsOrigin != loop.IsHole())
            {
                // There is no point in saving the loop index, because the error is a
                // property of the entire set of loops.  In general there is no way to
                // determine which ones are incorrect.
                error_inconsistent_loop_orientations_ = true;
#if s2debug
                // The s2debug validity checking usually happens in InitIndex(),
                // but this error is detected too late for that.
                MyDebug.Assert(IsValid(), "Always fails.");
#endif
            }
        }
    }

    // Initialize a polygon from a single loop.  Note that this method
    // automatically converts the special empty loop (see S2Loop) into an empty
    // polygon, unlike the vector-of-loops InitNested() method which does not
    // allow empty loops at all.
    private void Init(S2Loop loop)
    {
        // We don't allow empty loops in the other Init() methods because deleting
        // them changes the number of loops, which is awkward to handle.
        if (loop.IsEmpty())
        {
            InitLoopProperties();
        }
        else
        {
            loops_.Add(loop);
            InitOneLoop();
        }
    }

    // Releases ownership of and returns the loops of this polygon, and resets
    // the polygon to be empty.
    public List<S2Loop> Release()
    {
        // Reset the polygon to be empty.
        var loops = loops_.ToList();
        ClearLoops();
        NumVertices = 0;
        bound_ = S2LatLngRect.Empty;
        subregion_bound_ = S2LatLngRect.Empty;
        return loops;
    }

    // Makes a deep copy of the given source polygon.  The destination polygon
    // will be cleared if necessary.
    public void Copy(S2Polygon src)
    {
        ClearLoops();
        for (var i = 0; i < src.NumLoops(); ++i)
        {
            loops_.Add((S2Loop)src.Loop(i).CustomClone());
        }
        // Don't copy error_inconsistent_loop_orientations_, since this is not a
        // property of the polygon but only of the way the polygon was constructed.
        NumVertices = src.NumVertices;
        unindexed_contains_calls_ = 0;
        bound_ = src.bound_;
        subregion_bound_ = src.subregion_bound_;
        InitIndex();  // TODO(ericv): Copy the index efficiently.
    }

    // Destroys the polygon and frees its loops.
    //public   ~S2Polygon() override;

    // Returns true if this is a valid polygon (including checking whether all
    // the loops are themselves valid).  Note that validity is checked
    // automatically during initialization when --s2debug is enabled (true by
    // default in debug binaries).
    public bool IsValid() => !FindValidationError(out _);

    // Returns true if this is *not* a valid polygon and sets "error"
    // appropriately.  Otherwise returns false and leaves "error" unchanged.
    //
    // Note that in error messages, loops that represent holes have their edges
    // numbered in reverse order, starting from the last vertex of the loop.
    //
    // REQUIRES: error != null
    public bool FindValidationError(out S2Error error)
    {
        for (var i = 0; i < NumLoops(); ++i)
        {
            // Check for loop errors that don't require building an S2ShapeIndex.
            if (Loop(i).FindValidationErrorNoIndex(out error))
            {
                error = new(error.Code, $"Loop {i}: {error.Text}");
                return true;
            }
            // Check that no loop is empty, and that the full loop only appears in the
            // full polygon.
            if (Loop(i).IsEmpty())
            {
                error = new(S2ErrorCode.POLYGON_EMPTY_LOOP, $"Loop {i}: empty loops are not allowed");
                return true;
            }
            if (Loop(i).IsFull() && NumLoops() > 1)
            {
                error = new(S2ErrorCode.POLYGON_EXCESS_FULL_LOOP, $"Loop {i}: full loop appears in non-full polygon");
                return true;
            }
        }

        // Check for loop self-intersections and loop pairs that cross
        // (including duplicate edges and vertices).
        if (S2ShapeUtil.EdgePairs.FindSelfIntersection(Index, out error)) return true;

        // Check whether InitOriented detected inconsistent loop orientations.
        if (error_inconsistent_loop_orientations_)
        {
            error = new(S2ErrorCode.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS, "Inconsistent loop orientations detected");
            return true;
        }

        // Finally, verify the loop nesting hierarchy.
        return FindLoopNestingError(out error);
    }

    // Return true if this is the empty polygon (consisting of no loops).
    public bool IsEmpty() => loops_.Count==0;

    // Return true if this is the full polygon (consisting of a single loop that
    // encompasses the entire sphere).
    public bool IsFull() => NumLoops() == 1 && Loop(0).IsFull();

    // Return the number of loops in this polygon.
    public int NumLoops() { return loops_.Count; }

    // Return the loop at the given index.  Note that during initialization, the
    // given loops are reordered according to a preorder traversal of the loop
    // nesting hierarchy.  This implies that every loop is immediately followed
    // by its descendants.  This hierarchy can be traversed using the methods
    // GetParent, GetLastDescendant(), and S2Loop.depth.
    public S2Loop Loop(int k) { return loops_[k]; }
    public List<S2Loop> Loops() { return [.. loops_]; }

    // Return the index of the parent of loop k, or -1 if it has no parent.
    public int GetParent(int k)
    {
        var depth = Loop(k).Depth;
        if (depth == 0) return -1;  // Optimization.
        while (--k >= 0 && Loop(k).Depth >= depth) continue;
        return k;
    }

    // Return the index of the last loop that is contained within loop k.
    // Returns num_loops() - 1 if k < 0.  Note that loops are indexed according
    // to a preorder traversal of the nesting hierarchy, so the immediate
    // children of loop k can be found by iterating over loops
    // (k+1)..GetLastDescendant(k) and selecting those whose depth is equal to
    // (loop(k).depth + 1).
    public int GetLastDescendant(int k)
    {
        if (k < 0) return NumLoops() - 1;
        var depth = Loop(k).Depth;
        while (++k < NumLoops() && Loop(k).Depth > depth) continue;
        return k - 1;
    }

    // Return the area of the polygon interior, i.e. the region on the left side
    // of an odd number of loops.  The return value is between 0 and 4*Pi.
    public double GetArea()
    {
        double area = 0;
        for (var i = 0; i < NumLoops(); ++i)
        {
            area += Loop(i).Sign() * Loop(i).Area();
        }
        return area;
    }

    // Return the true centroid of the polygon multiplied by the area of the
    // polygon (see s2centroids.h for details on centroids).  The result is not
    // unit length, so you may want to normalize it.  Also note that in general,
    // the centroid may not be contained by the polygon.
    //
    // We prescale by the polygon area for two reasons: (1) it is cheaper to
    // compute this way, and (2) it makes it easier to compute the centroid of
    // more complicated shapes (by splitting them into disjoint regions and
    // adding their centroids).
    public S2Point GetCentroid()
    {
        var centroid = S2Point.Empty;
        for (var i = 0; i < NumLoops(); ++i)
        {
            centroid += Loop(i).Sign() * Loop(i).Centroid();
        }
        return centroid;
    }

    // If all of the polygon's vertices happen to be the centers of S2Cells at
    // some level, then return that level, otherwise return -1.  See also
    // InitToSnapped() and s2builderutil.S2CellIdSnapFunction.
    // Returns -1 if the polygon has no vertices.
    public int GetSnapLevel()
    {
        var snap_level = -1;
        foreach (var child in loops_)
        {
            for (var j = 0; j < child.NumVertices; ++j)
            {
                var level = S2.XYZtoFaceSiTi(child.Vertex(j), out var face, out var si, out var ti);
                if (level < 0) return level;  // Vertex is not a cell center.
                if (level != snap_level)
                {
                    if (snap_level < 0)
                    {
                        snap_level = level;  // First vertex.
                    }
                    else
                    {
                        return -1;  // Vertices at more than one cell level.
                    }
                }
            }
        }
        return snap_level;
    }

    // Return the distance from the given point to the polygon interior.  If the
    // polygon is empty, return S1Angle.Infinity().  "x" should be unit length.
    public S1Angle GetDistance(S2Point x)
    {
        // Note that S2Polygon.Contains(S2Point) is slightly more efficient than
        // the generic version used by S2ClosestEdgeQuery.
        if (Contains(x)) return S1Angle.Zero;
        return GetDistanceToBoundary(x);
    }

    // Return the distance from the given point to the polygon boundary.  If the
    // polygon is empty or full, return S1Angle.Infinity() (since the polygon
    // has no boundary).  "x" should be unit length.
    public S1Angle GetDistanceToBoundary(S2Point x)
    {
        S2ClosestEdgeQuery.Options options = new()
        {
            IncludeInteriors = false,
        };
        S2ClosestEdgeQuery.PointTarget t = new(x);
        return new S2ClosestEdgeQuery(Index, options).GetDistance(t).ToAngle();
    }

    // Return the overlap fractions between two polygons, i.e. the ratios of the
    // area of intersection to the area of each polygon.
    public static (double, double) GetOverlapFractions(S2Polygon a, S2Polygon b)
    {
        var intersection = new S2Polygon();
        intersection.InitToIntersection(a, b);
        var intersection_area = intersection.GetArea();
        var a_area = a.GetArea();
        var b_area = b.GetArea();
        return (
            intersection_area >= a_area ? 1 : intersection_area / a_area,
            intersection_area >= b_area ? 1 : intersection_area / b_area);
    }

    // If the given point is contained by the polygon, return it.  Otherwise
    // return the closest point on the polygon boundary.  If the polygon is
    // empty, return the input argument.  Note that the result may or may not be
    // contained by the polygon.  "x" should be unit length.
    public S2Point Project(S2Point x)
    {
        if (Contains(x)) return x;
        return ProjectToBoundary(x);
    }

    // Return the closest point on the polygon boundary to the given point.  If
    // the polygon is empty or full, return the input argument (since the
    // polygon has no boundary).  "x" should be unit length.
    public S2Point ProjectToBoundary(S2Point x)
    {
        S2ClosestEdgeQuery.Options? options = new()
        {
            IncludeInteriors = false,
        };
        S2ClosestEdgeQuery q = new(Index, options);
        S2ClosestEdgeQuery.PointTarget target = new(x);
        var edge = q.FindClosestEdge(target);
        return q.Project(x, edge);
    }

    // Return true if this polygon contains the given other polygon, i.e.
    // if polygon A contains all points contained by polygon B.
    public bool Contains(S2Polygon b)
    {
        // It's worth checking bounding rectangles, since they are precomputed.
        // Note that the first bound has been expanded to account for possible
        // numerical errors in the second bound.
        if (!subregion_bound_.Contains(b.bound_))
        {
            // It is possible that A contains B even though Bound(A) does not contain
            // Bound(B).  This can only happen when polygon B has at least two outer
            // shells and the union of the two bounds spans all longitudes.  For
            // example, suppose that B consists of two shells with a longitude gap
            // between them, while A consists of one shell that surrounds both shells
            // of B but goes the other way around the sphere (so that it does not
            // intersect the longitude gap).
            //
            // For convenience we just check whether B has at least two loops rather
            // than two outer shells.
            if (b.NumLoops() == 1 || !bound_.Lng.Union(b.bound_.Lng).IsFull())
            {
                return false;
            }
        }

        // The following case is not handled by S2BooleanOperation because it only
        // determines whether the boundary of the result is empty (which does not
        // distinguish between the full and empty polygons).
        if (IsEmpty() && b.IsFull()) return false;

        return S2BooleanOperation.Contains(Index, b.Index);
    }

    // Returns true if this polgyon (A) approximately contains the given other
    // polygon (B). This is true if it is possible to move the vertices of B
    // no further than "tolerance" such that A contains the modified B.
    //
    // For example, the empty polygon will contain any polygon whose maximum
    // width is no more than "tolerance".
    public bool ApproxContains(S2Polygon b, S1Angle tolerance)
    {
        var difference = new S2Polygon();
        difference.InitToDifference(b, this, new IdentitySnapFunction(tolerance));
        return difference.IsEmpty();
    }

    // Return true if this polygon intersects the given other polygon, i.e.
    // if there is a point that is contained by both polygons.
    public bool Intersects(S2Polygon b)
    {
        // It's worth checking bounding rectangles, since they are precomputed.
        if (!bound_.Intersects(b.bound_)) return false;

        // The following case is not handled by S2BooleanOperation because it only
        // determines whether the boundary of the result is empty (which does not
        // distinguish between the full and empty polygons).
        if (IsFull() && b.IsFull()) return true;

        return S2BooleanOperation.Intersects(Index, b.Index);
    }

    // Returns true if this polgyon (A) and the given polygon (B) are
    // approximately disjoint.  This is true if it is possible to ensure that A
    // and B do not intersect by moving their vertices no further than
    // "tolerance".  This implies that in borderline cases where A and B overlap
    // slightly, this method returns true (A and B are approximately disjoint).
    //
    // For example, any polygon is approximately disjoint from a polygon whose
    // maximum width is no more than "tolerance".
    public bool ApproxDisjoint(S2Polygon b, S1Angle tolerance)
    {
        var intersection = new S2Polygon();
        intersection.InitToIntersection(b, this, new IdentitySnapFunction(tolerance));
        return intersection.IsEmpty();
    }

    // Initialize this polygon to the intersection, union, difference (A - B),
    // or symmetric difference (XOR) of the given two polygons.
    //
    // "snap_function" allows you to specify a minimum spacing between output
    // vertices, and/or that the vertices should be snapped to a discrete set of
    // points (e.g. S2CellId centers or E7 lat/lng coordinates).  Any snap
    // function can be used, including the IdentitySnapFunction with a
    // snap_radius of zero (which preserves the input vertices exactly).
    //
    // The boundary of the output polygon before snapping is guaranteed to be
    // accurate to within S2EdgeCrossings.kIntersectionError of the exact result.
    // Snapping can move the boundary by an additional distance that depends on
    // the snap function.  Finally, any degenerate portions of the output
    // polygon are automatically removed (i.e., regions that do not contain any
    // points) since S2Polygon does not allow such regions.
    //
    // See S2Builder and s2builderutil for more details on snap functions.  For
    // example, you can snap to E7 coordinates by setting "snap_function" to
    // s2builderutil.IntLatLngSnapFunction(7).
    //
    // The default snap function is the IdentitySnapFunction with a snap radius
    // of S2EdgeCrossings.kIntersectionMergeRadius (equal to about 1.8e-15 radians
    // or 11 nanometers on the Earth's surface).  This means that vertices may
    // be positioned arbitrarily, but vertices that are extremely close together
    // can be merged together.  The reason for a non-zero default snap radius is
    // that it helps to eliminate narrow cracks and slivers when T-vertices are
    // present.  For example, adjacent S2Cells at different levels do not share
    // exactly the same boundary, so there can be a narrow crack between them.
    // If a polygon is intersected with those cells and the pieces are unioned
    // together, the result would have a narrow crack unless the snap radius is
    // set to a non-zero value.
    //
    // Note that if you want to encode the vertices in a lower-precision
    // representation (such as S2CellIds or E7), it is much better to use a
    // suitable SnapFunction rather than rounding the vertices yourself, because
    // this will create self-intersections unless you ensure that the vertices
    // and edges are sufficiently well-separated first.  In particular you need
    // to use a snap function whose min_edge_vertex_separation() is at least
    // twice the maximum distance that a vertex can move when rounded.
    //
    // The versions of these functions with an S2Error argument return true on
    // success and set "error" appropriately otherwise.  However note that these
    // functions should never return an error provided that both input polygons
    // are valid (i.e., IsValid returns true).
    public void InitToIntersection(S2Polygon a, S2Polygon b)
    {
        InitToIntersection(a, b, new IdentitySnapFunction(S2.kIntersectionMergeRadiusS1Angle));
    }
    public void InitToIntersection(S2Polygon a, S2Polygon b, SnapFunction snap_function)
    {
        if (!a.bound_.Intersects(b.bound_))
        {
            //InitNested({ });
            return;
        }
        InitToOperation(S2BooleanOperation.OpType.INTERSECTION,
                        snap_function, a, b);
    }
    public bool InitToIntersection(S2Polygon a, S2Polygon b, SnapFunction snap_function, out S2Error error)
    {
        error = S2Error.OK;
        if (!a.bound_.Intersects(b.bound_))
        {
            //InitNested({ });
            return true;  // Success.
        }
        return InitToOperation(S2BooleanOperation.OpType.INTERSECTION, snap_function, a, b, out error);
    }

    public void InitToUnion(S2Polygon a, S2Polygon b)
    {
        InitToUnion(a, b, new IdentitySnapFunction(S2.kIntersectionMergeRadiusS1Angle));
    }
    public void InitToUnion(S2Polygon a, S2Polygon b, SnapFunction snap_function)
    {
        InitToOperation(S2BooleanOperation.OpType.UNION, snap_function, a, b);
    }
    public bool InitToUnion(S2Polygon a, S2Polygon b, SnapFunction snap_function, out S2Error error)
    {
        return InitToOperation(S2BooleanOperation.OpType.UNION, snap_function, a, b, out error);
    }

    public void InitToDifference(S2Polygon a, S2Polygon b)
    {
        InitToDifference(a, b, new IdentitySnapFunction(S2.kIntersectionMergeRadiusS1Angle));
    }
    public void InitToDifference(S2Polygon a, S2Polygon b, SnapFunction snap_function)
    {
        InitToOperation(S2BooleanOperation.OpType.DIFFERENCE, snap_function, a, b);
    }
    public bool InitToDifference(S2Polygon a, S2Polygon b, SnapFunction snap_function, out S2Error error)
    {
        return InitToOperation(S2BooleanOperation.OpType.DIFFERENCE, snap_function, a, b, out error);
    }

    public void InitToSymmetricDifference(S2Polygon a, S2Polygon b)
    {
        InitToSymmetricDifference(a, b, new IdentitySnapFunction(S2.kIntersectionMergeRadiusS1Angle));
    }
    public void InitToSymmetricDifference(S2Polygon a, S2Polygon b, SnapFunction snap_function)
    {
        InitToOperation(S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, snap_function, a, b);
    }
    public bool InitToSymmetricDifference(S2Polygon a, S2Polygon b, SnapFunction snap_function, out S2Error error)
    {
        return InitToOperation(S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, snap_function, a, b, out error);
    }

    // Convenience function that snaps the vertices to S2CellId centers at the
    // given level (default level 30, which has S2CellId centers spaced about 1
    // centimeter apart).  Polygons can be efficiently encoded by Encode() after
    // they have been snapped.
    public void InitToSnapped(S2Polygon polygon, SnapFunction snapFunction)
    {
        S2Builder builder = new(new(snapFunction));
        InitFromBuilder(polygon, builder);
    }
    public void InitToSnapped(S2Polygon a, int snapLevel = S2.kMaxCellLevel)
    {
        InitToSnapped(a, new S2CellIdSnapFunction(snapLevel));
    }

    // Snaps the input polygon according to the given "snap_function" and
    // reduces the number of vertices if possible, while ensuring that no vertex
    // moves further than snap_function.snap_radius().
    //
    // Simplification works by replacing nearly straight chains of short edges
    // with longer edges, in a way that preserves the topology of the input
    // polygon up to the creation of degeneracies.  This means that loops or
    // portions of loops may become degenerate, in which case they are removed.
    // For example, if there is a very small island in the original polygon, it
    // may disappear completely.  (Even if there are dense islands, they could
    // all be removed rather than being replaced by a larger simplified island
    // if more area is covered by water than land.)
    public void InitToSimplified(S2Polygon a, SnapFunction snap_function)
    {
        var options = new S2Builder.Options(snap_function)
        {
            SimplifyEdgeChains = true
        };
        var builder = new S2Builder(options);
        InitFromBuilder(a, builder);
    }

    public void InitToSimplifiedInCell(S2Polygon a, S2Cell cell, S1Angle snap_radius)
    {
        InitToSimplifiedInCell(a, cell, snap_radius, S1Angle.FromRadians(S2.DoubleError));
    }

    // Like InitToSimplified, except that any vertices or edges on the boundary
    // of the given S2Cell are preserved if possible.  This method requires that
    // the polygon has already been clipped so that it does not extend outside
    // the cell by more than "boundary_tolerance".  In other words, it operates
    // on polygons that have already been intersected with a cell.
    //
    // Typically this method is used in geometry-processing pipelines that
    // intersect polygons with a collection of S2Cells and then process those
    // cells in parallel, where each cell generates some geometry that needs to
    // be simplified.  In contrast, if you just need to simplify the *input*
    // geometry then it is easier and faster to do the simplification before
    // computing the intersection with any S2Cells.
    //
    // "boundary_tolerance" specifies how close a vertex must be to the cell
    // boundary to be kept.  The default tolerance is large enough to handle any
    // reasonable way of interpolating points along the cell boundary, such as
    // S2EdgeCrossings.GetIntersection(), S2EdgeDistances.Interpolate(), or direct (u,v)
    // interpolation using S2Coords.FaceUVtoXYZ().  However, if the vertices have
    // been snapped to a lower-precision representation (e.g., S2CellId centers
    // or E7 coordinates) then you will need to set this tolerance explicitly.
    // For example, if the vertices were snapped to E7 coordinates then
    // "boundary_tolerance" should be set to
    //
    //   s2builderutil.IntLatLngSnapFunction.MinSnapRadiusForExponent(7)
    //
    // Degenerate portions of loops are always removed, so if a vertex on the
    // cell boundary belongs only to degenerate regions then it will not be
    // kept.  For example, if the input polygon is a narrow strip of width less
    // than "snap_radius" along one side of the cell, then the entire loop may
    // become degenerate and be removed.
    //
    // REQUIRES: all vertices of "a" are within "boundary_tolerance" of "cell".
    public void InitToSimplifiedInCell(S2Polygon a, S2Cell cell, S1Angle snap_radius, S1Angle boundary_tolerance)
    {
        // The polygon to be simplified consists of "boundary edges" that follow the
        // cell boundary and "interior edges" that do not.  We want to simplify the
        // interior edges while leaving the boundary edges unchanged.  It's not
        // sufficient to call S2Builder.ForceVertex() on all boundary vertices.
        // For example, suppose the polygon includes a triangle ABC where all three
        // vertices are on the cell boundary and B is a cell corner.  Then if
        // interior edge AC snaps to vertex B, this loop would become degenerate and
        // be removed.  Similarly, we don't want boundary edges to snap to interior
        // vertices, since this also would cause portions of the polygon along the
        // boundary to be removed.
        //
        // Instead we use a two-pass algorithm.  In the first pass, we simplify
        // *only* the interior edges, using ForceVertex() to ensure that any edge
        // endpoints on the cell boundary do not move.  In the second pass, we add
        // the boundary edges (which are guaranteed to still form loops with the
        // interior edges) and build the output polygon.
        //
        // Note that in theory, simplifying the interior edges could create an
        // intersection with one of the boundary edges, since if two interior edges
        // intersect very near the boundary then the intersection point could be
        // slightly outside the cell (by at most S2EdgeCrossings.kIntersectionError).
        // This is the *only* way that a self-intersection can be created, and it is
        // expected to be extremely rare.  Nevertheless we use a small snap radius
        // in the second pass in order to eliminate any such self-intersections.
        //
        // We also want to preserve the cyclic vertex order of loops, so that the
        // original polygon can be reconstructed when no simplification is possible
        // (i.e., idempotency).  In order to do this, we group the input edges into
        // a sequence of polylines.  Each polyline contains only one type of edge
        // (interior or boundary).  We use S2Builder to simplify the interior
        // polylines, while the boundary polylines are passed through unchanged.
        // Each interior polyline is in its own S2Builder layer in order to keep the
        // edges in sequence.  This lets us ensure that in the second pass, the
        // edges are added in their original order so that S2PolygonLayer can
        // reconstruct the original loops.

        // We want an upper bound on how much "u" or "v" can change when a point on
        // the boundary of the S2Cell is moved away by up to "boundary_tolerance".
        // Inverting this, instead we could compute a lower bound on how far a point
        // can move away from an S2Cell edge when "u" or "v" is changed by a given
        // amount.  The latter quantity is simply (S2Metrics.kMinWidth.deriv() / 2)
        // under the S2_LINEAR_PROJECTION model, where we divide by 2 because we
        // want the bound in terms of (u = 2 * s - 1) rather than "s" itself.
        // Consulting s2metrics.cc, this value is Math.Sqrt(2/3)/2 = Math.Sqrt(1/6).
        // Going back to the original problem, this gives:
        var boundary_tolerance_uv = Math.Sqrt(6) * boundary_tolerance.Radians;

        // The first pass yields a collection of simplified polylines that preserve
        // the original cyclic vertex order.
        var polylines = SimplifyEdgesInCell(a, cell, boundary_tolerance_uv, snap_radius);

        // The second pass eliminates any intersections between interior edges and
        // boundary edges, and then assembles the edges into a polygon.
        var options = new S2Builder.Options(new IdentitySnapFunction(S2.kIntersectionErrorS1Angle))
        {
            Idempotent = false  // Force snapping up to the given radius
        };
        var builder = new S2Builder(options);
        builder.StartLayer(new S2PolygonLayer(this));
        foreach (var polyline in polylines)
        {
            builder.AddPolyline(polyline);
        }
        if (!builder.Build(out var error))
        {
            throw new ApplicationException("Could not build polygon: " + error);
        }
        // If there are no loops, check whether the result should be the full
        // polygon rather than the empty one.  (See InitToIntersection.)
        if (NumLoops() == 0)
        {
            if (a.bound_.Area() > S2.M_2_PI && a.GetArea() > S2.M_2_PI) Invert();
        }
    }

    // Initialize this polygon to the complement of the given polygon.
    public void InitToComplement(S2Polygon a)
    {
        Copy(a);
        Invert();
    }

    // Initialize this polygon to the outline of the given cell union.
    // In principle this polygon should exactly contain the cell union and
    // this polygon's inverse should not intersect the cell union, but rounding
    // issues may cause this not to be the case.
    public void InitToCellUnionBorder(S2CellUnion cells)
    {
        // We use S2Builder to compute the union.  Due to rounding errors, we can't
        // compute an exact union - when a small cell is adjacent to a larger cell,
        // the shared edges can fail to line up exactly.  Two cell edges cannot come
        // closer then kMinWidth, so if we have S2Builder snap edges within half
        // that distance, then we should always merge shared edges without merging
        // different edges.
        var snap_radius = 0.5 * S2.kMinWidth.GetValue(S2.kMaxCellLevel);
        var builder = new S2Builder(new S2Builder.Options(new IdentitySnapFunction(S1Angle.FromRadians(snap_radius))));
        builder.StartLayer(new S2PolygonLayer(this));
        foreach (var id in cells)
        {
            builder.AddLoop(new S2Loop(new S2Cell(id)));
        }
        if (!builder.Build(out var error))
        {
            throw new ApplicationException("InitToCellUnionBorder failed: " + error);
        }
        // If there are no loops, check whether the result should be the full
        // polygon rather than the empty one.  There are only two ways that this can
        // happen: either the cell union is empty, or it consists of all six faces.
        if (NumLoops() == 0)
        {
            if (cells.IsEmpty()) return;
            MyDebug.Assert(6UL << (2 * S2.kMaxCellLevel) == cells.LeafCellsCovered());
            Invert();
        }
    }

    // Given that loops_ contains a single loop, initialize all other fields.
    //
    // This is an internal method that expects that loops_ has already been
    // initialized with a single non-empty loop.
    private void InitOneLoop()
    {
        MyDebug.Assert(1 == NumLoops());
        var loop = loops_[0];
        loop.Depth = 0;
        error_inconsistent_loop_orientations_ = false;
        NumVertices = loop.NumVertices;
        bound_ = loop.GetRectBound();
        subregion_bound_ = S2LatLngRectBounder.ExpandForSubregions(bound_);
        InitIndex();
    }

    // Compute num_vertices_, bound_, subregion_bound_.
    private void InitLoopProperties()
    {
        NumVertices = 0;
        bound_ = S2LatLngRect.Empty;
        for (var i = 0; i < NumLoops(); ++i)
        {
            if (Loop(i).Depth == 0)
            {
                bound_ = bound_.Union(Loop(i).GetRectBound());
            }
            NumVertices += Loop(i).NumVertices;
        }
        subregion_bound_ = S2LatLngRectBounder.ExpandForSubregions(bound_);
        InitIndex();
    }

    // Invert the polygon (replace it by its complement).
    public void Invert()
    {
        // Inverting any one loop will invert the polygon.  The best loop to invert
        // is the one whose area is largest, since this yields the smallest area
        // after inversion.  The loop with the largest area is always at depth 0.
        // The descendents of this loop all have their depth reduced by 1, while the
        // former siblings of this loop all have their depth increased by 1.

        // The empty and full polygons are handled specially.
        if (IsEmpty())
        {
            loops_.Add(S2Loop.KFull);
        }
        else if (IsFull())
        {
            ClearLoops();
        }
        else
        {
            // Find the loop whose area is largest (i.e., whose curvature is
            // smallest), minimizing calls to GetCurvature().  In particular, for
            // polygons with a single shell at level 0 there is not need to call
            // GetCurvature() at all.  (This method is relatively expensive.)
            var best = 0;
            var kNone = 10.0;  // Flag that means "not computed yet"
            var best_angle = kNone;
            for (var i = 1; i < NumLoops(); ++i)
            {
                if (Loop(i).Depth == 0)
                {
                    // We defer computing the curvature of loop 0 until we discover
                    // that the polygon has another top-level shell.
                    if (best_angle == kNone) best_angle = Loop(best).Curvature();
                    var angle = Loop(i).Curvature();
                    // We break ties deterministically in order to avoid having the output
                    // depend on the input order of the loops.
                    if (angle < best_angle ||
                        (angle == best_angle && CompareLoops(Loop(i), Loop(best)) < 0))
                    {
                        best = i;
                        best_angle = angle;
                    }
                }
            }
            // Build the new loops vector, starting with the inverted loop.
            Loop(best).Invert();
            var new_loops = new List<S2Loop>(NumLoops());
            // Add the former siblings of this loop as descendants.
            var last_best = GetLastDescendant(best);
            new_loops.Add(loops_[best]);
            for (var i = 0; i < NumLoops(); ++i)
            {
                if (i < best || i > last_best)
                {
                    Loop(i).Depth = Loop(i).Depth + 1;
                    new_loops.Add(loops_[i]);
                }
            }
            // Add the former children of this loop as siblings.
            for (var i = 0; i < NumLoops(); ++i)
            {
                if (i > best && i <= last_best)
                {
                    Loop(i).Depth = Loop(i).Depth - 1;
                    new_loops.Add(loops_[i]);
                }
            }
            var tmp = loops_.ToList();
            loops_.Clear();
            loops_.AddRange(new_loops);
            new_loops.Clear();
            new_loops.AddRange(tmp);
            MyDebug.Assert(new_loops.Count == NumLoops());
        }
        ClearIndex();
        InitLoopProperties();
    }

    // Return true if this polygon contains the given polyline.  This method
    // returns an exact result, according to the following model:
    //
    //  - All edges are geodesics (of course).
    //
    //  - Vertices are ignored for the purposes of defining containment.
    //    (This is because polygons often do not contain their vertices, in
    //    order to that when a set of polygons tiles the sphere then every point
    //    is contained by exactly one polygon.)
    //
    //  - Points that lie exactly on geodesic edges are resolved using symbolic
    //    perturbations (i.e., they are considered to be infinitesmally offset
    //    from the edge).
    //
    //  - If the polygon and polyline share an edge, it is handled as follows.
    //    First, the polygon edges are oriented so that the interior is always
    //    on the left.  Then the shared polyline edge is contained if and only
    //    if it is in the same direction as the corresponding polygon edge.
    //    (This model ensures that when a polyline is intersected with a polygon
    //    and its complement, the edge only appears in one of the two results.)
    //
    // TODO(ericv): Update the implementation to correspond to the model above.
    public bool Contains(S2Polyline b)
    {
        return ApproxContains(b, S2.kIntersectionMergeRadiusS1Angle);
    }

    // Returns true if this polgyon approximately contains the given polyline
    // This is true if it is possible to move the polyline vertices no further
    // than "tolerance" such that the polyline is now contained.
    public bool ApproxContains(S2Polyline b, S1Angle tolerance)
    {
        var difference = ApproxSubtractFromPolyline(b, tolerance);
        return difference.Count==0;
    }

    // Return true if this polygon intersects the given polyline.  This method
    // returns an exact result; see Contains(S2Polyline) for details.
    public bool Intersects(S2Polyline b)
    {
        return !ApproxDisjoint(b, S2.kIntersectionMergeRadiusS1Angle);
    }

    // Returns true if this polgyon is approximately disjoint from the given
    // polyline.  This is true if it is possible to avoid intersection by moving
    // their vertices no further than "tolerance".
    //
    // This implies that in borderline cases where there is a small overlap,
    // this method returns true (i.e., they are approximately disjoint).
    public bool ApproxDisjoint(S2Polyline b, S1Angle tolerance)
    {
        var intersection = ApproxIntersectWithPolyline(b, tolerance);
        return intersection.Count==0;
    }

    // Intersect this polygon with the polyline "in" and return the resulting
    // zero or more polylines.  The polylines are returned in the order they
    // would be encountered by traversing "in" from beginning to end.
    // Note that the output may include polylines with only one vertex,
    // but there will not be any zero-vertex polylines.
    //
    // This is equivalent to calling ApproxIntersectWithPolyline() with the
    // "snap_radius" set to S2EdgeCrossings.kIntersectionMergeRadius.
    public List<S2Polyline> IntersectWithPolyline(S2Polyline input)
    {
        return ApproxIntersectWithPolyline(input, S2.kIntersectionMergeRadiusS1Angle);
    }

    // Similar to IntersectWithPolyline(), except that vertices will be
    // dropped as necessary to ensure that all adjacent vertices in the
    // sequence obtained by concatenating the output polylines will be
    // farther than "snap_radius" apart.  Note that this can change
    // the number of output polylines and/or yield single-vertex polylines.
    public List<S2Polyline> ApproxIntersectWithPolyline(S2Polyline input, S1Angle snap_radius)
    {
        return IntersectWithPolyline(input, new IdentitySnapFunction(snap_radius));
    }

    // TODO(ericv): Update documentation.
    public List<S2Polyline> IntersectWithPolyline(S2Polyline input, SnapFunction snap_function)
    {
        return OperationWithPolyline(S2BooleanOperation.OpType.INTERSECTION, snap_function, input);
    }

    // Same as IntersectWithPolyline, but subtracts this polygon from
    // the given polyline.
    public List<S2Polyline> SubtractFromPolyline(S2Polyline input)
    {
        return ApproxSubtractFromPolyline(input, S2.kIntersectionMergeRadiusS1Angle);
    }

    // Same as ApproxIntersectWithPolyline, but subtracts this polygon
    // from the given polyline.
    public List<S2Polyline> ApproxSubtractFromPolyline(S2Polyline input, S1Angle snap_radius)
    {
        return SubtractFromPolyline(input, new IdentitySnapFunction(snap_radius));
    }

    public List<S2Polyline> SubtractFromPolyline(S2Polyline input, SnapFunction snap_function)
    {
        return OperationWithPolyline(S2BooleanOperation.OpType.DIFFERENCE, snap_function, input);
    }

    // Return a polygon which is the union of the given polygons.
    public static S2Polygon DestructiveUnion(IEnumerable<S2Polygon> polygons)
    {
        return DestructiveApproxUnion(polygons, S2.kIntersectionMergeRadiusS1Angle);
    }
    public static S2Polygon DestructiveApproxUnion(IEnumerable<S2Polygon> polygons, S1Angle snap_radius)
    {
        return DestructiveUnion(polygons, new IdentitySnapFunction(snap_radius));
    }
    public static S2Polygon DestructiveUnion(IEnumerable<S2Polygon> polygons, SnapFunction snap_function)
    {
        // Effectively create a priority queue of polygons in order of number of
        // vertices.  Repeatedly union the two smallest polygons and add the result
        // to the queue until we have a single polygon to return.
        // using QueueType = multimap<int, S2Polygon>;
        var queue = new SortedSet<(int, S2Polygon)>();  // Map from # of vertices to polygon.
        foreach (var polygon in polygons)
            queue.Add((polygon.NumVertices, polygon));

        while (queue.Count > 1)
        {
            // Pop two simplest polygons from queue.
            var smallest_it = queue.First();
            var a_size = smallest_it.Item1;
            var a_polygon = smallest_it.Item2;
            queue.Remove(smallest_it);
            smallest_it = queue.First();
            var b_size = smallest_it.Item1;
            var b_polygon = smallest_it.Item2;
            queue.Remove(smallest_it);

            // Union and add result back to queue.
            var union_polygon = new S2Polygon();
            union_polygon.InitToUnion(a_polygon, b_polygon, snap_function);
            queue.Add((a_size + b_size, union_polygon));
            // We assume that the number of vertices in the union polygon is the
            // sum of the number of vertices in the original polygons, which is not
            // always true, but will almost always be a decent approximation, and
            // faster than recomputing.
        }

        if (queue.Count==0)
            return new S2Polygon();
        else
            return queue.First().Item2;
    }

    // Return true if every loop of this polygon shares at most one vertex with
    // its parent loop.  Every polygon has a unique normalized form.  A polygon
    // can be normalized by passing it through S2Builder (with no snapping) in
    // order to reconstruct the polygon from its edges.
    //
    // Generally there is no reason to convert polygons to normalized form.  It
    // is mainly useful for testing in order to compare whether two polygons
    // have exactly the same interior, even when they have a different loop
    // structure.  For example, a diamond nested within a square (touching at
    // four points) could be represented as a square with a diamond-shaped hole,
    // or as four triangles.  Methods such as BoundaryApproxEquals() will report
    // these polygons as being different (because they have different
    // boundaries) even though they contain the same points.  However if they
    // are both converted to normalized form (the "four triangles" version) then
    // they can be compared more easily.
    //
    // Also see ApproxEquals(), which can determine whether two polygons contain
    // approximately the same set of points without any need for normalization.
    public bool IsNormalized()
    {
        // TODO(ericv): The condition tested here is insufficient.  The correct
        // condition is that each *connected component* of child loops can share at
        // most one vertex with their parent loop.  Example: suppose loop A has
        // children B, C, D, and the following pairs are connected: AB, BC, CD, DA.
        // Then the polygon is not normalized.
        var vertices = new SortedSet<S2Point>();
        S2Loop? last_parent = null;
        for (var i = 0; i < NumLoops(); ++i)
        {
            var child = Loop(i);
            if (child.Depth == 0) continue;
            var parent = Loop(GetParent(i));
            if (parent != last_parent)
            {
                vertices.Clear();
                for (var j = 0; j < parent.NumVertices; ++j)
                {
                    vertices.Add(parent.Vertex(j));
                }
                last_parent = parent;
            }
            var count = 0;
            for (var j = 0; j < child.NumVertices; ++j)
            {
                if (vertices.Contains(child.Vertex(j))) ++count;
            }
            if (count > 1) return false;
        }
        return true;
    }

    // Return true if two polygons are approximately equal to within the given
    // tolerance.  This is true if it is possible to move the vertices of the
    // two polygons so that they contain the same set of points.
    //
    // Note that according to this model, small regions less than "tolerance" in
    // width do not need to be considered, since these regions can be collapsed
    // into degenerate loops (which contain no points) by moving their vertices.
    //
    // This model is not as strict as using the Hausdorff distance would be, and
    // it is also not as strict as BoundaryNear (defined below).  However, it is
    // a good choice for comparing polygons that have been snapped, simplified,
    // unioned, etc, since these operations use a model similar to this one
    // (i.e., degenerate loops or portions of loops are automatically removed).
    public bool ApproxEquals(S2Polygon b, S1Angle tolerance)
    {
        // TODO(ericv): This can be implemented more cheaply with S2Builder, by
        // simply adding all the edges from one polygon, adding the reversed edge
        // from the other polygon, and turning on the options to split edges and
        // discard sibling pairs.  Then the polygons are approximately equal if the
        // output graph has no edges.
        var symmetric_difference = new S2Polygon();
        symmetric_difference.InitToSymmetricDifference(
            b, this, new IdentitySnapFunction(tolerance));
        return symmetric_difference.IsEmpty();
    }

    // Returns true if two polygons have the same boundary.  More precisely,
    // this method requires that both polygons have loops with the same cyclic
    // vertex order and the same nesting hierarchy.  (This implies that vertices
    // may be cyclically rotated between corresponding loops, and the loop
    // ordering may be different between the two polygons as long as the nesting
    // hierarchy is the same.)
    public bool BoundaryEquals(S2Polygon b)
    {
        if (NumLoops() != b.NumLoops()) return false;

        for (var i = 0; i < NumLoops(); ++i)
        {
            var a_loop = Loop(i);
            var success = false;
            for (var j = 0; j < NumLoops(); ++j)
            {
                var b_loop = b.Loop(j);
                if ((b_loop.Depth == a_loop.Depth) &&
                    b_loop.BoundaryEquals(a_loop))
                {
                    success = true;
                    break;
                }
            }
            if (!success) return false;
        }
        return true;
    }

    // Return true if two polygons have the same boundary except for vertex
    // perturbations.  Both polygons must have loops with the same cyclic vertex
    // order and the same nesting hierarchy, but the vertex locations are
    // allowed to differ by up to "max_error".
    public bool BoundaryApproxEquals(S2Polygon b)
    {
        return BoundaryApproxEquals(b, S1Angle.FromRadians(S2.DoubleError));
    }
    public bool BoundaryApproxEquals(S2Polygon b, S1Angle max_error)
    {
        if (NumLoops() != b.NumLoops()) return false;

        // For now, we assume that there is at most one candidate match for each
        // loop.  (So far this method is just used for testing.)

        for (var i = 0; i < NumLoops(); ++i)
        {
            var a_loop = Loop(i);
            var success = false;
            for (var j = 0; j < NumLoops(); ++j)
            {
                var b_loop = b.Loop(j);
                if (b_loop.Depth == a_loop.Depth &&
                    b_loop.BoundaryApproxEquals(a_loop, max_error))
                {
                    success = true;
                    break;
                }
            }
            if (!success) return false;
        }
        return true;
    }

    // Return true if two polygons have boundaries that are within "max_error"
    // of each other along their entire lengths.  More precisely, there must be
    // a bijection between the two sets of loops such that for each pair of
    // loops, "a_loop.BoundaryNear(b_loop)" is true.
    public bool BoundaryNear(S2Polygon b)
    {
        return BoundaryNear(b, S1Angle.FromRadians(S2.DoubleError));
    }
    public bool BoundaryNear(S2Polygon b, S1Angle max_error)
    {
        if (NumLoops() != b.NumLoops()) return false;

        // For now, we assume that there is at most one candidate match for each
        // loop.  (So far this method is just used for testing.)

        for (var i = 0; i < NumLoops(); ++i)
        {
            var a_loop = Loop(i);
            var success = false;
            for (var j = 0; j < NumLoops(); ++j)
            {
                var b_loop = b.Loop(j);
                if (b_loop.Depth == a_loop.Depth &&
                    b_loop.BoundaryNear(a_loop, max_error))
                {
                    success = true;
                    break;
                }
            }
            if (!success) return false;
        }
        return true;
    }

    // Deletes the contents of the loops_ vector and resets the polygon state.
    private void ClearLoops()
    {
        ClearIndex();
        loops_.Clear();
        error_inconsistent_loop_orientations_ = false;
    }

    // Return true if there is an error in the loop nesting hierarchy.
    private bool FindLoopNestingError(out S2Error error)
    {
        // First check that the loop depths make sense.
        var last_depth = -1;
        for (var i = 0; i < NumLoops(); ++i)
        {
            var depth = Loop(i).Depth;
            if (depth < 0 || depth > last_depth + 1)
            {
                error = new(S2ErrorCode.POLYGON_INVALID_LOOP_DEPTH, $"Loop {i}: invalid loop depth ({depth})");
                return true;
            }
            last_depth = depth;
        }
        // Then check that they correspond to the actual loop nesting.  This test
        // is quadratic in the number of loops but the cost per iteration is small.
        for (var i = 0; i < NumLoops(); ++i)
        {
            var last = GetLastDescendant(i);
            for (var j = 0; j < NumLoops(); ++j)
            {
                if (i == j) continue;
                var nested = (j >= i + 1) && (j <= last);
                var reverse_b = false;
                if (Loop(i).ContainsNonCrossingBoundary(Loop(j), reverse_b) != nested)
                {
                    error = new(S2ErrorCode.POLYGON_INVALID_LOOP_NESTING,
                                $"Invalid nesting: loop {i} should {(nested ? "" : "not ")}contain loop {j}");
                    return true;
                }
            }
        }
        error = S2Error.OK;
        return false;
    }

    private static void InsertLoop(S2Loop new_loop, S2Loop parent, LoopMap loop_map)
    {
        List<S2Loop> children;
        bool done;
        do
        {
            children = loop_map.GetOrCreate(parent, () => new List<S2Loop>());
            done = true;
            foreach (var child in children)
            {
                if (child.ContainsNested(new_loop))
                {
                    parent = child;
                    done = false;
                    break;
                }
            }
        } while (!done);

        // Some of the children of the parent loop may now be children of
        // the new loop.
        var new_children = loop_map.GetOrCreate(new_loop, () => new List<S2Loop>());
        for (var i = 0; i < children.Count;)
        {
            var child = children[i];
            if (new_loop.ContainsNested(child))
            {
                new_children.Add(child);
                children.RemoveAt(i);
            }
            else
            {
                ++i;
            }
        }
        children.Add(new_loop);
    }

    private void InitLoops(LoopMap loop_map)
    {
        var loop_stack = new List<S2Loop> { S2Loop.NullLoop() };
        var depth = -1;
        while (loop_stack.Count!=0)
        {
            var loop = loop_stack.First();
            loop_stack.RemoveAt(0);
            if (loop != S2Loop.NullLoop())
            {
                depth = loop.Depth;
                loops_.Add(loop);
            }
            if (loop_map.TryGetValue(loop, out List<S2Loop>? value))
            {
                var children = value;
                for (var i = children.Count - 1; i >= 0; --i)
                {
                    var child = children[i];
                    //Assert.True(child is not null);
                    child!.Depth = depth + 1;
                    loop_stack.Add(child);
                }
            }
        }
    }

    // Add the polygon's loops to the S2ShapeIndex.  (The actual work of
    // building the index only happens when the index is first used.)
    private void InitIndex()
    {
        MyDebug.Assert(Index.NumShapeIds() == 0);
        Index.Add(new Shape(this));
#if s2polygon_not_lazy_indexing
            index_.ForceBuild();
#endif
#if s2debug
        // Note that s2debug is false in optimized builds (by default).
        MyDebug.Assert(IsValid());
#endif
    }

    // When the loop is modified (Invert(), or Init() called again) then the
    // indexing structures need to be cleared since they become invalid.
    private void ClearIndex()
    {
        unindexed_contains_calls_ = 0;
        Index.Clear();
    }

    // Initializes the polygon to the result of the given boolean operation,
    // returning an error on failure.
    private bool InitToOperation(S2BooleanOperation.OpType op_type, SnapFunction snap_function, S2Polygon a, S2Polygon b, out S2Error error)
    {
        var options = new S2BooleanOperation.Options
        {
            SnapFunction_ = snap_function
        };
        var op = new S2BooleanOperation(op_type, new S2PolygonLayer(this), options);
        return op.Build(a.Index, b.Index, out error);
    }

    // Initializes the polygon to the result of the given boolean operation,
    // logging an error on failure (fatal in debug builds).
    private void InitToOperation(S2BooleanOperation.OpType op_type, SnapFunction snap_function, S2Polygon a, S2Polygon b)
    {
        if (!InitToOperation(op_type, snap_function, a, b, out var error))
        {
            throw new ApplicationException(op_type + " operation failed: " + error);
        }
    }

    // Initializes the polygon from input polygon "a" using the given S2Builder.
    // If the result has an empty boundary (no loops), also decides whether the
    // result should be the full polygon rather than the empty one based on the
    // area of the input polygon.  (See comments in InitToApproxIntersection.)
    private void InitFromBuilder(S2Polygon a, S2Builder builder)
    {
        builder.StartLayer(new S2PolygonLayer(this));
        builder.AddPolygon(a);
        if (!builder.Build(out var error))
        {
            throw new ApplicationException("Could not build polygon: " + error);
        }
        // If there are no loops, check whether the result should be the full
        // polygon rather than the empty one.  (See InitToIntersection.)
        if (NumLoops() == 0)
        {
            if (a.bound_.Area() > S2.M_2_PI && a.GetArea() > S2.M_2_PI) Invert();
        }
    }

    private List<S2Polyline> OperationWithPolyline(S2BooleanOperation.OpType op_type, SnapFunction snap_function, S2Polyline a)
    {
        var options = new S2BooleanOperation.Options
        {
            SnapFunction_ = snap_function
        };
        var result = new List<S2Polyline>();
        S2PolylineVectorLayer.Options layer_options = new()
        {
            PolylineType_ = S2Builder.Graph.PolylineType.WALK,
        };
        S2BooleanOperation op = new(op_type, new S2PolylineVectorLayer(result, layer_options), options);
        MutableS2ShapeIndex a_index =
        [
            new S2Polyline.Shape(a)
        ];
        if (!op.Build(a_index, Index, out var error))
        {
            throw new ApplicationException("Polyline " + op_type + " operation failed: " + error);
        }
        return result;
    }

    // See comments in InitToSimplifiedInCell.
    private static S2Polyline[] SimplifyEdgesInCell(S2Polygon a, S2Cell cell, double tolerance_uv, S1Angle snap_radius)
    {
        var options = new S2Builder.Options(new IdentitySnapFunction(snap_radius))
        {
            SimplifyEdgeChains = true
        };
        var builder = new S2Builder(options);
        // The output consists of a sequence of polylines.  Polylines consisting of
        // interior edges are simplified using S2Builder, while polylines consisting
        // of boundary edges are returned unchanged.
        var polylines = new List<S2Polyline>();
        for (var i = 0; i < a.NumLoops(); ++i)
        {
            var a_loop = a.Loop(i);
            var v0 = a_loop.OrientedVertex(0);
            var mask0 = GetCellEdgeIncidenceMask(cell, v0, tolerance_uv);
            var in_interior = false;  // Was the last edge an interior edge?
            for (var j = 1; j <= a_loop.NumVertices; ++j)
            {
                var v1 = a_loop.OrientedVertex(j);
                var mask1 = GetCellEdgeIncidenceMask(cell, v1, tolerance_uv);
                if ((mask0 & mask1) != 0)
                {
                    // This is an edge along the cell boundary.  Such edges do not get
                    // simplified; we add them directly to the output.  (We create a
                    // separate polyline for each edge to keep things simple.)  We call
                    // ForceVertex on all boundary vertices to ensure that they don't
                    // move, and so that nearby interior edges are snapped to them.
                    MyDebug.Assert(!in_interior);
                    builder.ForceVertex(v1);
                    polylines.Add(new S2Polyline(new S2Point[] { v0, v1 }));
                }
                else
                {
                    // This is an interior edge.  If this is the first edge of an interior
                    // chain, then start a new S2Builder layer.  Also ensure that any
                    // polyline vertices on the boundary do not move, so that they will
                    // still connect with any boundary edge(s) there.
                    if (!in_interior)
                    {
                        var polyline = new S2Polyline();
                        builder.StartLayer(new S2PolylineLayer(polyline));
                        polylines.Add(polyline);
                        in_interior = true;
                    }
                    builder.AddEdge(v0, v1);
                    if (mask1 != 0)
                    {
                        builder.ForceVertex(v1);
                        in_interior = false;  // Terminate this polyline.
                    }
                }
                v0 = v1;
                mask0 = mask1;
            }
        }
        if (!builder.Build(out var error))
        {
            throw new ApplicationException("InitToSimplifiedInCell failed: " + error);
        }
        return [.. polylines];
    }

    // Internal implementation of intersect/subtract polyline functions above.
    // private S2Polyline[] InternalClipPolyline(bool invert, S2Polyline a, S1Angle snap_radius);

    // Defines a total ordering on S2Loops that does not depend on the cyclic
    // order of loop vertices.  This function is used to choose which loop to
    // invert in the case where several loops have exactly the same area.
    //
    // TODO(ericv): Consider adding this to the S2Loop API.  (May also want an
    // undirected version (CompareDirected vs CompareUndirected); should they
    // return a sign, or have separate "<" and "==" methods?)
    private static int CompareLoops(S2Loop a, S2Loop b)
    {
        if (a.NumVertices != b.NumVertices)
        {
            return a.NumVertices - b.NumVertices;
        }
        var ao = a.GetCanonicalLoopOrder();
        var bo = b.GetCanonicalLoopOrder();
        if (ao.Dir != bo.Dir) return ao.Dir - bo.Dir;
        for (int n = a.NumVertices, ai = ao.First, bi = bo.First;
             --n >= 0; ai += ao.Dir, bi += bo.Dir)
        {
            if (a.Vertex(ai) < b.Vertex(bi)) return -1;
            if (a.Vertex(ai) > b.Vertex(bi)) return 1;
        }
        return 0;
    }

    // Given a point "p" inside an S2Cell or on its boundary, return a mask
    // indicating which of the S2Cell edges the point lies on.  All boundary
    // comparisons are to within a maximum "u" or "v" error of "tolerance_uv".
    // Bit "i" in the result is set if and only "p" is incident to the edge
    // corresponding to S2Cell.edge(i).
    private static byte GetCellEdgeIncidenceMask(S2Cell cell, S2Point p, double tolerance_uv)
    {
        byte mask = 0;
        if (S2.FaceXYZtoUV(cell.Face, p, out var uv))
        {
            var bound = cell.BoundUV;
#if s2debug
            MyDebug.Assert(bound.Expanded(tolerance_uv).Contains(uv));
#endif
            if (Math.Abs(uv[1] - bound[1][0]) <= tolerance_uv) mask |= 1;
            if (Math.Abs(uv[0] - bound[0][1]) <= tolerance_uv) mask |= 2;
            if (Math.Abs(uv[1] - bound[1][1]) <= tolerance_uv) mask |= 4;
            if (Math.Abs(uv[0] - bound[0][0]) <= tolerance_uv) mask |= 8;
        }
        return mask;
    }

    // Returns the total number of bytes used by the polygon.
    public int SpaceUsed()
    {
        var size = Marshal.SizeOf(this);
        for (var i = 0; i < NumLoops(); ++i)
        {
            size += Loop(i).SpaceUsed();
        }
        size += Index.SpaceUsed() - Marshal.SizeOf(Index);
        return size;
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    // GetRectBound() returns essentially tight results, while GetCapBound()
    // might have a lot of extra padding.  Both bounds are conservative in that
    // if the loop contains a point P, then the bound contains P also.

    public S2Cap GetCapBound()  // Cap surrounding rect bound.
    {
        return bound_.GetCapBound();
    }
    public S2LatLngRect GetRectBound() { return bound_; }
    public void GetCellUnionBound(List<S2CellId> cell_ids)
    {
        Index.MakeS2ShapeIndexRegion().GetCellUnionBound(cell_ids);
    }

    public bool Contains(S2Cell cell)
    {
        return Index.MakeS2ShapeIndexRegion().Contains(cell);
    }
    public bool MayIntersect(S2Cell cell)
    {
        return Index.MakeS2ShapeIndexRegion().MayIntersect(cell);
    }

    // The point 'p' does not need to be normalized.
    public bool Contains(S2Point p)
    {
        // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
        // worthwhile only when it might allow us to delay building the index.
        if (!Index.IsFresh() && !bound_.Contains(p)) return false;

        // For small polygons it is faster to just check all the crossings.
        // Otherwise we keep track of the number of calls to Contains() and only
        // build the index once enough calls have been made so that we think it is
        // worth the effort.  See S2Loop.Contains(S2Point) for detailed comments.
        const int kMaxBruteForceVertices = 32;
        const int kMaxUnindexedContainsCalls = 20;
        if (NumVertices <= kMaxBruteForceVertices ||
            (!Index.IsFresh() &&
             ++unindexed_contains_calls_ != kMaxUnindexedContainsCalls))
        {
            var inside = false;
            for (var i = 0; i < NumLoops(); ++i)
            {
                // Use brute force to avoid building the loop's S2ShapeIndex.
                inside ^= Loop(i).BruteForceContains(p);
            }
            return inside;
        }
        // Otherwise we look up the S2ShapeIndex cell containing this point.
        return Index.MakeS2ContainsPointQuery().Contains(p);
    }

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2Polygon to "encoder".
    //
    // The encoding uses about 4 bytes per vertex for typical polygons in
    // Google's geographic repository, assuming that most vertices have been
    // snapped to the centers of S2Cells at some fixed level (typically using
    // InitToSnapped). The remaining vertices are stored using 24 bytes.
    // Decoding a polygon encoded this way always returns the original polygon,
    // without any loss of precision.
    //
    // The snap level is chosen to be the one that has the most vertices snapped
    // to S2Cells at that level.  If most vertices need 24 bytes, then all
    // vertices are encoded this way (this method automatically chooses the
    // encoding that has the best chance of giving the smaller output size).
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        if (NumVertices == 0)
        {
            EncodeCompressed(encoder, null, S2.kMaxCellLevel);
            return;
        }
        // Converts all the polygon vertices to S2XYZFaceSiTi format.
        var all_vertices = new S2PointCompression.S2XYZFaceSiTi[NumVertices];
        var current_loop_vertices = 0;
        foreach (var loop in loops_)
        {
            loop.GetXYZFaceSiTiVertices(all_vertices, current_loop_vertices);
            current_loop_vertices += loop.NumVertices;
        }
        // Computes a histogram of the cell levels at which the vertices are snapped.
        // cell_level is -1 for unsnapped, or 0 through kMaxCellLevel if snapped,
        // so we add one to it to get a non-negative index.  (histogram[0] is the
        // number of unsnapped vertices, histogram[i] the number of vertices
        // snapped at level i-1).
        var histogram = new int[S2.kMaxCellLevel + 2].Fill(0);
        foreach (var v in all_vertices)
        {
            histogram[v.CellLevel + 1] += 1;
        }
        // Compute the level at which most of the vertices are snapped.
        // If multiple levels have the same maximum number of vertices
        // snapped to it, the first one (lowest level number / largest
        // area / smallest encoding length) will be chosen, so this
        // is desired.  Start with histogram[1] since histogram[0] is
        // the number of unsnapped vertices, which we don't care about.
        var max_item = histogram.Skip(1).Max();
        var max_iter = Array.IndexOf(histogram, max_item);
        // snap_level will be at position histogram[snap_level + 1], see above.
        var snap_level = max_iter - 1;
        var num_snapped = histogram[max_iter];
        // Choose an encoding format based on the number of unsnapped vertices and a
        // rough estimate of the encoded sizes.

        // The compressed encoding requires approximately 4 bytes per vertex plus
        // "exact_point_size" for each unsnapped vertex (encoded as an S2Point plus
        // the index at which it is located).
        var exact_point_size = SizeHelper.SizeOf(typeof(S2Point)) + 2;
        var num_unsnapped = NumVertices - num_snapped;
        var compressed_size = 4 * NumVertices + exact_point_size * num_unsnapped;
        var lossless_size = SizeHelper.SizeOf(typeof(S2Point)) * NumVertices;
        if (compressed_size < lossless_size)
        {
            EncodeCompressed(encoder, all_vertices, snap_level);
        }
        else
        {
            EncodeUncompressed(encoder);
        }
    }

    // Encodes the polygon's S2Points directly as three doubles using
    // (40 + 43 * num_loops + 24 * num_vertices) bytes.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void EncodeUncompressed(Encoder encoder)
    {
        encoder.Ensure(10);  // Sufficient
        encoder.Put8(kCurrentUncompressedEncodingVersionNumber);
        // This code used to write "owns_loops_", so write "true" for compatibility.
        encoder.Put8(1);
        // Encode obsolete "has_holes_" field for backwards compatibility.
        var has_holes = 0;
        for (var i = 0; i < NumLoops(); ++i)
        {
            if (Loop(i).IsHole()) has_holes = 1;
        }
        encoder.Put8((byte)has_holes);
        encoder.Put32(loops_.Count);
        MyDebug.Assert(encoder.Avail() >= 0);

        for (var i = 0; i < NumLoops(); ++i)
        {
            Loop(i).Encode(encoder, CodingHint.FAST);
        }
        bound_.Encode(encoder, CodingHint.FAST);
    }

    // Decodes a polygon encoded with Encode().  Returns true on success.
    public static (bool, S2Polygon?) Decode(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte)) return (false, null);
        var version = decoder.Get8();
        return version switch
        {
            kCurrentUncompressedEncodingVersionNumber => DecodeUncompressed(decoder),
            kCurrentCompressedEncodingVersionNumber => DecodeCompressed(decoder),
            _ => (false, null),
        };
    }

    // Decodes a polygon by pointing the S2Loop vertices directly into the
    // decoder's memory buffer (which needs to persist for the lifetime of the
    // decoded S2Polygon).  It is much faster than Decode(), but requires that
    // all the polygon vertices were encoded exactly using 24 bytes per vertex.
    // This essentially requires that the polygon was not snapped beforehand to
    // a given S2Cell level; otherwise this method falls back to Decode().
    //
    // Returns true on success.
    public static (bool, S2Polygon?) DecodeWithinScope(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte)) return (false, null);
        var version = decoder.Get8();
        return version switch
        {
            kCurrentUncompressedEncodingVersionNumber => DecodeUncompressed(decoder),
            kCurrentCompressedEncodingVersionNumber => DecodeCompressed(decoder),
            _ => (false, null),
        };
    }

    // Decode a polygon encoded with EncodeUncompressed().  Used by the Decode
    // and DecodeWithinScope methods above.  The within_scope parameter
    // specifies whether to call DecodeWithinScope on the loops.
    private static (bool, S2Polygon?) DecodeUncompressed(Decoder decoder)
    {
        if (decoder.Avail() < 2 * sizeof(byte) + sizeof(uint))
            return (false, null);

        decoder.Get8();  // Ignore irrelevant serialized owns_loops_ value.
        decoder.Get8();  // Ignore irrelevant serialized has_holes_ value.
                         // Polygons with no loops are explicitly allowed here: a newly created
                         // polygon has zero loops and such polygons encode and decode properly.
        var num_loops = (int)decoder.Get32();
        if (num_loops > s2polygon_decode_max_num_loops) return (false, null);
        var loops = new List<S2Loop>(num_loops);
        for (var i = 0; i < num_loops; ++i)
        {
            var (success, newLoop) = S2Loop.Decode(decoder);
            if (!success) return (false, null);
            loops.Add(newLoop!);
        }
        var (success2, bound) = S2LatLngRect.Decode(decoder);
        if (!success2) return (false, null);
        var subregion_bound = S2LatLngRectBounder.ExpandForSubregions(bound);

        return (true, new S2Polygon(loops, bound, subregion_bound));
    }

    // Encode the polygon's vertices using about 4 bytes / vertex plus 24 bytes /
    // unsnapped vertex. All the loop vertices must be converted first to the
    // S2XYZFaceSiTi format using S2Loop.GetXYZFaceSiTiVertices, and concatenated
    // in the all_vertices array.
    //
    // REQUIRES: snap_level >= 0.
    private void EncodeCompressed(Encoder encoder, S2PointCompression.S2XYZFaceSiTi[]? all_vertices, int snap_level)
    {
        MyDebug.Assert(snap_level >= 0);
        // Sufficient for what we write. Typically enough for a 4 vertex polygon.
        encoder.Ensure(40);
        encoder.Put8(kCurrentCompressedEncodingVersionNumber);
        encoder.Put8((byte)snap_level);
        encoder.PutVarUInt32((uint)NumLoops());
        MyDebug.Assert(encoder.Avail() >= 0);
        var current_loop_vertices = 0;
        for (var i = 0; i < NumLoops(); ++i)
        {
            //if NumLoops() != 0 => all_vertices is not null
            loops_[i].EncodeCompressed(encoder, all_vertices!, current_loop_vertices, snap_level);
            current_loop_vertices += loops_[i].NumVertices;
        }
        // Do not write the bound or num_vertices as they can be cheaply recomputed
        // by DecodeCompressed.  Microbenchmarks show the speed difference is
        // inconsequential.
    }

    // Decode a polygon encoded with EncodeCompressed().
    private static (bool success, S2Polygon? result) DecodeCompressed(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte)) return (false, null);
        S2Polygon pol = new();
        pol.ClearLoops();
        var snap_level = decoder.Get8();
        if (snap_level > S2.kMaxCellLevel) return (false, null);
        // Polygons with no loops are explicitly allowed here: a newly created
        // polygon has zero loops and such polygons encode and decode properly.
        if (!decoder.TryGetVarUInt32(out var num_loops)) return (false, null);
        if (num_loops > s2polygon_decode_max_num_loops) return (false, null);
        pol.loops_.Capacity = (int)num_loops;
        for (var i = 0; i < num_loops; ++i)
        {
            var (success, loop) = S2Loop.DecodeCompressed(decoder, snap_level);
            if (!success)
            {
                return (false, null);
            }
            pol.loops_.Add(loop!);
        }
        pol.InitLoopProperties();
        return (true, pol);
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone()
    {
        var result = new S2Polygon();
        result.Copy(this);
        return result;
    }

    #endregion

    #region IEquatable

    // Return true if two polygons have exactly the same loops.  The loops must
    // appear in the same order, and corresponding loops must have the same
    // linear vertex ordering (i.e., cyclic rotations are not allowed).
    public bool Equals(S2Polygon? b)
    {
        if (b is null) return false;

        return loops_.SequenceEqual(b.loops_);
    }
    public override int GetHashCode() => LinqUtils.GetSequenceHashCode(loops_);

    #endregion

    #region IDisposable

    public void Dispose()
    {
        ClearLoops();
    }

    #endregion

    // Wrapper class for indexing a polygon (see S2ShapeIndex).  Once this
    // object is inserted into an S2ShapeIndex it is owned by that index, and
    // will be automatically deleted when no longer needed by the index.  Note
    // that this class does not take ownership of the polygon itself (see
    // OwningShape below).  You can also subtype this class to store additional
    // data (see S2Shape for details).
    //
    // Note that unlike S2Polygon, the edges of S2Polygon.Shape are directed
    // such that the polygon interior is always on the left.
    public class Shape : S2Shape
    {
        #region Fields, Constants

        public const TypeTag kTypeTag = TypeTag.S2Polygon;

        // The loop that contained the edge returned by the previous call to the
        // edge() method.  This is used as a hint to speed up edge location when
        // there are many loops.  Note that this field does not take up any space
        // due to field packing with S2Shape::id_.
        private int prev_loop_ = 0;

        public S2Polygon Polygon { get; protected set; }

        // An array where element "i" is the total number of edges in loops 0..i-1.
        // This field is only used for polygons that have a large number of loops.
        private readonly int[]? loop_starts_;

        #endregion

        #region Constructors

        public Shape(S2Polygon polygon)
        {
            Polygon = polygon;
            //inlined Init()
            loop_starts_ = null;
            int offset = 0;
            if (!polygon.IsFull())
            {
                var kMaxLinearSearchLoops = 12;  // From benchmarks.
                var num_loops = polygon.NumLoops();
                if (num_loops > kMaxLinearSearchLoops)
                {
                    // Unlike make_unique<>, new T[] does not default-construct each element.
                    loop_starts_ = new int[num_loops + 1];
                }
                for (var i = 0; i < num_loops; ++i)
                {
                    if (loop_starts_ is not null) loop_starts_[i] = offset;
                    offset += polygon.Loop(i).NumVertices;
                }
                if (loop_starts_ is not null) loop_starts_[num_loops] = offset;
            }
        }

        #endregion

        #region IEncoder

        public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
        {
            if (hint == CodingHint.FAST)
            {
                Polygon.EncodeUncompressed(encoder);
            }
            else
            {
                Polygon.Encode(encoder, hint);
            }
        }
        // Decoding is defined only for S2Polygon::OwningShape below.

        #endregion

        #region S2Shape

        // S2Shape interface:

        public sealed override int NumEdges()
        {
            return (Polygon.NumVertices != 1)
                ? Polygon.NumVertices
                : Polygon.IsFull() ? 0 : 1;
        }

        public sealed override Edge GetEdge(int e)
        {
            // Method names are fully specified to enable inlining.
            ChainPosition pos = GetChainPositionInternal(e);
            return ChainEdgeInternal(pos.ChainId, pos.Offset);
        }
        public sealed override int Dimension() => 2;
        public sealed override ReferencePoint GetReferencePoint()
        {
            var contains_origin = false;
            for (var i = 0; i < Polygon.NumLoops(); ++i)
            {
                contains_origin ^= Polygon.Loop(i).ContainsOrigin;
            }
            return new ReferencePoint(S2.Origin, contains_origin);
        }
        public sealed override int NumChains() => Polygon.NumLoops();
        public sealed override Chain GetChain(int i)
        {
            MyDebug.Assert(i < NumChains());
            if (loop_starts_ is not null)
            {
                var start = loop_starts_[i];
                return new Chain(start, loop_starts_[i + 1] - start);
            }
            else
            {
                var e = 0;
                for (var j = 0; j < i; ++j) e += Polygon.Loop(j).NumVertices;
                // S2Polygon represents a full loop as a loop with one vertex, while
                // S2Shape represents a full loop as a chain with no vertices.
                var num_vertices = Polygon.Loop(i).NumVertices;
                return new Chain(e, (num_vertices == 1) ? 0 : num_vertices);
            }
        }
        public sealed override Edge ChainEdge(int i, int j) => ChainEdgeInternal(i, j);
        private Edge ChainEdgeInternal(int i, int j)
        {
            MyDebug.Assert(i < NumChains());
            var loop = Polygon.Loop(i);
            MyDebug.Assert(j < loop.NumVertices);
            return new Edge(loop.OrientedVertex(j), loop.OrientedVertex(j + 1));
        }
        public sealed override ChainPosition GetChainPosition(int e) => GetChainPositionInternal(e);
        private ChainPosition GetChainPositionInternal(int e)
        {
            MyDebug.Assert(e < NumEdges());
            int i;
            if (loop_starts_ is null)
            {
                // When the number of loops is small, linear search is faster.  Most often
                // there is exactly one loop and the code below executes zero times.
                for (i = 0; e >= Polygon.Loop(i).NumVertices; ++i)
                {
                    e -= Polygon.Loop(i).NumVertices;
                }
            }
            else
            {
                i = prev_loop_; //.load(std::memory_order_relaxed);
                if (e >= loop_starts_[i] && e < loop_starts_[i + 1])
                {
                    // This edge belongs to the same loop as the previous call.
                }
                else
                {
                    if (e == loop_starts_[i + 1])
                    {
                        // This edge immediately follows the loop from the previous call.
                        // Note that S2Polygon does not allow empty loops.
                        ++i;
                    }
                    else
                    {
                        // "upper_bound" finds the loop just beyond the one we want.
                        i = LinqUtils.GetUpperBound(loop_starts_, e, loop_starts_[1], loop_starts_[Polygon.NumLoops()]) - loop_starts_[1];
                    }
                    prev_loop_ = i; //.store(i, std::memory_order_relaxed);
                }

                e -= loop_starts_[i];
            }
            return new ChainPosition(i, e);
        }
        public override TypeTag GetTypeTag() => kTypeTag;

        #endregion
    }

    // Like Shape, except that the S2Polygon is automatically deleted when this
    // object is deleted by the S2ShapeIndex.  This is useful when an S2Polygon
    // is constructed solely for the purpose of indexing it.
    public class OwningShape(S2Polygon polygon) : Shape(polygon), IInitEncoder<OwningShape>
    {
        public S2Polygon OwnedPolygon { get; set; } = polygon;

        public static (bool, OwningShape?) Init(Decoder decoder)
        {
            var (success, polygon) = Decode(decoder);
            if (!success) return (false, null);

            return (true, new OwningShape(polygon!));
        }
    }
}
