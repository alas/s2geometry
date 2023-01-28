namespace S2Geometry;

using System.Runtime.InteropServices;
using S2BuilderUtil;
using static S2PointCompression;

// This struct represents a search state in the NearlyCovers algorithm
// below.  See the description of the algorithm for details.
// IComparable Needed for storing SearchStates in a SortedSet.  The ordering
// chosen has no special meaning.
using SearchState = Key3<int, int, bool>;

// An S2Polyline represents a sequence of zero or more vertices connected by
// straight edges (geodesics).  Edges of length 0 and 180 degrees are not
// allowed, i.e. adjacent vertices should not be identical or antipodal.
public record struct S2Polyline : IS2Region<S2Polyline>, IDecoder<S2Polyline>
{
    #region Fields, Constants

    public const byte kCurrentCompressedEncodingVersionNumber = 2;

    public S2Point[] Vertices { get; set; }

    // Allows overriding the automatic validity checking controlled by the
    // --s2debug flag.
    public S2Debug s2debug_override_ { get; set; } 

    #endregion

    #region Constructors

    public S2Polyline() : this(Array.Empty<S2Point>()) { }

    // Convenience constructor that call Init() with the given vertices.
    public S2Polyline(S2Point[] vertices, S2Debug s2debug_override = S2Debug.ALLOW)
    {
        Vertices = vertices;
        s2debug_override_ = s2debug_override;
        //Init();
    //}

    // Initialize a polyline that connects the given vertices. Empty polylines are
    // allowed.  Adjacent vertices should not be identical or antipodal.  All
    // vertices should be unit length.
    //public void Init()
    //{
#if s2debug
        if (s2debug_override_ == S2Debug.ALLOW)
        {
            // Note that s2debug is false in optimized builds (by default).
            MyDebug.Assert(IsValid());
        }
#endif
    }

    // Convenience constructor that call Init() with the given vertices.
    public S2Polyline(S2LatLng[] vertices)
        : this(vertices.Select(t => t.ToPoint()).ToArray()) { }

    private S2Polyline(S2Polyline src)
    {
        Vertices = (S2Point[])src.Vertices.Clone();
        s2debug_override_ = src.s2debug_override_;
        // TODO(user, ericv): Decide whether to use a canonical empty
        // representation.
        // If num_vertices_ == 0, then src.vertices_ will be null if src was default
        // constructed.
    }

    #endregion

    #region S2Polyline

    // Convenience function that snaps the vertices to S2CellId centers at the
    // given level (default level 30, which has S2CellId centers spaced about 1
    // centimeter apart)
    public void InitToSnapped(S2Polyline polyline, int snap_level = S2CellId.kMaxLevel)
    {
        S2Builder builder =new(new(new S2CellIdSnapFunction(snap_level)));
        InitFromBuilder(polyline, builder);
    }

    // Snaps the input polyline according to the given "snap_function" and
    // reduces the number of vertices if possible, while ensuring that no vertex
    // moves further than snap_function.snap_radius().
    //
    // Simplification works by replacing nearly straight chains of short edges
    // with longer edges, in a way that preserves the topology of the input
    // polyline up to the creation of degeneracies.
    public void InitToSimplified(S2Polyline polyline, SnapFunction snap_function)
    {
        S2Builder.Options options=new(snap_function);
        options.SimplifyEdgeChains = true;
        S2Builder builder=new(options);
        InitFromBuilder(polyline, builder);
    }

    // Initializes this polyline from input polyline using the given S2Builder.
    public void InitFromBuilder(S2Polyline polyline, S2Builder builder)
    {
        builder.StartLayer(new S2PolylineLayer(this));
        builder.AddPolyline(polyline);
        S2Error error;
        MyDebug.Assert(builder.Build(out error), "Could not build polyline: " + error);
    }

    // Return true if the given vertices form a valid polyline.
    public bool IsValid()
    {
        if (FindValidationError(out S2Error error))
        {
#if s2debug
            MyDebug.WriteLine($"IsValid: {error}");
#endif
            return false;
        }
        return true;
    }

    // Returns true if this is *not* a valid polyline and sets "error"
    // appropriately.  Otherwise returns false and leaves "error" unchanged.
    //
    // REQUIRES: error != null
    public bool FindValidationError(out S2Error error)
    {
        // All vertices must be unit length.
        for (int i = 0; i < Vertices.Length; ++i)
        {
            if (!Vertex(i).IsUnitLength())
            {
                error = new(S2ErrorCode.NOT_UNIT_LENGTH, $"Vertex {i} is not unit length");
                return true;
            }
        }
        // Adjacent vertices must not be identical or antipodal.
        for (int i = 1; i < Vertices.Length; ++i)
        {
            if (Vertex(i - 1) == Vertex(i))
            {
                error = new(S2ErrorCode.DUPLICATE_VERTICES, $"Vertices {i - 1} and {i} are identical");
                return true;
            }
            if (Vertex(i - 1) == -Vertex(i))
            {
                error = new(S2ErrorCode.ANTIPODAL_VERTICES, $"Vertices {i - 1} and {i} are antipodal");
                return true;
            }
        }
        error = S2Error.OK;
        return false;
    }

    public int NumVertices() => Vertices?.Length ?? 0;

    public S2Point[] VerticesClone() => (S2Point[])Vertices.Clone();

    public S2Point Vertex(int k)
    {
        MyDebug.Assert(k >= 0);
        MyDebug.Assert(k <= Vertices.Length);
        return Vertices[k];
    }

    // Returns an S2PointSpan containing the polyline's vertices.
    public S2PointSpan GetVerticesSpan() => new(Vertices);

    // Return the length of the polyline.
    public S1Angle GetLength() => S2PolylineMeasures.GetLength(GetVerticesSpan());

    // Return the true centroid of the polyline multiplied by the length of the
    // polyline (see s2centroids.h for details on centroids).  The result is not
    // unit length, so you may want to normalize it.
    //
    // Prescaling by the polyline length makes it easy to compute the centroid
    // of several polylines (by simply adding up their centroids).
    public S2Point GetCentroid() => S2PolylineMeasures.GetCentroid(GetVerticesSpan());

    // If all of the polyline's vertices happen to be the centers of S2Cells at
    // some level, then returns that level, otherwise returns -1.  See also
    // InitToSnapped() and s2builderutil::S2CellIdSnapFunction.
    // Returns -1 if the polyline has no vertices.
    public int GetSnapLevel()
    {
        int snap_level = -1;
        for (int i = 0; i < NumVertices(); ++i)
        {
            int face;
            uint si, ti;
            var level = S2.XYZtoFaceSiTi(Vertices[i], out face, out si, out ti);
            if (level< 0) return level;  // Vertex is not a cell center.
            if (level != snap_level)
            {
                if (snap_level< 0)
                {
                    snap_level = level;  // First vertex.
                }
                else
                {
                    return -1;  // Vertices at more than one cell level.
                }
            }
        }
        return snap_level;
    }


    // Return the point whose distance from vertex 0 along the polyline is the
    // given fraction of the polyline's total length.  Fractions less than zero
    // or greater than one are clamped.  The return value is unit length.  This
    // cost of this function is currently linear in the number of vertices.
    // The polyline must not be empty.
    public S2Point Interpolate(double fraction) => GetSuffix(fraction, out _);

    // Like Interpolate(), but also return the index of the next polyline
    // vertex after the interpolated point P.  This allows the caller to easily
    // construct a given suffix of the polyline by concatenating P with the
    // polyline vertices starting at "next_vertex".  Note that P is guaranteed
    // to be different than vertex(*next_vertex), so this will never result in
    // a duplicate vertex.
    //
    // The polyline must not be empty.  Note that if "fraction" >= 1.0, then
    // "next_vertex" will be set to Vertices.Length (indicating that no vertices
    // from the polyline need to be appended).  The value of "next_vertex" is
    // always between 1 and Vertices.Length.
    //
    // This method can also be used to construct a prefix of the polyline, by
    // taking the polyline vertices up to "next_vertex - 1" and appending the
    // returned point P if it is different from the last vertex (since in this
    // case there is no guarantee of distinctness).
    public S2Point GetSuffix(double fraction, out int next_vertex)
    {
        MyDebug.Assert(Vertices.Length > 0);
        // We intentionally let the (fraction >= 1) case fall through, since
        // we need to handle it in the loop below in any case because of
        // possible roundoff errors.
        if (fraction <= 0)
        {
            next_vertex = 1;
            return Vertex(0);
        }
        var length_sum = S1Angle.Zero;
        for (int i = 1; i < Vertices.Length; ++i)
        {
            length_sum += new S1Angle(Vertex(i - 1), Vertex(i));
        }
        S1Angle target = fraction * length_sum;
        for (int i = 1; i < Vertices.Length; ++i)
        {
            var length = new S1Angle(Vertex(i - 1), Vertex(i));
            if (target < length)
            {
                // This interpolates with respect to arc length rather than
                // straight-line distance, and produces a unit-length result.
                S2Point result = S2.GetPointOnLine(Vertex(i - 1), Vertex(i), target);
                // It is possible that (result == vertex(i)) due to rounding errors.
                next_vertex = (result == Vertex(i)) ? (i + 1) : i;
                return result;
            }
            target -= length;
        }
        next_vertex = Vertices.Length;
        return Vertex(Vertices.Length - 1);
    }

    // The inverse operation of GetSuffix/Interpolate.  Given a point on the
    // polyline, returns the ratio of the distance to the point from the
    // beginning of the polyline over the length of the polyline.  The return
    // value is always betwen 0 and 1 inclusive.  See GetSuffix() for the
    // meaning of "next_vertex".
    //
    // The polyline should not be empty.  If it has fewer than 2 vertices, the
    // return value is zero.
    public double UnInterpolate(S2Point point, int next_vertex)
    {
        MyDebug.Assert(Vertices.Length > 0);
        if (Vertices.Length < 2)
        {
            return 0;
        }
        var length_sum = S1Angle.Zero;
        for (int i = 1; i < next_vertex; ++i)
        {
            length_sum += new S1Angle(Vertex(i - 1), Vertex(i));
        }
        S1Angle length_to_point = length_sum + new S1Angle(Vertex(next_vertex - 1), point);
        for (int i = next_vertex; i < Vertices.Length; ++i)
        {
            length_sum += new S1Angle(Vertex(i - 1), Vertex(i));
        }
        // The ratio can be greater than 1.0 due to rounding errors or because the
        // point is not exactly on the polyline.
        return Math.Min(1.0, length_to_point / length_sum);
    }

    // Given a point, returns a point on the polyline that is closest to the given
    // point.  See GetSuffix() for the meaning of "next_vertex", which is chosen
    // here w.r.t. the projected point as opposed to the interpolated point in
    // GetSuffix().
    //
    // The polyline must be non-empty.
    public S2Point Project(S2Point point, out int next_vertex)
    {
        MyDebug.Assert(Vertices.Length > 0);

        if (Vertices.Length == 1)
        {
            // If there is only one vertex, it is always closest to any given point.
            next_vertex = 1;
            return Vertex(0);
        }

        // Initial value larger than any possible distance on the unit sphere.
        S1Angle min_distance = S1Angle.FromRadians(10);
        int min_index = -1;

        // Find the line segment in the polyline that is closest to the point given.
        for (int i = 1; i < Vertices.Length; ++i)
        {
            var distance_to_segment = S2.GetDistance(point, Vertex(i - 1), Vertex(i));
            if (distance_to_segment < min_distance)
            {
                min_distance = distance_to_segment;
                min_index = i;
            }
        }
        MyDebug.Assert(min_index != -1);

        // Compute the point on the segment found that is closest to the point given.
        var closest_point = S2.Project(point, Vertex(min_index - 1), Vertex(min_index));

        next_vertex = min_index + (closest_point == Vertex(min_index) ? 1 : 0);
        return closest_point;
    }

    // Returns true if the point given is on the right hand side of the polyline,
    // using a naive definition of "right-hand-sideness" where the point is on
    // the RHS of the polyline iff the point is on the RHS of the line segment in
    // the polyline which it is closest to.
    //
    // The polyline must have at least 2 vertices.
    public bool IsOnRight(S2Point point)
    {
        MyDebug.Assert(Vertices.Length >= 2);

        S2Point closest_point = Project(point, out int next_vertex);

        MyDebug.Assert(next_vertex >= 1);
        MyDebug.Assert(next_vertex <= Vertices.Length);

        // If the closest point C is an interior vertex of the polyline, let B and D
        // be the previous and next vertices.  The given point P is on the right of
        // the polyline (locally) if B, P, D are ordered CCW around vertex C.
        if (closest_point == Vertex(next_vertex - 1) && next_vertex > 1 && next_vertex < Vertices.Length)
        {
            if (point == Vertex(next_vertex - 1))
                return false;  // Polyline vertices are not on the RHS.
            return S2Pred.OrderedCCW(
                Vertex(next_vertex - 2), point,
                Vertex(next_vertex),
                Vertex(next_vertex - 1));
        }

        // Otherwise, the closest point C is incident to exactly one polyline edge.
        // We test the point P against that edge.
        if (next_vertex == Vertices.Length)
            --next_vertex;

        return S2Pred.Sign(point, Vertex(next_vertex), Vertex(next_vertex - 1)) > 0;
    }

    // Return true if this polyline intersects the given polyline. If the
    // polylines share a vertex they are considered to be intersecting. When a
    // polyline endpoint is the only intersection with the other polyline, the
    // function may return true or false arbitrarily.
    //
    // The running time is quadratic in the number of vertices.  (To intersect
    // polylines more efficiently, or compute the actual intersection geometry,
    // use S2BooleanOperation.)
    public bool Intersects(S2Polyline line)
    {
        if (Vertices.Length <= 0 || line.Vertices.Length <= 0)
        {
            return false;
        }

        if (!GetRectBound().Intersects(line.GetRectBound()))
        {
            return false;
        }

        // TODO(ericv): Use S2ShapeIndex here.
        for (int i = 1; i < Vertices.Length; ++i)
        {
            var crosser = new S2EdgeCrosser(Vertex(i - 1), Vertex(i), line.Vertex(0));
            for (int j = 1; j < line.Vertices.Length; ++j)
            {
                if (crosser.CrossingSign(line.Vertex(j)) >= 0)
                {
                    return true;
                }
            }
        }
        return false;
    }

    // Return a subsequence of vertex indices such that the polyline connecting
    // these vertices is never further than "tolerance" from the original
    // polyline.  Provided the first and last vertices are distinct, they are
    // always preserved; if they are not, the subsequence may contain only a
    // single index.
    //
    // Some useful properties of the algorithm:
    //
    //  - It runs in linear time.
    //
    //  - The output is always a valid polyline.  In particular, adjacent
    //    output vertices are never identical or antipodal.
    //
    //  - The method is not optimal, but it tends to produce 2-3% fewer
    //    vertices than the Douglas-Peucker algorithm with the same tolerance.
    //
    //  - The output is *parametrically* equivalent to the original polyline to
    //    within the given tolerance.  For example, if a polyline backtracks on
    //    itself and then proceeds onwards, the backtracking will be preserved
    //    (to within the given tolerance).  This is different than the
    //    Douglas-Peucker algorithm, which only guarantees geometric equivalence.
    //
    // See also S2PolylineSimplifier, which uses the same algorithm but is more
    // efficient and supports more features, and also S2Builder, which can
    // simplify polylines and polygons, supports snapping (e.g. to E7 lat/lng
    // coordinates or S2CellId centers), and can split polylines at intersection
    // points.
    public void SubsampleVertices(S1Angle tolerance, out int[]? indices)
    {
        var indicesTmp = new InputEdgeLoop();
        if (Vertices.Length == 0)
        {
            indices = null;
            return;
        }

        indicesTmp.Add(0);
        var clamped_tolerance = S1Angle.Max(tolerance, S1Angle.FromRadians(0));
        for (int index = 0; index + 1 < Vertices.Length;)
        {
            int next_index = FindEndVertex(this, clamped_tolerance, index);
            // Don't create duplicate adjacent vertices.
            if (Vertex(next_index) != Vertex(index))
            {
                indicesTmp.Add(next_index);
            }
            index = next_index;
        }
        indices = indicesTmp.ToArray();
    }

    // Given a polyline, a tolerance distance, and a start index, this function
    // returns the maximal end index such that the line segment between these two
    // vertices passes within "tolerance" of all interior vertices, in order.
    private static int FindEndVertex(S2Polyline polyline, S1Angle tolerance, int index)
    {
        MyDebug.Assert(tolerance.Radians >= 0);
        MyDebug.Assert((index + 1) < polyline.Vertices.Length);

        // The basic idea is to keep track of the "pie wedge" of angles from the
        // starting vertex such that a ray from the starting vertex at that angle
        // will pass through the discs of radius "tolerance" centered around all
        // vertices processed so far.

        // First we define a "coordinate frame" for the tangent and normal spaces
        // at the starting vertex.  Essentially this means picking three
        // orthonormal vectors X,Y,Z such that X and Y span the tangent plane at
        // the starting vertex, and Z is "up".  We use the coordinate frame to
        // define a mapping from 3D direction vectors to a one-dimensional "ray
        // angle" in the range (-Pi, Pi].  The angle of a direction vector is
        // computed by transforming it into the X,Y,Z basis, and then calculating
        // atan2(y,x).  This mapping allows us to represent a wedge of angles as a
        // 1D interval.  Since the interval wraps around, we represent it as an
        // S1Interval, i.e. an interval on the unit circle.
        var origin = polyline.Vertex(index);
        var frame = S2.GetFrame(origin);

        // As we go along, we keep track of the current wedge of angles and the
        // distance to the last vertex (which must be non-decreasing).
        S1Interval current_wedge = S1Interval.Full;
        double last_distance = 0;

        for (++index; index < polyline.Vertices.Length; index++)
        {
            var candidate = polyline.Vertex(index);
            double distance = origin.Angle(candidate);

            // We don't allow simplification to create edges longer than 90 degrees,
            // to avoid numeric instability as lengths approach 180 degrees.  (We do
            // need to allow for original edges longer than 90 degrees, though.)
            if (distance > S2.M_PI_2 && last_distance > 0) break;

            // Vertices must be in increasing order along the ray, except for the
            // initial disc around the origin.
            if (distance < last_distance && last_distance > tolerance.Radians) break;
            last_distance = distance;

            // Points that are within the tolerance distance of the origin do not
            // constrain the ray direction, so we can ignore them.
            if (distance <= tolerance.Radians) continue;

            // If the current wedge of angles does not contain the angle to this
            // vertex, then stop right now.  Note that the wedge of possible ray
            // angles is not necessarily empty yet, but we can't continue unless we
            // are willing to backtrack to the last vertex that was contained within
            // the wedge (since we don't create new vertices).  This would be more
            // complicated and also make the worst-case running time more than linear.
            S2Point direction = S2.ToFrame(frame, candidate);
            double center = Math.Atan2(direction.Y, direction.X);
            if (!current_wedge.Contains(center)) break;

            // To determine how this vertex constrains the possible ray angles,
            // consider the triangle ABC where A is the origin, B is the candidate
            // vertex, and C is one of the two tangent points between A and the
            // spherical cap of radius "tolerance" centered at B.  Then from the
            // spherical law of sines, sin(a)/sin(A) = sin(c)/sin(C), where "a" and
            // "c" are the lengths of the edges opposite A and C.  In our case C is a
            // 90 degree angle, therefore A = asin(sin(a) / sin(c)).  Angle A is the
            // half-angle of the allowable wedge.

            double half_angle = Math.Asin(Math.Sin(tolerance.Radians) / Math.Sin(distance));
            S1Interval target = S1Interval.FromPoint(center).Expanded(half_angle);
            current_wedge = current_wedge.Intersection(target);
            MyDebug.Assert(!current_wedge.IsEmpty());
        }
        // We break out of the loop when we reach a vertex index that can't be
        // included in the line segment, so back up by one vertex.
        return index - 1;
    }

    // Return true if two polylines have the same number of vertices, and
    // corresponding vertex pairs are separated by no more than "max_error".
    // (For testing purposes.)
    public bool ApproxEquals(S2Polyline b) => ApproxEquals(b, S1Angle.FromRadians(S2.DoubleError));

    public bool ApproxEquals(S2Polyline b, S1Angle max_error)
    {
        if (Vertices.Length != b.Vertices.Length) return false;
        for (int offset = 0; offset < Vertices.Length; ++offset)
        {
            if (!S2.ApproxEquals(Vertex(offset), b.Vertex(offset), max_error))
            {
                return false;
            }
        }
        return true;
    }

    // Return true if "covered" is within "max_error" of a contiguous subpath of
    // this polyline over its entire length.  Specifically, this method returns
    // true if this polyline has parameterization a:[0,1] . S^2, "covered" has
    // parameterization b:[0,1] . S^2, and there is a non-decreasing function
    // f:[0,1] . [0,1] such that distance(a(f(t)), b(t)) <= max_error for all t.
    //
    // You can think of this as testing whether it is possible to drive a car
    // along "covered" and a car along some subpath of this polyline such that no
    // car ever goes backward, and the cars are always within "max_error" of each
    // other.
    //
    // This function is well-defined for empty polylines:
    //    anything.covers(empty) = true
    //    empty.covers(nonempty) = false
    public bool NearlyCovers(S2Polyline covered, S1Angle max_error)
    {
        // NOTE: This algorithm is described assuming that adjacent vertices in a
        // polyline are never at the same point.  That is, the ith and i+1th vertices
        // of a polyline are never at the same point in space.  The implementation
        // does not make this assumption.

        // DEFINITIONS:
        //   - edge "i" of a polyline is the edge from the ith to i+1th vertex.
        //   - covered_j is a polyline consisting of edges 0 through j of "covered."
        //   - this_i is a polyline consisting of edges 0 through i of this polyline.
        //
        // A search state is represented as an (int, int, bool) tuple, (i, j,
        // i_in_progress).  Using the "drive a car" analogy from the header comment, a
        // search state signifies that you can drive one car along "covered" from its
        // first vertex through a point on its jth edge, and another car along this
        // polyline from some point on or before its ith edge to a to a point on its
        // ith edge, such that no car ever goes backward, and the cars are always
        // within "max_error" of each other.  If i_in_progress is true, it means that
        // you can definitely drive along "covered" through the jth vertex (beginning
        // of the jth edge). Otherwise, you can definitely drive along "covered"
        // through the point on the jth edge of "covered" closest to the ith vertex of
        // this polyline.
        //
        // The algorithm begins by finding all edges of this polyline that are within
        // "max_error" of the first vertex of "covered," and adding search states
        // representing all of these possible starting states to the stack of
        // "pending" states.
        //
        // The algorithm proceeds by popping the next pending state,
        // (i,j,i_in_progress), off of the stack.  First it checks to see if that
        // state represents finding a valid covering of "covered" and returns true if
        // so.  Next, if the state represents reaching the end of this polyline
        // without finding a successful covering, the algorithm moves on to the next
        // state in the stack.  Otherwise, if state (i+1,j,false) is valid, it is
        // added to the stack of pending states.  Same for state (i,j+1,true).
        //
        // We need the stack because when "i" and "j" can both be incremented,
        // sometimes only one choice leads to a solution.  We use a set to keep track
        // of visited states to avoid duplicating work.  With the set, the worst-case
        // number of states examined is O(n+m) where n = this.Vertices.Length and m =
        // covered.Vertices.Length.  Without it, the amount of work could be as high as
        // O((n*m)^2).  Using set, the running time is O((n*m) log (n*m)).
        //
        // TODO(user): Benchmark this, and see if the set is worth it.

        if (covered.Vertices.Length == 0) return true;
        if (this.Vertices.Length == 0) return false;

        var pending = new List<SearchState>();
        var done = new SortedSet<SearchState>();

        // Find all possible starting states.
        for (int i = 0, next_i = NextDistinctVertex(this, 0), next_next_i;
             next_i < this.Vertices.Length; i = next_i, next_i = next_next_i)
        {
            next_next_i = NextDistinctVertex(this, next_i);
            S2Point closest_point = S2.Project(
                covered.Vertex(0), Vertex(i), Vertex(next_i));

            // In order to avoid duplicate starting states, we exclude the end vertex
            // of each edge *except* for the last non-degenerate edge.
            if ((next_next_i == this.Vertices.Length ||
                 closest_point != Vertex(next_i)) &&
                new S1Angle(closest_point, covered.Vertex(0)) <= max_error)
            {
                pending.Add(new SearchState(i, 0, true));
            }
        }

        while (pending.Any())
        {
            SearchState state = pending.Last();
            pending.RemoveAt(pending.Count - 1);
            if (!done.Add(state)) continue;

            int next_i = NextDistinctVertex(this, state.Item1);
            int next_j = NextDistinctVertex(covered, state.Item2);
            if (next_j == covered.Vertices.Length)
            {
                return true;
            }
            else if (next_i == this.Vertices.Length)
            {
                continue;
            }

            S2Point i_begin, j_begin;
            if (state.Item3)
            {
                j_begin = covered.Vertex(state.Item2);
                i_begin = S2.Project(
                    j_begin, Vertex(state.Item1), Vertex(next_i));
            }
            else
            {
                i_begin = Vertex(state.Item1);
                j_begin = S2.Project(
                    i_begin, covered.Vertex(state.Item2), covered.Vertex(next_j));
            }

            if (S2.IsEdgeBNearEdgeA(j_begin, covered.Vertex(next_j),
                                     i_begin, Vertex(next_i), max_error))
            {
                pending.Add(new SearchState(next_i, state.Item2, false));
            }
            if (S2.IsEdgeBNearEdgeA(i_begin, Vertex(next_i),
                                     j_begin, covered.Vertex(next_j), max_error))
            {
                pending.Add(new SearchState(state.Item1, next_j, true));
            }
        }
        return false;
    }

    // Returns the total number of bytes used by the polyline.
    public int SpaceUsed() =>
        Marshal.SizeOf(typeof(S2Polyline)) + Vertices.Length * SizeHelper.SizeOf(typeof(S2Point));

    // Return the first i > "index" such that the ith vertex of "pline" is not at
    // the same point as the "index"th vertex.  Returns pline.num_vertices() if
    // there is no such value of i.
    private static int NextDistinctVertex(S2Polyline pline, int index)
    {
        S2Point initial = pline.Vertex(index);
        do
        {
            ++index;
        } while (index < pline.Vertices.Length && pline.Vertex(index) == initial);
        return index;
    }

    public void Reverse()
    {
        if (NumVertices() > 0)
        {
            Array.Reverse(Vertices);
        }
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound() => GetRectBound().GetCapBound();
    public S2LatLngRect GetRectBound()
    {
        var bounder = new S2LatLngRectBounder();
        for (int i = 0; i < Vertices.Length; ++i)
        {
            bounder.AddPoint(Vertex(i));
        }
        return bounder.GetBound();
    }
    public bool Contains(S2Cell cell) => false;
    public bool MayIntersect(S2Cell cell)
    {
        if (Vertices.Length == 0) return false;

        // We only need to check whether the cell contains vertex 0 for correctness,
        // but these tests are cheap compared to edge crossings so we might as well
        // check all the vertices.
        for (int i = 0; i < Vertices.Length; ++i)
        {
            if (cell.Contains(Vertex(i))) return true;
        }
        var cell_vertices = new S2Point[4];
        for (int i = 0; i < 4; ++i)
        {
            cell_vertices[i] = cell.Vertex(i);
        }
        for (int j = 0; j < 4; ++j)
        {
            var crosser = new S2EdgeCrosser(cell_vertices[j], cell_vertices[(j + 1) & 3], Vertex(0));
            for (int i = 1; i < Vertices.Length; ++i)
            {
                if (crosser.CrossingSign(Vertex(i)) >= 0)
                {
                    // There is a proper crossing, or two vertices were the same.
                    return true;
                }
            }
        }
        return false;
    }

    // Always return false, because "containment" is not numerically
    // well-defined except at the polyline vertices.
    public bool Contains(S2Point p) => false;

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2Polyline to "encoder".
    // Currently just calls EncodeUncompressed().
    //
    // TODO(b/128865764): After transition period, replace it with
    // EncodeMostCompact.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        EncodeUncompressed(encoder);
    }

    // Appends a serialized uncompressed representation of the S2Polyline to
    // "encoder".
    public void EncodeUncompressed(Encoder encoder)
    {
        encoder.Ensure(Vertices.Length * SizeHelper.SizeOf(typeof(S2Point)) + 10);  // sufficient

        encoder.Put8(S2.kCurrentLosslessEncodingVersionNumber);
        encoder.Put32(Vertices.Length);
        encoder.PutPoints(Vertices);

        MyDebug.Assert(encoder.Avail() >= 0);
    }

    // Encode the polylines's vertices using the most compact way: compressed or
    // uncompressed.
    public void EncodeMostCompact(Encoder encoder) 
    {
        if (NumVertices() == 0)
        {
            EncodeCompressed(encoder, Array.Empty<S2XYZFaceSiTi>(), S2.kMaxCellLevel);
            return;
        }
        
        // Convert S2Points to (face, si, ti) representation.
        var all_vertices = new S2XYZFaceSiTi[NumVertices()];
        for (int i = 0; i < NumVertices(); ++i) {
            var xyz = Vertices[i];
            var cell_level =  S2.XYZtoFaceSiTi(xyz, out var face, out var si, out var ti);
            all_vertices[i] = new(xyz, face, si, ti, cell_level);
        }

        // Computes a histogram of the cell levels at which the vertices are snapped.
        // cell_level is -1 for unsnapped, or 0 through kMaxCellLevel if snapped,
        // so we add one to it to get a non-negative index.  (histogram[0] is the
        // number of unsnapped vertices, histogram[i] the number of vertices
        // snapped at level i-1).
        var histogram = new int[S2.kMaxCellLevel + 2];
        histogram.Fill(0);
        for (int i = 0; i < NumVertices(); ++i)
        {
            histogram[all_vertices[i].CellLevel + 1] += 1;
        }
        // Compute the level at which most of the vertices are snapped.
        // If multiple levels have the same maximum number of vertices
        // snapped to it, the first one (lowest level number / largest
        // area / smallest encoding length) will be chosen, so this
        // is desired.  Start with histogram[1] since histogram[0] is
        // the number of unsnapped vertices, which we don't care about.
        var max_iter = histogram.Select((n, i) => (n, i)).Skip(1).Max().i;

        // snap_level will be at position histogram[snap_level + 1], see above.
        var snap_level = max_iter - (max_iter + 1);
        var num_snapped = histogram[max_iter];

        // The compressed encoding requires approximately 4 bytes per vertex plus
        // "exact_point_size" for each unsnapped vertex (encoded as an S2Point plus
        // the index at which it is located).
        int exact_point_size = SizeHelper.SizeOf(typeof(S2Point)) + 2;
        int num_unsnapped = NumVertices() - num_snapped;
        int compressed_size = 4 * NumVertices() + exact_point_size * num_unsnapped;
        int lossless_size = SizeHelper.SizeOf(typeof(S2Point)) * NumVertices();
        if (compressed_size < lossless_size)
        {
            EncodeCompressed(encoder, all_vertices, snap_level);
        }
        else
        {
            EncodeUncompressed(encoder);
        }
    }

    // Encode the polylines's vertices using S2EncodePointsCompressed().
    public void EncodeCompressed(Encoder encoder, S2XYZFaceSiTi[] all_vertices, int snap_level)
    {
        // Set version number.
        encoder.Ensure(2 + Encoder.kVarintMax32);
        encoder.Put8(kCurrentCompressedEncodingVersionNumber);
        encoder.Put8((sbyte)snap_level);
        encoder.PutVarInt32(NumVertices());
        S2EncodePointsCompressed(all_vertices, snap_level, encoder);
    }

    // Decodes an S2Polyline encoded with any of Encode*() methods. Returns true
    // on success.
    public static (bool, S2Polyline) Decode(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte) + sizeof(UInt32)) return (false, default);
        byte version = decoder.Get8();
        return version switch
        {
            S2.kCurrentLosslessEncodingVersionNumber => DecodeUncompressed(decoder),
            kCurrentCompressedEncodingVersionNumber => DecodeCompressed(decoder),
            _ => (false, default),
        };
    }

    // Decode a polyline encoded with EncodeUncompressed(). Used by the Decode
    // method above.
    public static (bool, S2Polyline) DecodeUncompressed(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(UInt32)) return (false, default);

        var count = (int)decoder.Get32();

        // Check the bytes available before allocating memory in case of
        // corrupt/malicious input.
        if (decoder.Avail() < count * SizeHelper.SizeOf(typeof(S2Point))) return (false, default);
        var verts = new S2Point[count];
        decoder.GetPoints(verts, 0, count);

        var pol = new S2Polyline(verts);
#if s2debug
        //if (s2debug_override_ == S2Debug.ALLOW)
        {
            MyDebug.Assert(pol.IsValid());
        }
#endif
        return (true, pol);
    }

    // Decodes a polyline encoded with EncodeCompressed(). Used by the Decode
    // method above.
    public static (bool, S2Polyline) DecodeCompressed(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte)) return (false, default);
        var snap_level = decoder.Get8();
        if (snap_level > S2.kMaxCellLevel) return (false, default);

        if (!decoder.TryGetVarInt32(out var num_vertices)) return (false, default);
        if (num_vertices == 0)
        {
            // Empty polylines are allowed.
            S2Polyline shape = new(Array.Empty<S2Point>());
            return (true, shape);
        }

        // TODO(b/209937354): Prevent large allocations like in DecodeUncompressed.
        // This is more complicated due to the compressed encoding, but perhaps the
        // minimum required size can be bounded.
        var points =new S2Point[num_vertices];
        if (!S2DecodePointsCompressed(decoder, snap_level, points, 0))
        {
            return (false, default);
        }

        {
            S2Polyline shape = new(points);
            return (true, shape);
        }
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2Polyline(this);

    #endregion

    #region IEquatable

    public bool Equals(S2Polyline b)
    {
        return Vertices.SequenceEqual(b.Vertices);
    }
    public override int GetHashCode() => LinqUtils.GetSequenceHashCode(Vertices);

    #endregion

    // Wrapper class for indexing a polyline (see S2ShapeIndex).  Once this
    // object is inserted into an S2ShapeIndex it is owned by that index, and
    // will be automatically deleted when no longer needed by the index.  Note
    // that this class does not take ownership of the polyline itself (see
    // OwningShape below).  You can also subtype this class to store additional
    // data (see S2Shape for details).
    public class Shape : S2Shape
    {
        #region Fields, Constants

        // Define as enum so we don't have to declare storage.
        // TODO(user, b/210097200): Use static constexpr when C++17 is
        // allowed in opensource.
        public const TypeTag kTypeTag = TypeTag.S2Polyline;

        public S2Polyline Polyline { get; private set; }

        #endregion

        #region Constructors

        /// <summary>
        /// Must call Init().
        /// </summary>
        public Shape() { }

        // Initialization.  Does not take ownership of "polyline".
        //
        // Note that a polyline with one vertex is defined to have no edges.  Use
        // S2LaxPolylineShape or S2LaxClosedPolylineShape if you want to define a
        // polyline consisting of a single degenerate edge.
        public Shape(S2Polyline polyline) { Init(polyline); }

        #endregion

        #region Shape

        public virtual void Init(S2Polyline polyline)
        {
            // if (polyline.Vertices.Length == 1) => "S2Polyline.Shape with one vertex has no edges";
            Polyline = polyline;
        }

        #endregion

        #region IEncoder

        // Encodes the polyline using S2Polyline.Encode().
        public override void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
        {
            if (hint == CodingHint.FAST)
            {
                Polyline.EncodeUncompressed(encoder);
            }
            else
            {
                Polyline.EncodeMostCompact(encoder);
            }
        }
        // Decoding is defined only for S2Polyline::OwningShape below.

        #endregion

        #region S2Shape

        // S2Shape interface:

        public sealed override int NumEdges() => NumEdgesStatic(Polyline);

        private static int NumEdgesStatic(S2Polyline pol) =>
            Math.Max(0, pol.Vertices.Length - 1);

        public sealed override Edge GetEdge(int e) =>
            new Edge(Polyline.Vertex(e), Polyline.Vertex(e + 1));

        public sealed override int Dimension() => 1;

        public sealed override ReferencePoint GetReferencePoint() =>
            ReferencePoint.FromContained(false);

        public sealed override int NumChains() =>
            Math.Min(1, NumEdgesStatic(Polyline));  // Avoid virtual call.

        public sealed override Chain GetChain(int i)
        {
            MyDebug.Assert(i == 0);
            return new Chain(0, NumEdgesStatic(Polyline));  // Avoid virtual call.
        }

        public sealed override Edge ChainEdge(int i, int j)
        {
            MyDebug.Assert(i == 0);
            return new Edge(Polyline.Vertex(j), Polyline.Vertex(j + 1));
        }

        public sealed override ChainPosition GetChainPosition(int e) => new(0, e);

        public override TypeTag GetTypeTag() => kTypeTag;

        #endregion
    }

    // Like Shape, except that the S2Polyline is automatically deleted when this
    // object is deleted by the S2ShapeIndex.  This is useful when an S2Polyline
    // is constructed solely for the purpose of indexing it.
    public class OwningShape : Shape, IInitEncoder<OwningShape>
    {
        private S2Polyline owned_polyline_;

        public OwningShape() { }  // Must call Init().
        public OwningShape(S2Polyline polyline) : base(polyline) => owned_polyline_ = polyline;

        public override void Init(S2Polyline polyline)
        {
            base.Init(polyline);
            owned_polyline_ = polyline;
        }

        public static (bool, OwningShape?) Init(Decoder decoder)
        {
            var (success, result) = Decode(decoder);
            if (!success) return (false, null);

            return (true, new OwningShape(result));
        }
    }
}

