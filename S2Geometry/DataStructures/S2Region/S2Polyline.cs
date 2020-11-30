using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

namespace S2Geometry
{
    // An S2Polyline represents a sequence of zero or more vertices connected by
    // straight edges (geodesics).  Edges of length 0 and 180 degrees are not
    // allowed, i.e. adjacent vertices should not be identical or antipodal.
    public sealed class S2Polyline : S2Region<S2Polyline>, IEquatable<S2Polyline>, ICoder
    {
        #region Fields, Constants

        private S2Point[] vertices_; 

        #endregion

        #region Constructors

        // Creates an empty S2Polyline that should be initialized by calling Init()
        // or Decode().
        public S2Polyline() { }

        /// <summary>
        /// must call Init
        /// </summary>
        /// <param name="other"></param>
        public S2Polyline(S2Polyline other)
        {
            vertices_ = other.vertices_;
            other.vertices_ = null;
        }

        // Convenience constructor that call Init() with the given vertices.
        public S2Polyline(S2Point[] vertices, bool assertIsValid = true)
        {
            Init(vertices, assertIsValid);
        }

        // Convenience constructor that call Init() with the given vertices.
        public S2Polyline(S2LatLng[] vertices)
        {
            Init(vertices);
        }

        #endregion

        #region S2Polyline

        // Initialize a polyline that connects the given vertices. Empty polylines are
        // allowed.  Adjacent vertices should not be identical or antipodal.  All
        // vertices should be unit length.
        public void Init(S2Point[] vertices, bool assertIsValid = true)
        {
            vertices_ = vertices;
#if s2debug
            if (assertIsValid)
            {
                // Note that s2debug is false in optimized builds (by default).
                Assert.True(IsValid);
            }
#endif
        }

        // Convenience initialization function that accepts latitude-longitude
        // coordinates rather than S2Points.
        public void Init(S2LatLng[] vertices)
        {
            Init(vertices.Select(t => t.ToPoint()).ToArray());
        }

        // Return true if the given vertices form a valid polyline.
        public bool IsValid
        {
            get
            {
                if (FindValidationError(out S2Error error))
                {
#if s2debug
                    System.Diagnostics.Debug.WriteLine($"IsValid: {error}");
#endif
                    return false;
                }
                return true;
            }
        }

        // Returns true if this is *not* a valid polyline and sets "error"
        // appropriately.  Otherwise returns false and leaves "error" unchanged.
        //
        // REQUIRES: error != null
        public bool FindValidationError(out S2Error error)
        {
            // All vertices must be unit length.
            for (int i = 0; i < vertices_.Length; ++i)
            {
                if (!Vertex(i).IsUnitLength)
                {
                    error = new(S2ErrorCode.NOT_UNIT_LENGTH, $"Vertex {i} is not unit length");
                    return true;
                }
            }
            // Adjacent vertices must not be identical or antipodal.
            for (int i = 1; i < vertices_.Length; ++i)
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

        public int NumVertices { get { return vertices_?.Length ?? 0; } }
        public S2Point[] VerticesClone { get { return (S2Point[])vertices_.Clone(); } }

        public S2Point Vertex(int k)
        {
            Assert.True(k >= 0);
            Assert.True(k <= vertices_.Length);
            return vertices_[k];
        }

        // Return the length of the polyline.
        public S1Angle GetLength()
        {
            return S2PolylineMeasures.GetLength(vertices_.Take(vertices_.Length).ToArray());
        }

        // Return the true centroid of the polyline multiplied by the length of the
        // polyline (see s2centroids.h for details on centroids).  The result is not
        // unit length, so you may want to normalize it.
        //
        // Prescaling by the polyline length makes it easy to compute the centroid
        // of several polylines (by simply adding up their centroids).
        public S2Point GetCentroid()
        {
            return S2PolylineMeasures.GetCentroid(vertices_.Take(vertices_.Length).ToArray());
        }

        // Return the point whose distance from vertex 0 along the polyline is the
        // given fraction of the polyline's total length.  Fractions less than zero
        // or greater than one are clamped.  The return value is unit length.  This
        // cost of this function is currently linear in the number of vertices.
        // The polyline must not be empty.
        public S2Point Interpolate(double fraction)
        {
            return GetSuffix(fraction, out _);
        }

        // Like Interpolate(), but also return the index of the next polyline
        // vertex after the interpolated point P.  This allows the caller to easily
        // construct a given suffix of the polyline by concatenating P with the
        // polyline vertices starting at "next_vertex".  Note that P is guaranteed
        // to be different than vertex(*next_vertex), so this will never result in
        // a duplicate vertex.
        //
        // The polyline must not be empty.  Note that if "fraction" >= 1.0, then
        // "next_vertex" will be set to vertices_.Length (indicating that no vertices
        // from the polyline need to be appended).  The value of "next_vertex" is
        // always between 1 and vertices_.Length.
        //
        // This method can also be used to construct a prefix of the polyline, by
        // taking the polyline vertices up to "next_vertex - 1" and appending the
        // returned point P if it is different from the last vertex (since in this
        // case there is no guarantee of distinctness).
        public S2Point GetSuffix(double fraction, out int next_vertex)
        {
            Assert.True(vertices_.Length > 0);
            // We intentionally let the (fraction >= 1) case fall through, since
            // we need to handle it in the loop below in any case because of
            // possible roundoff errors.
            if (fraction <= 0)
            {
                next_vertex = 1;
                return Vertex(0);
            }
            var length_sum = S1Angle.Zero;
            for (int i = 1; i < vertices_.Length; ++i)
            {
                length_sum += new S1Angle(Vertex(i - 1), Vertex(i));
            }
            S1Angle target = fraction * length_sum;
            for (int i = 1; i < vertices_.Length; ++i)
            {
                var length = new S1Angle(Vertex(i - 1), Vertex(i));
                if (target < length)
                {
                    // This interpolates with respect to arc length rather than
                    // straight-line distance, and produces a unit-length result.
                    S2Point result = S2EdgeDistances.InterpolateAtDistance(target, Vertex(i - 1),
                                                                       Vertex(i));
                    // It is possible that (result == vertex(i)) due to rounding errors.
                    next_vertex = (result == Vertex(i)) ? (i + 1) : i;
                    return result;
                }
                target -= length;
            }
            next_vertex = vertices_.Length;
            return Vertex(vertices_.Length - 1);
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
            Assert.True(vertices_.Length > 0);
            if (vertices_.Length < 2)
            {
                return 0;
            }
            var length_sum = S1Angle.Zero;
            for (int i = 1; i < next_vertex; ++i)
            {
                length_sum += new S1Angle(Vertex(i - 1), Vertex(i));
            }
            S1Angle length_to_point = length_sum + new S1Angle(Vertex(next_vertex - 1), point);
            for (int i = next_vertex; i < vertices_.Length; ++i)
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
            Assert.True(vertices_.Length > 0);

            if (vertices_.Length == 1)
            {
                // If there is only one vertex, it is always closest to any given point.
                next_vertex = 1;
                return Vertex(0);
            }

            // Initial value larger than any possible distance on the unit sphere.
            S1Angle min_distance = S1Angle.FromRadians(10);
            int min_index = -1;

            // Find the line segment in the polyline that is closest to the point given.
            for (int i = 1; i < vertices_.Length; ++i)
            {
                var distance_to_segment = S2EdgeDistances.GetDistance(point, Vertex(i - 1), Vertex(i));
                if (distance_to_segment < min_distance)
                {
                    min_distance = distance_to_segment;
                    min_index = i;
                }
            }
            Assert.True(min_index != -1);

            // Compute the point on the segment found that is closest to the point given.
            var closest_point = S2EdgeDistances.Project(point, Vertex(min_index - 1), Vertex(min_index));

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
            Assert.True(vertices_.Length >= 2);

            S2Point closest_point = Project(point, out int next_vertex);

            Assert.True(next_vertex >= 1);
            Assert.True(next_vertex <= vertices_.Length);

            // If the closest point C is an interior vertex of the polyline, let B and D
            // be the previous and next vertices.  The given point P is on the right of
            // the polyline (locally) if B, P, D are ordered CCW around vertex C.
            if (closest_point == Vertex(next_vertex - 1) && next_vertex > 1 && next_vertex < vertices_.Length)
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
            if (next_vertex == vertices_.Length)
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
            if (vertices_.Length <= 0 || line.vertices_.Length <= 0)
            {
                return false;
            }

            if (!GetRectBound().Intersects(line.GetRectBound()))
            {
                return false;
            }

            // TODO(ericv): Use S2ShapeIndex here.
            for (int i = 1; i < vertices_.Length; ++i)
            {
                var crosser = new S2EdgeCrosser(Vertex(i - 1), Vertex(i), line.Vertex(0));
                for (int j = 1; j < line.vertices_.Length; ++j)
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
        public void SubsampleVertices(S1Angle tolerance, out int[] indices)
        {
            indices = null;
            var indicesTmp = new List<int>();
            if (vertices_.Length == 0) return;

            indicesTmp.Add(0);
            var clamped_tolerance = S1Angle.Max(tolerance, S1Angle.FromRadians(0));
            for (int index = 0; index + 1 < vertices_.Length;)
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
            Assert.True(tolerance.Radians >= 0);
            Assert.True((index + 1) < polyline.vertices_.Length);

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
            var frame = S2PointUtil.GetFrame(origin);

            // As we go along, we keep track of the current wedge of angles and the
            // distance to the last vertex (which must be non-decreasing).
            S1Interval current_wedge = S1Interval.Full;
            double last_distance = 0;

            for (++index; index < polyline.vertices_.Length; index++)
            {
                var candidate = polyline.Vertex(index);
                double distance = origin.Angle(candidate);

                // We don't allow simplification to create edges longer than 90 degrees,
                // to avoid numeric instability as lengths approach 180 degrees.  (We do
                // need to allow for original edges longer than 90 degrees, though.)
                if (distance > S2Constants.M_PI_2 && last_distance > 0) break;

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
                S2Point direction = S2PointUtil.ToFrame(frame, candidate);
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
                Assert.True(!current_wedge.IsEmpty);
            }
            // We break out of the loop when we reach a vertex index that can't be
            // included in the line segment, so back up by one vertex.
            return index - 1;
        }

        // Return true if two polylines have the same number of vertices, and
        // corresponding vertex pairs are separated by no more than "max_error".
        // (For testing purposes.)
        public bool ApproxEquals(S2Polyline b) => ApproxEquals(b, S1Angle.FromRadians(S2Constants.DoubleError));

        public bool ApproxEquals(S2Polyline b, S1Angle max_error)
        {
            if (vertices_.Length != b.vertices_.Length) return false;
            for (int offset = 0; offset < vertices_.Length; ++offset)
            {
                if (!S2PointUtil.ApproxEquals(Vertex(offset), b.Vertex(offset), max_error))
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
            // number of states examined is O(n+m) where n = this.vertices_.Length and m =
            // covered.vertices_.Length.  Without it, the amount of work could be as high as
            // O((n*m)^2).  Using set, the running time is O((n*m) log (n*m)).
            //
            // TODO(user): Benchmark this, and see if the set is worth it.

            if (covered.vertices_.Length == 0) return true;
            if (this.vertices_.Length == 0) return false;

            var pending = new List<SearchState>();
            var done = new SortedSet<SearchState>(new SearchStateComparer());

            // Find all possible starting states.
            for (int i = 0, next_i = NextDistinctVertex(this, 0), next_next_i;
                 next_i < this.vertices_.Length; i = next_i, next_i = next_next_i)
            {
                next_next_i = NextDistinctVertex(this, next_i);
                S2Point closest_point = S2EdgeDistances.Project(
                    covered.Vertex(0), Vertex(i), Vertex(next_i));

                // In order to avoid duplicate starting states, we exclude the end vertex
                // of each edge *except* for the last non-degenerate edge.
                if ((next_next_i == this.vertices_.Length ||
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

                int next_i = NextDistinctVertex(this, state.I);
                int next_j = NextDistinctVertex(covered, state.J);
                if (next_j == covered.vertices_.Length)
                {
                    return true;
                }
                else if (next_i == this.vertices_.Length)
                {
                    continue;
                }

                S2Point i_begin, j_begin;
                if (state.IInProgress)
                {
                    j_begin = covered.Vertex(state.J);
                    i_begin = S2EdgeDistances.Project(
                        j_begin, Vertex(state.I), Vertex(next_i));
                }
                else
                {
                    i_begin = Vertex(state.I);
                    j_begin = S2EdgeDistances.Project(
                        i_begin, covered.Vertex(state.J), covered.Vertex(next_j));
                }

                if (S2EdgeDistances.IsEdgeBNearEdgeA(j_begin, covered.Vertex(next_j),
                                         i_begin, Vertex(next_i), max_error))
                {
                    pending.Add(new SearchState(next_i, state.J, false));
                }
                if (S2EdgeDistances.IsEdgeBNearEdgeA(i_begin, Vertex(next_i),
                                         j_begin, covered.Vertex(next_j), max_error))
                {
                    pending.Add(new SearchState(state.I, next_j, true));
                }
            }
            return false;
        }

        // Returns the total number of bytes used by the polyline.
        public int SpaceUsed()
        {
            return Marshal.SizeOf(typeof(S2Polyline)) + vertices_.Length * Marshal.SizeOf(typeof(S2Point));
        }

        // Return the first i > "index" such that the ith vertex of "pline" is not at
        // the same point as the "index"th vertex.  Returns pline.num_vertices() if
        // there is no such value of i.
        private static int NextDistinctVertex(S2Polyline pline, int index)
        {
            S2Point initial = pline.Vertex(index);
            do
            {
                ++index;
            } while (index < pline.vertices_.Length && pline.Vertex(index) == initial);
            return index;
        } 

        public void Reverse()
        {
            Array.Reverse(vertices_);
        }

        #endregion

        #region S2Region

        ////////////////////////////////////////////////////////////////////////
        // S2Region interface (see s2region.h for details):

        public override object Clone()
        {
            return new S2Polyline { vertices_ = (S2Point[])this.vertices_.Clone() };
        }
        public override S2Cap GetCapBound() => GetRectBound().GetCapBound();
        public override S2LatLngRect GetRectBound()
        {
            var bounder = new S2LatLngRectBounder();
            for (int i = 0; i < vertices_.Length; ++i)
            {
                bounder.AddPoint(Vertex(i));
            }
            return bounder.GetBound();
        }
        public override bool Contains(S2Cell cell) => false;
        public override bool MayIntersect(S2Cell cell)
        {
            if (vertices_.Length == 0) return false;

            // We only need to check whether the cell contains vertex 0 for correctness,
            // but these tests are cheap compared to edge crossings so we might as well
            // check all the vertices.
            for (int i = 0; i < vertices_.Length; ++i)
            {
                if (cell.Contains(Vertex(i))) return true;
            }
            var cell_vertices = new S2Point[4];
            for (int i = 0; i < 4; ++i)
            {
                cell_vertices[i] = cell.GetVertex(i);
            }
            for (int j = 0; j < 4; ++j)
            {
                var crosser = new S2EdgeCrosser(cell_vertices[j], cell_vertices[(j + 1) & 3], Vertex(0));
                for (int i = 1; i < vertices_.Length; ++i)
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
        public override bool Contains(S2Point p) => false;

        #endregion

        #region ICoder

        // Appends a serialized representation of the S2Polyline to "encoder".
        //
        // REQUIRES: "encoder" uses the defaultructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public void Encode(Encoder encoder)
        {
            encoder.Ensure(vertices_.Length * Marshal.SizeOf(typeof(S2Point)) + 10);  // sufficient

            encoder.Put8(S2Constants.kCurrentLosslessEncodingVersionNumber);
            encoder.Put32(vertices_.Length);
            encoder.PutPoints(vertices_);

            Assert.True(encoder.Avail() >= 0);
        }

        // Decodes an S2Polyline encoded with Encode().  Returns true on success.
        public static (bool success, S2Polyline result) DecodeStatic(Decoder decoder)
        {
            if (decoder.Avail() < sizeof(byte) + sizeof(UInt32)) return (false, null);
            byte version = decoder.Get8();
            if (version > S2Constants.kCurrentLosslessEncodingVersionNumber) return (false, null);

            var count = (int)decoder.Get32();
            if (decoder.Avail() < count) return (false, null);
            var verts = new S2Point[count];
            decoder.GetPoints(verts, 0, count);

            var pol = new S2Polyline(verts);
#if s2debug
            Assert.True(pol.IsValid);
#endif
            return (true, pol);
        }

        #endregion

        #region IEquatable

        // Return true if two polylines are exactly the same.
        public override bool Equals(object obj)
        {
            return obj is S2Polyline poly && Equals(poly);
        }
        public override bool Equals(S2Polyline b)
        {
            return vertices_.SequenceEqual(b.vertices_);
        }
        public override int GetHashCode()
        {
            return LinqUtils.GetSequenceHashCode(vertices_);
        }
        public static bool operator ==(S2Polyline left, S2Polyline right)
        {
            return Equals(left, right);
        }
        public static bool operator !=(S2Polyline left, S2Polyline right)
        {
            return !Equals(left, right);
        }

        #endregion

        // This struct represents a search state in the NearlyCovers algorithm
        // below.  See the description of the algorithm for details.
        public readonly struct SearchState
        {
            public readonly int I;
            public readonly int J;
            public readonly bool IInProgress;

            public SearchState(int i_val, int j_val, bool i_in_progress_val)
            {
                I = i_val;
                J = j_val;
                IInProgress = i_in_progress_val;
            }
        }

        // Needed for storing SearchStates in a SortedSet.  The ordering
        // chosen has no special meaning.
        public class SearchStateComparer : IComparer<SearchState>
        {
            public int Compare(SearchState a, SearchState b)
            {
                if (a.I.CompareTo(b.I) != 0)
                {
                    return a.I.CompareTo(b.I);
                }
                else if (a.J.CompareTo(b.J) != 0)
                {
                    return a.J.CompareTo(b.J);
                }
                else if (a.IInProgress.CompareTo(b.IInProgress) != 0)
                {
                    return a.IInProgress.CompareTo(b.IInProgress);
                }
                else
                {
                    return 0;
                }
            }
        }

        // Wrapper class for indexing a polyline (see S2ShapeIndex).  Once this
        // object is inserted into an S2ShapeIndex it is owned by that index, and
        // will be automatically deleted when no longer needed by the index.  Note
        // that this class does not take ownership of the polyline itself (see
        // OwningShape below).  You can also subtype this class to store additional
        // data (see S2Shape for details).
        public class Shape : S2Shape
        {
            #region Fields, Constants

            public S2Polyline Polyline { get; private set; }

            #endregion

            #region Constructors

            /// <summary>
            /// Must call Init().
            /// </summary>
            public Shape() { Polyline = null; }

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
                // if (polyline.vertices_.Length == 1) => "S2Polyline.Shape with one vertex has no edges";
                Polyline = polyline;
            }

            // Encodes the polyline using S2Polyline.Encode().
            public void Encode(Encoder encoder)
            {
                // TODO(geometry-library): Support compressed encodings.
                Polyline.Encode(encoder);
            }

            #endregion

            #region S2Shape

            // Decoding is defined only for S2Polyline.OwningShape below.

            // S2Shape interface:

            public sealed override int NumEdges => NumEdgesStatic(Polyline);
            private static int NumEdgesStatic(S2Polyline pol)
            {
                return Math.Max(0, pol.vertices_.Length - 1);
            }
            public sealed override Edge GetEdge(int e)
            {
                return new Edge(Polyline.Vertex(e), Polyline.Vertex(e + 1));
            }
            public sealed override int Dimension() { return 1; }
            public sealed override ReferencePoint GetReferencePoint()
            {
                return ReferencePoint.FromContained(false);
            }
            public sealed override int NumChains()
            {
                return Math.Min(1, NumEdgesStatic(Polyline));  // Avoid virtual call.
            }
            public sealed override Chain GetChain(int i)
            {
                Assert.True(i == 0);
                return new Chain(0, NumEdgesStatic(Polyline));  // Avoid virtual call.
            }
            public sealed override Edge ChainEdge(int i, int j)
            {
                Assert.True(i == 0);
                return new Edge(Polyline.Vertex(j), Polyline.Vertex(j + 1));
            }
            public sealed override ChainPosition GetChainPosition(int e)
            {
                return new ChainPosition(0, e);
            }
            public override TypeTag GetTypeTag() => TypeTag.S2Polyline; 

            #endregion
        }

        // Like Shape, except that the S2Polyline is automatically deleted when this
        // object is deleted by the S2ShapeIndex.  This is useful when an S2Polyline
        // is constructed solely for the purpose of indexing it.
        public class OwningShape : Shape, IEncodeInit
        {
            private S2Polyline owned_polyline_;

            public OwningShape() { }  // Must call Init().
            public OwningShape(S2Polyline polyline) : base(polyline) => owned_polyline_ = polyline;

            public override void Init(S2Polyline polyline)
            {
                base.Init(polyline);
                owned_polyline_ = polyline;
            }

            public bool Init(Decoder decoder)
            {
                var (success, result) = DecodeStatic(decoder);
                owned_polyline_ = result;
                if (!success) return false;
                base.Init(owned_polyline_);
                return true;
            }
        }
    }
}

