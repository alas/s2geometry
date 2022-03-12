// This class provides a way to expand an arbitrary collection of geometry by
// a fixed radius (an operation variously known as "buffering", "offsetting",
// or "Minkowski sum with a disc").  The output consists of a polygon
// (possibly with multiple shells) that contains all points within the given
// radius of the original geometry.
//
// The radius can also be negative, in which case the geometry is contracted.
// This causes the boundaries of polygons to shrink or disappear, and removes
// all points and polylines.
//
// The input consists of a sequence of layers.  Each layer may consist of any
// combination of points, polylines, and polygons, with the restriction that
// polygon interiors within each layer may not intersect any other geometry
// (including other polygon interiors).  The output is the union of the
// buffered input layers.  Note that only a single layer is allowed if the
// buffer radius is negative.
//
// This class may be used to compute polygon unions by setting the buffer
// radius to zero.  The union is computed using a single snapping operation.
//
// Note that if you only want to compute an S2CellId covering of the buffered
// geometry, it is much faster to use S2ShapeIndexBufferedRegion instead.
//
// Keywords: buffer, buffering, expand, expanding, offset, offsetting,
//           widen, contract, shrink, Minkowski sum


// The algorithm below essentially computes the offset curve of the original
// boundary, and uses this curve to divide the sphere into regions of constant
// winding number.  Since winding numbers on the sphere are relative rather
// than absolute (see s2winding_operation.h), we also need to keep track of
// the desired winding number at a fixed reference point.  The initial winding
// number for this point is the number of input shapes that contain it.  We
// then update it during the buffering process by imagining a "sweep edge"
// that extends from the current point A on the input boundary to the
// corresponding point B on the offset curve.  As we process an input loop and
// generate the corresponding offset curve, the sweep edge moves continuously
// and covers the entire buffer region (i.e., the region added to or
// subtracted from the input geometry).  We increase the winding number of the
// reference point by one whenever it crosses the sweep edge from left to
// right, and we decrease the winding number by one whenever it crosses the
// sweep edge from right to left.
//
// Concave vertices require special handling, because the corresponding offset
// curve can leave behind regions whose winding number is zero or negative.
// We handle this by splicing the concave vertex into the offset curve itself;
// this effectively terminates the current buffer region and starts a new one,
// such that the region of overlap is counted twice (i.e., its winding number
// increases by two).  The end result is the same as though we had computed
// the union of a sequence of buffered convex boundary segments.  This trick
// is described in the following paper: "Polygon Offsetting by Computing
// Winding Numbers" (Chen and McMains, Proceedings of IDETC/CIE 2005).
//
// TODO(ericv): The algorithm below is much faster than, say, computing the
// union of many buffered edges.  However further improvements are possible.
// In particular, there is an unimplemented optimization that would make it
// much faster to buffer concave boundaries when the buffer radius is large.


// The errors due to buffering can be categorized as follows:
//
//  1. Requested error.  This represents the error due to approximating the
//     buffered boundary as a sequence of line segments rather than a sequence
//     of circular arcs.  It is largely controlled by options.error_fraction(),
//     and can be bounded as
//
//       max(kMinRequestedError, error_fraction * buffer_radius)
//
//     where kMinRequestedError reflects the fact that S2Points do not have
//     infinite precision.  (For example, it makes no sense to attempt to
//     buffer geometry by 1e-100 radians because the spacing between
//     representable S2Points is only about 2e-16 radians in general.)
//
//  2. Relative interpolation errors.  These are numerical errors whose
//     magnitude is proportional to the buffer radius.  For such errors the
//     worst-case coefficient of proportionality turns out to be so tiny
//     compared to the smallest allowable error fraction (kMinErrorFraction)
//     that we simply ignore such errors.
//
//  3. Absolute interpolation errors.  These are numerical errors that are not
//     proportional to the buffer radius.  The main sources of such errors are
//     (1) calling S2.RobustCrossProd() to compute edge normals, and (2) calls
//     to S2.GetPointOnRay() to interpolate points along the buffered
//     boundary.  It is possible to show that this error is at most
//     kMaxAbsoluteInterpolationError as defined below.
//
// Putting these together, the final error bound looks like this:
//
//   max_error = kMaxAbsoluteInterpolationError +
//               max(kMinRequestedError,
//                   max(kMinErrorFraction, options.error_fraction()) *
//                   options.buffer_radius())

namespace S2Geometry;

public class S2BufferOperation
{
    // The maximum angular spacing between representable S2Points on the unit
    // sphere is roughly 2 * DBL_ERR.  We require the requested absolute error to
    // be at least this large because attempting to achieve a smaller error does
    // not increase the precision of the result and can increase the running time
    // and space requirements considerably.
    private static readonly S1Angle kMinRequestedError = S1Angle.FromRadians(2 * S2.DoubleError);

    // The maximum absolute error due to interpolating points on the buffered
    // boundary.  The following constant bounds the maximum additional error
    // perpendicular to the buffered boundary due to all steps of the calculation
    // (S2.RobustCrossProd, the two calls to GetPointOnRay, etc).
    //
    // This distance represents about 10 nanometers on the Earth's surface.  Note
    // that this is a conservative upper bound and that it is difficult to
    // construct inputs where the error is anywhere close to this large.
    private const double kMaxAbsoluteInterpolationError =
        S2.kGetPointOnLineError + S2.kGetPointOnRayPerpendicularError;

    private static readonly S1Angle kMaxAbsoluteInterpolationErrorS1Angle = S1Angle.FromRadians(kMaxAbsoluteInterpolationError);

    // Default constructor; requires Init() to be called.
    public S2BufferOperation() { }

    // Convenience constructor that calls Init().
    public S2BufferOperation(S2Builder.Layer result_layer, Options? options = null)
        => Init(result_layer, options ?? new());

    // Starts a buffer operation that sends the output polygon to the given
    // S2Builder layer.  This method may be called more than once.
    //
    // Note that buffering always yields a polygon, even if the input includes
    // polylines and points.  If the buffer radius is zero, points and polylines
    // will be converted into degenerate polygon loops; if the buffer radius is
    // negative, points and polylines will be removed.
    public void Init(S2Builder.Layer result_layer,
        Options? options = null)
    {
        Options_ = options ?? new();
        ref_point_ = S2.Origin;
        ref_winding_ = 0;
        have_input_start_ = false;
        have_offset_start_ = false;
        buffer_sign_ = Options_.buffer_radius_.Radians.Sgn();
        S1Angle abs_radius = S1Angle.Abs(Options_.buffer_radius_);
        S1Angle requested_error = S1Angle.Max(kMinRequestedError,
            Options_.error_fraction_ * abs_radius);
        S1Angle max_error = kMaxAbsoluteInterpolationErrorS1Angle + requested_error;
        if (abs_radius <= max_error)
        {
            // If the requested radius is smaller than the maximum error, buffering
            // could yield points on the wrong side of the original input boundary
            // (e.g., shrinking geometry slightly rather than expanding it).  Rather
            // than taking that risk, we set the buffer radius to zero when this
            // happens (which causes the original geometry to be returned).
            abs_radius_ = S1ChordAngle.Zero;
            buffer_sign_ = 0;
        }
        else if (abs_radius + max_error >= S1Angle.FromRadians(S2.M_PI))
        {
            // If the permissible range of buffer angles includes Pi then we might
            // as well take advantage of that.
            abs_radius_ = S1ChordAngle.Straight;
        }
        else
        {
            abs_radius_ = new S1ChordAngle(abs_radius);
            S1Angle vertex_step = GetMaxEdgeSpan(abs_radius, requested_error);
            vertex_step_ = new S1ChordAngle(vertex_step);

            // We take extra care to ensure that points are buffered as regular
            // polygons.  The step angle is adjusted up slightly to ensure that we
            // don't wind up with a tiny extra edge.
            point_step_ = S1ChordAngle.FromRadians(
                2 * S2.M_PI / Math.Ceiling(2 * S2.M_PI / vertex_step.Radians) + 1e-15);

            // Edges are buffered only if the buffer radius (including permissible
            // error) is less than 90 degrees.
            S1Angle edge_radius = S1Angle.FromRadians(S2.M_PI_2) - abs_radius;
            if (edge_radius > max_error)
            {
                edge_step_ = new S1ChordAngle(GetMaxEdgeSpan(edge_radius, requested_error));
            }
        }

        // The buffered output should include degeneracies (i.e., isolated points
        // and/or sibling edge pairs) only if (1) the user specified a non-negative
        // buffer radius, and (2) the adjusted buffer radius is zero.  The only
        // purpose of keeping degeneracies is to allow points/polylines in the input
        // geometry to be converted back to points/polylines in the output if the
        // client so desires.
        S2WindingOperation.Options winding_options = new(options.snap_function_);
        winding_options.include_degeneracies_ = (
            buffer_sign_ == 0 && Options_.buffer_radius_ >= S1Angle.Zero);
        winding_options.memory_tracker_ = options.memory_tracker_;
        op_.Init(result_layer, winding_options);
        tracker_.Init(options.memory_tracker_);
    }

    public Options Options_ { get; private set; }

    // Each call below represents a different input layer.  Note that if the
    // buffer radius is negative, then at most one input layer is allowed
    // (ignoring any layers that contain only points and polylines).

    // Adds an input layer containing a single point.
    public void AddPoint(S2Point point)
    {
        // If buffer_radius < 0, points are discarded.
        if (buffer_sign_ < 0) return;

        // Buffering by 180 degrees or more always yields the full polygon.
        // (We don't need to worry about buffering by 180 degrees yielding
        // a degenerate hole because error_fraction_ is always positive.
        if (abs_radius_ >= S1ChordAngle.Straight)
        {
            AddFullPolygon();
            return;
        }

        // If buffer_radius == 0, points are converted into degenerate loops.
        if (buffer_sign_ == 0)
        {
            if (!tracker_.AddSpace(path_, 1)) return;
            path_.Add(point);
        }
        else
        {
            // Since S1ChordAngle can only represent angles between 0 and 180 degrees,
            // we generate the circle in four 90 degree increments.
            SetInputVertex(point);
            S2Point start = S2.Ortho(point);
            S1ChordAngle angle = S1ChordAngle.Zero;
            for (int quadrant = 0; quadrant < 4; ++quadrant)
            {
                // Generate 90 degrees of the circular arc.  Normalize "rotate_dir" at
                // each iteration to avoid magnifying normalization errors in "point".
                S2Point rotate_dir = point.CrossProd(start).Normalize();
                for (; angle < S1ChordAngle.Right; angle += point_step_)
                {
                    S2Point dir = S2.GetPointOnRay(start, rotate_dir, angle);
                    AddOffsetVertex(S2.GetPointOnRay(point, dir, abs_radius_));
                }
                angle -= S1ChordAngle.Right;
                start = rotate_dir;
            }
            CloseBufferRegion();
        }
        OutputPath();
    }

    // Adds an input layer containing a polyline.  Note the following:
    //
    //  - Polylines with 0 or 1 vertices are considered to be empty.
    //  - A polyline with 2 identical vertices is equivalent to a point.
    //  - Polylines have end caps (see Options.end_cap_style).
    //  - One-sided polyline buffering is supported (see Options.polyline_side).
    public void AddPolyline(S2PointSpan polyline)
    {
        // Left-sided buffering is supported by reversing the polyline and then
        // buffering on the right.
        if (Options_.polyline_side_ == PolylineSide.LEFT)
        {
            polyline.Reverse();
        }

        // If buffer_radius < 0, polylines are discarded.
        if (buffer_sign_ < 0 || !tracker_.Ok()) return;

        // Polylines with 0 or 1 vertices are defined to have no edges.
        int n = polyline.Count;
        if (n <= 1) return;

        // Polylines with one degenerate edge are treated as points.
        if (n == 2 && polyline[0] == polyline[1])
        {
            AddPoint(polyline[0]);
            return;
        }

        // Buffering by 180 degrees or more always yields the full polygon.
        if (abs_radius_ >= S1ChordAngle.Straight)
        {
            AddFullPolygon();
            return;
        }

        // If buffer_radius == 0, polylines are converted into degenerate loops.
        if (buffer_sign_ == 0)
        {
            if (!tracker_.AddSpace(path_, 2 * (n - 1))) return;

            path_ = polyline.Take(polyline.Count - 1).ToList();
            path_.AddRange(polyline.Skip(1).Reverse());
        }
        else
        {
            // Otherwise we buffer each side of the polyline separately.
            SetInputVertex(polyline[0]);
            AddStartCap(polyline[0], polyline[1]);
            for (int i = 0; i < n - 2; ++i)
            {
                BufferEdgeAndVertex(polyline[i], polyline[i + 1], polyline[i + 2]);
            }
            AddEdgeArc(polyline[n - 2], polyline[n - 1]);
            AddEndCap(polyline[n - 2], polyline[n - 1]);

            if (Options_.polyline_side_ == PolylineSide.BOTH)
            {
                for (int i = n - 3; i >= 0; --i)
                {
                    BufferEdgeAndVertex(polyline[i + 2], polyline[i + 1], polyline[i]);
                }
                AddEdgeArc(polyline[1], polyline[0]);
                CloseBufferRegion();
            }
            else
            {
                // The other side of the polyline is not buffered.  Note that for
                // PolylineSide.LEFT, the polyline direction has been reversed.
                if (!tracker_.AddSpace(path_, n)) return;
                path_.AddRange(polyline.AsEnumerable().Reverse());
                // Don't call CloseBufferRegion() since the path has already been closed.
            }
        }
        OutputPath();
    }

    // Adds an input layer containing a loop.  Note the following:
    //
    //  - A loop with no vertices is empty.
    //  - A loop with 1 vertex is equivalent to a point.
    //  - The interior of the loop is on its left.
    //  - Buffering a self-intersecting loop produces undefined results.
    public void AddLoop(S2PointLoopSpan loop)
    {
        if (!loop.Any()) return;
        BufferLoop(loop);

        // The vertex copying below could be avoided by adding a version of
        // S2LaxLoopShape that doesn't own its vertices.
        if (!tracker_.Ok()) return;
        ref_winding_ += S2ShapeUtil.ContainsBruteForce(new S2LaxLoopShape(loop.ToArray()),
            ref_point_) ? 1 : 0;
        num_polygon_layers_ += 1;
    }

    // Adds an input layer containing the given shape.  Shapes are handled as
    // points, polylines, or polygons according to the rules above.  In addition
    // note the following:
    //
    //  - Polygons holes may be degenerate (e.g., consisting of a
    //    single vertex or entirely of sibling pairs such as ABCBDB).
    //  - Full polygons are supported.  Note that since full polygons do
    //    not have a boundary, they are not affected by buffering.
    public void AddShape(S2Shape shape)
    {
        BufferShape(shape);
        ref_winding_ += S2ShapeUtil.ContainsBruteForce(shape, ref_point_) ? 1 : 0;
        num_polygon_layers_ += (shape.Dimension() == 2) ? 1 : 0;
    }

    // Adds an input layer containing all of the shapes in the given index.
    //
    // REQUIRES: The interiors of polygons must be disjoint from all other
    //           indexed geometry, including other polygon interiors.
    //           (S2BooleanOperation also requires this.)
    public void AddShapeIndex(S2ShapeIndex index)
    {
        int max_dimension = -1;
        foreach (var shape in index)
        {
            if (shape == null) continue;

            max_dimension = Math.Max(max_dimension, shape.Dimension());
            BufferShape(shape);
        }
        ref_winding_ += index.MakeS2ContainsPointQuery().Contains(ref_point_) ? 1 : 0;
        num_polygon_layers_ += (max_dimension == 2) ? 1 : 0;
    }

    // Computes the union of the buffered input shapes and sends the output
    // polygon to the S2Builder layer specified in the constructor.  Returns
    // true on success and otherwise sets "error" appropriately.
    //
    // Note that if the buffer radius is negative, only a single input layer is
    // allowed (ignoring any layers that contain only points and polylines).
    public bool Build(out S2Error error)
    {
        if (buffer_sign_ < 0 && num_polygon_layers_ > 1)
        {
            error = new(S2ErrorCode.FAILED_PRECONDITION,
                        "Negative buffer radius requires at most one polygon layer");
            return false;
        }
        return op_.Build(ref_point_, ref_winding_,
            S2WindingOperation.WindingRule.POSITIVE, out error);
    }

    private S1Angle GetMaxEdgeSpan(S1Angle radius, S1Angle requested_error)
    {
        // If the allowable radius range spans Pi/2 then we can use edges as long as
        // we like, however we always use at least 3 edges to approximate a circle.
        S1Angle step = S1Angle.FromRadians(2 * S2.M_PI / 3 + 1e-15);
        S1Angle min_radius = radius - requested_error;
        System.Diagnostics.Debug.Assert(min_radius >= S1Angle.Zero);
        if (radius.Radians < S2.M_PI_2)
        {
            step = S1Angle.Min(step, S1Angle.FromRadians(2 * Math.Acos(min_radius.Tan() / radius.Tan())));
        }
        else if (min_radius.Radians > S2.M_PI_2)
        {
            step = S1Angle.Min(step, S1Angle.FromRadians(2 * Math.Acos(radius.Tan() / min_radius.Tan())));
        }
        return step;
    }

    // The sweep edge AB (see introduction) consists of one point on the input
    // boundary (A) and one point on the offset curve (B).  This function advances
    // the sweep edge by moving its first vertex A to "new_a" and updating the
    // winding number of the reference point if necessary.
    private void SetInputVertex(S2Point new_a)
    {
        if (have_input_start_)
        {
            System.Diagnostics.Debug.Assert(have_offset_start_);
            UpdateRefWinding(sweep_a_, sweep_b_, new_a);
        }
        else
        {
            input_start_ = new_a;
            have_input_start_ = true;
        }
        sweep_a_ = new_a;
    }

    // Adds the point "new_b" to the offset path.  Also advances the sweep edge AB
    // by moving its second vertex B to "new_b" and updating the winding number of
    // the reference point if necessary (see introduction).
    private void AddOffsetVertex(S2Point new_b)
    {
        if (!tracker_.AddSpace(path_, 1)) return;
        path_.Add(new_b);
        if (have_offset_start_)
        {
            System.Diagnostics.Debug.Assert(have_input_start_);
            UpdateRefWinding(sweep_a_, sweep_b_, new_b);
        }
        else
        {
            offset_start_ = new_b;
            have_offset_start_ = true;
        }
        sweep_b_ = new_b;
    }

    // Finishes buffering the current loop by advancing the sweep edge back to its
    // starting location, updating the winding number of the reference point if
    // necessary.
    private void CloseBufferRegion()
    {
        if (have_offset_start_ && have_input_start_)
        {
            UpdateRefWinding(sweep_a_, sweep_b_, input_start_);
            UpdateRefWinding(input_start_, sweep_b_, offset_start_);
        }
    }

    // Outputs the current buffered path (which is assumed to be a loop), and
    // resets the state to prepare for buffering a new loop.
    private void OutputPath()
    {
        op_.AddLoop(path_);
        path_.Clear();  // Does not change capacity.
        have_input_start_ = false;
        have_offset_start_ = false;
    }

    // Given a triangle ABC that has just been covered by the sweep edge AB,
    // updates the winding number of the reference point if necessary.
    private void UpdateRefWinding(S2Point a, S2Point b, S2Point c)
    {
        // TODO(ericv): This code could be made much faster by maintaining a
        // bounding plane that separates the current sweep edge from the reference
        // point.  Whenever the sweep_a_ or sweep_b_ is updated we would just need
        // to check that the new vertex is still on the opposite side of the
        // bounding plane (i.e., one dot product).  If not, we test the current
        // triangle using the code below and then compute a new bounding plane.
        //
        // Another optimization would be to choose the reference point to be 90
        // degrees away from the first input vertex, since then triangle tests would
        // not be needed unless the input geometry spans more than 90 degrees.  This
        // would involve adding a new flag have_ref_point_ rather than always
        // choosing the reference point to be S2.Origin().
        //
        // According to profiling these optimizations are not currently worthwhile,
        // but this is worth revisiting if and when other improvements are made.
        int sign = S2Pred.Sign(a, b, c);
        if (sign == 0) return;
        bool inside = S2.AngleContainsVertex(a, b, c) == (sign > 0);
        S2EdgeCrosser crosser = new(b, ref_point_);
        inside ^= crosser.EdgeOrVertexCrossing(a, b);
        inside ^= crosser.EdgeOrVertexCrossing(b, c);
        inside ^= crosser.EdgeOrVertexCrossing(c, a);
        if (inside) ref_winding_ += sign;
    }

    // Ensures that the output will be the full polygon.
    private void AddFullPolygon()
    {
        ref_winding_ += 1;
    }

    // Returns the edge normal for the given edge AB.  The sign is chosen such
    // that the normal is on the right of AB if buffer_sign_ > 0, and on the left
    // of AB if buffer_sign_ < 0.
    private S2Point GetEdgeAxis(S2Point a, S2Point b)
    {
        System.Diagnostics.Debug.Assert(buffer_sign_ != 0);
        return buffer_sign_ * S2.RobustCrossProd(b, a).Normalize();
    }

    // Adds a semi-open offset arc around vertex V.  The arc proceeds CCW from
    // "start" to "end" (both of which must be perpendicular to V).
    private void AddVertexArc(S2Point v, S2Point start, S2Point end)
    {
        // Make sure that we output at least one point even when span == 0.
        S2Point rotate_dir = buffer_sign_ * v.CrossProd(start).Normalize();
        S1ChordAngle angle = new(), span = new(start, end);
        do
        {
            S2Point dir = S2.GetPointOnRay(start, rotate_dir, angle);
            AddOffsetVertex(S2.GetPointOnRay(v, dir, abs_radius_));
        } while ((angle += vertex_step_) < span);
    }

    // Closes the semi-open arc generated by AddVertexArc().
    private void CloseVertexArc(S2Point v, S2Point end)
    {
        AddOffsetVertex(S2.GetPointOnRay(v, end, abs_radius_));
    }

    // Adds a semi-open offset arc for the given edge AB.
    private void AddEdgeArc(S2Point a, S2Point b)
    {
        S2Point ab_axis = GetEdgeAxis(a, b);
        if (edge_step_ == S1ChordAngle.Zero)
        {
            // If the buffer radius is more than 90 degrees, edges do not contribute to
            // the buffered boundary.  Instead we force the offset path to pass
            // through a vertex located at the edge normal.  This is similar to the
            // case of concave vertices (below) where it is necessary to route the
            // offset path through the concave vertex to ensure that the winding
            // numbers in all output regions have the correct sign.
            AddOffsetVertex(ab_axis);
        }
        else
        {
            // Make sure that we output at least one point even when span == 0.
            S2Point rotate_dir = buffer_sign_ * a.CrossProd(ab_axis).Normalize();
            S1ChordAngle angle=new(), span=new(a, b);
            do
            {
                S2Point p = S2.GetPointOnRay(a, rotate_dir, angle);
                AddOffsetVertex(S2.GetPointOnRay(p, ab_axis, abs_radius_));
            } while ((angle += edge_step_) < span);
        }
        SetInputVertex(b);
    }

    // Closes the semi-open arc generated by AddEdgeArc().
    private void CloseEdgeArc(S2Point a, S2Point b)
    {
        if (edge_step_ != S1ChordAngle.Zero)
        {
            AddOffsetVertex(S2.GetPointOnRay(b, GetEdgeAxis(a, b), abs_radius_));
        }
    }

    // Buffers the edge AB and the vertex B.  (The vertex C is used to determine
    // the range of angles that should be buffered at B.)
    //
    // TODO(ericv): Let A* denote the possible offset points of A with respect to
    // the edge AB for buffer radii in the range specified by "radius" and
    // "error_fraction".  Rather than requiring that the path so far terminates at
    // a point in A*, as you might expect, instead we only require that the path
    // terminates at a point X such that for any point Y in A*, the edge XY does
    // not leave the valid buffer zone of the previous edge and vertex.
    private void BufferEdgeAndVertex(S2Point a, S2Point b, S2Point c)
    {
        System.Diagnostics.Debug.Assert(a != b);
        System.Diagnostics.Debug.Assert(b != c);
        System.Diagnostics.Debug.Assert(buffer_sign_ != 0);
        if (!tracker_.Ok()) return;

        // For left (convex) turns we need to add an offset arc.  For right
        // (concave) turns we connect the end of the current offset path to the
        // vertex itself and then to the start of the offset path for the next edge.
        // Note that A == C is considered to represent a convex (left) turn.
        AddEdgeArc(a, b);
        if (buffer_sign_ * S2Pred.Sign(a, b, c) >= 0)
        {
            // The boundary makes a convex turn.  If there is no following edge arc
            // then we need to generate a closed vertex arc.
            S2Point start = GetEdgeAxis(a, b);
            S2Point end = GetEdgeAxis(b, c);
            AddVertexArc(b, start, end);
            if (edge_step_ == S1ChordAngle.Zero) CloseVertexArc(b, end);
        }
        else
        {
            // The boundary makes a concave turn.  It is tempting to simply connect
            // the end of the current offset path to the start of the offset path for
            // the next edge, however this can create output regions where the winding
            // number is incorrect.  A solution that always works is to terminate the
            // current offset path and start a new one by connecting the two offset
            // paths through the input vertex whenever it is concave.  We first need
            // to close the previous semi-open edge arc if necessary.
            CloseEdgeArc(a, b);
            AddOffsetVertex(b);  // Connect through the input vertex.
        }
    }

    // Given a polyline that starts with the edge AB, adds an end cap (as
    // specified by end_cap_style() and polyline_side()) for the vertex A.
    private void AddStartCap(S2Point a, S2Point b)
    {
        S2Point axis = GetEdgeAxis(a, b);
        if (Options_.end_cap_style_ == EndCapStyle.FLAT)
        {
            // One-sided flat end caps require no additional vertices since the
            // "offset curve" for the opposite side is simply the reversed polyline.
            if (Options_.polyline_side_ == PolylineSide.BOTH)
            {
                AddOffsetVertex(S2.GetPointOnRay(a, -axis, abs_radius_));
            }
        }
        else
        {
            System.Diagnostics.Debug.Assert(Options_.end_cap_style_ == EndCapStyle.ROUND);
            if (Options_.polyline_side_ == PolylineSide.BOTH)
            {
                // The end cap consists of a semicircle.
                AddVertexArc(a, -axis, axis);
            }
            else
            {
                // The end cap consists of a quarter circle.  Note that for
                // PolylineSide.LEFT, the polyline direction has been reversed.
                AddVertexArc(a, axis.CrossProd(a).Normalize(), axis);
            }
        }
    }

    // Given a polyline that ends with the edge AB, adds an end cap (as specified
    // by end_cap_style() and polyline_side()) for the vertex B.
    private void AddEndCap(S2Point a, S2Point b)
    {
        S2Point axis = GetEdgeAxis(a, b);
        if (Options_.end_cap_style_ == EndCapStyle.FLAT)
        {
            CloseEdgeArc(a, b);  // Close the previous semi-open edge arc if necessary.
        }
        else
        {
            System.Diagnostics.Debug.Assert(Options_.end_cap_style_ == EndCapStyle.ROUND);
            if (Options_.polyline_side_ == PolylineSide.BOTH)
            {
                // The end cap consists of a semicircle.
                AddVertexArc(b, axis, -axis);
            }
            else
            {
                // The end cap consists of a quarter circle.  We close the arc since it
                // will be followed by the reversed polyline vertices.  Note that for
                // PolylineSide.LEFT, the polyline direction has been reversed.
                S2Point end = b.CrossProd(axis).Normalize();
                AddVertexArc(b, axis, end);
                CloseVertexArc(b, end);
            }
        }
    }

    // Helper function that buffers the given loop.
    private void BufferLoop(S2PointLoopSpan loop)
    {
        // Empty loops always yield an empty path.
        if (!loop.Any() || !tracker_.Ok()) return;

        // Loops with one degenerate edge are treated as points.
        if (loop.Count == 1)
        {
            AddPoint(loop[0]);
            return;
        }

        // Buffering by 180 degrees or more always yields the full polygon.
        // Buffering by -180 degrees or more always yields the empty polygon.
        if (abs_radius_ >= S1ChordAngle.Straight)
        {
            if (buffer_sign_ > 0) AddFullPolygon();
            return;
        }

        // If buffer_radius == 0, the loop is passed through unchanged.
        if (buffer_sign_ == 0)
        {
            if (!tracker_.AddSpace(path_, loop.Count)) return;
            path_ = loop.ToList();
        }
        else
        {
            SetInputVertex(loop[0]);
            for (int i = 0; i < loop.Count; ++i)
            {
                BufferEdgeAndVertex(loop[i], loop[i + 1], loop[i + 2]);
            }
            CloseBufferRegion();
        }
        OutputPath();
    }

    private void BufferShape(S2Shape shape)
    {
        int dimension = shape.Dimension();
        int num_chains = shape.NumChains();
        for (int c = 0; c < num_chains; ++c)
        {
            S2Shape.Chain chain = shape.GetChain(c);
            if (chain.Length == 0) continue;
            if (dimension == 0)
            {
                AddPoint(shape.GetEdge(c).V0);
            }
            else
            {
                S2.GetChainVertices(shape, c, out tmp_vertices_);
                if (dimension == 1)
                {
                    AddPolyline(new S2PointSpan(tmp_vertices_));
                }
                else
                {
                    BufferLoop(new S2PointLoopSpan(tmp_vertices_));
                }
            }
        }
    }

    // The number of layers containing two-dimension geometry that have been
    // added so far.  This is used to enforce the requirement that negative
    // buffer radii allow only a single such layer.
    private int num_polygon_layers_ = 0;

    // Parameters for buffering vertices and edges.
    private int buffer_sign_;  // The sign of buffer_radius (-1, 0, or +1).
    private S1ChordAngle abs_radius_;
    private S1ChordAngle vertex_step_, edge_step_;

    // We go to extra effort to ensure that points are transformed into regular
    // polygons.  (We don't do this for arcs in general because we would rather
    // use the allowable error to reduce the complexity of the output rather
    // than increase its symmetry.)
    private S1ChordAngle point_step_;

    // Contains the buffered loops that have been accumulated so far.
    private S2WindingOperation op_;

    // The current offset path.  When each path is completed into a loop it is
    // added to op_ (the S2WindingOperation).
    private List<S2Point> path_;

    // As buffered loops are added we keep track of the winding number of a
    // fixed reference point.  This is used to derive the winding numbers of
    // every region in the spherical partition induced by the buffered loops.
    private S2Point ref_point_;

    // The winding number associated with ref_point_.
    private int ref_winding_;

    // The endpoints of the current sweep edge.  sweep_a_ is a vertex of the
    // original geometry and sweep_b_ is a vertex of the current offset path.
    private S2Point sweep_a_, sweep_b_;

    // The starting vertices of the current input loop and offset curve.  These
    // are used to close the buffer region when a loop is completed.
    private S2Point input_start_, offset_start_;
    private bool have_input_start_, have_offset_start_;

    // Used internally as a temporary to avoid excessive memory allocation.
    private S2Point[] tmp_vertices_;

    private S2MemoryTracker.Client tracker_;

    // For polylines, specifies whether the end caps should be round or flat.
    // See Options.set_end_cap_style() below.
    public enum EndCapStyle : byte { ROUND, FLAT };

    // Specifies whether polylines should be buffered only on the left, only on
    // the right, or on both sides.
    public enum PolylineSide : byte { LEFT, RIGHT, BOTH };

    public class Options
    {
        public Options() => snap_function_ = new S2BuilderUtil.IdentitySnapFunction(S1Angle.Zero);

        // Convenience constructor that calls set_buffer_radius().
        public Options(S1Angle buffer_radius)
            : this() => buffer_radius_ = buffer_radius;

        // Options may be assigned and copied.
        public Options(Options options)
        {
            buffer_radius_ = options.buffer_radius_;
            error_fraction_ = options.error_fraction_;
            end_cap_style_ = options.end_cap_style_;
            polyline_side_ = options.polyline_side_;
            snap_function_ = (S2BuilderUtil.SnapFunction)options.snap_function_.CustomClone();
            memory_tracker_ = options.memory_tracker_;
        }

        // If positive, specifies that all points within the given radius of the
        // input geometry should be added to the output.  If negative, specifies
        // that all points within the given radius of complement of the input
        // geometry should be subtracted from the output.  If the buffer radius
        // is zero then the input geometry is passed through to the output layer
        // after first converting points and polylines into degenerate loops.
        //
        // DEFAULT: S1Angle.Zero()
        public S1Angle buffer_radius_ { get; set; } = S1Angle.Zero;

        // Specifies the allowable error when buffering, expressed as a fraction
        // of buffer_radius().  The actual buffer distance will be in the range
        // [(1-f) * r - C, (1 + f) * r + C] where "f" is the error fraction, "r"
        // is the buffer radius, and "C" is S2BufferOperation.kAbsError.
        //
        // Be aware that the number of output edges increases proportionally to
        // (1 / sqrt(error_fraction)), so setting a small value may increase the
        // size of the output considerably.
        //
        // REQUIRES: error_fraction() >= kMinErrorFraction
        // REQUIRES: error_fraction() <= 1.0
        //
        // DEFAULT: 0.01  (i.e., maximum error of 1%)
        //
        public static double kMinErrorFraction = 1e-6;
        //public double error_fraction_ { get; set; } = 0.01;
        public double error_fraction_ { get; set; } = 0.02;

        // Returns the maximum error in the buffered result for the current
        // buffer_radius(), error_fraction(), and snap_function().  Note that the
        // error due to buffering consists of both relative errors (those
        // proportional to the buffer radius) and absolute errors.  The maximum
        // relative error is controlled by error_fraction(), while the maximum
        // absolute error is about 10 nanometers on the Earth's surface and is
        // defined internally.  The error due to snapping is defined by the
        // specified snap_function().
        public S1Angle max_error()
        {
            // See comments for kMinRequestedError above.
            S2Builder.Options builder_options = new(snap_function_);
            builder_options.SplitCrossingEdges = true;
            return S1Angle.Max(kMinRequestedError, error_fraction_ * S1Angle.Abs(buffer_radius_))
                + kMaxAbsoluteInterpolationErrorS1Angle + builder_options.MaxEdgeDeviation();
        }

        // Alternatively, error_fraction() may be specified as the number of
        // polyline segments used to approximate a planar circle.  These two
        // values are related according to the formula
        //
        //    error_fraction = (1 - cos(theta)) / (1 + cos(theta))
        //                  ~= 0.25 * (theta ** 2)
        //
        // where (theta == Pi / circle_segments), i.e. error decreases
        // quadratically with the number of circle segments.
        //
        // REQUIRES: circle_segments() >= 2.0
        // REQUIRES: circle_segments() <= kMaxCircleSegments
        //           (about 1570; corresponds to kMinErrorFraction)
        //
        // DEFAULT: about 15.76 (corresponding to  error_fraction() default value)
        public double kMaxCircleSegments = 1570.7968503979573;
        public double circle_segments_
        {
            get
            {
#if false
                // This formula assumes that vertices can be placed anywhere.  TODO(ericv).
                return S2Constants.M_PI / Math.Acos((1 - error_fraction_) / (1 + error_fraction_));
#else
                // This formula assumes that all vertices are placed on the midline.
                return S2.M_PI / Math.Acos(1 - error_fraction_);
#endif
            }
            set
            {
                System.Diagnostics.Debug.Assert(value >= 2.0);
                System.Diagnostics.Debug.Assert(value <= kMaxCircleSegments);
                value = Math.Max(2.0, Math.Min(kMaxCircleSegments, value));

                // We convert value to error_fraction using planar geometry,
                // because the number of segments required to approximate a circle on the
                // sphere to within a given tolerance is not constant.  Unlike in the plane,
                // the total curvature of a circle on the sphere decreases as the area
                // enclosed by the circle increases; great circles have no curvature at all.
                // We round up when converting to ensure that we won't generate any tiny
                // extra edges.
                //
#if false
                // Note that we take advantage of both positive and negative errors when
                // approximating circles (i.e., vertices are not necessarily on the midline)
                // and thus the relationships between value and error_fraction are
                //        e = (1 - cos(Pi/n)) / (1 + cos(Pi/n))
                //        n = Pi / acos((1 - e) / (1 + e))
                double r = Math.Cos(S2Constants.M_PI / value);
                error_fraction_ = (1 - r) / (1 + r) + 1e-15;
#else
                // When all vertices are on the midline, the relationships are
                //        e = 1 - cos(Pi/n)
                //        n = Pi / acos(1 - e)
                error_fraction_ = 1 - Math.Cos(S2.M_PI / value) + 1e-15;
#endif
            }
        }

        // For polylines, specifies whether the end caps should be round or flat.
        //
        // Note that with flat end caps, there is no buffering beyond the polyline
        // endpoints (unlike "square" end caps, which are not implemented).
        //
        // DEFAULT: EndCapStyle.ROUND
        public EndCapStyle end_cap_style_ { get; set; } = EndCapStyle.ROUND;

        // Specifies whether polylines should be buffered only on the left, only
        // on the right, or on both sides.  For one-sided buffering please note
        // the following:
        //
        //  - EndCapStyle.ROUND yields two quarter-circles, one at each end.
        //
        //  - To buffer by a different radius on each side of the polyline, you
        //    can use two S2BufferOperations and compute their union.  (Note that
        //    round end caps will yield two quarter-circles at each end of the
        //    polyline with different radii.)
        //
        //  - Polylines consisting of a single degenerate edge are always buffered
        //    identically to points, i.e. this option has no effect.
        //
        //  - When the polyline turns right by more than 90 degrees, buffering may
        //    or may not extend to the non-buffered side of the polyline.  For
        //    example if ABC makes a 170 degree right turn at B, it is unspecified
        //    whether the buffering of AB extends across edge BC and vice versa.
        //    Similarly if ABCD represents two right turns of 90 degrees where AB
        //    and CD are separated by less than the buffer radius, it is
        //    unspecified whether buffering of AB extends across CD and vice versa.
        //
        // DEFAULT: PolylineSide.BOTH
        public PolylineSide polyline_side_ { get; set; } = PolylineSide.BOTH;

        // Specifies the function used for snap rounding the output during the
        // call to Build().  Note that any errors due to snapping are in addition
        // to those specified by error_fraction().
        //
        // DEFAULT: s2builderutil.IdentitySnapFunction(S1Angle.Zero())
        public S2BuilderUtil.SnapFunction snap_function_
        { get => snap_function__; set => snap_function__ = (S2BuilderUtil.SnapFunction)value.CustomClone(); }
        private S2BuilderUtil.SnapFunction snap_function__;

        // Specifies that internal memory usage should be tracked using the given
        // S2MemoryTracker.  If a memory limit is specified and more more memory
        // than this is required then an error will be returned.  Example usage:
        //
        //   S2MemoryTracker tracker;
        //   tracker.set_limit(500 << 20);  // 500 MB
        //   S2BufferOperation.Options options;
        //   options.set_buffer_radius(S1Angle.Degrees(1e-5));
        //   options.set_memory_tracker(&tracker);
        //   S2BufferOperation op{options};
        //   ...
        //   S2Error error;
        //   if (!op.Build(&error)) {
        //     if (error.code() == S2Error.RESOURCE_EXHAUSTED) {
        //       S2_LOG(ERROR) << error;  // Memory limit exceeded
        //     }
        //   }
        //
        // CAVEATS:
        //
        //  - Memory allocated by the output S2Builder layer is not tracked.
        //
        //  - While memory tracking is reasonably complete and accurate, it does
        //    not account for every last byte.  It is intended only for the
        //    purpose of preventing clients from running out of memory.
        //
        // DEFAULT: nullptr (memory tracking disabled)
        public S2MemoryTracker memory_tracker_ { get; set; } = null;
    }
}
