// S2Builder is a tool for assembling polygonal geometry from edges.  Here are
// some of the things it is designed for:
//
// 1. Building polygons, polylines, and polygon meshes from unsorted
//    collections of edges.
//
// 2. Snapping geometry to discrete representations (such as S2CellId centers
//    or E7 lat/lng coordinates) while preserving the input topology and with
//    guaranteed error bounds.
//
// 3. Simplifying geometry (e.g. for indexing, display, or storage).
//
// 4. Importing geometry from other formats, including repairing geometry
//    that has errors.
//
// 5. As a tool for implementing more complex operations such as polygon
//    intersections and unions.
//
// The implementation is based on the framework of "snap rounding".  Unlike
// most snap rounding implementations, S2Builder defines edges as geodesics on
// the sphere (straight lines) and uses the topology of the sphere (i.e.,
// there are no "seams" at the poles or 180th meridian).  The algorithm is
// designed to be 100% robust for arbitrary input geometry.  It offers the
// following properties:
//
//   - Guaranteed bounds on how far input vertices and edges can move during
//     the snapping process (i.e., at most the given "snap_radius").
//
//   - Guaranteed minimum separation between edges and vertices other than
//     their endpoints (similar to the goals of Iterated Snap Rounding).  In
//     other words, edges that do not intersect in the output are guaranteed
//     to have a minimum separation between them.
//
//   - Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the
//     input already meets the output criteria then it will not be modified.
//
//   - Preservation of the input topology (up to the creation of
//     degeneracies).  This means that there exists a continuous deformation
//     from the input to the output such that no vertex crosses an edge.  In
//     other words, self-intersections won't be created, loops won't change
//     orientation, etc.
//
//   - The ability to snap to arbitrary discrete point sets (such as S2CellId
//     centers, E7 lat/lng points on the sphere, or simply a subset of the
//     input vertices), rather than being limited to an integer grid.
//
// Here are some of its other features:
//
//  - It can handle both directed and undirected edges.  Undirected edges can
//    be useful for importing data from other formats, e.g. where loops have
//    unspecified orientations.
//
//  - It can eliminate self-intersections by finding all edge pairs that cross
//    and adding a new vertex at each intersection point.
//
//  - It can simplify polygons to within a specified tolerance.  For example,
//    if two vertices are close enough they will be merged, and if an edge
//    passes nearby a vertex then it will be rerouted through that vertex.
//    Optionally, it can also detect nearly straight chains of short edges and
//    replace them with a single long edge, while maintaining the same
//    accuracy, separation, and topology guarantees ("simplify_edge_chains").
//
//  - It supports many different output types through the concept of "layers"
//    (polylines, polygons, polygon meshes, etc).  You can build multiple
//    layers at once in order to ensure that snapping does not create
//    intersections between different objects (for example, you can simplify a
//    set of contour lines without the risk of having them cross each other).
//
//  - It supports edge labels, which allow you to attach arbitrary information
//    to edges and have it preserved during the snapping process.  (This can
//    also be achieved using layers, at a coarser level of granularity.)
//
// Caveats:
//
//  - Because S2Builder only works with edges, it cannot distinguish between
//    the empty and full polygons.  If your application can generate both the
//    empty and full polygons, you must implement logic outside of this class.
//
// Example showing how to snap a polygon to E7 coordinates:
//
//  using s2builderutil.IntLatLngSnapFunction;
//  S2Builder builder(S2Builder.Options(IntLatLngSnapFunction(7)));
//  S2Polygon output;
//  builder.StartLayer(new s2builderutil.S2PolygonLayer(output));
//  builder.AddPolygon(input);
//  S2Error error;
//  if (!builder.Build(out error)) {
//    S2_LOG(ERROR) << error;
//    ...
//  }
//
// The algorithm is based on the idea of choosing a set of sites and computing
// their "limited radius Voronoi diagram", which is obtained by intersecting
// each Voronoi region with a disc of fixed radius (the "snap radius")
// centered around the corresponding site.
//
// For each input edge, we then determine the sequence of Voronoi regions
// crossed by that edge, and snap the edge to the corresponding sequence of
// sites.  (In other words, each input edge is replaced by an edge chain.)
//
// The sites are chosen by starting with the set of input vertices, optionally
// snapping them to discrete point set (such as S2CellId centers or lat/lng E7
// coordinates), and then choosing a subset such that no two sites are closer
// than the given "snap_radius".  Note that the sites do not need to be spaced
// regularly -- their positions are completely arbitrary.
//
// Rather than computing the full limited radius Voronoi diagram, instead we
// compute on demand the sequence of Voronoi regions intersected by each edge.
// We do this by first finding all the sites that are within "snap_radius" of
// the edge, sorting them by distance from the edge origin, and then using an
// incremental algorithm.
//
// We implement the minimum edge-vertex separation property by snapping all
// the input edges and checking whether any site (the "site to avoid") would
// then be too close to the edge.  If so we add another site (the "separation
// site") along the input edge, positioned so that the new snapped edge is
// guaranteed to be far enough away from the site to avoid.  We then find all
// the input edges that are within "snap_radius" of the new site, and resnap
// those edges.  (It is very rare that we need to add a separation site, even
// when sites are densely spaced.)
//
// Idempotency is implemented by explicitly checking whether the input
// geometry already meets the output criteria.  This is not as sophisticated
// as Stable Snap Rounding (Hershberger); I have worked out the details and it
// is possible to adapt those ideas here, but it would make the implementation
// significantly more complex.
//
// The only way that different output layers interact is in the choice of
// Voronoi sites:
//
//  - Vertices from all layers contribute to the initial selection of sites.
//
//  - Edges in any layer that pass too close to a site can cause new sites to
//    be added (which affects snapping in all layers).
//
//  - Simplification can be thought of as removing sites.  A site can be
//    removed only if the snapped edges stay within the error bounds of the
//    corresponding input edges in all layers.
//
// Otherwise all layers are processed independently.  For example, sibling
// edge pairs can only cancel each other within a single layer (if desired).

// Internal flag intended to be set from within a debugger.
// #define s2builder_verbose = false;

namespace S2Geometry;

using System.Runtime.InteropServices;
using InputVertexId = Int32;
using SiteId = Int32;
using S2BuilderUtil;
using ShapeEdge = S2ShapeUtil.ShapeEdge;
using System.Diagnostics;

public partial class S2Builder
{
    // For output layers that represent polygons, there is an ambiguity inherent
    // in spherical geometry that does not exist in planar geometry.  Namely, if
    // a polygon has no edges, does it represent the empty polygon (containing
    // no points) or the full polygon (containing all points)?  This ambiguity
    // also occurs for polygons that consist only of degeneracies, e.g. a
    // degenerate loop with only two edges could be either a degenerate shell in
    // the empty polygon or a degenerate hole in the full polygon.
    //
    // To resolve this ambiguity, an IsFullPolygonPredicate may be specified for
    // each output layer (see AddIsFullPolygonPredicate below).  If the output
    // after snapping consists only of degenerate edges and/or sibling pairs
    // (including the case where there are no edges at all), then the layer
    // implementation calls the given predicate to determine whether the polygon
    // is empty or full except for those degeneracies.  The predicate is given
    // an S2Builder.Graph containing the output edges, but note that in general
    // the predicate must also have knowledge of the input geometry in order to
    // determine the correct result.
    //
    // This predicate is only needed by layers that are assembled into polygons.
    // It is not used by other layer types.
    public delegate bool IsFullPolygonPredicate(Graph g, out S2Error error);

    // Convenience constructor that calls Init().  Note that to use the default
    // options, C++ syntax requires an extra layer of parentheses:
    //
    //   S2Builder builder{S2Builder.Options()};
    public S2Builder(Options options)
    {
        Options_ = options;
        var snap_function = options.SnapFunction;
        var snap_radius = snap_function.SnapRadius;
        Debug.Assert(snap_radius <= SnapFunction.kMaxSnapRadius);

        // Convert the snap radius to an S1ChordAngle.  This is the "true snap
        // radius" used when evaluating exact predicates (s2predicates.h).
        site_snap_radius_ca_ = new S1ChordAngle(snap_radius);

        // When intersection_tolerance() is non-zero we need to use a larger snap
        // radius for edges than for vertices to ensure that both edges are snapped
        // to the edge intersection location.  This is because the computed
        // intersection point is not exact; it may be up to intersection_tolerance()
        // away from its true position.  The computed intersection point might then
        // be snapped to some other vertex up to snap_radius away.  So to ensure
        // that both edges are snapped to a common vertex, we need to increase the
        // snap radius for edges to at least the sum of these two values (calculated
        // conservatively).
        S1Angle edge_snap_radius = options.EdgeSnapRadius();
        edge_snap_radius_ca_ = RoundUp(edge_snap_radius);
        snapping_requested_ = edge_snap_radius > S1Angle.Zero;

        // Compute the maximum distance that a vertex can be separated from an
        // edge while still affecting how that edge is snapped.
        max_edge_deviation_ = options.MaxEdgeDeviation();
        edge_site_query_radius_ca_ = new S1ChordAngle(
            max_edge_deviation_ + snap_function.MinEdgeVertexSeparation());

        // Compute the maximum edge length such that even if both endpoints move by
        // the maximum distance allowed (i.e., edge_snap_radius), the center of the
        // edge will still move by less than max_edge_deviation().  This saves us a
        // lot of work since then we don't need to check the actual deviation.
        if (!snapping_requested_)
        {
            min_edge_length_to_split_ca_ = S1ChordAngle.Infinity;
        }
        else
        {
            // This value varies between 30 and 50 degrees depending on the snap radius.
            min_edge_length_to_split_ca_ = S1ChordAngle.FromRadians(
                2 * Math.Acos(edge_snap_radius.Sin() / max_edge_deviation_.Sin()));
        }

        // In rare cases we may need to explicitly check that the input topology is
        // preserved, i.e. that edges do not cross vertices when snapped.  This is
        // only necessary (1) for vertices added using ForceVertex(), and (2) when the
        // snap radius is smaller than intersection_tolerance() (which is typically
        // either zero or S2::kIntersectionError, about 9e-16 radians).  This
        // condition arises because when a geodesic edge is snapped, the edge center
        // can move further than its endpoints.  This can cause an edge to pass on the
        // wrong side of an input vertex.  (Note that this could not happen in a
        // planar version of this algorithm.)  Usually we don't need to consider this
        // possibility explicitly, because if the snapped edge passes on the wrong
        // side of a vertex then it is also closer than min_edge_vertex_separation()
        // to that vertex, which will cause a separation site to be added.
        //
        // If the condition below is true then we need to check all sites (i.e.,
        // snapped input vertices) for topology changes.  However this is almost never
        // the case because
        //
        //            max_edge_deviation() == 1.1 * edge_snap_radius()
        //      and   min_edge_vertex_separation() >= 0.219 * snap_radius()
        //
        // for all currently implemented snap functions.  The condition below is
        // only true when intersection_tolerance() is non-zero (which causes
        // edge_snap_radius() to exceed snap_radius() by S2::kIntersectionError) and
        // snap_radius() is very small (at most S2::kIntersectionError / 1.19).
        check_all_site_crossings_ = (options.MaxEdgeDeviation() >
                                     options.EdgeSnapRadius() +
                                     snap_function.MinEdgeVertexSeparation());
        if (options.IntersectionTolerance <= S1Angle.Zero)
        {
            Debug.Assert(!check_all_site_crossings_);
        }

        // To implement idempotency, we check whether the input geometry could
        // possibly be the output of a previous S2Builder invocation.  This involves
        // testing whether any site/site or edge/site pairs are too close together.
        // This is done using exact predicates, which require converting the minimum
        // separation values to an S1ChordAngle.
        min_site_separation_ = snap_function.MinVertexSeparation();
        min_site_separation_ca_ = new(min_site_separation_);
        min_edge_site_separation_ca_ =
            new(snap_function.MinEdgeVertexSeparation());

        // This is an upper bound on the distance computed by S2ClosestPointQuery
        // where the true distance might be less than min_edge_site_separation_ca_.
        min_edge_site_separation_ca_limit_ =
            AddPointToEdgeError(min_edge_site_separation_ca_);

        // Compute the maximum possible distance between two sites whose Voronoi
        // regions touch.  (The maximum radius of each Voronoi region is
        // edge_snap_radius_.)  Then increase this bound to account for errors.
        max_adjacent_site_separation_ca_ =
            AddPointToPointError(RoundUp(2 * edge_snap_radius));

        // Finally, we also precompute sin^2(edge_snap_radius), which is simply the
        // squared distance between a vertex and an edge measured perpendicular to
        // the plane containing the edge, and increase this value by the maximum
        // error in the calculation to compare this distance against the bound.
        var d = edge_snap_radius.Sin();
        edge_snap_radius_sin2_ = d * d;
        edge_snap_radius_sin2_ += ((9.5 * d + 2.5 + 2 * Math.Sqrt(3)) * d + 9 * S2.DoubleEpsilon) * S2.DoubleEpsilon;

        // Initialize the current label set.
        label_set_id_ = IdSetLexicon.kEmptySetId;
        label_set_modified_ = false;

        // If snapping was requested, we try to determine whether the input geometry
        // already meets the output requirements.  This is necessary for
        // idempotency, and can also save work.  If we discover any reason that the
        // input geometry needs to be modified, snapping_needed_ is set to true.
        snapping_needed_ = false;

        tracker_.Init(options.MemoryTracker);
    }

    // Starts a new output layer.  This method must be called before adding any
    // edges to the S2Builder.  You may call this method multiple times to build
    // multiple geometric objects that are snapped to the same set of sites.
    //
    // For example, if you have a set of contour lines, then you could put each
    // contour line in a separate layer.  This keeps the contour lines separate
    // from each other, while also ensuring that no crossing edges are created
    // when they are snapped and/or simplified.  (This is not true if the
    // contour lines are snapped or simplified independently.)
    //
    // Similarly, if you have a set of polygons that share common boundaries
    // (e.g., countries), you can snap and/or simplify them at the same time by
    // putting them in different layers, while ensuring that their boundaries
    // remain consistent (i.e., no crossing edges or T-vertices are introduced).
    //
    // Ownership of the layer is transferred to the S2Builder.  Example usage:
    //
    // S2Polyline line1, line2;
    // builder.StartLayer(new s2builderutil.S2PolylineLayer(line1)));
    // ... Add edges using builder.AddEdge(), etc ...
    // builder.StartLayer(new s2builderutil.S2PolylineLayer(line2)));
    // ... Add edges using builder.AddEdge(), etc ...
    // S2Error error;
    // Debug.Assert(builder.Build(out error)) << error;  // Builds "line1" & "line2"
    public void StartLayer(Layer layer)
    {
        layer_options_.Add(layer.GraphOptions_());
        layer_begins_.Add(input_edges_.Count);
        layer_is_full_polygon_predicates_.Add(IsFullPolygon(false));
        layers_.Add(layer);
    }

    // Adds a degenerate edge (representing a point) to the current layer.
    public void AddPoint(S2Point v)
    {
        AddEdge(v, v);
    }

    // Adds the given edge to the current layer.
    public void AddEdge(S2Point v0, S2Point v1)
    {
        Debug.Assert(layers_.Any()); // Call StartLayer before adding any edges
        if (v0 == v1 && (layer_options_.Last().DegenerateEdges_ ==
                         GraphOptions.DegenerateEdges.DISCARD))
        {
            return;
        }
        var j0 = AddVertex(v0);
        var j1 = AddVertex(v1);
        if (!tracker_.AddSpace(input_edges_, 1)) return;
        input_edges_.Add(new(j0, j1));

        // If there are any labels, then attach them to this input edge.
        if (label_set_modified_)
        {
            if (!label_set_ids_.Any())
            {
                // Populate the missing entries with empty label sets.
                label_set_ids_.Fill(label_set_id_, input_edges_.Count - 1);
            }
            label_set_id_ = label_set_lexicon_.Add(label_set_);
            label_set_ids_.Add(label_set_id_);
            label_set_modified_ = false;
        }
        else if (label_set_ids_.Any())
        {
            label_set_ids_.Add(label_set_id_);
        }
    }

    // Adds the edges in the given polyline.  Note that polylines with 0 or 1
    // vertices are defined to have no edges.
    public void AddPolyline(S2PointSpan polyline)
    {
        for (int i = 1; i < polyline.Count; ++i)
        {
            AddEdge(polyline[i - 1], polyline[i]);
        }
    }

    public void AddPolyline(S2Polyline polyline)
    {
        var n = polyline.NumVertices();
        for (var i = 1; i < n; ++i)
        {
            AddEdge(polyline.Vertex(i - 1), polyline.Vertex(i));
        }
    }

    // Adds the edges in the given loop.  Note that a loop consisting of one
    // vertex adds a single degenerate edge.
    //
    // If the sign() of an S2Loop is negative (i.e. the loop represents a hole
    // within a polygon), the edge directions are automatically reversed to
    // ensure that the polygon interior is always to the left of every edge.
    public void AddLoop(S2PointLoopSpan loop)
    {
        for (int i = 0; i < loop.Count; ++i)
        {
            AddEdge(loop[i], loop.GetPoint(i + 1));
        }
    }

    public void AddLoop(S2Loop loop)
    {
        // Ignore loops that do not have a boundary.
        if (loop.IsEmptyOrFull()) return;

        // For loops that represent holes, we add the edge from vertex n-1 to vertex
        // n-2 first.  This is because these edges will be assembled into a
        // clockwise loop, which will eventually be normalized in S2Polygon by
        // calling S2Loop.Invert().  S2Loop.Invert() reverses the order of the
        // vertices, so to end up with the original vertex order (0, 1, ..., n-1) we
        // need to build a clockwise loop with vertex order (n-1, n-2, ..., 0).
        // This is done by adding the edge (n-1, n-2) first, and then ensuring that
        // Build() assembles loops starting from edges in the order they were added.
        var n = loop.NumVertices;
        for (var i = 0; i < n; ++i)
        {
            AddEdge(loop.OrientedVertex(i), loop.OrientedVertex(i + 1));
        }
    }

    // Adds the loops in the given polygon.  Loops representing holes have their
    // edge directions automatically reversed as described for AddLoop().  Note
    // that this method does not distinguish between the empty and full polygons,
    // i.e. adding a full polygon has the same effect as adding an empty one.
    public void AddPolygon(S2Polygon polygon)
    {
        for (var i = 0; i < polygon.NumLoops(); ++i)
        {
            AddLoop(polygon.Loop(i));
        }
    }

    // If "vertex" is the intersection point of two edges AB and CD (as computed
    // by S2::GetIntersection()), this method ensures that AB and CD snap to a
    // common vertex.  (Note that the common vertex may be different than
    // "vertex" in order to ensure that no pair of vertices is closer than the
    // given snap radius.)  Unlike Options::split_crossing_edges(), this method
    // may be used to split crossing edge pairs selectively.
    //
    // This method can also be used to tessellate edges using S2::GetPointOnLine()
    // or S2::Project() provided that a suitable intersection tolerance is
    // specified (see intersection_tolerance() for details).
    //
    // This method implicitly overrides the idempotent() option, since adding an
    // intersection point implies a desire to have nearby edges snapped to it
    // even if these edges already satsify the S2Builder output guarantees.
    // (Otherwise for example edges would never be snapped to nearby
    // intersection points when the snap radius is zero.)
    //
    // Note that unlike ForceVertex(), this method maintains all S2Builder
    // guarantees regarding minimum vertex-vertex separation, minimum
    // edge-vertex separation, and edge chain simplification.
    //
    // REQUIRES: options().intersection_tolerance() > S1Angle::Zero()
    // REQUIRES: "vertex" was computed by S2::GetIntersection() (in order to
    //           guarantee that both edges snap to a common vertex)
    public void AddIntersection(S2Point vertex)
    {
        // It is an error to call this method without first setting
        // intersection_tolerance() to a non-zero value.
        Debug.Assert(Options_.IntersectionTolerance > S1Angle.Zero);

        // Calling this method also overrides the idempotent() option.
        snapping_needed_ = true;

        AddVertex(vertex);
    }

    // Adds the edges of the given shape to the current layer.
    public void AddShape(S2Shape shape)
    {
        for (int e = 0, n = shape.NumEdges(); e < n; ++e)
        {
            var edge = shape.GetEdge(e);
            AddEdge(edge.V0, edge.V1);
        }
    }

    // For layers that are assembled into polygons, this method specifies a
    // predicate that is called when the output consists entirely of degenerate
    // edges and/or sibling pairs.  The predicate is given an S2Builder.Graph
    // containing the output edges (if any) and is responsible for deciding
    // whether this graph represents the empty polygon (possibly with degenerate
    // shells) or the full polygon (possibly with degenerate holes).  Note that
    // this cannot be determined from the output edges alone; it also requires
    // knowledge of the input geometry.  (Also see IsFullPolygonPredicate above.)
    //
    // This method should be called at most once per layer; additional calls
    // simply overwrite the previous value for the current layer.
    //
    // The default predicate simply returns false (i.e., degenerate polygons are
    // assumed to be empty).  Arguably it would better to return an error in
    // this case, but the fact is that relatively few clients need to be able to
    // construct full polygons, and it is unreasonable to expect all such
    // clients to supply an appropriate predicate.
    //
    // The reason for having a predicate rather than a boolean value is that the
    // predicate is responsible for determining whether the output polygon is
    // empty or full.  In general the input geometry is not degenerate, but
    // rather collapses into a degenerate configuration due to snapping and/or
    // simplification.
    //
    // TODO(ericv): Provide standard predicates to handle common cases,
    // e.g. valid input geometry that becomes degenerate due to snapping.
    public void AddIsFullPolygonPredicate(IsFullPolygonPredicate predicate)
    {
        layer_is_full_polygon_predicates_[^1] = predicate;
    }

    // A predicate that returns an error indicating that no polygon predicate
    // has been specified.
    public static bool IsFullPolygonUnspecified(out S2Error error)
    {
        error = new(S2ErrorCode.BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED, "A degenerate polygon was found, but no predicate was specified to determine whether the polygon is empty or full.  Call S2Builder.AddIsFullPolygonPredicate() to fix this problem.");
        // Assumes the polygon is empty.
        return false;
    }

    // Returns a predicate that returns a constant value (true or false);
    public static IsFullPolygonPredicate IsFullPolygon(bool is_full)
    {
        return (Graph g, out S2Error error) =>
        {
            error = S2Error.OK;
            return is_full;
        };
    }

    // Forces a vertex to be located at the given position.  This can be used to
    // prevent certain input vertices from moving.  However if you are trying to
    // preserve input edges, be aware that this option does not prevent edges from
    // being split by new vertices.
    //
    // Forced vertices are subject to the following limitations:
    //
    //  - Forced vertices are never snapped.  This is true even when the given
    //    position is not allowed by the given snap function (e.g. you can force
    //    a vertex at a non-S2CellId center when using S2CellIdSnapFunction).
    //    If you want to ensure that forced vertices obey the snap function
    //    restrictions, you must call snap_function().SnapPoint() explicitly.
    //
    //  - There is no guaranteed minimum separation between pairs of forced
    //    vertices, i.e. snap_function().min_vertex_separation() does not apply.
    //    (This must be true because forced vertices can be placed arbitrarily.)
    //
    //  - There is no guaranteed minimum separation between forced vertices and
    //    non-incident edges, i.e. snap_function().min_edge_vertex_separation()
    //    does not apply.
    //
    //  - Forced vertices are never simplified away (i.e. when simplification is
    //    requested using options().simplify_edge_chains()).
    //
    // All other guarantees continue to hold, e.g. the input topology will always
    // be preserved.
    public void ForceVertex(S2Point vertex)
    {
        if (!tracker_.AddSpace(sites_, 1)) return;
        sites_.Add(vertex);
    }

    // Every edge can have a set of non-negative integer labels attached to it.
    // When used with an appropriate layer type, you can then retrieve the
    // labels associated with each output edge.  This can be useful when merging
    // or combining data from several sources.  (Note that in many cases it is
    // easier to use separate output layers rather than labels.)
    //
    // Labels are 32-bit non-negative integers.  To support other label types,
    // you can use ValueLexicon to store the set of unique labels seen so far:
    //
    //   ValueLexicon<MyLabel> my_label_lexicon;
    //   builder.set_label(my_label_lexicon.Add(label));
    //
    // The current set of labels is represented as a stack.  This makes it easy
    // to add and remove labels hierarchically (e.g., polygon 5, loop 2).  Use
    // set_label() and clear_labels() if you need at most one label per edge.
    //

    // Clear the stack of labels.
    public void ClearLabels()
    {
        label_set_.Clear();
        label_set_modified_ = true;
    }

    // Add a label to the stack.
    // REQUIRES: label >= 0.
    public void PushLabel(int label)
    {
        Debug.Assert(label >= 0);
        label_set_.Add(label);
        label_set_modified_ = true;
    }

    // Remove a label from the stack.
    public void PopLabel()
    {
        label_set_.RemoveAt(label_set_.Count - 1);
        label_set_modified_ = true;
    }

    // Convenience function that clears the stack and adds a single label.
    // REQUIRES: label >= 0.
    public void SetLabel(int label)
    {
        Debug.Assert(label >= 0);
        label_set_.Capacity = 1;
        label_set_[0] = label;
        label_set_modified_ = true;
    }

    // Performs the requested edge splitting, snapping, simplification, etc, and
    // then assembles the resulting edges into the requested output layers.
    //
    // Returns true if all edges were assembled; otherwise sets "error"
    // appropriately.  Depending on the error, some or all output layers may
    // have been created.  Automatically resets the S2Builder state so that it
    // can be reused.
    //
    // REQUIRES: error != null.
    public bool Build(out S2Error error)
    {
        // Mark the end of the last layer.
        layer_begins_.Add(input_edges_.Count);

        // See the algorithm overview at the top of this file.
        if (snapping_requested_ && !Options_.Idempotent)
        {
            snapping_needed_ = true;
        }
        ChooseSites();
        BuildLayers();
        Reset();
        if (!tracker_.Ok()) error_ = tracker_.Error();
        error = error_;
        return error.IsOk();
    }

    // Clears all input data and resets the builder state.  Any options
    // specified are preserved.
    public void Reset()
    {
        // Note that these calls do not change vector capacities.
        input_vertices_.Clear();
        input_edges_.Clear();
        layers_.Clear();
        layer_options_.Clear();
        layer_begins_.Clear();
        layer_is_full_polygon_predicates_.Clear();
        label_set_ids_.Clear();
        label_set_lexicon_.Clear();
        label_set_.Clear();
        label_set_modified_ = false;
        sites_.Clear();
        edge_sites_.Clear();
        snapping_needed_ = false;
    }

    ///////////////////////////////////////////////////////////////////////////
    // The following methods may be called at any time, including from
    // S2Builder::Layer implementations.

    // Returns the number of input edges.
    public int NumInputEdges() => input_edges_.Count;

    // Returns the endpoints of the given input edge.
    //
    // REQUIRES: 0 <= input_edge_id < num_input_edges()
    public S2Shape.Edge InputEdge(int input_edge_id)
    {
        var edge = input_edges_[input_edge_id];
        return new S2Shape.Edge(input_vertices_[edge.Item1],
                             input_vertices_[edge.Item2]);
    }

    //////////////////////  Input Types  /////////////////////////
    // All types associated with the S2Builder inputs are prefixed with "Input".

    // Input vertices are stored in a vector, with some removal of duplicates.
    // Edges are represented as (VertexId, VertexId) pairs.  All edges are stored
    // in a single vector; each layer corresponds to a contiguous range.

    private int AddVertex(S2Point v)
    {
        // Remove duplicate vertices that follow the pattern AB, BC, CD.  If we want
        // to do anything more sophisticated, either use a ValueLexicon, or sort the
        // vertices once they have all been added, remove duplicates, and update the
        // edges.
        if (!input_vertices_.Any() || v != input_vertices_.Last())
        {
            if (!tracker_.AddSpace(input_vertices_, 1)) return -1;
            input_vertices_.Add(v);
        }
        return input_vertices_.Count - 1;
    }

    private void ChooseSites()
    {
        if (!tracker_.Ok() || !input_vertices_.Any()) return;

        // Note that although we always create an S2ShapeIndex, often it is not
        // actually built (because this happens lazily).  Therefore we only test
        // its memory usage at the places where it is used.
        MutableS2ShapeIndex input_edge_index = new();
        input_edge_index.MemoryTracker = tracker_.Tracker;
        input_edge_index.Add(new VertexIdEdgeVectorShape(input_edges_, input_vertices_));
        if (Options_.SplitCrossingEdges)
        {
            AddEdgeCrossings(input_edge_index);
        }

        if (snapping_requested_)
        {
            S2PointIndex<int> site_index = new();
            try
            {
                AddForcedSites(site_index);
                ChooseInitialSites(site_index);
                if (!tracker_.FixSiteIndexTally(site_index)) return;
                CollectSiteEdges(site_index);
            }
            finally
            {
                tracker_.DoneSiteIndex(/*site_index*/);
            }
        }
        if (snapping_needed_)
        {
            AddExtraSites(input_edge_index);
        }
        else
        {
            ChooseAllVerticesAsSites();
        }
    }
    private void ChooseAllVerticesAsSites()
    {
        // Sort the input vertices, discard duplicates, and use the result as the
        // list of sites.  (We sort in the same order used by ChooseInitialSites()
        // to avoid inconsistencies in tests.)  We also copy the result back to
        // input_vertices_ and update the input edges to use the new vertex
        // numbering (so that InputVertexId == SiteId).  This simplifies the
        // implementation of SnapEdge() for this case.
        sites_.Clear();
        if (!tracker_.AddSpaceExact(sites_, input_vertices_.Count)) return;
        Int64 kTempPerVertex = Marshal.SizeOf(typeof(InputVertexKey)) + sizeof(InputVertexId);
        if (!tracker_.TallyTemp(input_vertices_.Count * kTempPerVertex)) return;
        var sorted = SortInputVertices();
        var vmap = new int[input_vertices_.Count];
        for (var in_ = 0; in_ < sorted.Count;)
        {
            var site = input_vertices_[sorted[in_].Item2];
            vmap[sorted[in_].Item2] = sites_.Count;
            while (++in_ < sorted.Count && input_vertices_[sorted[in_].Item2] == site)
            {
                vmap[sorted[in_].Item2] = sites_.Count;
            }
            sites_.Add(site);
        }
        input_vertices_ = sites_;  // Does not change allocated size.
        for (var i = 0; i < input_edges_.Count; i++)
        {
            var e = input_edges_[i];
            input_edges_[i] = new(vmap[e.Item1], vmap[e.Item2]);
        }
    }
    private List<InputVertexKey> SortInputVertices()
    {
        // Sort all the input vertices in the order that we wish to consider them as
        // candidate Voronoi sites.  Any sort order will produce correct output, so
        // we have complete flexibility in choosing the sort key.  We could even
        // leave them unsorted, although this would have the disadvantage that
        // changing the order of the input edges could cause S2Builder to snap to a
        // different set of Voronoi sites.
        //
        // We have chosen to sort them primarily by S2CellId since this improves the
        // performance of many S2Builder phases (due to better spatial locality).
        // It also allows the possibility of replacing the current S2PointIndex
        // approach with a more efficient recursive divide-and-conquer algorithm.
        //
        // However, sorting by leaf S2CellId alone has two small disadvantages in
        // the case where the candidate sites are densely spaced relative to the
        // snap radius (e.g., when using the IdentitySnapFunction, or when snapping
        // to E6/E7 near the poles, or snapping to S2CellId/E6/E7 using a snap
        // radius larger than the minimum value required):
        //
        //  - First, it tends to bias the Voronoi site locations towards points that
        //    are earlier on the S2CellId Hilbert curve.  For example, suppose that
        //    there are two parallel rows of input vertices on opposite sides of the
        //    edge between two large S2Cells, and the rows are separated by less
        //    than the snap radius.  Then only vertices from the cell with the
        //    smaller S2CellId are selected, because they are considered first and
        //    prevent us from selecting the sites from the other cell (because they
        //    are closer than "snap_radius" to an existing site).
        //
        //  - Second, it tends to choose more Voronoi sites than necessary, because
        //    at each step we choose the first site along the Hilbert curve that is
        //    at least "snap_radius" away from all previously selected sites.  This
        //    tends to yield sites whose "coverage discs" overlap quite a bit,
        //    whereas it would be better to cover all the input vertices with a
        //    smaller set of coverage discs that don't overlap as much.  (This is
        //    the "geometric set cover problem", which is NP-hard.)
        //
        // It is not worth going to much trouble to fix these problems, because they
        // really aren't that important (and don't affect the guarantees made by the
        // algorithm), but here are a couple of heuristics that might help:
        //
        // 1. Sort the input vertices by S2CellId at a coarse level (down to cells
        // that are O(snap_radius) in size), and then sort by a fingerprint of the
        // S2Point coordinates (i.e., quasi-randomly).  This would retain most of
        // the advantages of S2CellId sorting, but makes it more likely that we will
        // select sites that are further apart.
        //
        // 2. Rather than choosing the first uncovered input vertex and snapping it
        // to obtain the next Voronoi site, instead look ahead through following
        // candidates in S2CellId order and choose the furthest candidate whose
        // snapped location covers all previous uncovered input vertices.
        //
        // TODO(ericv): Experiment with these approaches.

        List<InputVertexKey> keys = new(input_vertices_.Count);
        for (var i = 0; i < input_vertices_.Count; ++i)
        {
            keys.Add(new(new(input_vertices_[i]), i));
        }
        keys.Sort((InputVertexKey a, InputVertexKey b) =>
        {
            if (a.Item1 < b.Item1) return -1;
            if (b.Item1 < a.Item1) return 1;
            return input_vertices_[a.Item2].CompareTo(input_vertices_[b.Item2]);
        });
        return keys;
    }

    // Check all edge pairs for crossings, and add the corresponding intersection
    // points to input_vertices_.  (The intersection points will be snapped and
    // merged with the other vertices during site selection.)
    private void AddEdgeCrossings(MutableS2ShapeIndex input_edge_index)
    {
        input_edge_index.ForceBuild();
        if (!tracker_.Ok()) return;

        // We need to build a list of intersections and add them afterwards so that
        // we don't reallocate vertices_ during the VisitCrossings() call.
        List<S2Point> new_vertices = new();
        try
        {
            S2ShapeUtil.EdgePairs.VisitCrossingEdgePairs(input_edge_index, CrossingType.INTERIOR,
                (ShapeEdge a, ShapeEdge b, bool c) =>
                {
                    if (!tracker_.AddSpace(new_vertices, 1)) return false;
                    new_vertices.Add(S2.GetIntersection(a.V0, a.V1, b.V0, b.V1, null));
                    return true;  // Continue visiting.
                });
            if (!new_vertices.Any())
                return;

            snapping_needed_ = true;
            if (!tracker_.AddSpaceExact(input_vertices_, new_vertices.Count)) return;
            input_vertices_.InsertRange(input_vertices_.Count, new_vertices);
        }
        finally
        {
            tracker_.Untally(new_vertices);
        }
    }
    private void AddForcedSites(S2PointIndex<int> site_index)
    {
        // Sort the forced sites and remove duplicates.
        SortedSet<S2Point> tmp = new(sites_);
        sites_.Clear();
        sites_.AddRange(tmp);
        // Add the forced sites to the index.
        for (var id = 0; id < sites_.Count; ++id)
        {
            if (!tracker_.TallyIndexedSite()) return;
            site_index.Add(sites_[id], id);
        }
        num_forced_sites_ = sites_.Count;
    }
    private bool IsForced(int v) => v < num_forced_sites_;
    private void ChooseInitialSites(S2PointIndex<int> site_index)
    {
        // Prepare to find all points whose distance is <= min_site_separation_ca_.
        S2ClosestPointQuery<int>.Options? options = new();
        options.ConservativeMaxDistance = min_site_separation_ca_;
        S2ClosestPointQuery<int> site_query = new(site_index, options);
        List<S2ClosestPointQueryBase<S1ChordAngle, int>.Result> results = new();

        // Apply the snap_function() to each input vertex, then check whether any
        // existing site is closer than min_vertex_separation().  If not, then add a
        // new site.
        //
        // NOTE(ericv): There are actually two reasonable algorithms, which we call
        // "snap first" (the one above) and "snap last".  The latter checks for each
        // input vertex whether any existing site is closer than snap_radius(), and
        // only then applies the snap_function() and adds a new site.  "Snap last"
        // can yield slightly fewer sites in some cases, but it is also more
        // expensive and can produce surprising results.  For example, if you snap
        // the polyline "0:0, 0:0.7" using IntLatLngSnapFunction(0), the result is
        // "0:0, 0:0" rather than the expected "0:0, 0:1", because the snap radius
        // is approximately Math.Sqrt(2) degrees and therefore it is legal to snap both
        // input points to "0:0".  "Snap first" produces "0:0, 0:1" as expected.
        //
        // Track the memory used by SortInputVertices() before calling it.
        if (!tracker_.Tally(input_vertices_.Count * Marshal.SizeOf(typeof(InputVertexKey)))) return;
        var sorted_keys = SortInputVertices();
        try
        {
            foreach (var key in sorted_keys)
            {
                var vertex = input_vertices_[key.Item2];
                var site = SnapSite(vertex);
                // If any vertex moves when snapped, the output cannot be idempotent.
                snapping_needed_ = snapping_needed_ || site != vertex;

                var add_site = true;

                if (site_snap_radius_ca_ == S1ChordAngle.Zero)
                {
                    add_site = !sites_.Any() || site != sites_.Last();
                }
                else
                {
                    // FindClosestPoints() measures distances conservatively, so we need to
                    // recheck the distances using exact predicates.
                    //
                    // NOTE(ericv): When the snap radius is large compared to the average
                    // vertex spacing, we could possibly avoid the call the FindClosestPoints
                    // by checking whether sites_.Last() is close enough.
                    var target = new S2ClosestPointQuery<int>.PointTarget(site);
                    site_query.FindClosestPoints(target, results);
                    foreach (var result in results)
                    {
                        if (S2Pred.CompareDistance(site, result.Point,
                                                    min_site_separation_ca_) <= 0)
                        {
                            add_site = false;
                            // This pair of sites is too close.  If the sites are distinct, then
                            // the output cannot be idempotent.
                            snapping_needed_ = snapping_needed_ || site != result.Point;
                        }
                    }
                    if (add_site)
                    {
                        if (!tracker_.TallyIndexedSite()) return;
                        site_index.Add(site, sites_.Count);
                        if (!tracker_.AddSpace(sites_, 1)) return;
                        sites_.Add(site);
                        site_query.ReInit();
                    }
                }
            }
        }
        finally
        {
            tracker_.Untally(sorted_keys);
        }
    }
    private S2Point SnapSite(S2Point point)
    {
        if (!snapping_requested_) return point;
        var site = Options_.SnapFunction.SnapPoint(point);
        S1ChordAngle dist_moved = new(site, point);
        if (dist_moved > site_snap_radius_ca_)
        {
            error_ = new(S2ErrorCode.BUILDER_SNAP_RADIUS_TOO_SMALL, $"Snap function moved vertex ({point.X:g}, {point.Y:g}, {point.Z:g}) by {dist_moved.ToAngle().Radians:g}, which is more than the specified snap radius of {site_snap_radius_ca_.ToAngle().Radians:g}");
        }
        return site;
    }

    // For each edge, find all sites within edge_site_query_radius_ca_  and
    // store them in edge_sites_.  Also, to implement idempotency this method also
    // checks whether the input vertices and edges may already satisfy the output
    // criteria.  If any problems are found then snapping_needed_ is set to true.
    private void CollectSiteEdges(S2PointIndex<int> site_index)
    {
        // Find all points whose distance is <= edge_site_query_radius_ca_.
        //
        // Memory used by S2ClosestPointQuery is not tracked, but it is temporary,
        // typically insignificant, and does not affect the high water mark.
        S2ClosestPointQuery<int>.Options? options = new();
        options.ConservativeMaxDistance = edge_site_query_radius_ca_;
        var site_query = new S2ClosestPointQuery<int>(site_index, options);
        var results = new List<S2ClosestPointQueryBase<S1ChordAngle, int>.Result>();
        if (!tracker_.AddSpaceExact(edge_sites_, input_edges_.Count)) return;
        edge_sites_.Capacity = input_edges_.Count;  // Construct all elements.
        for (var e = 0; e < input_edges_.Count; ++e)
        {
            var edge = input_edges_[e];
            var v0 = input_vertices_[edge.Item1];
            var v1 = input_vertices_[edge.Item2];

            WriteS2Polyline(v0, v1);

            S2ClosestPointQuery<int>.EdgeTarget target = new(v0, v1);
            site_query.FindClosestPoints(target, results);
            var sites = edge_sites_[e];
            if (sites == null)
            {
                sites = new();
                edge_sites_[e] = sites;
            }
            sites.Capacity = Math.Max(results.Count, sites.Count);
            foreach (var result in results)
            {
                sites.Add(result.Data);
                if (!snapping_needed_ &&
                    result.Distance < min_edge_site_separation_ca_limit_ &&
                    result.Point != v0 && result.Point != v1 &&
                    S2Pred.CompareEdgeDistance(result.Point, v0, v1, min_edge_site_separation_ca_) < 0)
                {
                    snapping_needed_ = true;
                }
            }
            SortSitesByDistance(v0, sites);
            if (!tracker_.TallyEdgeSites(sites)) return;
        }
    }

    [Conditional("s2builder_verbose")]
    private static void WriteS2Polyline(S2Point v0, S2Point v1) =>
        Debug.WriteLine($"S2Polyline: {v0.ToDebugString()}, {v1.ToDebugString()}");

    // Sort sites in increasing order of distance to X.
    private void SortSitesByDistance(S2Point x, List<SiteId> sites) =>
        sites.Sort(new SiteIdsComp(x, sites_));

    // Like the above, but inserts "new_site_id" into an already-sorted list.
    private void InsertSiteByDistance(SiteId new_site_id, S2Point x,
                                     List<SiteId> sites)
    {
        if (!tracker_.ReserveEdgeSite(sites)) return;
        sites.Insert(sites.GetLowerBound(new_site_id, new SiteIdsComp(x, sites_)), new_site_id);
    }

    // There are two situatons where we need to add extra Voronoi sites in order to
    // ensure that the snapped edges meet the output requirements:
    //
    //  (1) If a snapped edge deviates from its input edge by more than
    //      max_edge_deviation(), we add a new site on the input edge near the
    //      middle of the snapped edge.  This causes the snapped edge to split
    //      into two pieces, so that it follows the input edge more closely.
    //
    //  (2) If a snapped edge is closer than min_edge_vertex_separation() to any
    //      nearby site (the "site to avoid") or passes on the wrong side of it
    //      relative to the input edge, then we add a new site (the "separation
    //      site") along the input edge near the site to avoid.  This causes the
    //      snapped edge to follow the input edge more closely, so that it is
    //      guaranteed to pass on the correct side of the site to avoid with a
    //      separation of at least the required distance.
    //
    // We check these conditions by snapping all the input edges to a chain of
    // Voronoi sites and then testing each edge in the chain.  If a site needs to
    // be added, we mark all nearby edges for re-snapping.
    private void AddExtraSites(MutableS2ShapeIndex input_edge_index)
    {
        // Note that we could save some work in AddSnappedEdges() by saving the
        // snapped edge chains in a vector, but currently this is not worthwhile
        // since SnapEdge() accounts for less than 5% of the runtime.

        // Note that we intentionally use dense_hash_set rather than flat_hash_set
        // in order to ensure that iteration is deterministic when debugging.
        Dictionary<InputEdgeId, InputEdgeId> edges_to_resnap = new(16 /*expected_max_elements*/);
        //edges_to_resnap.set_empty_key(-1);
        //edges_to_resnap.set_deleted_key(-2);

        List<SiteId> chain = new();  // Temporary storage.
        int num_edges_after_snapping = 0;

        // CheckEdge() defines the body of the loops below.
        var CheckEdge = (InputEdgeId e) =>
        {
            if (!tracker_.Ok()) return false;
            SnapEdge(e, chain);
            edges_to_resnap.Remove(e);
            num_edges_after_snapping += chain.Count;
            MaybeAddExtraSites(e, chain, input_edge_index, edges_to_resnap);
            return true;
        };

        // The first pass is different because we snap every edge.  In the following
        // passes we only snap edges that are near the extra sites that were added.
        Debug.WriteLine($"Before pass 0: sites={sites_.Count}");
        for (InputEdgeId e = 0; e < input_edges_.Count; ++e)
        {
            if (!CheckEdge(e)) return;
        }
        Debug.WriteLine($"Pass 0: edges snapped={input_edges_.Count}, output edges={num_edges_after_snapping}, sites={sites_.Count}");

        for (int num_passes = 1; edges_to_resnap.Any(); ++num_passes)
        {
            var edges_to_snap = edges_to_resnap;
            edges_to_resnap.Clear();
            num_edges_after_snapping = 0;
            foreach (InputEdgeId e in edges_to_snap.Values)
            {
                if (!CheckEdge(e)) return;
            }
            Debug.WriteLine($"Pass {num_passes}: edges snapped={edges_to_snap.Count}, output edges={num_edges_after_snapping}, sites={sites_.Count}");
        }
    }

    private void MaybeAddExtraSites(
        int edge_id, List<int> chain, 
        MutableS2ShapeIndex input_edge_index,
        Dictionary<InputEdgeId, InputEdgeId> edges_to_resnap)
    {
        // If the memory tracker has a periodic callback function, tally an amount
        // of memory proportional to the work being done so that the caller has an
        // opportunity to cancel the operation if necessary.
        if (!tracker_.TallyTemp(chain.Count * Marshal.SizeOf(chain[0]))) return;

        // If the input includes NaN vertices, snapping can produce an empty chain.
        if (!chain.Any()) return;

        // The snapped edge chain is always a subsequence of the nearby sites
        // (edge_sites_), so we walk through the two arrays in parallel looking for
        // sites that weren't snapped.  These are the "sites to avoid".  We also keep
        // track of the current snapped edge, since it is the only edge that can be
        // too close or pass on the wrong side of a site to avoid.  Vertices beyond
        // the chain endpoints in either direction can be ignored because only the
        // interiors of chain edges can be too close to a site to avoid.
        var edge = input_edges_[edge_id];
        var a0 = input_vertices_[edge.Item1];
        var a1 = input_vertices_[edge.Item2];
        var nearby_sites = edge_sites_[edge_id];
        for (int i = 0, j = 0; j < nearby_sites.Count; ++j)
        {
            SiteId id = nearby_sites[j];
            if (id == chain[i])
            {
                // This site is a vertex of the snapped edge chain.
                if (++i == chain.Count)
                {
                    break;  // Sites beyond the end of the snapped chain can be ignored.
                }


                // Check whether this snapped edge deviates too far from its original
                // position.  If so, we split the edge by adding an extra site.
                var v0 = sites_[chain[i - 1]];
                var v1 = sites_[chain[i]];
                if (new S1ChordAngle(v0, v1) < min_edge_length_to_split_ca_) continue;
                if (!S2.IsEdgeBNearEdgeA(a0, a1, v0, v1, max_edge_deviation_))
                {
                    // Add a new site on the input edge, positioned so that it splits the
                    // snapped edge into two approximately equal pieces.  Then we find all
                    // the edges near the new site (including this one) and add them to
                    // the snap queue.
                    //
                    // Note that with large snap radii, it is possible that the snapped
                    // edge wraps around the sphere the "wrong way".  To handle this we
                    // find the preferred split location by projecting both endpoints onto
                    // the input edge and taking their midpoint.
                    var mid = (S2.Project(v0, a0, a1) +
                                   S2.Project(v1, a0, a1)).Normalize();
                    var new_site = GetSeparationSite(mid, v0, v1, edge_id);
                    AddExtraSite(new_site, input_edge_index, edges_to_resnap);

                    // In the case where the edge wrapped around the sphere the "wrong
                    // way", it is not safe to continue checking this edge.  It will be
                    // marked for resnapping and we will come back to it in the next pass.
                    return;
                }
            }
            else
            {
                // This site is near the input edge but is not part of the snapped chain.
                if (i == 0)
                {
                    continue;  // Sites before the start of the chain can be ignored.
                }
                // We need to ensure that non-forced sites are separated by at least
                // min_edge_vertex_separation() from the snapped chain.  This happens
                // automatically as part of the algorithm except where there are portions
                // of the input edge that are not within edge_snap_radius() of any site.
                // These portions of the original edge are called "coverage gaps".
                // Therefore if we find that a site to avoid that is too close to the
                // snapped edge chain, we can fix the problem by adding a new site (the
                // "separation site") in the corresponding coverage gap located as closely
                // as possible to the site to avoid.  This technique is is guaranteed to
                // produce the required minimum separation, and the entire process of
                // adding separation sites is guaranteed to terminate.
                var site_to_avoid = sites_[id];
                var v0 = sites_[chain[i - 1]];
                var v1 = sites_[chain[i]];
                bool add_separation_site = false;
                if (!IsForced(id) &&
                    min_edge_site_separation_ca_ > S1ChordAngle.Zero &&
                    S2Pred.CompareEdgeDistance(
                        site_to_avoid, v0, v1, min_edge_site_separation_ca_) < 0)
                {
                    add_separation_site = true;
                }
                // Similarly, we also add a separation site whenever a snapped edge passes
                // on the wrong side of a site to avoid.  Normally we don't need to worry
                // about this, since if an edge passes on the wrong side of a nearby site
                // then it is also too close to it.  However if the snap radius is very
                // small and intersection_tolerance() is non-zero then we need to check
                // this condition explicitly (see the "check_all_site_crossings_" flag for
                // details).  We also need to check this condition explicitly for forced
                // vertices.  Again, we can solve this problem by adding a "separation
                // site" in the corresponding coverage gap located as closely as possible
                // to the site to avoid.
                //
                // It is possible to show that when all points are projected onto the
                // great circle through (a0, a1), no improper crossing occurs unless the
                // the site to avoid is located between a0 and a1, and also between v0
                // and v1.  TODO(ericv): Verify whether all these checks are necessary.
                if (!add_separation_site &&
                    (IsForced(id) || check_all_site_crossings_) &&
                    (S2Pred.Sign(a0, a1, site_to_avoid) !=
                     S2Pred.Sign(v0, v1, site_to_avoid)) &&
                    S2Pred.CompareEdgeDirections(a0, a1, a0, site_to_avoid) > 0 &&
                    S2Pred.CompareEdgeDirections(a0, a1, site_to_avoid, a1) > 0 &&
                    S2Pred.CompareEdgeDirections(a0, a1, v0, site_to_avoid) > 0 &&
                    S2Pred.CompareEdgeDirections(a0, a1, site_to_avoid, v1) > 0)
                {
                    add_separation_site = true;
                }
                if (add_separation_site)
                {
                    // We add a new site (the separation site) in the coverage gap along the
                    // input edge, located as closely as possible to the site to avoid.
                    // Then we find all the edges near the new site (including this one) and
                    // add them to the snap queue.
                    var new_site = GetSeparationSite(site_to_avoid, v0, v1, edge_id);
                    Debug.Assert(site_to_avoid != new_site);
                    AddExtraSite(new_site, input_edge_index, edges_to_resnap);

                    // Skip the remaining sites near this chain edge, and then continue
                    // scanning this chain.  Note that this is safe even though the call
                    // to AddExtraSite() above added a new site to "nearby_sites".
                    for (; nearby_sites[j + 1] != chain[i]; ++j) { }
                }
            }
        }
    }

    // Adds a new site, then updates "edge_sites"_ for all edges near the new site
    // and adds them to "edges_to_resnap" for resnapping.
    private void AddExtraSite(S2Point new_site, 
        MutableS2ShapeIndex input_edge_index,
        Dictionary<InputEdgeId, InputEdgeId> edges_to_resnap)
    {
        if (sites_.Any()) Debug.Assert(new_site != sites_.Last());
        if (!tracker_.AddSpace(sites_, 1)) return;
        var new_site_id = sites_.Count;
        sites_.Add(new_site);

        // Find all edges whose distance is <= edge_site_query_radius_ca_.
        S2ClosestEdgeQuery.Options? options = new();
        options.ConservativeMaxDistance = edge_site_query_radius_ca_;
        options.IncludeInteriors = (false);

        if (!input_edge_index.IsFresh()) input_edge_index.ForceBuild();
        if (!tracker_.Ok()) return;

        // Memory used by S2ClosestEdgeQuery is not tracked, but it is temporary,
        // typically insignificant, and does not affect the high water mark.
        S2ClosestEdgeQuery query = new(input_edge_index, options);
        S2ClosestEdgeQuery.PointTarget target = new(new_site);
        foreach (var result in query.FindClosestEdges(target))
        {
            var e = result.EdgeId;
            var v0 = input_vertices_[input_edges_[e].Item1];
            InsertSiteByDistance(new_site_id, v0, edge_sites_[e]);
            edges_to_resnap.Add(e, e);
        }
    }

    private S2Point GetSeparationSite(S2Point site_to_avoid, S2Point v0, S2Point v1, int input_edge_id)
    {
        // Define the "coverage disc" of a site S to be the disc centered at S with
        // radius "snap_radius".  Similarly, define the "coverage interval" of S for
        // an edge XY to be the intersection of XY with the coverage disc of S.  The
        // SnapFunction implementations guarantee that the only way that a snapped
        // edge can be closer than min_edge_vertex_separation() to a non-snapped
        // site (i.e., site_to_avoid) if is there is a gap in the coverage of XY
        // near this site.  We can fix this problem simply by adding a new site to
        // fill this gap, located as closely as possible to the site to avoid.
        //
        // To calculate the coverage gap, we look at the two snapped sites on
        // either side of site_to_avoid, and find the endpoints of their coverage
        // intervals.  The we place a new site in the gap, located as closely as
        // possible to the site to avoid.  Note that the new site may move when it
        // is snapped by the snap_function, but it is guaranteed not to move by
        // more than snap_radius and therefore its coverage interval will still
        // intersect the gap.
        var edge = input_edges_[input_edge_id];
        var x = input_vertices_[edge.Item1];
        var y = input_vertices_[edge.Item2];
        var xy_dir = y - x;
        var n = S2.RobustCrossProd(x, y);
        var new_site = S2.Project(site_to_avoid, x, y, n);
        var gap_min = GetCoverageEndpoint(v0/*, x, y*/, n);
        var gap_max = GetCoverageEndpoint(v1/*, y, x*/, -n);
        if ((new_site - gap_min).DotProd(xy_dir) < 0)
        {
            new_site = gap_min;
        }
        else if ((gap_max - new_site).DotProd(xy_dir) < 0)
        {
            new_site = gap_max;
        }
        new_site = SnapSite(new_site);
        Debug.Assert(v0 != new_site);
        Debug.Assert(v1 != new_site);
        return new_site;
    }

    // Given a site P and an edge XY with normal N, intersect XY with the disc of
    // radius snap_radius() around P, and return the intersection point that is
    // further along the edge XY toward Y.
    private S2Point GetCoverageEndpoint(S2Point p/*, S2Point x, S2Point y*/, S2Point n)
    {
        // Consider the plane perpendicular to P that cuts off a spherical cap of
        // radius snap_radius().  This plane intersects the plane through the edge
        // XY (perpendicular to N) along a line, and that line intersects the unit
        // sphere at two points Q and R, and we want to return the point R that is
        // further along the edge XY toward Y.
        //
        // Let M be the midpoint of QR.  This is the point along QR that is closest
        // to P.  We can now express R as the sum of two perpendicular vectors OM
        // and MR in the plane XY.  Vector MR is in the direction N x P, while
        // vector OM is in the direction (N x P) x N, where N = X x Y.
        //
        // The length of OM can be found using the Pythagorean theorem on triangle
        // OPM, and the length of MR can be found using the Pythagorean theorem on
        // triangle OMR.
        //
        // In the calculations below, we save some work by scaling all the vectors
        // by n.CrossProd(p).Norm2, and normalizing at the end.
        var n2 = n.Norm2();
        var nDp = n.DotProd(p);
        var nXp = n.CrossProd(p);
        var nXpXn = n2 * p - nDp * n;
        var om = Math.Sqrt(1 - edge_snap_radius_sin2_) * nXpXn;
        var mr2 = edge_snap_radius_sin2_ * n2 - nDp * nDp;

        // MR is constructed so that it points toward Y (rather than X).
        var mr = Math.Sqrt(Math.Max(0.0, mr2)) * nXp;
        return (om + mr).Normalize();
    }

    private void SnapEdge(int e, List<int> chain)
    {
        chain.Clear();
        var edge = input_edges_[e];
        if (!snapping_needed_)
        {
            // Note that the input vertices have been renumbered such that
            // InputVertexId and SiteId are the same (see ChooseAllVerticesAsSites).
            chain.Add(edge.Item1);
            chain.Add(edge.Item2);
            return;
        }

        var x = input_vertices_[edge.Item1];
        var y = input_vertices_[edge.Item2];

        // Optimization: if there is only one nearby site, return.
        // Optimization: if there are exactly two nearby sites, and one is close
        // enough to each vertex, then return.

        // Now iterate through the sites.  We keep track of the sequence of sites
        // that are visited.
        var candidates = edge_sites_[e];
        foreach (var site_id in candidates)
        {
            var c = sites_[site_id];
            // Skip any sites that are too far away.  (There will be some of these,
            // because we also keep track of "sites to avoid".)  Note that some sites
            // may be close enough to the line containing the edge, but not to the
            // edge itself, so we can just use the dot product with the edge normal.
            if (S2Pred.CompareEdgeDistance(c, x, y, edge_snap_radius_ca_) > 0)
            {
                continue;
            }
            // Check whether the new site C excludes the previous site B.  If so,
            // repeat with the previous site, and so on.
            var add_site_c = true;
            for (; chain.Any(); chain.RemoveAt(chain.Count - 1))
            {
                var b = sites_[chain.Last()];

                // First, check whether B and C are so far apart that their clipped
                // Voronoi regions can't intersect.
                var bc = new S1ChordAngle(b, c);
                if (bc >= max_adjacent_site_separation_ca_) break;

                // Otherwise, we want to check whether site C prevents the Voronoi
                // region of B from intersecting XY, or vice versa.  This can be
                // determined by computing the "coverage interval" (the segment of XY
                // intersected by the coverage disc of radius snap_radius) for each
                // site.  If the coverage interval of one site contains the coverage
                // interval of the other, then the contained site can be excluded.
                var result = S2Pred.GetVoronoiSiteExclusion(
                    b, c, x, y, edge_snap_radius_ca_);
                if (result == S2Pred.Excluded.FIRST) continue;  // Site B excluded by C
                if (result == S2Pred.Excluded.SECOND)
                {
                    add_site_c = false;  // Site C is excluded by B.
                    break;
                }
                Debug.Assert(S2Pred.Excluded.NEITHER == result);

                // Otherwise check whether the previous site A is close enough to B and
                // C that it might further clip the Voronoi region of B.
                if (chain.Count < 2) break;
                var a = sites_[chain[^2]];
                S1ChordAngle ac = new(a, c);
                if (ac >= max_adjacent_site_separation_ca_) break;

                // If triangles ABC and XYB have the same orientation, the circumcenter
                // Z of ABC is guaranteed to be on the same side of XY as B.
                var xyb = S2Pred.Sign(x, y, b);
                if (S2Pred.Sign(a, b, c) == xyb)
                {
                    break;  // The circumcenter is on the same side as B but further away.
                }
                // Other possible optimizations:
                //  - if AB > max_adjacent_site_separation_ca_ then keep B.
                //  - if d(B, XY) < 0.5 * Math.Min(AB, BC) then keep B.

                // If the circumcenter of ABC is on the same side of XY as B, then B is
                // excluded by A and C combined.  Otherwise B is needed and we can exit.
                if (S2Pred.EdgeCircumcenterSign(x, y, a, b, c) != xyb) break;
            }
            if (add_site_c)
            {
                chain.Add(site_id);
            }
        }
        Debug.Assert(chain.Any());

        WriteS2Polyline(edge, chain);
    }

    [Conditional("s2builder_verbose")]
    private static void WriteS2Polyline(InputEdge edge, List<int> chain)
    {
        System.Text.StringBuilder sb = new($"({edge.Item1},{edge.Item2}): ");
        foreach (var id in chain) sb.Append(id + " ");
        Debug.WriteLine(sb.ToString());
    }

    private const int kMinLayersForVertexFiltering = 10;

    private void BuildLayers()
    {
        if (!tracker_.Ok()) return;

        // Each output edge has an "input edge id set id" (an int) representing
        // the set of input edge ids that were snapped to this edge.  The actual
        // ints can be retrieved using "input_edge_id_set_lexicon".
        var layer_edges = new List<List<OutputEdge>>();
        var layer_input_edge_ids = new List<List<int>>();
        var input_edge_id_set_lexicon = new IdSetLexicon();
        List<List<S2Point>> layer_vertices = new();
        BuildLayerEdges(layer_edges, layer_input_edge_ids, input_edge_id_set_lexicon);

        try
        {
            // If there are a large number of layers, then we build a minimal subset of
            // vertices for each layer.  This ensures that layer types that iterate over
            // vertices will run in time proportional to the size of that layer rather
            // than the size of all layers combined.
            if (layers_.Count >= kMinLayersForVertexFiltering)
            {
                // Disable vertex filtering if it is disallowed by any layer.  (This could
                // be optimized, but in current applications either all layers allow
                // filtering or none of them do.)
                var allow_vertex_filtering = true;
                foreach (var options in layer_options_)
                {
                    allow_vertex_filtering &= options.AllowVertexFiltering;
                }
                if (allow_vertex_filtering)
                {
                    // Track the temporary memory used by FilterVertices().  Note that
                    // although vertex filtering can increase the number of vertices stored
                    // (i.e., if the same vertex is referred to by multiple layers), it
                    // never increases storage quadratically because there can be at most
                    // two filtered vertices per edge.
                    if (!tracker_.TallyFilterVertices(sites_.Count, layer_edges)) return;
                    try
                    {
                        layer_vertices.Capacity = layers_.Count;
                        List<VertexId> filter_tmp = new();  // Temporary used by FilterVertices.
                        for (var i = 0; i < layers_.Count; ++i)
                        {
                            layer_vertices[i] = Graph.FilterVertices(sites_, layer_edges[i], filter_tmp);
                            if (!tracker_.Tally(layer_vertices[i])) return;
                        }
                        tracker_.Clear(sites_);  // Release memory
                    }
                    finally
                    {
                        tracker_.DoneFilterVertices();
                    }
                }
            }
            if (!tracker_.Ok()) return;

            for (var i = 0; i < layers_.Count; ++i)
            {
                var vertices = (!layer_vertices.Any() ? sites_ : layer_vertices[i]);
                var graph = new Graph(layer_options_[i], vertices, layer_edges[i],
                layer_input_edge_ids[i], input_edge_id_set_lexicon,
                label_set_ids_, label_set_lexicon_,
                layer_is_full_polygon_predicates_[i]);
                layers_[i].Build(graph, out error_);
                // Don't free the layer data until all layers have been built, in order to
                // support building multiple layers at once (e.g. ClosedSetNormalizer).
            }
        }
        finally
        {
            for (int i = 0; i < layers_.Count; ++i)
            {
                tracker_.Untally(layer_edges[i]);
                tracker_.Untally(layer_input_edge_ids[i]);
                if (layer_vertices.Any()) tracker_.Untally(layer_vertices[i]);
            }
        }
    }

    // Snaps and possibly simplifies the edges for each layer, populating the
    // given output arguments.  The resulting edges can be used to construct an
    // S2Builder.Graph directly (no further processing is necessary).
    //
    // This method is not "const" because Graph.ProcessEdges can modify
    // layer_options_ in some cases (changing undirected edges to directed ones).
    private void BuildLayerEdges(List<List<OutputEdge>> layer_edges, List<List<int>> layer_input_edge_ids, IdSetLexicon input_edge_id_set_lexicon)
    {
        // Edge chains are simplified only when a non-zero snap radius is specified.
        // If so, we build a map from each site to the set of input vertices that
        // snapped to that site.  (Note that site_vertices is relatively small and
        // that its memory tracking is deferred until TallySimplifyEdgeChains.)
        List<List<int>> site_vertices = new();
        var simplify = snapping_needed_ && Options_.SimplifyEdgeChains;
        if (simplify) site_vertices.Capacity = sites_.Count;

        while (layer_edges.Count < layers_.Count) layer_edges.Add(new());
        while (layer_input_edge_ids.Count < layers_.Count) layer_input_edge_ids.Add(new List<int>());
        for (var i = 0; i < layers_.Count; ++i)
        {
            AddSnappedEdges(layer_begins_[i], layer_begins_[i + 1], layer_options_[i],
                            layer_edges[i], layer_input_edge_ids[i]/*,
                            input_edge_id_set_lexicon*/, site_vertices);
        }

        // We simplify edge chains before processing the per-layer GraphOptions
        // because simplification can create duplicate edges and/or sibling edge
        // pairs which may need to be removed.
        if (simplify)
        {
            SimplifyEdgeChains(site_vertices, layer_edges, layer_input_edge_ids,
                               input_edge_id_set_lexicon);
            site_vertices.Clear();
        }

        // At this point we have no further need for nearby site data, so we clear
        // it to save space.  We keep input_vertices_ and input_edges_ so that
        // S2Builder::Layer implementations can access them if desired.  (This is
        // useful for determining how snapping has changed the input geometry.)
        tracker_.ClearEdgeSites(edge_sites_);
        for (var i = 0; i < layers_.Count; ++i)
        {
            // The errors generated by ProcessEdges are really warnings, so we simply
            // record them and continue.
            Graph.ProcessEdges(layer_options_[i], layer_edges[i],
                                layer_input_edge_ids[i],
                                input_edge_id_set_lexicon, out error_, tracker_);

            if (!tracker_.Ok()) return;
        }
    }

    // Snaps all the input edges for a given layer, populating the given output
    // arguments.  If (*site_vertices) is non-empty then it is updated so that
    // (*site_vertices)[site] contains a list of all input vertices that were
    // snapped to that site.
    private void AddSnappedEdges(int begin, int end, GraphOptions options, List<OutputEdge> edges, List<int> input_edge_ids/*, IdSetLexicon input_edge_id_set_lexicon*/, List<List<int>> site_vertices)
    {
        var discard_degenerate_edges = options.DegenerateEdges_ == GraphOptions.DegenerateEdges.DISCARD;
        var chain = new List<int>();
        for (var e = begin; e < end; ++e)
        {
            var id = IdSetLexicon.AddSingleton(e);
            SnapEdge(e, chain);
            int num_snapped_edges = Math.Max(1, chain.Count - 1);
            if (options.EdgeType_ == EdgeType.UNDIRECTED) num_snapped_edges *= 2;
            if (!tracker_.AddSpace(edges, num_snapped_edges)) return;
            if (!tracker_.AddSpace(input_edge_ids, num_snapped_edges)) return;
            MaybeAddInputVertex(input_edges_[e].Item1, chain[0], site_vertices);
            if (chain.Count == 1)
            {
                if (discard_degenerate_edges) continue;
                AddSnappedEdge(chain[0], chain[0], id, options.EdgeType_,
                               edges, input_edge_ids);
            }
            else
            {
                MaybeAddInputVertex(input_edges_[e].Item2, chain.Last(), site_vertices);
                for (var i = 1; i < chain.Count; ++i)
                {
                    AddSnappedEdge(chain[i - 1], chain[i], id, options.EdgeType_,
                                   edges, input_edge_ids);
                }
            }
        }
        DumpEdges(edges, sites_);
    }

    [Conditional("s2builder_verbose")]
    private static void DumpEdges(List<OutputEdge> edges, List<S2Point> vertices)
    {
        foreach (var e in edges)
        {
            var v = new S2Point[2] { vertices[e.ShapeId], vertices[e.EdgeId] };
            Debug.WriteLine($"S2Polyline: {v.ToDebugString()}({e.ShapeId},{e.EdgeId})");
        }
    }

    // If "site_vertices" is non-empty, ensures that (*site_vertices)[id] contains
    // "v".  Duplicate entries are allowed.  The purpose of this function is to
    // build a map so that SimplifyEdgeChains() can quickly find all the input
    // vertices that snapped to a particular site.
    private static void MaybeAddInputVertex(int v, int id, List<List<int>> site_vertices)
    {
        if (!site_vertices.Any()) return;

        // Optimization: check if we just added this vertex.  This is worthwhile
        // because the input edges usually form a continuous chain, i.e. the
        // destination of one edge is the same as the source of the next edge.
        var vertices = site_vertices[id];
        if (!vertices.Any() || vertices.Last() != v)
        {
            // Memory tracking is deferred until SimplifyEdgeChains.
            vertices.Add(v);
        }
    }

    // Adds the given edge to "edges" and "input_edge_ids".  If undirected edges
    // are being used, also adds an edge in the opposite direction.
    private static void AddSnappedEdge(int src, int dst, int id, EdgeType edge_type, List<OutputEdge> edges, List<int> input_edge_ids)
    {
        edges.Add(new(src, dst));
        input_edge_ids.Add(id);
        if (edge_type == EdgeType.UNDIRECTED)
        {
            edges.Add(new(dst, src));
            // Automatically created edges do not have input edge ids or labels.  This
            // can be used to distinguish the original direction of the undirected edge.
            input_edge_ids.Add(IdSetLexicon.kEmptySetId);
        }
    }

    // Simplifies edge chains, updating its input/output arguments as necessary.
    private void SimplifyEdgeChains(List<List<int>> site_vertices, List<List<OutputEdge>> layer_edges, List<List<int>> layer_input_edge_ids, IdSetLexicon input_edge_id_set_lexicon)
    {
        if (!layers_.Any()) return;
        if (!tracker_.TallySimplifyEdgeChains(site_vertices, layer_edges)) return;

        // Merge the edges from all layers (in order to build a single graph).
        List<OutputEdge> merged_edges = new();
        List<int> merged_input_edge_ids = new();
        List<int> merged_edge_layers = new();
        MergeLayerEdges(layer_edges, layer_input_edge_ids,
                        merged_edges, merged_input_edge_ids, merged_edge_layers);

        // The following fields will be reconstructed by EdgeChainSimplifier.
        foreach (var edges in layer_edges) edges.Clear();
        foreach (var input_edge_ids in layer_input_edge_ids) input_edge_ids.Clear();

        // The graph options are irrelevant for edge chain simplification, but we
        // try to set them appropriately anyway.
        var graph_options = new GraphOptions(EdgeType.DIRECTED,
                                              GraphOptions.DegenerateEdges.KEEP,
                                              GraphOptions.DuplicateEdges.KEEP,
                                              GraphOptions.SiblingPairs.KEEP);
        var graph = new Graph(graph_options, sites_, merged_edges, merged_input_edge_ids,
          input_edge_id_set_lexicon, null, null, null);
        var simplifier = new EdgeChainSimplifier(
            this, graph, merged_edge_layers, site_vertices,
            layer_edges, layer_input_edge_ids, input_edge_id_set_lexicon);
        simplifier.Run();
    }

    // Merges the edges from all layers and sorts them in lexicographic order so
    // that we can construct a single graph.  The sort is stable, which means that
    // any duplicate edges within each layer will still be sorted by int.
    private static void MergeLayerEdges(List<List<OutputEdge>> layer_edges, List<List<int>> layer_input_edge_ids, List<OutputEdge> edges, List<int> input_edge_ids, List<int> edge_layers)
    {
        List<LayerEdgeId> order = new();
        for (var i = 0; i < layer_edges.Count; ++i)
        {
            for (var e = 0; e < layer_edges[i].Count; ++e)
            {
                order.Add(new(i, e));
            }
        }
        order.Sort((LayerEdgeId ai, LayerEdgeId bi) =>
        {
            return StableLessThan(layer_edges[ai.Item1][ai.Item2],
                                  layer_edges[bi.Item1][bi.Item2], ai, bi);
        });
        edges.Capacity = order.Count;
        input_edge_ids.Capacity = order.Count;
        edge_layers.Capacity = order.Count;
        foreach (var id in order)
        {
            edges.Add(layer_edges[id.Item1][id.Item2]);
            input_edge_ids.Add(layer_input_edge_ids[id.Item1][id.Item2]);
            edge_layers.Add(id.Item1);
        }
    }

    // A comparison function that allows stable sorting with sort (which is
    // fast but not stable).  It breaks ties between equal edges by comparing
    // their LayerEdgeIds.
    private static int StableLessThan(OutputEdge a, OutputEdge b, LayerEdgeId ai, LayerEdgeId bi)
    {
        // The compiler doesn't optimize this as well as it should:
        //   return (a, ai) < (b, bi);
        if (a.ShapeId.CompareTo(b.ShapeId) != 0) return a.ShapeId.CompareTo(b.ShapeId);
        if (a.EdgeId.CompareTo(b.EdgeId) != 0) return a.EdgeId.CompareTo(b.EdgeId);
        if (b.EdgeId.CompareTo(a.EdgeId) < 0) return 1;
        return ai.CompareTo(bi);  // Stable sort.
    }

    // Helper functions for computing error bounds:
    private static S1ChordAngle RoundUp(S1Angle a)
    {
        var ca = new S1ChordAngle(a);
        return ca.PlusError(ca.S1AngleConstructorMaxError);
    }
    private static S1ChordAngle AddPointToPointError(S1ChordAngle ca) => ca.PlusError(ca.GetS2PointConstructorMaxError());
    private static S1ChordAngle AddPointToEdgeError(S1ChordAngle ca) => ca.PlusError(S2.GetUpdateMinDistanceMaxError(ca));

    //////////// Parameters /////////////

    // S2Builder options.
    public Options Options_ { get; set; }

    // The maximum distance (inclusive) that a vertex can move when snapped,
    // equal to S1ChordAngle(options_.snap_function().snap_radius()).
    private readonly S1ChordAngle site_snap_radius_ca_;

    // The maximum distance (inclusive) that an edge can move when snapping to a
    // snap site.  It can be slightly larger than the site snap radius when
    // edges are being split at crossings.
    private readonly S1ChordAngle edge_snap_radius_ca_;

    // True if we need to check that snapping has not changed the input topology
    // around any vertex (i.e. Voronoi site).  Normally this is only necessary for
    // forced vertices, but if the snap radius is very small (e.g., zero) and
    // split_crossing_edges() is true then we need to do this for all vertices.
    // In all other situations, any snapped edge that crosses a vertex will also
    // be closer than min_edge_vertex_separation() to that vertex, which will
    // cause us to add a separation site anyway.
    private readonly bool check_all_site_crossings_;

    private readonly S1Angle max_edge_deviation_;
    private readonly S1ChordAngle edge_site_query_radius_ca_;
    private readonly S1ChordAngle min_edge_length_to_split_ca_;

    private readonly S1Angle min_site_separation_;
    private readonly S1ChordAngle min_site_separation_ca_;
    private readonly S1ChordAngle min_edge_site_separation_ca_;
    private readonly S1ChordAngle min_edge_site_separation_ca_limit_;

    private readonly S1ChordAngle max_adjacent_site_separation_ca_;

    // The squared sine of the edge snap radius.  This is equivalent to the snap
    // radius (squared) for distances measured through the interior of the
    // sphere to the plane containing an edge.  This value is used only when
    // interpolating new points along edges (see GetSeparationSite).
    private readonly double edge_snap_radius_sin2_;

    // A copy of the argument to Build().
    private S2Error error_;

    // True if snapping was requested.  This is true if either snap_radius() is
    // positive, or split_crossing_edges() is true (which implicitly requests
    // snapping to ensure that both crossing edges are snapped to the
    // intersection point).
    private readonly bool snapping_requested_;

    // Initially false, and set to true when it is discovered that at least one
    // input vertex or edge does not meet the output guarantees (e.g., that
    // vertices are separated by at least snap_function.min_vertex_separation).
    private bool snapping_needed_;

    //////////// Input Data /////////////

    // A flag indicating whether label_set_ has been modified since the last
    // time label_set_id_ was computed.
    private bool label_set_modified_;

    private List<S2Point> input_vertices_ = new();
    private readonly List<InputEdge> input_edges_ = new();

    private readonly List<Layer> layers_ = new();
    private readonly List<GraphOptions> layer_options_ = new();
    private readonly List<int> layer_begins_ = new();
    private readonly List<IsFullPolygonPredicate> layer_is_full_polygon_predicates_ = new();

    // Each input edge has "label set id" (an int) representing the set of
    // labels attached to that edge.  This vector is populated only if at least
    // one label is used.
    private readonly List<int> label_set_ids_ = new();
    private readonly IdSetLexicon label_set_lexicon_ = new();

    // The current set of labels (represented as a stack).
    private readonly List<int> label_set_ = new();

    // The LabelSetId corresponding to the current label set, computed on demand
    // (by adding it to label_set_lexicon()).
    private int label_set_id_;

    ////////////// Data for Snapping and Simplifying //////////////

    // The number of sites specified using ForceVertex().  These sites are
    // always at the beginning of the sites_ vector.
    private int num_forced_sites_;

    // The set of snapped vertex locations ("sites").
    private readonly List<S2Point> sites_ = new();

    // A map from each input edge to the set of sites "nearby" that edge,
    // defined as the set of sites that are candidates for snapping and/or
    // avoidance.  Note that compact_array will inline up to two sites, which
    // usually takes care of the vast majority of edges.  Sites are kept sorted
    // by increasing distance from the origin of the input edge.
    //
    // Once snapping is finished, this field is discarded unless edge chain
    // simplification was requested, in which case instead the sites are
    // filtered by removing the ones that each edge was snapped to, leaving only
    // the "sites to avoid" (needed for simplification).
    private readonly List<List<int>> edge_sites_ = new();

    // An object to track the memory usage of this class.
    private readonly MemoryTracker tracker_ = new();

    // Indicates whether the input edges are undirected.  Typically this is
    // specified for each output layer (e.g., s2builderutil.S2PolygonLayer).
    //
    // Directed edges are preferred, since otherwise the output is ambiguous.
    // For example, output polygons may be the *inverse* of the intended result
    // (e.g., a polygon intended to represent the world's oceans may instead
    // represent the world's land masses).  Directed edges are also somewhat
    // more efficient.
    //
    // However even with undirected edges, most S2Builder layer types try to
    // preserve the input edge direction whenever possible.  Generally, edges
    // are reversed only when it would yield a simpler output.  For example,
    // S2PolygonLayer assumes that polygons created from undirected edges should
    // cover at most half of the sphere.  Similarly, S2PolylineVectorLayer
    // assembles edges into as few polylines as possible, even if this means
    // reversing some of the "undirected" input edges.
    //
    // For shapes with interiors, directed edges should be oriented so that the
    // interior is to the left of all edges.  This means that for a polygon with
    // holes, the outer loops ("shells") should be directed counter-clockwise
    // while the inner loops ("holes") should be directed clockwise.  Note that
    // S2Builder.AddPolygon() follows this convention automatically.
    public enum EdgeType : byte { DIRECTED, UNDIRECTED };

    public class Options
    {
        public Options() => SnapFunction = new IdentitySnapFunction(S1Angle.Zero);

        // Convenience constructor that calls set_snap_function().
        public Options(SnapFunction snapFunction) => SnapFunction = (SnapFunction)snapFunction.CustomClone();

        public Options(Options options)
        {
            SnapFunction = (SnapFunction)options.SnapFunction.CustomClone();
            SplitCrossingEdges = options.SplitCrossingEdges;
            intersection_tolerance_ = options.intersection_tolerance_;
            SimplifyEdgeChains = options.SimplifyEdgeChains;
            Idempotent = options.Idempotent;
            MemoryTracker = options.MemoryTracker;
        }

        // The maximum distance from snapped edge vertices to the original edge.
        // This is the same as snap_function().snap_radius() except when
        // split_crossing_edges() is true (see below), in which case the edge snap
        // radius is increased by S2::kIntersectionError.
        public S1Angle EdgeSnapRadius() => SnapFunction.SnapRadius + IntersectionTolerance;

        // Sets the desired snap function.  The snap function is copied
        // internally, so you can safely pass a temporary object.
        //
        // Note that if your input data includes vertices that were created using
        // S2EdgeCrossings.GetIntersection(), then you should use a "snap_radius" of
        // at least S2.kIntersectionSnapRadius, e.g. by calling
        //
        //  options.set_snap_function(s2builderutil.IdentitySnapFunction(
        //      S2.kIntersectionSnapRadius));
        //
        // DEFAULT: s2builderutil.IdentitySnapFunction(S1Angle.Zero())
        // [This does no snapping and preserves all input vertices exactly.]
        public SnapFunction SnapFunction { get; set; }

        // If true, then detect all pairs of crossing edges and eliminate them by
        // adding a new vertex at their intersection point.  See also the
        // AddIntersection() method which allows intersection points to be added
        // selectively.
        //
        // When this option if true, intersection_tolerance() is automatically set
        // to a minimum of S2::kIntersectionError (see intersection_tolerance()
        // for why this is necessary).  Note that this means that edges can move
        // by up to S2::kIntersectionError even when the specified snap radius is
        // zero.  The exact distance that edges can move is always given by
        // max_edge_deviation() defined above.
        //
        // Undirected edges should always be used when the output is a polygon,
        // since splitting a directed loop at a self-intersection converts it into
        // two loops that don't define a consistent interior according to the
        // "interior is on the left" rule.  (On the other hand, it is fine to use
        // directed edges when defining a polygon *mesh* because in that case the
        // input consists of sibling edge pairs.)
        //
        // Self-intersections can also arise when importing data from a 2D
        // projection.  You can minimize this problem by subdividing the input
        // edges so that the S2 edges (which are geodesics) stay close to the
        // original projected edges (which are curves on the sphere).  This can
        // be done using S2EdgeTessellator, for example.
        //
        // DEFAULT: false
        public bool SplitCrossingEdges { get; set; } = false;

        // Specifes the maximum allowable distance between a vertex added by
        // AddIntersection() and the edge(s) that it is intended to snap to.  This
        // method must be called before AddIntersection() can be used.  It has the
        // effect of increasing the snap radius for edges (but not vertices) by
        // the given distance.
        //
        // The intersection tolerance should be set to the maximum error in the
        // intersection calculation used.  For example, if S2::GetIntersection()
        // is used then the error should be set to S2::kIntersectionError.  If
        // S2::GetPointOnLine() is used then the error should be set to
        // S2::kGetPointOnLineError.  If S2::Project() is used then the error
        // should be set to S2::kProjectPerpendicularError.  If more than one
        // method is used then the intersection tolerance should be set to the
        // maximum such error.
        //
        // The reason this option is necessary is that computed intersection
        // points are not exact.  For example, S2::GetIntersection(a, b, c, d)
        // returns a point up to S2::kIntersectionError away from the true
        // mathematical intersection of the edges AB and CD.  Furthermore such
        // intersection points are subject to further snapping in order to ensure
        // that no pair of vertices is closer than the specified snap radius.  For
        // example, suppose the computed intersection point X of edges AB and CD
        // is 1 nanonmeter away from both edges, and the snap radius is 1 meter.
        // In that case X might snap to another vertex Y exactly 1 meter away,
        // which would leave us with a vertex Y that could be up to 1.000000001
        // meters from the edges AB and/or CD.  This means that AB and/or CD might
        // not snap to Y leaving us with two edges that still cross each other.
        //
        // However if the intersection tolerance is set to 1 nanometer then the
        // snap radius for edges is increased to 1.000000001 meters ensuring that
        // both edges snap to a common vertex even in this worst case.  (Tthis
        // technique does not work if the vertex snap radius is increased as well;
        // it requires edges and vertices to be handled differently.)
        //
        // Note that this option allows edges to move by up to the given
        // intersection tolerance even when the snap radius is zero.  The exact
        // distance that edges can move is always given by max_edge_deviation()
        // defined above.
        //
        // When split_crossing_edges() is true, the intersection tolerance is
        // automatically set to a minimum of S2::kIntersectionError.  A larger
        // value can be specified by calling this method explicitly.
        //
        // DEFAULT: S1Angle::Zero()
        public S1Angle IntersectionTolerance
        { 
            get
            {
                if (!SplitCrossingEdges) return intersection_tolerance_;
                return S1Angle.Max(intersection_tolerance_, S2.kIntersectionErrorS1Angle);
            }
            set
            {
                Debug.Assert(value >= S1Angle.Zero);
                intersection_tolerance_ = value;
            }
        }
        private S1Angle intersection_tolerance_ = S1Angle.Zero;

        // Specifies that internal memory usage should be tracked using the given
        // S2MemoryTracker.  If a memory limit is specified and more more memory
        // than this is required then an error will be returned.  Example usage:
        //
        //   S2MemoryTracker tracker;
        //   tracker.set_limit(500 << 20);  // 500 MB
        //   S2Builder::Options options;
        //   options.set_memory_tracker(&tracker);
        //   S2Builder builder{options};
        //   ...
        //   S2Error error;
        //   if (!builder.Build(&error)) {
        //     if (error.code() == S2Error::RESOURCE_EXHAUSTED) {
        //       S2_LOG(ERROR) << error;  // Memory limit exceeded
        //     }
        //   }
        //
        // CAVEATS:
        //
        //  - Memory allocated by the output S2Builder layers is not tracked.
        //
        //  - While memory tracking is reasonably complete and accurate, it does
        //    not account for every last byte.  It is intended only for the
        //    purpose of preventing clients from running out of memory.
        //
        // DEFAULT: nullptr (memory tracking disabled)
        public S2MemoryTracker? MemoryTracker { get; set; } = null;

        // If true, then simplify the output geometry by replacing nearly straight
        // chains of short edges with a single long edge.
        //
        // The combined effect of snapping and simplifying will not change the
        // input by more than the guaranteed tolerances (see the list documented
        // with the SnapFunction class).  For example, simplified edges are
        // guaranteed to pass within snap_radius() of the *original* positions of
        // all vertices that were removed from that edge.  This is a much tighter
        // guarantee than can be achieved by snapping and simplifying separately.
        //
        // However, note that this option does not guarantee idempotency.  In
        // other words, simplifying geometry that has already been simplified once
        // may simplify it further.  (This is unavoidable, since tolerances are
        // measured with respect to the original geometry, which is no longer
        // available when the geometry is simplified a second time.)
        //
        // When the output consists of multiple layers, simplification is
        // guaranteed to be consistent: for example, edge chains are simplified in
        // the same way across layers, and simplification preserves topological
        // relationships between layers (e.g., no crossing edges will be created).
        // Note that edge chains in different layers do not need to be identical
        // (or even have the same number of vertices, etc) in order to be
        // simplified together.  All that is required is that they are close
        // enough together so that the same simplified edge can meet all of their
        // individual snapping guarantees.
        //
        // Note that edge chains are approximated as parametric curves rather than
        // point sets.  This means that if an edge chain backtracks on itself (for
        // example, ABCDEFEDCDEFGH) then such backtracking will be preserved to
        // within snap_radius() (for example, if the preceding point were all in a
        // straight line then the edge chain would be simplified to ACFCFH, noting
        // that C and F have degree > 2 and therefore can't be simplified away).
        //
        // Simplified edges are assigned all labels associated with the edges of
        // the simplified chain.
        //
        // For this option to have any effect, a SnapFunction with a non-zero
        // snap_radius() must be specified.  Also note that vertices specified
        // using ForceVertex are never simplified away.
        //
        // DEFAULT: false
        public bool SimplifyEdgeChains
        {
            get => _simplifyEdgeChains;
            set
            {
                _simplifyEdgeChains = value;

                // Simplification requires a non-zero snap radius, and while it might be
                // possible to do some simplifying without snapping, it is much simpler to
                // always snap (even if the input geometry already meets the other output
                // requirements).  We need to compute edge_sites_ in order to avoid
                // approaching non-incident vertices too closely, for example.
                Idempotent = false;
            }
        }
        private bool _simplifyEdgeChains = false;

        // If true, then snapping occurs only when the input geometry does not
        // already meet the S2Builder output guarantees (see the SnapFunction
        // class description for details).  This means that if all input vertices
        // are at snapped locations, all vertex pairs are separated by at least
        // min_vertex_separation(), and all edge-vertex pairs are separated by at
        // least min_edge_vertex_separation(), then no snapping is done.
        //
        // If false, then all vertex pairs and edge-vertex pairs closer than
        // "snap_radius" will be considered for snapping.  This can be useful, for
        // example, if you know that your geometry contains errors and you want to
        // make sure that features closer together than "snap_radius" are merged.
        //
        // This option is automatically turned off by simplify_edge_chains(),
        // since simplifying edge chains is never guaranteed to be idempotent.
        //
        // DEFAULT: true
        public bool Idempotent { get; set; } = true;

        // The maximum distance that any point along an edge can move when snapped.
        // It is slightly larger than edge_snap_radius() because when a geodesic
        // edge is snapped, the edge center moves further than its endpoints.
        // S2Builder ensures that this distance is at most 10% larger than
        // edge_snap_radius().
        public S1Angle MaxEdgeDeviation()
        {
            // We want max_edge_deviation() to be large enough compared to snap_radius()
            // such that edge splitting is rare.
            //
            // Using spherical trigonometry, if the endpoints of an edge of length L
            // move by at most a distance R, the center of the edge moves by at most
            // asin(sin(R) / cos(L / 2)).  Thus the (max_edge_deviation / snap_radius)
            // ratio increases with both the snap radius R and the edge length L.
            //
            // We arbitrarily limit the edge deviation to be at most 10% more than the
            // snap radius.  With the maximum allowed snap radius of 70 degrees, this
            // means that edges up to 30.6 degrees long are never split.  For smaller
            // snap radii, edges up to 49 degrees long are never split.  (Edges of any
            // length are not split unless their endpoints move far enough so that the
            // actual edge deviation exceeds the limit; in practice, splitting is rare
            // even with long edges.)  Note that it is always possible to split edges
            // when max_edge_deviation() is exceeded; see MaybeAddExtraSites().
            Debug.Assert(SnapFunction.SnapRadius <= SnapFunction.kMaxSnapRadius);
            double kMaxEdgeDeviationRatio = 1.1;
            return kMaxEdgeDeviationRatio * EdgeSnapRadius();
        }
    }

    // An S2Shape used to represent the entire collection of S2Builder input edges.
    // Vertices are specified as indices into a vertex vector to save space.
    public sealed class VertexIdEdgeVectorShape : S2Shape
    {
        // Requires that "edges" isant for the lifetime of this object.
        public VertexIdEdgeVectorShape(List<InputEdge> edges, List<S2Point> vertices)
        {
            edges_ = edges; vertices_ = vertices;
        }

        public S2Point Vertex0(int e) => Vertex(edges_[e].Item1);
        public S2Point Vertex1(int e) => Vertex(edges_[e].Item2);

        // S2Shape interface:
        public override int NumEdges() => edges_.Count;
        public override Edge GetEdge(int e) => new(vertices_[edges_[e].Item1], vertices_[edges_[e].Item2]);
        public override int Dimension() => 1;
        public override ReferencePoint GetReferencePoint() => ReferencePoint.FromContained(false);
        public override int NumChains() => edges_.Count;
        public override Chain GetChain(int i) => new(i, 1);
        public override Edge ChainEdge(int i, int j) => GetEdge(i);
        public override ChainPosition GetChainPosition(int e) => new(e, 0);

        private S2Point Vertex(int i) => vertices_[i];

        private readonly List<InputEdge> edges_;
        private readonly List<S2Point> vertices_;
    }

    //////////////////////  Internal Types  /////////////////////////

    // A class that encapsulates the state needed for simplifying edge chains.
    public class EdgeChainSimplifier
    {
        // The graph "g" contains all edges from all layers.  "edge_layers"
        // indicates the original layer for each edge.  "site_vertices" is a map
        // from SiteId to the set of InputVertexIds that were snapped to that site.
        // "layer_edges" and "layer_input_edge_ids" are output arguments where the
        // simplified edge chains will be placed.  The input and output edges are
        // not sorted.
        public EdgeChainSimplifier(
          S2Builder builder, Graph g,
          List<int> edge_layers,
          List<List<int>> site_vertices,
          List<List<OutputEdge>> layer_edges,
          List<List<int>> layer_input_edge_ids,
          IdSetLexicon input_edge_id_set_lexicon)
        {
            builder_ = builder;
            g_ = g;
            in_ = new Graph.VertexInMap(g);
            out_ = new Graph.VertexOutMap(g);
            edge_layers_ = edge_layers;
            site_vertices_ = site_vertices;
            layer_edges_ = layer_edges;
            layer_input_edge_ids_ = layer_input_edge_ids;
            input_edge_id_set_lexicon_ = input_edge_id_set_lexicon;
            layer_begins_ = builder_.layer_begins_;
            is_interior_ = new List<bool>(g.NumVertices);
            used_ = new List<bool>(g.NumEdges);
            tmp_vertex_set_ = new(16); /*expected_max_elements*/
            new_edges_.Capacity = g.NumEdges;
            new_input_edge_ids_.Capacity = g.NumEdges;
            new_edge_layers_.Capacity = g.NumEdges;
        }

        public void Run()
        {
            // Determine which vertices can be interior vertices of an edge chain.
            for (var v = 0; v < g_.NumVertices; ++v)
            {
                is_interior_[v] = IsInterior(v);
            }
            // Attempt to simplify all edge chains that start from a non-interior
            // vertex.  (This takes care of all chains except loops.)
            for (var e = 0; e < g_.NumEdges; ++e)
            {
                if (used_[e]) continue;
                var edge = g_.GetEdge(e);
                if (is_interior_[edge.ShapeId]) continue;
                if (!is_interior_[edge.EdgeId])
                {
                    OutputEdge(e);  // An edge between two non-interior vertices.
                }
                else
                {
                    SimplifyChain(edge.ShapeId, edge.EdgeId);
                }
            }
            // If there are any edges left, they form one or more disjoint loops where
            // all vertices are interior vertices.
            //
            // TODO(ericv): It would be better to start from the edge with the smallest
            // min_input_edge_id(), since that would make the output more predictable
            // for testing purposes.  It also means that we won't create an edge that
            // spans the start and end of a polyline if the polyline is snapped into a
            // loop.  (Unfortunately there are pathological examples that prevent us
            // from guaranteeing this in general, e.g. there could be two polylines in
            // different layers that snap to the same loop but start at different
            // positions.  In general we only consider input edge ids to be a hint
            // towards the preferred output ordering.)
            for (var e = 0; e < g_.NumEdges; ++e)
            {
                if (used_[e]) continue;
                var edge = g_.GetEdge(e);
                if (edge.ShapeId == edge.EdgeId)
                {
                    // Note that it is safe to output degenerate edges as we go along,
                    // because this vertex has at least one non-degenerate outgoing edge and
                    // therefore we will (or just did) start an edge chain here.
                    OutputEdge(e);
                }
                else
                {
                    SimplifyChain(edge.ShapeId, edge.EdgeId);
                }
            }
            // TODO(ericv): The graph is not needed past here, so we could save some
            // memory by clearing the underlying Edge and InputEdgeIdSetId vectors.

            // Finally, copy the output edges into the appropriate layers.  They don't
            // need to be sorted because the input edges were also unsorted.
            for (var e = 0; e < new_edges_.Count; ++e)
            {
                var layer = new_edge_layers_[e];
                layer_edges_[layer].Add(new_edges_[e]);
                layer_input_edge_ids_[layer].Add(new_input_edge_ids_[e]);
            }
        }

        // A helper class for determining whether a vertex can be an interior vertex
        // of a simplified edge chain.  Such a vertex must be adjacent to exactly two
        // vertices (across all layers combined), and in each layer the number of
        // incoming edges from one vertex must equal the number of outgoing edges to
        // the other vertex (in both directions).  Furthermore the vertex cannot have
        // any degenerate edges in a given layer unless it has at least one
        // non-degenerate edge in that layer as well.  (Note that usually there will
        // not be any degenerate edges at all, since most layer types discard them.)
        //
        // The last condition is necessary to prevent the following: suppose that one
        // layer has a chain ABC and another layer has a degenerate edge BB (with no
        // other edges incident to B).  Then we can't simplify ABC to AC because there
        // would be no suitable replacement for the edge BB (since the input edge that
        // mapped to BB can't be replaced by any of the edges AA, AC, or CC without
        // moving further than snap_radius).
        public class InteriorVertexMatcher
        {
            // Checks whether "v0" can be an interior vertex of an edge chain.
            public InteriorVertexMatcher(int v0)
            {
                v0_ = v0; v1_ = -1; v2_ = -1; n0_ = 0; n1_ = 0; n2_ = 0; excess_out_ = 0;
                too_many_endpoints_ = false;
            }

            // Starts analyzing the edges of a new layer.
            public void StartLayer()
            {
                excess_out_ = n0_ = n1_ = n2_ = 0;
            }

            // This method should be called for each edge incident to "v0" in a given
            // layer.  (For degenerate edges, it should be called twice.)
            public void Tally(int v, bool outgoing)
            {
                excess_out_ += outgoing ? 1 : -1;  // outdegree - indegree
                if (v == v0_)
                {
                    ++n0_;  // Counts both endpoints of each degenerate edge.
                }
                else
                {
                    // We keep track of the total number of edges (incoming or outgoing)
                    // connecting v0 to up to two adjacent vertices.
                    if (v1_ < 0) v1_ = v;
                    if (v1_ == v)
                    {
                        ++n1_;
                    }
                    else
                    {
                        if (v2_ < 0) v2_ = v;
                        if (v2_ == v)
                        {
                            ++n2_;
                        }
                        else
                        {
                            too_many_endpoints_ = true;
                        }
                    }
                }
            }

            // This method should be called after processing the edges for each layer.
            // It returns true if "v0" is an interior vertex based on the edges so far.
            public bool Matches()
            {
                // We check that there are the same number of incoming and outgoing edges
                // in each direction by verifying that (1) indegree(v0) == outdegree(v0)
                // and (2) the total number of edges (incoming and outgoing) to "v1" and
                // "v2" are equal.  We also check the condition on degenerate edges that
                // is documented above.
                return (!too_many_endpoints_ && excess_out_ == 0 && n1_ == n2_ &&
                        (n0_ == 0 || n1_ > 0));
            }

            private readonly int v0_;
            private int v1_;
            private int v2_;
            private int n0_, n1_, n2_;
            private int excess_out_;           // outdegree(v0) - indegree(v0)
            private bool too_many_endpoints_;  // Have we seen more than two adjacent vertices?
        }

        // Copies the given edge to the output and marks it as used.
        private void OutputEdge(int e)
        {
            new_edges_.Add(g_.GetEdge(e));
            new_input_edge_ids_.Add(g_.InputEdgeIdSetId(e));
            new_edge_layers_.Add(edge_layers_[e]);
            used_[e] = true;
        }
        // Returns the layer than a given input edge belongs to.
        private int InputEdgeLayer(int id)
        {
            // NOTE(ericv): If this method shows up in profiling, the result could be
            // stored with each edge (i.e., edge_layers_ and new_edge_layers_).
            Debug.Assert(id >= 0);
            return layer_begins_.GetUpperBound(id) - 1;
        }
        // Returns true if VertexId "v" can be an interior vertex of a simplified edge
        // chain.  (See the InteriorVertexMatcher class for what this implies.)
        private bool IsInterior(int v)
        {
            // Check a few simple prerequisites.
            if (out_.Degree(v) == 0) return false;
            if (out_.Degree(v) != in_.Degree(v)) return false;
            if (builder_.IsForced(v)) return false;  // Keep forced vertices.

            // Sort the edges so that they are grouped by layer.
            var edges = new List<int>();
            foreach (var e in out_.EdgeIds(v)) edges.Add(e);
            foreach (var e in in_.EdgeIds(v)) edges.Add(e);
            edges.Sort(new GraphEdgeComp(edge_layers_));
            // Now feed the edges in each layer to the InteriorVertexMatcher.
            var matcher = new InteriorVertexMatcher(v);
            for (var e = 0; e < edges.Count;)
            {
                var layer = edge_layers_[e];
                matcher.StartLayer();
                for (; e < edges.Count && edge_layers_[e] == layer; ++e)
                {
                    var edge = g_.GetEdge(edges[e]);
                    if (edge.ShapeId == v) matcher.Tally(edge.EdgeId, true /*outgoing*/);
                    if (edge.EdgeId == v) matcher.Tally(edge.ShapeId, false /*outgoing*/);
                }
                if (!matcher.Matches()) return false;
            }
            return true;
        }
        private class GraphEdgeComp : IComparer<int>
        {
            private readonly List<int> edge_layers_;
            public GraphEdgeComp(List<int> edge_layers) { edge_layers_ = edge_layers; }
            public int Compare(int x, int y)
            {
                return edge_layers_[x].CompareTo(edge_layers_[y]);
            }
        }
        // Follows the edge chain starting with (v0, v1) until either we find a
        // non-interior vertex or we return to the original vertex v0.  At each vertex
        // we simplify a subchain of edges that is as long as possible.
        private void SimplifyChain(int v0, int v1)
        {
            // Avoid allocating "chain" each time by reusing it.
            var chain = tmp_vertices_;
            // Contains the set of vertices that have either been avoided or added to
            // the chain so far.  This is necessary so that AvoidSites() doesn't try to
            // avoid vertices that have already been added to the chain.
            Dictionary<VertexId, VertexId> used_vertices = tmp_vertex_set_;
            S2PolylineSimplifier simplifier;
            var vstart = v0;
            bool done;
            do
            {
                // Simplify a subchain of edges starting with (v0, v1).
                chain.Add(v0);
                used_vertices.Add(v0, v0);
                simplifier = new(g_.Vertex(v0));
                // Note that if the first edge (v0, v1) is longer than the maximum length
                // allowed for simplification, then AvoidSites() will return false and we
                // exit the loop below after the first iteration.
                var simplify = AvoidSites(v0, v0, v1, used_vertices, simplifier);
                do
                {
                    chain.Add(v1);
                    used_vertices.Add(v1, v1);
                    done = !is_interior_[v1] || v1 == vstart;
                    if (done) break;

                    // Attempt to extend the chain to the next vertex.
                    var vprev = v0;
                    v0 = v1;
                    v1 = FollowChain(vprev, v0);
                } while (simplify && TargetInputVertices(v0, simplifier) &&
                         AvoidSites(chain[0], v0, v1, used_vertices, simplifier) &&
                         simplifier.Extend(g_.Vertex(v1)));

                if (chain.Count == 2)
                {
                    OutputAllEdges(chain[0], chain[1]);  // Could not simplify.
                }
                else
                {
                    MergeChain(chain);
                }
                // Note that any degenerate edges that were not merged into a chain are
                // output by EdgeChainSimplifier.Run().
                chain.Clear();
                used_vertices.Clear();
            } while (!done);
        }

        // Given an edge (v0, v1) where v1 is an interior vertex, returns the (unique)
        // next vertex in the edge chain.
        private int FollowChain(int v0, int v1)
        {
            Debug.Assert(is_interior_[v1]);
            foreach (var e in out_.EdgeIds(v1))
            {
                var v = g_.GetEdge(e).EdgeId;
                if (v != v0 && v != v1) return v;
            }
            throw new ApplicationException("Could not find next edge in edge chain");
        }

        // Copies all input edges between v0 and v1 (in both directions) to the output.
        private void OutputAllEdges(int v0, int v1)
        {
            foreach (var e in out_.EdgeIds(v0, v1)) OutputEdge(e);
            foreach (var e in out_.EdgeIds(v1, v0)) OutputEdge(e);
        }

        // Ensures that the simplified edge passes within "edge_snap_radius" of all
        // the *input* vertices that snapped to the given vertex "v".
        private bool TargetInputVertices(int v, S2PolylineSimplifier simplifier)
        {
            foreach (var i in site_vertices_[v])
            {
                if (!simplifier.TargetDisc(builder_.input_vertices_[i],
                                            builder_.edge_snap_radius_ca_))
                {
                    return false;
                }
            }
            return true;
        }

        // Given the starting vertex v0 and last edge (v1, v2) of an edge chain,
        // restricts the allowable range of angles in order to ensure that all sites
        // near the edge (v1, v2) are avoided by at least min_edge_vertex_separation.
        private bool AvoidSites(int v0, int v1, int v2,
            Dictionary<VertexId, VertexId> used_vertices,
            S2PolylineSimplifier simplifier)
        {
            var p0 = g_.Vertex(v0);
            var p1 = g_.Vertex(v1);
            var p2 = g_.Vertex(v2);
            var r1 = new S1ChordAngle(p0, p1);
            var r2 = new S1ChordAngle(p0, p2);

            // The distance from the start of the edge chain must increase monotonically
            // for each vertex, since we don't want to simplify chains that backtrack on
            // themselves (we want a parametric approximation, not a geometric one).
            if (r2 < r1) return false;

            // We also limit the maximum edge length in order to guarantee that the
            // simplified edge stays with max_edge_deviation() of all the input edges
            // that snap to it.
            if (r2 >= builder_.min_edge_length_to_split_ca_) return false;

            // Otherwise it is sufficient to consider the nearby sites (edge_sites_) for
            // a single input edge that snapped to (v1, v2) or (v2, v1).  This is
            // because each edge has a list of all sites within (max_edge_deviation +
            // min_edge_vertex_separation), and since the output edge is within
            // max_edge_deviation of all input edges, this list includes all sites
            // within min_edge_vertex_separation of the output edge.
            //
            // Usually there is only one edge to choose from, but it's not much more
            // effort to choose the edge with the shortest list of edge_sites_.
            var best = -1;
            var edge_sites = builder_.edge_sites_;
            foreach (var e in out_.EdgeIds(v1, v2))
            {
                foreach (var id in g_.InputEdgeIds(e))
                {
                    if (best < 0 || edge_sites[id].Count < edge_sites[best].Count)
                        best = id;
                }
            }
            foreach (var e in out_.EdgeIds(v2, v1))
            {
                foreach (var id in g_.InputEdgeIds(e))
                {
                    if (best < 0 || edge_sites[id].Count < edge_sites[best].Count)
                        best = id;
                }
            }
            Debug.Assert(best >= 0);  // Because there is at least one outgoing edge.

            foreach (var v in edge_sites[best])
            {
                // Sites whose distance from "p0" is at least "r2" are not relevant yet.
                var p = g_.Vertex(v);
                var r = new S1ChordAngle(p0, p);
                if (r >= r2) continue;

                // The following test prevents us from avoiding previous vertices of the
                // edge chain that also happen to be nearby the current edge.  (It also
                // happens to ensure that each vertex is avoided at most once, but this is
                // just an optimization.)
                if (used_vertices.ContainsKey(v)) continue;

                used_vertices.Add(v, v);

                // We need to figure out whether this site is to the left or right of the
                // edge chain.  For the first edge this is easy.  Otherwise, since we are
                // only considering sites in the radius range (r1, r2), we can do this by
                // checking whether the site is to the left of the wedge (p0, p1, p2).
                var disc_on_left = (v1 == v0) ? (S2Pred.Sign(p1, p2, p) > 0)
                                               : S2Pred.OrderedCCW(p0, p2, p, p1);
                if (!simplifier.AvoidDisc(p, builder_.min_edge_site_separation_ca_,
                                           disc_on_left))
                {
                    return false;
                }
            }
            return true;
        }

        // Given the vertices in a simplified edge chain, adds the corresponding
        // simplified edge(s) to the output.  Note that (1) the edge chain may exist
        // in multiple layers, (2) the edge chain may exist in both directions, (3)
        // there may be more than one copy of an edge chain (in either direction)
        // within a single layer.
        private void MergeChain(List<int> vertices)
        {
            // Suppose that all interior vertices have M outgoing edges and N incoming
            // edges.  Our goal is to group the edges into M outgoing chains and N
            // incoming chains, and then replace each chain by a single edge.
            var merged_input_ids = new List<List<int>>();
            var degenerate_ids = new List<int>();
            int num_out;  // Edge count in the outgoing direction.
            for (var i = 1; i < vertices.Count; ++i)
            {
                var v0 = vertices[i - 1];
                var v1 = vertices[i];
                var out_edges = out_.EdgeIds(v0, v1);
                var in_edges = out_.EdgeIds(v1, v0);
                if (i == 1)
                {
                    // Allocate space to store the input edge ids associated with each edge.
                    num_out = out_edges.Count;
                    merged_input_ids.Capacity = num_out + in_edges.Count;
                    foreach (var ids in merged_input_ids)
                    {
                        ids.Capacity = vertices.Count - 1;
                    }
                }
                else
                {
                    // For each interior vertex, we build a list of input edge ids
                    // associated with degenerate edges.  Each input edge ids will be
                    // assigned to one of the output edges later.  (Normally there are no
                    // degenerate edges at all since most layer types don't want them.)
                    Debug.Assert(is_interior_[v0]);
                    foreach (var e in out_.EdgeIds(v0, v0))
                    {
                        foreach (var id in g_.InputEdgeIds(e))
                        {
                            degenerate_ids.Add(id);
                        }
                        used_[e] = true;
                    }
                }
                // Because the edges were created in layer order, and all sorts used are
                // stable, the edges are still in layer order.  Therefore we can simply
                // merge together all the edges in the same relative position.
                var j = 0;
                foreach (var e in out_edges)
                {
                    foreach (var id in g_.InputEdgeIds(e))
                    {
                        merged_input_ids[j].Add(id);
                    }
                    used_[e] = true;
                    ++j;
                }
                foreach (var e in in_edges)
                {
                    foreach (var id in g_.InputEdgeIds(e))
                    {
                        merged_input_ids[j].Add(id);
                    }
                    used_[e] = true;
                    ++j;
                }
                Debug.Assert(merged_input_ids.Count == j);
            }
            if (degenerate_ids.Any())
            {
                degenerate_ids.Sort();
                AssignDegenerateEdges(degenerate_ids, merged_input_ids);
            }
            // Output the merged edges.
            int v0_ = vertices[0], v1_ = vertices[1], vb = vertices.Last();
            foreach (var e in out_.EdgeIds(v0_, v1_))
            {
                new_edges_.Add(new(v0_, vb));
                new_edge_layers_.Add(edge_layers_[e]);
            }
            foreach (var e in out_.EdgeIds(v1_, v0_))
            {
                new_edges_.Add(new(vb, v0_));
                new_edge_layers_.Add(edge_layers_[e]);
            }
            foreach (var ids in merged_input_ids)
            {
                new_input_edge_ids_.Add(input_edge_id_set_lexicon_.Add(ids));
            }
        }

        // Given a list of the input edge ids associated with degenerate edges in the
        // interior of an edge chain, assigns each input edge id to one of the output
        // edges.
        private void AssignDegenerateEdges(List<int> degenerate_ids, List<List<int>> merged_ids)
        {
            // Each degenerate edge is assigned to an output edge in the appropriate
            // layer.  If there is more than one candidate, we use heuristics so that if
            // the input consists of a chain of edges provided in consecutive order
            // (some of which became degenerate), then all those input edges are
            // assigned to the same output edge.  For example, suppose that one output
            // edge is labeled with input edges 3,4,7,8, while another output edge is
            // labeled with input edges 50,51,54,57.  Then if we encounter degenerate
            // edges labeled with input edges 5 and 6, we would prefer to assign them to
            // the first edge (yielding the continuous range 3,4,5,6,7,8).
            //
            // The heuristic below is only smart enough to handle the case where the
            // candidate output edges have non-overlapping ranges of input edges.
            // (Otherwise there is probably not a good heuristic for assigning the
            // degenerate edges in any case.)

            // Duplicate edge ids don't affect the heuristic below, so we don't bother
            // removing them.  (They will be removed by IdSetLexicon.Add.)
            foreach (var ids in merged_ids) ids.Sort();

            // Sort the output edges by their minimum input edge id.  This is sufficient
            // for the purpose of determining which layer they belong to.  With
            // EdgeType.UNDIRECTED, some edges might not have any input edge ids (i.e.,
            // if they consist entirely of siblings of input edges).  We simply remove
            // such edges from the lists of candidates.
            var order = new List<int>(merged_ids.Count);
            for (var i = 0; i < merged_ids.Count; ++i)
            {
                if (merged_ids[i].Any()) order.Add(i);
            }
            order.Sort((int i, int j) => merged_ids[i][0].CompareTo(merged_ids[j][0]));

            // Now determine where each degenerate edge should be assigned.
            var comp = new MergedIdsComp(merged_ids);
            foreach (var degenerate_id in degenerate_ids)
            {
                var layer = InputEdgeLayer(degenerate_id);

                // Find the first output edge whose range of input edge ids starts after
                // "degenerate_id".  If the previous edge belongs to the correct layer,
                // then we assign the degenerate edge to it.
                var index = order.GetUpperBound(degenerate_id, comp);
                if (index < order.Count)
                {
                    if (merged_ids[order[index - 1]][0] >= layer_begins_[layer]) --index;
                }
                Debug.Assert(layer == InputEdgeLayer(merged_ids[(int)order[index]][0]));
                merged_ids[order[index]].Add(degenerate_id);
            }
        }

        private class MergedIdsComp : IComparer<int>
        {
            private readonly List<List<int>> merged_input_ids_;
            public MergedIdsComp(List<List<int>> merged_input_ids) { merged_input_ids_ = merged_input_ids; }

            public int Compare(int x, int y)
            {
                return x.CompareTo(merged_input_ids_[y][0]);
            }
        }

        private readonly S2Builder builder_;
        private readonly Graph g_;
        private readonly Graph.VertexInMap in_;
        private readonly Graph.VertexOutMap out_;
        private readonly List<int> edge_layers_;
        private readonly List<List<int>> site_vertices_;
        private readonly List<List<OutputEdge>> layer_edges_;
        private readonly List<List<int>> layer_input_edge_ids_;
        private readonly IdSetLexicon input_edge_id_set_lexicon_;

        // Convenience member copied from builder_.
        private readonly List<int> layer_begins_;

        // is_interior_[v] indicates that VertexId "v" is eligible to be an interior
        // vertex of a simplified edge chain.  You can think of it as vertex whose
        // indegree and outdegree are both 1 (although the actual definition is a
        // bit more complicated because of duplicate edges and layers).
        private readonly List<bool> is_interior_;

        // used_[e] indicates that EdgeId "e" has already been processed.
        private readonly List<bool> used_;

        // Temporary objects declared here to avoid repeated allocation.
        private readonly List<int> tmp_vertices_ = new();

        private readonly Dictionary<VertexId, VertexId> tmp_vertex_set_ = new();

        // The output edges after simplification.
        private readonly List<OutputEdge> new_edges_ = new();
        private readonly List<int> new_input_edge_ids_ = new();
        private readonly List<int> new_edge_layers_ = new();
    }

    // This class is only needed by S2Builder.Layer implementations.  A layer is
    // responsible for assembling an S2Builder.Graph of snapped edges into the
    // desired output format (e.g., an S2Polygon).  The GraphOptions class allows
    // each Layer type to specify requirements on its input graph: for example, if
    // DegenerateEdges.DISCARD is specified, then S2Builder will ensure that all
    // degenerate edges are removed before passing the graph to the S2Layer.Build
    // method.
    public class GraphOptions : IEquatable<GraphOptions>
    {
        #region Fields, Constants

        // Specifies whether the S2Builder input edges should be treated as
        // undirected.  If true, then all input edges are duplicated into pairs
        // consisting of an edge and a sibling (reverse) edge.  Note that the
        // automatically created sibling edge has an empty set of labels and does
        // not have an associated InputEdgeId.
        //
        // The layer implementation is responsible for ensuring that exactly one
        // edge from each pair is used in the output, i.e. *only half* of the graph
        // edges will be used.  (Note that some values of the sibling_pairs() option
        // automatically take care of this issue by removing half of the edges and
        // changing edge_type() to DIRECTED.)
        //
        // DEFAULT: EdgeType.DIRECTED
        public EdgeType EdgeType_ { get; set; }

        public DegenerateEdges DegenerateEdges_ { get; set; }

        public DuplicateEdges DuplicateEdges_ { get; set; }

        public SiblingPairs SiblingPairs_ { get; set; }

        // This is a specialized option that is only needed by clients want to work
        // with the graphs for multiple layers at the same time (e.g., in order to
        // check whether the same edge is present in two different graphs).  [Note
        // that if you need to do this, usually it is easier just to build a single
        // graph with suitable edge labels.]
        //
        // When there are a large number of layers, then by default S2Builder builds
        // a minimal subgraph for each layer containing only the vertices needed by
        // the edges in that layer.  This ensures that layer types that iterate over
        // the vertices run in time proportional to the size of that layer rather
        // than the size of all layers combined.  (For example, if there are a
        // million layers with one edge each, then each layer would be passed a
        // graph with 2 vertices rather than 2 million vertices.)
        //
        // If this option is set to false, this optimization is disabled.  Instead
        // the graph passed to this layer will contain the full set of vertices.
        // (This is not recommended when the number of layers could be large.)
        //
        // DEFAULT: true
        public bool AllowVertexFiltering { get; set; }

        #endregion

        #region Constructors

        // All S2Builder.Layer subtypes should specify GraphOptions explicitly
        // using this constructor, rather than relying on default values.
        public GraphOptions(EdgeType edge_type, DegenerateEdges degenerate_edges,
                       DuplicateEdges duplicate_edges, SiblingPairs sibling_pairs)
        {
            EdgeType_ = edge_type; DegenerateEdges_ = degenerate_edges;
            DuplicateEdges_ = duplicate_edges; SiblingPairs_ = sibling_pairs;
            AllowVertexFiltering = true;
        }

        // The default options specify that all edges should be kept, since this
        // produces the least surprising output and makes it easier to diagnose the
        // problem when an option is left unspecified.
        public GraphOptions()
        {
            EdgeType_ = EdgeType.DIRECTED;
            DegenerateEdges_ = DegenerateEdges.KEEP;
            DuplicateEdges_ = DuplicateEdges.KEEP;
            SiblingPairs_ = SiblingPairs.KEEP;
            AllowVertexFiltering = true;
        }

        #endregion

        #region IEquatable

        public static bool operator ==(GraphOptions x, GraphOptions y) => Equals(x, y);

        public static bool operator !=(GraphOptions x, GraphOptions y) => !Equals(x, y);
        public override bool Equals(object? other) => other is GraphOptions go && Equals(go);
        public bool Equals(GraphOptions? other)
        {
            if (other is null) return false;

            return EdgeType_ == other.EdgeType_ &&
                DegenerateEdges_ == other.DegenerateEdges_ &&
                DuplicateEdges_ == other.DuplicateEdges_ &&
                SiblingPairs_ == other.SiblingPairs_ &&
                AllowVertexFiltering == other.AllowVertexFiltering;
        }
        public override int GetHashCode()
        {
            return HashCode.Combine(EdgeType_, DegenerateEdges_, DuplicateEdges_, SiblingPairs_, AllowVertexFiltering);
        }

        #endregion

        // Controls how degenerate edges (i.e., an edge from a vertex to itself) are
        // handled.  Such edges may be present in the input, or they may be created
        // when both endpoints of an edge are snapped to the same output vertex.
        // The options available are:
        //
        // DISCARD: Discards all degenerate edges.  This is useful for layers that
        //          do not support degeneracies, such as S2PolygonLayer.
        //
        // DISCARD_EXCESS: Discards all degenerate edges that are connected to
        //                 non-degenerate edges and merges any remaining duplicate
        //                 degenerate edges.  This is useful for simplifying
        //                 polygons while ensuring that loops that collapse to a
        //                 single point do not disappear.
        //
        // KEEP: Keeps all degenerate edges.  Be aware that this may create many
        //       redundant edges when simplifying geometry (e.g., a polyline of the
        //       form AABBBBBCCCCCCDDDD).  DegenerateEdges.KEEP is mainly useful
        //       for algorithms that require an output edge for every input edge.
        //
        // DEFAULT: DegenerateEdges.KEEP
        public enum DegenerateEdges : byte { DISCARD, DISCARD_EXCESS, KEEP }

        // Controls how duplicate edges (i.e., edges that are present multiple
        // times) are handled.  Such edges may be present in the input, or they can
        // be created when vertices are snapped together.  When several edges are
        // merged, the result is a single edge labelled with all of the original
        // input edge ids.
        //
        // DEFAULT: DuplicateEdges.KEEP
        public enum DuplicateEdges : byte { MERGE, KEEP }

        // Controls how sibling edge pairs (i.e., pairs consisting of an edge and
        // its reverse edge) are handled.  Layer types that define an interior
        // (e.g., polygons) normally discard such edge pairs since they do not
        // affect the result (i.e., they define a "loop" with no interior).  The
        // various options include:
        //
        // DISCARD: Discards all sibling edge pairs.
        //
        // DISCARD_EXCESS: Like DISCARD, except that a single sibling pair is kept
        //                 if the result would otherwise be empty.  This is useful
        //                 for polygons with degeneracies (S2LaxPolygonShape), and
        //                 for simplifying polylines while ensuring that they are
        //                 not split into multiple disconnected pieces.
        //
        // KEEP: Keeps sibling pairs.  This can be used to create polylines that
        //       double back on themselves, or degenerate loops (with a layer type
        //       such as S2LaxPolygonShape).
        //
        // REQUIRE: Requires that all edges have a sibling (and returns an error
        //          otherwise).  This is useful with layer types that create a
        //          collection of adjacent polygons (a polygon mesh).
        //
        // CREATE: Ensures that all edges have a sibling edge by creating them if
        //         necessary.  This is useful with polygon meshes where the input
        //         polygons do not cover the entire sphere.  Such edges always have
        //         an empty set of labels and do not have an associated InputEdgeId.
        //
        // If edge_type() is EdgeType.UNDIRECTED, a sibling edge pair is considered
        // to consist of four edges (two duplicate edges and their siblings), since
        // only two of these four edges will be used in the final output.
        //
        // Furthermore, since the options REQUIRE and CREATE guarantee that all
        // edges will have siblings, S2Builder implements these options for
        // undirected edges by discarding half of the edges in each direction and
        // changing the edge_type() to EdgeType.DIRECTED.  For example, two
        // undirected input edges between vertices A and B would first be converted
        // into two directed edges in each direction, and then one edge of each pair
        // would be discarded leaving only one edge in each direction.
        //
        // Degenerate edges are considered not to have siblings.  If such edges are
        // present, they are passed through unchanged by SiblingPairs.DISCARD.  For
        // SiblingPairs.REQUIRE or SiblingPairs.CREATE with undirected edges, the
        // number of copies of each degenerate edge is reduced by a factor of two.
        //
        // Any of the options that discard edges (DISCARD, DISCARD_EXCESS, and
        // REQUIRE/CREATE in the case of undirected edges) have the side effect that
        // when duplicate edges are present, all of the corresponding edge labels
        // are merged together and assigned to the remaining edges.  (This avoids
        // the problem of having to decide which edges are discarded.)  Note that
        // this merging takes place even when all copies of an edge are kept.  For
        // example, consider the graph {AB1, AB2, AB3, BA4, CD5, CD6} (where XYn
        // denotes an edge from X to Y with label "n").  With SiblingPairs::DISCARD,
        // we need to discard one of the copies of AB.  But which one?  Rather than
        // choosing arbitrarily, instead we merge the labels of all duplicate edges
        // (even ones where no sibling pairs were discarded), yielding {AB123,
        // AB123, CD45, CD45} (assuming that duplicate edges are being kept).
        // Notice that the labels of duplicate edges are merged even if no siblings
        // were discarded (such as CD5, CD6 in this example), and that this would
        // happen even with duplicate degenerate edges (e.g. the edges EE7, EE8).
        //
        // DEFAULT: SiblingPairs.KEEP
        public enum SiblingPairs : byte { DISCARD, DISCARD_EXCESS, KEEP, REQUIRE, CREATE }
    }

    private class SiteIdsComp : IComparer<SiteId>
    {
        private readonly S2Point X;
        private readonly List<S2Point> Sites;

        public SiteIdsComp(S2Point x, List<S2Point> sites) { X = x; Sites = sites; }

        public int Compare(SiteId i, SiteId j)
        {
            return S2Pred.CompareDistances(X, Sites[i], Sites[j]);
        }
    }

    // MemoryTracker is a helper class to measure S2Builder memory usage.  It is
    // based on a detailed analysis of the data structures used.  This approach
    // is fragile because the memory tracking code needs to be updated whenever
    // S2Builder is modified, however S2Builder has been quite stable and this
    // approach allows the memory usage to be measured quite accurately.
    //
    // CAVEATS:
    //
    //  - Does not track memory used by edge labels.  (It is tricky to do this
    //    accurately because they are stored in an IdSetLexicon, and labels
    //    are typically a tiny fraction of the total space used.)
    //
    //  - Does not track memory used to represent layers internally.  (The
    //    number of layers is typically small compared to the numbers of
    //    vertices and edges, and the amount of memory used by the Layer and
    //    IsFullPolygonPredicate objects is difficult to measure.)
    //
    //  - Does not track memory used by the output layer Build() methods.  (This
    //    includes both temporary space, e.g. due to calling S2Builder::Graph
    //    methods, and also any geometric objects created by these layers.)
    public class MemoryTracker : S2MemoryTracker.Client
    {
        public MemoryTracker() : base() { }
        public MemoryTracker(S2MemoryTracker tracker) : base(tracker) { }

        // Called to track memory used to store the set of sites near a given edge.
        public bool TallyEdgeSites(List<SiteId> sites)
        {
            var size = GetCompactArrayAllocBytes(sites);
            edge_sites_bytes_ += size;
            return Tally(size);
        }
        // Ensures that "sites" contains space for at least one more edge site.
        public bool ReserveEdgeSite(List<SiteId> sites)
        {
            var new_size = sites.Count + 1;
            if (new_size <= sites.Capacity) return true;
            var old_bytes = GetCompactArrayAllocBytes(sites);
            sites.Capacity = new_size;
            var added_bytes = GetCompactArrayAllocBytes(sites) - old_bytes;
            edge_sites_bytes_ += added_bytes;
            return Tally(added_bytes);
        }
        // Releases and tracks the memory used to store nearby edge sites.
        public bool ClearEdgeSites(List<List<SiteId>> edge_sites)
        {
            Tally(-edge_sites_bytes_);
            edge_sites_bytes_ = 0;
            return Clear(edge_sites);
        }

        // Called when a site is added to the S2PointIndex.
        public bool TallyIndexedSite()
        {
            // S2PointIndex stores its data in a btree.  In general btree nodes are only
            // guaranteed to be half full, but in our case all nodes are full except for
            // the rightmost node at each btree level because the values are added in
            // sorted order.
            var delta_bytes = GetBtreeMinBytesPerEntry<
                List<KeyData2<S2CellId, S2Point, SiteId>>>();
            site_index_bytes_ += delta_bytes;
            return Tally(delta_bytes);
        }
        // Corrects the approximate S2PointIndex memory tracking done above.
        public bool FixSiteIndexTally(S2PointIndex<SiteId> index)
        {
            var delta_bytes = index.SpaceUsed() - site_index_bytes_;
            site_index_bytes_ += delta_bytes;
            return Tally(delta_bytes);
        }
        // Tracks memory due to destroying the site index.
        public bool DoneSiteIndex(/*S2PointIndex<SiteId> index*/)
        {
            Tally(-site_index_bytes_);
            site_index_bytes_ = 0;
            return Ok();
        }

        // Called to indicate that edge simplification was requested.
        public bool TallySimplifyEdgeChains(
                        List<List<InputVertexId>> site_vertices,
                        List<List<Edge>> layer_edges)
        {
            if (!IsActive()) return true;

            // The simplify_edge_chains() option uses temporary memory per site
            // (output vertex) and per output edge, as outlined below.
            //
            // Per site:
            //  vector<compact_array<InputVertexId>> site_vertices;  // BuildLayerEdges
            //   - compact_array non-inlined space is tallied separately
            //  vector<bool> is_interior_;  // EdgeChainSimplifier
            //  Graph::VertexInMap in_;     // EdgeChainSimplifier
            //  Graph::VertexOutMap out_;   // EdgeChainSimplifier
            var kTempPerSite =
                Marshal.SizeOf(typeof(List<InputVertexId>)) + sizeof(bool) + 2 * sizeof(EdgeId);

            // Per output edge:
            //  vector<bool> used_;                              // EdgeChainSimplifier
            //  Graph::VertexInMap in_;                          // EdgeChainSimplifier
            //  vector<Edge> merged_edges;                       // SimplifyEdgeChains
            //  vector<InputEdgeIdSetId> merged_input_edge_ids;  // SimplifyEdgeChains
            //  vector<int> merged_edge_layers;                  // SimplifyEdgeChains
            //  vector<Edge> new_edges_;                         // EdgeChainSimplifier
            //  vector<InputEdgeIdSetId> new_input_edge_ids_;    // EdgeChainSimplifier
            //  vector<int> new_edge_layers_;                    // EdgeChainSimplifier
            //
            // Note that the temporary vector<LayerEdgeId> in MergeLayerEdges() does not
            // affect peak usage.
            var kTempPerEdge =
                sizeof(bool) + sizeof(EdgeId) + 2 * Marshal.SizeOf(typeof(Edge)) +
                2 * sizeof(InputEdgeIdSetId) + 2 * sizeof(int);
            var simplify_bytes = site_vertices.Count * kTempPerSite;
            foreach (var array in site_vertices) {
                simplify_bytes += GetCompactArrayAllocBytes(array);
            }
            foreach (var edges in layer_edges) {
                simplify_bytes += edges.Count * kTempPerEdge;
            }
            return TallyTemp(simplify_bytes);
        }

        // Tracks the temporary memory used by Graph::FilterVertices.
        public bool TallyFilterVertices(int num_sites,
                             List<List<Edge>> layer_edges)
        {
            if (!IsActive()) return true;

            // Vertex filtering (see BuildLayers) uses temporary space of one VertexId
            // per Voronoi site plus 2 VertexIds per layer edge, plus space for all the
            // vertices after filtering.
            //
            //  vector<VertexId> *tmp;      // Graph::FilterVertices
            //  vector<VertexId> used;      // Graph::FilterVertices
            const Int64 kTempPerSite = sizeof(VertexId);
            const Int64 kTempPerEdge = 2 * sizeof(VertexId);

            int max_layer_edges = 0;
            foreach (var edges in layer_edges) {
                max_layer_edges = Math.Max(max_layer_edges, edges.Count);
            }
            filter_vertices_bytes_ = (num_sites * kTempPerSite +
                                      max_layer_edges * kTempPerEdge);
            return Tally(filter_vertices_bytes_);
        }
        public bool DoneFilterVertices()
        {
            Tally(-filter_vertices_bytes_);
            filter_vertices_bytes_ = 0;
            return Ok();
        }

        // The amount of non-inline memory used to store edge sites.
        private Int64 edge_sites_bytes_ = 0;

        // The amount of memory used by the S2PointIndex for sites.
        private Int64 site_index_bytes_ = 0;

        // The amount of temporary memory used by Graph::FilterVertices().
        private Int64 filter_vertices_bytes_ = 0;
}
}
