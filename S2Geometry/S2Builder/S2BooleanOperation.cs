// Boolean operations are implemented by constructing the boundary of the
// result and then using S2Builder to assemble the edges.  The boundary is
// obtained by clipping each of the two input regions to the interior or
// exterior of the other region.  For example, to compute the union of A and
// B, we clip the boundary of A to the exterior of B and the boundary of B to
// the exterior of A; the resulting set of edges defines the union of the two
// regions.
//
// We use exact predicates, but inexactructions (e.g. computing the
// intersection point of two edges).  Nevertheless, the following algorithm is
// guaranteed to be 100% robust, in that the computed boundary stays within a
// small tolerance (snap_radius + S2EdgeCrossings.kIntersectionError) of the exact
// result, and also preserves the correct topology (i.e., no crossing edges).
//
// Unfortunately this robustness cannot quite be achieved using the strategy
// outlined above (clipping the two input regions and assembling the
// resulting edges).  Since computed intersection points are not exact, the
// input geometry passed to S2Builder might contain self-intersections, and
// these self-intersections cannot be eliminated reliably by snap rounding.
//
// So instead, we pass S2Builder the entire set of input edges where at least
// some portion of each edge belongs to the output boundary.  We allow
// S2Builder to compute the intersection points and snap round the edges
// (which it does in a way that is guaranteed to preserve the input topology).
// Then once this is finished, we remove the portions of each edge that would
// have been clipped if we had done the clipping first.  This step only
// involves deciding whether to keep or discard each edge in the output, since
// all intersection points have already been resolved, and therefore there is
// no risk of creating new self-intersections.
//
// This is implemented using the following classes:
//
//  - S2BooleanOperation.Impl: the top-level class that clips each of
//                              the two regions to the other region.
//
//  - CrossingProcessor: a class that processes edge crossings and maintains
//                       the necessary state in order to clip the boundary
//                       of one region to the interior or exterior of the
//                       other region.
//
//  - EdgeClippingLayer: an S2Builder.Layer that removes graph edges that
//                       correspond to clipped portions of input edges, and
//                       passes the result to another layer for assembly.
//
//  - GraphEdgeClipper: a helper class that does the actual work of the
//                      EdgeClippingLayer.


// This class implements boolean operations (intersection, union, difference,
// and symmetric difference) for regions whose boundaries are defined by
// geodesic edges.
//
// S2BooleanOperation operates on exactly two input regions at a time.  Each
// region is represented as an S2ShapeIndex and may contain any number of
// points, polylines, and polygons.  The region is essentially the union of
// these objects, except that polygon interiors must be disjoint from all
// other geometry (including other polygon interiors).  If the input geometry
// for a region does not meet this condition, it can be normalized by
// computing its union first.  Duplicate polygon edges are not allowed (even
// among different polygons), however polylines may have duplicate edges and
// may even be self-intersecting.  Note that points or polylines are allowed
// to coincide with the boundaries of polygons.
//
// Degeneracies are fully supported.  Supported degeneracy types include the
// following:
//
//  - Point polylines consisting of a single degenerate edge AA.
//
//  - Point loops consisting of a single vertex A.  Such loops may represent
//    either shells or holes according to whether the loop adds to or
//    subtracts from the surrounding region of the polygon.
//
//  - Sibling edge pairs of the form {AB, BA}.  Such sibling pairs may
//    represent either shells or holes according to whether they add to or
//    subtract from the surrounding region.  The edges of a sibling pair may
//    belong to the same polygon loop (e.g. a loop AB) or to different polygon
//    loops or polygons (e.g. the polygons {ABC, CBD}).
//
// A big advantage of degeneracy support is that geometry may be simplified
// without completely losing small details.  For example, if a polygon
// representing a land area with many lakes and rivers is simplified using a
// tolerance of 1 kilometer, every water feature in the input is guaranteed to
// be within 1 kilometer of some water feature in the input (even if some
// lakes and rivers are merged and/or reduced to degenerate point or sibling
// edge pair holes).  Mathematically speaking, degeneracy support allows
// geometry to be simplified while guaranteeing that the Hausdorff distance
// betweeen the boundaries of the original and simplified geometries is at
// most the simplification tolerance.  It also allows geometry to be
// simplified without changing its dimension, thus preserving boundary
// semantics.  (Note that the boundary of a polyline ABCD is {A,D}, whereas
// the boundary of a degenerate shell ABCDCB is its entire set of vertices and
// edges.)
//
// Points and polyline edges are treated as multisets: if the same point or
// polyline edge appears multiple times in the input, it will appear multiple
// times in the output.  For example, the union of a point with an identical
// point consists of two points.  This feature is useful for modeling large
// sets of points or polylines as a single region while maintaining their
// distinct identities, even when the points or polylines intersect each
// other.  It is also useful for reconstructing polylines that loop back on
// themselves (e.g., time series such as GPS tracks).  If duplicate geometry
// is not desired, it can easily be removed by choosing the appropriate
// S2Builder output layer options.
//
// Self-intersecting polylines can be manipulated without materializing new
// vertices at the self-intersection points.  This feature is important when
// processing polylines with large numbers of self-intersections such as GPS
// tracks (e.g., consider the path of a race car in the Indy 500).
//
// Polylines are always considered to be directed.  Polyline edges between the
// same pair of vertices are defined to intersect even if the two edges are in
// opposite directions.  (Undirected polylines can be modeled by specifying
// GraphOptions.EdgeType.UNDIRECTED in the S2Builder output layer.)
//
// The output of each operation is sent to an S2Builder.Layer provided by the
// client.  This allows clients to build any representation of the geometry
// they choose.  It also allows the client to do additional postprocessing of
// the output before building data structures; for example, the client can
// easily discard degeneracies or convert them to another data type.
//
// The boundaries of polygons and polylines can be modeled as open, semi-open,
// or closed.  Polyline boundaries are controlled by the PolylineModel class,
// whose options are as follows:
//
//  - In the OPEN model, polylines do not contain their first or last vertex
//    except for one special case: namely, if the polyline forms a loop and
//    the polyline_loops_have_boundaries() option is set to false, then the
//    first/last vertex is contained.
//
//  - In the SEMI_OPEN model, polylines contain all vertices except the last.
//    Therefore if one polyline starts where another polyline stops, the two
//    polylines do not intersect.
//
//  - In the CLOSED model, polylines contain all of their vertices.
//
// When multiple polylines are present, they are processed independently and
// have no effect on each other.  For example, in the OPEN boundary model the
// polyline ABC contains the vertex B, while set of polylines {AB, BC} does
// not.  (If you want to treat the polylines as a union instead, with
// boundaries merged according to the "mod 2" rule, this can be achieved by
// reassembling the edges into maximal polylines using S2PolylineVectorLayer
// with EdgeType.UNDIRECTED, DuplicateEdges.MERGE, and PolylineType.WALK.)
//
// Polygon boundaries are controlled by the PolygonModel class, which has the
// following options:
//
//  - In the OPEN model, polygons do not contain their vertices or edges.
//    This implies that a polyline that follows the boundary of a polygon will
//    not intersect it.
//
//  - In the SEMI_OPEN model, polygon point containment is defined such that
//    if several polygons tile the region around a vertex, then exactly one of
//    those polygons contains that vertex.  Similarly polygons contain all of
//    their edges, but none of their reversed edges.  This implies that a
//    polyline and polygon edge with the same endpoints intersect if and only
//    if they are in the same direction.  (This rule ensures that if a
//    polyline is intersected with a polygon and its complement, the two
//    resulting polylines do not have any edges in common.)
//
//  - In the CLOSED model, polygons contain all their vertices, edges, and
//    reversed edges.  This implies that a polyline that shares an edge (in
//    either direction) with a polygon is defined to intersect it.  Similarly,
//    this is the only model where polygons that touch at a vertex or along an
//    edge intersect.
//
// Note that PolylineModel and PolygonModel are defined as separate classes in
// order to allow for possible future extensions.
//
// Operations between geometry of different dimensions are defined as follows:
//
//  - For UNION, the higher-dimensional shape always wins.  For example the
//    union of a closed polygon A with a polyline B that coincides with the
//    boundary of A consists only of the polygon A.
//
//  - For INTERSECTION, the lower-dimensional shape always wins.  For example,
//    the intersection of a closed polygon A with a point B that coincides
//    with a vertex of A consists only of the point B.
//
//  - For DIFFERENCE, higher-dimensional shapes are not affected by
//    subtracting lower-dimensional shapes.  For example, subtracting a point
//    or polyline from a polygon A yields the original polygon A.  This rule
//    exists because in general, it is impossible to represent the output
//    using the specified boundary model(s).  (Consider subtracting one vertex
//    from a PolylineModel.CLOSED polyline, or subtracting one edge from a
//    PolygonModel.CLOSED polygon.)  If you want to perform operations like
//    this, consider representing all boundaries explicitly (topological
//    boundaries) using OPEN boundary models.  Another option for polygons is
//    to subtract a degenerate loop, which yields a polygon with a degenerate
//    hole (see S2LaxPolygonShape).
//
// Note that in the case of Precision.EXACT operations, the above remarks
// only apply to the output before snapping.  Snapping may cause nearby
// distinct edges to become coincident, e.g. a polyline may become coincident
// with a polygon boundary.  However also note that S2BooleanOperation is
// perfectly happy to accept such geometry as input.
//
// Note the following differences between S2BooleanOperation and the similar
// S2MultiBooleanOperation class:
//
//  - S2BooleanOperation operates on exactly two regions at a time, whereas
//    S2MultiBooleanOperation operates on any number of regions.
//
//  - S2BooleanOperation is potentially much faster when the input is already
//    represented as S2ShapeIndexes.  The algorithm is output sensitive and is
//    often sublinear in the input size.  This can be a big advantage if, say,
//
//  - S2BooleanOperation supports exact predicates and the corresponding
//    exact operations (i.e., operations that are equivalent to computing the
//    exact result and then snap rounding it).
//
//  - S2MultiBooleanOperation has better error guarantees when there are many
//    regions, since it requires only one snapping operation for any number of
//    input regions.
//
// Example usage:
//   S2ShapeIndex a, b;  // Input geometry, e.g. containing polygons.
//   S2Polygon polygon;  // Output geometry.
//   S2BooleanOperation.Options options;
//   options.set_snap_function(snap_function);
//   S2BooleanOperation op(S2BooleanOperation.OpType.INTERSECTION,
//                         new S2PolygonLayer(polygon),
//                         options);
//   S2Error error;
//   if (!op.Build(a, b, error)) {
//     S2_LOG(ERROR) << error;
//     ...
//   }
//
// If the output includes objects of different dimensions, they can be
// assembled into different layers with code like this:
//
//   S2Point[] points;
//   S2Polyline[] polylines;
//   S2Polygon polygon;
//   S2BooleanOperation op(
//       S2BooleanOperation.OpType.UNION,
//       new s2builderutil.PointVectorLayer(points),
//       new s2builderutil.S2PolylineVectorLayer(polylines),
//       new S2PolygonLayer(polygon));

namespace S2Geometry;

using System.Runtime.InteropServices;
using S2BuilderUtil;
using static S2Builder;
using static S2Builder.GraphOptions;
using ShapeEdge = S2ShapeUtil.ShapeEdge;
using ShapeEdgeId = S2ShapeUtil.ShapeEdgeId;
using SourceIdMap = Dictionary<SourceId, Int32>; // gtl.btree_map

public class S2BooleanOperation
{
    // A collection of special Int32s that allow the GraphEdgeClipper state
    // modifications to be inserted into the list of edge crossings.
    public const InputEdgeId kSetInside = -1;
    public const InputEdgeId kSetInvertB = -2;
    public const InputEdgeId kSetReverseA = -3;

    private readonly Options Options_;

    public OpType OpType_ { get; }

    // The input regions.
    private readonly S2ShapeIndex[] regions_ = new S2ShapeIndex[2];

    // The output consists either of zero layers, one layer, or three layers.
    private readonly List<Layer> layers_;

    // The following field is set if and only if there are no output layers.
    private readonly Action<bool>? result_empty_;

    /// <summary>
    /// Specifies that the output boundary edges should be sent to a single
    /// S2Builder layer.  This version can be used when the dimension of the
    /// output geometry is known (e.g., intersecting two polygons to yield a
    /// third polygon).
    /// </summary>
    public S2BooleanOperation(OpType op_type, Layer layer, Options? options = null)
    {
        Options_ = options ?? new Options();
        OpType_ = op_type;
        result_empty_ = null;
        layers_ = new List<Layer>() { layer };
    }

    // Specifies that the output boundary edges should be sent to three
    // different layers according to their dimension.  Points (represented by
    // degenerate edges) are sent to layer 0, polyline edges are sent to
    // layer 1, and polygon edges are sent to layer 2.
    //
    // The dimension of an edge is defined as the minimum dimension of the two
    // input edges that produced it.  For example, the intersection of two
    // crossing polyline edges is a considered to be a degenerate polyline
    // rather than a point, so it is sent to layer 1.  Clients can easily
    // reclassify such polylines as points if desired, but this rule makes it
    // easier for clients that want to process point, polyline, and polygon
    // inputs differently.
    //
    // The layers are always built in the order 0, 1, 2, and all arguments to
    // the Build() calls are guaranteed to be valid until the last call returns.
    // All Graph objects have the same set of vertices and the same lexicon
    // objects, in order to make it easier to write classes that process all the
    // edges in parallel.
    public S2BooleanOperation(OpType op_type, List<Layer> layers, Options? options = null)
    {
        Options_ = options ?? new Options();
        OpType_ = op_type;
        layers_ = layers; result_empty_ = null;
    }

    // Specifies that "result_empty" should be set to indicate whether the exact
    // result of the operation is empty.  This constructor is used to efficiently
    // test boolean relationships (see IsEmpty above).
    private S2BooleanOperation(OpType op_type, Action<bool> result_empty, Options? options = null)
    {
        Options_ = options ?? new Options();
        OpType_ = op_type;
        result_empty_ = result_empty;
    }

    // Executes the given operation.  Returns true on success, and otherwise
    // sets "error" appropriately.  (This class does not generate any errors
    // itself, but the S2Builder.Layer might.)
    public bool Build(S2ShapeIndex a, S2ShapeIndex b, out S2Error error)
    {
        regions_[0] = a;
        regions_[1] = b;
        return new Impl(this).Build(out error);
    }

    // Convenience method that returns true if the result of the given operation
    // is empty.
    public static bool IsEmpty(OpType op_type, S2ShapeIndex a, S2ShapeIndex b, Options? options = null)
    {
        bool result_empty = true;
        var op = new S2BooleanOperation(op_type, (res) => result_empty = res, options ?? new Options());
        op.Build(a, b, out S2Error error);
        System.Diagnostics.Debug.Assert(error.IsOk());
        return result_empty;
    }

    // Convenience method that returns true if A intersects B.
    public static bool Intersects(S2ShapeIndex a, S2ShapeIndex b, Options? options = null)
    {
        return !IsEmpty(OpType.INTERSECTION, b, a, options ?? new Options());
    }

    // Convenience method that returns true if A contains B, i.e., if the
    // difference (B - A) is empty.
    public static bool Contains(S2ShapeIndex a, S2ShapeIndex b, Options? options = null)
    {
        return IsEmpty(OpType.DIFFERENCE, b, a, options ?? new Options());
    }

    // Convenience method that returns true if the symmetric difference of A and
    // B is empty.  (Note that A and B may still not be identical, e.g. A may
    // contain two copies of a polyline while B contains one.)
    public static bool Equals(S2ShapeIndex a, S2ShapeIndex b, Options? options = null)
    {
        return IsEmpty(OpType.SYMMETRIC_DIFFERENCE, b, a, options ?? new Options());
    }

    // The supported operation types.
    public enum OpType : byte
    {
        None,
        UNION,                // Contained by either region.
        INTERSECTION,         // Contained by both regions.
        DIFFERENCE,           // Contained by the first region but not the second.
        SYMMETRIC_DIFFERENCE  // Contained by one region but not the other.
    }

    // Defines whether polygons are considered to contain their vertices and/or
    // edges (see definitions above).
    public enum PolygonModel : byte { OPEN, SEMI_OPEN, CLOSED }

    // Defines whether polylines are considered to contain their endpoints
    // (see definitions above).
    public enum PolylineModel : byte { OPEN, SEMI_OPEN, CLOSED }

    // With Precision.EXACT, the operation is evaluated using the exact input
    // geometry.  Predicates that use this option will produce exact results;
    // for example, they can distinguish between a polyline that barely
    // intersects a polygon from one that barely misses it.  Constructive
    // operations (ones that yield new geometry, as opposed to predicates) are
    // implemented by computing the exact result and then snap rounding it
    // according to the given snap_function() (see below).  This is as close as
    // it is possible to get to the exact result while requiring that vertex
    // coordinates have type "double".
    //
    // With Precision.SNAPPED, the input regions are snapped together *before*
    // the operation is evaluated.  So for example, two polygons that overlap
    // slightly will be treated as though they share a common boundary, and
    // similarly two polygons that are slightly separated from each other will
    // be treated as though they share a common boundary.  Snapped results are
    // useful for dealing with points, since in S2 the only points that lie
    // exactly on a polyline or polygon edge are the endpoints of that edge.
    //
    // Conceptually, the difference between these two options is that with
    // Precision.SNAPPED, the inputs are snap rounded (together), whereas with
    // Precision.EXACT only the result is snap rounded.
    public enum Precision : byte { EXACT, SNAPPED }

    public class Options
    {
        public Options()
        {
            SnapFunction_ = new IdentitySnapFunction(S1Angle.Zero);
        }

        // Convenience constructor that calls set_snap_function().
        public Options(SnapFunction snap_function)
        {
            SnapFunction_ = (SnapFunction)snap_function.CustomClone();
        }

        public Options(Options options)
        {
            SnapFunction_ = (SnapFunction)options.SnapFunction_.CustomClone();
            PolygonModel_ = options.PolygonModel_;
            PolylineModel_ = options.PolylineModel_;
            PolylineLoopsHaveBoundaries = options.PolylineLoopsHaveBoundaries;
            split_all_crossing_polyline_edges_ = options.split_all_crossing_polyline_edges_;
            Precision = options.Precision;
            ConservativeOutput = options.ConservativeOutput;
            SourceIdLexicon_ = options.SourceIdLexicon_;
            memory_tracker_ = options.memory_tracker_;
        }

        // Specifies the function to be used for snap rounding.
        //
        // DEFAULT: s2builderutil.IdentitySnapFunction(S1Angle.Zero())
        //  - This does no snapping and preserves all input vertices exactly unless
        //    there are crossing edges, in which case the snap radius is increased
        //    to the maximum intersection point error (S2EdgeCrossings.kIntersectionError).
        public SnapFunction SnapFunction_ { get; set; }

        // Defines whether polygons are considered to contain their vertices
        // and/or edges (see comments above).
        //
        // DEFAULT: PolygonModel.SEMI_OPEN
        public PolygonModel PolygonModel_ { get; set; } = PolygonModel.SEMI_OPEN;

        // Defines whether polylines are considered to contain their vertices (see
        // comments above).
        //
        // DEFAULT: PolylineModel.CLOSED
        public PolylineModel PolylineModel_ { get; set; } = PolylineModel.CLOSED;

        // Specifies whether a polyline loop is considered to have a non-empty
        // boundary.  By default this option is true, meaning that even if the
        // first and last vertices of a polyline are the same, the polyline is
        // considered to have a well-defined "start" and "end".  For example, if
        // the polyline boundary model is OPEN then the polyline loop would not
        // include the start/end vertices.  These are the best semantics for most
        // applications, such as GPS tracks or road network segments.
        //
        // If the polyline forms a loop and this option is set to false, then
        // instead the first and last vertices are considered to represent a
        // single vertex in the interior of the polyline.  In this case the
        // boundary of the polyline is empty, meaning that the first/last vertex
        // will be contained by the polyline even if the boundary model is OPEN.
        // (Note that this option also has a small effect on the CLOSED boundary
        // model, because the first/last vertices of a polyline loop are
        // considered to represent one vertex rather than two.)
        //
        // The main reason for this option is to implement the "mod 2 union"
        // boundary semantics of the OpenGIS Simple Features spec.  This can be
        // achieved by making sure that all polylines are constructed using
        // S2Builder.Graph.PolylineType.WALK (which ensures that all polylines
        // are as long as possible), and then setting this option to false.
        //
        // DEFAULT: true
        public bool PolylineLoopsHaveBoundaries { get; set; } = true;

        // Specifies that a new vertex should be added whenever a polyline edge
        // crosses another polyline edge.  Note that this can cause the size of
        // polylines with many self-intersections to increase quadratically.
        //
        // If false, new vertices are added only when a polyline from one input
        // region cross a polyline from the other input region.  This allows
        // self-intersecting input polylines to be modified as little as possible.
        //
        // DEFAULT: false
        public bool split_all_crossing_polyline_edges_ { get; set; } = false;

        // Specifies whether the operation should use the exact input geometry
        // (Precision.EXACT), or whether the two input regions should be snapped
        // together first (Precision.SNAPPED).
        //
        // DEFAULT: Precision.EXACT
        public Precision Precision { get; } = Precision.EXACT;
        // void set_precision(Precision precision);

        // If true, the input geometry is interpreted as representing nearby
        // geometry that has been snapped or simplified.  It then outputs a
        // conservative result based on the value of polygon_model() and
        // polyline_model().  For the most part, this only affects the handling of
        // degeneracies.
        //
        // - If the model is OPEN, the result is as open as possible.  For
        //   example, the intersection of two identical degenerate shells is empty
        //   under PolygonModel.OPEN because they could have been disjoint before
        //   snapping.  Similarly, two identical degenerate polylines have an
        //   empty intersection under PolylineModel.OPEN.
        //
        // - If the model is CLOSED, the result is as closed as possible.  In the
        //   case of the DIFFERENCE operation, this is equivalent to evaluating
        //   A - B as Closure(A) - Interior(B).  For other operations, it affects
        //   only the handling of degeneracies.  For example, the union of two
        //   identical degenerate holes is empty under PolygonModel.CLOSED
        //   (i.e., the hole disappears) because the holes could have been
        //   disjoint before snapping.
        //
        // - If the model is SEMI_OPEN, the result is as degenerate as possible.
        //   New degeneracies will not be created, but all degeneracies that
        //   coincide with the opposite region's boundary are retained unless this
        //   would cause a duplicate polygon edge to be created.  This model is
        //   is very useful for working with input data that has both positive and
        //   negative degeneracies (i.e., degenerate shells and holes).
        //
        // DEFAULT: false
        public bool ConservativeOutput { get; } = false;
        // void set_conservative_output(bool conservative);

        // Specifies that internal memory usage should be tracked using the given
        // S2MemoryTracker.  If a memory limit is specified and more more memory
        // than this is required then an error will be returned.  Example usage:
        //
        //   S2MemoryTracker tracker;
        //   tracker.set_limit(500 << 20);  // 500 MB
        //   S2BooleanOperation.Options options;
        //   options.set_memory_tracker(&tracker);
        //   S2BooleanOperation op(..., options);
        //   ...
        //   S2Error error;
        //   if (!op.Build(..., &error)) {
        //     if (error.code() == S2Error.RESOURCE_EXHAUSTED) {
        //       S2_LOG(ERROR) << error;  // Memory limit exceeded
        //     }
        //   }
        //
        // CAVEATS:
        //
        //  - Memory used by the input S2ShapeIndexes and the output S2Builder
        //    layers is not counted towards the total.
        //
        //  - While memory tracking is reasonably complete and accurate, it does
        //    not account for every last byte.  It is intended only for the
        //    purpose of preventing clients from running out of memory.
        //
        // DEFAULT: nullptr (memory tracking disabled)
        public S2MemoryTracker? memory_tracker_ { get; set; } = null;

        // If specified, then each output edge will be labelled with one or more
        // SourceIds indicating which input edge(s) it corresponds to.  This
        // can be useful if your input geometry has additional data that needs to
        // be propagated from the input to the output (e.g., elevations).
        //
        // You can access the labels by using an S2Builder.Layer type that
        // supports labels, such as S2PolygonLayer.  The layer outputs a
        // "label_set_lexicon" and an "label_set_id" for each edge.  You can then
        // look up the source information for each edge like this:
        //
        // for (Int32 label : label_set_lexicon.id_set(label_set_id)) {
        //   SourceId src = source_id_lexicon.value(label);
        //   // region_id() specifies which S2ShapeIndex the edge is from (0 or 1).
        //   DoSomething(src.region_id(), src.shape_id(), src.edge_id());
        // }
        //
        // DEFAULT: null
        public SourceIdLexicon? SourceIdLexicon_ { get; } = null;
    }

    private class Impl
    {
        public Impl(S2BooleanOperation op)
        {
            op_ = op;
            index_crossings_first_region_id_ = -1;
            tracker_ = new(op.Options_.memory_tracker_);
        }

        public void DoBuild(out S2Error error)
        {
            if (!tracker_.Ok())
            {
                error = S2Error.OK;
                return;
            }
            builder_options_ = new(op_.Options_.SnapFunction_)
            {
                IntersectionTolerance = S2.kIntersectionErrorS1Angle,
                MemoryTracker = tracker_.Tracker,
            };
            if (op_.Options_.split_all_crossing_polyline_edges_)
            {
                builder_options_.SplitCrossingEdges = true;
            }
            // TODO(ericv): Ideally idempotent() should be true, but existing clients
            // expect vertices closer than the full "snap_radius" to be snapped.
            builder_options_.Idempotent = false;

            error = S2Error.OK;
            if (IsBooleanOutput())
            {
                // BuildOpType() returns true if and only if the result has no edges.
                /*var g = new Graph();*/  // Unused by IsFullPolygonResult() implementation.
                op_.result_empty_(BuildOpType(op_.OpType_) && !IsFullPolygonResult(/*g,*/ out error));
                return;
            }
            builder_ = new S2Builder(builder_options_);
            builder_.StartLayer(new EdgeClippingLayer(
                op_.layers_, input_dimensions_, input_crossings_, tracker_));

            // Add a predicate that decides whether a result with no polygon edges should
            // be interpreted as the empty polygon or the full polygon.
            builder_.AddIsFullPolygonPredicate((Graph g, out S2Error error) => IsFullPolygonResult(/*g,*/ out error));
            BuildOpType(op_.OpType_);

            // Release memory that is no longer needed.
            if (!tracker_.Clear(index_crossings_)) return;
            builder_.Build(out error);
        }

        public bool Build(out S2Error error)
        {
            // This wrapper ensures that memory tracking errors are reported.
            DoBuild(out error);
            if (!tracker_.Ok()) error = tracker_.Error();
            return error.IsOk();
        }

        // An IndexCrossing represents a pair of intersecting S2ShapeIndex edges
        // ("a_edge" and "b_edge").  We store all such intersections because the
        // algorithm needs them twice, once when processing the boundary of region A
        // and once when processing the boundary of region B.
        public record struct IndexCrossing : IComparable<IndexCrossing>
        {
            public Edge A { get; set; }
            public Edge B { get; set; }

            // True if S2EdgeCrossings.CrossingSign(a_edge, b_edge) > 0.
            public bool IsInteriorCrossing { get; set; }

            // True if "a_edge" crosses "b_edge" from left to right.  Undefined if
            // is_interior_crossing is false.
            public bool LeftToRight { get; set; }

            // Equal to S2EdgeCrossings.VertexCrossing(a_edge, b_edge).  Undefined if "a_edge" and
            // "b_edge" do not share exactly one vertex or either edge is degenerate.
            public bool IsVertexCrossing { get; set; }

            // All flags are "false" by default.
            public IndexCrossing(Edge a, Edge b)
            {
                A = a; B = b; IsInteriorCrossing = false;
                LeftToRight = false; IsVertexCrossing = false;
            }

            public override int GetHashCode() => HashCode.Combine(A, B);
            public bool Equals(IndexCrossing other) => A == other.A && B == other.B;

            public int CompareTo(IndexCrossing other)
            {
                var c = A.CompareTo(other.A);
                if (c != 0)
                    return c;

                return B.CompareTo(other.B);
            }
            public static bool operator <(IndexCrossing x, IndexCrossing y) => x.CompareTo(y) < 0;
            public static bool operator >(IndexCrossing x, IndexCrossing y) => x.CompareTo(y) > 0;
            public static bool operator <=(IndexCrossing x, IndexCrossing y) => x.CompareTo(y) <= 0;
            public static bool operator >=(IndexCrossing x, IndexCrossing y) => x.CompareTo(y) >= 0;
        }

        public class MemoryTracker : S2MemoryTracker.Client
        {
            public MemoryTracker(S2MemoryTracker tracker) : base(tracker) { }

            // Used to track memory used by CrossingProcessor.source_id_map_.  (The
            // type is a template parameter so that SourceIdMap can be private.)
            public bool TallySourceIdMap<T>(int num_entries) where T : SourceIdMap
            {
                Int64 delta_bytes = num_entries * GetBtreeMinBytesPerEntry<T>();
                source_id_map_bytes_ += delta_bytes;
                return Tally(delta_bytes);
            }

            // Used to clear CrossingProcessor.source_id_map_ and update the tracked
            // memory usage accordingly.
            public bool ClearSourceIdMap<T>(T source_id_map) where T : SourceIdMap
            {
                source_id_map.Clear();
                Tally(-source_id_map_bytes_);
                source_id_map_bytes_ = 0;
                return Ok();
            }

            // The amount of memory used by CrossingProcessor.source_id_map_.
            private Int64 source_id_map_bytes_ = 0;
        }

        // CrossingProcessor is a helper class that processes all the edges from one
        // region that cross a specific edge of the other region.  It outputs the
        // appropriate edges to an S2Builder, and outputs other information required
        // by GraphEdgeClipper to the given vectors.
        private class CrossingProcessor
        {
            // Prepares to build output for the given polygon and polyline boundary
            // models.  Edges are emitted to "builder", while other auxiliary data is
            // appended to the given vectors.
            //
            // If a predicate is being evaluated (i.e., we do not need to construct the
            // actual result), then "builder" and the various output vectors should all
            // be null.
            public CrossingProcessor(PolygonModel polygon_model, PolylineModel polyline_model,
                bool polyline_loops_have_boundaries, S2Builder builder,
                List<sbyte> input_dimensions, InputEdgeCrossings input_crossings,
                MemoryTracker tracker)
            {
                polygon_model_ = polygon_model; polyline_model_ = polyline_model;
                polyline_loops_have_boundaries_ = polyline_loops_have_boundaries;
                builder_ = builder; input_dimensions_ = input_dimensions;
                input_crossings_ = input_crossings; tracker_ = tracker;
                prev_inside_ = false;
            }

            private bool is_degenerate(ShapeEdgeId a_id)
            {
                return is_degenerate_hole_.ContainsKey(a_id);
            }

            // Starts processing edges from the given region.  "invert_a", "invert_b",
            // and "invert_result" indicate whether region A, region B, and/or the
            // result should be inverted, which allows operations such as union and
            // difference to be implemented.  For example, union is ~(~A & ~B).
            //
            // This method should be called in pairs, once to process the edges from
            // region A and once to process the edges from region B.
            public void StartBoundary(int a_region_id, bool invert_a, bool invert_b, bool invert_result)
            {
                a_region_id_ = a_region_id;
                b_region_id_ = 1 - a_region_id;
                invert_a_ = invert_a;
                invert_b_ = invert_b;
                invert_result_ = invert_result;
                is_union_ = invert_b && invert_result;

                // Specify to GraphEdgeClipper how these edges should be clipped.
                SetClippingState(kSetReverseA, invert_a != invert_result);
                SetClippingState(kSetInvertB, invert_b);
            }

            // Starts processing edges from the given shape.
            public void StartShape(S2Shape a_shape)
            {
                a_shape_ = a_shape;
                a_dimension_ = a_shape.Dimension();
            }

            // Starts processing edges from the given chain.
            public void StartChain(int chain_id, S2Shape.Chain chain, bool inside)
            {
                chain_id_ = chain_id;
                chain_start_ = chain.Start;
                chain_limit_ = chain.Start + chain.Length;
                Inside = inside;
                v0_emitted_max_edge_id_ = chain.Start - 1;  // No edges emitted yet.
                chain_v0_emitted_ = false;
            }

            // Processes the given edge "a_id".  "it" should be positioned to the set of
            // edges from the other region that cross "a_id" (if any).
            //
            // Supports "early exit" in the case of boolean results by returning false
            // as soon as the result is known to be non-empty.
            public bool ProcessEdge(Edge a_id, CrossingIterator it)
            {
                // chain_edge() is faster than edge() when there are multiple chains.
                var a = a_shape_.ChainEdge(chain_id_, a_id.EdgeId - chain_start_);
                if (a_dimension_ == 0)
                {
                    return ProcessEdge0(a_id, a, it);
                }
                else if (a_dimension_ == 1)
                {
                    return ProcessEdge1(a_id, a, it);
                }
                else
                {
                    System.Diagnostics.Debug.Assert(2 == a_dimension_);
                    return ProcessEdge2(a_id, a, it);
                }
            }

            // This method should be called after each pair of calls to StartBoundary.
            // (The only operation that processes more than one pair of boundaries is
            // SYMMETRIC_DIFFERENCE, which computes the union of A-B and B-A.)
            //
            // Resets the state of the CrossingProcessor.
            //
            // Translates the temporary representation of crossing edges (SourceId) into
            // the format expected by EdgeClippingLayer (Int32).
            public void DoneBoundaryPair()
            {
                tracker_.AddSpaceExact(input_crossings_, source_edge_crossings_.Count);
                if (!tracker_.TallySourceIdMap<SourceIdMap>(3)) return;

                // Add entries that translate the "special" crossings.
                source_id_map_[new SourceId(kSetInside)] = kSetInside;
                source_id_map_[new SourceId(kSetInvertB)] = kSetInvertB;
                source_id_map_[new SourceId(kSetReverseA)] = kSetReverseA;
                foreach (var tmp in source_edge_crossings_)
                {
                    System.Diagnostics.Debug.Assert(source_id_map_.ContainsKey(tmp.Item2.Item1));
                    var it = source_id_map_[tmp.Item2.Item1];
                    input_crossings_.Add((tmp.Item1, new CrossingInputEdge(it, tmp.Item2.Item2)));
                }
                tracker_.Clear(source_edge_crossings_);
                tracker_.ClearSourceIdMap(source_id_map_);
            }

            // Indicates whether the point being processed along the current edge chain
            // is in the polygonal interior of the opposite region, using semi-open
            // boundaries.  If "invert_b_" is true then this field is inverted.
            //
            // This value along with the set of incident edges can be used to compute
            // whether the opposite region contains this point under any of the
            // supported boundary models (PolylineModel.CLOSED, etc).
            public bool Inside { get; private set; }

            // PointCrossingResult describes the relationship between a point from region A
            // and a set of crossing edges from region B.  For example, "matches_polygon"
            // indicates whether a polygon vertex from region B matches the given point.
            private struct PointCrossingResult
            {
                public static readonly PointCrossingResult Zero = new(false, false, false);
                public PointCrossingResult(bool matches_point, bool matches_polyline, bool matches_polygon)
                {
                    MatchesPoint = matches_point; MatchesPolyline = matches_polyline; MatchesPolygon = matches_polygon;
                }
                // Note that "matches_polyline" is true only if the point matches a polyline
                // vertex of B *and* the polyline contains that vertex, whereas
                // "matches_polygon" is true if the point matches any polygon vertex.
                public bool MatchesPoint;     // Matches point.
                public bool MatchesPolyline;  // Matches contained polyline vertex.
                public bool MatchesPolygon;   // Matches polygon vertex.
            }
            // EdgeCrossingResult describes the relationship between an edge (a0, a1) from
            // region A and a set of crossing edges from region B.  For example,
            // "matches_polygon" indicates whether (a0, a1) matches a polygon edge from
            // region B.
            private struct EdgeCrossingResult
            {
                // These fields indicate that (a0, a1) exactly matches an edge of B.
                public bool matches_polyline = false;  // Matches polyline edge (either direction).

                // These fields indicate that a B polyline contains the degenerate polyline
                // (a0, a0) or (a1, a1).  (This is identical to whether the B polyline
                // contains the point a0 or a1 except when the B polyline is degenerate,
                // since a degenerate polyline VV contains itself in all boundary models but
                // contains the point V only in the CLOSED polyline model.)
                public bool A0MatchesPolyline = false;  // B polyline contains (a0, a0)
                public bool A1MatchesPolyline = false;  // B polyline contains (a1, a1)

                // These fields indicate that a vertex of (a0, a1) matches a polygon vertex
                // of B.  (Unlike with polylines, the polygon may not contain that vertex.)
                public bool a0_matches_polygon = false;   // a0 matches polygon vertex.
                public bool a1_matches_polygon = false;   // a1 matches polygon vertex.

                // When a0 != a1, the first two fields identify any B polygon edge that
                // exactly matches (a0, a1) or the sibling edge (a1, a0).  The third field
                // identifies any B polygon edge that exactly matches (a0, a0).
                public ShapeEdgeId polygon_match_id;  // B polygon edge that matches (a0, a1).
                public ShapeEdgeId sibling_match_id;  // B polygon edge that matches (a1, a0).
                public ShapeEdgeId a0_loop_match_id;  // B polygon edge that matches (a0, a0).

                // Convenience functions to test whether a matching edge was found.
                public bool matches_polygon() { return polygon_match_id.EdgeId >= 0; }
                public bool matches_sibling() { return sibling_match_id.EdgeId >= 0; }
                public bool loop_matches_a0() { return a0_loop_match_id.EdgeId >= 0; }

                // These fields count the number of edge crossings at a0, a1, and the
                // interior of (a0, a1).
                public int a0_crossings = 0;        // Count of polygon crossings at a0.
                public int a1_crossings = 0;        // Count of polygon crossings at a1.
                public int interior_crossings = 0;  // Count of polygon crossings in edge interior.

                public EdgeCrossingResult(bool matches_polyline, bool a0MatchesPolyline, bool a1MatchesPolyline, bool a0_matches_polygon, bool a1_matches_polygon, ShapeEdgeId polygon_match_id, ShapeEdgeId sibling_match_id, ShapeEdgeId a0_loop_match_id, int a0_crossings, int a1_crossings, int interior_crossings)
                {
                    this.matches_polyline = matches_polyline;
                    A0MatchesPolyline = a0MatchesPolyline;
                    A1MatchesPolyline = a1MatchesPolyline;
                    this.a0_matches_polygon = a0_matches_polygon;
                    this.a1_matches_polygon = a1_matches_polygon;
                    this.polygon_match_id = polygon_match_id;
                    this.sibling_match_id = sibling_match_id;
                    this.a0_loop_match_id = a0_loop_match_id;
                    this.a0_crossings = a0_crossings;
                    this.a1_crossings = a1_crossings;
                    this.interior_crossings = interior_crossings;
                }
            }

            private InputEdgeId InputEdgeId() => input_dimensions_.Count;

            // Returns true if the edges on either side of the first vertex of the
            // current edge have not been emitted.
            //
            // REQUIRES: This method is called just after updating "inside_" for "v0".
            private bool IsV0Isolated(Edge a_id)
            {
                return !Inside && v0_emitted_max_edge_id_ < a_id.EdgeId;
            }

            // Returns true if "a_id" is the last edge of the current chain, and the
            // edges on either side of the last vertex have not been emitted (including
            // the possibility that the chain forms a loop).
            private bool IsChainLastVertexIsolated(Edge a_id)
            {
                return (a_id.EdgeId == chain_limit_ - 1 && !chain_v0_emitted_ &&
                        v0_emitted_max_edge_id_ <= a_id.EdgeId);
            }

            // Returns true if the given polyline edge contains "v0", taking into
            // account the specified PolylineModel.
            private bool PolylineContainsV0(int edge_id, int chain_start)
            {
                return (polyline_model_ != PolylineModel.OPEN || edge_id > chain_start);
            }

            private void AddCrossing((SourceId, bool) crossing)
            {
                if (!tracker_.AddSpace(source_edge_crossings_, 1)) return;
                source_edge_crossings_.Add((InputEdgeId(), crossing));
            }

            private void AddInteriorCrossing(SourceEdgeCrossing crossing)
            {
                // Crossing edges are queued until the S2Builder edge that they are
                // supposed to be associated with is created (see AddEdge() and
                // pending_source_edge_crossings_ for details).
                pending_source_edge_crossings_.Add(crossing);
            }

            private void SetClippingState(InputEdgeId parameter, bool state) =>
                AddCrossing((new SourceId(parameter), state));

            // Supports "early exit" in the case of boolean results by returning false
            // as soon as the result is known to be non-empty.
            private bool AddEdge(Edge a_id, S2Shape.Edge a, int dimension, int interior_crossings)
            {
                if (builder_ == null) return false;  // Boolean output.
                if (interior_crossings > 0)
                {
                    // Add the edges that cross this edge to the output so that
                    // GraphEdgeClipper can find them.
                    if (!tracker_.AddSpace(source_edge_crossings_,
                                            pending_source_edge_crossings_.Count))
                    {
                        return false;
                    }
                    foreach (var crossing in pending_source_edge_crossings_)
                    {
                        source_edge_crossings_.Add(new(InputEdgeId(), crossing));
                    }

                    // Build a map that translates temporary edge ids (SourceId) to
                    // the representation used by EdgeClippingLayer (InputEdgeId).
                    if (!tracker_.TallySourceIdMap<SourceIdMap>(1))
                    {
                        return false;
                    }
                    var src_id = new SourceId(a_region_id_, a_id.ShapeId, a_id.EdgeId);
                    source_id_map_[src_id] = InputEdgeId();
                }
                // Set the GraphEdgeClipper's "inside" state to match ours.
                if (Inside != prev_inside_) SetClippingState(kSetInside, Inside);
                if (!tracker_.AddSpace(input_dimensions_, 1)) return false;
                input_dimensions_.Add((sbyte)dimension);
                builder_.AddEdge(a.V0, a.V1);
                Inside ^= (interior_crossings & 1) != 0;
                prev_inside_ = Inside;
                return tracker_.Ok();
            }

            // Supports "early exit" in the case of boolean results by returning false
            // as soon as the result is known to be non-empty.
            private bool AddPointEdge(S2Point p, int dimension)
            {
                if (builder_ == null) return false;  // Boolean output.
                if (!prev_inside_) SetClippingState(kSetInside, true);
                if (!tracker_.AddSpace(input_dimensions_, 1)) return false;
                input_dimensions_.Add((sbyte)dimension);
                builder_.AddEdge(p, p);
                prev_inside_ = true;
                return tracker_.Ok();
            }

            // Processes an edge of dimension 0 (i.e., a point) from region A.
            //
            // Supports "early exit" in the case of boolean results by returning false
            // as soon as the result is known to be non-empty.
            private bool ProcessEdge0(Edge a_id, S2Shape.Edge a, CrossingIterator it)
            {
                System.Diagnostics.Debug.Assert(a.V0 == a.V1);
                // When a region is inverted, all points and polylines are discarded.
                if (invert_a_ != invert_result_)
                {
                    SkipCrossings(a_id, it);
                    return true;
                }
                var r = ProcessPointCrossings(a_id, a.V0, it);

                // "contained" indicates whether the current point is inside the polygonal
                // interior of the opposite region, using semi-open boundaries.
                bool contained = Inside ^ invert_b_;
                if (r.MatchesPolygon && polygon_model_ != PolygonModel.SEMI_OPEN)
                {
                    contained = (polygon_model_ == PolygonModel.CLOSED);
                }
                if (r.MatchesPolyline) contained = true;

                // The output of UNION includes duplicate values, so ensure that points are
                // not suppressed by other points.
                if (r.MatchesPoint && !is_union_) contained = true;

                // Test whether the point is contained after region B is inverted.
                if (contained == invert_b_) return true;  // Don't exit early.
                return AddPointEdge(a.V0, 0);
            }
            // Processes an edge of dimension 1 (i.e., a polyline edge) from region A.
            //
            // Supports "early exit" in the case of boolean results by returning false
            // as soon as the result is known to be non-empty.
            private bool ProcessEdge1(Edge a_id, S2Shape.Edge a, CrossingIterator it)
            {
                // When a region is inverted, all points and polylines are discarded.
                if (invert_a_ != invert_result_)
                {
                    SkipCrossings(a_id, it);
                    return true;
                }
                // Evaluate whether the start vertex should belong to the output, in case it
                // needs to be emitted as an isolated vertex.
                var r = ProcessEdgeCrossings(a_id, a, it);
                bool a0_inside = IsPolylineVertexInside(r.A0MatchesPolyline, r.a0_matches_polygon);

                // Test whether the entire polyline edge should be emitted (or not emitted)
                // because it matches a polyline or polygon edge.
                bool is_degenerate = (a.V0 == a.V1);
                Inside ^= (r.a0_crossings & 1) != 0;
                if (Inside != IsPolylineEdgeInside(r, is_degenerate))
                {
                    Inside ^= true;   // Invert the inside_ state.
                    ++r.a1_crossings;  // Restore the correct (semi-open) state later.
                }

                // If neither edge adjacent to v0 was emitted, and this polyline contains
                // v0, and the other region contains v0, then emit an isolated vertex.
                if (!polyline_loops_have_boundaries_ && a_id.EdgeId == chain_start_ &&
                    a.V0 == a_shape_.ChainEdge(chain_id_,
                        chain_limit_ - chain_start_ - 1).V1)
                {
                    // This is the first vertex of a polyline loop, so we can't decide if it
                    // is isolated until we process the last polyline edge.
                    chain_v0_emitted_ = Inside;
                }
                else if (IsV0Isolated(a_id) && !is_degenerate &&
                         PolylineContainsV0(a_id.EdgeId, chain_start_) && a0_inside)
                {
                    if (!AddPointEdge(a.V0, 1)) return false;
                }

                // Test whether the entire edge or any part of it belongs to the output.
                if (Inside || r.interior_crossings > 0)
                {
                    // Note: updates "inside_" to correspond to the state just before a1.
                    if (!AddEdge(a_id, a, 1 /*dimension*/, r.interior_crossings))
                    {
                        return false;
                    }
                }
                // Remember whether the edge portion just before "a1" was emitted, so that
                // we can decide whether "a1" need to be emitted as an isolated vertex.
                if (Inside) v0_emitted_max_edge_id_ = a_id.EdgeId + 1;

                // Verify that edge crossings are being counted correctly.
                Inside ^= (r.a1_crossings & 1) != 0;
                if (it.CrossingsComplete())
                {
                    System.Diagnostics.Debug.Assert(it.BIndex().MakeS2ContainsPointQuery().Contains(a.V1) == Inside ^ invert_b_);
                }

                // Special case to test whether the last vertex of a polyline should be
                // emitted as an isolated vertex.
                if (it.CrossingsComplete() && !is_degenerate &&
                    IsChainLastVertexIsolated(a_id) &&
                    (polyline_model_ == PolylineModel.CLOSED ||
                     (!polyline_loops_have_boundaries_ &&
                      a.V1 == a_shape_.ChainEdge(chain_id_, chain_start_).V0)) &&
                    IsPolylineVertexInside(r.A1MatchesPolyline, r.a1_matches_polygon))
                {
                    if (!AddPointEdge(a.V1, 1)) return false;
                }
                return true;
            }
            // Processes an edge of dimension 2 (i.e., a polygon edge) from region A.
            //
            // Supports "early exit" in the case of boolean results by returning false
            // as soon as the result is known to be non-empty.
            private bool ProcessEdge2(Edge a_id, S2Shape.Edge a, CrossingIterator it)
            {
                // Whenever the two regions contain the same edge, or opposite edges of a
                // sibling pair, or one region contains a point loop while the other
                // contains a matching vertex, then in general the result depends on whether
                // one or both sides represent a degenerate shell or hole.
                //
                // In each pass it is easy to determine whether edges in region B represent
                // degenerate geometry, and if so whether they represent a shell or hole,
                // since this can be determined from the inside_ state and the
                // matches_polygon() / matches_sibling() methods of EdgeCrossingResult.
                // However this information is not readily available for region A.
                //
                // We handle this by saving the shell/hole status of each degenerate loop in
                // region B during the first pass, and deferring the processing of any edges
                // that meet the criteria above until the second pass.  (Note that regions
                // A,B correspond to regions 0,1 respectively in the first pass whereas they
                // refer to regions 1,0 respectively in the second pass.)
                //
                // The first pass ignores:
                //  - degenerate edges of A that are incident to any edge of B
                //  - non-degenerate edges of A that match or are siblings to an edge of B
                //
                // The first pass also records the shell/hole status of:
                //  - degenerate edges of B that are incident to any edge of A
                //  - sibling pairs of B where either edge matches an edge of A
                //
                // The second pass processes and perhaps outputs:
                //  - degenerate edges of B that are incident to any edge of A
                //  - non-degenerate edges of B that match or are siblings to an edge of A
                //
                // The following flag indicates that we are in the second pass described
                // above, i.e. that we are emitting any necessary edges that were ignored by
                // the first pass.
                bool emit_shared = (a_region_id_ == 1);

                // Degeneracies such as isolated vertices and sibling pairs can only be
                // created by intersecting CLOSED polygons or unioning OPEN polygons.
                bool create_degen =
                    (polygon_model_ == PolygonModel.CLOSED && !invert_a_ && !invert_b_) ||
                    (polygon_model_ == PolygonModel.OPEN && invert_a_ && invert_b_);

                // In addition, existing degeneracies are kept when an open boundary is
                // subtracted.  Note that "keep_degen_b" is only defined for completeness.
                // It is needed to ensure that the "reverse subtraction operator" (B - A)
                // preserves degeneracies correctly, however in practice this operator is
                // only used internally to implement symmetric difference, and in that
                // situation the preserved degeneracy is always removed from the final
                // result because it overlaps other geometry.
                bool keep_degen_a = (polygon_model_ == PolygonModel.OPEN && invert_b_);
                bool keep_degen_b = (polygon_model_ == PolygonModel.OPEN && invert_a_);

                var r = ProcessEdgeCrossings(a_id, a, it);
                System.Diagnostics.Debug.Assert(!r.matches_polyline);

                // If only one region is inverted, matching/sibling relations are reversed.
                if (invert_a_ != invert_b_)
                {
                    (r.sibling_match_id, r.polygon_match_id) = (r.polygon_match_id, r.sibling_match_id);
                }

                bool is_point = (a.V0 == a.V1);
                if (!emit_shared)
                {
                    // Remember the shell/hole status of degenerate B edges that are incident
                    // to any edge of A.  (We don't need to do this for vertex a1 since it is
                    // the same as vertex a0 of the following A loop edge.)
                    if (r.loop_matches_a0())
                    {
                        is_degenerate_hole_[r.a0_loop_match_id] = Inside;
                        if (is_point) return true;
                    }

                    // Point loops are handled identically to points in the semi-open model,
                    // and are easier to process in the first pass (since otherwise in the
                    // r.a0_matches_polygon case we would need to remember the containment
                    // status of the matching vertex).  Otherwise we defer processing such
                    // loops to the second pass so that we can distinguish whether the
                    // degenerate edge represents a hole or shell.
                    if (polygon_model_ != PolygonModel.SEMI_OPEN)
                    {
                        if (is_point && r.a0_matches_polygon) return true;
                    }
                }
                Inside ^= (r.a0_crossings & 1) == 1;
                if (!emit_shared)
                {
                    // Defer processing A edges that match or are siblings to an edge of B.
                    if (r.matches_polygon() || r.matches_sibling())
                    {
                        // For sibling pairs, also remember their shell/hole status.
                        if (r.matches_polygon() && r.matches_sibling())
                        {
                            is_degenerate_hole_[r.polygon_match_id] = Inside;
                            is_degenerate_hole_[r.sibling_match_id] = Inside;
                        }
                        System.Diagnostics.Debug.Assert(r.interior_crossings == 0);
                        Inside ^= (r.a1_crossings & 1) == 1;
                        return true;
                    }
                }

                // Remember whether the B geometry represents a sibling pair hole.
                bool is_b_hole = r.matches_polygon() && r.matches_sibling() && Inside;

                // At this point, "inside_" indicates whether the initial part of the A edge
                // is contained by the B geometry using semi-open rules.  The following code
                // implements the various other polygon boundary rules by changing the value
                // of "inside_" when necessary to indicate whether the current A edge should
                // be emitted to the output or not.  "semi_open_inside" remembers the true
                // value of "inside_" so that it can be restored later.
                bool semi_open_inside = Inside;
                if (is_point)
                {
                    if (r.loop_matches_a0())
                    {
                        // Both sides are point loops.  The edge is kept only:
                        //  - for closed intersection, open union, and open difference;
                        //  - if A and B are both holes or both shells.
                        Inside = create_degen || keep_degen_a ||
                                  (Inside == is_degenerate_hole_[r.a0_loop_match_id]);
                    }
                    else if (r.a0_matches_polygon)
                    {
                        // A point loop in A matches a polygon vertex in B.  Note that this code
                        // can emit an extra isolated vertex if A represents a point hole, but
                        // this doesn't matter (see comments on the call to AddPointEdge below).
                        if (polygon_model_ != PolygonModel.SEMI_OPEN)
                        {
                            Inside = create_degen || keep_degen_a;
                        }
                    }
                }
                else if (r.matches_polygon())
                {
                    if (is_degenerate(a_id))
                    {
                        // The A edge has a sibling.  The edge is kept only:
                        //  - for closed intersection, open union, and open difference;
                        //  - if the A sibling pair is a hole and the B edge has no sibling; or
                        //  - if the B geometry is also a sibling pair and A and B are both
                        //    holes or both shells.
                        Inside = create_degen || keep_degen_a ||
                                  (!r.matches_sibling() || Inside) == is_degenerate_hole_[a_id];
                    }
                    else
                    {
                        // Matching edges are kept unless the B geometry is a sibling pair, in
                        // which case it is kept only for closed intersection, open union, and
                        // open difference.
                        if (!r.matches_sibling() || create_degen || keep_degen_b) Inside = true;
                    }
                }
                else if (r.matches_sibling())
                {
                    if (is_degenerate(a_id))
                    {
                        // The A edge has a sibling.  The edge is kept only if A is a sibling
                        // pair shell and the operation is closed intersection, open union, or
                        // open difference.
                        Inside = (create_degen || keep_degen_a) && !is_degenerate_hole_[a_id];
                    }
                    else
                    {
                        Inside = create_degen;
                    }
                }
                if (Inside != semi_open_inside)
                {
                    ++r.a1_crossings;  // Restores the correct (semi-open) state later.
                }

                // Test whether the first vertex of this edge should be emitted as an
                // isolated degenerate vertex.  This is only needed in the second pass when:
                //  - a0 matches a vertex of the B polygon;
                //  - the initial part of the A edge will not be emitted; and
                //  - the operation is closed intersection or open union, or open difference
                //    and the B geometry is a point loop.
                //
                // The logic does not attempt to avoid redundant extra vertices (e.g. the
                // extra code in ProcessEdge1() that checks whether the vertex is the
                // endpoint of the preceding emitted edge) since these these will be removed
                // during S2Builder.Graph creation by DegenerateEdges.DISCARD or
                // DISCARD_EXCESS (which are necessary in any case due to snapping).
                if (emit_shared && r.a0_matches_polygon && !Inside &&
                    (create_degen || (keep_degen_b && r.loop_matches_a0())))
                {
                    if (!AddPointEdge(a.V0, 2)) return false;
                }

                // Since we skipped edges in the first pass that only had a sibling pair
                // match in the B geometry, we sometimes need to emit the sibling pair of an
                // edge in the second pass.  This happens only if:
                //  - the operation is closed intersection, open union, or open difference;
                //  - the A geometry is not a sibling pair (since otherwise we will process
                //    that edge as well); and
                //  - the B geometry is not a sibling pair hole (since then only one edge
                //    should be emitted).
                if (r.matches_sibling() && (create_degen || keep_degen_b) &&
                    !is_degenerate(a_id) && !is_b_hole)
                {
                    S2Shape.Edge sibling = new(a.V1, a.V0);
                    if (!AddEdge(r.sibling_match_id, sibling, 2 /*dimension*/, 0))
                    {
                        return false;
                    }
                }

                // Test whether the entire edge or any part of it belongs to the output.
                if (Inside || r.interior_crossings > 0)
                {
                    // Note: updates "inside_" to correspond to the state just before a1.
                    if (!AddEdge(a_id, a, 2 /*dimension*/, r.interior_crossings))
                    {
                        return false;
                    }
                }

                Inside ^= (r.a1_crossings & 1) != 0;

                // Verify that edge crossings are being counted correctly.
                if (it.CrossingsComplete())
                {
                    System.Diagnostics.Debug.Assert(it.BIndex().MakeS2ContainsPointQuery().Contains(a.V1) == Inside ^ invert_b_);
                }

                return true;
            }

            // Skip any crossings that were not needed to determine the result.
            private static void SkipCrossings(Edge a_id, CrossingIterator it)
            {
                while (!it.Done(a_id)) it.Next();
            }
            // Returns a summary of the relationship between a point from region A and
            // a set of crossing edges from region B (see PointCrossingResult).
            private PointCrossingResult ProcessPointCrossings(Edge a_id, S2Point a0, CrossingIterator it)
            {
                var r = new PointCrossingResult();
                for (; !it.Done(a_id); it.Next())
                {
                    if (it.BDimension() == 0)
                    {
                        r.MatchesPoint = true;
                    }
                    else if (it.BDimension() == 1)
                    {
                        if (PolylineEdgeContainsVertex(a0, it, 0))
                        {
                            r.MatchesPolyline = true;
                        }
                    }
                    else
                    {
                        r.MatchesPolygon = true;
                    }
                }
                return r;
            }
            // Returns a summary of the relationship between a test edge from region A and
            // a set of crossing edges from region B (see EdgeCrossingResult).
            //
            // NOTE(ericv): We could save a bit of work when matching polygon vertices by
            // passing in a flag saying whether this information is needed.  For example
            // it is only needed in ProcessEdge2 when (emit_shared && create_degenerate).
            private EdgeCrossingResult ProcessEdgeCrossings(Edge a_id, S2Shape.Edge a, CrossingIterator it)
            {
                pending_source_edge_crossings_.Clear();
                var r = new EdgeCrossingResult();
                if (it.Done(a_id)) return r;

                for (; !it.Done(a_id); it.Next())
                {
                    // Polyline and polygon "inside" states are not affected by point geometry.
                    if (it.BDimension() == 0) continue;
                    var b = it.BEdge();
                    if (it.IsInteriorCrossing())
                    {
                        // The crossing occurs in the edge interior.  The condition below says
                        // that (1) polyline crossings don't affect the polygon "inside" state,
                        // and (2) subtracting a crossing polyline from a polyline does not
                        // affect its "inside" state.  (Note that vertices are still created at
                        // the intersection points.)
                        if (a_dimension_ <= it.BDimension() &&
                            !(invert_b_ != invert_result_ && it.BDimension() == 1))
                        {
                            var src_id = new SourceId(b_region_id_, it.BShapeId(), it.BEdgeId());
                            AddInteriorCrossing((src_id, it.LeftToRight()));
                        }
                        r.interior_crossings += (it.BDimension() == 1) ? 2 : 1;
                    }
                    else if (it.BDimension() == 1)
                    {
                        // The polygon "inside" state is not affected by polyline geometry.
                        if (a_dimension_ == 2) continue;
                        if ((a.V0 == b.V0 && a.V1 == b.V1) || (a.V0 == b.V1 && a.V1 == b.V0))
                        {
                            r.matches_polyline = true;
                        }
                        if ((a.V0 == b.V0 || a.V0 == b.V1) &&
                            PolylineEdgeContainsVertex(a.V0, it, 1))
                        {
                            r.A0MatchesPolyline = true;
                        }
                        if ((a.V1 == b.V0 || a.V1 == b.V1) &&
                            PolylineEdgeContainsVertex(a.V1, it, 1))
                        {
                            r.A1MatchesPolyline = true;
                        }
                    }
                    else
                    {
                        System.Diagnostics.Debug.Assert(2 == it.BDimension());
                        if (a.V0 == a.V1 || b.V0 == b.V1)
                        {
                            // There are no edge crossings since at least one edge is degenerate.
                            if (a.V0 == b.V0 && a.V0 == b.V1)
                            {
                                r.a0_loop_match_id = it.BId();
                            }
                        }
                        else if (a.V0 == b.V0 && a.V1 == b.V1)
                        {
                            ++r.a0_crossings;
                            r.polygon_match_id = it.BId();
                        }
                        else if (a.V0 == b.V1 && a.V1 == b.V0)
                        {
                            ++r.a0_crossings;
                            r.sibling_match_id = it.BId();
                        }
                        else if (it.IsVertexCrossing())
                        {
                            if (a.V0 == b.V0 || a.V0 == b.V1)
                            {
                                ++r.a0_crossings;
                            }
                            else
                            {
                                ++r.a1_crossings;
                            }
                        }
                        if (a.V0 == b.V0 || a.V0 == b.V1)
                        {
                            r.a0_matches_polygon = true;
                        }
                        if (a.V1 == b.V0 || a.V1 == b.V1)
                        {
                            r.a1_matches_polygon = true;
                        }
                    }
                }
                return r;
            }

            // Returns true if the current point being processed (which must be a polyline
            // vertex) is contained by the opposite region (after inversion if "invert_b_"
            // is true).  "matches_polyline" and "matches_polygon" indicate whether the
            // vertex matches a polyline/polygon vertex of the opposite region.
            private bool IsPolylineVertexInside(bool matches_polyline, bool matches_polygon)
            {
                // Initially "contained" indicates whether the current point is inside the
                // polygonal interior of region B using semi-open boundaries.
                bool contained = Inside ^ invert_b_;

                // For UNION the output includes duplicate polylines.  The test below
                // ensures that isolated polyline vertices are not suppressed by other
                // polyline vertices in the output.
                if (matches_polyline && !is_union_)
                {
                    contained = true;
                }
                else if (matches_polygon && polygon_model_ != PolygonModel.SEMI_OPEN)
                {
                    contained = (polygon_model_ == PolygonModel.CLOSED);
                }
                // Finally, invert the result if the opposite region should be inverted.
                return contained ^ invert_b_;
            }
            // Returns true if the current polyline edge is contained by the opposite
            // region (after inversion if "invert_b_" is true).
            private bool IsPolylineEdgeInside(EdgeCrossingResult r, bool is_degenerate)
            {
                bool contained = Inside ^ invert_b_;

                // Note that if r.matches_polyline and is_union_ is true, then "contained"
                // will be false (unless there is also a matching polygon edge) since
                // polyline edges are not allowed in the interior of B.  In this case we
                // leave "contained" as false since it causes both matching edges to be
                // emitted.
                if (r.matches_polyline && !is_union_)
                {
                    contained = true;
                }
                else if (is_degenerate)
                {
                    // First allow the polygon boundary model to override the semi-open rules.
                    // Note that a polygon vertex (dimension 2) is considered to completely
                    // contain degenerate OPEN and SEMI_OPEN polylines (dimension 1) even
                    // though the latter do not contain any points.  This is because dimension
                    // 2 points are considered to be a strict superset of dimension 1 points.
                    if (polygon_model_ != PolygonModel.SEMI_OPEN && r.a0_matches_polygon)
                    {
                        contained = (polygon_model_ == PolygonModel.CLOSED);
                    }
                    // Note that r.a0_matches_polyline is true if and only if some B polyline
                    // contains the degenerate polyline (a0, a0).
                    if (r.A0MatchesPolyline && !is_union_) contained = true;
                }
                else if (r.matches_polygon())
                {
                    // In the SEMI_OPEN model, polygon sibling pairs cancel each other and
                    // have no effect on point or edge containment.
                    if (!(polygon_model_ == PolygonModel.SEMI_OPEN && r.matches_sibling()))
                    {
                        contained = (polygon_model_ != PolygonModel.OPEN);
                    }
                }
                else if (r.matches_sibling())
                {
                    contained = (polygon_model_ == PolygonModel.CLOSED);
                }
                // Finally, invert the result if the opposite region should be inverted.
                return contained ^ invert_b_;
            }
            // Returns true if the vertex "v" is contained by the polyline edge referred
            // to by the CrossingIterator "it", taking into account the PolylineModel.
            // "dimension" is 0 or 1 according to whether "v" should be modeled as a point
            // or as a degenerate polyline.  (This only makes a difference when the
            // containing polyline is degenerate, since the polyline AA contains itself in
            // all boundary models but contains the point A only in the CLOSED model.)
            //
            // REQUIRES: it.b_dimension() == 1
            // REQUIRES: "v" is an endpoint of it.b_edge()
            private bool PolylineEdgeContainsVertex(S2Point v, CrossingIterator it, int dimension)
            {
                System.Diagnostics.Debug.Assert(1 == it.BDimension());
                System.Diagnostics.Debug.Assert(it.BEdge().V0 == v || it.BEdge().V1 == v);
                System.Diagnostics.Debug.Assert(dimension == 0 || dimension == 1);

                // Closed polylines contain all their vertices.
                if (polyline_model_ == PolylineModel.CLOSED) return true;

                // Note that the code below is structured so that it.b_edge() is not usually
                // needed (since accessing the edge can be relatively expensive).
                var (ChainId, Start, Limit) = it.BChainInfo();
                int b_edge_id = it.BEdgeId();

                // A polyline contains its last vertex only when the polyline is degenerate
                // (v0 == v1) and "v" is modeled as a degenerate polyline (dimension == 1).
                // This corresponds to the fact that the polyline AA contains itself in all
                // boundary models, but contains the point A only in the CLOSED model.
                if (b_edge_id == Limit - 1 && v == it.BEdge().V1 &&
                    (dimension == 0 || b_edge_id > 0 || v != it.BEdge().V0))
                {
                    return false;
                }

                // Otherwise all interior vertices are contained.  The first polyline
                // vertex is contained if either the polyline model is not OPEN, or the
                // polyline forms a loop and polyline_loops_have_boundaries_ is false.
                if (PolylineContainsV0(b_edge_id, Start)) return true;
                if (v != it.BEdge().V0) return true;
                if (polyline_loops_have_boundaries_) return false;
                return v == it.BShape().ChainEdge(ChainId, Limit - Start - 1).V1;
            }

            // Constructor parameters:

            private readonly PolygonModel polygon_model_;
            private readonly PolylineModel polyline_model_;
            private readonly bool polyline_loops_have_boundaries_;

            // The output of the CrossingProcessor consists of a subset of the input
            // edges that are emitted to "builder_", and some auxiliary information
            // that allows GraphEdgeClipper to determine which segments of those input
            // edges belong to the output.  The auxiliary information consists of the
            // dimension of each input edge, and set of input edges from the other
            // region that cross each input input edge.
            private readonly S2Builder? builder_;  // (null if boolean output was requested)
            private readonly List<sbyte> input_dimensions_;
            private readonly InputEdgeCrossings input_crossings_;
            private readonly MemoryTracker tracker_;

            // Fields set by StartBoundary:

            private int a_region_id_, b_region_id_;
            private bool invert_a_, invert_b_, invert_result_;
            private bool is_union_;  // True if this is a UNION operation.

            // Fields set by StartShape:

            private S2Shape a_shape_;
            private int a_dimension_;

            // Fields set by StartChain:

            private int chain_id_;
            private int chain_start_;
            private int chain_limit_;

            // Fields updated by ProcessEdge:

            // A temporary representation of input_crossings_ that is used internally
            // until all necessary edges from *both* polygons have been emitted to the
            // S2Builder.  This field is then converted by DoneBoundaryPair() into
            // the InputEdgeCrossings format expected by GraphEdgeClipper.
            //
            // The reason that we can'truct input_crossings_ directly is that it
            // uses Int32s to identify the edges from both polygons, and when we
            // are processing edges from the first polygon, Int32s have not yet
            // been assigned to the second polygon.  So instead this field identifies
            // edges from the first polygon using an InputEdgeId, and edges from the
            // second polygon using a (region_id, shape_id, edge_id) tuple (i.e., a
            // SourceId).
            //
            // All crossings are represented twice, once to indicate that an edge from
            // polygon 0 is crossed by an edge from polygon 1, and once to indicate that
            // an edge from polygon 1 is crossed by an edge from polygon 0.  The entries
            // are sorted lexicographically by their eventual InputEdgeIds except for
            // GraphEdgeClipper state modifications, which are sorted by the first
            // InputEdgeId only.
            private readonly SourceEdgeCrossings source_edge_crossings_ = new();

            // A set of edges that cross the current edge being processed by
            // ProcessEdge() but that have not yet been associated with a particular
            // S2Builder edge.  This is necessary because ProcessEdge can create up to
            // three S2Builder edges per input edge: one to represent the edge interior,
            // and up to two more to represent an isolated start and/or end vertex.  The
            // crossing edges must be associated with the S2Builder edge that represents
            // the edge interior, and they are stored here until that edge is created.
            private readonly List<SourceEdgeCrossing> pending_source_edge_crossings_ = new();

            // A map that translates from SourceId (the (region_id, shape_id, edge_id)
            // triple that identifies an S2ShapeIndex edge) to InputEdgeId (the
            // sequentially increasing numbers assigned to input edges by S2Builder).
            private readonly SourceIdMap source_id_map_ = new();

            // For each edge in region B that defines a degenerate loop (either a point
            // loop or a sibling pair), indicates whether that loop represents a shell
            // or a hole.  This information is used during the second pass of
            // AddBoundaryPair() to determine the output for degenerate edges.
            private readonly Dictionary<ShapeEdgeId, bool> is_degenerate_hole_ = new();

            // The value of that "inside_" would have just before the end of the
            // previous edge added to S2Builder.  This value is used to determine
            // whether the GraphEdgeClipper state needs to be updated when jumping from
            // one edge chain to another.
            private bool prev_inside_;

            // The maximum edge id of any edge in the current chain whose v0 vertex has
            // already been emitted.  This is used to determine when an isolated vertex
            // needs to be emitted, e.g. when two closed polygons share only a vertex.
            private int v0_emitted_max_edge_id_;

            // True if the first vertex of the current chain has been emitted.  This is
            // used when processing loops in order to determine whether the first/last
            // vertex of the loop should be emitted as an isolated vertex.
            private bool chain_v0_emitted_;
        }

        // A helper class for iterating through the edges from region B that cross a
        // particular edge from region A.  It caches information from the current
        // shape, chain, and edge so that it doesn't need to be looked up repeatedly.
        // Typical usage:
        //
        //  void SomeFunction(S2ShapeUtil.Edge a_id, CrossingIterator *it) {
        //    // Iterate through the edges that cross edge "a_id".
        //    for (; !it.Done(a_id); it.Next()) {
        //      ... use it.b_shape(), it.b_edge(), etc ...
        //    }
        private class CrossingIterator
        {
            // Creates an iterator over crossing edge pairs (a, b) where "b" is an edge
            // from "b_index".  "crossings_complete" indicates that "crossings" contains
            // all edge crossings between the two regions (rather than a subset).
            public CrossingIterator(S2ShapeIndex b_index, List<IndexCrossing> crossings, bool crossings_complete)
            {
                b_index_ = b_index;
                crossings_ = crossings;
                index_ = 0;
                b_shape_id_ = -1;
                crossings_complete_ = crossings_complete;
                Update();
            }
            public void Next()
            {
                ++index_;
                Update();
            }
            public bool Done(Edge id) { return AId() != id; }

            // True if all edge crossings are available (see above).
            public bool CrossingsComplete() { return crossings_complete_; }

            // True if this crossing occurs at a point interior to both edges.
            public bool IsInteriorCrossing() { return crossings_[index_].IsInteriorCrossing; }

            // Equal to S2EdgeCrossings.VertexCrossing(a_edge, b_edge), provided that a_edge and
            // b_edge have exactly one vertex in common and neither edge is degenerate.
            public bool IsVertexCrossing() { return crossings_[index_].IsVertexCrossing; }

            // True if a_edge crosses b_edge from left to right (for interior crossings).
            public bool LeftToRight() { return crossings_[index_].LeftToRight; }

            public Edge AId() { return crossings_[index_].A; }
            public Edge BId() { return crossings_[index_].B; }
            public S2ShapeIndex BIndex() { return b_index_; }
            public S2Shape BShape() { return b_shape_; }
            public int BDimension() { return b_dimension_; }
            public int BShapeId() { return b_shape_id_; }
            public int BEdgeId() { return BId().EdgeId; }

            public S2Shape.Edge BEdge()
            {
                return b_shape_.GetEdge(BEdgeId());  // Opportunity to cache this.
            }

            // Returns a description of the chain to which the current B edge belongs.
            public (int ChainId, int Start, int Limit) BChainInfo()
            {
                if (b_info_.ChainId < 0)
                {
                    var chainId = BShape().GetChainPosition(BEdgeId()).ChainId;
                    var chain = BShape().GetChain(b_info_.ChainId);
                    var start = chain.Start;
                    var limit = chain.Start + chain.Length;
                    return (chainId, start, limit);
                }
                return b_info_;
            }

            // Updates information about the B shape whenever it changes.
            private void Update()
            {
                if (AId() != kSentinel && BId().ShapeId != b_shape_id_)
                {
                    b_shape_id_ = BId().ShapeId;
                    b_shape_ = b_index_.Shape(b_shape_id_);
                    b_dimension_ = b_shape_.Dimension();
                    b_info_ = (-1, -1, -1);  // Computed on demand.
                }
            }

            private readonly S2ShapeIndex b_index_;
            private S2Shape b_shape_;
            private int b_shape_id_;
            private int b_dimension_;
            private (int ChainId, int Start, int Limit) b_info_;  // Computed on demand.
            private readonly bool crossings_complete_;
            private readonly List<IndexCrossing> crossings_;
            private int index_;
        }

        private bool IsBooleanOutput() => op_.result_empty_ != null;

        // All of the methods below support "early exit" in the case of boolean
        // results by returning "false" as soon as the result is known to be
        // non-empty.
        //
        // Clips the boundary of A to the interior of the opposite region B and adds
        // the resulting edges to the output.  Optionally, any combination of region
        // A, region B, and the result may be inverted, which allows operations such
        // as union and difference to be implemented.
        //
        // Note that when an input region is inverted with respect to the output
        // (e.g., invert_a != invert_result), all polygon edges are reversed and all
        // points and polylines are discarded, since the complement of such objects
        // cannot be represented.  (If you want to compute the complement of points
        // or polylines, you can use S2LaxPolygonShape to represent your geometry as
        // degenerate polygons instead.)
        //
        // This method must be called an even number of times (first to clip A to B
        // and then to clip B to A), calling DoneBoundaryPair() after each pair.
        //
        // Supports "early exit" in the case of boolean results by returning false
        // as soon as the result is known to be non-empty.
        private bool AddBoundary(int a_region_id, bool invert_a, bool invert_b, bool invert_result, List<Edge> a_chain_starts, CrossingProcessor cp)
        {
            var a_index = op_.regions_[a_region_id];
            var b_index = op_.regions_[1 - a_region_id];
            if (!GetIndexCrossings(a_region_id)) return false;
            cp.StartBoundary(a_region_id, invert_a, invert_b, invert_result);

            // Walk the boundary of region A and build a list of all edge crossings.
            // We also keep track of whether the current vertex is inside region B.
            var next_start = 0;
            CrossingIterator next_crossing = new(b_index, index_crossings_, true /*crossings_complete*/);
            var next_id = new[] { a_chain_starts[next_start], next_crossing.AId() }.Min();
            while (next_id != kSentinel)
            {
                int a_shape_id = next_id.ShapeId;
                var a_shape = a_index.Shape(a_shape_id);
                cp.StartShape(a_shape);
                while (next_id.ShapeId == a_shape_id)
                {
                    // TODO(ericv): Special handling of dimension 0?  Can omit most of this
                    // code, including the loop, since all chains are of length 1.
                    int edge_id = next_id.EdgeId;
                    var chain_position = a_shape.GetChainPosition(edge_id);
                    int chain_id = chain_position.ChainId;
                    var chain = a_shape.GetChain(chain_id);
                    bool start_inside = next_id == a_chain_starts[next_start];
                    if (start_inside) ++next_start;
                    cp.StartChain(chain_id, chain, start_inside);
                    int chain_limit = chain.Start + chain.Length;
                    while (edge_id < chain_limit)
                    {
                        Edge a_id = new(a_shape_id, edge_id);
                        System.Diagnostics.Debug.Assert(cp.Inside || next_crossing.AId() == a_id);
                        if (!cp.ProcessEdge(a_id, next_crossing))
                        {
                            return false;
                        }
                        if (cp.Inside)
                        {
                            ++edge_id;
                        }
                        else if (next_crossing.AId().ShapeId == a_shape_id && next_crossing.AId().EdgeId < chain_limit)
                        {
                            edge_id = next_crossing.AId().EdgeId;
                        }
                        else
                        {
                            break;
                        }
                    }
                    next_id = new[] { a_chain_starts[next_start], next_crossing.AId() }.Min();
                }
            }
            return true;
        }
        // Returns the first edge of each edge chain from "a_region_id" whose first
        // vertex is contained by opposite region's polygons (using the semi-open
        // boundary model).  Each input region and the result region are inverted as
        // specified (invert_a, invert_b, and invert_result) before testing for
        // containment.  The algorithm uses these "chain starts" in order to clip the
        // boundary of A to the interior of B in an output-senstive way.
        //
        // This method supports "early exit" in the case where a boolean predicate is
        // being evaluated and the algorithm discovers that the result region will be
        // non-empty.
        private bool GetChainStarts(int a_region_id, bool invert_a, bool invert_b, bool invert_result, CrossingProcessor cp, List<Edge> chain_starts)
        {
            var a_index = op_.regions_[a_region_id];
            var b_index = op_.regions_[1 - a_region_id];

            if (IsBooleanOutput())
            {
                // If boolean output is requested, then we use the CrossingProcessor to
                // determine whether the first edge of each chain will be emitted to the
                // output region.  This lets us terminate the operation early in many
                // cases.
                cp.StartBoundary(a_region_id, invert_a, invert_b, invert_result);
            }

            // If region B has no two-dimensional shapes and is not inverted, then by
            // definition no chain starts are contained.  However if boolean output is
            // requested then we check for containment anyway, since as a side effect we
            // may discover that the result region is non-empty and terminate the entire
            // operation early.
            bool b_has_interior = HasInterior(b_index);
            if (b_has_interior || invert_b || IsBooleanOutput())
            {
                var query = b_index.MakeS2ContainsPointQuery();
                int num_shape_ids = a_index.NumShapeIds();
                for (int shape_id = 0; shape_id < num_shape_ids; ++shape_id)
                {
                    var a_shape = a_index.Shape(shape_id);
                    if (a_shape == null) continue;

                    // If region A is being subtracted from region B, points and polylines
                    // in region A can be ignored since these shapes never contribute to the
                    // output (they can only remove edges from region B).
                    if (invert_a != invert_result && a_shape.Dimension() < 2) continue;

                    if (IsBooleanOutput()) cp.StartShape(a_shape);
                    int num_chains = a_shape.NumChains();
                    for (int chain_id = 0; chain_id < num_chains; ++chain_id)
                    {
                        var chain = a_shape.GetChain(chain_id);
                        if (chain.Length == 0) continue;
                        ShapeEdge a = new(shape_id, chain.Start, a_shape.ChainEdge(chain_id, 0));
                        bool inside = (b_has_interior && query.Contains(a.V0)) != invert_b;
                        if (inside)
                        {
                            if (!tracker_.AddSpace(chain_starts, 1)) return false;
                            chain_starts.Add(new(shape_id, chain.Start));
                        }
                        if (IsBooleanOutput())
                        {
                            cp.StartChain(chain_id, chain, inside);
                            if (!ProcessIncidentEdges(a, query, cp)) return false;
                        }
                    }
                }
            }
            if (!tracker_.AddSpace(chain_starts, 1)) return false;
            chain_starts.Add(kSentinel);
            return true;
        }
        private bool ProcessIncidentEdges(ShapeEdge a, S2ContainsPointQuery<S2ShapeIndex> query, CrossingProcessor cp)
        {
            tmp_crossings_.Clear();
            query.VisitIncidentEdges(a.V0, (ShapeEdge b) => AddIndexCrossing(a, b, false /*is_interior*/, tmp_crossings_));
            // Fast path for the common case where there are no incident edges.  We
            // return false (terminating early) if the first chain edge will be emitted.
            if (!tmp_crossings_.Any())
            {
                return !cp.Inside;
            }
            // Otherwise we invoke the full CrossingProcessor logic to determine whether
            // the first chain edge will be emitted.
            if (tmp_crossings_.Count > 1)
            {
                tmp_crossings_.Sort();
                // VisitIncidentEdges() should not generate any duplicate values.
                var anyDup = false;
                for (var i = 1; i < tmp_crossings_.Count && !anyDup; i++)
                    anyDup = tmp_crossings_[i - 1] == tmp_crossings_[i];
                System.Diagnostics.Debug.Assert(!anyDup);
            }
            tmp_crossings_.Add(new IndexCrossing(kSentinel, kSentinel));
            var next_crossing = new CrossingIterator(query.Index, tmp_crossings_, false /*crossings_complete*/);
            return cp.ProcessEdge(a.Id, next_crossing);
        }
        private static bool HasInterior(S2ShapeIndex index)
        {
            for (var s = index.NumShapeIds(); --s >= 0;)
            {
                var shape = index.Shape(s);
                if (shape != null && shape.Dimension() == 2) return true;
            }
            return false;
        }
        private bool AddIndexCrossing(ShapeEdge a, ShapeEdge b, bool is_interior, List<IndexCrossing> crossings)
        {
            if (!tracker_.AddSpace(crossings, 1)) return false;
            crossings.Add(new IndexCrossing(a.Id, b.Id));
            var crossing = crossings.Last();
            if (is_interior)
            {
                crossing.IsInteriorCrossing = true;
                if (S2Pred.Sign(a.V0, a.V1, b.V0) > 0)
                {
                    crossing.LeftToRight = true;
                }
                builder_.AddIntersection(
                    S2.GetIntersection(a.V0, a.V1, b.V0, b.V1, null));
            }
            else
            {
                // TODO(ericv): This field isn't used unless one shape is a polygon and
                // the other is a polyline or polygon, but we don't have the shape
                // dimension information readily available here.
                if (S2.VertexCrossing(a.V0, a.V1, b.V0, b.V1))
                {
                    crossing.IsVertexCrossing = true;
                }
            }
            return true;  // Continue visiting.
        }
        // Initialize index_crossings_ to the set of crossing edge pairs such that the
        // first element of each pair is an edge from "region_id".
        //
        // Supports "early exit" in the case of boolean results by returning false
        // as soon as the result is known to be non-empty.
        private bool GetIndexCrossings(int region_id)
        {
            if (region_id == index_crossings_first_region_id_) return true;
            if (index_crossings_first_region_id_ < 0)
            {
                System.Diagnostics.Debug.Assert(region_id == 0);  // For efficiency, not correctness.
                // TODO(ericv): This would be more efficient if VisitCrossingEdgePairs()
                // returned the sign (+1 or -1) of the interior crossing, i.e.
                // "int interior_crossing_sign" rather than "bool is_interior".
                if (!S2ShapeUtil.EdgePairs.VisitCrossingEdgePairs(op_.regions_[0], op_.regions_[1], CrossingType.ALL,
                    (ShapeEdge a, ShapeEdge b, bool is_interior) =>
                    {
                        // For all supported operations (union, intersection, and
                        // difference), if the input edges have an interior crossing
                        // then the output is guaranteed to have at least one edge.
                        if (is_interior && IsBooleanOutput()) return false;
                        return AddIndexCrossing(a, b, is_interior, index_crossings_);
                    }))
                {
                    return false;
                }
                if (index_crossings_.Count > 1)
                {
                    var ss = new SortedSet<IndexCrossing>(index_crossings_);
                    index_crossings_.Clear();
                    index_crossings_.AddRange(ss);
                }
                // Add a sentinel value to simplify the loop logic.
                tracker_.AddSpace(index_crossings_, 1);
                index_crossings_.Add(new IndexCrossing(kSentinel, kSentinel));
                index_crossings_first_region_id_ = 0;
            }
            if (region_id != index_crossings_first_region_id_)
            {
                for (var i = 0; i < index_crossings_.Count; i++)
                {
                    var crossing = index_crossings_[i];
                    (crossing.B, crossing.A) = (crossing.A, crossing.B);
                    // The following predicates get inverted when the edges are swapped.
                    crossing.LeftToRight ^= true;
                    crossing.IsVertexCrossing ^= true;
                }
                index_crossings_.Sort();
                index_crossings_first_region_id_ = region_id;
            }
            return tracker_.Ok();
        }
        // Supports "early exit" in the case of boolean results by returning false
        // as soon as the result is known to be non-empty.
        private bool AddBoundaryPair(bool invert_a, bool invert_b, bool invert_result, CrossingProcessor cp)
        {
            // Optimization: if the operation is DIFFERENCE or SYMMETRIC_DIFFERENCE,
            // it is worthwhile checking whether the two regions are identical (in which
            // case the output is empty).
            var type = op_.OpType_;
            if (type == OpType.DIFFERENCE || type == OpType.SYMMETRIC_DIFFERENCE)
            {
                if (AreRegionsIdentical()) return true;
            }
            else if (IsBooleanOutput())
            {
                // TODO(ericv): When boolean output is requested there are other quick
                // checks that could be done here, such as checking whether a full cell from
                // one S2ShapeIndex intersects a non-empty cell of the other S2ShapeIndex.
            }
            List<Edge> a_starts = new(), b_starts = new();
            try
            {
                if (!GetChainStarts(0, invert_a, invert_b, invert_result, cp, a_starts) ||
                    !GetChainStarts(1, invert_b, invert_a, invert_result, cp, b_starts) ||
                    !AddBoundary(0, invert_a, invert_b, invert_result, a_starts, cp) ||
                    !AddBoundary(1, invert_b, invert_a, invert_result, b_starts, cp))
                {
                    return false;
                }
                if (!IsBooleanOutput()) cp.DoneBoundaryPair();
                return tracker_.Ok();
            }
            finally
            {
                tracker_.Untally(a_starts);
                tracker_.Untally(b_starts);
            }
        }
        // When subtracting regions, we can save a lot of work by detecting the
        // relatively common case where the two regions are identical.
        private bool AreRegionsIdentical()
        {
            var a = op_.regions_[0];
            var b = op_.regions_[1];
            if (a == b) return true;

            // If the regions are not identical, we would like to detect that fact as
            // quickly as possible.  In particular we would like to avoid fully decoding
            // both shapes if they are represented as encoded shape types.
            //
            // First we test whether the two geometries have the same dimensions and
            // chain structure.  This can be done without decoding any S2Points.
            int num_shape_ids = a.NumShapeIds();
            if (num_shape_ids != b.NumShapeIds()) return false;
            for (int s = 0; s < num_shape_ids; ++s)
            {
                var a_shape = a.Shape(s);
                var b_shape = b.Shape(s);
                int dimension = a_shape.Dimension();
                if (dimension != b_shape.Dimension()) return false;
                int num_chains = a_shape.NumChains();
                if (num_chains != b_shape.NumChains()) return false;
                int num_edges = a_shape.NumEdges();
                if (num_edges != b_shape.NumEdges()) return false;
                if (dimension == 0)
                {
                    System.Diagnostics.Debug.Assert(num_edges == num_chains);  // All chains are of length 1.
                    continue;
                }
                for (int c = 0; c < num_chains; ++c)
                {
                    var a_chain = a_shape.GetChain(c);
                    var b_chain = b_shape.GetChain(c);
                    System.Diagnostics.Debug.Assert(a_chain.Start == b_chain.Start);
                    if (a_chain.Length != b_chain.Length) return false;
                }
            }
            // Next we test whether both geometries have the same vertex positions.
            for (int s = 0; s < num_shape_ids; ++s)
            {
                var a_shape = a.Shape(s);
                var b_shape = b.Shape(s);
                int num_chains = a_shape.NumChains();
                for (int c = 0; c < num_chains; ++c)
                {
                    var a_chain = a_shape.GetChain(c);
                    for (int i = 0; i < a_chain.Length; ++i)
                    {
                        var a_edge = a_shape.ChainEdge(c, i);
                        var b_edge = b_shape.ChainEdge(c, i);
                        if (a_edge.V0 != b_edge.V0) return false;
                        if (a_edge.V1 != b_edge.V1) return false;
                    }
                }
                // Note that we don't need to test whether both shapes have the same
                // GetReferencePoint(), because S2Shape requires that the full geometry of
                // the shape (including its interior) must be derivable from its chains
                // and edges.  This is why the "full loop" exists; see s2shape.h.
            }
            return true;
        }
        // Supports "early exit" in the case of boolean results by returning false
        // as soon as the result is known to be non-empty.
        private bool BuildOpType(OpType op_type)
        {
            // CrossingProcessor does the real work of emitting the output edges.
            var cp = new CrossingProcessor(op_.Options_.PolygonModel_,
                                 op_.Options_.PolylineModel_,
                                 op_.Options_.PolylineLoopsHaveBoundaries,
                                 builder_, input_dimensions_, input_crossings_, tracker_);
            return op_type switch
            {
                OpType.UNION => AddBoundaryPair(true, true, true, cp),    // A | B == ~(~A & ~B)
                OpType.INTERSECTION => AddBoundaryPair(false, false, false, cp), // A & B

                // A - B = A & ~B
                //
                // Note that degeneracies are implemented such that the symmetric
                // operation (-B + A) also produces correct results.  This can be tested
                // by swapping op_->regions[0, 1] and calling AddBoundaryPair(true,
                // false, false), which computes (~B & A).
                OpType.DIFFERENCE => AddBoundaryPair(false, true, false, cp),
                OpType.SYMMETRIC_DIFFERENCE => (AddBoundaryPair(false, true, false, cp) &&
                                                AddBoundaryPair(true, false, false, cp)), // Compute the union of (A - B) and (B - A).
                _ => throw new ApplicationException("Invalid S2BooleanOperation.OpType"),
            };
        }
        // Given a polygon edge graph containing only degenerate edges and sibling edge
        // pairs, the purpose of this function is to decide whether the polygon is empty
        // or full except for the degeneracies, i.e. whether the degeneracies represent
        // shells or holes.
        private bool IsFullPolygonResult(/*S2Builder.Graph g,*/ out S2Error error)
        {
            error = S2Error.OK;

            // If there are no edges of dimension 2, the result could be either the
            // empty polygon or the full polygon.  Note that this is harder to determine
            // than you might think due to snapping.  For example, the union of two
            // non-empty polygons can be empty, because both polygons consist of tiny
            // loops that are eliminated by snapping.  Similarly, even if two polygons
            // both contain a common point their intersection can still be empty.
            //
            // We distinguish empty from full results using two heuristics:
            //
            //  1. We compute a bit mask representing the subset of the six S2 cube faces
            //     intersected by each input geometry, and use this to determine if only
            //     one of the two results is possible.  (This test is very fast.)  Note
            //     that snapping will never cause the result to cover an entire extra
            //     cube face because the maximum allowed snap radius is too small.
            System.Diagnostics.Debug.Assert(SnapFunction.kMaxSnapRadius.GetDegrees() <= 70);
            //
            //  2. We compute the area of each input geometry, and use this to bound the
            //     minimum and maximum area of the result.  If only one of {0, 4*Pi} is
            //     possible then we are done.  If neither is possible then we choose the
            //     one that is closest to being possible (since snapping can change the
            //     result area).  Both results are possible only when computing the
            //     symmetric difference of two regions of area 2*Pi each, in which case we
            //     must resort to additional heuristics (see below).
            //
            // TODO(ericv): Implement a predicate that uses the results of edge snapping
            // directly, rather than computing areas.  This would not only be much faster
            // but would also allows all cases to be handled 100% robustly.
            var a = op_.regions_[0];
            var b = op_.regions_[1];
            switch (op_.OpType_)
            {
                case OpType.UNION: return IsFullPolygonUnion(a, b);

                case OpType.INTERSECTION: return IsFullPolygonIntersection(a, b);

                case OpType.DIFFERENCE: return IsFullPolygonDifference(a, b);

                case OpType.SYMMETRIC_DIFFERENCE: return IsFullPolygonSymmetricDifference(a, b);

                default:
                    {
                        error = new(S2ErrorCode.INVALID_ARGUMENT, "Invalid S2BooleanOperation.OpType");
                        throw new ApplicationException("Invalid S2BooleanOperation.OpType");
                    }
            }
        }
        private static bool IsFullPolygonUnion(S2ShapeIndex a, S2ShapeIndex b)
        {
            // See comments in IsFullPolygonResult().  The most common case is that
            // neither input polygon is empty but the result is empty due to snapping.

            // The result can be full only if the union of the two input geometries
            // intersects all six faces of the S2 cube.  This test is fast.
            if ((GetFaceMask(a) | GetFaceMask(b)) != kAllFacesMask) return false;

            // The union area satisfies:
            //
            //   Math.Max(A, B) <= Union(A, B) <= Math.Min(4*Pi, A + B)
            //
            // where A, B can refer to a polygon or its area.  We then choose the result
            // that assumes the smallest amount of error.
            double a_area = S2ShapeIndexMeasures.GetArea(a), b_area = S2ShapeIndexMeasures.GetArea(b);
            double min_area = Math.Max(a_area, b_area);
            double max_area = Math.Min(S2.M_4_PI, a_area + b_area);
            return min_area > S2.M_4_PI - max_area;
        }
        private static bool IsFullPolygonIntersection(S2ShapeIndex a, S2ShapeIndex b)
        {
            // See comments in IsFullPolygonResult().  By far the most common case is
            // that the result is empty.

            // The result can be full only if each of the two input geometries
            // intersects all six faces of the S2 cube.  This test is fast.
            if ((GetFaceMask(a) & GetFaceMask(b)) != kAllFacesMask) return false;

            // The intersection area satisfies:
            //
            //   Math.Max(0, A + B - 4*Pi) <= Intersection(A, B) <= Math.Min(A, B)
            //
            // where A, B can refer to a polygon or its area.  We then choose the result
            // that assumes the smallest amount of error.
            double a_area = S2ShapeIndexMeasures.GetArea(a), b_area = S2ShapeIndexMeasures.GetArea(b);
            double min_area = Math.Max(0.0, a_area + b_area - S2.M_4_PI);
            double max_area = Math.Min(a_area, b_area);
            return min_area > S2.M_4_PI - max_area;
        }
        private static bool IsFullPolygonDifference(S2ShapeIndex a, S2ShapeIndex b)
        {
            // See comments in IsFullPolygonResult().  By far the most common case is
            // that the result is empty.

            // The result can be full only if each cube face is intersected by the first
            // geometry.  (The second geometry is irrelevant, since for example it could
            // consist of a tiny loop on each S2 cube face.)  This test is fast.
            if (GetFaceMask(a) != kAllFacesMask) return false;

            // The difference area satisfies:
            //
            //   Math.Max(0, A - B) <= Difference(A, B) <= Math.Min(A, 4*Pi - B)
            //
            // where A, B can refer to a polygon or its area.  We then choose the result
            // that assumes the smallest amount of error.
            double a_area = S2ShapeIndexMeasures.GetArea(a), b_area = S2ShapeIndexMeasures.GetArea(b);
            double min_area = Math.Max(0.0, a_area - b_area);
            double max_area = Math.Min(a_area, S2.M_4_PI - b_area);
            return min_area > S2.M_4_PI - max_area;
        }
        private bool IsFullPolygonSymmetricDifference(S2ShapeIndex a, S2ShapeIndex b)
        {
            // See comments in IsFullPolygonResult().  By far the most common case is
            // that the result is empty.

            // The result can be full only if the union of the two input geometries
            // intersects all six faces of the S2 cube.  This test is fast.
            byte a_mask = GetFaceMask(a);
            byte b_mask = GetFaceMask(b);
            if ((a_mask | b_mask) != kAllFacesMask) return false;

            // The symmetric difference area satisfies:
            //
            //   |A - B| <= SymmetricDifference(A, B) <= 4*Pi - |4*Pi - (A + B)|
            //
            // where A, B can refer to a polygon or its area.
            double a_area = S2ShapeIndexMeasures.GetArea(a), b_area = S2ShapeIndexMeasures.GetArea(b);
            double min_area = Math.Abs(a_area - b_area);
            double max_area = S2.M_4_PI - Math.Abs(S2.M_4_PI - (a_area + b_area));

            // Now we choose the result that assumes the smallest amount of error
            // (min_area in the empty case, and (4*Pi - max_area) in the full case).
            // However in the case of symmetric difference these two errors may be equal,
            // meaning that the result is ambiguous.  This happens when both polygons have
            // area 2*Pi.  Furthermore, this can happen even when the areas are not
            // exactly 2*Pi due to snapping and area calculation errors.
            //
            // To determine whether the result is ambiguous, we compute a rough estimate
            // of the maximum expected area error (including errors due to snapping),
            // using the worst-case error bound for a hemisphere defined by 4 vertices.
            var edge_snap_radius = builder_options_.EdgeSnapRadius();
            double hemisphere_area_error = S2.M_2_PI * edge_snap_radius.Radians +
                                           40 * S2.DoubleEpsilon;  // GetCurvatureMaxError

            // The following sign is the difference between the error needed for an empty
            // result and the error needed for a full result.  It is negative if an
            // empty result is possible, positive if a full result is possible, and zero
            // if both results are possible.
            double error_sign = min_area - (S2.M_4_PI - max_area);
            if (Math.Abs(error_sign) <= hemisphere_area_error)
            {
                // Handling the ambiguous case correctly requires a more sophisticated
                // algorithm (see below), but we can at least handle the simple cases by
                // testing whether both input geometries intersect all 6 cube faces.  If
                // not, then the result is definitely full.
                if ((a_mask & b_mask) != kAllFacesMask) return true;

                // Otherwise both regions have area 2*Pi and intersect all 6 cube faces.
                // We choose "empty" in this case under the assumption that it is more
                // likely that the user is computing the difference between two nearly
                // identical polygons.
                //
                // TODO(ericv): Implement a robust algorithm based on examining the edge
                // snapping results directly, or alternatively add another heuristic (such
                // as testing containment of random points, or using a larger bit mask in
                // the tests above, e.g. a 24-bit mask representing all level 1 cells).
                return false;
            }
            return error_sign > 0;
        }

        // Returns a bit mask indicating which of the 6 S2 cube faces intersect the
        // index contents.
        private static byte GetFaceMask(S2ShapeIndex index)
        {
            byte mask = 0;
            var pos = 0;
            var count = index.GetEnumerableCount();
            while (pos < count)
            {
                var cellId = index.GetCellId(pos)!.Value;
                int face = (int)cellId.Face();
                mask |= (byte)(1 << face);
                var (pos2, _) = index.SeekCell(S2CellId.FromFace(face + 1).RangeMin());
                pos = pos2;
            }
            return mask;
        }

        // A bit mask representing all six faces of the S2 cube.
        private const byte kAllFacesMask = 0x3f;

        private readonly S2BooleanOperation op_;

        // The S2Builder options used to construct the output.
        S2Builder.Options builder_options_;

        // The S2Builder used to construct the output.  Note that the S2Builder
        // object is created only when is_boolean_output() is false.
        private S2Builder builder_;

        // A vector specifying the dimension of each edge added to S2Builder.
        private readonly List<sbyte> input_dimensions_ = new();

        // The set of all input edge crossings, which is used by EdgeClippingLayer
        // to construct the clipped output polygon.
        private readonly InputEdgeCrossings input_crossings_ = new();

        // kSentinel is a sentinel value used to mark the end of vectors.
        private static readonly Edge kSentinel = new(int.MaxValue, 0);

        // A vector containing all pairs of crossing edges from the two input
        // regions (including edge pairs that share a common vertex).  The first
        // element of each pair is an edge from "index_crossings_first_region_id_",
        // while the second element of each pair is an edge from the other region.
        private readonly List<IndexCrossing> index_crossings_ = new();

        // Indicates that the first element of each crossing edge pair in
        // "index_crossings_" corresponds to an edge from the given region.
        // This field is negative if index_crossings_ has not been computed yet.
        private int index_crossings_first_region_id_;

        // Temporary storage used in GetChainStarts(), declared here to avoid
        // repeatedly allocating memory.
        private readonly List<IndexCrossing> tmp_crossings_ = new();

        // An object to track the memory usage of this class.
        private readonly MemoryTracker tracker_;
    }
}

// CrossingInputEdge represents an input edge B that crosses some other input
// edge A.  It stores the input edge id of edge B and also whether it crosses
// edge A from left to right (or vice versa).
//
// Constructor:
// Indicates that input edge "input_id" crosses another edge (from left to
// right if "left_to_right" is true).
public readonly record struct CrossingInputEdge(InputEdgeId InputId, bool LeftToRight)
    : IComparable<CrossingInputEdge>, IComparable<InputEdgeId>
{
    public int CompareTo(CrossingInputEdge other) => InputId.CompareTo(other.InputId);

    public static bool operator <(CrossingInputEdge x, CrossingInputEdge y) => x.InputId.CompareTo(y.InputId) < 0;
    public static bool operator >(CrossingInputEdge x, CrossingInputEdge y) => x.InputId.CompareTo(y.InputId) > 0;
    public static bool operator <=(CrossingInputEdge x, CrossingInputEdge y) => x.InputId.CompareTo(y.InputId) <= 0;
    public static bool operator >=(CrossingInputEdge x, CrossingInputEdge y) => x.InputId.CompareTo(y.InputId) >= 0;

    public int CompareTo(InputEdgeId other) => InputId.CompareTo(other);

    public static bool operator <(CrossingInputEdge x, InputEdgeId y) => x.InputId.CompareTo(y) < 0;
    public static bool operator >(CrossingInputEdge x, InputEdgeId y) => x.InputId.CompareTo(y) > 0;
    public static bool operator <=(CrossingInputEdge x, InputEdgeId y) => x.InputId.CompareTo(y) <= 0;
    public static bool operator >=(CrossingInputEdge x, InputEdgeId y) => x.InputId.CompareTo(y) >= 0;
}

// Given two input edges A and B that intersect, suppose that A maps to a
// chain of snapped edges A_0, A_1, ..., A_m and B maps to a chain of snapped
// edges B_0, B_1, ..., B_n.  CrossingGraphEdge represents an edge from chain
// B that shares a vertex with chain A.  It is used as a temporary data
// representation while processing chain A.  The arguments are:
//
//   "id" - the Int32 of an edge from chain B.
//   "a_index" - the index of the vertex (A_i) that is shared with chain A.
//   "outgoing" - true if the shared vertex is the first vertex of the B edge.
//   "dst" - the Int32 of the vertex that is not shared with chain A.
//
// Note that if an edge from the B chain shares both vertices with the A
// chain, there will be two entries: an outgoing edge that treats its first
// vertex as being shared, and an incoming edge that treats its second vertex
// as being shared.
public readonly struct CrossingGraphEdge
{
    public readonly int Id;
    public readonly int AIndex;
    public readonly bool Outgoing;
    public readonly int Dst;
    public CrossingGraphEdge(int id, int a_index, bool outgoing, int dst)
    { Id = id; AIndex = a_index; Outgoing = outgoing; Dst = dst; }
}

// Given a set of clipping instructions encoded as a set of InputEdgeCrossings,
// GraphEdgeClipper determines which graph edges correspond to clipped
// portions of input edges and removes them.
//
// The clipping model is as follows.  The input consists of edge chains.  The
// clipper maintains an "inside" boolean state as it clips each chain, and
// toggles this state whenever an input edge is crossed.  Any edges that are
// deemed to be "outside" after clipping are removed.
//
// The "inside" state can be reset when necessary (e.g., when jumping to the
// start of a new chain) by adding a special crossing marked kSetInside.
// There are also two other special "crossings" that modify the clipping
// parameters: kSetInvertB specifies that edges should be clipped to the
// exterior of the other region, and kSetReverseA specifies that edges should
// be reversed before emitting them (which is needed to implement difference
// operations).
public class GraphEdgeClipper
{
    // "input_dimensions" is a vector specifying the dimension of each input
    // edge (0, 1, or 2).  "input_crossings" is the set of all crossings to be
    // used when clipping the edges of "g", sorted in lexicographic order.
    //
    // The clipped set of edges and their corresponding set of input edge ids
    // are returned in "new_edges" and "new_input_edge_ids".  (These can be used
    // to construct a new S2Builder.Graph.)
    public GraphEdgeClipper(Graph g, List<sbyte> input_dimensions, InputEdgeCrossings input_crossings, List<Edge> new_edges, List<int> new_input_edge_ids)
    {
        g_ = g; in_ = new Graph.VertexInMap(g); out_ = new Graph.VertexOutMap(g);
        input_dimensions_ = input_dimensions;
        input_crossings_ = input_crossings;
        new_edges_ = new_edges;
        new_input_edge_ids_ = new_input_edge_ids;
        input_ids_ = g.InputEdgeIdSetIds;
        order_ = GetInputEdgeChainOrder(g_, input_ids_);
        rank_ = new int[order_.Count];
        for (int i = 0; i < order_.Count; ++i)
        {
            rank_[order_[i]] = i;
        }
        // new_edges_ is obtained by filtering the graph edges and therefore the
        // number of graph edges is an upper bound on its size.
        new_edges_.EnsureCapacity(g_.NumEdges);
        new_input_edge_ids_.EnsureCapacity(g_.NumEdges);
    }


    // Returns a vector of EdgeIds sorted by input edge id.  When more than one
    // output edge has the same input edge id (i.e., the input edge snapped to a
    // chain of edges), the edges are sorted so that they form a directed edge
    // chain.
    //
    // This function could possibily be moved to S2Builder.Graph, but note that
    // it has special requirements.  Namely, duplicate edges and sibling pairs
    // must be kept in order to ensure that every output edge corresponds to
    // exactly one input edge.  (See also S2Builder.Graph.GetInputEdgeOrder.)
    private static List<int> GetInputEdgeChainOrder(Graph g, List<InputEdgeId> input_ids)
    {
        System.Diagnostics.Debug.Assert(g.Options.EdgeType_ == EdgeType.DIRECTED);
        System.Diagnostics.Debug.Assert(g.Options.DuplicateEdges_ == DuplicateEdges.KEEP);
        System.Diagnostics.Debug.Assert(g.Options.SiblingPairs_ == SiblingPairs.KEEP);

        // First, sort the edges so that the edges corresponding to each input edge
        // are consecutive.  (Each input edge was snapped to a chain of output
        // edges, or two chains in the case of undirected input edges.)
        var order = Graph.GetInputEdgeOrder(input_ids);

        // Now sort the group of edges corresponding to each input edge in edge
        // chain order (e.g.  AB, BC, CD).
        var vmap = new List<Int32Int32>();     // Map from source vertex to edge id.
        var indegree = new int[g.NumVertices];  // Restricted to current input edge.
        for (int end, begin = 0; begin < order.Count; begin = end)
        {
            // Gather the edges that came from a single input edge.
            var input_id = input_ids[order[begin]];
            for (end = begin; end < order.Count; ++end)
            {
                if (input_ids[order[end]] != input_id) break;
            }
            if (end - begin == 1) continue;

            // Build a map from the source vertex of each edge to its edge id,
            // and also compute the indegree at each vertex considering only the edges
            // that came from the current input edge.
            for (int i = begin; i < end; ++i)
            {
                var e = order[i];
                vmap.Add(new Int32Int32(g.GetEdge(e).ShapeId, e));
                indegree[g.GetEdge(e).EdgeId] += 1;
            }
            vmap.Sort();

            // Find the starting edge for building the edge chain.
            var next = g.NumEdges;
            for (int i = begin; i < end; ++i)
            {
                int e = order[i];
                if (indegree[g.GetEdge(e).ShapeId] == 0) next = e;
            }
            // Build the edge chain.
            for (int i = begin; ;)
            {
                order[i] = next;
                int v = g.GetEdge(next).EdgeId;
                indegree[v] = 0;  // Clear as we go along.
                if (++i == end) break;
                var index = vmap.GetLowerBound(new(v, 0));
                var output = vmap[index];
                System.Diagnostics.Debug.Assert(v == output.Item1);
                next = output.Item2;
            }
            vmap.Clear();
        }
        return order;
    }

    public void Run()
    {
        // Declare vectors here and reuse them to avoid reallocation.
        var a_vertices = new List<int>();
        var a_num_crossings = new List<int>();
        var a_isolated = new List<bool>();
        var b_input_edges = new List<CrossingInputEdge>();
        var b_edges = new List<CrossingGraphEdgeVector>();

        bool inside = false;
        bool invert_b = false;
        bool reverse_a = false;
        var next = 0;
        for (int i = 0; i < order_.Count; ++i)
        {
            // For each input edge (the "A" input edge), gather all the input edges
            // that cross it (the "B" input edges).
            int a_input_id = input_ids_[order_[i]];
            var edge0 = g_.GetEdge(order_[i]);
            b_input_edges.Clear();
            for (; next < input_crossings_.Count; ++next)
            {
                if (input_crossings_[next].Item1 != a_input_id) break;
                if (input_crossings_[next].Item2.InputId >= 0)
                {
                    b_input_edges.Add(input_crossings_[next].Item2);
                }
                else if (input_crossings_[next].Item2.InputId == S2BooleanOperation.kSetInside)
                {
                    inside = input_crossings_[next].Item2.LeftToRight;
                }
                else if (input_crossings_[next].Item2.InputId == S2BooleanOperation.kSetInvertB)
                {
                    invert_b = input_crossings_[next].Item2.LeftToRight;
                }
                else
                {
                    System.Diagnostics.Debug.Assert(input_crossings_[next].Item2.InputId == S2BooleanOperation.kSetReverseA);
                    reverse_a = input_crossings_[next].Item2.LeftToRight;
                }
            }
            // Optimization for degenerate edges.
            // TODO(ericv): If the output layer for this edge dimension specifies
            // DegenerateEdges.DISCARD, then remove the edge here.
            if (edge0.ShapeId == edge0.EdgeId)
            {
                inside ^= (b_input_edges.Count & 1) != 0;
                AddEdge(edge0, a_input_id);
                continue;
            }
            // Optimization for the case where there are no crossings.
            if (!b_input_edges.Any())
            {
                // In general the caller only passes edges that are part of the output
                // (i.e., we could System.Diagnostics.Debug.Assert(inside) here).  The one exception is for
                // polyline/polygon operations, where the polygon edges are needed to
                // compute the polyline output but are not emitted themselves.
                if (inside)
                {
                    AddEdge(reverse_a ? Graph.Reverse(edge0) : edge0, a_input_id);
                }
                continue;
            }
            // Walk along the chain of snapped edges for input edge A, and at each
            // vertex collect all the incident edges that belong to one of the
            // crossing edge chains (the "B" input edges).
            a_vertices.Clear();
            a_vertices.Add(edge0.ShapeId);
            b_edges.Clear();
            b_edges.Capacity = b_input_edges.Count;
            GatherIncidentEdges(a_vertices, 0, b_input_edges, b_edges);
            for (; i < order_.Count && input_ids_[order_[i]] == a_input_id; ++i)
            {
                a_vertices.Add(g_.GetEdge(order_[i]).EdgeId);
                GatherIncidentEdges(a_vertices, a_vertices.Count - 1, b_input_edges, b_edges);
            }
            --i;
#if s2builder_verbose
            System.Diagnostics.Debug.Write($"input edge {a_input_id} (inside={inside}):");
            foreach (var id in a_vertices) System.Diagnostics.Debug.Write(" " + id);
#endif
            // Now for each B edge chain, decide which vertex of the A chain it
            // crosses, and keep track of the number of signed crossings at each A
            // vertex.  The sign of a crossing depends on whether the other edge
            // crosses from left to right or right to left.
            //
            // This would not be necessary if all calculations were done in exact
            // arithmetic, because crossings would have strictly alternating signs.
            // But because we have already snapped the result, some crossing locations
            // are ambiguous, and GetCrossedVertexIndex() handles this by choosing a
            // candidate vertex arbitrarily.  The end result is that rarely, we may
            // see two crossings in a row with the same sign.  We correct for this by
            // adding extra output edges that essentially link up the crossings in the
            // correct (alternating sign) order.  Compared to the "correct" behavior,
            // the only difference is that we have added some extra sibling pairs
            // (consisting of an edge and its corresponding reverse edge) which do not
            // affect the result.
            a_num_crossings.Clear();
            a_num_crossings.Capacity = a_vertices.Count;
            a_isolated.Clear();
            a_isolated.Capacity = a_vertices.Count;
            for (int bi = 0; bi < b_input_edges.Count; ++bi)
            {
                bool left_to_right = b_input_edges[bi].LeftToRight;
                int a_index = GetCrossedVertexIndex(a_vertices, b_edges[bi],
                                                    left_to_right);
                if (a_index >= 0)
                {
#if s2builder_verbose
                    System.Diagnostics.Debug.Write($"{Environment.NewLine}  b input edge {b_input_edges[bi].InputId} (l2r={left_to_right}, crossing={a_vertices[a_index]})");

                    foreach (var x in b_edges[bi])
                    {
                        var e = g_.GetEdge(x.Id);
                        System.Diagnostics.Debug.Write($" ({e.ShapeId}, {e.EdgeId})");
                    }
#endif
                    // Keep track of the number of signed crossings (see above).
                    bool is_line = input_dimensions_[b_input_edges[bi].InputId] == 1;
                    int sign = is_line ? 0 : (left_to_right == invert_b) ? -1 : 1;
                    a_num_crossings[a_index] += sign;

                    // Any polyline or polygon vertex that has at least one crossing but no
                    // adjacent emitted edge may be emitted as an isolated vertex.
                    a_isolated[a_index] = true;
                }
                else
                {
                    // TODO(b/112043775): fix this condition.
                    throw new ApplicationException("Failed to get crossed vertex index.");
                }
            }
#if s2builder_verbose
            System.Diagnostics.Debug.WriteLine("");
#endif

            // Finally, we iterate through the A edge chain, keeping track of the
            // number of signed crossings as we go along.  The "multiplicity" is
            // defined as the cumulative number of signed crossings, and indicates how
            // many edges should be output (and in which direction) in order to link
            // up the edge crossings in the correct order.  (The multiplicity is
            // almost always either 0 or 1 except in very rare cases.)
            int multiplicity = a_num_crossings[0] + (inside ? 1 : 0);
            for (int ai = 1; ai < a_vertices.Count; ++ai)
            {
                if (multiplicity != 0)
                {
                    a_isolated[ai - 1] = a_isolated[ai] = false;
                }
                int edge_count = reverse_a ? -multiplicity : multiplicity;
                // Output any forward edges required.
                for (var i0 = 0; i0 < edge_count; ++i0)
                {
                    AddEdge(new Edge(a_vertices[ai - 1], a_vertices[ai]), a_input_id);
                }
                // Output any reverse edges required.
                for (var i1 = edge_count; i1 < 0; ++i1)
                {
                    AddEdge(new Edge(a_vertices[ai], a_vertices[ai - 1]), a_input_id);
                }
                multiplicity += a_num_crossings[ai];
            }
            // Multiplicities other than 0 or 1 can only occur in the edge interior.
            System.Diagnostics.Debug.Assert(multiplicity == 0 || multiplicity == 1);
            inside = (multiplicity != 0);

            // Output any isolated polyline vertices.
            // TODO(ericv): Only do this if an output layer wants degenerate edges.
            if (input_dimensions_[a_input_id] != 0)
            {
                for (int ai = 0; ai < a_vertices.Count; ++ai)
                {
                    if (a_isolated[ai])
                    {
                        AddEdge(new Edge(a_vertices[ai], a_vertices[ai]), a_input_id);
                    }
                }
            }
        }
    }

    private void AddEdge(Edge edge, InputEdgeId input_edge_id)
    {
        new_edges_.Add(edge);
        new_input_edge_ids_.Add(input_edge_id);
    }

    // Given the vertices of the snapped edge chain for an input edge A and the
    // set of input edges B that cross input edge A, this method gathers all of
    // the snapped edges of B that are incident to a given snapped vertex of A.
    // The incident edges for each input edge of B are appended to a separate
    // output vector.  (A and B can refer to either the input edge or the
    // corresponding snapped edge chain.)
    private void GatherIncidentEdges(List<int> a, int ai, List<CrossingInputEdge> b_input_edges, List<CrossingGraphEdgeVector> b_edges)
    {
        // Examine all of the edges incident to the given vertex of A.  If any edge
        // comes from a B input edge, append it to the appropriate vector.
        System.Diagnostics.Debug.Assert(b_input_edges.Count == b_edges.Count);
        foreach (var e in in_.EdgeIds(a[ai]))
        {
            int id = input_ids_[e];
            var it = b_input_edges.GetLowerBound(new CrossingInputEdge(id, false));
            if (it < b_input_edges.Count && b_input_edges[it].InputId == id)
            {
                var edges = b_edges[it];
                edges.Add(new CrossingGraphEdge(e, ai, false, g_.GetEdge(e).ShapeId));
            }
        }
        foreach (var e in out_.EdgeIds(a[ai]))
        {
            var id = input_ids_[e];
            var it = b_input_edges.GetLowerBound(new CrossingInputEdge(id, false));
            if (it < b_input_edges.Count && b_input_edges[it].InputId == id)
            {
                var edges = b_edges[it];
                edges.Add(new CrossingGraphEdge(e, ai, true, g_.GetEdge(e).EdgeId));
            }
        }
    }
    // Given an edge chain A that is crossed by another edge chain B (where
    // "left_to_right" indicates whether B crosses A from left to right), this
    // method decides which vertex of A the crossing takes place at.  The
    // parameters are the vertices of the A chain ("a") and the set of edges in
    // the B chain ("b") that are incident to vertices of A.  The B chain edges
    // are sorted in increasing order of (a_index, outgoing) tuple.
    private int GetCrossedVertexIndex(List<int> a, CrossingGraphEdgeVector b, bool left_to_right)
    {
        if (!a.Any() || !b.Any())
        {
            System.Diagnostics.Debug.WriteLine(
                $"GraphEdgeClipper.GetCrossedVertexIndex called with {a.Count} vertex ids and {b.Count} crossing graph edges.");
            return -1;
        }

        // The reason this calculation is tricky is that after snapping, the A and B
        // chains may meet and separate several times.  For example, if B crosses A
        // from left to right, then B may touch A, make an excursion to the left of
        // A, come back to A, then make an excursion to the right of A and come back
        // to A again, like this:
        //
        //  *--B--*-\             /-*-\
        //           B-\       /-B     B-\      6     7     8     9
        //  *--A--*--A--*-A,B-*--A--*--A--*-A,B-*--A--*--A--*-A,B-*
        //  0     1     2     3     4     5      \-B     B-/
        //                                          \-*-/
        //
        // (where "*" is a vertex, and "A" and "B" are edge labels).  Note that B
        // may also follow A for one or more edges whenever they touch (e.g. between
        // vertices 2 and 3).  In this case the only vertices of A where the
        // crossing could take place are 5 and 6, i.e. after all excursions of B to
        // the left of A, and before all excursions of B to the right of A.
        //
        // Other factors to consider are that the portion of B before and/or after
        // the crossing may be degenerate, and some or all of the B edges may be
        // reversed relative to the A edges.

        // First, check whether edge A is degenerate.
        int n = a.Count;
        if (n == 1) return 0;

        // If edge chain B is incident to only one vertex of A, we're done.
        if (b[0].AIndex == b.Last().AIndex) return b[0].AIndex;

        // Determine whether the B chain visits the first and last vertices that it
        // shares with the A chain in the same order or the reverse order.  This is
        // only needed to implement one special case (see below).
        bool b_reversed = GetVertexRank(b[0]) > GetVertexRank(b.Last());

        // Examine each incident B edge and use it to narrow the range of positions
        // where the crossing could occur in the B chain.  Vertex positions are
        // represented as a range [lo, hi] of vertex ranks in the B chain (see
        // GetVertexRank).
        //
        // Note that if an edge of B is incident to the first or last vertex of A,
        // we can't test which side of the A chain it is on.  (An S2Pred.Sign test
        // doesn't work; e.g. if the B edge is XY and the first edge of A is YZ,
        // then snapping can change the sign of XYZ while maintaining topological
        // guarantees.)  There can be up to 4 such edges (one incoming and one
        // outgoing edge at each endpoint of A).  Two of these edges logically
        // extend past the end of the A chain and place no restrictions on the
        // crossing vertex.  The other two edges define the ends of the subchain
        // where B shares vertices with A.  We save these edges in order to handle a
        // special case (see below).
        int lo = -1, hi = order_.Count;   // Vertex ranks of acceptable crossings
        var b_first = -1;
        var b_last = -1;  // "b" subchain connecting "a" endpoints
        foreach (var e in b)
        {
            int ai = e.AIndex;
            if (ai == 0)
            {
                if (e.Outgoing != b_reversed && e.Dst != a[1]) b_first = e.Id;
            }
            else if (ai == n - 1)
            {
                if (e.Outgoing == b_reversed && e.Dst != a[n - 2]) b_last = e.Id;
            }
            else
            {
                // This B edge is incident to an interior vertex of the A chain.  First
                // check whether this edge is identical (or reversed) to an edge in the
                // A chain, in which case it does not create any restrictions.
                if (e.Dst == a[ai - 1] || e.Dst == a[ai + 1]) continue;

                // Otherwise we can test which side of the A chain the edge lies on.
                bool on_left = S2Pred.OrderedCCW(g_.Vertex(a[ai + 1]), g_.Vertex(e.Dst),
                                                  g_.Vertex(a[ai - 1]), g_.Vertex(a[ai]));

                // Every B edge that is incident to an interior vertex of the A chain
                // places some restriction on where the crossing vertex could be.
                if (left_to_right == on_left)
                {
                    // This is a pre-crossing edge, so the crossing cannot be before the
                    // destination vertex of this edge.  (For example, the input B edge
                    // crosses the input A edge from left to right and this edge of the B
                    // chain is to the left of the A chain.)
                    lo = Math.Max(lo, rank_[e.Id] + 1);
                }
                else
                {
                    // This is a post-crossing edge, so the crossing cannot be after the
                    // source vertex of this edge.
                    hi = Math.Min(hi, rank_[e.Id]);
                }
            }
        }
        // There is one special case.  If a subchain of B connects the first and
        // last vertices of A, then together with the edges of A this forms a loop
        // whose orientation can be tested to determine whether B is on the left or
        // right side of A.  This is only possible (and only necessary) if the B
        // subchain does not include any interior vertices of A, since otherwise the
        // B chain might cross from one side of A to the other.
        //
        // Note that it would be possible to avoid this test in some situations by
        // checking whether either endpoint of the A chain has two incident B edges,
        // in which case we could check which side of the B chain the A edge is on
        // and use this to limit the possible crossing locations.
        if (b_first >= 0 && b_last >= 0)
        {
            // Swap the edges if necessary so that they are in B chain order.
            if (b_reversed)
            {
                (b_last, b_first) = (b_first, b_last);
            }

            // The B subchain connects the first and last vertices of A.  We test
            // whether the chain includes any interior vertices of A by iterating
            // through the incident B edges again, looking for ones that belong to
            // the B subchain and are not incident to the first or last vertex of A.
            bool has_interior_vertex = false;
            foreach (var e in b)
            {
                if (e.AIndex > 0 && e.AIndex < n - 1 &&
                    rank_[e.Id] >= rank_[b_first] && rank_[e.Id] <= rank_[b_last])
                {
                    has_interior_vertex = true;
                    break;
                }
            }
            if (!has_interior_vertex)
            {
                // The B subchain is not incident to any interior vertex of A.
                bool on_left = EdgeChainOnLeft(a, b_first, b_last);
                if (left_to_right == on_left)
                {
                    lo = Math.Max(lo, rank_[b_last] + 1);
                }
                else
                {
                    hi = Math.Min(hi, rank_[b_first]);
                }
            }
        }

        // Otherwise we choose the smallest shared Int32 in the acceptable range,
        // in order to ensure that both chains choose the same crossing vertex.
        int best = -1;
        System.Diagnostics.Debug.Assert(lo <= hi);
        foreach (var e in b)
        {
            int ai = e.AIndex;
            int vrank = GetVertexRank(e);
            if (vrank >= lo && vrank <= hi && (best < 0 || a[ai] < a[best]))
            {
                best = ai;
            }
        }
        return best;
    }

    // Returns the "vertex rank" of the shared vertex associated with the given
    // CrossingGraphEdge.  Recall that graph edges are sorted in input edge order,
    // and that the rank of an edge is its position in this order (rank_[e]).
    // VertexRank(e) is defined such that VertexRank(e.src) == rank_[e] and
    // VertexRank(e.dst) == rank_[e] + 1.  Note that the concept of "vertex rank"
    // is only defined within a single edge chain (since different edge chains can
    // have overlapping vertex ranks).
    private int GetVertexRank(CrossingGraphEdge e)
    {
        return rank_[e.Id] + (e.Outgoing ? 0 : 1);
    }
    // Given edge chains A and B that form a loop (after possibly reversing the
    // direction of chain B), returns true if chain B is to the left of chain A.
    // Chain A is given as a sequence of vertices, while chain B is specified as
    // the first and last edges of the chain.
    private bool EdgeChainOnLeft(List<int> a, EdgeId b_first, EdgeId b_last)
    {
        // Gather all the interior vertices of the B subchain.
        var loop = new List<int>();
        for (int i = rank_[b_first]; i < rank_[b_last]; ++i)
        {
            loop.Add(g_.GetEdge(order_[i]).EdgeId);
        }
        // Possibly reverse the chain so that it forms a loop when "a" is appended.
        if (g_.GetEdge(b_last).EdgeId != a[0]) loop.Reverse();
        loop.AddRange(a);
        // Duplicate the first two vertices to simplify vertex indexing.
        loop.Add(loop[0]);
        loop.Add(loop[1]);
        // Now B is to the left of A if and only if the loop is counterclockwise.
        double sum = 0;
        for (int i = 2; i < loop.Count; ++i)
        {
            sum += S2.TurnAngle(g_.Vertex(loop[i - 2]), g_.Vertex(loop[i - 1]),
                                 g_.Vertex(loop[i]));
        }
        return sum > 0;
    }

    private readonly Graph g_;
    private readonly Graph.VertexInMap in_;
    private readonly Graph.VertexOutMap out_;
    private readonly List<sbyte> input_dimensions_;
    private readonly InputEdgeCrossings input_crossings_;
    private readonly List<Edge> new_edges_;
    private readonly List<int> new_input_edge_ids_;

    // Every graph edge is associated with exactly one input edge in our case,
    // which means that we can declare g_.input_edge_id_set_ids() as a vector of
    // Int32s rather than a vector of Int32s.  (This also takes
    // advantage of the fact that IdSetLexicon represents a singleton set as the
    // value of its single element.)
    private readonly List<InputEdgeId> input_ids_;
    private readonly List<EdgeId> order_;  // Graph edges sorted in input edge id order.
    private readonly int[] rank_;      // The rank of each graph edge within order_.
}

// Given a set of clipping instructions encoded as a set of intersections
// between input edges, EdgeClippingLayer determines which graph edges
// correspond to clipped portions of input edges and removes them.  It
// assembles the remaining edges into a new S2Builder.Graph and passes the
// result to the given output layer for assembly.
public class EdgeClippingLayer : Layer
{
    public EdgeClippingLayer(
        List<Layer> layers,
        List<sbyte> input_dimensions,
        InputEdgeCrossings input_crossings,
        S2MemoryTracker.Client tracker)
    {
        layers_ = layers;
        input_dimensions_ = input_dimensions;
        input_crossings_ = input_crossings;
        tracker_ = tracker;
    }

    // Layer interface:
    public override GraphOptions GraphOptions_()
    {
        // We keep all edges, including degenerate ones, so that we can figure out
        // the correspondence between input edge crossings and output edge
        // crossings.
        return new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
                                DuplicateEdges.KEEP, SiblingPairs.KEEP);
    }
    public override void Build(Graph g, out S2Error error)
    {
        // Data per graph edge:
        //   vector<EdgeId> order_;
        //   vector<int> rank_;
        //   vector<Graph.Edge> new_edges;
        //   vector<InputEdgeIdSetId> new_input_edge_ids;
        // Data per graph vertex:
        //   Graph.VertexInMap in_;
        //   Graph.VertexOutMap out_;
        //
        // The first and last two vectors above are freed upon GraphEdgeClipper
        // destruction.  There is also a temporary vector "indegree" in
        // GetInputEdgeChainOrder() but this does not affect peak memory usage.
        Int64 tmp_bytes = g.NumEdges * (sizeof(EdgeId) + sizeof(int)) +
                          g.NumVertices * (2 * sizeof(EdgeId));
        Int64 final_bytes = g.NumEdges * (Marshal.SizeOf(typeof(Edge)) +
                                             Marshal.SizeOf(typeof(InputEdgeIdSetId)));

        // The order of the calls below is important.  Note that all memory tracked
        // through this client is automatically untallied upon object destruction.
        if (!tracker_.Tally(final_bytes) || !tracker_.TallyTemp(tmp_bytes))
        {
            error = S2Error.OK;
            // We don't need to copy memory tracking errors to "error" because this
            // is already done for us in S2BooleanOperation.Impl.Build().
            return;
        }
        error = S2Error.OK;
        // The bulk of the work is handled by GraphEdgeClipper.
        var new_edges = new List<Edge>();
        var new_input_edge_ids = new List<EdgeId>();
        // Destroy the GraphEdgeClipper immediately to save memory.
        new GraphEdgeClipper(g, input_dimensions_, input_crossings_,
            new_edges, new_input_edge_ids).Run();
#if s2builder_verbose
        System.Diagnostics.Debug.WriteLine("Edges after clipping: ");
        for (int i = 0; i < new_edges.Count; ++i)
        {
            System.Diagnostics.Debug.WriteLine(
                $"  {new_input_edge_ids[i]} ({new_edges[i].ShapeId}, {new_edges[i].EdgeId})");
        }
#endif
        // Construct one or more subgraphs from the clipped edges and pass them to
        // the given output layer(s).  We start with a copy of the input graph's
        // IdSetLexicon because this is necessary in general, even though in this
        // case it is guaranteed to be empty because no edges have been merged.
        var new_input_edge_id_set_lexicon = g.InputEdgeIdSetLexicon;
        if (layers_.Count == 1)
        {
            Graph new_graph = g.MakeSubgraph(
                layers_[0].GraphOptions_(), new_edges, new_input_edge_ids,
                new_input_edge_id_set_lexicon, g.IsFullPolygonPredicate(),
                out error, tracker_)!;
            if (tracker_.Ok()) layers_[0].Build(new_graph, out error);
            tracker_.Untally(new_edges);
            tracker_.Untally(new_input_edge_ids);
        }
        else
        {
            // The Graph objects must be valid until the last Build() call completes,
            // so we store all of the graph data in arrays with 3 elements.
            System.Diagnostics.Debug.Assert(3 == layers_.Count);
            List<Edge>[] layer_edges = { new List<Edge>(), new List<Edge>(), new List<Edge>(), };
            List<EdgeId>[] layer_input_edge_ids = { new List<EdgeId>(), new List<EdgeId>(), new List<EdgeId>(), };
            // Separate the edges according to their dimension.
            for (int i = 0; i < new_edges.Count; ++i)
            {
                int d = input_dimensions_[new_input_edge_ids[i]];
                if (!tracker_.AddSpace(layer_edges[d], 1)) return;
                if (!tracker_.AddSpace(layer_input_edge_ids[d], 1)) return; layer_edges[d].Add(new_edges[i]);
                layer_input_edge_ids[d].Add(new_input_edge_ids[i]);
            }
            // Clear variables to save space.
            new_edges.Clear();
            new_input_edge_ids.Clear();

            var layer_graphs = new List<Graph>(3);

            for (int d = 0; d < 3; ++d)
            {
                layer_graphs.Add(g.MakeSubgraph(
                    layers_[d].GraphOptions_(), layer_edges[d],
                    layer_input_edge_ids[d], new_input_edge_id_set_lexicon,
                    g.IsFullPolygonPredicate(), out error, tracker_)!);
                if (tracker_.Ok()) layers_[d].Build(layer_graphs[d], out error);
            }
            for (int d = 0; d < 3; ++d)
            {
                tracker_.Untally(layer_edges[d]);
                tracker_.Untally(layer_input_edge_ids[d]);
            }
        }
    }

    private readonly List<Layer> layers_;
    private readonly List<sbyte> input_dimensions_;
    private readonly InputEdgeCrossings input_crossings_;
    private readonly S2MemoryTracker.Client tracker_;
}
