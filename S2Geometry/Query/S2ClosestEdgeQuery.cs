using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace S2Geometry
{
    // S2MinDistance is a thin wrapper around S1ChordAngle that implements the
    // Distance concept required by S2ClosestPointQueryBase.
    using Distance = S1ChordAngle; // S2MinDistance
    // Each "Result" object represents a closest edge.  It has the following
    // fields:
    //
    //   S1ChordAngle distance;  // The distance from the target to this edge.
    //   Int32 shape_id;         // Identifies an indexed shape.
    //   Int32 edge_id;          // Identifies an edge within the shape.
    using Result = S2ClosestEdgeQueryBase<S1ChordAngle>.Result;
    // "Target" represents the geometry to which the distance is measured.
    // There are subtypes for measuring the distance to a point, an edge, an
    // S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
    // using S2MinDistance = S1ChordAngle;
    // using S2MinDistanceTarget = S2DistanceTarget<S2MinDistance>;
    using Target = S2DistanceTarget<S1ChordAngle>; // S2MinDistanceTarget

    // S2ClosestEdgeQuery is a helper class for finding the closest edge(s) to a
    // given point, edge, S2Cell, or geometry collection.  For example, given a
    // set of polylines, the following code efficiently finds the closest 5 edges
    // to a query point:
    //
    // void Test(S2Polyline[] polylines, S2Point point) {
    //   MutableS2ShapeIndex index;
    //   for (S2Polyline polyline : polylines) {
    //     index.Add(new S2Polyline.Shape(polyline));
    //   }
    //   S2ClosestEdgeQuery query(out index);
    //   query.Options().set_max_results(5);
    //   S2ClosestEdgeQuery.PointTarget target(point);
    //   foreach (var& result in query.FindClosestEdges(out target)) {
    //     // The Result struct contains the following fields:
    //     //   "distance" is the distance to the edge.
    //     //   "shape_id" identifies the S2Shape containing the edge.
    //     //   "edge_id" identifies the edge with the given shape.
    //     // The following convenience methods may also be useful:
    //     //   query.GetEdge(result) returns the endpoints of the edge.
    //     //   query.Project(point, result) computes the closest point on the
    //     //       result edge to the given target point.
    //     int polyline_index = result.shape_id;
    //     int edge_index = result.edge_id;
    //     S1ChordAngle distance = result.distance;  // Use ToAngle() for S1Angle.
    //     S2Shape.Edge edge = query.GetEdge(result);
    //     S2Point closest_point = query.Project(point, result);
    //   }
    // }
    //
    // You can find either the k closest edges, or all edges within a given
    // radius, or both (i.e., the k closest edges up to a given maximum radius).
    // E.g. to find all the edges within 5 kilometers, call
    //
    //   query.Options().set_max_distance(
    //       S2Earth.ToAngle(util.units.Kilometers(5)));
    //
    // By default *all* edges are returned, so you should always specify either
    // max_results() or max_distance() or both.  There is also a FindClosestEdge()
    // convenience method that returns only the closest edge.
    //
    // Note that by default, distances are measured to the boundary and interior
    // of polygons.  For example, if a point is inside a polygon then its distance
    // is zero.  To change this behavior, call set_include_interiors(false).
    //
    // If you only need to test whether the distance is above or below a given
    // threshold (e.g., 10 km), you can use the IsDistanceLess() method.  This is
    // much faster than actually calculating the distance with FindClosestEdge(),
    // since the implementation can stop as soon as it can prove that the minimum
    // distance is either above or below the threshold.
    //
    // To find the closest edges to a query edge rather than a point, use:
    //
    //   S2ClosestEdgeQuery.EdgeTarget target(v0, v1);
    //   query.FindClosestEdges(out target);
    //
    // Similarly you can find the closest edges to an S2Cell by using an
    // S2ClosestEdgeQuery.CellTarget, and you can find the closest edges to an
    // arbitrary collection of points, polylines, and polygons by using an
    // S2ClosestEdgeQuery.ShapeIndexTarget.
    //
    // The implementation is designed to be fast for both simple and complex
    // geometric objects.
    public class S2ClosestEdgeQuery
    {
        // See S2ClosestEdgeQueryBase for full documentation.

        // Options that control the set of edges returned.  Note that by default
        // *all* edges are returned, so you will always want to set either the
        // max_results() option or the max_distance() option (or both).
        //
        // See S2ClosestEdgeQueryBase.Options for the full set of options.
        public class Options : S2ClosestEdgeQueryBase<Distance>.Options
        {
            // Like set_max_distance(), except that edges whose distance is exactly
            // equal to "max_distance" are also returned.  Equivalent to calling
            // set_max_distance(max_distance.Successor()).
            public Distance InclusiveMaxDistance { set => MaxDistance = value.Successor; }

            // Like set_inclusive_max_distance(), except that "max_distance" is also
            // increased by the maximum error in the distance calculation.  This
            // ensures that all edges whose true distance is less than or equal to
            // "max_distance" will be returned (along with some edges whose true
            // distance is slightly greater).
            //
            // Algorithms that need to do exact distance comparisons can use this
            // option to find a set of candidate edges that can then be filtered
            // further (e.g., using S2Pred.CompareDistance).
            public Distance ConservativeMaxDistance { set => MaxDistance = value.PlusError(
                S2EdgeDistances.GetUpdateMinDistanceMaxError(value)).Successor; }
        }

        // Target subtype that computes the closest distance to a point.
        public sealed class PointTarget : S2MinDistancePointTarget
        {
            public PointTarget(S2Point point) : base(point) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindClosest (which finds the single closest edge), the
                // break-even points are approximately 80, 100, and 250 edges for point
                // cloud, fractal, and regular loop geometry respectively.
                120;
        }

        // Target subtype that computes the closest distance to an edge.
        public sealed class EdgeTarget : S2MinDistanceEdgeTarget
        {
            public EdgeTarget(S2Point a, S2Point b) : base(a, b) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindClosestToEdge (which finds the single closest edge), the
                // break-even points are approximately 40, 50, and 100 edges for point
                // cloud, fractal, and regular loop geometry respectively.
                60;
        }

        // Target subtype that computes the closest distance to an S2Cell
        // (including the interior of the cell).
        public sealed class CellTarget : S2MinDistanceCellTarget
        {
            public CellTarget(S2Cell cell) : base(cell) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindClosestToCell (which finds the single closest edge), the
                // break-even points are approximately 20, 25, and 40 edges for point cloud,
                // fractal, and regular loop geometry respectively.
                30;
        }

        // Target subtype that computes the closest distance to an S2ShapeIndex
        // (an arbitrary collection of points, polylines, and/or polygons).
        //
        // By default, distances are measured to the boundary and interior of
        // polygons in the S2ShapeIndex rather than to polygon boundaries only.
        // If you wish to change this behavior, you may call
        //
        //   target.set_include_interiors(false);
        //
        // (see S2MinDistanceShapeIndexTarget for details).
        public sealed class ShapeIndexTarget : S2MinDistanceShapeIndexTarget
        {
            public ShapeIndexTarget(S2ShapeIndex index) : base(index) { }
            public override int MaxBruteForceIndexSize =>
                // For BM_FindClosestToSameSizeAbuttingIndex (which uses two nearby indexes
                // with similar edge counts), the break-even points are approximately 20,
                // 30, and 40 edges for point cloud, fractal, and regular loop geometry
                // respectively.
                25;
        };

        // Convenience constructor that calls Init().  Options may be specified here
        // or changed at any time using the Options() accessor method.
        public S2ClosestEdgeQuery(S2ShapeIndex index, Options options = null)
        {
            Init(index, options ?? new Options());
        }

        // Default constructor; requires Init() to be called.
        public S2ClosestEdgeQuery() { }

        // Initializes the query.  Options may be specified here or changed at any
        // time using the Options() accessor method.
        //
        // REQUIRES: "index" must persist for the lifetime of this object.
        // REQUIRES: ReInit() must be called if "index" is modified.
        public void Init(S2ShapeIndex index, Options options = null)
        {
            Options_ = options ?? new Options();
            base_.Init(index);
        }

        // Reinitializes the query.  This method must be called whenever the
        // underlying S2ShapeIndex is modified.
        public void ReInit()
        {
            base_.ReInit();
        }

        // Returns a reference to the underlying S2ShapeIndex.
        public S2ShapeIndex Index()
        {
            return base_.Index;
        }

        // Options_ can be modified between queries.
        public Options Options_ { get; private set; }

        // Returns the closest edges to the given target that satisfy the current
        // options.  This method may be called multiple times.
        //
        // Note that if Options_.include_interiors() is true, the result vector may
        // include some entries with edge_id == -1.  This indicates that the target
        // intersects the indexed polygon with the given shape_id.
        public List<Result> FindClosestEdges(Target target)
        {
            return base_.FindClosestEdges(target, Options_);
        }

        // This version can be more efficient when this method is called many times,
        // since it does not require allocating a new vector on each call.
        public void FindClosestEdges(Target target, List<Result> results)
        {
            base_.FindClosestEdges(target, Options_, results);
        }

        //////////////////////// Convenience Methods ////////////////////////

        // Returns the closest edge to the target.  If no edge satisfies the search
        // criteria, then the Result object will have distance == Infinity(),
        // IsEmpty == true, and shape_id == edge_id == -1.
        //
        // Note that if options.include_interiors() is true, edge_id == -1 is also
        // used to indicate that the target intersects an indexed polygon (but in
        // that case distance == Zero() and shape_id >= 0).
        public Result FindClosestEdge(Target target)
        {
            // Assert.True(Marshal.SizeOf(Options_) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            return base_.FindClosestEdge(target, tmp_options);
        }

        // Returns the minimum distance to the target.  If the index or target is
        // empty, returns S1ChordAngle.Infinity.
        //
        // Use IsDistanceLess() if you only want to compare the distance against a
        // threshold value, since it is often much faster.
        public Distance GetDistance(Target target)
        {
            return FindClosestEdge(target).Distance;
        }

        // Returns true if the distance to "target" is less than "limit".
        //
        // This method is usually much faster than GetDistance(), since it is much
        // less work to determine whether the minimum distance is above or below a
        // threshold than it is to calculate the actual minimum distance.
        public bool IsDistanceLess(Target target, Distance limit)
        {
            Assert.True(Marshal.SizeOf(Options_) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.MaxDistance = (limit);
            tmp_options.MaxError = (Distance.Straight);
            return !base_.FindClosestEdge(target, tmp_options).IsEmpty;
        }

        // Like IsDistanceLess(), but also returns true if the distance to "target"
        // is exactly equal to "limit".
        public bool IsDistanceLessOrEqual(Target target, Distance limit)
        {
            Assert.True(Marshal.SizeOf(Options_) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.InclusiveMaxDistance = (limit);
            tmp_options.MaxError = (Distance.Straight);
            return !base_.FindClosestEdge(target, tmp_options).IsEmpty;
        }

        // Like IsDistanceLessOrEqual(), except that "limit" is increased by the
        // maximum error in the distance calculation.  This ensures that this
        // function returns true whenever the true, exact distance is less than
        // or equal to "limit".
        //
        // For example, suppose that we want to test whether two geometries might
        // intersect each other after they are snapped together using S2Builder
        // (using the IdentitySnapFunction with a given "snap_radius").  Since
        // S2Builder uses exact distance predicates (s2predicates.h), we need to
        // measure the distance between the two geometries conservatively.  If the
        // distance is definitely greater than "snap_radius", then the geometries
        // are guaranteed to not intersect after snapping.
        public bool IsConservativeDistanceLessOrEqual(Target target, Distance limit)
        {
            Assert.True(Marshal.SizeOf(Options_) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.ConservativeMaxDistance = (limit);
            tmp_options.MaxError = (Distance.Straight);
            return !base_.FindClosestEdge(target, tmp_options).IsEmpty;
        }

        // Returns the endpoints of the given result edge.
        //
        // CAVEAT: If Options_.include_interiors() is true, then clients must not
        // pass this method any Result objects that correspond to shape interiors,
        // i.e. those where result.edge_id < 0.
        //
        // REQUIRES: result.edge_id >= 0
        public S2Shape.Edge GetEdge(Result result)
        {
            return Index().Shape(result.ShapeId).GetEdge(result.EdgeId);
        }

        // Returns the point on given result edge that is closest to "point".
        public S2Point Project(S2Point point, Result result)
        {
            if (result.EdgeId< 0) return point;
            var edge = GetEdge(result);
            return S2EdgeDistances.Project(point, edge.V0, edge.V1);
        }

        private readonly S2ClosestEdgeQueryBase<Distance> base_ = new S2ClosestEdgeQueryBase<Distance>();
    }
}
