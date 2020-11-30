using System;
using Distance = S2Geometry.S2MaxDistance;
using Base = S2Geometry.S2ClosestEdgeQueryBase<S2Geometry.S2MaxDistance>;

// "Target" represents the geometry to which the distance is measured.
// There are subtypes for measuring the distance to a point, an edge, an
// S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).  Note
// that S2DistanceTarget<Distance> is equivalent to S2MaxDistanceTarget in
// s2max_distance_targets, which the following subtypes
// (e.g. S2MaxDistancePointTarget) extend.
using Target = S2Geometry.S2DistanceTarget<S2Geometry.S2MaxDistance>;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Diagnostics.CodeAnalysis;

namespace S2Geometry
{
    // S2FurthestEdgeQuery is a helper class for searching within an S2ShapeIndex
    // for the furthest edge(s) to a given query point, edge, S2Cell, or geometry
    // collection.  The furthest edge is defined as the one which maximizes the
    // distance from any point on that edge to any point on the target geometry.
    //
    // As an example, given a set of polylines, the following code efficiently
    // finds the furthest 5 edges to a query point:
    //
    // void Test(S2Polyline[] polylines, S2Point point) {
    //   MutableS2ShapeIndex index;
    //   for (S2Polyline* polyline : polylines) {
    //     index.Add(new S2Polyline.Shape(polyline));
    //   }
    //   S2FurthestEdgeQuery query(out index);
    //   query.Options().set_max_results(5);
    //   S2FurthestEdgeQuery.PointTarget target(point);
    //   foreach (var& result in query.FindFurthestEdges(out target)) {
    //     // The Result struct contains the following accessors:
    //     //   "distance()" is the distance to the edge.
    //     //   "shape_id()" identifies the S2Shape containing the edge.
    //     //   "edge_id()" identifies the edge with the given shape.
    //     // The following convenience methods may also be useful:
    //     //   query.GetEdge(result) returns the endpoints of the edge.
    //     int polyline_index = result.shape_id();
    //     int edge_index = result.edge_id();
    //     S1ChordAngle distance = result.distance();
    //     S2Shape.Edge edge = query.GetEdge(result);
    //   }
    // }
    //
    // You can find either the k furthest edges, or all edges no closer than a
    // given radius, or both (i.e., the k furthest edges no closer than a given
    // minimum radius).
    // E.g. to find all the edges further than 5 kilometers, call
    //
    //   query.Options().set_min_distance(
    //       S2Earth.ToAngle(util.units.Kilometers(5)));
    //
    // By default *all* edges are returned, so you should always specify either
    // max_results() or min_distance() or both.  Setting min distance may not be
    // very restrictive, so strongly consider using max_results().  There is also a
    // FindFurthestEdge() convenience method that returns only the single furthest
    // edge.
    //
    // Note that by default, distances are measured to the boundary and interior
    // of polygons.  To change this behavior, call set_include_interiors(false).
    //
    // If you only need to test whether the distance is above or below a given
    // threshold (e.g., 10 km), you can use the IsDistanceGreater() method.  This
    // is much faster than actually calculating the distance with
    // FindFurthestEdge(), since the implementation can stop as soon as it can
    // prove that the maximum distance is either above or below the threshold.
    //
    // To find the furthest edges to a query edge rather than a point, use:
    //
    //   S2FurthestEdgeQuery.EdgeTarget target(v0, v1);
    //   query.FindFurthestEdges(out target);
    //
    // Similarly you can find the furthest edges to an S2Cell by using an
    // S2FurthestEdgeQuery.CellTarget, and you can find the furthest edges to an
    // arbitrary collection of points, polylines, and polygons by using an
    // S2FurthestEdgeQuery.ShapeIndexTarget.
    //
    // The implementation is designed to be fast for both simple and complex
    // geometric objects.
    //
    // See S2FurthestEdgeQueryBase for full documentation.
    public class S2FurthestEdgeQuery
    {
        // Options that control the set of edges returned.  Note that by default
        // *all* edges are returned, so you will always want to set either the
        // max_results() option or the min_distance() option (or both).
        public class Options : Base.Options
        {
            // The following methods are like the corresponding max_distance() methods
            // in S2ClosestEdgeQuery, except that they specify the minimum distance to
            // the target.
            public S1ChordAngle MinDistance { get => MaxDistance.Distance; set => MaxDistance = new(value); }

            public S1ChordAngle InclusiveMinDistance { set => MinDistance = value.Predecessor; }

            public S1ChordAngle ConservativeMinDistance { set => MaxDistance = (new(value.PlusError(
                    -S2EdgeDistances.GetUpdateMinDistanceMaxError(value)).Predecessor)); }
        }


        // Target subtype that computes the furthest distance to a point.
        public sealed class PointTarget : S2MaxDistancePointTarget
        {
            public PointTarget(S2Point point) : base(point) { }
            // See s2closest_edge_query.cc for justifications of
            // max_brute_force_index_size() for that query.
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindFurthest (which finds the single furthest edge), the
                // break-even points are approximately 100, 400, and 600 edges for point
                // cloud, fractal, and regular loop geometry respectively.
                300;
        }

        // Target subtype that computes the furthest distance to an edge.
        public sealed class EdgeTarget : S2MaxDistanceEdgeTarget
        {
            public EdgeTarget(S2Point a, S2Point b) : base(a, b) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindFurthestToEdge (which finds the single furthest edge), the
                // break-even points are approximately 80, 100, and 230 edges for point
                // cloud, fractal, and regular loop geometry respectively.
                110;
        }

        // Target subtype that computes the furthest distance to an S2Cell
        // (including the interior of the cell).
        public sealed class CellTarget : S2MaxDistanceCellTarget
        {
            public CellTarget(S2Cell cell) : base(cell) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindFurthestToCell (which finds the single furthest edge), the
                // break-even points are approximately 70, 100, and 170 edges for point
                // cloud, fractal, and regular loop geometry respectively.
                100;
        }

        // Target subtype that computes the furthest distance to an S2ShapeIndex
        // (an arbitrary collection of points, polylines, and/or polygons).
        //
        // By default, distances are measured to the boundary and interior of
        // polygons in the S2ShapeIndex rather than to polygon boundaries only.
        // If you wish to change this behavior, you may call
        //
        //   target.set_include_interiors(false);
        //
        // (see S2MaxDistanceShapeIndexTarget for details).
        public sealed class ShapeIndexTarget : S2MaxDistanceShapeIndexTarget
        {
            public ShapeIndexTarget(S2ShapeIndex index) : base(index) { }
            public override int MaxBruteForceIndexSize =>
                // For BM_FindFurthestToSameSizeAbuttingIndex (which uses two nearby indexes
                // with similar edge counts), the break-even points are approximately 30,
                // 100, and 130 edges for point cloud, fractal, and regular loop geometry
                // respectively.
                70;
        }

        // Result class to pass back to user.  We choose to pass back this result
        // type, which has an S1ChordAngle as its distance, rather than the
        // Base.Result returned from the query which uses S2MaxDistance.
        public class Result : IEquatable<Result>, IComparable<Result>
        {
            // The defaultructor, which yields an invalid result.
            public Result() : this(S1ChordAngle.Negative, -1, -1) { }

            // Construct a Result from a Base.Result.
            public Result(Base.Result base_)
                : this(base_.Distance.Distance, base_.ShapeId, base_.EdgeId) { }

            // Constructs a Result object for the given edge with the given distance.
            public Result(S1ChordAngle distance, Int32 _shape_id, Int32 _edge_id)
            {
                distance_ = distance; shape_id_ = _shape_id; edge_id_ = _edge_id;
            }

            // Returns true if this Result object does not refer to any edge.
            // (The only case where an empty Result is returned is when the
            // FindFurthestEdges() method does not find any edges that meet the
            // specified criteria.)
            public bool IsEmpty => edge_id_ < 0;

            // The distance from the target to this point.
            public S1ChordAngle Distance() { return distance_; }

            // The edge identifiers.
            public Int32 ShapeId() { return shape_id_; }
            public Int32 EdgeId() { return edge_id_; }

            public override bool Equals(object obj)
            {
                return obj is Result res && Equals(res);
            }
            public bool Equals(Result other)
            {
                return distance_ == other.distance_ &&
                    shape_id_ == other.shape_id_ &&
                    edge_id_ == other.edge_id_;
            }
            public override int GetHashCode()
            {
                return HashCode.Combine(distance_, shape_id_, edge_id_);
            }

            // Returns true if two Result objects are identical.
            public static bool operator ==(Result x, Result y)
            {
                return Equals(x, y);
            }
            public static bool operator !=(Result x, Result y)
            {
                return !Equals(x, y);
            }

            public int CompareTo([AllowNull] Result other)
            {
                return this < other ? -1 : this > other ? 1 : 0;
            }

            // Compares edges first by distance, then by (shape_id, edge_id).
            public static bool operator <(Result x, Result y)
            {
                if (x.distance_ < y.distance_) return true;
                if (y.distance_ < x.distance_) return false;
                if (x.shape_id_ < y.shape_id_) return true;
                if (y.shape_id_ < x.shape_id_) return false;
                return x.edge_id_ < y.edge_id_;
            }
            public static bool operator >(Result x, Result y)
            {
                return !(x < y || Equals(x, y));
            }

            private readonly S1ChordAngle distance_;
            private readonly Int32 shape_id_;     // Identifies an indexed shape.
            private readonly Int32 edge_id_;      // Identifies an edge within the shape.
        }

        // Convenience constructor that calls Init().  Options may be specified here
        // or changed at any time using the Options() accessor method.
        public S2FurthestEdgeQuery(S2ShapeIndex index, Options options = null)
        {
            Init(index, options ?? new Options());
        }

        // Default constructor; requires Init() to be called.
        public S2FurthestEdgeQuery() { }

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

        // Returns the query options.  Options can be modified between queries.
        public Options Options_ { get; set; }

        // Returns the furthest edges to the given target that satisfy the given
        // options.  This method may be called multiple times.
        public List<Result> FindFurthestEdges(Target target)
        {
            var results = new List<Result>();
            FindFurthestEdges(target, results);
            return results;
        }

        // This version can be more efficient when this method is called many times,
        // since it does not require allocating a new vector on each call.
        public void FindFurthestEdges(Target target, List<Result> results)
        {
            results.Clear();
            foreach (var result in base_.FindClosestEdges(target, Options_))
            {
                results.Add(new S2FurthestEdgeQuery.Result(result));
            }
        }

        //////////////////////// Convenience Methods ////////////////////////

        // Returns the furthest edge to the target.  If no edge satisfies the search
        // criteria, then the Result object will have
        // distance == S1ChordAngle.Negative and shape_id == edge_id == -1.
        public Result FindFurthestEdge(Target target)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            Base.Result base_result = base_.FindClosestEdge(target, tmp_options);
            return new S2FurthestEdgeQuery.Result(base_result);
        }

        // Returns the maximum distance to the target.  If the index or target is
        // empty, returns S1ChordAngle.Negative.
        //
        // Use IsDistanceGreater() if you only want to compare the distance against a
        // threshold value, since it is often much faster.
        public S1ChordAngle GetDistance(Target target)
        {
            return FindFurthestEdge(target).Distance();
        }

        // Returns true if the distance to "target" is greater than "limit".
        //
        // This method is usually much faster than GetDistance(), since it is much
        // less work to determine whether the maximum distance is above or below a
        // threshold than it is to calculate the actual maximum distance.
        public bool IsDistanceGreater(Target target, S1ChordAngle limit)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.MinDistance = (limit);
            tmp_options.MaxError = (S1ChordAngle.Straight);
            return base_.FindClosestEdge(target, tmp_options).ShapeId>= 0;
        }

        // Like IsDistanceGreater(), but also returns true if the distance to
        // "target" is exactly equal to "limit".
        public bool IsDistanceGreaterOrEqual(Target target, S1ChordAngle limit)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.InclusiveMinDistance = (limit);
            tmp_options.MaxError = (S1ChordAngle.Straight);
            return base_.FindClosestEdge(target, tmp_options).ShapeId>= 0;
        }

        // Like IsDistanceGreaterOrEqual(), except that "limit" is decreased by the
        // maximum error in the distance calculation.  This ensures that this
        // function returns true whenever the true, exact distance is greater than
        // or equal to "limit".
        public bool IsConservativeDistanceGreaterOrEqual(Target target, S1ChordAngle limit)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.ConservativeMinDistance = (limit);
            tmp_options.MaxError = (S1ChordAngle.Straight);
            return base_.FindClosestEdge(target, tmp_options).ShapeId>= 0;
        }

        // Returns the endpoints of the given result edge.
        // REQUIRES: result.edge_id >= 0
        public S2Shape.Edge GetEdge(Result result)
        {
            return Index().Shape(result.ShapeId()).GetEdge(result.EdgeId());
        }

        private readonly Base base_ = new Base();
    }
}
