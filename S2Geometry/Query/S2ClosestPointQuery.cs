using System;
using System.Collections.Generic;
using Distance = S2Geometry.S1ChordAngle;
using Target = S2Geometry.S2DistanceTarget<S2Geometry.S1ChordAngle>;

namespace S2Geometry
{
    // Given a set of points stored in an S2PointIndex, S2ClosestPointQuery
    // provides methods that find the closest point(s) to a given query point
    // or query edge.  Example usage:
    //
    // void Test(S2Point[] index_points,
    //           S2Point[] target_points) {
    //   // The template argument allows auxiliary data to be attached to each
    //   // point (in this case, the array index).
    //   S2PointIndex<int> index;
    //   for (int i = 0; i < index_points.size(); ++i) {
    //     index.Add(index_points[i], i);
    //   }
    //   S2ClosestPointQuery<int> query(out index);
    //   query.Options().set_max_results(5);
    //   for (S2Point target_point : target_points) {
    //     S2ClosestPointQueryPointTarget target(target_point);
    //     foreach (var& result in query.FindClosestPoints(out target)) {
    //       // The Result class contains the following methods:
    //       //   distance() is the distance to the target.
    //       //   point() is the indexed point.
    //       //   data() is the auxiliary data.
    //       DoSomething(target_point, result);
    //     }
    //   }
    // }
    //
    // You can find either the k closest points, or all points within a given
    // radius, or both (i.e., the k closest points up to a given maximum radius).
    // E.g. to find all the points within 5 kilometers, call
    //
    //   query.Options().set_max_distance(
    //       S2Earth.ToAngle(util.units.Kilometers(5)));
    //
    // By default *all* points are returned, so you should always specify either
    // max_results() or max_distance() or both.  There is also a FindClosestPoint()
    // convenience method that returns only the closest point.
    //
    // You can restrict the results to an arbitrary S2Region, for example:
    //
    //   S2LatLngRect rect(...);
    //   query.Options().set_region(out rect);  // Does *not* take ownership.
    //
    // To find the closest points to a query edge rather than a point, use:
    //
    //   S2ClosestPointQueryEdgeTarget target(v0, v1);
    //   query.FindClosestPoints(out target);
    //
    // Similarly you can find the closest points to an S2Cell by using an
    // S2ClosestPointQuery.CellTarget, and you can find the closest points to an
    // arbitrary collection of points, polylines, and polygons by using an
    // S2ClosestPointQuery.ShapeIndexTarget.
    //
    // The implementation is designed to be fast for both small and large
    // point sets.
    //
    // See S2ClosestPointQueryBase for full documentation.
    public class S2ClosestPointQuery<Data> where Data : IComparable<Data>
    {
        // Options that control the set of points returned.  Note that by default
        // *all* points are returned, so you will always want to set either the
        // max_results() option or the max_distance() option (or both).
        //
        // This class is also available as S2ClosestPointQuery<Data>.Options.
        // (It is defined here to avoid depending on the "Data" template argument.)
        public class Options : S2ClosestPointQueryBase<Distance, Data>.Options
        {
            // See S2ClosestPointQueryBaseOptions for the full set of options.

            // Specifies that only points whose distance to the target is less than
            // "max_distance" should be returned.
            //
            // Note that points whose distance is exactly equal to "max_distance" are
            // not returned.  Normally this doesn't matter, because distances are not
            // computed exactly in the first place, but if such points are needed then
            // see set_inclusive_max_distance() below.
            //
            // DEFAULT: Distance.Infinity()
            //public S1ChordAngle MaxDistance { set => base.MaxDistance = value; }

            // Like set_max_distance(), except that points whose distance is exactly
            // equal to "max_distance" are also returned.  Equivalent to calling
            // set_max_distance(max_distance.Successor()).
            public S1ChordAngle InclusiveMaxDistance { set => MaxDistance = (value.Successor); }

            // Like set_inclusive_max_distance(), except that "max_distance" is also
            // increased by the maximum error in the distance calculation.  This ensures
            // that all points whose true distance is less than or equal to
            // "max_distance" will be returned (along with some points whose true
            // distance is slightly greater).
            //
            // Algorithms that need to do exact distance comparisons can use this
            // option to find a set of candidate points that can then be filtered
            // further (e.g., using S2Pred.CompareDistance).
            public S1ChordAngle ConservativeMaxDistance { set => MaxDistance = (value.PlusError(
                S2EdgeDistances.GetUpdateMinDistanceMaxError(value)).Successor); }
        }

        // Target subtype that computes the closest distance to a point.
        //
        // This class is also available as S2ClosestPointQuery<Data>.PointTarget.
        // (It is defined here to avoid depending on the "Data" template argument.)
        public sealed class PointTarget : S2MinDistancePointTarget
        {
            public PointTarget(S2Point point)
            : base(point) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindClosest (which finds the single closest point), the
                // break-even points are approximately X, Y, and Z points for grid,
                // fractal, and regular loop geometry respectively.
                //
                // TODO(ericv): Adjust using benchmarks.
                150;
        }

        // Target subtype that computes the closest distance to an edge.
        //
        // This class is also available as S2ClosestPointQuery<Data>.EdgeTarget.
        // (It is defined here to avoid depending on the "Data" template argument.)
        public sealed class EdgeTarget : S2MinDistanceEdgeTarget
        {
            public EdgeTarget(S2Point a, S2Point b) : base(a, b) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindClosestToEdge (which finds the single closest point), the
                // break-even points are approximately X, Y, and Z points for grid,
                // fractal, and regular loop geometry respectively.
                //
                // TODO(ericv): Adjust using benchmarks.
                100;
        }

        // Target subtype that computes the closest distance to an S2Cell
        // (including the interior of the cell).
        //
        // This class is also available as S2ClosestPointQuery<Data>.CellTarget.
        // (It is defined here to avoid depending on the "Data" template argument.)
        public sealed class CellTarget : S2MinDistanceCellTarget
        {
            public CellTarget(S2Cell cell) : base(cell) { }
            public override int MaxBruteForceIndexSize =>
                // Using BM_FindClosestToCell (which finds the single closest point), the
                // break-even points are approximately X, Y, and Z points for grid,
                // fractal, and regular loop geometry respectively.
                //
                // TODO(ericv): Adjust using benchmarks.
                50;
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
        //
        // This class is also available as S2ClosestPointQuery<Data>.ShapeIndexTarget.
        // (It is defined here to avoid depending on the "Data" template argument.)
        public sealed class ShapeIndexTarget : S2MinDistanceShapeIndexTarget
        {
            public ShapeIndexTarget(S2ShapeIndex index) : base(index) { }
            public override int MaxBruteForceIndexSize =>
                // For BM_FindClosestToSameSizeAbuttingIndex (which uses a nearby
                // S2ShapeIndex target of similar complexity), the break-even points are
                // approximately X, Y, and Z points for grid, fractal, and regular loop
                // geometry respectively.
                //
                // TODO(ericv): Adjust using benchmarks.
                30;
        }

        // Convenience constructor that calls Init().  Options may be specified here
        // or changed at any time using the Options() accessor method.
        public S2ClosestPointQuery(S2PointIndex<Data> index, Options options = null)
        {
            Init(index, options ?? new Options());
        }

        // Default constructor; requires Init() to be called.
        public S2ClosestPointQuery() { }

        // Initializes the query.  Options may be specified here or changed at any
        // time using the Options() accessor method.
        //
        // REQUIRES: "index" must persist for the lifetime of this object.
        // REQUIRES: ReInit() must be called if "index" is modified.
        public void Init(S2PointIndex<Data> index, Options options = null)
        {
            Options_ = options ?? new Options();
            base_.Init(index);
        }

        // Reinitializes the query.  This method must be called whenever the
        // underlying index is modified.
        public void ReInit()
        {
            base_.ReInit();
        }

        // Returns a reference to the underlying S2PointIndex.
        public S2PointIndex<Data> Index()
        {
            return base_.Index;
        }

        // Returns the query options.  Options can be modifed between queries.
        public Options Options_;

        // Returns the closest points to the given target that satisfy the current
        // options.  This method may be called multiple times.
        public List<S2ClosestPointQueryBase<Distance, Data>.Result> FindClosestPoints(Target target)
        {
            return base_.FindClosestPoints(target, Options_);
        }

        // This version can be more efficient when this method is called many times,
        // since it does not require allocating a new vector on each call.
        public void FindClosestPoints(Target target, List<S2ClosestPointQueryBase<Distance, Data>.Result> results)
        {
            base_.FindClosestPoints(target, Options_, results);
        }

        //////////////////////// Convenience Methods ////////////////////////

        // Returns the closest point to the target.  If no point satisfies the search
        // criteria, then a Result object with distance() == Infinity() and
        // IsEmpty == true is returned.
        public S2ClosestPointQueryBase<S1ChordAngle, Data>.Result FindClosestPoint(Target target)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            return base_.FindClosestPoint(target, tmp_options);
        }

        // Returns the minimum distance to the target.  If the index or target is
        // empty, returns S1ChordAngle.Infinity.
        //
        // Use IsDistanceLess() if you only want to compare the distance against a
        // threshold value, since it is often much faster.
        public Distance GetDistance(Target target)
        {
            return (Distance)FindClosestPoint(target).Distance;
        }

        // Returns true if the distance to "target" is less than "limit".
        //
        // This method is usually much faster than GetDistance(), since it is much
        // less work to determine whether the minimum distance is above or below a
        // threshold than it is to calculate the actual minimum distance.
        public bool IsDistanceLess(Target target, S1ChordAngle limit)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.MaxDistance = (limit);
            tmp_options.MaxError = (S1ChordAngle.Straight);
            return !base_.FindClosestPoint(target, tmp_options).IsEmpty;
        }

        // Like IsDistanceLess(), but also returns true if the distance to "target"
        // is exactly equal to "limit".
        public bool IsDistanceLessOrEqual(Target target, S1ChordAngle limit)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.InclusiveMaxDistance = (limit);
            tmp_options.MaxError = (S1ChordAngle.Straight);
            return !base_.FindClosestPoint(target, tmp_options).IsEmpty;
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
        public bool IsConservativeDistanceLessOrEqual(Target target, S1ChordAngle limit)
        {
            // Assert.True(Marshal.SizeOf(typeof(Options)) <= 32); // Consider not copying Options here
            Options tmp_options = Options_;
            tmp_options.MaxResults = (1);
            tmp_options.ConservativeMaxDistance = (limit);
            tmp_options.MaxError = (S1ChordAngle.Straight);
            return !base_.FindClosestPoint(target, tmp_options).IsEmpty;
        }

        private readonly S2ClosestPointQueryBase<Distance, Data> base_ = new();
    }
}
