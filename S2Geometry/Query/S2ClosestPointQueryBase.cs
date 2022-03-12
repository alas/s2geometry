// S2ClosestPointQueryBase is a templatized class for finding the closest
// point(s) to a given target.  It is not intended to be used directly, but
// rather to serve as the implementation of various specialized classes with
// more convenient APIs (such as S2ClosestPointQuery).  It is flexible enough
// so that it can be adapted to compute maximum distances and even potentially
// Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between a point collection A and a target B.
//  - Find all points in collection A that are within a distance D of target B.
//  - Find the k points of collection A that are closest to a given point P.
//
// The target is any class that implements the S2DistanceTarget interface.
// There are predefined targets for points, edges, S2Cells, and S2ShapeIndexes
// (arbitrary collctions of points, polylines, and polygons).
//
// The Distance template argument is used to represent distances.  Usually it
// is a thin wrapper around S1ChordAngle, but another distance type may be
// used as long as it implements the Distance concept described in
// s2distance_targets.h.  For example this can be used to measure maximum
// distances, to get more accuracy, or to measure non-spheroidal distances.

using System.Diagnostics.CodeAnalysis;

namespace S2Geometry;

using Delta = S1ChordAngle;

public class S2ClosestPointQueryBase<Distance, Data> where Distance : IDistance where Data : IComparable<Data>
{
    #region Fields, Constants

    private const int kMaxMaxResults = int.MaxValue;
    private static readonly IDistance Infinity = (IDistance)typeof(Distance).GetField("Infinity")!.GetValue(null)!;
    private static readonly IDistance Zero = (IDistance)typeof(Distance).GetField("Zero")!.GetValue(null)!;
    // The minimum number of points that a cell must contain to enqueue it
    // rather than processing its contents immediately.
    private const int kMinPointsToEnqueue = 13;

    // Return a reference to the underlying S2PointIndex.
    public S2PointIndex<Data> Index { get; private set; }
    private Options Options_ { get; set; }
    private readonly SortedSet<QueueEntry> queue_ = new();

    // Temporaries, defined here to avoid multiple allocations / initializations.
    private readonly List<S2CellId> intersection_with_region_ = new();
    private readonly List<S2CellId> intersection_with_max_distance_ = new();
    private readonly (S2Point, Data)[] tmp_point_data_ = new (S2Point, Data)[kMinPointsToEnqueue - 1];

    private S2DistanceTarget<Distance> target_;

    // True if max_error() must be subtracted from priority queue cell distances
    // in order to ensure that such distances are measured conservatively.  This
    // is true only if the target takes advantage of max_error() in order to
    // return faster results, and 0 < max_error() < distance_limit_.
    private bool use_conservative_cell_distance_;

    // For the optimized algorihm we precompute the top-level S2CellIds that
    // will be added to the priority queue.  There can be at most 6 of these
    // cells.  Essentially this is just a covering of the indexed points.
    private readonly List<S2CellId> index_covering_ = new();

    // The distance beyond which we can safely ignore further candidate points.
    // (Candidates that are exactly at the limit are ignored; this is more
    // efficient for UpdateMinDistance() and should not affect clients since
    // distance measurements have a small amount of error anyway.)
    //
    // Initially this is the same as the maximum distance specified by the user,
    // but it can also be updated by the algorithm (see MaybeAddResult).
    private Distance distance_limit_;

    // The current result set is stored in one of three ways:
    //
    //  - If max_results() == 1, the best result is kept in result_singleton_.
    //
    //  - If max_results() == "infinity", results are appended to result_vector_
    //    and sorted/uniqued at the end.
    //
    //  - Otherwise results are kept in a priority queue so that we can
    //    progressively reduce the distance limit once max_results() results
    //    have been found.
    private Result result_singleton_;
    private readonly List<Result> result_vector_ = new();
    private readonly SortedSet<Result> result_set_ = new();

    #endregion

    #region Constructors

    // Default constructor; requires Init() to be called.
    public S2ClosestPointQueryBase() { }

    // Convenience constructor that calls Init().
    public S2ClosestPointQueryBase(S2PointIndex<Data> index)
    {
        Init(index);
    }

    #endregion

    // Initializes the query.
    // REQUIRES: ReInit() must be called if "index" is modified.
    public void Init(S2PointIndex<Data> index)
    {
        Index = index;
        ReInit();
    }

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    public void ReInit()
    {
        //iter_ = new S2PointIndex<Data>.Iterator(index_);
        index_covering_.Clear();
    }

    // Returns the closest points to the given target that satisfy the given
    // options.  This method may be called multiple times.
    public List<Result> FindClosestPoints(S2DistanceTarget<Distance> target, Options options)
    {
        var results = new List<Result>();
        FindClosestPoints(target, options, results);
        return results;
    }

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    public void FindClosestPoints(S2DistanceTarget<Distance> target, Options options, List<Result> results)
    {
        FindClosestPointsInternal(target, options);
        results.Clear();
        if (options.MaxResults == 1)
        {
            if (!result_singleton_.IsEmpty)
            {
                results.Add(result_singleton_);
            }
        }
        else if (options.MaxResults == kMaxMaxResults)
        {
            result_vector_.Sort();
            results.Clear();
            results.AddRange(result_vector_);
            result_vector_.Clear();
        }
        else
        {
            results.Capacity = result_set_.Count;
            for (; result_set_.Any(); result_set_.Remove(result_set_.First()))
            {
                results.Add(result_set_.First());
            }
            // The priority queue returns the largest elements first.
            results.Reverse();
            System.Diagnostics.Debug.Assert(results.IsSorted());
        }
    }

    // Convenience method that returns exactly one point.  If no points satisfy
    // the given search criteria, then a Result with distance() == Infinity()
    // and IsEmpty == true is returned.
    //
    // REQUIRES: options.max_results() == 1
    public Result FindClosestPoint(S2DistanceTarget<Distance> target, Options options)
    {
        System.Diagnostics.Debug.Assert(options.MaxResults == 1);
        FindClosestPointsInternal(target, options);
        return result_singleton_;
    }

    private void FindClosestPointsInternal(S2DistanceTarget<Distance> target, Options options)
    {
        target_ = target;
        Options_ = options;

        distance_limit_ = (Distance)options.MaxDistance;
        result_singleton_ = new Result();
        System.Diagnostics.Debug.Assert(!result_vector_.Any());
        System.Diagnostics.Debug.Assert(!result_set_.Any());
        System.Diagnostics.Debug.Assert(target.MaxBruteForceIndexSize >= 0);
        if (Equals(distance_limit_, Zero)) return;

        // If max_error() > 0 and the target takes advantage of this, then we may
        // need to adjust the distance estimates to the priority queue cells to
        // ensure that they are always a lower bound on the true distance.  For
        // example, suppose max_distance == 100, max_error == 30, and we compute the
        // distance to the target from some cell C0 as d(C0) == 80.  Then because
        // the target takes advantage of max_error(), the true distance could be as
        // low as 50.  In order not to miss edges contained by such cells, we need
        // to subtract max_error() from the distance estimates.  This behavior is
        // controlled by the use_conservative_cell_distance_ flag.
        //
        // However there is one important case where this adjustment is not
        // necessary, namely when max_distance() < max_error().  This is because
        // max_error() only affects the algorithm once at least max_results() edges
        // have been found that satisfy the given distance limit.  At that point,
        // max_error() is subtracted from distance_limit_ in order to ensure that
        // any further matches are closer by at least that amount.  But when
        // max_distance() < max_error(), this reduces the distance limit to 0,
        // i.e. all remaining candidate cells and edges can safely be discarded.
        // (Note that this is how IsDistanceLess() and friends are implemented.)
        //
        // Note that Distance.Delta only supports operator==.
        bool target_uses_max_error = false;
        if (options.MaxError != Delta.Zero)
        {
            target_.MaxError = options.MaxError;
            target_uses_max_error = true;
        }

        // Note that we can't compare max_error() and distance_limit_ directly
        // because one is a Delta and one is a Distance.  Instead we subtract them.
        use_conservative_cell_distance_ = target_uses_max_error &&
            (Equals(distance_limit_, Infinity) ||
                Zero.IsLessThan(distance_limit_.SubstractChord(options.MaxError)));

        // Note that given point is processed only once (unlike S2ClosestEdgeQuery),
        // and therefore we don't need to worry about the possibility of having
        // duplicate points in the results.
        if (options.UseBruteForce ||
            Index.NumPoints() <= target_.MaxBruteForceIndexSize)
        {
            FindClosestPointsBruteForce();
        }
        else
        {
            FindClosestPointsOptimized();
        }
    }
    private void FindClosestPointsBruteForce()
    {
        for (var pos = 0; pos < Index.NumPoints(); pos++)
        {
            var it = Index.Get(pos);
            MaybeAddResult(it!.Value.Point, it.Value.Data);
        }
    }
    private void FindClosestPointsOptimized()
    {
        InitQueue();
        while (queue_.Any())
        {
            // We need to copy the top entry before removing it, and we need to remove
            // it before adding any new entries to the queue.
            var entry = queue_.First();
            queue_.Remove(entry);
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            IDistance distance = entry.Distance;
            if (!(distance.CompareTo(distance_limit_) < 0))
            {
                queue_.Clear();  // Clear any remaining entries.
                break;
            }
            S2CellId child = entry.Id.ChildBegin();
            // We already know that it has too many points, so process its children.
            // Each child may either be processed directly or enqueued again.  The
            // loop is optimized so that we don't seek unnecessarily.
            bool seek = true;
            var pos = 0;
            for (int i = 0; i < 4; ++i, child = child.Next())
            {
                seek = ProcessOrEnqueue(child, seek, ref pos);
            }
        }
    }
    private void InitQueue()
    {
        System.Diagnostics.Debug.Assert(!queue_.Any());

        // Optimization: rather than starting with the entire index, see if we can
        // limit the search region to a small disc.  Then we can find a covering for
        // that disc and intersect it with the covering for the index.  This can
        // save a lot of work when the search region is small.
        S2Cap cap = target_.GetCapBound();
        if (cap.IsEmpty()) return;  // Empty target.
        if (Options_.MaxResults == 1)
        {
            // If the user is searching for just the closest point, we can compute an
            // upper bound on search radius by seeking to the center of the target's
            // bounding cap and looking at the adjacent index points (in S2CellId
            // order).  The minimum distance to either of these points is an upper
            // bound on the search radius.
            //
            // TODO(ericv): The same strategy would also work for small values of
            // max_results() > 1, e.g. max_results() == 20, except that we would need to
            // examine more neighbors (at least 20, and preferably 20 in each
            // direction).  It's not clear whether this is a common case, though, and
            // also this would require extending MaybeAddResult() so that it can
            // remove duplicate entries.  (The points added here may be re-added by
            // ProcessOrEnqueue(), but this is okay when max_results() == 1.)
            var i = Index.Seek(new S2CellId(cap.Center));
            if (i < Index.NumPoints())
            {
                var it_ = Index.Get(i);
                MaybeAddResult(it_!.Value.Point, it_.Value.Data);
            }
            var it = Index.Get(--i);
            if (it != null)
            {
                MaybeAddResult(it.Value.Point, it.Value.Data);
            }
            // Skip the rest of the algorithm if we found a matching point.
            if (Equals(distance_limit_, Zero)) return;
        }
        // We start with a covering of the set of indexed points, then intersect it
        // with the given region (if any) and maximum search radius disc (if any).
        if (!index_covering_.Any()) InitCovering();
        var initial_cells = index_covering_;
        var region = Options_.Region;
        if (region != null)
        {
            var coverer = new S2RegionCoverer();
            coverer.Options_.MaxCells = 4;
            coverer.GetCovering(region, out var region_covering_);
            S2CellUnion.GetIntersection(index_covering_, region_covering_, intersection_with_region_);
            initial_cells = intersection_with_region_;
        }
        if (distance_limit_.IsLessThan(Infinity.ToS1ChordAngle()))
        {
            var coverer = new S2RegionCoverer();
            coverer.Options_.MaxCells = 4;
            Delta radius = cap.Radius + distance_limit_.GetChordAngleBound();
            var search_cap = new S2Cap(cap.Center, radius);
            var max_distance_covering_ = new List<S2CellId>();
            coverer.GetFastCovering(search_cap, max_distance_covering_);
            S2CellUnion.GetIntersection(initial_cells, max_distance_covering_, intersection_with_max_distance_);
            initial_cells = intersection_with_max_distance_;
        }
        var pos = 0;
        for (int i = 0; i < initial_cells.Count && pos < Index.NumPoints(); ++i)
        {
            S2CellId id = initial_cells[i];
            ProcessOrEnqueue(id, id.RangeMin() > Index.Get(pos)?.Id /*seek*/, ref pos);
        }
    }
    private void InitCovering()
    {
        // Compute the "index covering", which is a small number of S2CellIds that
        // cover the indexed points.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one covering cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a covering cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the covering cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.
        index_covering_.Capacity = 6;
        if (Index.NumPoints() == 0) return;

        var pos = 0;
        var index_last_id = Index.Get(Index.NumPoints() - 1)!.Value.Id;
        var index_first_id = Index.Get(0)!.Value.Id;
        if (index_first_id != index_last_id)
        {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            int level = index_first_id.CommonAncestorLevel(index_last_id) + 1;

            // Visit each potential covering cell except the last (handled below).
            var last_id = index_last_id.Parent(level);
            var id = index_first_id.Parent(level);
            for (; id != last_id; id = id.Next())
            {
                var current = Index.Get(pos)!.Value.Id;

                // Skip any covering cells that don't contain any index cells.
                if (id.RangeMax() < current) continue;

                // Find the range of index cells contained by this covering cell and
                // then shrink the cell if necessary so that it just covers them.
                var cell_first_id = current;
                pos = Index.Seek(id.RangeMax().Next()) - 1;
                current = Index.Get(pos)!.Value.Id;
                S2CellId cell_last_id = current;
                pos++;
                AddInitialRange(cell_first_id, cell_last_id);
            }
        }
        AddInitialRange(Index.Get(pos)!.Value.Id, index_last_id);
    }

    // Adds a cell to index_covering_ that covers the given inclusive range.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private void AddInitialRange(S2CellId first_id, S2CellId last_id)
    {
        // Add the lowest common ancestor of the given range.
        int level = first_id.CommonAncestorLevel(last_id);
        System.Diagnostics.Debug.Assert(level >= 0);
        index_covering_.Add(first_id.Parent(level));
    }
    private void MaybeAddResult(S2Point point, Data data)
    {
        var distance = distance_limit_;
        if (!target_.UpdateMinDistance(point, ref distance)) return;

        var region = Options_.Region;
        if (region != null && !region.Contains(point)) return;

        var result = new Result(distance, point, data);
        if (Options_.MaxResults == 1)
        {
            // Optimization for the common case where only the closest point is wanted.
            result_singleton_ = result;
            distance_limit_ = (Distance)result.Distance.Substract(Options_.MaxError);
        }
        else if (Options_.MaxResults == kMaxMaxResults)
        {
            result_vector_.Add(result);  // Sort/unique at end.
        }
        else
        {
            // Add this point to result_set_.  Note that with the current algorithm
            // each candidate point is considered at most once (except for one special
            // case where max_results() == 1, see InitQueue for details), so we don't
            // need to worry about possibly adding a duplicate entry here.
            if (result_set_.Count >= Options_.MaxResults)
            {
                result_set_.Remove(result_set_.First());  // Replace the furthest result point.
            }
            result_set_.Add(result);
            if (result_set_.Count >= Options_.MaxResults)
            {
                distance_limit_ = (Distance)result_set_.First().Distance.Substract(Options_.MaxError);
            }
        }
    }

    // Either process the contents of the given cell immediately, or add it to the
    // queue to be subdivided.  If "seek" is false, then "iter" must already be
    // positioned at the first indexed point within or after this cell.
    //
    // Returns "true" if the cell was added to the queue, and "false" if it was
    // processed immediately, in which case "iter" is left positioned at the next
    // cell in S2CellId order.
    private bool ProcessOrEnqueue(S2CellId id, bool seek, ref int pos)
    {
        if (seek) pos = Index.Seek(id.RangeMin());
        if (id.IsLeaf())
        {
            TreeNode<Data>? cur;
            // Leaf cells can't be subdivided.
            for (; pos < Index.NumPoints() && (cur = Index.Get(pos))?.Id == id; pos++)
            {
                MaybeAddResult(cur.Value.Point, cur.Value.Data);
            }
            return false;  // No need to seek to next child.
        }
        S2CellId last = id.RangeMax();
        int num_points = 0;
        TreeNode<Data>? item;
        for (; pos < Index.NumPoints() && (item = Index.Get(pos))?.Id <= last; pos++)
        {
            if (num_points == kMinPointsToEnqueue - 1)
            {
                // This cell has too many points (including this one), so enqueue it.
                var cell = new S2Cell(id);
                var distance = distance_limit_;
                // We check "region_" second because it may be relatively expensive.
                if (target_.UpdateMinDistance(cell, ref distance) &&
                    (Options_.Region == null || Options_.Region.MayIntersect(cell)))
                {
                    if (use_conservative_cell_distance_)
                    {
                        // Ensure that "distance" is a lower bound on distance to the cell.
                        distance = (Distance)distance.Substract(Options_.MaxError);
                    }
                    queue_.Add(new QueueEntry(distance, id));
                }
                return true;  // Seek to next child.
            }
            tmp_point_data_[num_points++] = (item.Value.Point, item.Value.Data);
        }
        // There were few enough points that we might as well process them now.
        for (int i = 0; i < num_points; ++i)
        {
            var item2 = tmp_point_data_[i];
            MaybeAddResult(item2.Item1, item2.Item2);
        }
        return false;  // No need to seek to next child.
    }

    // The Target class represents the geometry to which the distance is
    // measured.  For example, there can be subtypes for measuring the distance
    // to a point, an edge, or to an S2ShapeIndex (an arbitrary collection of
    // geometry).
    //
    // Implementations do *not* need to be thread-safe.  They may cache data or
    // allocate temporary data structures in order to improve performance.
    // using Target = S2DistanceTarget<Distance>;

    // Options that control the set of points returned.  Note that by default
    // *all* points are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    public class Options
    {
        // Specifies that at most "max_results" points should be returned.
        //
        // REQUIRES: max_results >= 1
        // DEFAULT: numeric_limits<int>.max()
        public int MaxResults
        {
            get => _maxResults;
            set
            {
                System.Diagnostics.Debug.Assert(value >= 1);
                _maxResults = value;
            }
        }
        private int _maxResults = kMaxMaxResults;

        // Specifies that only points whose distance to the target is less than
        // "max_distance" should be returned.
        //
        // Note that points whose distance is exactly equal to "max_distance" are
        // not returned.  In most cases this doesn't matter (since distances are
        // not computed exactly in the first place), but if such points are needed
        // then you can retrieve them by specifying "max_distance" as the next
        // largest representable Distance.  For example, if Distance is an
        // S1ChordAngle then you can specify max_distance.Successor().
        //
        // DEFAULT: Distance.Infinity
        public IDistance MaxDistance { get; set; } = Infinity;

        // Specifies that points up to max_error() further away than the true
        // closest points may be substituted in the result set, as long as such
        // points satisfy all the remaining search criteria (such as max_distance).
        // This option only has an effect if max_results() is also specified;
        // otherwise all points closer than max_distance() will always be returned.
        //
        // Note that this does not affect how the distance between points is
        // computed; it simply gives the algorithm permission to stop the search
        // early as soon as the best possible improvement drops below max_error().
        //
        // This can be used to implement distance predicates efficiently.  For
        // example, to determine whether the minimum distance is less than D, the
        // IsDistanceLess() method sets max_results() == 1 and max_distance() ==
        // max_error() == D.  This causes the algorithm to terminate as soon as it
        // finds any point whose distance is less than D, rather than continuing to
        // search for a point that is even closer.
        //
        // DEFAULT: Distance.Delta.Zero
        public Delta MaxError { get; set; } = Delta.Zero;

        // Specifies that points must be contained by the given S2Region.  "region"
        // is owned by the caller and must persist during the lifetime of this
        // object.  The value may be changed between calls to FindClosestPoints(),
        // or reset by calling set_region(null).
        //
        // Note that if you want to set the region to a disc around a target point,
        // it is faster to use a PointTarget with set_max_distance() instead.  You
        // can also call both methods, e.g. to set a maximum distance and also
        // require that points lie within a given rectangle.
        public IS2Region? Region { get; set; } = null;

        // Specifies that distances should be computed by examining every point
        // rather than using the S2ShapeIndex.  This is useful for testing,
        // benchmarking, and debugging.
        //
        // DEFAULT: false
        public bool UseBruteForce { get; set; } = false;
    }

    // Each "Result" object represents a closest point.
    public readonly record struct Result(
                    // The distance from the target to this point.
                    IDistance Distance,
                    // The point itself.
                    S2Point Point, 
                    // The client-specified data associated with this point.
                    Data Data) : IComparable<Result>
    {
        #region Constructors

        // The default constructor creates an "empty" result, with a distance() of
        // Infinity() and non-dereferencable point() and data() values.
        public Result() : this(Infinity, S2Point.Empty, default) { }

        #endregion

        #region Result

        // Returns true if this Result object does not refer to any data point.
        // (The only case where an empty Result is returned is when the
        // FindClosestPoint() method does not find any points that meet the
        // specified criteria.)
        public bool IsEmpty => Point == S2Point.Empty;

        #endregion

        #region IComparable

        // Compares two Result objects first by distance, then by point_data().
        public int CompareTo(Result other)
        {
            if (Distance.CompareTo(other.Distance) != 0)
                return Distance.CompareTo(other.Distance);

            if (Point.CompareTo(other.Point) != 0)
                return Point.CompareTo(other.Point);

            return Data.CompareTo(other.Data);
        }

        public static bool operator <(Result x, Result y) => x.CompareTo(y) < 0;
        public static bool operator >(Result x, Result y) => x.CompareTo(y) > 0;
        public static bool operator <=(Result x, Result y) => x.CompareTo(y) <= 0;
        public static bool operator >=(Result x, Result y) => x.CompareTo(y) >= 0;

        #endregion
    }

    // The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
    // in increasing order of distance from the target.
    public readonly struct QueueEntry
    {
        // A lower bound on the distance from the target to "id".  This is the key
        // of the priority queue.
        public readonly IDistance Distance;

        // The cell being queued.
        public readonly S2CellId Id;

        public QueueEntry(IDistance _distance, S2CellId _id)
        {
            Distance = _distance; Id = _id;
        }

        public static bool operator <(QueueEntry a, QueueEntry b)
        {
            // The priority queue returns the largest elements first, so we want the
            // "largest" entry to have the smallest distance.
            return b.Distance.CompareTo(a.Distance) < 0;
        }
        public static bool operator >(QueueEntry a, QueueEntry b)
        {
            // The priority queue returns the largest elements first, so we want the
            // "largest" entry to have the smallest distance.
            return b.Distance.CompareTo(a.Distance) > 0;
        }
    }
}
