namespace S2Geometry;

using ContentsEnumerator = S2CellIndex.ContentsEnumerator;
using NonEmptyRangeEnumerator = S2CellIndex.NonEmptyRangeEnumerator;
using LabelledCell = S2CellIndex.LabelledCell;

// S2ClosestCellQueryBase is a templatized class for finding the closest
// (cell_id, label) pairs in an S2CellIndex to a given target.  It is not
// intended to be used directly, but rather to serve as the implementation of
// various specialized classes with more convenient APIs (such as
// S2ClosestCellQuery).  It is flexible enough so that it can be adapted to
// compute maximum distances and even potentially Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between a cell collection A and a target B.
//  - Find all cells in collection A that are within a distance D of target B.
//  - Find the k cells of collection A that are closest to a given point P.
//
// The target is any class that implements the S2DistanceTarget interface.
// There are predefined targets for points, edges, S2Cells, S2CellUnions, and
// S2ShapeIndexes (arbitrary collctions of points, polylines, and polygons).
//
// The IDistance template argument is used to represent distances.  Usually it
// is a thin wrapper around S1ChordAngle, but another distance type may be
// used as long as it implements the IDistance concept described in
// s2distance_targets.h.  For example this can be used to measure maximum
// distances, to get more accuracy, or to measure non-spheroidal distances.
public class S2ClosestCellQueryBase<Distance> where Distance : IEquatable<Distance>, IComparable<Distance>, IDistance
{
    private static readonly Distance Infinity = (Distance)typeof(Distance).GetField("Infinity").GetValue(null);
    private static readonly Distance Zero = (Distance)typeof(Distance).GetField("Zero").GetValue(null);

    // Options that control the set of cells returned.  Note that by default
    // *all* cells are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    //
    // This class is also available as S2ClosestCellQueryBase<Data>.Options.
    //
    // The IDistance template argument is described below.
    public class Options
    {
        public Options() { }

        // Specifies that at most "max_results" cells should be returned.
        //
        // REQUIRES: max_results >= 1
        // DEFAULT: kMaxMaxResults
        public int MaxResults
        {
            get => maxResults;
            set
            {
                Assert.True(value >= 1);
                maxResults = value;
            }
        }
        private int maxResults = kMaxMaxResults;

        public const int kMaxMaxResults = int.MaxValue;

        // Specifies that only cells whose distance to the target is less than
        // "max_distance" should be returned.
        //
        // Note that cells whose distance is exactly equal to "max_distance" are
        // not returned.  In most cases this doesn't matter (since distances are
        // not computed exactly in the first place), but if such cells are needed
        // then you can retrieve them by specifying "max_distance" as the next
        // largest representable IDistance.  For example, if IDistance is an
        // S1ChordAngle then you can specify max_distance.Successor().
        //
        // DEFAULT: IDistance.Infinity()
        public Distance MaxDistance { get; set; } = Infinity;

        // Specifies that cells up to max_error() further away than the true
        // closest cells may be substituted in the result set, as long as such
        // cells satisfy all the remaining search criteria (such as max_distance).
        // This option only has an effect if max_results() is also specified;
        // otherwise all cells closer than max_distance() will always be returned.
        //
        // Note that this does not affect how the distance between cells is
        // computed; it simply gives the algorithm permission to stop the search
        // early as soon as the best possible improvement drops below max_error().
        //
        // This can be used to implement distance predicates efficiently.  For
        // example, to determine whether the minimum distance is less than D, the
        // IsIDistanceLess() method sets max_results() == 1 and max_distance() ==
        // max_error() == D.  This causes the algorithm to terminate as soon as it
        // finds any cell whose distance is less than D, rather than continuing to
        // search for a cell that is even closer.
        //
        // DEFAULT: S1ChordAngle.Zero
        public S1ChordAngle MaxError { get; set; } = S1ChordAngle.Zero;

        // Specifies that cells must intersect the given S2Region.  "region" is
        // owned by the caller and must persist during the lifetime of this
        // object.  The value may be changed between calls to FindClosestPoints(),
        // or reset by calling set_region(null).
        //
        // Note that if you want to set the region to a disc around a target
        // point, it is faster to use a PointTarget with set_max_distance()
        // instead.  You can also call both methods, e.g. to set a maximum
        // distance and also require that cells lie within a given rectangle.
        public IS2Region Region { get; set; } = null;

        // Specifies that distances should be computed by examining every cell
        // rather than using the S2ShapeIndex.  This is useful for testing,
        // benchmarking, and debugging.
        //
        // DEFAULT: false
        public bool UseBruteForce { get; set; } = false;
    }

    // Each "Result" object represents a closest (cell_id, label) pair.
    public readonly record struct Result(
                    // The distance from the target to this cell.
                    Distance Distance,
                    // The cell itself.
                    S2CellId CellId,
                    // The label associated with this S2CellId.
                    Int32 Label
        ) : IComparable<Result>
    {
        // The default constructor yields an empty result, with a distance() of
        // Infinity() and invalid cell_id() and label() values.
        public Result()
            : this(Infinity, S2CellId.None, -1)
        {
        }

        // Returns true if this Result object does not refer to any cell.
        // (The only case where an empty Result is returned is when the
        // FindClosestCell() method does not find any cells that meet the
        // specified criteria.)
        public bool IsEmpty() => CellId == S2CellId.None;

        // Returns true if two Result objects are identical.
        public bool Equals(Result other) =>
            Equals(Distance, other.Distance) &&
            CellId == other.CellId &&
            Label == other.Label;

        public override int GetHashCode() => HashCode.Combine(Distance, CellId, Label);

        // Compares two Result objects first by distance, then by cell_id and
        // finally by label.
        public int CompareTo(Result other) //=> this < other ? -1 : this > other ? 1 : 0;
        {
            var c = Distance.CompareTo(other.Distance);
            if (c != 0) return c;

            c = CellId.CompareTo(other.CellId);
            if (c != 0) return c;

            return Label.CompareTo(other.Label);
        }

        public static bool operator <(Result x, Result y) => x.CompareTo(y) < 0;
        public static bool operator >(Result x, Result y) => x.CompareTo(y) > 0;
        public static bool operator <=(Result x, Result y) => x.CompareTo(y) <= 0;
        public static bool operator >=(Result x, Result y) => x.CompareTo(y) >= 0;
    }

    // The minimum number of ranges that a cell must contain to enqueue it
    // rather than processing its contents immediately.
    private const int kMinRangesToEnqueue = 6;

    // Default constructor; requires Init() to be called.
    public S2ClosestCellQueryBase() { }

    // Convenience constructor that calls Init().
    public S2ClosestCellQueryBase(S2CellIndex index)
    {
        Init(index);
    }

    // Initializes the query.
    // REQUIRES: ReInit() must be called if "index" is modified.
    public void Init(S2CellIndex index)
    {
        Index = index;
        contents_it_ = new ContentsEnumerator(index);
        ReInit();
    }

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    public void ReInit()
    {
        index_covering_.Clear();
    }

    // Return a reference to the underlying S2CellIndex.
    public S2CellIndex Index { get; private set; }

    // Returns the closest (cell_id, label) pairs to the given target that
    // satisfy the given options.  This method may be called multiple times.
    public List<Result> FindClosestCells(S2DistanceTarget<Distance> target, Options options)
    {
        var results = new List<Result>();
        FindClosestCells(target, options, results);
        return results;
    }

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    public void FindClosestCells(S2DistanceTarget<Distance> target, Options options, List<Result> results)
    {
        FindClosestCellsInternal(target, options);
        results.Clear();
        if (options.MaxResults == 1)
        {
            if (!result_singleton_.IsEmpty())
            {
                results.Add(result_singleton_);
            }
        }
        else if (options.MaxResults == Options.kMaxMaxResults)
        {
            result_vector_.Sort();
            results.AddRange(result_vector_);
            result_vector_.Clear();
        }
        else
        {
            results.AddRange(result_set_);
            result_set_.Clear();
        }
    }

    // Convenience method that returns exactly one (cell_id, label) pair.  If no
    // cells satisfy the given search criteria, then a Result with
    // distance() == Infinity() and IsEmpty == true is returned.
    //
    // REQUIRES: options.max_results() == 1
    public Result FindClosestCell(S2DistanceTarget<Distance> target, Options options)
    {
        Assert.True(options.MaxResults == 1);
        FindClosestCellsInternal(target, options);
        return result_singleton_;
    }

    private Options Options_() { return options_; }
    private void FindClosestCellsInternal(S2DistanceTarget<Distance> target, Options options)
    {
        target_ = target;
        options_ = options;

        tested_cells_.Clear();
        contents_it_.Reset();
        distance_limit_ = options.MaxDistance;
        result_singleton_ = new Result();
        Assert.True(!result_vector_.Any());
        Assert.True(!result_set_.Any());
        Assert.True(target.MaxBruteForceIndexSize >= 0);
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
        // max_error() only affects the algorithm once at least max_edges() edges
        // have been found that satisfy the given distance limit.  At that point,
        // max_error() is subtracted from distance_limit_ in order to ensure that
        // any further matches are closer by at least that amount.  But when
        // max_distance() < max_error(), this reduces the distance limit to 0,
        // i.e. all remaining candidate cells and edges can safely be discarded.
        // (Note that this is how IsIDistanceLess() and friends are implemented.)
        bool target_uses_max_error = false;
        if (!Equals(options.MaxError, Zero))
        {
            target_.MaxError = options.MaxError;
            target_uses_max_error = true;
        }

        // Note that we can't compare max_error() and distance_limit_ directly
        // because one is a S1ChordAngle and one is a IDistance.  Instead we subtract them.
        use_conservative_cell_distance_ = target_uses_max_error &&
            (Equals(distance_limit_, Infinity) ||
             S1ChordAngle.Zero.IsLessThan(distance_limit_.Substract(options.MaxError)));

        // Use the brute force algorithm if the index is small enough.
        if (options.UseBruteForce ||
            Index.NumCells() <= target_.MaxBruteForceIndexSize)
        {
            avoid_duplicates_ = false;
            FindClosestCellsBruteForce();
        }
        else
        {
            // If the target takes advantage of max_error() then we need to avoid
            // duplicate edges explicitly.  (Otherwise it happens automatically.)
            avoid_duplicates_ = (target_uses_max_error && options.MaxResults > 1);
            FindClosestCellsOptimized();
        }
    }
    private void FindClosestCellsBruteForce()
    {
        foreach (var cell in Index.GetCellEnumerable())
        {
            MaybeAddResult(cell.CellId, cell.Label);
        }
    }
    private void FindClosestCellsOptimized()
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
            IDistance distance = entry.Item1;
            if (!(distance.IsLessThan(distance_limit_)))
            {
                queue_.Clear();
                break;
            }
            S2CellId child = entry.Item2.ChildBegin();
            // We already know that it has too many cells, so process its children.
            // Each child may either be processed directly or enqueued again.  The
            // loop is optimized so that we don't seek unnecessarily.
            bool seek = true;
            var range = Index.GetNERNEnum();
            for (int i = 0; i < 4; ++i, child = child.Next())
            {
                seek = ProcessOrEnqueue(child, range, seek);
            }
        }
    }
    private void InitQueue()
    {
        Assert.True(!queue_.Any());

        // Optimization: rather than starting with the entire index, see if we can
        // limit the search region to a small disc.  Then we can find a covering for
        // that disc and intersect it with the covering for the index.  This can
        // save a lot of work when the search region is small.
        S2Cap cap = target_.GetCapBound();
        if (cap.IsEmpty()) return;  // Empty target.
        if (Options_().MaxResults == 1)
        {
            // If the user is searching for just the closest cell, we can compute an
            // upper bound on search radius by seeking to the center of the target's
            // bounding cap and looking at the contents of that leaf cell range.  If
            // the range intersects any cells, then the distance is zero.  Otherwise
            // we can still look at the two neighboring ranges, and use the minimum
            // distance to any cell in those ranges as an upper bound on the search
            // radius.  These cells may wind up being processed twice, but in general
            // this is still faster.
            //
            // First check the range containing or immediately following "center".
            var range = Index.GetNERNEnum();
            var target = new S2CellId(cap.Center);
            range.Seek(target);
            AddRange(range);
            if (Equals(distance_limit_, Zero)) return;

            // If the range immediately follows "center" (rather than containing it),
            // then check the previous non-empty range as well.
            if (range.Current.StartId > target && range.MovePrevious())
            {
                AddRange(range);
                if (Equals(distance_limit_, Zero)) return;
            }
        }

        // We start with a covering of the set of indexed cells, then intersect it
        // with the maximum search radius disc (if any).
        //
        // Note that unlike S2ClosestPointQuery, we can't also intersect with the
        // given region (if any).  This is because the index cells in the result are
        // only required to intersect the region.  This means that an index cell that
        // intersects the region's covering may be much closer to the target than the
        // covering itself, which means that we cannot use the region's covering to
        // restrict the search.
        //
        // TODO(ericv): If this feature becomes important, this could be fixed by
        // (1) computing a covering of the region, (2) looking up any index cells
        // that contain each covering cell by seeking to covering_cell.RangeMin,
        // (3) replacing each covering cell by the largest such cell (if any), and
        // (4) normalizing the result.
        if (!index_covering_.Any()) InitCovering();
        var initial_cells = index_covering_;
        if (distance_limit_.IsLessThan(Infinity))
        {
            var coverer = new S2RegionCoverer();
            coverer.Options_.MaxCells = 4;
            var radius = cap.Radius + distance_limit_.GetChordAngleBound();
            var search_cap = new S2Cap(cap.Center, radius);
            coverer.GetFastCovering(search_cap, max_distance_covering_);
            S2CellUnion.GetIntersection(initial_cells, max_distance_covering_, intersection_with_max_distance_);
            initial_cells = intersection_with_max_distance_;
        }
        var range2 = Index.GetNERNEnum();
        for (int i = 0; i < initial_cells.Count; ++i)
        {
            var id = initial_cells[i];
            bool seek = (i == 0) || id.RangeMin() >= range2.GetLimitId();
            ProcessOrEnqueue(id, range2, seek);
            if (range2.Done()) break;
        }
    }
    private void InitCovering()
    {
        // Compute the "index covering", which is a small number of S2CellIds that
        // cover the indexed cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one covering cell
        // per spanned face, just big enough to cover the indexed cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the indexed cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a covering cell that is big enough to just fit
        // those indexed cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the covering cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.
        index_covering_.Capacity = 6;
        var it = Index.GetNERNEnum();
        var last = Index.GetNERNEnum();
        last.Finish();
        if (!last.MovePrevious()) return; // Empty index.

        it.MoveNext();
        var index_last_id = last.GetLimitId().Prev();
        if (it.Current.StartId != last.Current.StartId)
        {
            // The index contains at least two distinct S2CellIds (because otherwise
            // there would only be one non-empty range).  Choose a level such that the
            // entire index can be spanned with at most 6 cells (if the index spans
            // multiple faces) or 4 cells (it the index spans a single face).
            int level = it.Current.StartId.CommonAncestorLevel(index_last_id) + 1;

            // Visit each potential covering cell except the last (handled below).
            var start_id = it.Current.StartId.Parent(level);
            var last_id = index_last_id.Parent(level);
            for (var id = start_id; id != last_id; id = id.Next())
            {
                // Skip any covering cells that don't contain an indexed range.
                if (id.RangeMax() < it.Current.StartId) continue;

                // Find the indexed range contained by this covering cell and then
                // shrink the cell if necessary so that it just covers this range.
                var cell_first_id = it.Current.StartId;
                it.Seek(id.RangeMax().Next());
                // Find the last leaf cell covered by the previous non-empty range.
                last = it;
                last.MovePrevious();
                AddInitialRange(cell_first_id, last.GetLimitId().Prev());
            }
        }
        AddInitialRange(it.Current.StartId, index_last_id);
    }
    private void AddInitialRange(S2CellId first_id, S2CellId last_id)
    {
        // Add the lowest common ancestor of the given range.
        int level = first_id.CommonAncestorLevel(last_id);
        Assert.True(level >= 0);
        index_covering_.Add(first_id.Parent(level));
    }
    private void MaybeAddResult(S2CellId cell_id, Int32 label)
    {
        // TODO(ericv): Consider having this method return false when distance_limit_
        // is reduced to zero, and terminating any calling loops early.

        var hash = new LabelledCell(cell_id, label).GetHashCode();
        if (avoid_duplicates_ && tested_cells_.ContainsKey(hash))
        {
            return;
        }

        tested_cells_.Add(hash, 1);
        // TODO(ericv): It may be relatively common to add the same S2CellId
        // multiple times with different labels.  This could be optimized by
        // remembering the last "cell_id" argument and its distance.  However this
        // may not be beneficial when Options.max_results() == 1, for example.
        var cell = new S2Cell(cell_id);
        Distance distance = distance_limit_;
        if (!target_.UpdateMinDistance(cell, ref distance)) return;

        var region = Options_().Region;
        if (region != null && !region.MayIntersect(cell)) return;

        var result = new Result(distance, cell_id, label);
        if (Options_().MaxResults == 1)
        {
            // Optimization for the common case where only the closest cell is wanted.
            result_singleton_ = result;
            distance_limit_ = (Distance)result.Distance.Substract(Options_().MaxError);
        }
        else if (Options_().MaxResults == Options.kMaxMaxResults)
        {
            result_vector_.Add(result);  // Sort/unique at end.
        }
        else
        {
            // Add this cell to result_set_.  Note that even if we already have enough
            // edges, we can't erase an element before insertion because the "new"
            // edge might in fact be a duplicate.
            result_set_.Add(result);
            int size = result_set_.Count;
            if (size >= Options_().MaxResults)
            {
                if (size > Options_().MaxResults)
                {
                    result_set_.Remove(result_set_.Last());
                }
                distance_limit_ = (Distance)result_set_.Last().Distance.Substract(Options_().MaxError);
            }
        }
    }
    // Either process the contents of the given cell immediately, or add it to the
    // queue to be subdivided.  If "seek" is false, then "iter" must be positioned
    // at the first non-empty range (if any) with start_id() >= id.RangeMin.
    //
    // Returns "true" if the cell was added to the queue, and "false" if it was
    // processed immediately, in which case "iter" is positioned at the first
    // non-empty range (if any) with start_id() > id.RangeMax.
    private bool ProcessOrEnqueue(S2CellId id, NonEmptyRangeEnumerator iter, bool seek)
    {
        if (seek) iter.Seek(id.RangeMin());
        var last = id.RangeMax();
        if (iter.Current.StartId > last)
        {
            return false;  // No need to seek to next child.
        }
        // If this cell intersects at least "kMinRangesToEnqueue" leaf cell ranges
        // (including ranges whose contents are empty), then enqueue it.  We test
        // this by advancing (n - 1) ranges and checking whether that range also
        // intersects this cell.
        var max_it = (NonEmptyRangeEnumerator)iter.CustomClone();
        if (max_it.Advance(kMinRangesToEnqueue - 1) && max_it.Current.StartId <= last)
        {
            // This cell intersects at least kMinRangesToEnqueue ranges, so enqueue it.
            var cell = new S2Cell(id);
            var distance = distance_limit_;
            // We check "region_" second because it may be relatively expensive.
            if (target_.UpdateMinDistance(cell, ref distance) &&
                (Options_().Region != null || Options_().Region.MayIntersect(cell)))
            {
                if (use_conservative_cell_distance_)
                {
                    // Ensure that "distance" is a lower bound on distance to the cell.
                    distance = (Distance)distance.Substract(Options_().MaxError);
                }
                queue_.Add(new ReverseKeyData<IDistance, S2CellId>(distance, id));
            }
            return true;  // Seek to next child.
        }
        // There were few enough ranges that we might as well process them now.
        for (; iter.Current.StartId <= last; iter.MoveNext())
        {
            AddRange(iter);
        }
        return false;  // No need to seek to next child.
    }
    private void AddRange(NonEmptyRangeEnumerator it)
    {
        for (contents_it_.StartUnion(it.Current);
             contents_it_.MoveNext();)
        {
            MaybeAddResult(contents_it_.Current.CellId, contents_it_.Current.Label);
        }
    }

    private Options options_;
    private S2DistanceTarget<Distance> target_;

    // True if max_error() must be subtracted from priority queue cell distances
    // in order to ensure that such distances are measured conservatively.  This
    // is true only if the target takes advantage of max_error() in order to
    // return faster results, and 0 < max_error() < distance_limit_.
    private bool use_conservative_cell_distance_;

    // For the optimized algorithm we precompute the top-level S2CellIds that
    // will be added to the priority queue.  There can be at most 6 of these
    // cells.  Essentially this is just a covering of the indexed cells.
    private readonly List<S2CellId> index_covering_ = new();

    // The distance beyond which we can safely ignore further candidate cells.
    // (Candidates that are exactly at the limit are ignored; this is more
    // efficient for UpdateMinIDistance() and should not affect clients since
    // distance measurements have a small amount of error anyway.)
    //
    // Initially this is the same as the maximum distance specified by the user,
    // but it can also be updated by the algorithm (see MaybeAddResult).
    private Distance distance_limit_;

    // The current result set is stored in one of three ways:
    //
    //  - If max_results() == 1, the best result is kept in result_singleton_.
    //
    //  - If max_results() == kMaxMaxResults, results are appended to
    //    result_vector_ and sorted/uniqued at the end.
    //
    //  - Otherwise results are kept in a btree_set so that we can progressively
    //    reduce the distance limit once max_results() results have been found.
    //    (A priority queue is not sufficient because we need to be able to
    //    check whether a candidate cell is already in the result set.)
    //
    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.
    private Result result_singleton_;
    private readonly List<Result> result_vector_ = new();
    private readonly SortedSet<Result> result_set_ = new(); // absl.btree_set<Result>

    // When the results are stored in a btree_set (see above), usually
    // duplicates can be removed simply by inserting candidate cells in the
    // current result set.  However this is not true if Options.max_error() > 0
    // and the Target subtype takes advantage of this by returning suboptimal
    // distances.  This is because when UpdateMinIDistance() is called with
    // different "min_dist" parameters (i.e., the distance to beat), the
    // implementation may return a different distance for the same cell.  Since
    // the btree_set is keyed by (distance, cell_id, label) this can create
    // duplicate results.
    //
    // The flag below is true when duplicates must be avoided explicitly.  This
    // is achieved by maintaining a separate set keyed by (cell_id, label) only,
    // and checking whether each edge is in that set before computing the
    // distance to it.
    //
    // TODO(ericv): Check whether it is faster to avoid duplicates by default
    // (even when Options.max_results() == 1), rather than just when we need to.
    private bool avoid_duplicates_;
    private readonly Dictionary<int, int> tested_cells_ = new();

    // Priority queue of unprocessed S2CellIds, sorted
    // in increasing order of distance from the target.
    private readonly SortedSet<ReverseKeyData<IDistance, S2CellId>> queue_
               = new();

    // The Target class represents the geometry to which the distance is
    // measured.  For example, there can be subtypes for measuring the distance
    // to a point, an edge, or to an S2ShapeIndex (an arbitrary collection of
    // geometry).
    //
    // Implementations do *not* need to be thread-safe.  They may cache data or
    // allocate temporary data structures in order to improve performance.
    // using Target = S2DistanceTarget<Distance>;

    // Used to iterate over the contents of an S2CellIndex range.  It is defined
    // here to take advantage of the fact that when multiple ranges are visited
    // in increasing order, duplicates can automatically be eliminated.
    private ContentsEnumerator contents_it_;

    // Temporaries, defined here to avoid multiple allocations / initializations.
    private readonly List<S2CellId> max_distance_covering_ = new();
    private readonly List<S2CellId> intersection_with_max_distance_ = new();
}
