// The Target class represents the geometry to which the distance is
// measured.  For example, there can be subtypes for measuring the distance
// to a point, an edge, or to an S2ShapeIndex (an arbitrary collection of
// geometry).
//
// Implementations do *not* need to be thread-safe.  They may cache data or
// allocate temporary data structures in order to improve performance.
// using Target = S2DistanceTarget<Distance>;

// S2ClosestEdgeQueryBase is a templatized class for finding the closest
// edge(s) between two geometries.  It is not intended to be used directly,
// but rather to serve as the implementation of various specialized classes
// with more convenient APIs (such as S2ClosestEdgeQuery).  It is flexible
// enough so that it can be adapted to compute maximum distances and even
// potentially Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between two geometries A and B.
//  - Find all edges of geometry A that are within a distance D of geometry B.
//  - Find the k edges of geometry A that are closest to a given point P.
//
// You can also specify whether polygons should include their interiors (i.e.,
// if a point is contained by a polygon, should the distance be zero or should
// it be measured to the polygon boundary?)
//
// The input geometries may consist of any number of points, polylines, and
// polygons (collectively referred to as "shapes").  Shapes do not need to be
// disjoint; they may overlap or intersect arbitrarily.  The implementation is
// designed to be fast for both simple and complex geometries.
//
// The Distance template argument is used to represent distances.  Usually it
// is a thin wrapper around S1ChordAngle, but another distance type may be
// used as long as it implements the Distance concept described in
// s2distance_target.h.  For example this can be used to measure maximum
// distances, to get more accuracy, or to measure non-spheroidal distances.

using System;

namespace S2Geometry;

public class S2ClosestEdgeQueryBase<Distance> where Distance : IEquatable<Distance>, IComparable<Distance>, IDistance<Distance>
{
    // Returns a reference to the underlying S2ShapeIndex.
    public S2ShapeIndex Index { get; private set; }
    public Options Options_ { get; private set; }

    private S2ShapeIndex.Enumerator iter_;

    private S2DistanceTarget<Distance> target_;

    // True if max_error() must be subtracted from priority queue cell distances
    // in order to ensure that such distances are measured conservatively.  This
    // is true only if the target takes advantage of max_error() in order to
    // return faster results, and 0 < max_error() < distance_limit_.
    private bool use_conservative_cell_distance_;

    // For the optimized algorihm we precompute the top-level S2CellIds that
    // will be added to the priority queue.  There can be at most 6 of these
    // cells.  Essentially this is just a covering of the indexed edges, except
    // that we also store pointers to the corresponding S2ShapeIndexCells to
    // reduce the number of index seeks required.
    //
    // The covering needs to be stored in a vector so that we can use
    // S2CellUnion.GetIntersection().
    private readonly List<S2CellId> index_covering_ = [];
    private readonly Array6<S2ShapeIndexCell?> index_cells_ = new(() => new S2ShapeIndexCell());

    // The decision about whether to use the brute force algorithm is based on
    // counting the total number of edges in the index.  However if the index
    // contains a large number of shapes, this in itself might take too long.
    // So instead we only count edges up to (max_brute_force_index_size() + 1)
    // for the current target type (stored as index_num_edges_limit_).
    private int index_num_edges_;
    private int index_num_edges_limit_;

    // The distance beyond which we can safely ignore further candidate edges.
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
    //  - Otherwise results are kept in a btree_set so that we can progressively
    //    reduce the distance limit once max_results() results have been found.
    //    (A priority queue is not sufficient because we need to be able to
    //    check whether a candidate edge is already in the result set.)
    //
    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.
    private Result result_singleton_;
    private readonly List<Result> result_vector_ = [];
    private readonly SortedSet<Result> result_set_ = []; // absl.btree_set

    // When the result edges are stored in a btree_set (see above), usually
    // duplicates can be removed simply by inserting candidate edges in the
    // current set.  However this is not true if Options.max_error() > 0 and
    // the Target subtype takes advantage of this by returning suboptimal
    // distances.  This is because when UpdateMinDistance() is called with
    // different "min_dist" parameters (i.e., the distance to beat), the
    // implementation may return a different distance for the same edge.  Since
    // the btree_set is keyed by (distance, shape_id, edge_id) this can create
    // duplicate edges in the results.
    //
    // The flag below is true when duplicates must be avoided explicitly.  This
    // is achieved by maintaining a separate set keyed by (shape_id, edge_id)
    // only, and checking whether each edge is in that set before computing the
    // distance to it.
    //
    // TODO(ericv): Check whether it is faster to avoid duplicates by default
    // (even when Options.max_results() == 1), rather than just when we need to.
    private bool avoid_duplicates_;
    private readonly Dictionary<int, int> tested_edges_ = []; // gtl.dense_hash_set

    private readonly SortedSet<QueueEntry> queue_ = [];

    // Temporaries, defined here to avoid multiple allocations / initializations.

    private readonly List<S2CellId> max_distance_covering_ = [];
    private readonly List<S2CellId> initial_cells_ = [];

    #region Constructors

    // Convenience constructor.
    //
    // Initializes the query.
    //
    // REQUIRES: ReInit() must be called if "index" is modified.
    public S2ClosestEdgeQueryBase(S2ShapeIndex index)
    {
        Index = index;
        ReInit();
    }

    #endregion

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    public void ReInit()
    {
        index_num_edges_ = 0;
        index_num_edges_limit_ = 0;
        index_covering_.Clear();
        index_cells_.Clear();
    }

    // Returns the closest edges to the given target that satisfy the given
    // options.  This method may be called multiple times.
    //
    // Note that if options().include_interiors() is true, the result vector may
    // include some entries with edge_id == -1.  This indicates that the target
    // intersects the indexed polygon with the given shape_id.
    public List<Result> FindClosestEdges(S2DistanceTarget<Distance> target, Options options)
    {
        var results = new List<Result>();
        FindClosestEdges(target, options, results);
        return results;
    }

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    public void FindClosestEdges(S2DistanceTarget<Distance> target, Options options, List<Result> results)
    {
        FindClosestEdgesInternal(target, options);
        results.Clear();
        if (options.MaxResults == 1)
        {
            if (result_singleton_.ShapeId >= 0)
            {
                results.Add(result_singleton_);
            }
        }
        else if (options.MaxResults == Options.kMaxMaxResults)
        {
            var st = new SortedSet<Result>(result_vector_);
            results.AddRange(st);
            result_vector_.Clear();
        }
        else
        {
            results.AddRange(result_set_);
            result_set_.Clear();
        }
    }

    // Convenience method that returns exactly one edge.  If no edges satisfy
    // the given search criteria, then a Result with distance == Infinity() and
    // shape_id == edge_id == -1 is returned.
    //
    // Note that if options.include_interiors() is true, edge_id == -1 is also
    // used to indicate that the target intersects an indexed polygon (but in
    // that case distance == Zero() and shape_id >= 0).
    //
    // REQUIRES: options.max_results() == 1
    public Result FindClosestEdge(S2DistanceTarget<Distance> target, Options options)
    {
        MyDebug.Assert(options.MaxResults == 1);
        FindClosestEdgesInternal(target, options);
        return result_singleton_;
    }

    private void FindClosestEdgesInternal(S2DistanceTarget<Distance> target, Options options)
    {
        target_ = target;
        Options_ = options;

        tested_edges_.Clear();
        distance_limit_ = options.MaxDistance;
        result_singleton_ = new Result();
        MyDebug.Assert(result_vector_.Count==0);
        MyDebug.Assert(result_set_.Count==0);
        MyDebug.Assert(target.MaxBruteForceIndexSize >= 0);
        if (Equals(distance_limit_, Distance.Zero)) return;

        if (options.IncludeInteriors)
        {
            var shape_ids = new SortedSet<Int32>(); // absl.btree_set
            target.VisitContainingShapes(
                Index, (S2Shape containing_shape, S2Point target_point) =>
                {
                    shape_ids.Add(containing_shape.Id);
                    return shape_ids.Count < options.MaxResults;
                });
            foreach (int shape_id in shape_ids)
            {
                AddResult(new Result(Distance.Zero, shape_id, -1));
            }
            if (Equals(distance_limit_, Distance.Zero)) return;
        }

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
        if (options.MaxError != S1ChordAngle.Zero)
        {
            target_.MaxError = options.MaxError;
            target_uses_max_error = true;
        }

        // Note that we can't compare max_error() and distance_limit_ directly
        // because one is a Delta and one is a Distance.  Instead we subtract them.
        use_conservative_cell_distance_ = target_uses_max_error &&
            (Equals(distance_limit_, Distance.Infinity) ||
                Distance.Zero < (distance_limit_ - options.MaxError));

        // Use the brute force algorithm if the index is small enough.  To avoid
        // spending too much time counting edges when there are many shapes, we stop
        // counting once there are too many edges.  We may need to recount the edges
        // if we later see a target with a larger brute force edge threshold.
        int min_optimized_edges = target_.MaxBruteForceIndexSize + 1;
        if (min_optimized_edges > index_num_edges_limit_ &&
            index_num_edges_ >= index_num_edges_limit_)
        {
            index_num_edges_ = Index.GetCountEdgesUpTo(min_optimized_edges);
            index_num_edges_limit_ = min_optimized_edges;
        }

        if (options.UseBruteForce || index_num_edges_ < min_optimized_edges)
        {
            // The brute force algorithm considers each edge exactly once.
            avoid_duplicates_ = false;
            FindClosestEdgesBruteForce();
        }
        else
        {
            // If the target takes advantage of max_error() then we need to avoid
            // duplicate edges explicitly.  (Otherwise it happens automatically.)
            avoid_duplicates_ = target_uses_max_error && options.MaxResults > 1;
            FindClosestEdgesOptimized();
        }
    }

    private void FindClosestEdgesBruteForce()
    {
        foreach (var shape in Index)
        {
            if (shape is null) continue;
            int num_edges = shape.NumEdges();
            for (int e = 0; e < num_edges; ++e)
            {
                MaybeAddResult(shape, e);
            }
        }
    }

    private void FindClosestEdgesOptimized()
    {
        InitQueue();
        // Repeatedly find the closest S2Cell to "target" and either split it into
        // its four children or process all of its edges.
        while (queue_.Count!=0)
        {
            // We need to copy the top entry before removing it, and we need to
            // remove it before adding any new entries to the queue.
            var entry = queue_.First();
            queue_.Remove(entry);
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            Distance distance = entry.Distance;
            if (!(distance < distance_limit_))
            {
                queue_.Clear();
                break;
            }
            // If this is already known to be an index cell, just process it.
            if (entry.IndexCell is not null)
            {
                ProcessEdges(entry);
                continue;
            }
            // Otherwise split the cell into its four children.  Before adding a
            // child back to the queue, we first check whether it is empty.  We do
            // this in two seek operations rather than four by seeking to the key
            // between children 0 and 1 and to the key between children 2 and 3.
            S2CellId id = entry.Id!.Value;
            iter_.Seek(id.Child(1).RangeMin());
            if (!iter_.Done() && iter_.Id <= id.Child(1).RangeMax())
            {
                ProcessOrEnqueue(id.Child(1));
            }
            if (iter_.MovePrevious() && iter_.Id >= id.RangeMin())
            {
                ProcessOrEnqueue(id.Child(0));
            }
            iter_.Seek(id.Child(3).RangeMin());
            if (!iter_.Done() && iter_.Id <= id.RangeMax())
            {
                ProcessOrEnqueue(id.Child(3));
            }
            if (iter_.MovePrevious() && iter_.Id >= id.Child(2).RangeMin())
            {
                ProcessOrEnqueue(id.Child(2));
            }
        }
    }

    private void InitQueue()
    {
        MyDebug.Assert(queue_.Count==0);
        if (index_covering_.Count==0)
        {
            // We delay iterator initialization until now to make queries on very
            // small indexes a bit faster (i.e., where brute force is used).
            iter_ = new(Index, S2ShapeIndex.InitialPosition.UNPOSITIONED);
        }

        // Optimization: if the user is searching for just the closest edge, and the
        // center of the target's bounding cap happens to intersect an index cell,
        // then we try to limit the search region to a small disc by first
        // processing the edges in that cell.  This sets distance_limit_ based on
        // the closest edge in that cell, which we can then use to limit the search
        // area.  This means that the cell containing "target" will be processed
        // twice, but in general this is still faster.
        //
        // TODO(ericv): Even if the cap center is not contained, we could still
        // process one or both of the adjacent index cells in S2CellId order,
        // provided that those cells are closer than distance_limit_.
        S2Cap cap = target_.GetCapBound();
        if (cap.IsEmpty()) return;  // Empty target.
        if (Options_.MaxResults == 1 && iter_.Locate(cap.Center))
        {
            ProcessEdges(new QueueEntry(Distance.Zero, iter_.Id, iter_.Cell));
            // Skip the rest of the algorithm if we found an intersecting edge.
            if (Equals(distance_limit_, Distance.Zero)) return;
        }
        if (index_covering_.Count==0) InitCovering();
        if (Equals(distance_limit_, Distance.Infinity))
        {
            // Start with the precomputed index covering.
            for (int i = 0; i < index_covering_.Count; ++i)
            {
                ProcessOrEnqueue(index_covering_[i], index_cells_[i]);
            }
        }
        else
        {
            // Compute a covering of the search disc and intersect it with the
            // precomputed index covering.
            var coverer = new S2RegionCoverer();
            coverer.Options_.MaxCells = 4;
            var radius = cap.Radius + distance_limit_.GetChordAngleBound();
            var search_cap = new S2Cap(cap.Center, radius);
            coverer.GetFastCovering(search_cap, max_distance_covering_);
            S2CellUnion.GetIntersection(index_covering_, max_distance_covering_, initial_cells_);

            // Now we need to clean up the initial cells to ensure that they all
            // contain at least one cell of the S2ShapeIndex.  (Some may not intersect
            // the index at all, while other may be descendants of an index cell.)
            for (int i = 0, j = 0; i < initial_cells_.Count;)
            {
                S2CellId id_i = initial_cells_[i];
                // Find the top-level cell that contains this initial cell.
                while (index_covering_[j].RangeMax() < id_i) ++j;
                S2CellId id_j = index_covering_[j];
                if (id_i == id_j)
                {
                    // This initial cell is one of the top-level cells.  Use the
                    // precomputed S2ShapeIndexCell pointer to avoid an index seek.
                    ProcessOrEnqueue(id_j, index_cells_[j]);
                    ++i; ++j;
                }
                else
                {
                    // This initial cell is a proper descendant of a top-level cell.
                    // Check how it is related to the cells of the S2ShapeIndex.
                    var r = iter_.Locate(id_i);
                    if (r == S2CellRelation.INDEXED)
                    {
                        // This cell is a descendant of an index cell.  Enqueue it and skip
                        // any other initial cells that are also descendants of this cell.
                        ProcessOrEnqueue(iter_.Id, iter_.Cell);
                        S2CellId last_id = iter_.Id.RangeMax();
                        while (++i < initial_cells_.Count && initial_cells_[i] <= last_id)
                            continue;
                    }
                    else
                    {
                        // Enqueue the cell only if it contains at least one index cell.
                        if (r == S2CellRelation.SUBDIVIDED) ProcessOrEnqueue(id_i, null);
                        ++i;
                    }
                }
            }
        }
    }

    private void InitCovering()
    {
        // Find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  These are
        // the "top-level" cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one top-level cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a top-level cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the top-level cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.

        // Don't need to reserve index_cells_ since it is an InlinedVector (ArrayN).
        index_covering_.Capacity = 6;

        // TODO(ericv): Use a single iterator (iter_) below and save position
        // information using (S2CellId, S2ShapeIndexCell) type.
        S2ShapeIndex.Enumerator next= new(Index, S2ShapeIndex.InitialPosition.BEGIN);
        S2ShapeIndex.Enumerator last= new(Index, S2ShapeIndex.InitialPosition.END);
        last.MovePrevious();
        if (next.Id != last.Id)
        {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            int level = next.Id.CommonAncestorLevel(last.Id) + 1;

            // Visit each potential top-level cell except the last (handled below).
            var last_id = last.Id.Parent(level);
            for (var id = next.Id.Parent(level); id != last_id; id = id.Next())
            {
                // Skip any top-level cells that don't contain any index cells.
                if (id.RangeMax() < next.Id) continue;

                // Find the range of index cells contained by this top-level cell and
                // then shrink the cell if necessary so that it just covers them.
                var cell_first = next;
                next.Seek(id.RangeMax().Next());
                var cell_last = next;
                AddInitialRange(cell_first, cell_last);
            }
        }
        AddInitialRange(next, last);
    }

    // Add an entry to index_covering_ and index_cells_ that covers the given
    // inclusive range of cells.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private void AddInitialRange(S2ShapeIndex.Enumerator first, S2ShapeIndex.Enumerator last)
    {
        if (first.Id == last.Id)
        {
            // The range consists of a single index cell.
            index_covering_.Add(first.Id);
            index_cells_.Add(first.Cell);
        }
        else
        {
            // Add the lowest common ancestor of the given range.
            int level = first.Id.CommonAncestorLevel(last.Id);
            MyDebug.Assert(level >= 0);
            index_covering_.Add(first.Id.Parent(level));
            index_cells_.Add(null);
        }
    }

    private void MaybeAddResult(S2Shape shape, int edge_id)
    {
        var hash = (shape.Id, edge_id).GetHashCode();
        if (avoid_duplicates_ && tested_edges_.ContainsKey(hash)) return;

        tested_edges_.Add(hash, 1);
        var edge = shape.GetEdge(edge_id);
        Distance distance = distance_limit_;
        if (target_.UpdateMinDistance(edge.V0, edge.V1, ref distance))
        {
            AddResult(new Result(distance, shape.Id, edge_id));
        }
    }

    private void AddResult(Result result)
    {
        if (Options_.MaxResults == 1)
        {
            // Optimization for the common case where only the closest edge is wanted.
            result_singleton_ = result;
            distance_limit_ = result.Distance - Options_.MaxError;
        }
        else if (Options_.MaxResults == Options.kMaxMaxResults)
        {
            result_vector_.Add(result);  // Sort/unique at end.
        }
        else
        {
            // Add this edge to result_set_.  Note that even if we already have enough
            // edges, we can't erase an element before insertion because the "new"
            // edge might in fact be a duplicate.
            result_set_.Add(result);
            int size = result_set_.Count;
            if (size >= Options_.MaxResults)
            {
                if (size > Options_.MaxResults)
                {
                    result_set_.Remove(result_set_.Last());
                }
                distance_limit_ = (Distance)result_set_.Last().Distance - Options_.MaxError;
            }
        }
    }

    // Process all the edges of the given index cell.
    private void ProcessEdges(QueueEntry entry)
    {
        var index_cell = entry.IndexCell;
        for (int s = 0; s < index_cell.NumClipped(); ++s)
        {
            var clipped = index_cell.Clipped(s);
            var shape = Index.Shape(clipped.ShapeId);
            for (int j = 0; j < clipped.NumEdges; ++j)
            {
                MaybeAddResult(shape!, clipped.Edges[j]);
            }
        }
    }

    // Enqueue the given cell id.
    // REQUIRES: iter_ is positioned at a cell contained by "id".

    private void ProcessOrEnqueue(S2CellId id)
    {
        MyDebug.Assert(id.Contains(iter_.Id));
        if (iter_.Id == id)
        {
            ProcessOrEnqueue(id, iter_.Cell);
        }
        else
        {
            ProcessOrEnqueue(id, null);
        }
    }

    // Add the given cell id to the queue.  "index_cell" is the corresponding
    // S2ShapeIndexCell, or null if "id" is not an index cell.
    //
    // This version is called directly only by InitQueue().
    private void ProcessOrEnqueue(S2CellId id, S2ShapeIndexCell? index_cell)
    {
        if (index_cell is not null)
        {
            // If this index cell has only a few edges, then it is faster to check
            // them directly rather than computing the minimum distance to the S2Cell
            // and inserting it into the queue.
            const int kMinEdgesToEnqueue = 10;
            int num_edges = CountEdges(index_cell);
            if (num_edges == 0) return;
            if (num_edges < kMinEdgesToEnqueue)
            {
                // Set "distance" to zero to avoid the expense of computing it.
                ProcessEdges(new QueueEntry(Distance.Zero, id, index_cell));
                return;
            }
        }
        // Otherwise compute the minimum distance to any point in the cell and add
        // it to the priority queue.
        var cell = new S2Cell(id);
        Distance distance = distance_limit_;
        if (!target_.UpdateMinDistance(cell, ref distance)) return;
        if (use_conservative_cell_distance_)
        {
            // Ensure that "distance" is a lower bound on the true distance to the cell.
            distance -= Options_.MaxError;
        }
        queue_.Add(new QueueEntry(distance, id, index_cell));
    }

    // Return the number of edges in the given index cell.
    private static int CountEdges(S2ShapeIndexCell cell)
    {
        int count = 0;
        for (int s = 0; s < cell.NumClipped(); ++s)
        {
            count += cell.Clipped(s).NumEdges;
        }
        return count;
    }

    // Options that control the set of edges returned.  Note that by default
    // *all* edges are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    public class Options
    {
        // Specifies that at most "max_results" edges should be returned.
        //
        // REQUIRES: max_results >= 1
        // DEFAULT: kMaxMaxResults
        public int MaxResults
        {
            get => _maxResults;
            set
            {
                MyDebug.Assert(value >= 1);
                _maxResults = value;
            }
        }
        private int _maxResults = kMaxMaxResults;

        public const int kMaxMaxResults = int.MaxValue;

        // Specifies that only edges whose distance to the target is less than
        // "max_distance" should be returned.
        //
        // Note that edges whose distance is exactly equal to "max_distance" are
        // not returned.  In most cases this doesn't matter (since distances are
        // not computed exactly in the first place), but if such edges are needed
        // then you can retrieve them by specifying "max_distance" as the next
        // largest representable Distance.  For example, if Distance is an
        // S1ChordAngle then you can specify max_distance.Successor().
        //
        // DEFAULT: Distance.Infinity()
        public Distance MaxDistance { get; set; } = Distance.Infinity;

        // Specifies that edges up to max_error() further away than the true
        // closest edges may be substituted in the result set, as long as such
        // edges satisfy all the remaining search criteria (such as max_distance).
        // This option only has an effect if max_results() is also specified;
        // otherwise all edges closer than max_distance() will always be returned.
        //
        // Note that this does not affect how the distance between edges is
        // computed; it simply gives the algorithm permission to stop the search
        // early as soon as the best possible improvement drops below max_error().
        //
        // This can be used to implement distance predicates efficiently.  For
        // example, to determine whether the minimum distance is less than D, set
        // max_results() == 1 and max_distance() == max_error() == D.  This causes
        // the algorithm to terminate as soon as it finds any edge whose distance
        // is less than D, rather than continuing to search for an edge that is
        // even closer.
        //
        // DEFAULT: Distance.Delta.Zero()
        public S1ChordAngle MaxError { get; set; } = S1ChordAngle.Zero;

        // Specifies that polygon interiors should be included when measuring
        // distances.  In other words, polygons that contain the target should
        // have a distance of zero.  (For targets consisting of multiple connected
        // components, the distance is zero if any component is contained.)  This
        // is indicated in the results by returning a (shape_id, edge_id) pair
        // with edge_id == -1, i.e. this value denotes the polygons's interior.
        //
        // Note that for efficiency, any polygon that intersects the target may or
        // may not have an (edge_id == -1) result.  Such results are optional
        // because in that case the distance to the polygon is already zero.
        //
        // DEFAULT: true
        public bool IncludeInteriors { get; set; } = true;

        // Specifies that distances should be computed by examining every edge
        // rather than using the S2ShapeIndex.  This is useful for testing,
        // benchmarking, and debugging.
        //
        // DEFAULT: false
        public bool UseBruteForce { get; set; } = false;
    }

    // The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
    // in increasing order of distance from the target.
    private readonly record struct QueueEntry(
                    Distance Distance,
                    S2CellId? Id,
                    // If "id" belongs to the index, this field stores the corresponding
                    // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
                    // S2ShapeIndexCells and this field stores null.  The purpose of this
                    // field is to avoid an extra Seek() when the queue entry is processed.
                    S2ShapeIndexCell? IndexCell) : IComparable<QueueEntry>
    {
        #region IEquatable

        public bool Equals(QueueEntry other) => Equals(Distance, other.Distance);
        public override int GetHashCode() => Distance.GetHashCode();

        #endregion

        #region IComparable

        public int CompareTo(QueueEntry other) =>
            // The priority queue returns the largest elements first, so we want the
            // "largest" entry to have the smallest distance.
            other.Distance.CompareTo(Distance);

        public static bool operator <(QueueEntry x, QueueEntry y) => x.CompareTo(y) < 0;
        public static bool operator >(QueueEntry x, QueueEntry y) => x.CompareTo(y) > 0;
        public static bool operator <=(QueueEntry x, QueueEntry y) => x.CompareTo(y) <= 0;
        public static bool operator >=(QueueEntry x, QueueEntry y) => x.CompareTo(y) >= 0;

        #endregion
    }

    // Each "Result" object represents a closest edge.  Note the following
    // special cases:
    //
    //  - (shape_id() >= 0) && (edge_id() < 0) represents the interior of a shape.
    //    Such results may be returned when options.include_interiors() is true.
    //    Such results can be identified using the is_interior() method.
    //
    //  - (shape_id() < 0) && (edge_id() < 0) is returned by FindClosestEdge()
    //    to indicate that no edge satisfies the given query options.  Such
    //    results can be identified using IsEmpty method.
    public readonly record struct Result(
                    Distance Distance, // The distance from the target to this edge.
                    Int32 ShapeId, // Identifies an indexed shape.
                    Int32 EdgeId) : IComparable<Result> // Identifies an edge within the shape.
    {
        #region Constructors

        // The default constructor yields an empty result, with a distance() of
        // Infinity() and shape_id == edge_id == -1.
        public Result() : this(Distance.Infinity, -1, -1) { }

        #endregion

        #region Result

        // Returns true if this Result object represents the interior of a shape.
        // (Such results may be returned when options.include_interiors() is true.)
        public bool IsInterior() => ShapeId >= 0 && EdgeId < 0;

        // Returns true if this Result object indicates that no edge satisfies the
        // given query options.  (This result is only returned in one special
        // case, namely when FindClosestEdge() does not find any suitable edges.
        // It is never returned by methods that return a vector of results.)
        public bool IsEmpty() => ShapeId < 0;

        #endregion

        #region IComparable

        public int CompareTo(Result other)
        {
            var c = Distance.CompareTo(other.Distance);
            if (c != 0) return c;

            c = ShapeId.CompareTo(other.ShapeId);
            if (c != 0) return c;

            return EdgeId.CompareTo(other.EdgeId);
        }

        public static bool operator <(Result x, Result y) => x.CompareTo(y) < 0;
        public static bool operator >(Result x, Result y) => x.CompareTo(y) > 0;
        public static bool operator <=(Result x, Result y) => x.CompareTo(y) <= 0;
        public static bool operator >=(Result x, Result y) => x.CompareTo(y) >= 0;

        #endregion
    }
}
