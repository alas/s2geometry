// An S2RegionCoverer is a class that allows arbitrary regions to be
// approximated as unions of cells (S2CellUnion).  This is useful for
// implementing various sorts of search and precomputation operations.
//
// Typical usage:
//
// S2RegionCoverer.Options options;
// options.set_max_cells(5);
// S2RegionCoverer coverer(options);
// S2Cap cap(center, radius);
// S2CellUnion covering = coverer.GetCovering(cap);
//
// This yields a vector of at most 5 cells that is guaranteed to cover the
// given cap (a disc-shaped region on the sphere).
//
// The approximation algorithm is not optimal but does a pretty good job in
// practice.  The output does not always use the maximum number of cells
// allowed, both because this would not always yield a better approximation,
// and because max_cells() is a limit on how much work is done exploring the
// possible covering as well as a limit on the final output size.
//
// Because it is an approximation algorithm, one should not rely on the
// stability of the output.  In particular, the output of the covering algorithm
// may change across different versions of the library.
//
// One can also generate interior coverings, which are sets of cells which
// are entirely contained within a region.  Interior coverings can be
// empty, even for non-empty regions, if there are no cells that satisfy
// the providedraints and are contained by the region.  Note that for
// performance reasons, it is wise to specify a max_level when computing
// interior coverings - otherwise for regions with small or zero area, the
// algorithm may spend a lot of time subdividing cells all the way to leaf
// level to try to find contained cells.

namespace S2Geometry;

public class S2RegionCoverer
{
    // Constructs an S2RegionCoverer with the given options.
    public S2RegionCoverer(Options options)
    {
        Options_ = options;

        MyDebug.Assert(options.MinLevel <= options.MaxLevel);
    }

    // Default constructor.  Options can be set using Options().
    public S2RegionCoverer() { Options_ = new Options(); }

    // Returns the current options.  Options can be modifed between calls.
    public Options Options_ { get; set; }

    // Returns an S2CellUnion that covers (GetCovering) or is contained within
    // (GetInteriorCovering) the given region and satisfies the current options.
    //
    // Note that if options().min_level() > 0 or options().level_mod() > 1, then
    // by definition the S2CellUnion may not be normalized, i.e. there may be
    // groups of four child cells that can be replaced by their parent cell.
    public S2CellUnion GetCovering(IS2Region region)
    {
        interior_covering_ = false;
        var result = GetCoveringInternal(region);
        return S2CellUnion.FromVerbatim(result);
    }

    public S2CellUnion GetInteriorCovering(IS2Region region)
    {
        interior_covering_ = true;
        var result = GetCoveringInternal(region);
        return S2CellUnion.FromVerbatim(result);
    }

    // Like the methods above, but works directly with a vector of S2CellIds.
    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    public void GetCovering(IS2Region region, out List<S2CellId> covering)
    {
        interior_covering_ = false;
        var result = GetCoveringInternal(region);
        covering = result;
    }

    public void GetInteriorCovering(IS2Region region, out List<S2CellId> interior)
    {
        interior_covering_ = true;
        var result = GetCoveringInternal(region);
        interior = result;
    }

    // Like GetCovering(), except that this method is much faster and the
    // coverings are not as tight.  All of the usual parameters are respected
    // (max_cells, min_level, max_level, and level_mod), except that the
    // implementation makes no attempt to take advantage of large values of
    // max_cells().  (A small number of cells will always be returned.)
    //
    // This function is useful as a starting point for algorithms that
    // recursively subdivide cells.
    public void GetFastCovering(IS2Region region, List<S2CellId> covering)
    {
        region.GetCellUnionBound(covering);
        CanonicalizeCovering(covering);
    }

    // Given a connected region and a starting point on the boundary or inside the
    // region, returns a set of cells at the given level that cover the region.
    // The output cells are returned in arbitrary order.
    //
    // Note that this method is *not* faster than the regular GetCovering()
    // method for most region types, such as S2Cap or S2Polygon, and in fact it
    // can be much slower when the output consists of a large number of cells.
    // Currently it can be faster at generating coverings of long narrow regions
    // such as polylines, but this may change in the future, in which case this
    // method will most likely be removed.
    public static void GetSimpleCovering(IS2Region region, S2Point start, int level, List<S2CellId> output)
    {
        FloodFill(region, new S2CellId(start).Parent(level), output);
    }

    // Like GetSimpleCovering(), but accepts a starting S2CellId rather than a
    // starting point and cell level.  Returns all edge-connected cells at the
    // same level as "start" that intersect "region", in arbitrary order.
    public static void FloodFill(IS2Region region, S2CellId start, List<S2CellId> output)
    {
        var all = new List<(S2CellId, int)>();
        var frontier = new List<S2CellId>();
        output.Clear();
        all.Add((start, start.GetHashCode()));
        frontier.Add(start);
        while (frontier.Count!=0)
        {
            S2CellId id = frontier.Last();
            frontier.RemoveAt(frontier.Count - 1);
            if (!region.MayIntersect(new S2Cell(id))) continue;
            output.Add(id);

            var neighbors = new S2CellId[4];
            id.EdgeNeighbors(neighbors);
            for (int edge = 0; edge < 4; ++edge)
            {
                var nbr = neighbors[edge];
                var hash = nbr.GetHashCode();
                all.Add((nbr, hash));
                if (hash != 0)
                {
                    frontier.Add(nbr);
                }
            }
        }
    }

    // Returns true if the given S2CellId vector represents a valid covering
    // that conforms to the current covering parameters.  In particular:
    //
    //  - All S2CellIds must be valid.
    //
    //  - S2CellIds must be sorted and non-overlapping.
    //
    //  - S2CellId levels must satisfy min_level(), max_level(), and level_mod().
    //
    //  - If covering.size() > max_cells(), there must be no two cells with
    //    a common ancestor at min_level() or higher.
    //
    //  - There must be no sequence of cells that could be replaced by an
    //    ancestor (i.e. with level_mod() == 1, the 4 child cells of a parent).
    public bool IsCanonical(S2CellUnion covering)
    {
        return IsCanonical(covering.CellIds);
    }
    public bool IsCanonical(List<S2CellId> covering)
    {
        // We check this on each call because of Options().
        MyDebug.Assert(Options_.MinLevel <= Options_.MaxLevel);

        int min_level = Options_.MinLevel;
        int max_level = Options_.TrueMaxLevel;
        int level_mod = Options_.LevelMod;
        bool too_many_cells = covering.Count > Options_.MaxCells;
        int same_parent_count = 1;
        var prev_id = S2CellId.None;
        foreach (var id in covering)
        {
            if (!id.IsValid()) return false;

            // Check that the S2CellId level is acceptable.
            int level = id.Level();
            if (level < min_level || level > max_level) return false;
            if (level_mod > 1 && (level - min_level) % level_mod != 0) return false;

            if (prev_id != S2CellId.None)
            {
                // Check that cells are sorted and non-overlapping.
                if (prev_id.RangeMax() >= id.RangeMin()) return false;

                // If there are too many cells, check that no pair of adjacent cells
                // could be replaced by an ancestor.
                if (too_many_cells && id.CommonAncestorLevel(prev_id) >= min_level)
                {
                    return false;
                }

                // Check that there are no sequences of (4 ** level_mod) cells that all
                // have the same parent (considering only multiples of "level_mod").
                int plevel = level - level_mod;
                if (plevel < min_level || level != prev_id.Level() ||
                    id.Parent(plevel) != prev_id.Parent(plevel))
                {
                    same_parent_count = 1;
                }
                else if (++same_parent_count == (1 << (2 * level_mod)))
                {
                    return false;
                }
            }
            prev_id = id;
        }
        return true;
    }

    // Modify "covering" if necessary so that it conforms to the current
    // covering parameters (max_cells, min_level, max_level, and level_mod).
    // There are no restrictions on the input S2CellIds (they may be unsorted,
    // overlapping, etc).
    public S2CellUnion CanonicalizeCovering(S2CellUnion covering)
    {
        var ids = covering.CellIds;
        CanonicalizeCovering(ids);
        return new S2CellUnion(ids);
    }
    public void CanonicalizeCovering(List<S2CellId> covering)
    {
        // We check this on each call because of Options().
        MyDebug.Assert(Options_.MinLevel <= Options_.MaxLevel);

        // Note that when the covering parameters have their default values, almost
        // all of the code in this function is skipped.

        // If any cells are too small, or don't satisfy level_mod(), then replace
        // them with ancestors.
        if (Options_.MaxLevel < S2.kMaxCellLevel || Options_.LevelMod > 1)
        {
            for (var i = 0; i < covering.Count; i++)
            {
                var id = covering[i];
                int level = id.Level();
                int new_level = AdjustLevel(Math.Min(level, Options_.MaxLevel));
                if (new_level != level)
                {
                    covering[i] = id.Parent(new_level);
                }
            }
        }

        // Sort the cells and simplify them.
        S2CellUnion.Normalize(covering);

        // Make sure that the covering satisfies min_level() and level_mod(),
        // possibly at the expense of satisfying max_cells().
        if (Options_.MinLevel > 0 || Options_.LevelMod > 1)
        {
            var tmp = new List<S2CellId>();
            S2CellUnion.Denormalize(covering, Options_.MinLevel, Options_.LevelMod, tmp);
            covering.Clear();
            covering.AddRange(tmp);
        }

        // If there are too many cells and the covering is very large, use the
        // S2RegionCoverer to compute a new covering.  (This avoids possible O(n^2)
        // behavior of the simpler algorithm below.)
        var excess = covering.Count - Options_.MaxCells;
        if (excess <= 0 || IsCanonical(covering))
        {
            return;
        }
        if (excess * covering.Count > 10000)
        {
            GetCovering(new S2CellUnion(covering), out covering);
        }
        else
        {
            // Repeatedly replace two adjacent cells in S2CellId order by their lowest
            // common ancestor until the number of cells is acceptable.
            while (covering.Count > Options_.MaxCells)
            {
                int best_index = -1, best_level = -1;
                for (int i = 0; i + 1 < covering.Count; ++i)
                {
                    int level = covering[i].CommonAncestorLevel(covering[i + 1]);
                    level = AdjustLevel(level);
                    if (level > best_level)
                    {
                        best_level = level;
                        best_index = i;
                    }
                }
                if (best_level < Options_.MinLevel) break;

                // Replace all cells contained by the new ancestor cell.
                S2CellId id = covering[best_index].Parent(best_level);
                ReplaceCellsWithAncestor(covering, id);

                // Now repeatedly check whether all children of the parent cell are
                // present, in which case we can replace those cells with their parent.
                while (best_level > Options_.MinLevel)
                {
                    best_level -= Options_.LevelMod;
                    id = id.Parent(best_level);
                    if (!ContainsAllChildren(covering, id)) break;
                    ReplaceCellsWithAncestor(covering, id);
                }
            }
        }
        MyDebug.Assert(IsCanonical(covering));
    }

    public class Candidate(S2Cell cell, int maxChildren)
    {
        public S2Cell Cell { get; } = cell;
        public bool IsTerminal { get; set; } = maxChildren == 0;
        public int NumChildren { get; set; } = 0;
        public Candidate[] Children { get; } = new Candidate[maxChildren].Fill(() => default);
    }

    // If the cell intersects the given region, return a new candidate with no
    // children, otherwise return null.  Also marks the candidate as "terminal"
    // if it should not be expanded further.
    private Candidate? NewCandidate(S2Cell cell)
    {
        if (!region_.MayIntersect(cell)) return null;

        bool is_terminal = false;
        if (cell.Level >= Options_.MinLevel)
        {
            if (interior_covering_)
            {
                if (region_.Contains(cell))
                {
                    is_terminal = true;
                }
                else if (cell.Level + Options_.LevelMod > Options_.MaxLevel)
                {
                    return null;
                }
            }
            else
            {
                if (cell.Level + Options_.LevelMod > Options_.MaxLevel ||
                    region_.Contains(cell))
                {
                    is_terminal = true;
                }
            }
        }
        //++candidates_created_counter_;
        var max_children = is_terminal ? 0 : (1 << MaxChildrenShift());
        return new Candidate(cell, max_children);
    }

    // Returns the log base 2 of the maximum number of children of a candidate.
    private int MaxChildrenShift() { return 2 * Options_.LevelMod; }

    // Frees the memory associated with a candidate.
    private static void DeleteCandidate(Candidate candidate, bool delete_children)
    {
        if (delete_children)
        {
            for (int i = 0; i < candidate.NumChildren; ++i)
                DeleteCandidate(candidate.Children[i], true);
        }
        //delete candidate;
    }

    // Processes a candidate by either adding it to the result vector or
    // expanding its children and inserting it into the priority queue.
    // Passing an argument of null does nothing.
    private void AddCandidate(Candidate candidate, List<S2CellId> result)
    {
        if (candidate is null) return;

        if (candidate.IsTerminal)
        {
            result.Add(candidate.Cell.Id);
            DeleteCandidate(candidate, true);
            return;
        }
        MyDebug.Assert(candidate.NumChildren == 0);

        // Expand one level at a time until we hit min_level() to ensure that we
        // don't skip over it.
        int num_levels = (candidate.Cell.Level < Options_.MinLevel) ?
                          1 : Options_.LevelMod;
        int num_terminals = ExpandChildren(candidate, candidate.Cell, num_levels);

        if (candidate.NumChildren == 0)
        {
            DeleteCandidate(candidate, false);

        }
        else if (!interior_covering_ &&
                 num_terminals == 1 << MaxChildrenShift() &&
                 candidate.Cell.Level >= Options_.MinLevel)
        {
            // Optimization: add the parent cell rather than all of its children.
            // We can't do this for interior coverings, since the children just
            // intersect the region, but may not be contained by it - we need to
            // subdivide them further.
            candidate.IsTerminal = true;
            AddCandidate(candidate, result);

        }
        else
        {
            // We negate the priority so that smaller absolute priorities are returned
            // first.  The heuristic is designed to refine the largest cells first,
            // since those are where we have the largest potential gain.  Among cells
            // of the same size, we prefer the cells with the fewest children.
            // Finally, among cells with equal numbers of children we prefer those
            // with the smallest number of children that cannot be refined further.
            int priority = -((((candidate.Cell.Level << MaxChildrenShift())
                               + candidate.NumChildren) << MaxChildrenShift())
                             + num_terminals);
            pq_.Add(new KeyData<int, Candidate>(priority, candidate));
        }
    }

    // Populates the children of "candidate" by expanding the given number of
    // levels from the given cell.  Returns the number of children that were
    // marked "terminal".
    private int ExpandChildren(Candidate candidate, S2Cell cell, int num_levels)
    {
        num_levels--;
        var child_cells = new S2Cell[4];
        cell.Subdivide(child_cells);
        int num_terminals = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (num_levels > 0)
            {
                if (region_.MayIntersect(child_cells[i]))
                {
                    num_terminals += ExpandChildren(candidate, child_cells[i], num_levels);
                }
                continue;
            }
            var child = NewCandidate(child_cells[i]);
            if (child is not null)
            {
                candidate.Children[candidate.NumChildren++] = child;
                if (child.IsTerminal) ++num_terminals;
            }
        }
        return num_terminals;
    }

    // Computes a set of initial candidates that cover the given region.
    private void GetInitialCandidates(List<S2CellId> result)
    {
        // Optimization: start with a small (usually 4 cell) covering of the
        // region's bounding cap.
        var tmp_coverer = new S2RegionCoverer();
        tmp_coverer.Options_.MaxCells = Math.Min(4, Options_.MaxCells);
        tmp_coverer.Options_.MaxLevel = Options_.MaxLevel;
        var cells = new List<S2CellId>();
        tmp_coverer.GetFastCovering(region_, cells);
        AdjustCellLevels(cells);
        foreach (var cell_id in cells)
        {
            AddCandidate(NewCandidate(new S2Cell(cell_id)), result);
        }
    }

    // Generates a covering.
    private List<S2CellId> GetCoveringInternal(IS2Region region)
    {
        // We check this on each call because of Options().
        MyDebug.Assert(Options_.MinLevel <= Options_.MaxLevel);

        // Strategy: Start with the 6 faces of the cube.  Discard any
        // that do not intersect the shape.  Then repeatedly choose the
        // largest cell that intersects the shape and subdivide it.
        //
        // result_ contains the cells that will be part of the output, while pq_
        // contains cells that we may still subdivide further.  Cells that are
        // entirely contained within the region are immediately added to the output,
        // while cells that do not intersect the region are immediately discarded.
        // Therefore pq_ only contains cells that partially intersect the region.
        // Candidates are prioritized first according to cell size (larger cells
        // first), then by the number of intersecting children they have (fewest
        // children first), and then by the number of fully contained children
        // (fewest children first).

        MyDebug.Assert(pq_.Count==0);
        var result = new List<S2CellId>();
        region_ = region;
        //candidates_created_counter_ = 0;

        GetInitialCandidates(result);
        while (pq_.Count!=0 && (!interior_covering_ || result.Count < Options_.MaxCells))
        {
            var top = pq_.First();
            var candidate = top.Item2;
            pq_.Remove(top);
            // For interior coverings we keep subdividing no matter how many children
            // the candidate has.  If we reach max_cells() before expanding all
            // children, we will just use some of them.  For exterior coverings we
            // cannot do this, because the result has to cover the whole region, so
            // all children have to be used.  The (candidate.num_children == 1) case
            // takes care of the situation when we already have more than max_cells()
            // in results (min_level is too high).  Subdividing the candidate with one
            // child does no harm in this case.
            if (interior_covering_ ||
                candidate.Cell.Level < Options_.MinLevel ||
                candidate.NumChildren == 1 ||
                (result.Count + pq_.Count + candidate.NumChildren <=
                 Options_.MaxCells))
            {
                // Expand this candidate into its children.
                for (int i = 0; i < candidate.NumChildren; ++i)
                {
                    if (interior_covering_ && result.Count >= Options_.MaxCells)
                    {
                        DeleteCandidate(candidate.Children[i], true);
                    }
                    else
                    {
                        AddCandidate(candidate.Children[i], result);
                    }
                }
                DeleteCandidate(candidate, false);
            }
            else
            {
                candidate.IsTerminal = true;
                AddCandidate(candidate, result);
            }
        }
        while (pq_.Count!=0)
        {
            var top = pq_.First();
            DeleteCandidate(top.Item2, true);
            pq_.Remove(top);
        }
        region_ = null;

        // Rather than just returning the raw list of cell ids, we construct a cell
        // union and then denormalize it.  This has the effect of replacing four
        // child cells with their parent whenever this does not violate the covering
        // parameters specified (min_level, level_mod, etc).  This significantly
        // reduces the number of cells returned in many cases, and it is cheap
        // compared to computing the covering in the first place.
        S2CellUnion.Normalize(result);
        if (Options_.MinLevel > 0 || Options_.LevelMod > 1)
        {
            var result_copy = result.ToList();
            S2CellUnion.Denormalize(result_copy, Options_.MinLevel, Options_.LevelMod, result);
        }
        MyDebug.Assert(IsCanonical(result));
        return result;
    }

    // If level > min_level(), then reduces "level" if necessary so that it also
    // satisfies level_mod().  Levels smaller than min_level() are not affected
    // (since cells at these levels are eventually expanded).
    private int AdjustLevel(int level)
    {
        if (Options_.LevelMod > 1 && level > Options_.MinLevel)
        {
            level -= (level - Options_.MinLevel) % Options_.LevelMod;
        }
        return level;
    }

    // Ensures that all cells with level > min_level() also satisfy level_mod(),
    // by replacing them with an ancestor if necessary.  Cell levels smaller
    // than min_level() are not modified (see AdjustLevel).  The output is
    // then normalized to ensure that no redundant cells are present.
    private void AdjustCellLevels(List<S2CellId> cells)
    {
        MyDebug.Assert(cells.IsSorted());
        if (Options_.LevelMod == 1) return;

        int output = 0;
        foreach (var id in cells)
        {
            int level = id.Level();
            int new_level = AdjustLevel(level);
            var id2 = id;
            if (new_level != level) id2 = id.Parent(new_level);
            if (output > 0 && cells[output - 1].Contains(id2)) continue;
            while (output > 0 && id2.Contains(cells[output - 1])) --output;
            cells[output++] = id2;
        }
        cells.Resize(output, default(S2CellId));
    }

    // Returns true if "covering" contains all children of "id" at level
    // (id.Level + options_.level_mod()).
    private bool ContainsAllChildren(List<S2CellId> covering, S2CellId id)
    {
        var it = covering.GetLowerBound(id.RangeMin());
        int level = id.Level() + Options_.LevelMod;
        var limit = id.ChildEnd(level);
        for (var child = id.ChildBegin(level); child != limit; child = child.Next())
        {
            if (it == covering.Count || covering[it] != child) return false;
            it++;
        }
        return true;
    }

    // Replaces all descendants of "id" in "covering" with "id".
    // REQUIRES: "covering" contains at least one descendant of "id".
    private static void ReplaceCellsWithAncestor(List<S2CellId> covering, S2CellId id)
    {
        var begin = covering.GetLowerBound(id.RangeMin());
        var end = covering.GetUpperBound(id.RangeMax());
        MyDebug.Assert(begin != end);
        covering.RemoveRange(begin + 1, end);
        covering[begin] = id;
    }

    // We save a temporary copy of the pointer passed to GetCovering() in order
    // to avoid passing this parameter around internally.  It is only used (and
    // only valid) for the duration of a single GetCovering() call.
    private IS2Region? region_ = null;

    public class Options
    {
        public int MaxCells { get; set; } = kDefaultMaxCells;

        // Sets the minimum and maximum cell levels to be used.  The default is to
        // use all cell levels.
        //
        // To find the cell level corresponding to a given physical distance, use
        // the S2Cell metrics defined in s2metrics.h.  For example, to find the
        // cell level that corresponds to an average edge length of 10km, use:
        //
        //   int level =
        //       S2.kAvgEdge.GetClosestLevel(S2Earth.KmToRadians(length_km));
        //
        // Note that min_level() takes priority over max_cells(), i.e. cells below
        // the given level will never be used even if this causes a large number
        // of cells to be returned.  (This doesn't apply to interior coverings,
        // since interior coverings make no completeness guarantees -- the result
        // is simply a set of cells that covers as much of the interior as
        // possible while satisfying the given restrictions.)
        //
        // REQUIRES: min_level() <= max_level()
        // DEFAULT: 0
        public int MinLevel
        {
            get => _min_level;
            set
            {
                MyDebug.Assert(value >= 0);
                MyDebug.Assert(value <= S2.kMaxCellLevel);
                // min_level() <= max_level() is checked by S2RegionCoverer.
                _min_level = Math.Max(0, Math.Min(S2.kMaxCellLevel, value));
            }
        }
        private int _min_level = 0;

        // REQUIRES: min_level() <= max_level()
        // DEFAULT: S2Constants.kMaxCellLevel
        public int MaxLevel
        {
            get => max_level_;
            set
            {
                MyDebug.Assert(value >= 0);
                MyDebug.Assert(value <= S2.kMaxCellLevel);
                // min_level() <= max_level() is checked by S2RegionCoverer.
                max_level_ = Math.Max(0, Math.Min(S2.kMaxCellLevel, value));
            }
        }
        private int max_level_ = S2.kMaxCellLevel;

        // Convenience function that sets both the maximum and minimum cell levels.
        public int FixedLevel
        {
            set
            {
                MinLevel = value;
                MaxLevel = value;
            }
        }

        // If specified, then only cells where (level - min_level) is a multiple
        // of "level_mod" will be used (default 1).  This effectively allows the
        // branching factor of the S2CellId hierarchy to be increased.  Currently
        // the only parameter values allowed are 1, 2, or 3, corresponding to
        // branching factors of 4, 16, and 64 respectively.
        //
        // DEFAULT: 1
        public int LevelMod
        {
            get => level_mod_;
            set
            {
                MyDebug.Assert(value >= 1);
                MyDebug.Assert(value <= 3);
                level_mod_ = Math.Max(1, Math.Min(3, value));
            }
        }
        private int level_mod_ = 1;

        // Convenience function that returns the maximum level such that
        //
        //   (level <= max_level()) && (level - min_level()) % level_mod() == 0.
        //
        // This is the maximum level that will actually be used in coverings.
        public int TrueMaxLevel
        {
            get
            {
                if (level_mod_ == 1) return max_level_;
                return max_level_ - (max_level_ - MinLevel) % level_mod_;
            }
        }

        // Sets the desired maximum number of cells in the approximation.  Note
        // the following:
        //
        //  - For any setting of max_cells(), up to 6 cells may be returned if
        //    that is the minimum number required (e.g. if the region intersects
        //    all six cube faces).  Even for very tiny regions, up to 3 cells may
        //    be returned if they happen to be located at the intersection of
        //    three cube faces.
        //
        //  - min_level() takes priority over max_cells(), i.e. cells below the
        //    given level will never be used even if this causes a large number of
        //    cells to be returned.
        //
        //  - If max_cells() is less than 4, the area of the covering may be
        //    arbitrarily large compared to the area of the original region even
        //    if the region is convex (e.g. an S2Cap or S2LatLngRect).
        //
        // Accuracy is measured by dividing the area of the covering by the area
        // of the original region.  The following table shows the median and worst
        // case values for this area ratio on a test case consisting of 100,000
        // spherical caps of random size (generated using s2region_coverer_test):
        //
        //   max_cells:        3      4     5     6     8    12    20   100   1000
        //   median ratio:  5.33   3.32  2.73  2.34  1.98  1.66  1.42  1.11   1.01
        //   worst case:  215518  14.41  9.72  5.26  3.91  2.75  1.92  1.20   1.02
        //
        // The default value of 8 gives a reasonable tradeoff between the number
        // of cells used and the accuracy of the approximation.
        //
        // DEFAULT: kDefaultMaxCells
        public const int kDefaultMaxCells = 8;
    }

    // We keep the candidates in a priority queue.  We specify a vector to hold
    // the queue entries since for some reason priority_queue<> uses a deque by
    // default.  We define our own own comparison function on QueueEntries in
    // order to make the results deterministic.  (Using the default
    // less<QueueEntry>, entries of equal priority would be sorted according to
    // the memory address of the candidate.)

    private readonly SortedSet<KeyData<int, Candidate>> pq_ = [];

    // True if we're computing an interior covering.
    private bool interior_covering_;
}
