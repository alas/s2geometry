// An S2CellUnion is a region consisting of cells of various sizes.  Typically
// a cell union is used to approximate some other shape.  There is a tradeoff
// between the accuracy of the approximation and how many cells are used.
// Unlike polygons, cells have a fixed hierarchical structure.  This makes
// them more suitable for optimizations based on preprocessing.
//
// An S2CellUnion is represented as a vector of sorted, non-overlapping
// S2CellIds.  By default the vector is also "normalized", meaning that groups
// of 4 child cells have been replaced by their parent cell whenever possible.
// S2CellUnions are not required to be normalized, but certain operations will
// return different results if they are not (e.g., Contains(S2CellUnion).)
//
// S2CellUnion is movable and copyable.

namespace S2Geometry;

public record class S2CellUnion(List<S2CellId> CellIds) : IS2Region<S2CellUnion>, IEnumerable<S2CellId>, IDecoder<S2CellUnion>
{
    #region Fields, Constants

    // Constructs a cell union for the whole sphere.
    public static S2CellUnion WholeSphere() => new(new int[] { 0, 1, 2, 3, 4, 5 }
        .Select(t => S2CellId.FromFace(t)).ToList());

    #endregion

    #region Constructors

    // Creates an empty cell union.
    public S2CellUnion()
        : this(new List<S2CellId>()) { }

    // Constructs a cell union with the given S2CellIds, then calls Normalize()
    // to sort them, remove duplicates, and merge cells when possible.  (See
    // FromNormalized if your vector is already normalized.)
    //
    // The argument is passed by value, so if you are passing a named variable
    // and have no further use for it, consider using move().
    //
    // A cell union containing a single S2CellId may be constructed like this:
    //
    //     S2CellUnion example({cell_id});
    public S2CellUnion(List<S2CellId> cell_ids, bool normalize = true)
        : this(cell_ids)
    {
        if (normalize)
        {
            Normalize();
        }
    }

    public S2CellUnion(UInt64[] cell_ids, bool checkValidity = true)
        : this(cell_ids.Select(t => new S2CellId(t)).ToList(), checkValidity) { }

    #endregion

    #region Factories

    // Constructs a cell union from S2CellIds that have already been normalized
    // (typically because they were extracted from another S2CellUnion).
    //
    // The argument is passed by value, so if you are passing a named variable
    // and have no further use for it, consider using move().
    //
    // REQUIRES: "cell_ids" satisfies the requirements of IsNormalized().
    public static S2CellUnion FromNormalized(List<S2CellId> cell_ids)
    {
        var result = new S2CellUnion(cell_ids);
        Assert.True(result.IsNormalized());
        return result;
    }

    // Constructs a cell union from a vector of sorted, non-overlapping
    // S2CellIds.  Unlike the other constructors, FromVerbatim does not require
    // that groups of 4 child cells have been replaced by their parent cell.  In
    // other words, "cell_ids" must satisfy the requirements of IsValid but
    // not necessarily IsNormalized().
    //
    // Note that if the cell union is not normalized, certain operations may
    // return different results (e.g., Contains(S2CellUnion)).
    //
    // REQUIRES: "cell_ids" satisfies the requirements of IsValid.
    public static S2CellUnion FromVerbatim(List<S2CellId> cell_ids)
    {
        var result = new S2CellUnion(cell_ids);
#if s2debug
        Assert.True(result.IsValid());
#endif
        return result;
    }

#if s2debug
    public static S2CellUnion FromVerbatimNoCheck(List<S2CellId> cell_ids)
    {
        return new S2CellUnion(cell_ids);
    }
#endif

    // Constructs a cell union that corresponds to a continuous range of cell
    // ids.  The output is a normalized collection of cell ids that covers the
    // leaf cells between "min_id" and "max_id" inclusive.
    //
    // REQUIRES: min_id.IsLeaf, max_id.IsLeaf, min_id <= max_id.
    public static S2CellUnion FromMinMax(S2CellId min_id, S2CellId max_id)
    {
        var result = new S2CellUnion();
        result.InitFromMinMax(min_id, max_id);
        return result;
    }

    // Like FromMinMax() except that the union covers the range of leaf cells
    // from "begin" (inclusive) to "end" (exclusive), as with Python ranges or
    // STL iterator ranges.  If (begin == end) the result is empty.
    //
    // REQUIRES: begin.IsLeaf, end.IsLeaf, begin <= end.
    public static S2CellUnion FromBeginEnd(S2CellId begin, S2CellId end)
    {
        var result = new S2CellUnion();
        result.InitFromBeginEnd(begin, end);
        return result;
    }

    #endregion

    #region S2CellUnion

    public void InitFromMinMax(S2CellId min_id, S2CellId max_id)
    {
        Assert.True(max_id.IsValid());
        InitFromBeginEnd(min_id, max_id.Next());
    }
    public void InitFromBeginEnd(S2CellId begin, S2CellId end)
    {
        Assert.True(begin.IsLeaf());
        Assert.True(end.IsLeaf());

        var kLeafEnd = S2CellId.End(S2CellId.kMaxLevel);
        Assert.True(begin.IsValid() || begin == kLeafEnd);
        Assert.True(end.IsValid() || end == kLeafEnd);

        Assert.True(begin <= end);

        // We repeatedly add the largest cell we can.
        CellIds.Clear();
        for (var id = begin.MaximumTile(end); id != end; id = id.Next().MaximumTile(end))
        {
            CellIds.Add(id);
        }
        // The output is already normalized.
        Assert.True(IsNormalized());
    }

    // Clears the contents of the cell union and minimizes memory usage.
    public void Clear()
    {
        CellIds.Clear();
    }

    // Gives ownership of the vector data to the client without copying, and
    // clears the content of the cell union.  The original data in cell_ids
    // is lost if there was any.
    /*public List<S2CellId> Release()
    {
        // vector's rvalue reference constructor does not necessarily leave
        // moved-from value in empty state, so swap instead.
        var tmp = CellIds;
        CellIds = null;
        return tmp;
    }*/

    // Convenience methods for accessing the individual cell ids.
    public S2CellId CellId(int i) => CellIds[i];

    // Vector-like methods for accessing the individual cell ids.
    public int Size() => CellIds.Count;

    public bool IsEmpty() => !CellIds.Any();

    // Returns true if the cell union is valid, meaning that the S2CellIds are
    // valid, non-overlapping, and sorted in increasing order.
    public bool IsValid()
    {
        if (CellIds.Count > 0 && !CellId(0).IsValid()) return false;
        for (int i = 1; i < CellIds.Count; ++i)
        {
            if (!CellId(i).IsValid()) return false;
            if (CellId(i - 1).RangeMax() >= CellId(i).RangeMin()) return false;
        }
        return true;
    }

    // Returns true if the cell union is normalized, meaning that it is
    // satisfies IsValid and that no four cells have a common parent.
    // Certain operations such as Contains(S2CellUnion) will return a different
    // result if the cell union is not normalized.
    public bool IsNormalized()
    {
        if (CellIds.Count > 0 && !CellId(0).IsValid()) return false;
        for (int i = 1; i < CellIds.Count; ++i)
        {
            if (!CellId(i).IsValid()) return false;
            if (CellId(i - 1).RangeMax() >= CellId(i).RangeMin()) return false;
            if (i >= 3 && AreSiblings(CellId(i - 3), CellId(i - 2),
                                      CellId(i - 1), CellId(i)))
            {
                return false;
            }
        }
        return true;
    }

    // Normalizes the cell union by discarding cells that are contained by other
    // cells, replacing groups of 4 child cells by their parent cell whenever
    // possible, and sorting all the cell ids in increasing order.
    //
    // Returns true if the number of cells was reduced.
    // TODO(ericv): Change this method to return void.
    public void Normalize()
    {
        Normalize(CellIds);
    }

    // Replaces "output" with an expanded version of the cell union where any
    // cells whose level is less than "min_level" or where (level - min_level)
    // is not a multiple of "level_mod" are replaced by their children, until
    // either both of these conditions are satisfied or the maximum level is
    // reached.
    //
    // This method allows a covering generated by S2RegionCoverer using
    // min_level() or level_mod()raints to be stored as a normalized cell
    // union (which allows various geometric computations to be done) and then
    // converted back to the original list of cell ids that satisfies the
    // desiredraints.
    public void Denormalize(int min_level, int level_mod, List<S2CellId> output)
    {
        Denormalize(CellIds, min_level, level_mod, output);
    }

    // If there are more than "excess" elements of the cell_ids() vector that
    // are allocated but unused, reallocates the array to eliminate the excess
    // space.  This reduces memory usage when many cell unions need to be held
    // in memory at once.
    public void Pack()
    {
        CellIds.TrimExcess();
    }

    // Returns true if "a" lies entirely before "b" on the Hilbert curve.  Note that
    // this is not even a weak ordering, since incomparability is not transitive.
    // Nevertheless, given a sorted vector of disjoint S2CellIds it can be used be
    // used to find the first element that might intersect a given target S2CellId:
    //
    //   auto it = std::lower_bound(v.begin(), v.end(), target, EntirelyPrecedes);
    //
    // This works because std::lower_bound() only requires that the elements are
    // partitioned with respect to the given predicate (which is true as long as the
    // S2CellIds are sorted and disjoint).
    public static bool EntirelyPrecedes(S2CellId a, S2CellId b)
    {
        return a.RangeMax() < b.RangeMin();
    }

    // Returns true if the cell union contains the given cell id.  Containment
    // is defined with respect to regions, e.g. a cell contains its 4 children.
    // This is a fast operation (logarithmic in the size of the cell union).
    //
    // CAVEAT: If you have constructed a non-normalized S2CellUnion using
    // FromVerbatim, note that groups of 4 child cells are *not* considered to
    // contain their parent cell.  To get this behavior you must use one of the
    // other constructors or call Normalize() explicitly.
    public bool Contains(S2CellId id)
    {
        // This is an exact test.  Each cell occupies a linear span of the S2
        // space-filling curve, and the cell id is simply the position at the center
        // of this span.  The cell union ids are sorted in increasing order along
        // the space-filling curve.  So we simply find the first cell that might
        // intersect (is not entirely before) the target (using binary search).
        // There is containment if and only if this cell id contains the target id.
        Assert.True(id.IsValid());

        var i = CellIds.GetLowerBound(id, new EntirelyPrecedes_Comparer());
        return i != CellIds.Count && CellIds[i].Contains(id);
    }

    public class EntirelyPrecedes_Comparer : Comparer<S2CellId>
    {
        public override int Compare(S2CellId x, S2CellId y)
        {
            return EntirelyPrecedes(x, y) ? -1 : 1;
        }
    }

    // Returns true if the cell union intersects the given cell id.
    // This is a fast operation (logarithmic in the size of the cell union).
    public bool Intersects(S2CellId id)
    {
        // This is an exact test; see the comments for Contains() above.
        Assert.True(id.IsValid());

        var i = CellIds.GetLowerBound(id, new EntirelyPrecedes_Comparer());
        return i != CellIds.Count && CellIds[i].Intersects(id);
    }

    // Returns true if this cell union contains the given other cell union.
    //
    // CAVEAT: If you have constructed a non-normalized S2CellUnion using
    // FromVerbatim, note that groups of 4 child cells are *not* considered to
    // contain their parent cell.  To get this behavior you must use one of the
    // other constructors or call Normalize() explicitly.
    public bool Contains(S2CellUnion y)
    {
        if (y.IsEmpty()) return true;
        if (IsEmpty()) return false;
        var ii = 0;

        foreach (var y_id in y)
        {
            var i = CellIds[ii];
            // If our first cell ends before the one we need to contain, advance
            // where we start searching.
            if (EntirelyPrecedes(i, y_id))
            {
                ii = CellIds.GetLowerBound(y_id, ii + 1, CellIds.Count, new EntirelyPrecedes_Comparer());
                i = CellIds[ii];
                // If we're at the end, we don't contain the current y_id.
                if (ii == CellIds.Count) return false;
            }
            if (!i.Contains(y_id)) return false;
        }
        return true;
    }

    // Returns true if this cell union intersects the given other cell union.
    public bool Intersects(S2CellUnion y)
    {
        // Walk along the two sorted vectors, looking for overlap.
        for (int ii = 0, jj = 0; ii != CellIds.Count && jj != y.CellIds.Count;)
        {
            var i = CellIds[ii];
            var j = CellIds[jj];
            if (EntirelyPrecedes(i, j))
            {
                // Advance "i" to the first cell that might overlap *j.
                ii = CellIds.GetLowerBound(j, ii + 1, CellIds.Count, new EntirelyPrecedes_Comparer());
                i = CellIds[ii];
                continue;
            }
            if (EntirelyPrecedes(j, i))
            {
                // Advance "j" to the first cell that might overlap *i.
                jj = y.CellIds.GetLowerBound(i, jj + 1, y.CellIds.Count, new EntirelyPrecedes_Comparer());
                j = CellIds[jj];
                continue;
            }
            // Neither cell is to the left of the other, so they must intersect.
            Assert.True(i.Intersects(j));
            return true;
        }
        return false;
    }

    // Returns the union of the two given cell unions.
    public S2CellUnion Union(S2CellUnion y)
    {
        var cell_ids = new List<S2CellId>(CellIds.Count + y.CellIds.Count);
        cell_ids.AddRange(CellIds);
        cell_ids.AddRange(y.CellIds);
        return new S2CellUnion(cell_ids);
    }

    // Specialized version of GetIntersection() that returns the intersection of
    // a cell union with an S2CellId.  This can be useful for splitting a cell
    // union into pieces.
    public S2CellUnion Intersection(S2CellId id)
    {
        Assert.True(id.IsValid());
        S2CellUnion result = new();
        if (Contains(id))
        {
            result.CellIds.Add(id);
        }
        else
        {
            var i = CellIds.GetLowerBound(id.RangeMin());
            S2CellId id_max = id.RangeMax();
            while (i != CellIds.Count && CellIds[i] <= id_max)
                result.CellIds.Add(CellIds[i++]);
        }
        Assert.True(result.IsNormalized() || !IsNormalized());
        return result;
    }

    // Returns the intersection of the two given cell unions.
    public S2CellUnion Intersection(S2CellUnion y)
    {
        var result = new S2CellUnion();
        GetIntersection(CellIds, y.CellIds, result.CellIds);
        // The output is normalized as long as both inputs are normalized.
        Assert.True(result.IsNormalized() || !IsNormalized() || !y.IsNormalized());
        return result;
    }

    // Returns the difference of the two given cell unions.
    public S2CellUnion Difference(S2CellUnion y)
    {
        // TODO(ericv): this is approximately O(N*log(N)), but could probably
        // use similar techniques as GetIntersection() to be more efficient.
        var ret = new List<S2CellId>();
        foreach (S2CellId id in this)
        {
            GetDifferenceInternal(id, y, ret);
        }
        var result = new S2CellUnion(ret);
        // The output is normalized as long as the first argument is normalized.
        Assert.True(result.IsNormalized() || !IsNormalized());
        return result;
    }

    // Expands the cell union by adding a buffer of cells at "expand_level"
    // around the union boundary.
    //
    // For each cell "c" in the union, we add all neighboring cells at level
    // "expand_level" that are adjacent to "c".  Note that there can be many
    // such cells if "c" is large compared to "expand_level".  If "c" is smaller
    // than "expand_level", we first add the parent of "c" at "expand_level" and
    // then add all the neighbors of that cell.
    //
    // Note that the size of the output is exponential in "expand_level".  For
    // example, if expand_level == 20 and the input has a cell at level 10,
    // there will be on the order of 4000 adjacent cells in the output.  For
    // most applications the Expand(min_radius, max_level_diff) method below is
    // easier to use.
    public S2CellUnion Expand(int expand_level)
    {
        var output = new List<S2CellId>();
        var level_lsb = S2CellId.LowestOnBitForLevel(expand_level);
        for (int i = CellIds.Count; --i >= 0;)
        {
            S2CellId id = CellId(i);
            if (id.LowestOnBit() < level_lsb)
            {
                id = id.Parent(expand_level);
                // Optimization: skip over any cells contained by this one.  This is
                // especially important when very small regions are being expanded.
                while (i > 0 && id.Contains(CellId(i - 1))) --i;
            }
            output.Add(id);
            id.AppendAllNeighbors(expand_level, output);
        }
        return new S2CellUnion(output);
    }

    // Expands the cell union such that it contains all points whose distance to
    // the cell union is at most "min_radius", but do not use cells that are
    // more than "max_level_diff" levels higher than the largest cell in the
    // input.  The second parameter controls the tradeoff between accuracy and
    // output size when a large region is being expanded by a small amount
    // (e.g. expanding Canada by 1km).  For example, if max_level_diff == 4 the
    // region will always be expanded by approximately 1/16 the width of its
    // largest cell.  Note that in the worst case, the number of cells in the
    // output can be up to 4 * (1 + 2 ** max_level_diff) times larger than the
    // number of cells in the input.
    public S2CellUnion Expand(S1Angle min_radius, int max_level_diff)
    {
        var res = this;
        int min_level = S2.kMaxCellLevel;
        foreach (S2CellId id in this)
        {
            min_level = Math.Min(min_level, id.Level());
        }
        // Find the maximum level such that all cells are at least "min_radius" wide.
        int radius_level = S2.kMinWidth.GetLevelForMinValue(min_radius.Radians);
        if (radius_level == 0 && min_radius.Radians > S2.kMinWidth.GetValue(0))
        {
            // The requested expansion is greater than the width of a face cell.
            // The easiest way to handle this is to expand twice.
            res = Expand(0);
        }
        return res.Expand(Math.Min(min_level + max_level_diff, radius_level));
    }

    // The number of leaf cells covered by the union.
    // This will be no more than 6*2^60 for the whole sphere.
    public UInt64 LeafCellsCovered()
    {
        UInt64 num_leaves = 0;
        foreach (S2CellId id in this)
        {
            int inverted_level = S2.kMaxCellLevel - id.Level();
            num_leaves += (1UL << (inverted_level << 1));
        }
        return num_leaves;
    }

    // Approximates this cell union's area in steradians by summing the average
    // area of each contained cell's average area, using the AverageArea method
    // from the S2Cell class.  This is equivalent to the number of leaves covered,
    // multiplied by the average area of a leaf.  Note that AverageArea does not
    // take into account distortion of cell, and thus may be off by up to a
    // factor of up to 1.7.
    //
    // NOTE: Since this is proportional to LeafCellsCovered(), it is
    // always better to use that function if all you care about is
    // the relative average area between objects.
    public double AverageBasedArea()
    {
        return S2Cell.AverageArea(S2.kMaxCellLevel) * LeafCellsCovered();
    }

    // Calculates this cell union's area in steradians by summing the approximate
    // area for each contained cell, using the ApproxArea method from the S2Cell
    // class.
    public double ApproxArea()
    {
        double area = 0;
        foreach (S2CellId id in this)
        {
            area += new S2Cell(id).ApproxArea();
        }
        return area;
    }

    // Calculates this cell union's area in steradians by summing the exact area
    // for each contained cell, using the Exact method from the S2Cell class.
    public double ExactArea()
    {
        double area = 0;
        foreach (S2CellId id in this)
        {
            area += new S2Cell(id).ExactArea();
        }
        return area;
    }

    #region Static Versions

    ////////////////////////////////////////////////////////////////////////
    // Static methods intended for high-performance clients that prefer to
    // manage their own storage.

    // Like Normalize(), but works with a vector of S2CellIds.
    // Equivalent to:
    //   *cell_ids = S2CellUnion(move(*cell_ids)).Release();
    public static void Normalize(List<S2CellId> ids)
    {
        // Optimize the representation by discarding cells contained by other cells,
        // and looking for cases where all subcells of a parent cell are present.

        ids.Sort();
        int output = 0;
        for (var i = 0; i < ids.Count;)
        {
            // var skipInc = false;
            var id = ids[i];
            Assert.True(id.IsValid());

            // Check whether this cell is contained by the previous cell.
            if (output <= 0 || !ids[output - 1].Contains(id))
            {
                // Discard any previous cells contained by this cell.
                while (output > 0 && id.Contains(ids[output - 1])) --output;

                // Check whether the last 3 elements plus "id" can be collapsed into a
                // single parent cell.
                while (output >= 3 && AreSiblings(ids[output - 3], ids[output - 2], ids[output - 1], id))
                {
                    // Replace four children by their parent cell.
                    id = id.Parent();
                    output -= 3;
                    // skipInc = true;
                }
                ids[output++] = id;
            }

            i++;
        }
        if (ids.Count != output)
            ids.RemoveRange(output, ids.Count - output);
    }

    // Like Denormalize(), but works with a vector of S2CellIds.
    // REQUIRES: output != &input
    public static void Denormalize(List<S2CellId> input, int min_level, int level_mod, List<S2CellId> output)
    {
        Assert.True(min_level >= 0);
        Assert.True(min_level <= S2.kMaxCellLevel);
        Assert.True(level_mod >= 1);
        Assert.True(level_mod <= 3);
        Assert.True(!ReferenceEquals(input, output));

        output.Clear();
        foreach (var id in input)
        {
            int level = id.Level();
            int new_level = Math.Max(min_level, level);
            if (level_mod > 1)
            {
                // Round up so that (new_level - min_level) is a multiple of level_mod.
                // (Note that S2Constants.kMaxCellLevel is a multiple of 1, 2, and 3.)
                new_level += (S2.kMaxCellLevel - (new_level - min_level)) % level_mod;
                new_level = Math.Min(S2.kMaxCellLevel, new_level);
            }
            if (new_level == level)
            {
                output.Add(id);
            }
            else
            {
                S2CellId end = id.ChildEnd(new_level);
                for (var id2 = id.ChildBegin(new_level); id2 != end; id2 = id2.Next())
                {
                    output.Add(id2);
                }
            }
        }
    }

    // Like GetIntersection(), but works directly with vectors of S2CellIds,
    // Equivalent to:
    //
    //    *output = S2CellUnion(x).Intersection(S2CellUnion(y)).Release()
    //
    // except that this method has slightly more relaxed normalization
    // requirements: the input vectors may contain groups of 4 child cells that
    // all have the same parent.  (In a normalized S2CellUnion, such groups are
    // always replaced by the parent cell.)
    public static void GetIntersection(List<S2CellId> x, List<S2CellId> y, List<S2CellId> output)
    {
        Assert.True(x.IsSorted());
        Assert.True(y.IsSorted());

        // This is a fairly efficient calculation that uses binary search to skip
        // over sections of both input vectors.  It takes logarithmic time if all the
        // cells of "x" come before or after all the cells of "y" in S2CellId order.

        output.Clear();
        var i = 0;
        var j = 0;
        while (i != x.Count && j != y.Count)
        {
            S2CellId imin = x[i].RangeMin();
            S2CellId jmin = y[j].RangeMin();
            if (imin > jmin)
            {
                // Either j->contains(*i) or the two cells are disjoint.
                if (x[i] <= y[j].RangeMax())
                {
                    output.Add(x[i++]);
                }
                else
                {
                    // Advance "j" to the first cell that might overlap *i.
                    j = y.GetLowerBound(x[i], j + 1, y.Count, new EntirelyPrecedes_Comparer());
                }
            }
            else if (jmin > imin)
            {
                // Identical to the code above with "i" and "j" reversed.
                if (y[j] <= x[i].RangeMax())
                {
                    output.Add(y[j++]);
                }
                else
                {
                    i = y.GetLowerBound(y[j], i + 1, x.Count, new EntirelyPrecedes_Comparer());
                }
            }
            else
            {
                // "i" and "j" have the same RangeMin, so one contains the other.
                if (x[i] < y[j])
                    output.Add(x[i++]);
                else
                    output.Add(y[j++]);
            }
        }
        // The output is generated in sorted order.
        Assert.True(output.IsSorted());
    }

    // Returns true if the given four cells have a common parent.
    // REQUIRES: The four cells are distinct.
    private static bool AreSiblings(S2CellId a, S2CellId b, S2CellId c, S2CellId d)
    {
        // A necessary (but not sufficient) condition is that the XOR of the
        // four cells must be zero.  This is also very fast to test.
        if ((a.Id ^ b.Id ^ c.Id) != d.Id) return false;

        // Now we do a slightly more expensive but exact test.  First, compute a
        // mask that blocks out the two bits that encode the child position of
        // "id" with respect to its parent, then check that the other three
        // children all agree with "mask".
        UInt64 mask = d.LowestOnBit() << 1;
        mask = ~(mask + (mask << 1));
        UInt64 id_masked = (d.Id & mask);
        return ((a.Id & mask) == id_masked &&
                (b.Id & mask) == id_masked &&
                (c.Id & mask) == id_masked &&
                !d.IsFace());
    }

    private static void GetDifferenceInternal(S2CellId cell, S2CellUnion y, List<S2CellId> cell_ids)
    {
        // Add the difference between cell and y to cell_ids.
        // If they intersect but the difference is non-empty, divide and conquer.
        if (!y.Intersects(cell))
        {
            cell_ids.Add(cell);
        }
        else if (!y.Contains(cell))
        {
            var child = cell.ChildBegin();
            for (int i = 0; ; ++i)
            {
                GetDifferenceInternal(child, y, cell_ids);
                if (i == 3) break;  // Avoid unnecessary Next computation.
                child = child.Next();
            }
        }
    }

    #endregion

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        // Compute the approximate centroid of the region.  This won't produce the
        // bounding cap of minimal area, but it should be close enough.
        if (!CellIds.Any()) return S2Cap.Empty;

        var centroid = new S2Point(0, 0, 0);
        foreach (S2CellId id in this)
        {
            double area = S2Cell.AverageArea(id.Level());
            centroid += area * id.ToPoint();
        }
        if (centroid == S2Point.Empty)
        {
            centroid = new S2Point(1, 0, 0);
        }
        else
        {
            centroid = centroid.Normalize();
        }

        // Use the centroid as the cap axis, and expand the cap angle so that it
        // contains the bounding caps of all the individual cells.  Note that it is
        // *not* sufficient to just bound all the cell vertices because the bounding
        // cap may be concave (i.e. cover more than one hemisphere).
        S2Cap cap = S2Cap.FromPoint(centroid);
        foreach (S2CellId id in this)
        {
            cap.AddCap(new S2Cell(id).GetCapBound());
        }
        return cap;
    }
    public S2LatLngRect GetRectBound()
    {
        S2LatLngRect bound = S2LatLngRect.Empty;
        foreach (S2CellId id in this)
        {
            bound = bound.Union(new S2Cell(id).GetRectBound());
        }
        return bound;
    }

    // This is a fast operation (logarithmic in the size of the cell union).
    public bool Contains(S2Cell cell)
    {
        return Contains(cell.Id);
    }

    // This is a fast operation (logarithmic in the size of the cell union).
    public bool MayIntersect(S2Cell cell)
    {
        return Intersects(cell.Id);
    }

    // The point 'p' does not need to be normalized.
    // This is a fast operation (logarithmic in the size of the cell union).
    public bool Contains(S2Point p)
    {
        return Contains(new S2CellId(p));
    }

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2CellUnion to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        // Unsigned char for version number, and N+1 UInt64's for N cell_ids
        // (1 for vector length, N for the ids).
        encoder.Ensure(sizeof(byte) + sizeof(UInt64) * (1 + CellIds.Count));

        encoder.Put8(S2.kCurrentLosslessEncodingVersionNumber);
        encoder.Put64((ulong)CellIds.Count);
        foreach (S2CellId cell_id in CellIds)
        {
            cell_id.Encode(encoder, hint);
        }
    }

    // Decodes an S2CellUnion encoded with Encode().  Returns true on success.
    public static (bool, S2CellUnion?) Decode(Decoder decoder)
    {
        // Should contain at least version and vector length.
        if (decoder.Avail() < sizeof(byte) + sizeof(UInt64))
            return (false, null);

        byte version = decoder.Get8();
        if (version > S2.kCurrentLosslessEncodingVersionNumber)
            return (false, null);

        var num_cells = (int)decoder.Get64();
        if (num_cells > Union_decode_max_num_cells)
        {
            return (false, null);
        }

        List<S2CellId> temp_cell_ids = new();
        for (var i = 0; i < num_cells; i++)
        {
            var (success, id) = S2CellId.Decode(decoder);
            if (!success)
                return (false, null);
            temp_cell_ids.Add(id!.Value);
        }
        var result = new S2CellUnion(temp_cell_ids);
        return (true, result);
    }

    public const int Union_decode_max_num_cells = 1000000; // The maximum number of cells allowed by S2CellUnion.Decode

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2CellUnion(CellIds, false);

    #endregion

    #region IEnumerable

    public IEnumerator<S2CellId> GetEnumerator()
    {
        return CellIds.GetEnumerator();
    }

    System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    #endregion

    #region IEquatable

    public override int GetHashCode()
    {
        return LinqUtils.GetSequenceHashCode(CellIds);
    }

    public virtual bool Equals(S2CellUnion? other)
    {
        if (other is null) return false;

        return CellIds.SequenceEqual(other.CellIds);
    }

    #endregion

    #region Object

    // Returns a human-readable string describing the S2CellUnion, consisting of
    // the number of cells and the list of S2CellIds in S2CellId::ToToken()
    // format (limited to at most 500 cells).
    public override string ToString()
    {
        const int kMaxCount = 500;
        System.Text.StringBuilder output = new();
        output.AppendFormat("Size:{0} S2CellIds:", CellIds.Count);
        for (int i = 0, limit = Math.Min(kMaxCount, CellIds.Count); i < limit; ++i)
        {
            if (i > 0) output.Append(",");
            output.Append(CellId(i).ToToken());
        }
        if (CellIds.Count > kMaxCount) output.Append(",...");
        return output.ToString();
    } 

    #endregion
}
