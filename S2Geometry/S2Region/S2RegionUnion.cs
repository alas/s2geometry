namespace S2Geometry;

// An S2RegionUnion represents a union of possibly overlapping regions.
// It is convenient for computing a covering of a set of regions.  However, note
// that currently, using S2RegionCoverer to compute coverings of S2RegionUnions
// may produce coverings with considerably less than the requested number of
// cells in cases of overlapping or tiling regions.  This occurs because the
// current S2RegionUnion.Contains implementation for S2Cells only returns
// true if the cell is fully contained by one of the regions.  So, cells along
// internal boundaries in the region union will be subdivided by the coverer
// even though this is unnecessary, using up the maxSize cell budget.  Then,
// when the coverer normalizes the covering, groups of 4 sibling cells along
// these internal borders will be replaced by parents, resulting in coverings
// that may have significantly fewer than maxSize cells, and so are less
// accurate.  This is not a concern for unions of disjoint regions.
//
public sealed record class S2RegionUnion : IS2Region<S2RegionUnion>
{
    #region Fields, Constants

    private List<IS2Region> Regions_ { get; set; }

    #endregion

    #region Constructors

    // Create a region representing the union of the given regions.
    public S2RegionUnion(List<IS2Region> regions, bool check = true)
    {
        Regions_ = regions;
        if (check)
        {
            Debug.Assert(!Regions_.Any());
        }
    }

    #endregion

    #region S2RegionUnion

    // Releases ownership of the regions of this union and returns them,
    // leaving this region empty.
    /*public List<IS2Region> Release()
    {
        var result = regions_;
        regions_ = null;
        return result;
    }*/

    // Add the given region to the union.  This method can be called repeatedly
    // as an alternative to Init().
    public void Add(S2RegionUnion region)
    {
        Regions_.Add(region);
    }

    // Accessor methods.
    public int Count() => Regions_.Count;

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        // TODO(ericv): This could be optimized to return a tighter bound,
        // but doesn't seem worth it unless profiling shows otherwise.
        return GetRectBound().GetCapBound();
    }

    public S2LatLngRect GetRectBound()
    {
        var result = S2LatLngRect.Empty;
        foreach (var r in Regions_)
        {
            result = result.Union(r.GetRectBound());
        }
        return result;
    }
    public bool Contains(S2Point p)
    {
        return Regions_.Any(t => t.Contains(p));
    }

    // The current implementation only returns true if one of the regions in the
    // union fully contains the cell.
    public bool Contains(S2Cell cell)
    {
        // Note that this method is allowed to return false even if the cell
        // is contained by the region.
        return Regions_.Any(t => t.Contains(cell));
    }
    public bool MayIntersect(S2Cell cell)
    {
        return Regions_.Any(t => t.MayIntersect(cell));
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2RegionUnion(Regions_.DeepCustomClone(), false);

    #endregion

    #region IEquatable

    public bool Equals(S2RegionUnion? other)
    {
        if (other is null) return false;

        return Regions_.SequenceEqual(other!.Regions_);
    }
    public override int GetHashCode() => LinqUtils.GetSequenceHashCode(Regions_);

    #endregion
}
