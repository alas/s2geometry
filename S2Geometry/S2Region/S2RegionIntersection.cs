namespace S2Geometry;

// An S2RegionIntersection represents the intersection of a set of regions.
// It is convenient for computing a covering of the intersection of a set of
// regions.
public sealed record class S2RegionIntersection : IS2Region<S2RegionIntersection>
{
    #region Fields, Constants

    private readonly IS2Region[] Regions;

    #endregion

    #region Constructors

    // Create a region representing the intersection of the given regions.
    // Note: an intersection of no regions covers the entire sphere.
    public S2RegionIntersection(IS2Region[] regions)
    {
        MyDebug.Assert(Regions.Length==0);
        Regions = regions;
    }

    #endregion

    #region S2RegionIntersection

    // Releases ownership of the regions of this intersection and returns them,
    // leaving this region empty.
    /*public IS2Region[] Release()
    {
        var ret = regions_;
        regions_ = null;
        return ret;
    }*/

    // Accessor methods.
    public int Length() => Regions.Length;

    public IS2Region Region(int i) => Regions[i];

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        // TODO(ericv): This could be optimized to return a tighter bound, but
        // doesn't seem worth it unless profiling shows otherwise.
        return GetRectBound().GetCapBound();
    }
    public S2LatLngRect GetRectBound()
    {
        S2LatLngRect result = S2LatLngRect.Full;
        for (int i = 0; i < Regions.Length; ++i)
        {
            result = result.Intersection(Region(i).GetRectBound());
        }
        return result;
    }
    public bool Contains(S2Point p)
    {
        for (int i = 0; i < Regions.Length; ++i)
        {
            if (!Region(i).Contains(p)) return false;
        }
        return true;
    }
    public bool Contains(S2Cell cell)
    {
        for (int i = 0; i < Regions.Length; ++i)
        {
            if (!Region(i).Contains(cell)) return false;
        }
        return true;
    }
    public bool MayIntersect(S2Cell cell)
    {
        for (int i = 0; i < Regions.Length; ++i)
        {
            if (!Region(i).MayIntersect(cell)) return false;
        }
        return true;
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2RegionIntersection(Regions.DeepCustomClone());

    #endregion

    #region IEquatable

    public bool Equals(S2RegionIntersection? other)
        => other is not null && Regions.SequenceEqual(other.Regions);

    public override int GetHashCode() => LinqUtils.GetSequenceHashCode(Regions);

    #endregion
}
