using System;
using System.Linq;
namespace S2Geometry
{
    // An S2RegionIntersection represents the intersection of a set of regions.
    // It is convenient for computing a covering of the intersection of a set of
    // regions.
    public sealed class S2RegionIntersection : S2Region<S2RegionIntersection>, IEquatable<S2RegionIntersection>
    {
        #region Fields, Constants

        private IS2Region[] regions_;

        #endregion

        #region Constructors

        // Creates an empty intersection that should be initialized by calling Init().
        // Note: an intersection of no regions covers the entire sphere.
        public S2RegionIntersection() { }

        // Create a region representing the intersection of the given regions.
        public S2RegionIntersection(IS2Region[] regions) => Init(regions);

        #endregion

        #region S2RegionIntersection

        // Initialize region by taking ownership of the given regions.
        public void Init(IS2Region[] regions)
        {
            Assert.True(!regions_.Any());
            regions_ = regions;
        }

        // Releases ownership of the regions of this intersection and returns them,
        // leaving this region empty.
        public IS2Region[] Release()
        {
            var ret = regions_;
            regions_ = null;
            return ret;
        }

        // Accessor methods.
        public int Length => regions_.Length;
        public IS2Region Region(int i) { return regions_[i]; } 
        
        #endregion

        #region S2Region

        ////////////////////////////////////////////////////////////////////////
        // S2Region interface (see s2region.h for details):

        public override object Clone()
        {
            return new S2RegionIntersection { regions_ = regions_.DeepCopy() };
        }
        public override S2Cap GetCapBound()
        {
            // TODO(ericv): This could be optimized to return a tighter bound, but
            // doesn't seem worth it unless profiling shows otherwise.
            return GetRectBound().GetCapBound();
        }
        public override S2LatLngRect GetRectBound()
        {
            S2LatLngRect result = S2LatLngRect.Full;
            for (int i = 0; i < regions_.Length; ++i)
            {
                result = result.Intersection(Region(i).GetRectBound());
            }
            return result;
        }
        public override bool Contains(S2Point p)
        {
            for (int i = 0; i < regions_.Length; ++i)
            {
                if (!Region(i).Contains(p)) return false;
            }
            return true;
        }
        public override bool Contains(S2Cell cell)
        {
            for (int i = 0; i < regions_.Length; ++i)
            {
                if (!Region(i).Contains(cell)) return false;
            }
            return true;
        }
        public override bool MayIntersect(S2Cell cell)
        {
            for (int i = 0; i < regions_.Length; ++i)
            {
                if (!Region(i).MayIntersect(cell)) return false;
            }
            return true;
        }

        #endregion

        #region IEquatable

        public override bool Equals(S2RegionIntersection other)
        {
            return regions_.SequenceEqual(other.regions_);
        }
        public override bool Equals(object other)
        {
            return base.Equals(other);
        }
        public static bool operator ==(S2RegionIntersection left, S2RegionIntersection right)
        {
            return Equals(left, right);
        }
        public static bool operator !=(S2RegionIntersection left, S2RegionIntersection right)
        {
            return !Equals(left, right);
        }
        public override int GetHashCode()
        {
            return LinqUtils.GetSequenceHashCode(regions_);
        } 

        #endregion
    }
}
