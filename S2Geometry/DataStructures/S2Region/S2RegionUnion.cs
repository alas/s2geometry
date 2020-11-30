using System;
using System.Collections.Generic;
using System.Linq;

namespace S2Geometry
{
	// An S2RegionUnion represents a union of possibly overlapping regions.
	// It is convenient for computing a covering of a set of regions.
	public sealed class S2RegionUnion : S2Region<S2RegionUnion>
	{
        #region Fields, Constants

        private List<IS2Region> regions_; 

        #endregion

        #region Constructors

        // Create an empty region.  Can be made non-empty by calling Init() or Add().
        public S2RegionUnion() { }

        // Create a region representing the union of the given regions.
        public S2RegionUnion(List<IS2Region> regions)
        {
            Init(regions);
        }

        #endregion

        #region S2RegionUnion

        // Initialize region by taking ownership of the given regions.
        public void Init(List<IS2Region> regions)
        {
            Assert.True(!regions_.Any());
            regions_ = regions;
        }

        // Releases ownership of the regions of this union and returns them,
        // leaving this region empty.
        public List<IS2Region> Release()
        {
            var result = regions_;
            regions_ = null;
            return result;
        }

        // Add the given region to the union.  This method can be called repeatedly
        // as an alternative to Init().
        public void Add(S2RegionUnion region)
        {
            regions_.Add(region);
        }

        // Accessor methods.
        public int Count => regions_.Count; 

        #endregion

        #region S2Region

        ////////////////////////////////////////////////////////////////////////
        // S2Region interface (see s2region.h for details):

        public override object Clone()
        {
            return new S2RegionUnion { regions_ = regions_.DeepCopy() };
        }

        public override S2Cap GetCapBound()
        {
            // TODO(ericv): This could be optimized to return a tighter bound,
            // but doesn't seem worth it unless profiling shows otherwise.
            return GetRectBound().GetCapBound();
        }

        public override S2LatLngRect GetRectBound()
        {
            var result = S2LatLngRect.Empty;
            foreach (var r in regions_)
            {
                result = result.Union(r.GetRectBound());
            }
            return result;
        }
        public override bool Contains(S2Point p)
        {
            return regions_.Any(t => t.Contains(p));
        }
        public override bool Contains(S2Cell cell)
        {
            // Note that this method is allowed to return false even if the cell
            // is contained by the region.
            return regions_.Any(t => t.Contains(cell));
        }
        public override bool MayIntersect(S2Cell cell)
        {
            return regions_.Any(t => t.MayIntersect(cell));
        }

        #endregion

        #region IEquatable

        public override bool Equals(S2RegionUnion other)
        {
            return regions_.SequenceEqual(other.regions_);
        }
        public override bool Equals(object other)
        {
            return base.Equals(other);
        }
        public static bool operator ==(S2RegionUnion left, S2RegionUnion right)
        {
            return Equals(left, right);
        }
        public static bool operator !=(S2RegionUnion left, S2RegionUnion right)
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
