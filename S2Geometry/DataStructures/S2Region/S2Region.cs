using System;
using System.Collections.Generic;

namespace S2Geometry
{
	// An S2Region represents a two-dimensional region over the unit sphere.
	// It is an abstract interface with various concrete subtypes.
	//
	// The main purpose of this interface is to allow complex regions to be
	// approximated as simpler regions.  So rather than having a wide variety
	// of virtual methods that are implemented by all subtypes, the interface
	// is restricted to methods that are useful for computing approximations.
	public abstract class S2Region<T> : IEquatable<T>, IS2Region where T : IS2Region
	{
		// Returns a bounding spherical cap that contains the region.  The bound may
		// not be tight.
		public abstract S2Cap GetCapBound();

		// Returns a bounding latitude-longitude rectangle that contains the region.
		// The bound may not be tight.
		public abstract S2LatLngRect GetRectBound();

		// Returns a small collection of S2CellIds whose union covers the region.
		// The cells are not sorted, may have redundancies (such as cells that
		// contain other cells), and may cover much more area than necessary.
		//
		// This method is not intended for direct use by client code.  Clients
		// should typically use S2RegionCoverer.GetCovering, which has options to
		// control the size and accuracy of the covering.  Alternatively, if you
		// want a fast covering and don't care about accuracy, consider calling
		// S2RegionCoverer.GetFastCovering (which returns a cleaned-up version of
		// the covering computed by this method).
		//
		// GetCellUnionBound() implementations should attempt to return a small
		// covering (ideally 4 cells or fewer) that covers the region and can be
		// computed quickly.  The result is used by S2RegionCoverer as a starting
		// point for further refinement.
		//
		// TODO(ericv): Remove the default implementation.
		public virtual void GetCellUnionBound(List<S2CellId> cell_ids)
		{
			GetCapBound().GetCellUnionBound(cell_ids);
		}

		// Returns true if the region completely contains the given cell, otherwise
		// returns false.
		public abstract bool Contains(S2Cell cell);

		// If this method returns false, the region does not intersect the given
		// cell.  Otherwise, either region intersects the cell, or the intersection
		// relationship could not be determined.
		//
		// Note that there is currently exactly one implementation of this method
		// (S2LatLngRect.MayIntersect) that takes advantage of the semantics above
		// to be more efficient.  For all other S2Region subtypes, this method
		// returns true if the region intersect the cell and false otherwise.
		public abstract bool MayIntersect(S2Cell cell);

		// Returns true if and only if the given point is contained by the region.
		// The point 'p' is generally required to be unit length, although some
		// subtypes may relax this restriction.
		public abstract bool Contains(S2Point p);

		//////////////////////////////////////////////////////////////////////////
		// Many S2Region subtypes also define the following methods.
		//////////////////////////////////////////////////////////////////////////

		// public interface IS2RegionCodec<T> where T : S2Region<T>
		// {
		//     void Encode(Encoder encoder);
		//     static bool Decode(Decoder decoder, out T output);
		// }

		// Appends a serialized representation of the region to "encoder".
		//
		// The representation chosen is left up to the sub-classes but it should
		// satisfy the following constraints:
		// - It should encode a version number.
		// - It should be deserializable using the corresponding Decode method.
		// - Performance, not space, should be the chief consideration. Encode() and
		//   Decode() should be implemented such that the combination is equivalent
		//   to calling Clone().
		//
		// REQUIRES: "encoder" uses the defaultructor, so that its buffer
		//           can be enlarged as necessary by calling Ensure(int).
		//
		// public abstract void Encode(Encoder encoder);

		// Decodes an S2Region encoded with Encode().  Note that this method
		// requires that an S2Region object of the appropriate concrete type has
		// already been constructed.  It is not possible to decode regions of
		// unknown type.
		//
		// Whenever the Decode method is changed to deal with new serialized
		// representations, it should be done so in a manner that allows for
		// older versions to be decoded i.e. the version number in the serialized
		// representation should be used to decide how to decode the data.
		//
		// Returns true on success.
		//
		// public abstract bool Decode(Decoder decoder);

		// Provides the same functionality as Decode, except that decoded regions
		// are allowed to point directly into the Decoder's memory buffer rather
		// than copying the data.  This method can be much faster for regions that
		// have a lot of data (such as polygons), but the decoded region is only
		// valid within the scope (lifetime) of the Decoder's memory buffer.
		//
		// bool DecodeWithinScope(Decoder decoder);

		public static bool operator ==(S2Region<T> left, S2Region<T> right)
		{
			return Equals(left, right);
		}
		public static bool operator !=(S2Region<T> left, S2Region<T> right)
		{
			return !Equals(left, right);
		}

		public abstract bool Equals(T other);

		public override bool Equals(object other)
		{
			return other is T ot && Equals(ot);
		}

		public override abstract int GetHashCode();

		public abstract object Clone();
	}

	public interface IS2Region : ICloneable
	{
		bool Contains(S2Cell cell);
		bool MayIntersect(S2Cell cell);
		bool Contains(S2Point p);
		S2LatLngRect GetRectBound();
		void GetCellUnionBound(List<S2CellId> cell_ids);
	}
}