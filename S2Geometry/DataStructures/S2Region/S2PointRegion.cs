using System;
namespace S2Geometry
{
    // An S2PointRegion is a region that contains a single point.  It is more
    // expensive than the raw S2Point type and is useful mainly for completeness.
    //
    // This class is intended to be copied by value as desired.  It uses
    // the default copy constructor and assignment operator.
    public sealed class S2PointRegion : S2Region<S2PointRegion>, IEquatable<S2PointRegion>, ICoder
    {
        #region Fields, Constants

        public S2Point Point { get; } 
        
        #endregion

        #region Constructors

        // Create a region containing the given point, which must be unit length.
        public S2PointRegion(S2Point point)
        {
            Point = point;
            Assert.True(point.IsUnitLength);
        }

        #endregion

        #region S2Region

        ////////////////////////////////////////////////////////////////////////
        // S2Region interface (see s2region.h for details):

        public override object Clone()
        {
            return new S2PointRegion(Point);
        }
        public override S2Cap GetCapBound()
        {
            return S2Cap.FromPoint(Point);
        }
        public override S2LatLngRect GetRectBound()
        {
            var ll = new S2LatLng(Point);
            return new S2LatLngRect(ll, ll);
        }
        public override bool Contains(S2Cell cell) { return false; }
        public override bool MayIntersect(S2Cell cell)
        {
            return cell.Contains(Point);
        }
        public override bool Contains(S2Point p) { return Point == p; } 
        
        #endregion

        #region ICoder

        // Appends a serialized representation of the S2Point to "encoder".
        //
        // REQUIRES: "encoder" uses the defaultructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public void Encode(Encoder encoder)
        {
            encoder.Ensure(30);  // sufficient

            encoder.Put8(S2Constants.kCurrentLosslessEncodingVersionNumber);
            for (int i = 0; i < 3; ++i)
            {
                encoder.PutDouble(Point[i]);
            }
            Assert.True(encoder.Avail() >= 0);
        }

        // Decodes an S2Point encoded with Encode().  Returns true on success.
        // (Returns false if the encoded point is not unit length.)
        public static (bool success, S2PointRegion result) DecodeStatic(Decoder decoder)
        {
            if (decoder.Avail() < sizeof(byte) + 3 * sizeof(double))
                return (false, null);

            byte version = decoder.Get8();
            if (version > S2Constants.kCurrentLosslessEncodingVersionNumber)
                return (false, null);

            var coords = new double[3];
            for (int i = 0; i < 3; ++i)
            {
                coords[i] = decoder.GetDouble();
            }
            var p = new S2Point(coords);
            if (!p.IsUnitLength) return (false, null);

            return (true, new S2PointRegion(p));
        }
        
        #endregion

        #region IEquatable

        public static bool operator ==(S2PointRegion left, S2PointRegion right)
        {
            return Equals(left, right);
        }
        public static bool operator !=(S2PointRegion left, S2PointRegion right)
        {
            return !Equals(left, right);
        }
        public override bool Equals(object other)
        {
            return base.Equals(other);
        }
        public override bool Equals(S2PointRegion r)
        {
            return Point == r.Point;
        }
        public override int GetHashCode()
        {
            return Point.GetHashCode();
        }

        #endregion
    }
}
