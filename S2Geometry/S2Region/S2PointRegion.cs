// An S2PointRegion is a region that contains a single point.  It is more
// expensive than the raw S2Point type and is useful mainly for completeness.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.

namespace S2Geometry;

public readonly record struct S2PointRegion : IS2Region<S2PointRegion>, IDecoder<S2PointRegion>
{
    #region Fields, Constants

    public S2Point Point { get; init; }

    #endregion

    #region Constructors

    // Create a region containing the given point, which must be unit length.
    public S2PointRegion(S2Point point)
    {
        System.Diagnostics.Debug.Assert(point.IsUnitLength());
        Point = point;
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound() => S2Cap.FromPoint(Point);
    public S2LatLngRect GetRectBound()
    {
        var ll = new S2LatLng(Point);
        return new S2LatLngRect(ll, ll);
    }
    public bool Contains(S2Cell cell) => false;
    public bool MayIntersect(S2Cell cell) => cell.Contains(Point);
    public bool Contains(S2Point p) => Point == p;

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2Point to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(30); // sufficient

        encoder.Put8(S2.kCurrentLosslessEncodingVersionNumber);
        for (int i = 0; i < 3; ++i)
        {
            encoder.PutDouble(Point[i]);
        }
        System.Diagnostics.Debug.Assert(encoder.Avail() >= 0);
    }

    // Decodes an S2Point encoded with Encode().  Returns true on success.
    // (Returns false if the encoded point is not unit length.)
    public static (bool, S2PointRegion) Decode(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte) + 3 * sizeof(double))
            return (false, default);

        byte version = decoder.Get8();
        if (version > S2.kCurrentLosslessEncodingVersionNumber)
            return (false, default);

        var coords = new double[3];
        for (int i = 0; i < 3; ++i)
        {
            coords[i] = decoder.GetDouble();
        }
        var p = new S2Point(coords);
        if (!p.IsUnitLength()) return (false, default);

        return (true, new S2PointRegion(p));
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2PointRegion(Point);

    #endregion
}
