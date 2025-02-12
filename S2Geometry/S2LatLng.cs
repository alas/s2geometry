﻿// This class represents a point on the unit sphere as a pair
// of latitude-longitude coordinates.  Like the rest of the "geometry"
// package, the intent is to represent spherical geometry as a mathematical
// abstraction, so functions that are specifically related to the Earth's
// geometry (e.g. easting/northing conversions) should be put elsewhere.
// 
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.

namespace S2Geometry;

public readonly record struct S2LatLng : IComparable<S2LatLng>
{
    #region Fields, Constants

    /// <summary>
    /// Returns the latitude of this point as radians.
    /// </summary>
    public readonly double LatRadians { get; init; }

    /// <summary>
    /// Returns the longitude of this point as radians.
    /// </summary>
    public readonly double LngRadians { get; init; }

    /// <summary>
    /// The center point the lat/lng coordinate system.
    /// </summary>
    public static readonly S2LatLng Center = new(0.0, 0.0);

    /// <summary>
    /// Returns an S2LatLng for which is_valid() will return false.
    /// </summary>
    /// <remarks>
    /// These coordinates are outside the bounds allowed by IsValid.
    /// </remarks>
    public static readonly S2LatLng Invalid = new(Math.PI, S2.M_2_PI);

    #endregion

    #region Constructors

    /// <summary>
    /// This is internal to avoid ambiguity about which units are expected.
    /// </summary>
    private S2LatLng(double latRadians, double lngRadians)
    {
        LatRadians = latRadians;
        LngRadians = lngRadians;
    }

    public S2LatLng(R2Point p)
        : this(p.X, p.Y)
    {
    }

    /// <summary>
    /// Constructor.  The latitude and longitude are allowed to be outside
    /// the is_valid() range.  However, note that most methods that accept
    /// S2LatLngs expect them to be normalized (see Normalized() below).
    /// </summary>
    public S2LatLng(S1Angle lat, S1Angle lng)
        : this(lat.Radians, lng.Radians)
    {
    }

    /// <summary>
    /// Convert a direction vector (not necessarily unit length) to an S2LatLng.
    /// </summary>
    public S2LatLng(S2Point p)
        : this(Latitude(p).Radians, Longitude(p).Radians)
    {
        // The latitude and longitude are already normalized.
        if (!IsValid()) MyDebug.WriteLine($"Invalid S2LatLng in constructor {this}");
    }

    #endregion

    #region Factories

    public static S2LatLng FromRadians(double latRadians, double lngRadians) => new(latRadians, lngRadians);

    public static S2LatLng FromDegrees(double latDegrees, double lngDegrees) => new(S1Angle.FromDegrees(latDegrees), S1Angle.FromDegrees(lngDegrees));

    public static S2LatLng FromE5(Int32 latE5, Int32 lngE5) => new(S1Angle.FromE5(latE5), S1Angle.FromE5(lngE5));

    public static S2LatLng FromE5(long latE5, long lngE5) => new(S1Angle.FromE5(latE5), S1Angle.FromE5(lngE5));

    public static S2LatLng FromE6(Int32 latE6, Int32 lngE6) => new(S1Angle.FromE6(latE6), S1Angle.FromE6(lngE6));

    public static S2LatLng FromE6(long latE6, long lngE6) => new(S1Angle.FromE6(latE6), S1Angle.FromE6(lngE6));

    public static S2LatLng FromE7(Int32 latE7, Int32 lngE7) => new(S1Angle.FromE7(latE7), S1Angle.FromE7(lngE7));

    public static S2LatLng FromE7(long latE7, long lngE7) => new(S1Angle.FromE7(latE7), S1Angle.FromE7(lngE7));

    // Convenience functions -- to use when args have been fixed32s in protos.
    //
    // The arguments are static_cast into int32, so very large unsigned values
    // are treated as negative numbers.

    public static S2LatLng FromUnsignedE6(uint latE6, uint lngE6) => new(S1Angle.FromUnsignedE6(latE6), S1Angle.FromUnsignedE6(lngE6));

    public static S2LatLng FromUnsignedE7(uint latE7, uint lngE7) => new(S1Angle.FromUnsignedE7(latE7), S1Angle.FromUnsignedE7(lngE7));

    /// <summary>
    /// Convert a direction vector (not necessarily unit length) to an S2LatLng.
    /// </summary>
    /// <param name="p"></param>
    public static S2LatLng FromPoint(S2Point p)
    {
        var a = Math.Sqrt(p.X * p.X + p.Y * p.Y);
        var b = Math.Atan2(p.Z, a);
        var c = Math.Atan2(p.Y, p.X);
        return new S2LatLng(b, c);
    }

    #endregion

    #region S2LatLng

    /// <summary>
    /// Converts to an S2Point (equivalent unit-length vector).
    /// The maximum error in the result is 1.5 * S2Constants.DoubleEpsilon.  (This does not
    /// include the error of converting degrees, E5, E6, or E7 to radians.)
    /// </summary>
    public S2Point ToPoint()
    {
        if (!IsValid()) MyDebug.WriteLine($"Invalid S2LatLng in ToPoint {this}");
        var phi = LatRadians;
        var theta = LngRadians;
        var cosphi = Math.Cos(phi);
        return new S2Point(Math.Cos(theta) * cosphi, Math.Sin(theta) * cosphi, Math.Sin(phi));
    }

    public S1Angle Lat() => S1Angle.FromRadians(LatRadians);

    /// <summary>
    /// Returns the latitude of this point as degrees.
    /// </summary>
    public double LatDegrees() => LatRadians * S2.M_180_PI;

    /// <summary>
    /// Returns the longitude of this point as a new S1Angle.
    /// </summary>
    public S1Angle Lng() => S1Angle.FromRadians(LngRadians);

    /// <summary>
    /// Returns the longitude of this point as degrees.
    /// </summary>
    public double LngDegrees() => LngRadians * S2.M_180_PI;

    #region Methods to compute the latitude and longitude of a point separately. 

    public static S1Angle Latitude(S2Point p)
    {
        // We use atan2 rather than asin because the input vector is not necessarily
        // unit length, and atan2 is much more accurate than asin near the poles.
        // The "+ 0.0" is to ensure that points with coordinates of -0.0 and +0.0
        // (which compare equal) are converted to identical S2LatLng values, since
        // even though -0.0 == +0.0 they can be formatted differently.
        var sqrt = Math.Sqrt(p[0] * p[0] + p[1] * p[1]);
        var atan = Math.Atan2(p[2] + 0.0, sqrt);
        return S1Angle.FromRadians(atan);
    }

    public static S1Angle Longitude(S2Point p) =>
        // The "+ 0.0" is to ensure that points with coordinates of -0.0 and +0.0
        // (which compare equal) are converted to identical S2LatLng values, since
        // even though -0.0 == +0.0 and -180 == 180 degrees, they can be formatted
        // differently.  Also note that atan2(0, 0) is defined to be zero.
        S1Angle.FromRadians(Math.Atan2(p[1] + 0.0, p[0] + 0.0));

    #endregion

    /// <summary>
    /// Returns the distance (measured along the surface of the sphere) to the
    /// given S2LatLng, implemented using the Haversine formula.  This is
    /// equivalent to
    /// 
    ///   S1Angle(ToPoint(), o.ToPoint())
    /// 
    /// except that this function is slightly faster, and is also somewhat less
    /// accurate for distances approaching 180 degrees (see s1angle.h for
    /// details).  Both S2LatLngs must be normalized.
    /// </summary>
    public S1Angle GetDistance(S2LatLng o)
    {
        // This implements the Haversine formula, which is numerically stable for
        // small distances but only gets about 8 digits of precision for very large
        // distances (e.g. antipodal points).  Note that 8 digits is still accurate
        // to within about 10cm for a sphere the size of the Earth.
        //
        // This could be fixed with another sin() and cos() below, but at that point
        // you might as well just convert both arguments to S2Points and compute the
        // distance that way (which gives about 15 digits of accuracy for all
        // distances).

        if (!IsValid()) MyDebug.WriteLine($"Invalid S2LatLng in GetDistance {this}");
        if (!o.IsValid()) MyDebug.WriteLine($"Invalid S2LatLng in GetDistance {o}");

        var lat1 = LatRadians;
        var lat2 = o.LatRadians;
        var lng1 = LngRadians;
        var lng2 = o.LngRadians;
        var dlat = Math.Sin(0.5 * (lat2 - lat1));
        var dlng = Math.Sin(0.5 * (lng2 - lng1));
        var x = dlat * dlat + dlng * dlng * Math.Cos(lat1) * Math.Cos(lat2);
        return S1Angle.FromRadians(2 * Math.Asin(Math.Sqrt(Math.Min(1.0, x))));
    }

    /// <summary>
    /// Returns the surface distance to the given point assuming a constant radius.
    /// </summary>
    public double GetDistance(S2LatLng o, double radius) => GetDistance(o).Radians * radius;

    // TODO(dbeaumont): Maybe check that radius >= 0 ?

    /// <summary>
    /// Return true if the latitude is between -90 and 90 degrees inclusive
    /// and the longitude is between -180 and 180 degrees inclusive.
    /// </summary>
    public bool IsValid() => Math.Abs(LatRadians) <= S2.M_PI_2 && Math.Abs(LngRadians) <= Math.PI;

    /// <summary>
    /// Clamps the latitude to the range [-90, 90] degrees, and adds or subtracts
    /// a multiple of 360 degrees to the longitude if necessary to reduce it to
    /// the range [-180, 180].
    /// </summary>
    public S2LatLng Normalized()
    {
        // remainder(x, S2Constants.M_2_PI) reduces its argument to the range [-Math.PI, Math.PI]
        // inclusive, which is what we want here.
        var min = Math.Min(S2.M_PI_2, LatRadians);
        var max = Math.Max(-S2.M_PI_2, min);
        var rem = Math.IEEERemainder(LngRadians, S2.M_2_PI);
        return new S2LatLng(max, rem);
    }

    /// <summary>
    /// Returns true if both the latitude and longitude of the given point are
    /// within {@code maxError} radians of this point.
    /// </summary>
    public bool ApproxEquals(S2LatLng o, double maxError = S2.DoubleError)
        => (Math.Abs(LatRadians - o.LatRadians) < maxError)
        && (Math.Abs(LngRadians - o.LngRadians) < maxError);

    /// <summary>
    /// Returns true if the given point is within {@code 1e-9} radians of this
    /// point. This corresponds to a distance of less than {@code 1cm} at the
    /// surface of the Earth.
    /// </summary>
    public bool ApproxEquals1cm(S2LatLng o) => ApproxEquals(o, 1e-9);

    #endregion

    #region IComparable

    public int CompareTo(S2LatLng other)
    {
        if (LatRadians < other.LatRadians) return -1;
        if (LatRadians > other.LatRadians) return 1;
        if (LngRadians < other.LngRadians) return -1;
        if (LngRadians > other.LngRadians) return 1;
        return 0;
    }

    public static bool operator <(S2LatLng left, S2LatLng right) => left.LatRadians < right.LatRadians && left.LngRadians < right.LngRadians;

    public static bool operator >(S2LatLng left, S2LatLng right) => left.LatRadians > right.LatRadians && left.LngRadians > right.LngRadians;

    public static bool operator <=(S2LatLng left, S2LatLng right) => left.LatRadians <= right.LatRadians && left.LngRadians <= right.LngRadians;

    public static bool operator >=(S2LatLng left, S2LatLng right) => left.LatRadians >= right.LatRadians && left.LngRadians >= right.LngRadians;

    #endregion

    #region Operators

    // Simple arithmetic operations for manipulating latitude-longitude pairs.
    // The results are not normalized (see Normalized()).

    /// <summary>
    /// Adds the given point to this point.
    /// Note that there is no guarantee that the new point will be <em>valid</em>.
    /// </summary>
    public static S2LatLng operator +(S2LatLng x, S2LatLng y) => new(x.LatRadians + y.LatRadians, x.LngRadians + y.LngRadians);

    /// <summary>
    /// Subtracts the given point from this point.
    /// Note that there is no guarantee that the new point will be <em>valid</em>.
    /// </summary>
    public static S2LatLng operator -(S2LatLng x, S2LatLng y) => new(x.LatRadians - y.LatRadians, x.LngRadians - y.LngRadians);

    public static S2LatLng operator *(double m, S2LatLng a) => new(a.LatRadians * m, a.LngRadians * m);

    public static S2LatLng operator *(S2LatLng x, double m) => new(x.LatRadians * m, x.LngRadians * m);

    #endregion

    #region Object

    public override string ToString() => $"[{Lat()}, {Lng()}]";

    #endregion
}
