// This class represents a one-dimensional angle (as opposed to a
// two-dimensional solid angle).  It has methods for converting angles to
// or from radians, degrees, and the E5/E6/E7 representations (i.e. degrees
// multiplied by 1e5/1e6/1e7 and rounded to the nearest integer).
// 
// The internal representation is a double-precision value in radians, so
// conversion to and from radians is exact.  Conversions between E5, E6, E7,
// and Degrees are not always exact; for example, Degrees(3.1) is different
// from E6(3100000) or E7(310000000).  However, the following properties are
// guaranteed for any integer "n", provided that "n" is in the input range of
// both functions:
// 
//     Degrees(n) == E6(1000000 * n)
//     Degrees(n) == E7(10000000 * n)
//          E6(n) == E7(10 * n)
// 
// The corresponding properties are *not* true for E5, so if you use E5 then
// don't test for exact equality when comparing to other formats such as
// Degrees or E7.
// 
// The following conversions between degrees and radians are exact:
// 
//          Degrees(180) == Radians(Math.PI)
//       Degrees(45 * k) == Radians(k * Math.PI / 4)  for k == 0..8
// 
// These identities also hold when the arguments are scaled up or down by any
// power of 2.  Some similar identities are also true, for example,
// Degrees(60) == Radians(Math.PI / 3), but be aware that this type of identity
// does not hold in general.  For example, Degrees(3) != Radians(Math.PI / 60).
// 
// Similarly, the conversion to radians means that Angle.Degrees(x).degrees()
// does not always equal "x".  For example,
// 
//         Degrees(45 * k).degrees() == 45 * k      for k == 0..8
//   but       Degrees(60).degrees() != 60.
// 
// This means that when testing for equality, you should allow for numerical
// errors (Assert.True.Near) or convert to discrete E5/E6/E7 values first.
// 
// CAVEAT: All of the above properties depend on "double" being the usual
// 64-bit IEEE 754 type (which is true on almost all modern platforms).
// 
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.

namespace S2Geometry;

public readonly record struct S1Angle : IComparable<S1Angle>
{
    #region Fields, Constants

    public readonly double Radians { get; init; }

    public static readonly S1Angle Infinity = new(double.PositiveInfinity);
    public static readonly S1Angle Zero = new(0);

    #endregion

    #region Constructors

    private S1Angle(double radians) => Radians = radians;

    /// <summary>
    /// Return the angle between two points, which is also equal to the distance
    /// between these points on the unit sphere.  The points do not need to be
    /// normalized.  This function has a maximum error of 3.25 * S2Constants.DBL_EPSILON (or
    /// 2.5 * S2Constants.DBL_EPSILON for angles up to 1 radian). If either point is
    /// zero-length (e.g. an uninitialized S2Point), or almost zero-length, the
    /// resulting angle will be zero.
    /// </summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    public S1Angle(S2Point x, S2Point y)
        : this(x.Angle(y))
    {
    }

    /// <summary>
    /// Like the constructor above, but return the angle (i.e., distance) between
    /// two S2LatLng points.  This function has about 15 digits of accuracy for
    /// small distances but only about 8 digits of accuracy as the distance
    /// approaches 180 degrees (i.e., nearly-antipodal points).
    /// </summary>
    public S1Angle(S2LatLng x, S2LatLng y)
        : this(x.GetDistance(y).Radians)
    {
    }

    #endregion

    #region Factories

    // These methods construct S1Angle objects from their measure in radians
    // or degrees.

    public static S1Angle FromRadians(double radians) => new(radians);

    public static S1Angle FromDegrees(double degrees) => new(S2.M_PI_180 * degrees);

    public static S1Angle FromE5(long e5) => FromDegrees(e5 * 1e-5);

    public static S1Angle FromE6(long e6) =>
        // Multiplying by 1e-6 isn't quite as accurate as dividing by 1e6,
        // but it's about 10 times faster and more than accurate enough.
        FromDegrees(e6 * 1e-6);

    public static S1Angle FromE7(long e7) => FromDegrees(e7 * 1e-7);

    // Convenience functions -- to use when args have been fixed32s in protos.
    //
    // The arguments are static_cast into Int32, so very large unsigned values
    // are treated as negative numbers.

    public static S1Angle FromUnsignedE6(ulong e6) => FromE6((long)e6);

    public static S1Angle FromUnsignedE7(ulong e7) => FromE7((long)e7);

    #endregion

    #region S1Angle

    public double GetDegrees() => S2.M_180_PI * Radians;

    // Note that the E5, E6, and E7 conversion involve two multiplications rather
    // than one.  This is mainly for backwards compatibility (changing this would
    // break many tests), but it does have the nice side effect that conversions
    // between Degrees, E6, and E7 are exact when the arguments are integers.

    public long E5() => (long)Math.Round(GetDegrees() * 1e5);

    public long E6() => (long)Math.Round(GetDegrees() * 1e6);

    public long E7() => (long)Math.Round(GetDegrees() * 1e7);

    public static S1Angle Max(S1Angle left, S1Angle right) => right > left ? right : left;
    public static S1Angle Min(S1Angle left, S1Angle right) => right > left ? left : right;

    public double Sin() => Math.Sin(Radians);

    public double Cos() => Math.Cos(Radians);

    public double Tan() => Math.Tan(Radians);

    public double Abs() => Math.Abs(Radians);

    public static S1Angle Abs(S1Angle a) => new(Math.Abs(a.Radians));
    public S1Angle Normalize()
    {
        var rads = Math.IEEERemainder(Radians, S2.M_2_PI);
        if (rads <= -Math.PI) rads = Math.PI;

        return new S1Angle(rads);
    }

    #endregion

    #region IComparable

    public int CompareTo(S1Angle other)
    {
        if (Radians < other.Radians) return -1;
        if (Radians > other.Radians) return 1;
        return 0;
    }

    public static bool operator <(S1Angle x, S1Angle y) => x.Radians < y.Radians;

    public static bool operator >(S1Angle x, S1Angle y) => x.Radians > y.Radians;

    public static bool operator <=(S1Angle x, S1Angle y) => x.Radians <= y.Radians;

    public static bool operator >=(S1Angle x, S1Angle y) => x.Radians >= y.Radians;

    #endregion

    #region Operators

    public static S1Angle operator -(S1Angle x, S1Angle y) => FromRadians(x.Radians - y.Radians);
    public static S1Angle operator +(S1Angle a, S1Angle b) => FromRadians(a.Radians + b.Radians);
    public static S1Angle operator -(S1Angle a) => new(-a.Radians);
    public static S1Angle operator *(double m, S1Angle a) => FromRadians(m * a.Radians);
    public static S1Angle operator *(S1Angle a, double m) => FromRadians(m * a.Radians);
    public static S1Angle operator /(S1Angle a, double m) => FromRadians(a.Radians / m);
    public static double operator /(S1Angle a, S1Angle b) => a.Radians / b.Radians;

    #endregion

    #region Object

    /// </summary>
    /// Writes the angle in degrees with a "d" suffix, e.g. "17.3745d". By default
    /// 6 digits are printed; this can be changed using setprecision(). Up to 17
    /// digits are required to distinguish one angle from another.
    /// </summary>
    public override string ToString() => $"{GetDegrees()}d";

    #endregion
}
