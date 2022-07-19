/// <summary>
/// An S2Point represents a point on the unit sphere as a 3D vector.  Usually
/// points are normalized to be unit length, but some methods do not require
/// this.  See util/math/vector.h for the methods available.  Among other
/// things, there are overloaded operators that make it convenient to write
/// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
/// </summary>
global using S2Point = S2Geometry.Vector3<double>;

namespace S2Geometry;

#if s2debug
public static class S2PointExtensions
{
    public static string ToStringDegrees(this S2Point p)
    {
        var s2LatLng = new S2LatLng(p);
        return $"({s2LatLng.LatDegrees():g}d, {s2LatLng.LngDegrees():g}d)";
    }

    /// <summary>
    /// Return true if the given point is approximately unit length
    /// (this is mainly useful for assertions).
    /// </summary>
    public static bool IsUnitLength(this S2Point p)
    {
        // Normalize() is guaranteed to return a vector whose L2-norm differs from 1
        // by less than 2 * DBL_EPSILON.  Thus the squared L2-norm differs by less
        // than 4 * DBL_EPSILON.  The actual calculated Norm2() can have up to 1.5 *
        // DBL_EPSILON of additional error.  The total error of 5.5 * DBL_EPSILON
        // can then be rounded down since the result must be a representable
        // double-precision value.
        return Math.Abs(p.Norm2() - 1) <= 5 * S2.DoubleEpsilon;  // About 1.11e-15
    }
}
#endif

/// <summary>
/// todo
/// </summary>
public static class LongDoubleExtensions
{
    public static double ToLD(this double d) => d;
}

/// <summary>
/// todo
/// </summary>
public static class ExactFloat
{
    /// <summary>
    /// The maximum exponent supported.  If a value has an exponent larger than
    /// this, it is replaced by infinity (with the appropriate sign).
    /// </summary>
    /// <remarks>
    /// About 10**(60 million)
    /// </remarks>
    public const int kMaxExp = 200 * 1000 * 1000;

    /// <summary>
    /// The minimum exponent supported.  If a value has an exponent less than
    /// this, it is replaced by zero (with the appropriate sign).
    /// </summary>
    /// <remarks>
    /// About 10**(-60 million)
    /// </remarks>
    public const int kMinExp = -kMaxExp;

    /// <summary>
    ///  The maximum number of mantissa bits supported.  If a value has more
    ///  mantissa bits than this, it is replaced with NaN.  (It is expected that
    ///  users of this class will never want this much precision.)
    /// </summary>
    /// <remarks>
    /// About 20 million digits
    /// </remarks>
    public const int kMaxPrec = 64 << 20;

    public static double FromDouble(double p) => p;

    public static double ToDouble(this double that) => that;

    /// <summary>
    /// Return +1 if this ExactFloat is positive, -1 if it is negative, and 0
    /// if it is zero or NaN.  Note that unlike sign_bit(), sgn() returns 0 for
    /// both positive and negative zero.
    /// </summary>
    public static int Sgn(this double d)
    {
        if (double.IsNaN(d) || d == 0 || d == -0 || d == 0.0 || d == -0.0) return 0;

        return d > 0 ? 1 : -1;
    }

    /// <summary>
    /// Return true if this value is a normal floating-point number.  Non-normal
    /// values (zero, infinity, and NaN) often need to be handled separately
    /// because they are represented using special exponent values and their
    /// mantissa is not defined.
    /// </summary>
    public static bool IsNormal(this double xf) =>
        xf != 0.0 && xf != 0 && !double.IsInfinity(xf) && !double.IsNaN(xf);

    /// <summary>
    /// Return the exponent of this ExactFloat given that the mantissa is in the
    /// range [0.5, 1).  It is an error to call this method if the value is zero,
    /// infinity, or NaN.
    /// </summary>
    public static int Exp(this double xf) =>
        1 + (int)Math.Floor(Math.Log10(Math.Abs(xf)));
}
