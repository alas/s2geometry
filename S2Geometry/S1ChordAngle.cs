// S1ChordAngle represents the angle subtended by a chord (i.e., the straight
// line segment connecting two points on the sphere).  Its representation
// makes it very efficient for computing and comparing distances, but unlike
// S1Angle it is only capable of representing angles between 0 and Pi radians.
// S1ChordAngle is intended for applications where many angles need to be
// computed and compared, otherwise it is simpler to use S1Angle.
// 
// S1ChordAngle also loses some accuracy as the angle approaches Pi radians.
// There are several different ways to measure this error, including the
// representational error (i.e., how accurately S1ChordAngle can represent
// angles near Pi radians), the conversion error (i.e., how much precision is
// lost when an S1Angle is converted to an S1ChordAngle), and the measurement
// error (i.e., how accurate the S1ChordAngle(a, b) constructor is when the
// points A and B are separated by angles close to Pi radians).  All of these
// errors differ by a small constant factor.
//
// For the measurement error (which is the largest of these errors and also
// the most important in practice), let the angle between A and B be (Pi - x)
// radians, i.e. A and B are within "x" radians of being antipodal.  The
// corresponding chord length is
//
//    r = 2 * sin((Pi - x) / 2) = 2 * cos(x / 2) .
//
// For values of x not close to Pi the relative error in the squared chord
// length is at most 4.5 * DBL_EPSILON (see GetS2PointConstructorMaxError).
// The relative error in "r" is thus at most 2.25 * DBL_EPSILON ~= 5e-16.  To
// convert this error into an equivalent angle, we have
//
//    |dr / dx| = sin(x / 2)
//
// and therefore
//
//    |dx| = dr / sin(x / 2)
//         = 5e-16 * (2 * cos(x / 2)) / sin(x / 2)
//         = 1e-15 / tan(x / 2)
//
// The maximum error is attained when
//
//    x  = |dx|
//       = 1e-15 / tan(x / 2)
//      ~= 1e-15 / (x / 2)
//      ~= sqrt(2e-15)
//
// In summary, the measurement error for an angle (Pi - x) is at most
//
//    dx  = min(1e-15 / tan(x / 2), sqrt(2e-15))
//      (~= min(2e-15 / x, sqrt(2e-15)) when x is small).
//
// On the Earth's surface (assuming a radius of 6371km), this corresponds to
// the following worst-case measurement errors:
//
//     Accuracy:             Unless antipodal to within:
//     ---------             ---------------------------
//     6.4 nanometers        10,000 km (90 degrees)
//     1 micrometer          81.2 kilometers
//     1 millimeter          81.2 meters
//     1 centimeter          8.12 meters
//     28.5 centimeters      28.5 centimeters
//
// The representational and conversion errors referred to earlier are somewhat
// smaller than this.  For example, maximum distance between adjacent
// representable S1ChordAngle values is only 13.5 cm rather than 28.5 cm.  To
// see this, observe that the closest representable value to r^2 = 4 is
// r^2 =  4 * (1 - DBL_EPSILON / 2).  Thus r = 2 * (1 - DBL_EPSILON / 4) and
// the angle between these two representable values is
//
//    x  = 2 * acos(r / 2)
//       = 2 * acos(1 - DBL_EPSILON / 4)
//      ~= 2 * asin(sqrt(DBL_EPSILON / 2)
//      ~= sqrt(2 * DBL_EPSILON)
//      ~= 2.1e-8
//
// which is 13.5 cm on the Earth's surface.
//
// The worst case rounding error occurs when the value halfway between these
// two representable values is rounded up to 4.  This halfway value is
// r^2 = (4 * (1 - DBL_EPSILON / 4)), thus r = 2 * (1 - DBL_EPSILON / 8) and
// the worst case rounding error is
//
//    x  = 2 * acos(r / 2)
//       = 2 * acos(1 - DBL_EPSILON / 8)
//      ~= 2 * asin(sqrt(DBL_EPSILON / 4)
//      ~= sqrt(DBL_EPSILON)
//      ~= 1.5e-8
//
// which is 9.5 cm on the Earth's surface.
// 
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
// </summary>
// <Remarks>
// When S1ChordAngle is used as a key in one of the btree container types
// (util/btree), indicate that linear rather than binary search should be
// used.  This is much faster when the comparison function is cheap.

namespace S2Geometry;

public readonly record struct S1ChordAngle : IComparable<S1ChordAngle>, IDistance<S1ChordAngle>
{
    #region Fields, Constants

    /// <summary>
    /// The squared length of the chord.  (Most clients will not need this.)
    /// </summary>
    public readonly double Length2 { get; init; }

    private const double kMaxLength2 = 4.0;

    public static S1ChordAngle Zero { get; } = new(0);

    /// <summary>
    /// Return a chord angle of 90 degrees (a "right angle").
    /// </summary>
    public static readonly S1ChordAngle Right = new(2);

    /// <summary>
    /// Return a chord angle of 180 degrees (a "straight angle").  This is the
    /// maximum finite chord angle.
    /// </summary>
    public static readonly S1ChordAngle Straight = new(kMaxLength2);

    /// <summary>
    /// Return a chord angle larger than any finite chord angle.  The only valid
    /// operations on Infinity() are comparisons, S1Angle conversions, and
    /// Successor() / Predecessor().
    /// </summary>
    public static S1ChordAngle Infinity { get; } = new(S1Angle.Infinity);

    /// <summary>
    /// Return a chord angle smaller than Zero().  The only valid operations on
    /// Negative() are comparisons, S1Angle conversions, and Successor() /
    /// Predecessor().
    /// </summary>
    public static readonly S1ChordAngle Negative = new(-1);

    #endregion

    #region Constructors

    /// <summary>
    /// Construct the S1ChordAngle corresponding to the distance between the two
    /// given points.  The points must be unit length.
    /// </summary>
    public S1ChordAngle(S2Point x, S2Point y)
    {
        System.Diagnostics.Debug.Assert(x.IsUnitLength());
        System.Diagnostics.Debug.Assert(y.IsUnitLength());
        // The squared distance may slightly exceed kMaxLength2 due to roundoff errors.
        // The maximum error in the result is 2 * S2Constants.DoubleEpsilon * length2_.
        Length2 = Math.Min(kMaxLength2, (x - y).Norm2());
        System.Diagnostics.Debug.Assert(IsValid());
    }

    /// <summary>
    /// Conversion from an S1Angle.  Angles outside the range [0, Pi] are handled
    /// as follows: Infinity() is mapped to Infinity(), negative angles are
    /// mapped to Negative(), and finite angles larger than Pi are mapped to
    /// Straight().
    /// 
    /// Note that this operation is relatively expensive and should be avoided.
    /// To use S1ChordAngle effectively, you should structure your code so that
    /// input arguments are converted to S1ChordAngles at the beginning of your
    /// algorithm, and results are converted back to S1Angles only at the end.
    /// </summary>
    public S1ChordAngle(S1Angle angle)
    {
        if (angle.Radians < 0.0)
        {
            Length2 = -1.0;
        }
        else if (angle == S1Angle.Infinity)
        {
            Length2 = double.PositiveInfinity;
        }
        else
        {
            var ang = angle.Radians;
            // NOTE(Alas): angle can be NaN from dividing by zero and in C++ min(M_PI, NaN) == M_PI, in C# Math.Min(Math.PI, double.NaN) == double.NaN
            if (double.IsNaN(ang)) ang = double.PositiveInfinity;

            // The chord length is 2 * sin(angle / 2).
            var length = 2 * Math.Sin(0.5 * Math.Min(Math.PI, ang));
            // NOTE(Alas): sqrt(2) * sqrt(2) != 2 for System.Double
            if (length == Math.Sqrt(2))
                Length2 = 2;
            else
                Length2 = length * length;
        }
        System.Diagnostics.Debug.Assert(IsValid());
    }

    /// <summary>
    /// S1ChordAngles are represented by the squared chord length, which can
    /// range from 0 to 4.  Infinity() uses an infinite squared length.
    /// </summary>
    /// <param name="length2"></param>
    private S1ChordAngle(double length2)
    {
        Length2 = length2;
        System.Diagnostics.Debug.Assert(IsValid());
    }

    #endregion

    #region Factories

    /// <summary>
    /// Construct an S1ChordAngle from the squared chord length.  Note that the
    /// argument is automatically clamped to a maximum of kMaxLength2 to handle possible
    /// roundoff errors.  The argument must be non-negative.
    /// </summary>
    public static S1ChordAngle FromLength2(double length2) => new(Math.Min(kMaxLength2, length2));

    public static S1ChordAngle FromRadians(double radians) => new(S1Angle.FromRadians(radians));

    public static S1ChordAngle FromDegrees(double degrees) => new(S1Angle.FromDegrees(degrees));

    public static S1ChordAngle FromE5(long e5) => new(S1Angle.FromE5(e5));

    public static S1ChordAngle FromE6(long e6) => new(S1Angle.FromE6(e6));

    public static S1ChordAngle FromE7(long e7) => new(S1Angle.FromE7(e7));

    /// <summary>
    /// Construct an S1ChordAngle that is an upper bound on the given S1Angle,
    /// i.e. such that FastUpperBoundFrom(x).ToAngle() >= x.  Unlike the S1Angle
    /// constructor, this method is very fast, and the bound is accurate to
    /// within 1% for distances up to about 3100km on the Earth's surface.
    /// </summary>
    public static S1ChordAngle FastUpperBoundFrom(S1Angle angle) =>
        // This method uses the distance along the surface of the sphere as an upper
        // bound on the distance through the sphere's interior.
        FromLength2(angle.Radians * angle.Radians);

    #endregion

    #region S1ChordAngle

    /// <summary>
    /// Return true if the internal representation is valid.  Negative() and
    /// Infinity() are both considered valid.
    /// </summary>
    public bool IsValid() => (Length2 >= 0 && Length2 <= kMaxLength2) || IsSpecial();

    /// <summary>
    /// Converts to an S1Angle.
    /// 
    /// Infinity() is converted to S1Angle.Infinity(), and Negative() is
    /// converted to an unspecified negative S1Angle.
    /// 
    /// Note that the conversion uses trigonometric functions and therefore
    /// should be avoided in inner loops.
    /// </summary>
    public S1Angle ToAngle()
    {
        if (IsNegative()) return S1Angle.FromRadians(-1);
        if (IsInfinity()) return S1Angle.Infinity;
        return S1Angle.FromRadians(2 * Math.Asin(0.5 * Math.Sqrt(Length2)));
    }

    // Convenience methods implemented by calling ToAngle() first.  Note that
    // because of the S1Angle conversion these methods are relatively expensive
    // (despite their lowercase names), so the results should be cached if they
    // are needed inside loops.

    public double Radians() => ToAngle().Radians;

    public double Degrees() => ToAngle().GetDegrees();

    public long E5() => ToAngle().E5();

    public long E6() => ToAngle().E6();

    public long E7() => ToAngle().E7();

    public bool IsZero() => Length2 == 0.0;

    // TODO(ericv): Consider stricter check here -- only allow Negative().
    public bool IsNegative() => Length2 < 0.0;

    public bool IsInfinity() => Length2 == double.PositiveInfinity;

    /// <summary>
    /// Negative or infinity
    /// </summary>
    public bool IsSpecial() => IsNegative() || IsInfinity();

    public static S1ChordAngle Max(S1ChordAngle left, S1ChordAngle right) =>
        right > left ? right : left;

    public static S1ChordAngle Min(S1ChordAngle left, S1ChordAngle right) =>
        right > left ? left : right;

    /// Returns the smallest representable S1ChordAngle larger than this object.
    /// This can be used to convert a "<" comparison to a "<=" comparison.  For
    /// example:
    ///
    ///   S2ClosestEdgeQuery query(...);
    ///   S1ChordAngle limit = ...;
    ///   if (query.IsDistanceLess(target, limit.Successor())) {
    ///     // Distance to "target" is less than or equal to "limit".
    ///   }
    ///
    /// Note the following special cases:
    ///   Negative().Successor() == Zero()
    ///   Straight().Successor() == Infinity()
    ///   Infinity().Successor() == Infinity()
    public S1ChordAngle Successor()
    {
        if (Length2 >= kMaxLength2) return Infinity;
        if (Length2 < 0.0) return Zero;
        return new S1ChordAngle(MathUtils.NextAfter(Length2, 10.0));
    }

    /// Like Successor(), but returns the largest representable S1ChordAngle less
    /// than this object.
    ///
    /// Note the following special cases:
    ///   Infinity().Predecessor() == Straight()
    ///   Zero().Predecessor() == Negative()
    ///   Negative().Predecessor() == Negative()
    public S1ChordAngle Predecessor()
    {
        if (Length2 <= 0.0) return Negative;
        if (Length2 > kMaxLength2) return Straight;
        return new S1ChordAngle(MathUtils.NextAfter(Length2, -10.0));
    }

    /// <summary>
    /// Returns a new S1ChordAngle that has been adjusted by the given error
    /// bound (which can be positive or negative).  "error" should be the value
    /// returned by one of the error bound methods below.  For example:
    ///    S1ChordAngle a(x, y);
    ///    S1ChordAngle a1 = a.PlusError(a.GetS2PointConstructorMaxError());
    /// </summary>
    public S1ChordAngle PlusError(double error)
    {
        // If angle is Negative() or Infinity(), don't change it.
        // Otherwise clamp it to the valid range.
        if (IsSpecial()) return this;
        return new S1ChordAngle(Math.Max(0.0, Math.Min(kMaxLength2, Length2 + error)));
    }

    /// <summary>
    /// Return the maximum error in length2() for the S1ChordAngle(x, y)
    /// constructor, assuming that "x" and "y" are normalized to within the
    /// bounds guaranteed by S2Point.Normalize().  (The error is defined with
    /// respect to the true distance after the points are projected to lie
    /// exactly on the sphere.)
    /// </summary>
    public double GetS2PointConstructorMaxError() =>
        // There is a relative error of 2.5 * S2Constants.DoubleEpsilon when computing the squared
        // distance, plus a relative error of 2 * S2Constants.DoubleEpsilon and an absolute error
        // of (16 * S2Constants.DoubleEpsilon**2) because the lengths of the input points may
        // differ from 1 by up to (2 * S2Constants.DoubleEpsilon) each.  (This is the maximum
        // length error in S2Point.Normalize.)
        4.5 * S2.DoubleEpsilon * Length2 + 16 * S2.DoubleEpsilon * S2.DoubleEpsilon;

    #region Trigonmetric functions.

    // Trigonmetric functions.  It is more accurate and efficient to call these
    // rather than first converting to an S1Angle.

    /// <summary>
    /// Returns sin(a)**2, but computed more efficiently.
    /// </summary>
    public double Sin2()
    {
        System.Diagnostics.Debug.Assert(!IsSpecial());
        // Let "a" be the (non-squared) chord length, and let A be the corresponding
        // half-angle (a = 2*sin(A)).  The formula below can be derived from:
        //   sin(2*A) = 2 * sin(A) * cos(A)
        //   cos^2(A) = 1 - sin^2(A)
        // This is much faster than converting to an angle and computing its sine.
        return Length2 * (1 - 0.25 * Length2);
    }

    public double Sin() => Math.Sqrt(Sin2());

    public double Cos()
    {
        // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
        System.Diagnostics.Debug.Assert(!IsSpecial());
        return 1 - 0.5 * Length2;
    }

    public double Tan() => Sin() / Cos();

    /// <summary>
    /// Return the maximum error in length2() for the S1Angle constructor.
    /// </summary>
    /// <remarks>
    /// Assuming that an accurate math library is being used, the sin() call and
    /// the multiply each have a relative error of 0.5 * S2Constants.DoubleEpsilon.  However
    /// the sin() error is squared.
    /// </remarks>
    public double S1AngleConstructorMaxError => 1.5 * S2.DoubleEpsilon * Length2;

    #endregion

    #region Operators

    // Only addition and subtraction of S1ChordAngles is supported.  These
    // methods add or subtract the corresponding S1Angles, and clamp the result
    // to the range [0, Pi].  Both arguments must be non-negative and
    // non-infinite.
    //
    // REQUIRES: !a.is_special() && !b.is_special()

    public static S1ChordAngle operator +(S1ChordAngle a, S1ChordAngle b)
    {
        // Note that this method is much more efficient than converting the chord
        // angles to S1Angles and adding those.  It requires only one square root
        // plus a few additions and multiplications.
        System.Diagnostics.Debug.Assert(!a.IsSpecial());
        System.Diagnostics.Debug.Assert(!b.IsSpecial());

        // Optimization for the common case where "b" is an error tolerance
        // parameter that happens to be set to zero.
        var a2 = a.Length2;
        var b2 = b.Length2;
        if (b2 == 0) return a;

        // Clamp the angle sum to at most 180 degrees.
        if (a2 + b2 >= kMaxLength2) return Straight;

        // Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
        // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
        // Then the formula below can be derived from c = 2 * sin(A+B) and the
        // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
        //                 cos(X) = Math.Sqrt(1 - sin^2(X)) .

        var x = a2 * (1 - 0.25 * b2);  // is_valid() => non-negative
        var y = b2 * (1 - 0.25 * a2);  // is_valid() => non-negative
        return new S1ChordAngle(Math.Min(kMaxLength2, x + y + 2 * Math.Sqrt(x * y)));
    }

    public static S1ChordAngle operator -(S1ChordAngle a, S1ChordAngle b)
    {
        // See comments in operator+().
        System.Diagnostics.Debug.Assert(!a.IsSpecial());
        System.Diagnostics.Debug.Assert(!b.IsSpecial());
        var a2 = a.Length2;
        var b2 = b.Length2;
        if (b2 == 0) return a;
        if (a2 <= b2) return Zero;
        var x = a2 * (1 - 0.25 * b2);
        var y = b2 * (1 - 0.25 * a2);

        // The calculation below is formulated differently (with two square roots
        // rather than one) to avoid excessive cancellation error when two nearly
        // equal values are subtracted.
        var c = Math.Max(0.0, Math.Sqrt(x) - Math.Sqrt(y));
        return new S1ChordAngle(c * c);
    }

    #endregion

    #endregion

    #region IDistance

    public S1ChordAngle GetChordAngleBound() => PlusError(S1AngleConstructorMaxError);

    // If (dist < that), updates that and returns true (used internally).
    public static bool UpdateMin(S1ChordAngle dist, ref S1ChordAngle that)
    {
        if (dist < that)
        {
            that = dist;
            return true;
        }
        return false;
    }

    #endregion

    #region IComparable

    public int CompareTo(S1ChordAngle other)
    {
        if (Length2 < other.Length2) return -1;
        if (Length2 > other.Length2) return 1;
        return 0;
    }

    public static bool operator <(S1ChordAngle x, S1ChordAngle y) => x.Length2 < y.Length2;

    public static bool operator >(S1ChordAngle x, S1ChordAngle y) => x.Length2 > y.Length2;

    public static bool operator <=(S1ChordAngle x, S1ChordAngle y) => x.Length2 <= y.Length2;

    public static bool operator >=(S1ChordAngle x, S1ChordAngle y) => x.Length2 >= y.Length2;

    #endregion

    #region Object

    public override string ToString() => ToAngle().ToString();

    #endregion
}

public interface IDistance<T> where T : IEquatable<T>, IComparable<T>, IDistance<T>
{
    static abstract T operator -(T me, S1ChordAngle other);
    static abstract T operator -(T me, T other);
    static abstract bool operator <(T me, S1ChordAngle other);
    static abstract bool operator >(T me, S1ChordAngle other);
    static abstract bool operator <=(T me, S1ChordAngle other);
    static abstract bool operator >=(T me, S1ChordAngle other);
    static abstract bool operator <(T me, T other);
    static abstract bool operator >(T me, T other);
    static abstract bool operator <=(T me, T other);
    static abstract bool operator >=(T me, T other);
    S1ChordAngle GetChordAngleBound();
    static abstract T Zero { get; }
    static abstract T Infinity { get; }
}
