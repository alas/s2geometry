/// <summary>
/// An S2Point represents a point on the unit sphere as a 3D vector.  Usually
/// points are normalized to be unit length, but some methods do not require
/// this.  See util/math/vector.h for the methods available.  Among other
/// things, there are overloaded operators that make it convenient to write
/// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
/// </summary>
global using S2Point = S2Geometry.Vector3<double>;

namespace S2Geometry;

using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Numerics;

#if s2debug
public static class S2PointDebugExtensions
{
    public static string ToStringDegrees<T>(this Vector3<T> vec) where T : INumber<T>, IFloatingPointIeee754<T>
    {
        if (vec is S2Point p)
        {
            var s2LatLng = new S2LatLng(p);
            return $"({s2LatLng.LatDegrees():g}d, {s2LatLng.LngDegrees():g}d)";
        }

        return vec.ToString();
    }
}
#endif

public static class S2PointExtentions
{
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

    public static Vector3<ExactFloat> ToExact(this S2Point that) => new(
        ExactFloat.FromDouble(that.X),
        ExactFloat.FromDouble(that.Y),
        ExactFloat.FromDouble(that.Z));
}

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
public struct ExactFloat : INumber<ExactFloat>, IFloatingPointIeee754<ExactFloat>
{
    #region Constants

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

    #endregion

    #region Properties

    public required decimal Value { get; init; }

    #endregion

    #region ExactFloat Methods

    public static ExactFloat FromDouble(double p) => new() { Value = (decimal)p };

    public readonly double ToDouble() => decimal.ToDouble(Value);

    /// <summary>
    /// Return true if this value is a normal floating-point number.  Non-normal
    /// values (zero, infinity, and NaN) often need to be handled separately
    /// because they are represented using special exponent values and their
    /// mantissa is not defined.
    /// </summary>
    public readonly bool IsNormal() =>
        Value != 0 && Value != 0.0M /*&& !decimal.IsInfinity(Value) && !decimal.IsNaN(Value)*/;

    /// <summary>
    /// Return the exponent of this ExactFloat given that the mantissa is in the
    /// range [0.5, 1).  It is an error to call this method if the value is zero,
    /// infinity, or NaN.
    /// </summary>
    public readonly int Exp() =>
        1 + (int)decimal.Floor(MathM.Log10(decimal.Abs(Value)));

    #endregion

    #region INumber and IFloatingPointIeee754

    public static ExactFloat One => new() { Value = decimal.One };

    public static int Radix => 10;

    public static ExactFloat Zero => new() { Value = decimal.Zero };

    public static ExactFloat AdditiveIdentity => new() { Value = decimal.Zero };

    public static ExactFloat MultiplicativeIdentity => new() { Value = decimal.One };

    static ExactFloat IFloatingPointIeee754<ExactFloat>.Epsilon => throw new NotImplementedException();

    static ExactFloat IFloatingPointIeee754<ExactFloat>.NaN => throw new NotImplementedException(); //new () { Value = decimal.Zero / ExactFloat.Zero.Value };

    static ExactFloat IFloatingPointIeee754<ExactFloat>.NegativeInfinity => throw new NotImplementedException(); //new() { Value = decimal.MinusOne / decimal.Zero };

    static ExactFloat IFloatingPointIeee754<ExactFloat>.NegativeZero => new() { Value = -0.0M };

    static ExactFloat IFloatingPointIeee754<ExactFloat>.PositiveInfinity => throw new NotImplementedException(); //new() { Value = decimal.One / decimal.Zero };

    static ExactFloat ISignedNumber<ExactFloat>.NegativeOne => new() { Value = decimal.MinusOne };

    static ExactFloat IFloatingPointConstants<ExactFloat>.E => new() { Value = MathM.E };

    static ExactFloat IFloatingPointConstants<ExactFloat>.Pi => new() { Value = MathM.PI };

    static ExactFloat IFloatingPointConstants<ExactFloat>.Tau => new() { Value = MathM.TAU };

    public readonly int CompareTo(object? obj)
    {
        return obj is not ExactFloat ef
            ? throw new NotImplementedException("CompareTo unssuported type: " + obj?.GetType().FullName)
            : CompareTo(ef);
    }

    public readonly int CompareTo(ExactFloat other) => Value.CompareTo(other.Value);

    public static ExactFloat Abs(ExactFloat value) => new() { Value = MathM.Abs(value.Value) };

    public static bool IsCanonical(ExactFloat value) => true;

    public static bool IsComplexNumber(ExactFloat value) => false;

    public static bool IsEvenInteger(ExactFloat value) => IsInteger(value) && (MathM.Abs(value.Value % 2) == 0);

    public static bool IsFinite(ExactFloat value)
    {
        long bits = BitConverter.DoubleToInt64Bits(value.ToDouble());
        return (bits & 0x7FFFFFFFFFFFFFFF) < 0x7FF0000000000000;
    }

    public static bool IsImaginaryNumber(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsInfinity(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsInteger(ExactFloat value) => IsFinite(value) && (value == Truncate(value));
    public static ExactFloat Truncate(ExactFloat x) => new() { Value = MathM.Truncate(x.Value) };

    public static bool IsNaN(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsNegative(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsNegativeInfinity(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsNormal(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsOddInteger(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsPositive(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsPositiveInfinity(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsRealNumber(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsSubnormal(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static bool IsZero(ExactFloat value)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat MaxMagnitude(ExactFloat x, ExactFloat y)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat MaxMagnitudeNumber(ExactFloat x, ExactFloat y)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat MinMagnitude(ExactFloat x, ExactFloat y)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat MinMagnitudeNumber(ExactFloat x, ExactFloat y)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        throw new NotImplementedException();
    }

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, [MaybeNullWhen(false)] out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    public static bool TryParse([NotNullWhen(true)] string? s, NumberStyles style, IFormatProvider? provider, [MaybeNullWhen(false)] out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    public bool Equals(ExactFloat other)
    {
        throw new NotImplementedException();
    }

    public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format, IFormatProvider? provider)
    {
        throw new NotImplementedException();
    }

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat Parse(ReadOnlySpan<char> s, IFormatProvider? provider)
    {
        throw new NotImplementedException();
    }

    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, [MaybeNullWhen(false)] out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    public static ExactFloat Parse(string s, IFormatProvider? provider)
    {
        throw new NotImplementedException();
    }

    public static bool TryParse([NotNullWhen(true)] string? s, IFormatProvider? provider, [MaybeNullWhen(false)] out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    static bool INumberBase<ExactFloat>.TryConvertFromChecked<TOther>(TOther value, out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    static bool INumberBase<ExactFloat>.TryConvertFromSaturating<TOther>(TOther value, out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    static bool INumberBase<ExactFloat>.TryConvertFromTruncating<TOther>(TOther value, out ExactFloat result)
    {
        throw new NotImplementedException();
    }

    static bool INumberBase<ExactFloat>.TryConvertToChecked<TOther>(ExactFloat value, out TOther result)
    {
        throw new NotImplementedException();
    }

    static bool INumberBase<ExactFloat>.TryConvertToSaturating<TOther>(ExactFloat value, out TOther result)
    {
        throw new NotImplementedException();
    }

    static bool INumberBase<ExactFloat>.TryConvertToTruncating<TOther>(ExactFloat value, out TOther result)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPointIeee754<ExactFloat>.Atan2(ExactFloat y, ExactFloat x) => new() { Value = MathM.Atan2(y.Value, x.Value) };

    static ExactFloat IFloatingPointIeee754<ExactFloat>.Atan2Pi(ExactFloat y, ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPointIeee754<ExactFloat>.BitDecrement(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPointIeee754<ExactFloat>.BitIncrement(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPointIeee754<ExactFloat>.FusedMultiplyAdd(ExactFloat left, ExactFloat right, ExactFloat addend)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPointIeee754<ExactFloat>.Ieee754Remainder(ExactFloat left, ExactFloat right)
    {
        throw new NotImplementedException();
    }

    static int IFloatingPointIeee754<ExactFloat>.ILogB(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPointIeee754<ExactFloat>.ScaleB(ExactFloat x, int n)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IExponentialFunctions<ExactFloat>.Exp(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IExponentialFunctions<ExactFloat>.Exp10(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IExponentialFunctions<ExactFloat>.Exp2(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    int IFloatingPoint<ExactFloat>.GetExponentByteCount()
    {
        throw new NotImplementedException();
    }

    int IFloatingPoint<ExactFloat>.GetExponentShortestBitLength()
    {
        throw new NotImplementedException();
    }

    int IFloatingPoint<ExactFloat>.GetSignificandBitLength()
    {
        throw new NotImplementedException();
    }

    int IFloatingPoint<ExactFloat>.GetSignificandByteCount()
    {
        throw new NotImplementedException();
    }

    static ExactFloat IFloatingPoint<ExactFloat>.Round(ExactFloat x, int digits, MidpointRounding mode)
    {
        throw new NotImplementedException();
    }

    bool IFloatingPoint<ExactFloat>.TryWriteExponentBigEndian(Span<byte> destination, out int bytesWritten)
    {
        throw new NotImplementedException();
    }

    bool IFloatingPoint<ExactFloat>.TryWriteExponentLittleEndian(Span<byte> destination, out int bytesWritten)
    {
        throw new NotImplementedException();
    }

    bool IFloatingPoint<ExactFloat>.TryWriteSignificandBigEndian(Span<byte> destination, out int bytesWritten)
    {
        throw new NotImplementedException();
    }

    bool IFloatingPoint<ExactFloat>.TryWriteSignificandLittleEndian(Span<byte> destination, out int bytesWritten)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IHyperbolicFunctions<ExactFloat>.Acosh(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IHyperbolicFunctions<ExactFloat>.Asinh(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IHyperbolicFunctions<ExactFloat>.Atanh(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IHyperbolicFunctions<ExactFloat>.Cosh(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IHyperbolicFunctions<ExactFloat>.Sinh(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IHyperbolicFunctions<ExactFloat>.Tanh(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ILogarithmicFunctions<ExactFloat>.Log(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ILogarithmicFunctions<ExactFloat>.Log(ExactFloat x, ExactFloat newBase)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ILogarithmicFunctions<ExactFloat>.Log10(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ILogarithmicFunctions<ExactFloat>.Log2(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IPowerFunctions<ExactFloat>.Pow(ExactFloat x, ExactFloat y)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IRootFunctions<ExactFloat>.Cbrt(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IRootFunctions<ExactFloat>.Hypot(ExactFloat x, ExactFloat y)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IRootFunctions<ExactFloat>.RootN(ExactFloat x, int n)
    {
        throw new NotImplementedException();
    }

    static ExactFloat IRootFunctions<ExactFloat>.Sqrt(ExactFloat x) => new() { Value = MathM.Sqrt(x.Value) };

    static ExactFloat ITrigonometricFunctions<ExactFloat>.Acos(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.AcosPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.Asin(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.AsinPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.Atan(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.AtanPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.Cos(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.CosPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.Sin(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static (ExactFloat Sin, ExactFloat Cos) ITrigonometricFunctions<ExactFloat>.SinCos(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static (ExactFloat SinPi, ExactFloat CosPi) ITrigonometricFunctions<ExactFloat>.SinCosPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.SinPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.Tan(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    static ExactFloat ITrigonometricFunctions<ExactFloat>.TanPi(ExactFloat x)
    {
        throw new NotImplementedException();
    }

    public static bool operator >(ExactFloat left, ExactFloat right) => left.Value > right.Value;

    public static bool operator >=(ExactFloat left, ExactFloat right) => left.Value >= right.Value;

    public static bool operator <(ExactFloat left, ExactFloat right) => left.Value < right.Value;

    public static bool operator <=(ExactFloat left, ExactFloat right) => left.Value <= right.Value;

    public static ExactFloat operator %(ExactFloat left, ExactFloat right) => new() { Value = left.Value % right.Value };

    public static ExactFloat operator +(ExactFloat left, ExactFloat right) => new() { Value = left.Value + right.Value };

    public static ExactFloat operator --(ExactFloat value) => new() { Value = value.Value - 1 };

    public static ExactFloat operator /(ExactFloat left, ExactFloat right) => new() { Value = left.Value / right.Value };

    public static bool operator ==(ExactFloat left, ExactFloat right) => left.Value == right.Value;

    public static bool operator !=(ExactFloat left, ExactFloat right) => left.Value != right.Value;

    public static ExactFloat operator ++(ExactFloat value) => new() { Value = value.Value + 1 };

    public static ExactFloat operator *(ExactFloat left, ExactFloat right) => new() { Value = left.Value * right.Value };

    public static ExactFloat operator -(ExactFloat left, ExactFloat right) => new() { Value = left.Value - right.Value };

    public static ExactFloat operator -(ExactFloat value) => new() { Value = -value.Value };

    public static ExactFloat operator +(ExactFloat value) => value;

    public static int Sign(ExactFloat value) => MathM.Sign(value.Value);

    public override readonly bool Equals(object? obj) => obj is ExactFloat other && this == other;

    public override readonly int GetHashCode() => Value.GetHashCode();

    #endregion
}
