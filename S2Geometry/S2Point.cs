namespace S2Geometry;

/// <summary>
/// An S2Point represents a point on the unit sphere as a 3D vector.  Usually
/// points are normalized to be unit length, but some methods do not require
/// this.  See util/math/vector.h for the methods available.  Among other
/// things, there are overloaded operators that make it convenient to write
/// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
/// </summary>
[System.Diagnostics.DebuggerDisplay("{ToStringDegrees()}")]
public readonly record struct S2Point(double X, double Y, double Z) : IComparable<S2Point>
{
    // todo(Alas): add support for more floating point types and make S2Point generic
    // i.e.: S2Point<T> where T in (float, double, long double (16 bytes floating point), exact float, etc)

    #region Fields, Constants

    public static readonly S2Point Empty = new(0, 0, 0);

    #endregion

    #region Constructors

    public S2Point(IList<double> coords)
        : this(coords[0], coords[1], coords[2])
    {
    }

    public S2Point(double[] array, int offset)
        : this(array[offset], array[offset + 1], array[offset + 2])
    {
    }

    #endregion

    #region S2Point

    public double this[int axis] => axis switch
    {
        0 => X,
        1 => Y,
        2 => Z,
        _ => throw new ArgumentOutOfRangeException(nameof(axis)),
    };

    public S2Point SetAxis(int axis, double value) => axis switch
    {
        0 => new S2Point(value, Y, Z),
        1 => new S2Point(X, value, Z),
        2 => new S2Point(X, Y, value),
        _ => throw new ArgumentOutOfRangeException(nameof(axis)),
    };

    public double DotProd(S2Point that) => X * that.X + Y * that.Y + Z * that.Z;

    public S2Point CrossProd(S2Point p) => new(
            Y * p.Z - Z * p.Y,
            Z * p.X - X * p.Z,
            X * p.Y - Y * p.X);

    /// <summary>
    /// Returns a unit vector orthogonal to this one.
    /// </summary>
    public S2Point Ortho()
    {
        var temp = new[]
        {
                new[] { 0.0, 0.0, 1.0 },
                new[] { 1.0, 0.0, 0.0 },
                new[] { 0.0, 1.0, 0.0 },
            };

        return CrossProd(new S2Point(temp[LargestAbsComponent()])).Normalize();
    }

    /// <summary>
    /// return the index of the largest component (fabs)
    /// </summary>
    public int LargestAbsComponent()
    {
        var temp = Fabs();
        if (temp.X > temp.Y)
        {
            if (temp.X > temp.Z)
            {
                return 0;
            }
            else
            {
                return 2;
            }
        }
        else
        {
            if (temp.Y > temp.Z)
            {
                return 1;
            }
            else
            {
                return 2;
            }
        }
    }

    /// <summary>
    /// return a vector orthogonal to this one
    /// </summary>
    public S2Point Fabs()
    {
        return new(Math.Abs(X), Math.Abs(Y), Math.Abs(Z));
    }

    /// <summary>
    /// Squared Euclidean norm (the dot product with itself).
    /// </summary>
    /// <returns></returns>
    public double Norm2() => DotProd(this);

    /// <summary>
    /// Euclidean norm. For integer T, correct only if Norm2 does not overflow.
    /// </summary>
    /// <returns></returns>
    public double Norm() => Math.Sqrt(Norm2());

    /// <summary>
    /// Normalized vector if the norm is nonzero. Not for integer types.
    /// </summary>
    public S2Point Normalize()
    {
        var norm = Norm();
        if (norm != 0.0)
        {
            norm = 1.0 / norm;
        }
        return this * norm;
    }

    /// <summary>
    /// Return the angle between two vectors in radians
    /// </summary>
    public double Angle(S2Point va) => Math.Atan2(CrossProd(va).Norm(), DotProd(va));

    /// <summary>
    /// return the index of the smallest, median ,largest component of the vector
    /// </summary>
    public int[] ComponentOrder() => new[]
        {
            (0, X),
            (1, Y),
            (2, Z),
        }
        .OrderBy(t => t.Item2)
        .Select(t => t.Item1)
        .ToArray();

    /// <summary>
    /// Return true if the given point is approximately unit length
    /// (this is mainly useful for assertions).
    /// 
    /// Normalize() is guaranteed to return a vector whose L2-norm differs from 1
    /// by less than 2 * S2Constants.DoubleEpsilon.  Thus the squared L2-norm differs by less
    /// than 4 * S2Constants.DoubleEpsilon.  The actual calculated Norm2() can have up to 1.5 *
    /// S2Constants.DoubleEpsilon of additional error.  The total error of 5.5 * S2Constants.DoubleEpsilon
    /// can then be rounded down since the result must be a representable
    /// double-precision value.
    /// </summary>
    public bool IsUnitLength()
    {
        var na = Norm2() - 1;
        var abs = Math.Abs(na);

        const double err = 5 * S2.DoubleEpsilon;  // About 1.11e-15
        var isLess = abs <= err;
        return isLess;
    }

    /// <summary>
    /// Compare two vectors, return true if all their components are within a
    /// difference of margin.
    /// </summary>
    public bool ApproxEquals(S2Point that, double margin) =>
            Math.Abs(X - that.X) < margin &&
            Math.Abs(Y - that.Y) < margin &&
            Math.Abs(Z - that.Z) < margin;

    public void IntoArray(double[] array, int offset)
    {
        array[offset] = X;
        array[offset + 1] = Y;
        array[offset + 2] = Z;
    }

    public double[] ToArray()
    {
        var array = new double[3];
        array[0] = X;
        array[1] = Y;
        array[2] = Z;
        return array;
    }

    #endregion

    #region ExactFloat

    /// <summary>
    /// No ExactFloat support
    /// </summary>
    public S2Point ToExact() => this;

    public static S2Point FromExact(S2Point p) => p;

    #endregion

    #region Long Double

    /// <summary>
    /// No Long Double (16 bytes floating point) support
    /// </summary>
    /// <returns></returns>
    public S2Point ToLD() => this;

    public static S2Point FromLD(S2Point p) => p;

    #endregion

    #region Operators

    public static S2Point operator -(S2Point p1, S2Point p2) => new(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z);

    public static S2Point operator -(S2Point p) => new(-p.X, -p.Y, -p.Z);

    public static S2Point operator +(S2Point p1, S2Point p2) => new(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z);

    public static S2Point operator *(S2Point p, double m) => new(m * p.X, m * p.Y, m * p.Z);

    public static S2Point operator *(double m, S2Point p) => new(m * p.X, m * p.Y, m * p.Z);

    public static S2Point operator /(S2Point p, double m) => new(p.X / m, p.Y / m, p.Z / m);

    #endregion

    #region IComparable

    public int CompareTo(S2Point other)
    {
        if (X < other.X) return -1;
        if (X > other.X) return 1;
        if (Y < other.Y) return -1;
        if (Y > other.Y) return 1;
        if (Z < other.Z) return -1;
        if (Z > other.Z) return 1;
        return 0;
    }

    public static bool operator <(S2Point x, S2Point y) => x.CompareTo(y) < 0;
    public static bool operator >(S2Point x, S2Point y) => x.CompareTo(y) > 0;
    public static bool operator <=(S2Point x, S2Point y) => x.CompareTo(y) <= 0;
    public static bool operator >=(S2Point x, S2Point y) => x.CompareTo(y) >= 0;

    #endregion

    #region Object

    public override string ToString() => $"[{X:g}, {Y:g}, {Z:g}]";

#if s2debug
    public string ToStringDegrees()
    {
        var s2LatLng = new S2LatLng(this);
        return $"({s2LatLng.LatDegrees():g}d, {s2LatLng.LngDegrees():g}d)";
    }
#endif

    #endregion
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
