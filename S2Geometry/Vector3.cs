namespace S2Geometry;

using System.Numerics;

[DebuggerDisplay("{S2PointDebugExtensions.ToStringDegrees(this)}")]
public readonly record struct Vector3<T>(T X, T Y, T Z) : IComparable<Vector3<T>> where T : INumber<T>, IFloatingPointIeee754<T>
{
    #region Fields, Constants

    public static readonly Vector3<T> Empty = new(T.Zero, T.Zero, T.Zero);

    #endregion

    #region Constructors

    public Vector3(IList<T> coords)
        : this(coords[0], coords[1], coords[2])
    {
    }

    public Vector3(T[] array, int offset)
        : this(array[offset], array[offset + 1], array[offset + 2])
    {
    }

    #endregion

    #region Vector3

    public T this[int axis] => axis switch
    {
        0 => X,
        1 => Y,
        2 => Z,
        _ => throw new ArgumentOutOfRangeException(nameof(axis)),
    };

    public Vector3<T> SetAxis(int axis, T value) => axis switch
    {
        0 => new Vector3<T>(value, Y, Z),
        1 => new Vector3<T>(X, value, Z),
        2 => new Vector3<T>(X, Y, value),
        _ => throw new ArgumentOutOfRangeException(nameof(axis)),
    };

    public T DotProd(Vector3<T> that) => X * that.X + Y * that.Y + Z * that.Z;

    public Vector3<T> CrossProd(Vector3<T> p) =>
        new(Y * p.Z - Z * p.Y,
            Z * p.X - X * p.Z,
            X * p.Y - Y * p.X);

    /// <summary>
    /// Returns a unit vector orthogonal to this one.
    /// </summary>
    public Vector3<T> Ortho()
    {
        var temp = new T[][]
        {
            [T.Zero, T.Zero, T.One],
            [T.One,  T.Zero, T.Zero],
            [T.Zero, T.One,  T.Zero],
        };

        return CrossProd(new Vector3<T>(temp[LargestAbsComponent()])).Normalize();
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
    public Vector3<T> Fabs()
    {
        return new(T.Abs(X), T.Abs(Y), T.Abs(Z));
    }

    /// <summary>
    /// Squared Euclidean norm (the dot product with itself).
    /// </summary>
    /// <returns></returns>
    public T Norm2() => DotProd(this);

    /// <summary>
    /// Euclidean norm. For integer T, correct only if Norm2 does not overflow.
    /// </summary>
    /// <returns></returns>
    public T Norm() => T.Sqrt(Norm2());

    /// <summary>
    /// Normalized vector if the norm is nonzero. Not for integer types.
    /// </summary>
    public Vector3<T> Normalize()
    {
        var norm = Norm();
        if (norm != T.Zero)
        {
            norm = T.One / norm;
        }
        return this * norm;
    }

    /// <summary>
    /// Return the angle between two vectors in radians
    /// </summary>
    public T Angle(Vector3<T> va) => T.Atan2(CrossProd(va).Norm(), DotProd(va));

    /// <summary>
    /// return the index of the smallest, median ,largest component of the vector
    /// </summary>
    public int[] ComponentOrder() => [.. new[]
        {
            (0, X),
            (1, Y),
            (2, Z),
        }
        .OrderBy(t => t.Item2)
        .Select(t => t.Item1)];

    /// <summary>
    /// Compare two vectors, return true if all their components are within a
    /// difference of margin.
    /// </summary>
    public bool ApproxEquals(Vector3<T> that, T margin) =>
        T.Abs(X - that.X) < margin &&
        T.Abs(Y - that.Y) < margin &&
        T.Abs(Z - that.Z) < margin;

    public void IntoArray(T[] array, int offset)
    {
        array[offset] = X;
        array[offset + 1] = Y;
        array[offset + 2] = Z;
    }

    public T[] ToArray()
    {
        var array = new T[3];
        array[0] = X;
        array[1] = Y;
        array[2] = Z;
        return array;
    }

    public Vector3<T> MulComponents(Vector3<T> b) =>
        new(X * b.X, Y * b.Y, Z * b.Z);

    public Vector3<T> DivComponents(Vector3<T> b) =>
        new(X / b.X, Y / b.Y, Z / b.Z);

    public Vector3<T> Sqrt() =>
        new(T.Sqrt(X), T.Sqrt(Y), T.Sqrt(Z));

    public Vector3<T> Floor() =>
        new(T.Floor(X), T.Floor(Y), T.Floor(Z));

    public Vector3<T> Ceiling() =>
        new(T.Ceiling(X), T.Ceiling(Y), T.Ceiling(Z));

    public Vector3<T> Round() =>
        new(T.Round(X), T.Round(Y), T.Round(Z));

    #endregion

    #region Long Double

    /// <summary>
    /// No Long Double (16 bytes floating point) support
    /// </summary>
    /// <returns></returns>
    public Vector3<T> ToLD() => this;

    public static Vector3<T> FromLD(Vector3<T> p) => p;

    #endregion

    #region Operators

    public static Vector3<T> operator -(Vector3<T> p1, Vector3<T> p2) => new(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z);

    public static Vector3<T> operator -(Vector3<T> p) => new(-p.X, -p.Y, -p.Z);

    public static Vector3<T> operator +(Vector3<T> p1, Vector3<T> p2) => new(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z);

    public static Vector3<T> operator *(Vector3<T> p, T m) => new(m * p.X, m * p.Y, m * p.Z);

    public static Vector3<T> operator *(T m, Vector3<T> p) => new(m * p.X, m * p.Y, m * p.Z);

    public static Vector3<T> operator /(Vector3<T> p, T m) => new(p.X / m, p.Y / m, p.Z / m);

    #endregion

    #region IComparable

    public int CompareTo(Vector3<T> other)
    {
        if (X < other.X) return -1;
        if (X > other.X) return 1;
        if (Y < other.Y) return -1;
        if (Y > other.Y) return 1;
        if (Z < other.Z) return -1;
        if (Z > other.Z) return 1;
        return 0;
    }

    public static bool operator <(Vector3<T> x, Vector3<T> y) => x.CompareTo(y) < 0;
    public static bool operator >(Vector3<T> x, Vector3<T> y) => x.CompareTo(y) > 0;
    public static bool operator <=(Vector3<T> x, Vector3<T> y) => x.CompareTo(y) <= 0;
    public static bool operator >=(Vector3<T> x, Vector3<T> y) => x.CompareTo(y) >= 0;

    #endregion

    #region Object

    public override string ToString() => $"[{X:g}, {Y:g}, {Z:g}]";

    #endregion
}
