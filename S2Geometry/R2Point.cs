namespace S2Geometry;

public readonly record struct R2Point(double X, double Y) : IComparable<R2Point>
{
    #region Fields, Constants

    public static readonly R2Point Zero = new(0, 0);

    #endregion

    #region Factories

    /// <summary>
    /// Point as a list of 2; x is index 0, y is index 1
    /// </summary>
    /// <param name="coord"></param>
    public static R2Point FromCoords(IList<double> coord)
    {
        if (coord.Count != 2)
        {
            throw new ArgumentException("Points must have exactly 2 coordinates", nameof(coord));
        }
        return new R2Point(coord[0], coord[1]);
    }

    #endregion

    #region R2Point

    public double this[int index] => index switch
    {
        0 => X,
        1 => Y,
        _ => throw new ArgumentOutOfRangeException(nameof(index)),
    };

    public R2Point SetAxis(int index, double value) => index switch
    {
        0 => new R2Point(value, Y),
        1 => new R2Point(X, value),
        _ => throw new ArgumentOutOfRangeException(nameof(index)),
    };

    /// <summary>
    /// Squared Euclidean norm (the dot product with itself).
    /// </summary>
    /// <returns></returns>
    public double GetNorm2() => DotProd(this);

    /// <summary>
    /// Euclidean norm. For integer T, correct only if Norm2 does not overflow.
    /// </summary>
    /// <returns></returns>
    public double GetNorm() => Math.Sqrt(GetNorm2());

    public double DotProd(R2Point that) => DotProd(this, that);

    public static double DotProd(R2Point p1, R2Point p2) => (p1.X * p2.X) + (p1.Y * p2.Y);

    public static double CrossProd(R2Point p1, R2Point p2) => p1.X * p2.Y - p1.Y * p2.X;

    /// <summary>
    /// Returns a unit vector orthogonal to this one.
    /// </summary>
    public R2Point GetOrtho() => new(-Y, X);

    /// <summary>
    /// return a vector orthogonal to this one
    /// </summary>
    public static R2Point Fabs(R2Point p) => new(Math.Abs(p.X), Math.Abs(p.Y));

    /// <summary>
    /// Normalized vector if the norm is nonzero. Not for integer types.
    /// </summary>
    public static R2Point Normalize(R2Point p)
    {
        var norm = p.GetNorm();
        if (norm != 0.0)
        {
            norm = 1.0 / norm;
        }
        return p * norm;
    }

    /// <summary>
    /// Return the angle between two vectors in radians
    /// </summary>
    public double Angle(R2Point va) => Math.Atan2(CrossProd(this, va), DotProd(va));

    /// <summary>
    /// Compare two vectors, return true if all their components are within a
    /// difference of margin.
    /// </summary>
    public bool ApproxEquals(R2Point that, double margin) => (Math.Abs(X - that.X) < margin) && (Math.Abs(Y - that.Y) < margin);

    /// <summary>
    /// return the index of the smallest, median ,largest component of the vector
    /// </summary>
    public int[] ComponentOrder =>
        new[]
        {
            (0, this[0]),
            (1, this[1]),
            (2, this[2]),
        }
        .OrderBy(t => t.Item2)
        .Select(t => t.Item1)
        .ToArray();

    #endregion

    #region Operators

    public static R2Point operator -(R2Point p1, R2Point p2) => new(p1.X - p2.X, p1.Y - p2.Y);

    public static R2Point operator -(R2Point p) => new(-p.X, -p.Y);

    public static R2Point operator +(R2Point p1, R2Point p2) => new(p1.X + p2.X, p1.Y + p2.Y);

    public static R2Point operator *(R2Point p, double m) => new(m * p.X, m * p.Y);

    public static R2Point operator *(double m, R2Point p) => p * m;

    public static R2Point operator /(R2Point p, double m) => new(p.X / m, p.Y / m);

    #endregion

    #region IComparable

    public int CompareTo(R2Point other)
    {
        if (X < other.X) return -1;
        if (X > other.X) return 1;
        if (Y < other.Y) return -1;
        if (Y > other.Y) return 1;
        return 0;
    }

    public static bool operator <(R2Point x, R2Point y) => x.CompareTo(y) < 0;
    public static bool operator >(R2Point x, R2Point y) => x.CompareTo(y) > 0;
    public static bool operator <=(R2Point x, R2Point y) => x.CompareTo(y) <= 0;
    public static bool operator >=(R2Point x, R2Point y) => x.CompareTo(y) >= 0;

    #endregion

    #region Object

    public override string ToString() => $"[{X:g}, {Y:g}]";

    #endregion
}
