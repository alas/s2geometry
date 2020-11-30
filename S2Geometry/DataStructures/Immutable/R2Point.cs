using System;
using System.Collections.Generic;
using System.Linq;

namespace S2Geometry
{
    public readonly struct R2Point : IEquatable<R2Point>, IComparable<R2Point>
    {
        #region Fields, Constants

        public readonly double X;
        public readonly double Y;

        public static readonly R2Point Zero = new R2Point(0, 0);

        #endregion

        #region Constructors

        public R2Point(double x, double y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        /// Point as a list of 2; x is index 0, y is index 1
        /// </summary>
        /// <param name="coord"></param>
        public R2Point(IList<double> coord)
        {
            if (coord.Count != 2)
            {
                throw new ArgumentException("Points must have exactly 2 coordinates", nameof(coord));
            }
            X = coord[0];
            Y = coord[1];
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
    public double Norm2 => DotProd(this);

        /// <summary>
        /// Euclidean norm. For integer T, correct only if Norm2 does not overflow.
        /// </summary>
        /// <returns></returns>
        public double Norm => Math.Sqrt(Norm2);

        public double DotProd(R2Point that) => DotProd(this, that);

        public static double DotProd(R2Point p1, R2Point p2) => (p1.X * p2.X) + (p1.Y * p2.Y);

        public static double CrossProd(R2Point p1, R2Point p2) => p1.X * p2.Y - p1.Y * p2.X;

        /// <summary>
        /// Returns a unit vector orthogonal to this one.
        /// </summary>
        public R2Point Ortho => new R2Point(-Y, X);

        /// <summary>
        /// return a vector orthogonal to this one
        /// </summary>
        public static R2Point Fabs(R2Point p) => new R2Point(Math.Abs(p.X), Math.Abs(p.Y));

        /// <summary>
        /// Normalized vector if the norm is nonzero. Not for integer types.
        /// </summary>
        public static R2Point Normalize(R2Point p)
        {
            var norm = p.Norm;
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
        public int[] ComponentOrder => new[]
            {
                (0, this[0]),
                (1, this[1]),
                (2, this[2]),
            }
            .OrderBy(t => t.Item2)
            .Select(t => t.Item1)
            .ToArray();

#if s2debug
        public double[] Bounds => new double[2] { X, Y };
#endif

        #endregion

        #region Operators

        public static R2Point operator -(R2Point p1, R2Point p2) => new R2Point(p1.X - p2.X, p1.Y - p2.Y);

        public static R2Point operator -(R2Point p) => new R2Point(-p.X, -p.Y);

        public static R2Point operator +(R2Point p1, R2Point p2) => new R2Point(p1.X + p2.X, p1.Y + p2.Y);

        public static R2Point operator *(R2Point p, double m) => new R2Point(m * p.X, m * p.Y);

        public static R2Point operator *(double m, R2Point p) => p * m;

        public static R2Point operator /(R2Point p, double m) => new R2Point(p.X / m, p.Y / m);

        #endregion

        #region IEquatable

        public bool Equals(R2Point other) => X == other.X && Y == other.Y;

        public override bool Equals(object obj) => obj is R2Point point && Equals(point);

        /// <summary>
        /// Calcualates hashcode based on stored coordinates. Since we want +0.0 and
        /// -0.0 to be treated the same, we ignore the sign of the coordinates.
        /// </summary>
        public override int GetHashCode() => HashCode.Combine(Math.Abs(X), Math.Abs(Y));

        public static bool operator ==(R2Point left, R2Point right) => Equals(left, right);

        public static bool operator !=(R2Point left, R2Point right) => !Equals(left, right);

        #endregion

        #region IComparable

        public static bool operator <(R2Point x, R2Point y)
        {
            if (x.X < y.X)
            {
                return true;
            }
            if (y.X < x.X)
            {
                return false;
            }
            if (x.Y < y.Y)
            {
                return true;
            }
            return false;
        }

        public static bool operator >(R2Point x, R2Point y)
        {
            if (x.X > y.X)
            {
                return true;
            }
            if (y.X > x.X)
            {
                return false;
            }
            if (x.Y > y.Y)
            {
                return true;
            }
            return false;
        }
        public int CompareTo(R2Point other)
        {
            if (X.CompareTo(other.X) != 0)
                return X.CompareTo(other.X);

            return Y.CompareTo(other.Y);
        }

        #endregion

        #region Object

        public override string ToString() => $"[{X:g}, {Y:g}]";
        public string ToString2() => $"({X}, {Y})";

        #endregion   
    }
}
