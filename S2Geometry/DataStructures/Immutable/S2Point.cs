using System;
using System.Linq;

namespace S2Geometry
{
    /// <summary>
    /// An S2Point represents a point on the unit sphere as a 3D vector.  Usually
    /// points are normalized to be unit length, but some methods do not require
    /// this.  See util/math/vector.h for the methods available.  Among other
    /// things, there are overloaded operators that make it convenient to write
    /// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
    /// </summary>
    [System.Diagnostics.DebuggerDisplay("{ToStringDegrees()}")]
	public readonly struct S2Point : IEquatable<S2Point>, IComparable<S2Point>
    {
        // todo(Alas): add support for more floating point types and make S2Point generic
        // i.e.: S2Point<T> where T in (float, double, long double (16 bytes floating point), exact float, etc)
        // maybe move S2Point and S2EdgeCrossings to an assembly in F#

        #region Fields, Constants

        public readonly double X;
        public readonly double Y;
        public readonly double Z;

        public static readonly S2Point Empty = new S2Point(0, 0, 0);

        #endregion

        #region Constructors

        public S2Point(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public S2Point(double[] coords)
        {
            X = coords[0];
            Y = coords[1];
            Z = coords[2];
        }

        public S2Point(double[] array, int offset)
        {
            X = array[offset];
            Y = array[offset + 1];
            Z = array[offset + 2];
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

        public S2Point CrossProd(S2Point p) => new S2Point(
                Y * p.Z - Z * p.Y,
                Z * p.X - X * p.Z,
                X * p.Y - Y * p.X);

        /// <summary>
        /// Returns a unit vector orthogonal to this one.
        /// </summary>
        public S2Point Ortho
        {
            get
            {
                var temp = new[]
                {
                    new[] { 0.0, 0.0, 1.0 },
                    new[] { 1.0, 0.0, 0.0 },
                    new[] { 0.0, 1.0, 0.0 },
                };
                return CrossProd(new S2Point(temp[LargestAbsComponent])).Normalized;
            }
        }

        /// <summary>
        /// return the index of the largest component (fabs)
        /// </summary>
        public int LargestAbsComponent
        {
            get
            {
                var temp = Fabs;
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
        }

        /// <summary>
        /// return a vector orthogonal to this one
        /// </summary>
        public S2Point Fabs => new S2Point(Math.Abs(X), Math.Abs(Y), Math.Abs(Z));

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

        /// <summary>
        /// Normalized vector if the norm is nonzero. Not for integer types.
        /// </summary>
        public S2Point Normalized
        {
            get
            {
                var norm = Norm;
                if (norm != 0.0)
                {
                    norm = 1.0 / norm;
                }
                return this * norm;
            }
        }

        /// <summary>
        /// Return the angle between two vectors in radians
        /// </summary>
        public double Angle(S2Point va) => Math.Atan2(CrossProd(va).Norm, DotProd(va));

        /// <summary>
        /// return the index of the smallest, median ,largest component of the vector
        /// </summary>
        public int[] ComponentOrder => new[]
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
        public bool IsUnitLength
        {
            get
            {
                var na = Norm2 - 1;
                var abs = Math.Abs(na);

                const double err = 5 * S2Constants.DoubleEpsilon;  // About 1.11e-15
                var isLess = abs <= err;
                return isLess;
            }
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

        #endregion

        #region ExactFloat

        /// <summary>
        /// No ExactFloat support
        /// </summary>
        internal S2Point ToExact() => this;

        internal static S2Point FromExact(S2Point p) => p;

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

        public static S2Point operator -(S2Point p1, S2Point p2) => new S2Point(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z);

        public static S2Point operator -(S2Point p) => new S2Point(-p.X, -p.Y, -p.Z);

        public static S2Point operator +(S2Point p1, S2Point p2) => new S2Point(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z);

        public static S2Point operator *(S2Point p, double m) => new S2Point(m * p.X, m * p.Y, m * p.Z);

        public static S2Point operator *(double m, S2Point p) => new S2Point(m * p.X, m * p.Y, m * p.Z);

        public static S2Point operator /(S2Point p, double m) => new S2Point(p.X / m, p.Y / m, p.Z / m);

        #endregion

        #region IEquatable

        /// <summary>
        /// Calcualates hashcode based on stored coordinates. Since we want +0.0 and
        /// -0.0 to be treated the same, we ignore the sign of the coordinates.
        /// </summary>
        public override int GetHashCode() => HashCode.Combine(Math.Abs(X), Math.Abs(Y), Math.Abs(Z));

        public bool Equals(S2Point other) => X == other.X && Y == other.Y && Z == other.Z;

        public override bool Equals(object obj) => obj is S2Point point && Equals(point);

        public static bool operator ==(S2Point left, S2Point right) => Equals(left, right);

        public static bool operator !=(S2Point left, S2Point right) => !Equals(left, right);

        #endregion

        #region IComparable

        public int CompareTo(S2Point other) => this < other ? -1 : (Equals(other) ? 0 : 1);

        public static bool operator <(S2Point x, S2Point y)
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
            if (y.Y < x.Y)
            {
                return false;
            }
            if (x.Z < y.Z)
            {
                return true;
            }
            return false;
        }

        public static bool operator >(S2Point x, S2Point y)
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
            if (y.Y > x.Y)
            {
                return false;
            }
            if (x.Z > y.Z)
            {
                return true;
            }
            return false;
        }

        public static bool operator <=(S2Point x, S2Point y) => x < y || x == y;

        public static bool operator >=(S2Point x, S2Point y) => x > y || x == y;

        #endregion

        #region Object

        public override string ToString() => $"[{X:g}, {Y:g}, {Z:g}]";

        public string ToStringDegrees()
        {
            var s2LatLng = new S2LatLng(this);
            return $"({s2LatLng.LatDegrees:g}d, {s2LatLng.LngDegrees:g}d)";
        }

        #endregion
    }

    /// <summary>
    /// todo
    /// </summary>
    public static class ExactFloat
    {
        public static double FromDouble(double p) => p;

        public static double ToDouble(this double that) => that;
    }
}
