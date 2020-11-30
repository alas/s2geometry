using System;

namespace S2Geometry
{
    [System.Diagnostics.DebuggerDisplay("{ToDebugString()}")]
    public readonly struct R1Interval : IEquatable<R1Interval>
    {
        #region Fields, Constants

        public readonly double Lo;

        public readonly double Hi;

        /// <summary>
        ///     Returns an empty interval. (Any interval where lo > hi is considered empty.)
        /// </summary>
        public static readonly R1Interval Empty = new R1Interval(1, 0);

        #endregion

        #region Constructors

        /// <summary>
        /// Constructor.  If lo > hi, the interval is empty.
        /// </summary>
        public R1Interval(double[] lo_hi)
        {
            Lo = lo_hi[0];
            Hi = lo_hi[1];
        }

        /// <summary>
        /// Constructor.  If lo > hi, the interval is empty.
        /// </summary>
        public R1Interval(double lo, double hi)
        {
            Lo = lo;
            Hi = hi;
        }

        #endregion

        #region Factories

        /// <summary>
        /// Convenience method to construct an interval containing a single point.
        /// </summary>
        public static R1Interval FromPoint(double p) => new R1Interval(p, p);

        /// <summary>
        /// Convenience method to construct the minimal interval containing the two
        /// given points. This is equivalent to starting with an empty interval and
        /// calling AddPoint() twice, but it is more efficient.
        /// </summary>
        public static R1Interval FromPointPair(double p1, double p2)
        {
            if (p1 <= p2)
            {
                return new R1Interval(p1, p2);
            }
            else
            {
                return new R1Interval(p2, p1);
            }
        }

        #endregion

        #region R1Interval

        public double this[int index]
        {
            get
            {
                return index switch
                {
                    0 => Lo,
                    1 => Hi,
                    _ => throw new ArgumentOutOfRangeException(nameof(index)),
                };
            }
        }

        /// <summary>
        /// Return the center of the interval.  For empty intervals,
        /// the result is arbitrary.
        /// </summary>
        public double Center => 0.5 * (Lo + Hi);

        /// <summary>
        /// Return the length of the interval.  The length of an empty interval
        /// is negative.
        /// </summary>
        public double Length => Hi - Lo;

        /// <summary>
        /// Return true if the interval is empty, i.e. it contains no points.
        /// </summary>
        public bool IsEmpty => Lo > Hi;

        public bool Contains(double p) => p >= Lo && p <= Hi;

        public bool InteriorContains(double p) => p > Lo && p < Hi;

        /// <summary>
        /// Return true if this interval contains the interval 'y'.
        /// </summary>
        public bool Contains(R1Interval y) => y.IsEmpty || (y.Lo >= Lo && y.Hi <= Hi);

        /// <summary>
        /// Return true if the interior of this interval contains the entire
        /// interval 'y' (including its boundary).
        /// </summary>
        public bool InteriorContains(R1Interval y) => y.IsEmpty || (y.Lo > Lo && y.Hi < Hi);

        /// <summary>
        /// Return true if this interval intersects the given interval,
        /// i.e. if they have any points in common.
        /// </summary>
        public bool Intersects(R1Interval y)
        {
            if (Lo <= y.Lo)
            {
                return y.Lo <= Hi && y.Lo <= y.Hi;
            }
            else
            {
                return Lo <= y.Hi && Lo <= Hi;
            }
        }

        /// <summary>
        /// Return true if the interior of this interval intersects
        /// any point of the given interval (including its boundary).
        /// </summary>
        public bool InteriorIntersects(R1Interval y) => y.Lo < Hi && Lo < y.Hi && Lo < Hi && y.Lo <= y.Hi;

        /// <summary>
        /// Return the Hausdorff distance to the given interval 'y'. For two
        /// R1Intervals x and y, this distance is defined as
        ///     h(x, y) = max_{p in x} min_{q in y} d(p, q).
        /// </summary>
        public double GetDirectedHausdorffDistance(R1Interval y)
        {
            if (IsEmpty) return 0.0;

            if (y.IsEmpty) return double.MaxValue;

            return Math.Max(0.0, Math.Max(Hi - y.Hi, y.Lo - Lo));
        }

        /// <summary>
        /// Expand the interval so that it contains the given point "p".
        /// </summary>
        public static R1Interval AddPoint(R1Interval i, double p)
        {
            if (i.IsEmpty)
            {
                return FromPoint(p);
            }
            else if (p < i.Lo)
            {
                return new R1Interval(p, i.Hi);
            }
            else if (p > i.Hi)
            {
                return new R1Interval(i.Lo, p);
            }
            else
            {
                return i;
            }
        }

        /// <summary>
        /// Expand the interval so that it contains the given interval "y".
        /// </summary>
        public static R1Interval AddInterval(R1Interval x, R1Interval y)
        {
            if (y.IsEmpty) return x;
            if (x.IsEmpty) return y;

            return new R1Interval(
                (y.Lo < x.Lo) ? y.Lo : x.Lo,
                (y.Hi > x.Hi) ? y.Hi : x.Hi);
        }

        /// <summary>
        /// Return the closest point in the interval to the given point "p".
        /// The interval must be non-empty.
        /// </summary>
        public double Project(double p)
        {
            Assert.True(!IsEmpty);
            return Math.Max(Lo, Math.Min(Hi, p));
        }

        /// <summary>
        /// Return an interval that has been expanded on each side by the given
        /// distance "margin".  If "margin" is negative, then shrink the interval on
        /// each side by "margin" instead.  The resulting interval may be empty.  Any
        /// expansion of an empty interval remains empty.
        /// </summary>
        public R1Interval Expanded(double margin)
        {
            if (IsEmpty)
            {
                return this;
            }
            return new R1Interval(Lo - margin, Hi + margin);
        }

        /// <summary>
        /// Return the smallest interval that contains this interval and the
        /// given interval "y".
        /// </summary>
        public R1Interval Union(R1Interval y)
        {
            if (IsEmpty)
            {
                return y;
            }
            if (y.IsEmpty)
            {
                return this;
            }
            return new R1Interval(Math.Min(Lo, y.Lo), Math.Max(Hi, y.Hi));
        }

        /// <summary>
        /// Return the intersection of this interval with the given interval.
        /// Empty intervals do not need to be special-cased.
        /// </summary>
        public R1Interval Intersection(R1Interval y) => new R1Interval(Math.Max(Lo, y.Lo), Math.Min(Hi, y.Hi));

        /// <summary>
        /// Return true if this interval can be transformed into the given interval
        /// by moving each endpoint by at most "max_error".  The empty interval is
        /// considered to be positioned arbitrarily on the real line, thus any
        /// interval with (length <= 2*max_error) matches the empty interval.
        /// </summary>
        public bool ApproxEquals(R1Interval y, double maxError = S2Constants.DoubleError)
        {
            if (IsEmpty)
            {
                return y.Length <= 2 * maxError;
            }
            if (y.IsEmpty)
            {
                return Length <= 2 * maxError;
            }
            return Math.Abs(y.Lo - Lo) <= maxError
                && Math.Abs(y.Hi - Hi) <= maxError;
        }

        #endregion

        #region IEquatable

        public bool Equals(R1Interval other) => (Hi == other.Hi && Lo == other.Lo) || (IsEmpty && other.IsEmpty);

        public override bool Equals(object obj) => obj is R1Interval interval && Equals(interval);

        public override int GetHashCode() => HashCode.Combine(Hi, Lo);

        /// <summary>
        /// Return true if two intervals contain the same set of points.
        /// </summary>
        public static bool operator ==(R1Interval left, R1Interval right) => Equals(left, right);

        /// <summary>
        /// Return true if two intervals do not contain the same set of points.
        /// </summary>
        public static bool operator !=(R1Interval left, R1Interval right) => !Equals(left, right);

        #endregion

        #region Object

        public override string ToString() => $"[{Lo}, {Hi}]";
        public string ToDebugString() => IsEmpty ? "[Empty]" : $"{S1Angle.FromRadians(Lo)}, {S1Angle.FromRadians(Hi)}";

        #endregion
    }
}
