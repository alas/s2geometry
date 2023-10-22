// An S1Interval represents a closed interval on a unit circle (also known
// as a 1-dimensional sphere).  It is capable of representing the empty
// interval (containing no points), the full interval (containing all
// points), and zero-length intervals (containing a single point).
//
// Points are represented by the angle they make with the positive x-axis in
// the range [-Pi, Pi].  An interval is represented by its lower and upper
// bounds (both inclusive, since the interval is closed).  The lower bound may
// be greater than the upper bound, in which case the interval is "inverted"
// (i.e. it passes through the point (-1, 0)).
//
// Note that the point (-1, 0) has two valid representations, Pi and -Pi.
// The normalized representation of this point internally is Pi, so that
// endpoints of normal intervals are in the range (-Pi, Pi].  However, we
// take advantage of the point -Pi to construct two special intervals:
// the Full() interval is [-Pi, Pi], and the Empty() interval is [Pi, -Pi].
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.

namespace S2Geometry;

[DebuggerDisplay("{ToDebugString()}")]
public readonly record struct S1Interval
{
    #region Fields, Constants

    public readonly double Lo { get; init; }
    public readonly double Hi { get; init; }

    public static readonly S1Interval Empty = new(Math.PI, -Math.PI);
    public static readonly S1Interval Full = new(-Math.PI, Math.PI);

    #endregion

    #region Constructors

    // Constructor.  Both endpoints must be in the range -Pi to Pi inclusive.
    // The value -Pi is converted internally to Pi except for the Full()
    // and Empty() intervals.
    public S1Interval(double lo, double hi)
    {
        var lo_ = lo;
        var hi_ = hi;
        if (lo == -Math.PI && hi != Math.PI) lo_ = Math.PI;
        if (hi == -Math.PI && lo != Math.PI) hi_ = Math.PI;
        Lo = lo_;
        Hi = hi_;
        MyDebug.Assert(IsValid());
    }

    #endregion

    #region Factories

    /// <summary>
    /// Convenience method to construct an interval containing a single point.
    /// </summary>
    public static S1Interval FromPoint(double p)
    {
        if (p == -Math.PI) p = Math.PI;
        return new S1Interval(p, p);
    }

    /// <summary>
    /// Convenience method to construct the minimal interval containing
    /// the two given points.This is equivalent to starting with an empty
    /// interval and calling AddPoint() twice, but it is more efficient.
    /// </summary>
    public static S1Interval FromPointPair(double p1, double p2)
    {
        MyDebug.Assert(Math.Abs(p1) <= Math.PI);
        MyDebug.Assert(Math.Abs(p2) <= Math.PI);
        if (p1 == -Math.PI) p1 = Math.PI;
        if (p2 == -Math.PI) p2 = Math.PI;
        if (PositiveDistance(p1, p2) <= Math.PI)
        {
            return new S1Interval(p1, p2);
        }
        else
        {
            return new S1Interval(p2, p1);
        }
    }

    #endregion

    #region S1Interval

    public double this[int index] => index switch
    {
        0 => Lo,
        1 => Hi,
        _ => throw new ArgumentOutOfRangeException(nameof(index)),
    };

    // An interval is valid if neither bound exceeds Pi in absolute value,
    // and the value -Pi appears only in the Empty() and Full() intervals.
    public bool IsValid() => Math.Abs(Lo) <= Math.PI && Math.Abs(Hi) <= Math.PI &&
        !(Lo == -Math.PI && Hi != Math.PI) &&
        !(Hi == -Math.PI && Lo != Math.PI);

    // Return true if the interval contains all points on the unit circle.
    public bool IsFull() => Lo == -Math.PI && Hi == Math.PI;

    // Return true if the interval is empty, i.e. it contains no points.
    public bool IsEmpty() => Lo == Math.PI && Hi == -Math.PI;

    // Return true if Lo > Hi.  (This is true for empty intervals.)
    public bool IsInverted() => Lo > Hi;

    // Return the midpoint of the interval.  For full and empty intervals,
    // the result is arbitrary.
    public double GetCenter()
    {
        var center = 0.5 * (Lo + Hi);
        if (!IsInverted()) return center;
        // Return the center in the range (-Pi, Pi].
        return (center <= 0) ? (center + Math.PI) : (center - Math.PI);
    }

    // Return the length of the interval.  The length of an empty interval
    // is negative.
    public double GetLength()
    {
        var length = Hi - Lo;
        if (length >= 0) return length;
        length += S2.M_2_PI;
        // Empty intervals have a negative length.
        return (length > 0) ? length : -1;
    }

    // Return the complement of the interior of the interval.  An interval and
    // its complement have the same boundary but do not share any interior
    // values.  The complement operator is not a bijection, since the complement
    // of a singleton interval (containing a single value) is the same as the
    // complement of an empty interval.
    public S1Interval Complement()
    {
        if (Lo == Hi) return Full;   // Singleton.

        // The complement of(lo, hi), is (hi, lo)
        return new S1Interval(Hi, Lo);  // Handles empty and full.
    }

    // Return the midpoint of the complement of the interval. For full and empty
    // intervals, the result is arbitrary. For a singleton interval (containing a
    // single point), the result is its antipodal point on S1.
    public double ComplementCenter()
    {
        if (Lo != Hi)
        {
            return Complement().GetCenter();
        }
        else
        {  // Singleton.
            return (Hi <= 0) ? (Hi + Math.PI) : (Hi - Math.PI);
        }
    }

    // Return true if the interval (which is closed) contains the point 'p'.
    public bool Contains(double p)
    {
        // Works for empty, full, and singleton intervals.
        MyDebug.Assert(Math.Abs(p) <= Math.PI);
        if (p == -Math.PI) p = Math.PI;
        return FastContains(p);
    }

    // Return true if the interior of the interval contains the point 'p'.
    public bool InteriorContains(double p)
    {
        // Works for empty, full, and singleton intervals.
        MyDebug.Assert(Math.Abs(p) <= Math.PI);
        if (p == -Math.PI) p = Math.PI;

        if (IsInverted())
        {
            return p > Lo || p < Hi;
        }
        else
        {
            return (p > Lo && p < Hi) || IsFull();
        }
    }

    // Return true if the interval contains the given interval 'y'.
    // Works for empty, full, and singleton intervals.
    public bool Contains(S1Interval y)
    {
        // It might be helpful to compare the structure of these tests to
        // the simpler Contains(double) method above.

        if (IsInverted())
        {
            if (y.IsInverted()) return y.Lo >= Lo && y.Hi <= Hi;
            return (y.Lo >= Lo || y.Hi <= Hi) && !IsEmpty();
        }
        else
        {
            if (y.IsInverted()) return IsFull() || y.IsEmpty();
            return y.Lo >= Lo && y.Hi <= Hi;
        }
    }

    // Returns true if the interior of this interval contains the entire
    // interval 'y'.  Note that x.InteriorContains(x) is true only when
    // x is the empty or full interval, and x.InteriorContains(S1Interval(p,p))
    // is equivalent to x.InteriorContains(p).
    public bool InteriorContains(S1Interval y)
    {
        if (IsInverted())
        {
            if (!y.IsInverted()) return y.Lo > Lo || y.Hi < Hi;
            return (y.Lo > Lo && y.Hi < Hi) || y.IsEmpty();
        }
        else
        {
            if (y.IsInverted()) return IsFull() || y.IsEmpty();
            return (y.Lo > Lo && y.Hi < Hi) || IsFull();
        }
    }

    // Return true if the two intervals contain any points in common.
    // Note that the point +/-Pi has two representations, so the intervals
    // [-Pi,-3] and [2,Pi] intersect, for example.
    public bool Intersects(S1Interval y)
    {
        if (IsEmpty() || y.IsEmpty()) return false;
        if (IsInverted())
        {
            // Every non-empty inverted interval contains Pi.
            return y.IsInverted() || y.Lo <= Hi || y.Hi >= Lo;
        }
        else
        {
            if (y.IsInverted()) return y.Lo <= Hi || y.Hi >= Lo;
            return y.Lo <= Hi && y.Hi >= Lo;
        }
    }

    // Return true if the interior of this interval contains any point of the
    // interval 'y' (including its boundary).  Works for empty, full, and
    // singleton intervals.
    public bool InteriorIntersects(S1Interval y)
    {
        if (IsEmpty() || y.IsEmpty() || Lo == Hi) return false;
        if (IsInverted())
        {
            return y.IsInverted() || y.Lo < Hi || y.Hi > Lo;
        }
        else
        {
            if (y.IsInverted()) return y.Lo < Hi || y.Hi > Lo;
            return (y.Lo < Hi && y.Hi > Lo) || IsFull();
        }
    }

    // Return the Hausdorff distance to the given interval 'y'. For two
    // S1Intervals x and y, this distance is defined by
    //     h(x, y) = max_{p in x} min_{q in y} d(p, q),
    // where d(.,.) is measured along S1.
    public double GetDirectedHausdorffDistance(S1Interval y)
    {
        if (y.Contains(this)) return 0.0;  // this includes the case *this is empty
        if (y.IsEmpty()) return Math.PI;  // maximum possible distance on S1

        var y_complement_center = y.ComplementCenter();
        if (Contains(y_complement_center))
        {
            return PositiveDistance(y.Hi, y_complement_center);
        }
        else
        {
            // The Hausdorff distance is realized by either two Hi endpoints or two
            // Lo endpoints, whichever is farther apart.
            var hi_hi = new S1Interval(y.Hi, y_complement_center).Contains(Hi)
                ? PositiveDistance(y.Hi, Hi) : 0;
            var lo_lo = new S1Interval(y_complement_center, y.Lo).Contains(Lo)
                ? PositiveDistance(Lo, y.Lo) : 0;
            MyDebug.Assert(hi_hi > 0 || lo_lo > 0);
            return Math.Max(hi_hi, lo_lo);
        }
    }

    public static S1Interval AddPoints(S1Interval interval, double p1, double p2)
        => AddPoint(AddPoint(interval, p1), p2);

    // Expand the interval by the minimum amount necessary so that it
    // contains the given point "p" (an angle in the range [-Pi, Pi]).
    public static S1Interval AddPoint(S1Interval interval, double p)
    {
        MyDebug.Assert(Math.Abs(p) <= Math.PI);
        if (p == -Math.PI) p = Math.PI;

        if (!interval.FastContains(p))
        {
            if (interval.IsEmpty())
            {
                return new S1Interval(p, p);
            }
            else
            {
                // Compute distance from p to each endpoint.
                var dlo = PositiveDistance(p, interval.Lo);
                var dhi = PositiveDistance(interval.Hi, p);
                if (dlo < dhi)
                {
                    return new S1Interval(p, interval.Hi);
                }
                else
                {
                    return new S1Interval(interval.Lo, p);
                }
                // Adding a point can never turn a non-full interval into a full one.
            }
        }
        return interval;
    }

    // Return the closest point in the interval to the given point "p".
    // The interval must be non-empty.
    public double Project(double p)
    {
        MyDebug.Assert(!IsEmpty());
        MyDebug.Assert(Math.Abs(p) <= Math.PI);
        if (p == -Math.PI) p = Math.PI;
        if (FastContains(p)) return p;
        // Compute distance from p to each endpoint.
        var dlo = PositiveDistance(p, Lo);
        var dhi = PositiveDistance(Hi, p);
        return (dlo < dhi) ? Lo : Hi;
    }

    // Return an interval that has been expanded on each side by the given
    // distance "margin".  If "margin" is negative, then shrink the interval on
    // each side by "margin" instead.  The resulting interval may be empty or
    // full.  Any expansion (positive or negative) of a full interval remains
    // full, and any expansion of an empty interval remains empty.
    public S1Interval Expanded(double margin)
    {
        if (margin >= 0)
        {
            if (IsEmpty()) return this;
            // Check whether this interval will be full after expansion, allowing
            // for a 1-bit rounding error when computing each endpoint.
            if (GetLength() + 2 * margin + 2 * S2.DoubleEpsilon >= S2.M_2_PI) return Full;
        }
        else
        {
            if (IsFull()) return this;
            // Check whether this interval will be empty after expansion, allowing
            // for a 1-bit rounding error when computing each endpoint.
            if (GetLength() + 2 * margin - 2 * S2.DoubleEpsilon <= 0) return Empty;
        }
        var lo = Math.IEEERemainder(Lo - margin, S2.M_2_PI);
        var hi = Math.IEEERemainder(Hi + margin, S2.M_2_PI);
        if (lo <= -Math.PI) lo = Math.PI;

        return new S1Interval(lo, hi);
    }

    // Return the smallest interval that contains this interval and the
    // given interval "y".
    public S1Interval Union(S1Interval y)
    {
        // The y.IsFull case is handled correctly in all cases by the code
        // below, but can follow three separate code paths depending on whether
        // this interval is inverted, is non-inverted but contains Pi, or neither.

        if (y.IsEmpty()) return this;
        if (FastContains(y.Lo))
        {
            if (FastContains(y.Hi))
            {
                // Either this interval contains y, or the union of the two
                // intervals is the Full() interval.
                if (Contains(y)) return this;  // IsFull code path
                return Full;
            }
            return new S1Interval(Lo, y.Hi);
        }
        if (FastContains(y.Hi)) return new S1Interval(y.Lo, Hi);

        // This interval contains neither endpoint of y.  This means that either y
        // contains all of this interval, or the two intervals are disjoint.
        if (IsEmpty() || y.FastContains(Lo)) return y;

        // Check which pair of endpoints are closer together.
        var dlo = PositiveDistance(y.Hi, Lo);
        var dhi = PositiveDistance(Hi, y.Lo);
        if (dlo < dhi)
        {
            return new S1Interval(y.Lo, Hi);
        }
        else
        {
            return new S1Interval(Lo, y.Hi);
        }
    }

    // Return the smallest interval that contains the intersection of this
    // interval with "y".  Note that the region of intersection may
    // consist of two disjoint intervals.
    public S1Interval Intersection(S1Interval y)
    {
        // The y.IsFull case is handled correctly in all cases by the code
        // below, but can follow three separate code paths depending on whether
        // this interval is inverted, is non-inverted but contains Pi, or neither.

        if (y.IsEmpty()) return Empty;
        if (FastContains(y.Lo))
        {
            if (FastContains(y.Hi))
            {
                // Either this interval contains y, or the region of intersection
                // consists of two disjoint subintervals.  In either case, we want
                // to return the shorter of the two original intervals.
                if (y.GetLength() < GetLength()) return y;  // IsFull code path
                return this;
            }
            return new S1Interval(y.Lo, Hi);
        }
        if (FastContains(y.Hi)) return new S1Interval(Lo, y.Hi);

        // This interval contains neither endpoint of y.  This means that either y
        // contains all of this interval, or the two intervals are disjoint.

        if (y.FastContains(Lo)) return this;  // IsEmpty okay here
        MyDebug.Assert(!Intersects(y));
        return Empty;
    }

    // Return true if this interval can be transformed into the given interval by
    // moving each endpoint by at most "maxError" (and without the endpoints
    // crossing, which would invert the interval).  Empty and full intervals are
    // considered to start at an arbitrary point on the unit circle, thus any
    // interval with (length <= 2*maxError) matches the empty interval, and any
    // interval with (length >= 2*Pi - 2*maxError) matches the full interval.
    public bool ApproxEquals(S1Interval y, double maxError = S2.DoubleError)
    {
        // Full and empty intervals require special cases because the "endpoints"
        // are considered to be positioned arbitrarily.
        if (IsEmpty()) return y.GetLength() <= 2 * maxError;
        if (y.IsEmpty()) return GetLength() <= 2 * maxError;
        if (IsFull()) return y.GetLength() >= 2 * (Math.PI - maxError);
        if (y.IsFull()) return GetLength() >= 2 * (Math.PI - maxError);

        // The purpose of the last test below is to verify that moving the endpoints
        // does not invert the interval, e.g. [-1e20, 1e20] vs. [1e20, -1e20].
        return 
            Math.Abs(Math.IEEERemainder(y.Lo - Lo, S2.M_2_PI)) <= maxError &&
            Math.Abs(Math.IEEERemainder(y.Hi - Hi, S2.M_2_PI)) <= maxError &&
            Math.Abs(GetLength() - y.GetLength()) <= 2 * maxError;
    }

    // Return true if the interval (which is closed) contains the point 'p'.
    // Skips the normalization of 'p' from -Pi to Pi.
    private bool FastContains(double p)
    {
        if (IsInverted())
        {
            return (p >= Lo || p <= Hi) && !IsEmpty();
        }
        else
        {
            return p >= Lo && p <= Hi;
        }
    }

    private static double PositiveDistance(double a, double b)
    {
        // Compute the distance from "a" to "b" in the range [0, 2*Pi).
        // This is equivalent to (remainder(b - a - Math.PI, S2Constants.M_2_PI) + Math.PI),
        // except that it is more numerically stable (it does not lose
        // precision for very small positive distances).
        var d = b - a;
        if (d >= 0) return d;
        // We want to ensure that if b == Pi and a == (-Pi + eps),
        // the return result is approximately 2*Pi and not zero.
        return b + Math.PI - (a - Math.PI);
    }

    #endregion

    #region Object

    public override string ToString() => $"[{Lo}, {Hi}]";
    public string ToDebugString() => $"[{S1Angle.FromRadians(Lo)}, {S1Angle.FromRadians(Hi)}]";

    #endregion
}
