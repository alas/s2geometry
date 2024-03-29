﻿// S2MaxDistance is a class that allows maximum distances to be computed using
// a minimum distance algorithm.  Specifically, S2MaxDistance(x) represents the
// supplementary distance (Pi - x).  This has the effect of inverting the sort
// order, i.e.
//
//  (S2MaxDistance(x) < S2MaxDistance(y))  <=>  (Pi - x < Pi - y)  <=>  (x > y)
//
// All other operations are implemented similarly (using the supplementary
// distance Pi - x).  For example, S2MaxDistance(x) - S2MaxDistance(y) ==
// S2MaxDistance(x + y).

namespace S2Geometry;

public readonly record struct S2MaxDistance(S1ChordAngle Distance) : IComparable<S2MaxDistance>, IDistance<S2MaxDistance>
{
    #region Fields, Constants

    public static S2MaxDistance Zero { get; } = new(S1ChordAngle.Straight);
    public static S2MaxDistance Infinity { get; } = new(S1ChordAngle.Negative);
    public static S2MaxDistance Negative { get; } = new(S1ChordAngle.Infinity);
    public static S2MaxDistance DefaultValue { get; } = new(S1ChordAngle.Zero);

    #endregion

    #region S2MaxDistance

    public bool IsNegative() => Distance == S1ChordAngle.Infinity;

    // TODO(Alas): Consider stricter check here -- only allow Negative().
    public bool IsInfinity() => Distance.Length2 < 0;

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
        return S1Angle.FromRadians(2 * Math.Asin(0.5 * Math.Sqrt(Distance.Length2)));
    }

    #endregion

    #region IDistance

    public S1ChordAngle GetChordAngleBound() => S1ChordAngle.Straight - Distance;

    // If (dist < that), updates that and returns true (used internally).
    public static bool UpdateMin(S2MaxDistance dist, ref S2MaxDistance that)
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

    // All of these functions have inverted sort order

    public int CompareTo(S2MaxDistance other)
    {
        if (other.Distance < Distance) return -1;
        if (other.Distance > Distance) return 1;
        return 0;
    }

    public static bool operator <(S2MaxDistance x, S2MaxDistance y) => x.Distance > y.Distance;
    public static bool operator >(S2MaxDistance x, S2MaxDistance y) => x.Distance < y.Distance;
    public static bool operator <=(S2MaxDistance x, S2MaxDistance y) => x.Distance >= y.Distance;
    public static bool operator >=(S2MaxDistance x, S2MaxDistance y) => x.Distance <= y.Distance;

    public static bool operator <(S2MaxDistance x, S1ChordAngle y) => x.Distance > y;
    public static bool operator >(S2MaxDistance x, S1ChordAngle y) => x.Distance < y;
    public static bool operator <=(S2MaxDistance x, S1ChordAngle y) => x.Distance >= y;
    public static bool operator >=(S2MaxDistance x, S1ChordAngle y) => x.Distance <= y;

    #endregion

    #region Operators

    public static S2MaxDistance operator -(S2MaxDistance x, S2MaxDistance delta) => new(x.Distance + delta.Distance);
    public static S2MaxDistance operator -(S2MaxDistance x, S1ChordAngle delta) => new(x.Distance + delta);

    #endregion

    #region Object

    public override string ToString() => ToAngle().ToString();

    #endregion
}
