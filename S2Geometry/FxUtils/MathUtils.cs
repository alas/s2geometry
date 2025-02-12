﻿namespace S2Geometry;

using System.Numerics;

public static class MathUtils
{
    #region NextAfter

    /// <summary>
    /// Gets the floating-point number that is next after <paramref name="fromNumber"/> in the direction of <paramref name="towardNumber"/>.
    /// </summary>
    /// <param name="fromNumber">A floating-point number.</param>
    /// <param name="towardNumber">A floating-point number.</param>
    /// <returns>The floating-point number that is next after <paramref name="fromNumber"/> in the direction of <paramref name="towardNumber"/>.</returns>
    /// <remarks>
    /// <para>
    /// IEC 60559 recommends that <paramref name="fromNumber"/> be returned whenever <c><paramref name="fromNumber"/> == <paramref name="towardNumber"/></c>.
    /// These functions return <paramref name="towardNumber"/> instead, which makes the behavior around zero consistent: <c><see cref="math.nextafter(double, double)">nextafter</see>(-0.0, +0.0)</c>
    /// returns <c>+0.0</c> and <c><see cref="math.nextafter(double, double)">nextafter</see>(+0.0, -0.0)</c> returns <c>–0.0</c>.
    /// </para>
    /// <para>
    /// See <a href="http://en.cppreference.com/w/c/numeric/math/nextafter">nextafter</a> in the C standard documentation.
    /// </para>
    /// from: https://github.com/MachineCognitis/C.math.NET/blob/master/C.math/math.cs
    /// </remarks>
    /// <example>
    /// <code language="C#">
    /// Assert.IsTrue(math.nextafter(0D, 0D) == 0D);
    /// Assert.IsTrue(math.nextafter(-0D, 0D) == 0D;
    /// Assert.IsTrue(math.nextafter(0D, -0D) == -0D);
    /// 
    /// Assert.IsTrue(math.nextafter(math.DBL_MIN, 0D) == math.DBL_DENORM_MAX);
    /// Assert.IsTrue(math.nextafter(math.DBL_DENORM_MIN, 0D) == 0D);
    /// Assert.IsTrue(math.nextafter(math.DBL_MIN, -0D) == math.DBL_DENORM_MAX);
    /// Assert.IsTrue(math.nextafter(math.DBL_DENORM_MIN, -0D) == 0D);
    /// 
    /// Assert.IsTrue(math.nextafter(0D, System.Double.PositiveInfinity) == math.DBL_DENORM_MIN);
    /// Assert.IsTrue(math.nextafter(-0D, System.Double.NegativeInfinity) == -math.DBL_DENORM_MIN);
    /// </code> 
    /// <code language="VB.NET">
    /// Assert.IsTrue(math.nextafter(0D, 0D) = 0D);
    /// Assert.IsTrue(math.nextafter(-0D, 0D) = 0D);
    /// Assert.IsTrue(math.nextafter(0D, -0D) = -0D);
    /// 
    /// Assert.IsTrue(math.nextafter(math.DBL_MIN, 0D) = math.DBL_DENORM_MAX);
    /// Assert.IsTrue(math.nextafter(math.DBL_DENORM_MIN, 0D) = 0D);
    /// Assert.IsTrue(math.nextafter(math.DBL_MIN, -0D) = math.DBL_DENORM_MAX);
    /// Assert.IsTrue(math.nextafter(math.DBL_DENORM_MIN, -0D) = 0D);
    /// 
    /// Assert.IsTrue(math.nextafter(0D, System.Double.PositiveInfinity) = math.DBL_DENORM_MIN);
    /// Assert.IsTrue(math.nextafter(-0D, System.Double.NegativeInfinity) = -math.DBL_DENORM_MIN);
    /// </code> 
    /// </example>
    public static double NextAfter(double fromNumber, double towardNumber)
    {
        // If either fromNumber or towardNumber is NaN, return NaN.
        if (double.IsNaN(towardNumber) || double.IsNaN(fromNumber))
        {
            return double.NaN;
        }
        // If no direction.
        if (fromNumber == towardNumber)
        {
            return towardNumber;
        }
        // If fromNumber is zero, return smallest subnormal.
        if (fromNumber == 0)
        {
            return (towardNumber > 0) ? double.Epsilon : -double.Epsilon;
        }
        // All other cases are handled by incrementing or decrementing the bits value.
        // Transitions to infinity, to subnormal, and to zero are all taken care of this way.
        long bits = BitConverter.DoubleToInt64Bits(fromNumber);
        // A xor here avoids nesting conditionals. We have to increment if fromValue lies between 0 and toValue.
        if ((fromNumber > 0) ^ (fromNumber > towardNumber))
        {
            bits += 1;
        }
        else
        {
            bits -= 1;
        }
        return BitConverter.Int64BitsToDouble(bits);
    }

    #endregion

    public static double Ldexp(double x, int exp) => x * Math.Pow(2, exp);
    public static decimal Ldexp(decimal x, int exp) => x * MathM.Pow(2, exp);
    public static ExactFloat Ldexp(ExactFloat x, int exp) => new() { Value = x.Value * MathM.Pow(2, exp) };

    // Computes v^i, where i is a non-negative integer.
    // When T is a floating point type, this has the same semantics as pow(), but
    // is much faster.
    // T can also be any integral type, in which case computations will be
    // performed in the value domain of this integral type, and overflow semantics
    // will be those of T.
    // You can also use any type for which operator*= is defined.
    public static double IPow(double @base, int exp)
    {
        MyDebug.Assert(exp >= 0);
        var uexp = (uint)exp;

        if (uexp < 16)
        {
            var result = ((uexp & 1) != 0) ? @base : 1.0;
            if (uexp >= 2)
            {
                @base *= @base;
                if ((uexp & 2) != 0)
                {
                    result *= @base;
                }
                if (uexp >= 4)
                {
                    @base *= @base;
                    if ((uexp & 4) != 0)
                    {
                        result *= @base;
                    }
                    if (uexp >= 8)
                    {
                        @base *= @base;
                        result *= @base;
                    }
                }
            }
            return result;
        }

        {
            var result = @base;
            var count = 31 ^ BitsUtils.Log2FloorNonZero(uexp);

            uexp <<= count;
            count ^= 31;

            while (count-- > 0)
            {
                uexp <<= 1;
                result *= result;
                if (uexp >= 0x80000000)
                {
                    result *= @base;
                }
            }

            return result;
        }
    }

    public static int Sgn<T>(T x) where T : IFloatingPointIeee754<T>
    {
        return x == T.Zero ? 0 : (x < T.Zero ? -1 : 1);
    }
}
