//https://github.com/McNeight/MathM/blob/master/src/MathM.cs

// <copyright file="MathM.cs" company="Neil McNeight">
// Copyright © 2015-2019 Nathan P Jones.
// Copyright © 2019 Neil McNeight. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// full license information.
// </copyright>

namespace S2Geometry;

using System.Globalization;
using System.Runtime.CompilerServices;

/// <summary>
/// Provides constants and static methods for trigonometric, logarithmic,
/// and other common mathematical functions.
/// </summary>
/// <remarks>
/// The static fields and methods of the <see cref="MathM"/> class
/// correspond to those of the <c>System.MathF</c> classes, except that their
/// parameters are of type <see cref="decimal"/> rather than
/// <see cref="float"/>, and they return <see cref="decimal"/> rather
/// than <see cref="float"/> values.
/// </remarks>
public static class MathM
{
    /// <summary>
    /// Represents the natural logarithmic base, specified by the constant, e.
    /// </summary>
    /// <remarks>
    /// According to <see href="https://oeis.org/A001113">The On-Line Encyclopedia of Integer Sequences</see>,
    /// the value of e (also called Euler's number or Napier's constant) to 100 decimal places is
    /// 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274.
    /// </remarks>
    public const decimal E = 2.7182818284590452353602874714m;

    /// <summary>
    /// Represents the ratio of the circumference of a circle to its
    /// diameter, specified by the constant, π.
    /// </summary>
    /// <remarks>
    /// According to <see href="https://oeis.org/A000796">The On-Line Encyclopedia of Integer Sequences</see>,
    /// the value of π (also called Archimedes's constant) to 100 decimal places is
    /// 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679.
    /// </remarks>
    public const decimal PI = 3.1415926535897932384626433833m;

    /// <summary>
    /// Represents the ratio of the circumference of a circle to its
    /// radius, specified by the constant, τ.
    /// </summary>
    /// <remarks>
    /// According to <see href="https://tauday.com/tau-manifesto">The Tau Manifesto</see>,
    /// the value of τ to 100 decimal places is
    /// 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359.
    /// </remarks>
    public const decimal TAU = 6.2831853071795864769252867666m;

    /// <summary>
    /// π/2 - in radians is equivalent to 90 degrees.
    /// </summary>
    private const decimal PiHalf = 1.5707963267948966192313216916m;

    /// <summary>
    /// π/4 - in radians is equivalent to 45 degrees.
    /// </summary>
    private const decimal PiQuarter = 0.7853981633974483096156608458m;

    /// <summary>
    /// π/12 - in radians is equivalent to 15 degrees.
    /// </summary>
    private const decimal PiTwelfth = 0.2617993877991494365385536153m;

    /// <summary>
    /// The value of the natural logarithm of 10.
    /// </summary>
    /// <remarks>
    /// According to <see href="https://oeis.org/A002392">The On-Line Encyclopedia of Integer Sequences</see>,
    /// the value of Ln10 to 86 decimal places is
    /// 2.30258509299404568401799145468436420760110148862877297603332790096757260967735248023599.
    /// </remarks>
    private const decimal Ln10 = 2.3025850929940456840179914547m;

    /// <summary>
    /// The value of the natural logarithm of 2.
    /// </summary>
    /// <remarks>
    /// According to <see href="https://oeis.org/A002162">The On-Line Encyclopedia of Integer Sequences</see>,
    /// the value of Ln2 to 99 decimal places is
    /// 0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687.
    /// </remarks>
    private const decimal Ln2 = 0.6931471805599453094172321215m;

    private const int MaxRoundingDigits = 28;

    private const decimal DecimalRoundLimit = 1E28m;

    // Sign mask for the flags field. A value of zero in this bit indicates a
    // positive Decimal value, and a value of one in this bit indicates a
    // negative Decimal value.
    //
    // Look at OleAut's DECIMAL_NEG constant to check for negative values
    // in native code.
    private const int SignMask = unchecked((int)0x80000000);

    // Scale mask for the flags field. This byte in the flags field contains
    // the power of 10 to divide the Decimal value by. The scale byte must
    // contain a value between 0 and 28 inclusive.
    private const int ScaleMask = 0x00FF0000;

    // Number of bits scale is shifted by.
    private const int ScaleShift = 16;

    /// <summary>
    /// Smallest non-zero decimal value.
    /// </summary>
    /// <remarks>
    /// <code>new decimal(1, 0, 0, false, 28);</code> or 1e-28m.
    /// </remarks>
    private const decimal SmallestNonZeroDec = 0.0000000000000000000000000001m;

    // This table is required for the Round function which can specify the number of digits to round to
    private static readonly decimal[] roundPower10Decimal = new decimal[]
    {
        1E0m,  1E1m,  1E2m,  1E3m,  1E4m,  1E5m,  1E6m,  1E7m,  1E8m,  1E9m,
        1E10m, 1E11m, 1E12m, 1E13m, 1E14m, 1E15m, 1E16m, 1E17m, 1E18m, 1E19m,
        1E20m, 1E21m, 1E22m, 1E23m, 1E24m, 1E25m, 1E26m, 1E27m, 1E28m,
    };

    /// <summary>
    /// Returns the absolute value of a <see cref="decimal"/> floating-point number.
    /// </summary>
    /// <param name="m">A number that is greater than or equal to <see cref="decimal.MinValue"/>, but less than or equal to <see cref="decimal.MaxValue"/>.</param>
    /// <returns>A number, x, such that 0 ≤ x ≤ <see cref="decimal.MaxValue"/>.</returns>
    /// <remarks>
    /// The absolute value of a <see cref="decimal"/> is its numeric value without its sign. For example, the absolute value of both 1.2 and -1.2 is 1.2.
    /// </remarks>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Abs(decimal m)
    {
        // return new decimal(in d, d.flags & ~SignMask);
        //
        // Since decimal.Abs() is inaccessable, perform the operation
        // the old fashioned way.
        if (m < 0m)
        {
            m = -m;
        }

        return m;
    }

    /// <summary>
    /// Returns the angle whose cosine is the specified number.
    /// </summary>
    /// <param name="m">A number representing a cosine, where -1 ≤ <paramref name="m"/> ≤ 1.</param>
    /// <returns>An angle, θ, measured in radians, such that 0 ≤ θ ≤ π.</returns>
    /// <remarks>
    /// Multiply the return value by (180 / <see cref="PI"/>) to convert from radians to degrees.
    /// </remarks>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Acos(decimal m)
    {
        // See http://en.wikipedia.org/wiki/Inverse_trigonometric_function
        // and http://mathworld.wolfram.com/InverseCosine.html
        if (m < -1m || m > 1m)
        {
            throw new ArgumentOutOfRangeException(nameof(m), string.Format(CultureInfo.CurrentCulture, SR.ArgumentOutOfRange_Range, -1m, 1m));
        }

        // Special cases
        if (m == -1m)
        {
            return PI;
        }

        if (m == 0m)
        {
            return PiHalf;
        }

        if (m == 1m)
        {
            return 0m;
        }

        return 2m * Atan(Sqrt(1m - (m * m)) / (1m + m));
    }

    /// <summary>
    /// Returns the angle whose hyperbolic cosine is the specified number.
    /// </summary>
    /// <param name="m">A number representing a hyperbolic cosine, where <paramref name="m"/> must be greater than or equal to 1, but less than or equal to ∞.</param>
    /// <returns>An angle, θ, measured in radians, such that 0 ≤ θ ≤ ∞.</returns>
    /// <remarks>
    /// Multiply the return value by (180 / <see cref="PI"/>) to convert from radians to degrees.
    /// </remarks>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Acosh(decimal m) { throw null; }

    /// <summary>
    /// Returns the angle whose sine is the specified number.
    /// </summary>
    /// <param name="m">A number representing a sine, where -1 ≤ m ≤ 1.</param>
    /// <remarks>
    /// See http://en.wikipedia.org/wiki/Inverse_trigonometric_function
    /// and http://mathworld.wolfram.com/InverseSine.html
    /// I originally used the Taylor series for ASin, but it was extremely slow
    /// around -1 and 1 (millions of iterations) and still ends up being less
    /// accurate than deriving from the ATan function.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Asin(decimal m)
    {
        if (m < -1m || m > 1m)
        {
            throw new ArgumentOutOfRangeException(nameof(m), string.Format(CultureInfo.CurrentCulture, SR.ArgumentOutOfRange_Range, -1m, 1m));
        }

        // Special cases
        if (m == -1m)
        {
            return -PiHalf;
        }

        if (m == 0m)
        {
            return 0m;
        }

        if (m == 1m)
        {
            return PiHalf;
        }

        return 2m * Atan(m / (1m + Sqrt(1m - (m * m))));
    }

    /// <summary>
    /// Returns the angle whose hyperbolic sine is the specified number.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Asinh(decimal m) { throw null; }

    /// <summary>
    /// Returns the angle whose tangent is the specified number.
    /// </summary>
    /// <param name="m">A number representing a tangent.</param>
    /// <remarks>
    /// See http://mathworld.wolfram.com/InverseTangent.html for faster converging
    /// series from Euler that was used here.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Atan(decimal m)
    {
        // Special cases
        if (m == -1)
        {
            return -PiQuarter;
        }
        else if (m == 0)
        {
            return 0;
        }
        else if (m == 1)
        {
            return PiQuarter;
        }

        if (m < -1)
        {
            // Force down to -1 to 1 interval for faster convergence
            return -PiHalf - Atan(1 / m);
        }
        else if (m > 1)
        {
            // Force down to -1 to 1 interval for faster convergence
            return PiHalf - Atan(1 / m);
        }

        var result = 0m;
        var doubleIteration = 0; // current iteration * 2
        var y = (m * m) / (1 + (m * m));
        var nextAdd = 0m;

        while (true)
        {
            if (doubleIteration == 0)
            {
                nextAdd = m / (1 + (m * m));  // is = y / x  but this is better for very small numbers where y = 9
            }
            else
            {
                // We multiply by -1 each time so that the sign of the component
                // changes each time. The first item is positive and it
                // alternates back and forth after that.
                // Following is equivalent to: nextAdd *= y * (iteration * 2) / (iteration * 2 + 1);
                nextAdd *= y * doubleIteration / (doubleIteration + 1);
            }

            if (nextAdd == 0)
            {
                break;
            }

            result += nextAdd;

            doubleIteration += 2;
        }

        return result;
    }

    /// <summary>
    /// Returns the angle whose tangent is the quotient of two specified numbers.
    /// </summary>
    /// <param name="y">The y coordinate of a point.</param>
    /// <param name="x">The x coordinate of a point.</param>
    /// <returns>
    /// An angle, θ, measured in radians, such that -π ≤ θ ≤ π, and tan(θ) = y / x,
    /// where (x, y) is a point in the Cartesian plane. Observe the following:
    /// For (x, y) in quadrant 1, 0 &lt; θ &lt; π/2.
    /// For (x, y) in quadrant 2, π/2 &lt; θ ≤ π.
    /// For (x, y) in quadrant 3, -π &lt; θ &lt; -π/2.
    /// For (x, y) in quadrant 4, -π/2 &lt; θ &lt; 0.
    /// </returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Atan2(decimal y, decimal x)
    {
        if (x == 0 && y == 0)
        {
            return 0;                   // X0, Y0
        }

        if (x == 0)
        {
            return y > 0
                       ? PiHalf         // X0, Y+
                       : -PiHalf;       // X0, Y-
        }

        if (y == 0)
        {
            return x > 0
                       ? 0              // X+, Y0
                       : PI;            // X-, Y0
        }

        var aTan = Atan(y / x);

        if (x > 0)
        {
            return aTan;         // Q1&4: X+, Y+-
        }

        return y > 0
                   ? aTan + PI          // Q2: X-, Y+
                   : aTan - PI;         //   Q3: X-, Y-
    }

    /// <summary>
    /// Returns the angle whose hyperbolic tangent is the specified number.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Atanh(decimal m) { throw null; }

    /// <summary>
    /// Returns the cube root of a specified number.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Cbrt(decimal m) { throw null; }

    /// <summary>
    /// Returns the smallest integral value that is greater than or equal
    /// to the specified <see cref="decimal"/> floating-point number.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Ceiling(decimal m)
    {
        return decimal.Ceiling(m);
    }

    /// <summary>
    /// Returns the cosine of the specified angle.
    /// </summary>
    /// <param name="m">An angle, measured in radians.</param>
    /// <remarks>
    /// Uses a Taylor series to calculate sine. See
    /// http://en.wikipedia.org/wiki/Trigonometric_functions for details.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Cos(decimal m)
    {
        // Normalize to between -2Pi <= m <= 2Pi
        m = Remainder(m, TAU);

        if (m == 0 || m == TAU)
        {
            return 1m;
        }

        if (m == PI)
        {
            return -1m;
        }

        if (m == PiHalf || m == PI + PiHalf)
        {
            return 0m;
        }

        var result = 0m;
        var doubleIteration = 0; // current iteration * 2
        var xSquared = m * m;
        var nextAdd = 0m;

        while (true)
        {
            if (doubleIteration == 0)
            {
                nextAdd = 1m;
            }
            else
            {
                // We multiply by -1 each time so that the sign of the component
                // changes each time. The first item is positive and it
                // alternates back and forth after that.
                // Following is equivalent to: nextAdd *= -1 * x * x / ((2 * iteration - 1) * (2 * iteration));
                nextAdd *= -1 * xSquared / ((doubleIteration * doubleIteration) - doubleIteration);
            }

            if (nextAdd == 0m)
            {
                break;
            }

            result += nextAdd;

            doubleIteration += 2;
        }

        return result;
    }

    /// <summary>
    /// Returns the hyperbolic cosine of the specified angle.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
    public static decimal Cosh(decimal m) { throw null; }

    /// <summary>
    /// Returns <see cref="E">e</see> raised to the specified power.
    /// </summary>
    /// <param name="m">A number specifying a power.</param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Exp(decimal m)
    {
        decimal result;
        decimal nextAdd;
        int iteration;
        bool reciprocal;
        decimal t;

        reciprocal = m < 0;
        m = Abs(m);

        t = Truncate(m);

        if (m == 0)
        {
            result = 1;
        }
        else if (m == 1)
        {
            result = E;
        }
        else if (Abs(m) > 1 && t != m)
        {
            // Split up into integer and fractional
            result = Exp(t) * Exp(m - t);
        }
        else if (m == t)
        {
            // Integer power
            result = ExpBySquaring(E, m);
        }
        else
        {
            // Fractional power < 1
            // See http://mathworld.wolfram.com/ExponentialFunction.html
            iteration = 0;
            nextAdd = 0;
            result = 0;

            while (true)
            {
                if (iteration == 0)
                {
                    nextAdd = 1;               // == Pow(d, 0) / Factorial(0) == 1 / 1 == 1
                }
                else
                {
                    nextAdd *= m / iteration;  // == Pow(d, iteration) / Factorial(iteration)
                }

                if (nextAdd == 0)
                {
                    break;
                }

                result += nextAdd;

                iteration += 1;
            }
        }

        // Take reciprocal if this was a negative power
        // Note that result will never be zero at this point.
        if (reciprocal)
        {
            result = 1 / result;
        }

        return result;
    }

    /// <summary>
    /// Returns the largest integral value less than or equal to the
    /// specified <see cref="decimal"/> floating-point number.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Floor(decimal m)
    {
        return decimal.Floor(m);
    }

    /// <summary>
    /// Returns the remainder resulting from the division of a specified number by another specified number.
    /// </summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal IEEERemainder(decimal x, decimal y) { throw null; }

    /// <summary>
    /// Returns the natural (base e) logarithm of a specified number.
    /// </summary>
    /// <param name="m">A number whose logarithm is to be found.</param>
    /// <remarks>
    /// I'm still not satisfied with the speed. I tried several different
    /// algorithms that you can find in a historical version of this
    /// source file. The one I settled on was the best of mediocrity.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Log(decimal m)
    {
        if (m < 0)
        {
            throw new ArgumentException("Natural logarithm is a complex number for values less than zero!", nameof(m));
        }

        if (m == 0)
        {
            throw new OverflowException("Natural logarithm is defined as negative infinity at zero which the Decimal data type can't represent!");
        }

        if (m == 1)
        {
            return 0;
        }

        if (m >= 1)
        {
            decimal power = 0m;
            System.Math.Log(0d);

            decimal x = m;
            while (x > 1)
            {
                x /= 10;
                power += 1;
            }

            return Log(x) + (power * Ln10);
        }

        // See http://en.wikipedia.org/wiki/Natural_logarithm#Numerical_value
        // for more information on this faster-converging series.

        decimal y;
        decimal ySquared;

        decimal iteration = 0;
        decimal exponent = 0m;
        decimal nextAdd = 0m;
        decimal result = 0m;

        y = (m - 1) / (m + 1);
        ySquared = y * y;

        while (true)
        {
            if (iteration == 0)
            {
                exponent = 2 * y;
            }
            else
            {
                exponent *= ySquared;
            }

            nextAdd = exponent / ((2 * iteration) + 1);

            if (nextAdd == 0)
            {
                break;
            }

            result += nextAdd;

            iteration += 1;
        }

        return result;
    }

    /// <summary>
    /// Returns the logarithm of a specified number in a specified base.
    /// </summary>
    /// <param name="m">A number whose logarithm is to be found.</param>
    /// <param name="newBase">The base of the logarithm.</param>
    /// <remarks>
    /// This is a relatively naive implementation that simply divides the
    /// natural log of <paramref name="m"/> by the natural log of the base.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Log(decimal m, decimal newBase)
    {
        // Substituting null for double.NaN
        // if (!m.HasValue)
        // {
        //    return m;
        // }

        // Substituting null for double.NaN
        // if (!newBase.HasValue)
        // {
        //    return newBase;
        // }

        // Substituting null for double.NaN
        if (newBase == 1)
        {
            double y = 13.0d / 0.0d;
            // throw new InvalidOperationException("Logarithm for base 1 is undefined.");
            // return null;
        }

        // if ((a != 1) && ((newBase == 0) || double.IsPositiveInfinity(newBase)))
        if ((m != 1) && (newBase == 0))
        {
            // return null;
        }

        if (m < 0)
        {
            throw new ArgumentException("Logarithm is a complex number for values less than zero!", nameof(m));
        }

        if (m == 0)
        {
            throw new OverflowException("Logarithm is defined as negative infinity at zero which the Decimal data type can't represent!");
        }

        if (newBase < 0)
        {
            throw new ArgumentException("Logarithm base would be a complex number for values less than zero!", nameof(newBase));
        }

        if (newBase == 0)
        {
            throw new OverflowException("Logarithm base would be negative infinity at zero which the Decimal data type can't represent!");
        }

        // Short circuit the checks below if m is 1 because
        // that will yield 0 in the numerator below and give us
        // 0 for any base, even ones that would yield infinity.
        if (m == 1)
        {
            return 0m;
        }

        return Log(m) / Log(newBase);
    }

    /// <summary>
    /// Returns the base 10 logarithm of a specified number.
    /// </summary>
    /// <param name="m">A number whose logarithm is to be found.</param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Log10(decimal m)
    {
        // Substituting null for double.NaN
        // if (!m.HasValue)
        // {
        //    return null;
        // }

        if (m < 0)
        {
            throw new ArgumentException("Logarithm is a complex number for values less than zero!", nameof(m));
        }

        if (m == 0)
        {
            throw new OverflowException("Logarithm is defined as negative infinity at zero which the Decimal data type can't represent!");
        }

        return Log(m) / Ln10;
    }

    /// <summary>
    /// Returns the base 2 logarithm of a specified number.
    /// </summary>
    /// <param name="m">A number whose logarithm is to be found.</param>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    private static decimal Log2(decimal m)
    {
        if (m < 0)
        {
            throw new ArgumentException("Logarithm is a complex number for values less than zero!", nameof(m));
        }

        if (m == 0)
        {
            throw new OverflowException("Logarithm is defined as negative infinity at zero which the Decimal data type can't represent!");
        }

        return Log(m) / Ln2;
    }

    /// <summary>
    /// Returns the larger of two single-precision floating-point numbers.
    /// </summary>
    /// <param name="m1"></param>
    /// <param name="m2"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Max(decimal m1, decimal m2)
    {
        // return decimal.Max(ref val1, ref val2);
        //
        // Since decimal.Max() is inaccessable, perform the comparison
        // the old fashioned way.
        return (m1 >= m2) ? m1 : m2;
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="m1"></param>
    /// <param name="m2"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Min(decimal m1, decimal m2)
    {
        // return decimal.Min(ref val1, ref val2);
        //
        // Since decimal.Min() is inaccessable, perform the comparison
        // the old fashioned way.
        return (m1 <= m2) ? m1 : m2;
    }

    /// <summary>
    /// Returns a specified number <paramref name="x"/> raised to the specified power <paramref name="y"/>.
    /// </summary>
    /// <param name="x">A number to be raised to a power.</param>
    /// <param name="y">A number that specifies a power.</param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Pow(decimal x, decimal y)
    {
        decimal result;
        var isNegativeExponent = false;

        // Handle negative exponents
        if (y < 0)
        {
            isNegativeExponent = true;
            y = Abs(y);
        }

        if (y == 0)
        {
            result = 1;
        }
        else if (y == 1)
        {
            result = x;
        }
        else
        {
            var t = decimal.Truncate(y);

            if (y == t)
            {
                // Integer powers
                result = ExpBySquaring(x, y);
            }
            else
            {
                // Fractional power < 1
                // See http://en.wikipedia.org/wiki/Exponent#Real_powers
                // The next line is an optimization of Exp(y * Log(x)) for better precision
                result = ExpBySquaring(x, t) * Exp((y - t) * Log(x));
            }
        }

        if (isNegativeExponent)
        {
            // Note, for IEEE floats this would be Infinity and not an exception...
            if (result == 0)
            {
                throw new OverflowException("Negative power of 0 is undefined!");
            }

            result = 1 / result;
        }

        return result;
    }

    /// <summary>
    /// Rounds a <see cref="decimal"/> value to a given number of decimal places.
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Round(decimal m)
    {
#if NETFRAMEWORK || NETSTANDARD2_0 || NETSTANDARD2_1 || NETCOREAPP2_0 || NETCOREAPP2_1 || NETCOREAPP2_2 || NETCOREAPP3_0
        return decimal.Round(m);
#else
        return Round(m, 0, MidpointRounding.ToEven);
#endif
    }

    /// <summary>
    /// Rounds a <see cref="decimal"/> value to a given number of decimal places.
    /// </summary>
    /// <param name="m"></param>
    /// <param name="digits"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Round(decimal m, int digits)
    {
#if NETFRAMEWORK || NETSTANDARD2_0 || NETSTANDARD2_1 || NETCOREAPP2_0 || NETCOREAPP2_1 || NETCOREAPP2_2 || NETCOREAPP3_0
        return decimal.Round(m, digits);
#else
        return Round(m, digits, MidpointRounding.ToEven);
#endif
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="m"></param>
    /// <param name="mode"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Round(decimal m, MidpointRounding mode)
    {
#if NETFRAMEWORK || NETSTANDARD2_0 || NETSTANDARD2_1 || NETCOREAPP2_0 || NETCOREAPP2_1 || NETCOREAPP2_2 || NETCOREAPP3_0
        return decimal.Round(m, mode);
#else
        return Round(m, 0, mode);
#endif
    }

#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Round(decimal m, int digits, MidpointRounding mode)
    {
#if NETFRAMEWORK || NETSTANDARD2_0 || NETSTANDARD2_1 || NETCOREAPP2_0 || NETCOREAPP2_1 || NETCOREAPP2_2 || NETCOREAPP3_0
        return decimal.Round(m, digits, mode);
#else
        if ((digits < 0) || (digits > MaxRoundingDigits))
        {
            throw new ArgumentOutOfRangeException(nameof(digits), string.Format(CultureInfo.CurrentCulture, SR.ArgumentOutOfRange_RoundingDigits, MaxRoundingDigits));
        }

        // if (mode < MidpointRounding.ToEven || mode > MidpointRounding.ToPositiveInfinity)
        if (mode < MidpointRounding.ToEven)
        {
            throw new ArgumentException(string.Format(CultureInfo.CurrentCulture, SR.Argument_InvalidEnumValue, mode, nameof(MidpointRounding)), nameof(mode));
        }

        var bits = decimal.GetBits(m);
        var lo = bits[0];
        var mid = bits[1];
        var hi = bits[2];
        var flags = bits[3];

        int Scale = (byte)(flags >> ScaleShift);
        // int scale = m.Scale - digits;

        if (Abs(m) < DecimalRoundLimit)
        {
            decimal power10 = roundPower10Decimal[digits];

            m *= power10;

            switch (mode)
            {
                // Rounds to the nearest value; if the number falls midway,
                // it is rounded to the nearest value with an even least significant digit
                case MidpointRounding.ToEven:
                    {
                        m = Round(m);
                        break;
                    }

                // Rounds to the nearest value; if the number falls midway,
                // it is rounded to the nearest value above (for positive numbers) or below (for negative numbers)
                case MidpointRounding.AwayFromZero:
                    {
                        decimal fraction = 5; // ModF(x, &x);

                        if (Abs(fraction) >= 0.5m)
                        {
                            m += Sign(fraction);
                        }

                        break;
                    }

                default:
                    {
                        throw new ArgumentException(string.Format(CultureInfo.CurrentCulture, SR.Argument_InvalidEnumValue, mode, nameof(MidpointRounding)), nameof(mode));
                    }
            }

            m /= power10;
        }

        return m;
#endif
    }

    /// <summary>
    /// Returns a value indicating the sign of a <see cref="decimal"/> number.
    /// </summary>
    /// <param name="m">A signed <see cref="decimal"/> number.</param>
    /// <returns>
    /// A number that indicates the sign of <c>m</c>, as shown in the following table.
    /// <list type="table">
    ///     <listheader>
    ///         <term>Return value</term>
    ///         <description>Meaning</description>
    ///     </listheader>
    ///     <item>
    ///         <term>-1</term>
    ///         <description><c>m</c> is less than zero.</description>
    ///     </item>
    ///     <item>
    ///         <term>0</term>
    ///         <description><c>m</c> is equal to zero.</description>
    ///     </item>
    ///     <item>
    ///         <term>1</term>
    ///         <description><c>m</c> is greater than zero.</description>
    ///     </item>
    /// </list>
    /// </returns>
    /// <exception cref="ArithmeticException">
    /// <c>m</c> is equal to <see cref="null"/>.
    /// </exception>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static int Sign(decimal m)
    {
        // return decimal.Sign(ref value);
        //
        // Since decimal.Sign() is inaccessable, perform the comparison
        // the old fashioned way.
        if (m == null)
        {
            throw new ArithmeticException(SR.Arithmetic_Null);
        }
        else if (m < 0)
        {
            return -1;
        }
        else if (m > 0)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    /// <summary>
    /// Returns the sine of the specified angle.
    /// </summary>
    /// <param name="m">An angle, measured in radians.</param>
    /// <remarks>
    /// Uses a Taylor series to calculate sine. See
    /// http://en.wikipedia.org/wiki/Trigonometric_functions for details.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Sin(decimal m)
    {
        // Normalize to between -2Pi <= m <= 2Pi
        m = Remainder(m, TAU);

        if (m == 0 || m == PI || m == TAU)
        {
            return 0;
        }

        if (m == PiHalf)
        {
            return 1;
        }

        if (m == PI + PiHalf)
        {
            return -1;
        }

        var result = 0m;
        var doubleIteration = 0; // current iteration * 2
        var mSquared = m * m;
        var nextAdd = 0m;

        while (true)
        {
            if (doubleIteration == 0)
            {
                nextAdd = m;
            }
            else
            {
                // We multiply by -1 each time so that the sign of the component
                // changes each time. The first item is positive and it
                // alternates back and forth after that.
                // Following is equivalent to: nextAdd *= -1 * m * m / ((2 * iteration) * (2 * iteration + 1));
                nextAdd *= -1 * mSquared / ((doubleIteration * doubleIteration) + doubleIteration);
            }

            // Debug.WriteLine("{0:000}:{1,33:+0.0000000000000000000000000000;-0.0000000000000000000000000000} ->{2,33:+0.0000000000000000000000000000;-0.0000000000000000000000000000}",
            //    doubleIteration / 2, nextAdd, result + nextAdd);
            if (nextAdd == 0)
            {
                break;
            }

            result += nextAdd;

            doubleIteration += 2;
        }

        return result;
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Sinh(decimal m) { throw null; }

    /// <summary>
    /// Returns the square root of a given number.
    /// </summary>
    /// <param name="m">A non-negative number.</param>
    /// <remarks>
    /// Uses an implementation of the "Babylonian Method".
    /// See http://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Babylonian_method
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Sqrt(decimal m)
    {
        if (m < 0m)
        {
            throw new ArgumentException("Square root not defined for Decimal data type when less than zero!", nameof(m));
        }

        // Prevent divide-by-zero errors below. Dividing either
        // of the numbers below will yield a recurring 0 value
        // for halfS eventually converging on zero.
        if (m == 0m || m == SmallestNonZeroDec)
        {
            return 0m;
        }

        decimal x;
        var halfS = m / 2m;
        var lastX = -1m;
        decimal nextX;

        // Begin with an estimate for the square root.
        // Use hardware to get us there quickly.
        x = (decimal)Math.Sqrt(decimal.ToDouble(m));

        while (true)
        {
            nextX = (x / 2m) + (halfS / x);

            // The next check effectively sees if we've ran out of
            // precision for our data type.
            if (nextX == x || nextX == lastX)
            {
                break;
            }

            lastX = x;
            x = nextX;
        }

        return nextX;
    }

    /// <summary>
    /// Returns the tangent of the specified angle.
    /// </summary>
    /// <param name="m">An angle, measured in radians.</param>
    /// <remarks>
    /// Uses a Taylor series to calculate sine. See
    /// http://en.wikipedia.org/wiki/Trigonometric_functions for details.
    /// </remarks>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Tan(decimal m)
    {
        try
        {
            return Sin(m) / Cos(m);
        }
        catch (DivideByZeroException)
        {
            throw new Exception("Tangent is undefined at this angle!");
        }
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="m"></param>
    /// <returns></returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Tanh(decimal m) { throw null; }

    /// <summary>
    /// Calculates the integral part of a specified <see cref="decimal"/> number.
    /// </summary>
    /// <param name="m">A number to truncate.</param>
    /// <returns>
    /// The integral part of m; that is, the number that remains after
    /// any fractional digits have been discarded.
    /// </returns>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    public static decimal Truncate(decimal m)
    {
        return decimal.Truncate(m);
    }

    /// <summary>
    /// Raises one number to an integral power.
    /// </summary>
    /// <remarks>
    /// See http://en.wikipedia.org/wiki/Exponentiation_by_squaring
    /// </remarks>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    private static decimal ExpBySquaring(decimal x, decimal y)
    {
        Debug.Assert(y >= 0 && decimal.Truncate(y) == y, "Only non-negative, integer powers supported.");
        if (y < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(y), "Negative exponents not supported!");
        }

        if (decimal.Truncate(y) != y)
        {
            throw new ArgumentException("Exponent must be an integer!", nameof(y));
        }

        var result = 1m;
        var multiplier = x;

        while (y > 0)
        {
            if ((y % 2) == 1)
            {
                result *= multiplier;
                y -= 1;
                if (y == 0)
                {
                    break;
                }
            }

            multiplier *= multiplier;
            y /= 2;
        }

        return result;
    }

    /// <summary>
    /// Gets the number of decimal places in a decimal value.
    /// </summary>
    /// <remarks>
    /// Started with something found here: http://stackoverflow.com/a/6092298/856595
    /// </remarks>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    private static int GetDecimalPlaces(decimal m, bool countTrailingZeros)
    {
        const int signMask = unchecked((int)0x80000000);
        const int scaleMask = 0x00FF0000;
        const int scaleShift = 16;

        var bits = decimal.GetBits(m);
        var result = (bits[3] & scaleMask) >> scaleShift;  // extract exponent

        // Return immediately for values without a fractional portion or if we're counting trailing zeros
        if (countTrailingZeros || (result == 0))
        {
            return result;
        }

        // Get a raw version of the decimal's integer
        bits[3] = bits[3] & ~unchecked(signMask | scaleMask); // clear out exponent and negative bit
        var rawValue = new decimal(bits);

        // Account for trailing zeros
        while ((result > 0) && ((rawValue % 10) == 0))
        {
            result--;
            rawValue /= 10;
        }

        return result;
    }

    /// <summary>
    /// Gets the remainder of one number divided by another number in such a way as to retain maximum precision.
    /// </summary>
#if !NET20 && !NET35 && !NET40
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
    private static decimal Remainder(decimal m1, decimal m2)
    {
        if (Abs(m1) < Abs(m2))
        {
            return m1;
        }

        var timesInto = Truncate(m1 / m2);
        var shiftingNumber = m2;
        var sign = Sign(m1);

        for (var i = 0; i <= GetDecimalPlaces(m2, true); i++)
        {
            // Note that first "digit" will be the integer portion of d2
            var digit = Truncate(shiftingNumber);

            m1 -= timesInto * (digit / roundPower10Decimal[i]);

            shiftingNumber = (shiftingNumber - digit) * 10m; // remove used digit and shift for next iteration
            if (shiftingNumber == 0m)
            {
                break;
            }
        }

        // If we've crossed zero because of the precision mismatch,
        // we need to add a whole d2 to get a correct result.
        if (m1 != 0 && Sign(m1) != sign)
        {
            m1 = Sign(m2) == sign
                     ? m1 + m2
                     : m1 - m2;
        }

        return m1;
    }

    private static class SR
    {
        public const string ArgumentOutOfRange_Range = "Valid values are between {0} and {1}, inclusive.";
        public const string Argument_InvalidEnumValue = "The value '{0}' is not valid for this usage of the type {1}.";
        public const string ArgumentOutOfRange_RoundingDigits = "Rounding digits must be between 0 and {0}, inclusive.";
        public const string Arithmetic_Null = "Function does not accept decimal null values.";
    }
}

