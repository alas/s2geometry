using System;
using System.Globalization;

namespace S2Geometry
{
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

        public static double Ldexp(double x, int exp)
        {
            return x * Math.Pow(2, exp);
        }

        /// <summary>
        /// from: https://stackoverflow.com/questions/374316/round-a-double-to-x-significant-figures (RoundToSignificantDigits)
        /// </summary>
        public static string Round(double d)
        {
            if (d == 0) return "0";

            var scale = Math.Pow(10, Math.Floor(Math.Log10(Math.Abs(d))) + 1);
            var t = scale * Math.Round(d / scale, 1);

            return string.Format(CultureInfo.InvariantCulture, "{0:g15}", t);
        }

        public static double Scale(double d, int decimalPlaces)
        {
            if (d == 0) return 0;
            double exp = Math.Floor(Math.Log10(Math.Abs(d)));
            return d / (Math.Pow(10, exp) * Math.Pow(10, decimalPlaces));
        }
    }
}
