using System;
using Xunit;

namespace S2Geometry
{
    public class S1AngleTests
    {
        [Fact]
        public void Test_S1Angle_DefaultConstructor()
        {
            // Check that the default constructor returns an angle of 0.
            S1Angle a = S1Angle.FromRadians(0);
            Assert.Equal(0, a.Radians);
        }

        [Fact]
        public void Test_S1Angle_Infinity()
        {
            Assert.True(S1Angle.FromRadians(1e30) < S1Angle.Infinity);
            Assert.True(-S1Angle.Infinity < S1Angle.Zero);
            Assert.True(S1Angle.Infinity == S1Angle.FromRadians(double.PositiveInfinity));
        }

        [Fact]
        public void Test_S1Angle_Zero()
        {
            Assert.Equal(S1Angle.FromRadians(0), S1Angle.Zero);
        }

        [Fact]
        public void Test_S1Angle_PiRadiansExactly180Degrees()
        {
            // Check that the conversion between Pi radians and 180 degrees is exact.
            Assert.Equal(Math.PI, S1Angle.FromRadians(Math.PI).Radians);
            Assert.Equal(180.0, S1Angle.FromRadians(Math.PI).GetDegrees());
            Assert.Equal(Math.PI, S1Angle.FromDegrees(180).Radians);
            Assert.Equal(180.0, S1Angle.FromDegrees(180).GetDegrees());

            Assert.Equal(90.0, S1Angle.FromRadians(S2.M_PI_2).GetDegrees());

            // Check negative angles.
            Assert.Equal(-90.0, S1Angle.FromRadians(-S2.M_PI_2).GetDegrees());
            Assert.Equal(-S2.M_PI_4, S1Angle.FromDegrees(-45).Radians);
        }

        [Fact]
        public void Test_S1Angle_E5E6E7Representations()
        {
            // Check that E5/E6/E7 representations work as expected.
            Assert2.Near(S1Angle.FromDegrees(-45).Radians,
                             S1Angle.FromE5(-4500000).Radians);
            Assert2.Near(S1Angle.FromDegrees(-60).Radians,
                             S1Angle.FromE6(-60000000).Radians);
            Assert2.Near(S1Angle.FromDegrees(75).Radians,
                   S1Angle.FromE7(750000000).Radians);
            Assert.Equal(-17256123, S1Angle.FromDegrees(-172.56123).E5());
            Assert.Equal(12345678, S1Angle.FromDegrees(12.345678).E6());
            Assert.Equal(-123456789, S1Angle.FromDegrees(-12.3456789).E7());
        }

        [Fact]
        public void Test_S1Angle_E6E7RepresentationsUnsigned()
        {
            // Check that unsigned E6/E7 representations work as expected.
            Assert2.Near(
                S1Angle.FromDegrees(60).Radians,
                S1Angle.FromUnsignedE6((UInt32)60000000).Radians);
            Assert2.Near(
                S1Angle.FromDegrees(-60).Radians,
                S1Angle.FromUnsignedE6(0xFFFFFFFFFC6C7900U).Radians); // (UInt32)-60000000
            Assert2.Near(
                S1Angle.FromDegrees(75).Radians,
                S1Angle.FromUnsignedE7((UInt32)750000000).Radians);
            Assert2.Near(
                S1Angle.FromDegrees(-75).Radians,
                S1Angle.FromUnsignedE7(0xFFFFFFFFD34BE880U).Radians); // (UInt32)-750000000
        }

        [Fact]
        public void Test_S1Angle_NormalizeCorrectlyCanonicalizesAngles()
        {
            Assert2.Near(0.0, S1Angle.FromDegrees(360.0).Normalize().GetDegrees());
            Assert2.Near(-90.0, S1Angle.FromDegrees(-90.0).Normalize().GetDegrees());
            Assert2.Near(180.0, S1Angle.FromDegrees(-180.0).Normalize().GetDegrees());
            Assert2.Near(180.0, S1Angle.FromDegrees(180.0).Normalize().GetDegrees());
            Assert2.Near(180.0, S1Angle.FromDegrees(540.0).Normalize().GetDegrees());
            Assert2.Near(90.0, S1Angle.FromDegrees(-270.0).Normalize().GetDegrees());
        }

        [Fact]
        public void Test_S1Angle_ArithmeticOperationsOnAngles()
        {
            Assert2.Near(0.3, S1Angle.FromRadians(-0.3).Abs());
            Assert2.Near(-0.1, (-S1Angle.FromRadians(0.1)).Radians);
            Assert2.Near(0.4, (S1Angle.FromRadians(0.1) + S1Angle.FromRadians(0.3)).Radians);
            Assert2.Near(-0.2, (S1Angle.FromRadians(0.1) - S1Angle.FromRadians(0.3)).Radians);
            Assert2.Near(0.6, (2 * S1Angle.FromRadians(0.3)).Radians);
            Assert2.Near(0.6, (S1Angle.FromRadians(0.3) * 2).Radians);
            Assert2.Near(0.15, (S1Angle.FromRadians(0.3) / 2).Radians);
            Assert2.Near(0.5, (S1Angle.FromRadians(0.3) / S1Angle.FromRadians(0.6)));

            S1Angle tmp = S1Angle.FromRadians(1.0);
            tmp += S1Angle.FromRadians(0.5);
            Assert2.Near(1.5, tmp.Radians);
            tmp -= S1Angle.FromRadians(1.0);
            Assert2.Near(0.5, tmp.Radians);
            tmp *= 5;
            Assert2.Near(2.5, tmp.Radians);
            tmp /= 2;
            Assert2.Near(1.25, tmp.Radians);
        }

        [Fact]
        public void Test_S1Angle_Trigonometry()
        {
            // Spot check a few angles to ensure that the correct function is called.
            Assert2.Near(1, S1Angle.FromDegrees(0).Cos());
            Assert2.Near(1, S1Angle.FromDegrees(90).Sin());
            Assert2.Near(1, S1Angle.FromDegrees(45).Tan());
        }

        [Fact]
        public void Test_S1Angle_ConstructorsThatMeasureAngles()
        {
            Assert2.Near(S2.M_PI_2, new S1Angle(new S2Point(1, 0, 0), new S2Point(0, 0, 2)).Radians);
            Assert2.Near(0.0, new S1Angle(new S2Point(1, 0, 0), new S2Point(1, 0, 0)).Radians);
            Assert2.Near(50.0, new S1Angle(S2LatLng.FromDegrees(20, 20), S2LatLng.FromDegrees(70, 20)).GetDegrees(), 1e-13);
        }

        [Fact]
        public void Test_S1Angle_TestFormatting()
        {
            Assert.Equal("180d", S1Angle.FromDegrees(180).ToString()); // 180.0000000
        }

        // The current implementation guarantees exact conversions between
        // Degrees() and E6() when the Degrees() argument is an integer.
        [Fact]
        public void Test_S1Angle_DegreesVsE6()
        {
            for (int i = 0; i <= 180; ++i)
            {
                Assert.Equal(S1Angle.FromDegrees(i), S1Angle.FromE6(1000000 * i));
            }
        }

        // The current implementation guarantees exact conversions between
        // Degrees() and E7() when the Degrees() argument is an integer.
        [Fact]
        public void Test_S1Angle_DegreesVsE7()
        {
            for (int i = 0; i <= 180; ++i)
            {
                Assert.Equal(S1Angle.FromDegrees(i), S1Angle.FromE7(10000000 * i));
            }
        }

        // The current implementation guarantees exact conversions between
        // E6() and E7() when the E6() argument is an integer.
        [Fact]
        public void Test_S1Angle_E6VsE7()
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int iter = 0; iter < 1000; ++iter)
            {
                int i = S2Testing.Random.Uniform(180000000);
                Assert.Equal(S1Angle.FromE6(i), S1Angle.FromE7(10 * i));
            }
        }

        // The current implementation guarantees certain exact conversions between
        // degrees and radians (see the header file for details).
        [Fact]
        public void Test_S1Angle_DegreesVsRadians()
        {
            for (int k = -8; k <= 8; ++k)
            {
                Assert.Equal(S1Angle.FromDegrees(45 * k), S1Angle.FromRadians(k * S2.M_PI_4));
                Assert.Equal(45 * k, S1Angle.FromDegrees(45 * k).GetDegrees());
            }
            for (int k = 0; k <= 30; ++k)
            {
                int n = 1 << k;
                Assert.Equal(S1Angle.FromDegrees(180.0 / n), S1Angle.FromRadians(Math.PI / n));
                Assert.Equal(S1Angle.FromDegrees(60.0 / n), S1Angle.FromRadians(Math.PI / (3.0 * n)));
                Assert.Equal(S1Angle.FromDegrees(36.0 / n), S1Angle.FromRadians(Math.PI / (5.0 * n)));
                Assert.Equal(S1Angle.FromDegrees(20.0 / n), S1Angle.FromRadians(Math.PI / (9.0 * n)));
                Assert.Equal(S1Angle.FromDegrees(4.0 / n), S1Angle.FromRadians(Math.PI / (45.0 * n)));
            }
            // We also spot check a couple of non-identities.
            Assert.NotEqual(S1Angle.FromDegrees(3), S1Angle.FromRadians(Math.PI / 60));
            Assert.NotEqual(60, S1Angle.FromDegrees(60).GetDegrees());
        }
    }
}
