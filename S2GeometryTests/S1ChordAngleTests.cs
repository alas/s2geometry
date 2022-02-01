using System;
using Xunit;

namespace S2Geometry
{
    public class S1ChordAngleTests
    {
        [Fact]
        public void Test_S1ChordAngle_DefaultConstructor()
        {
            // Check that the default constructor returns an angle of 0.
            S1ChordAngle a = S1ChordAngle.Zero;
            Assert.Equal(S1ChordAngle.Zero, a);
        }

        [Fact]
        public void Test_S1ChordAngle_TwoPointConstructor()
        {
            for (int iter = 0; iter < 100; ++iter)
            {
                S2Testing.GetRandomFrame(out var x, out var y, out var z);
                Assert.Equal(S1Angle.Zero, new S1ChordAngle(z, z).ToAngle());
                Assert2.Near(Math.PI, new S1ChordAngle(-z, z).Radians(), 1e-7);
                Assert2.Near(S2.M_PI_2, new S1ChordAngle(x, z).Radians(), 0.03);
                S2Point w = (y + z).Normalize();
                Assert2.Near(S2.M_PI_4, new S1ChordAngle(w, z).Radians());
            }
        }

        [Fact]
        public void Test_S1ChordAngle_FromLength2()
        {
            Assert.Equal(0, S1ChordAngle.FromLength2(0).Degrees());
            Assert2.Near(60, S1ChordAngle.FromLength2(1).Degrees());
            Assert2.Near(90, S1ChordAngle.FromLength2(2).Degrees());
            Assert.Equal(180, S1ChordAngle.FromLength2(4).Degrees());
            Assert.Equal(180, S1ChordAngle.FromLength2(5).Degrees());
        }

        [Fact]
        public void Test_S1ChordAngle_Zero()
        {
            Assert.Equal(S1Angle.Zero, S1ChordAngle.Zero.ToAngle());
        }

        [Fact]
        public void Test_S1ChordAngle_Right()
        {
            Assert2.Near(90, S1ChordAngle.Right.Degrees());
        }

        [Fact]
        public void Test_S1ChordAngle_Straight()
        {
            Assert.True(S1Angle.FromDegrees(180) == S1ChordAngle.Straight.ToAngle());
        }

        [Fact]
        public void Test_S1ChordAngle_Infinity()
        {
            Assert.True(S1ChordAngle.Straight < S1ChordAngle.Infinity);
            Assert.True(S1ChordAngle.Infinity == new S1ChordAngle(S1Angle.Infinity));
            Assert.True(S1Angle.Infinity == S1ChordAngle.Infinity.ToAngle());
        }

        [Fact]
        public void Test_S1ChordAngle_Negative()
        {
            Assert.True(S1ChordAngle.Negative < S1ChordAngle.Zero);
            Assert.True(S1ChordAngle.Negative == S1ChordAngle.FromLength2(-1));
            Assert.True(S1ChordAngle.Negative.ToAngle() < S1Angle.Zero);
        }

        [Fact]
        public void Test_S1ChordAngle_Predicates()
        {
            Assert.True(S1ChordAngle.Zero.IsZero());
            Assert.False(S1ChordAngle.Zero.IsNegative());
            Assert.False(S1ChordAngle.Zero.IsSpecial());
            Assert.False(S1ChordAngle.Straight.IsSpecial());
            Assert.True(S1ChordAngle.Negative.IsNegative());
            Assert.True(S1ChordAngle.Negative.IsSpecial());
            Assert.True(S1ChordAngle.Infinity.IsInfinity());
            Assert.True(S1ChordAngle.Infinity.IsSpecial());
        }

        [Fact]
        public void Test_S1ChordAngle_ToFromS1Angle()
        {
            Assert.Equal(0, new S1ChordAngle(S1Angle.Zero).Radians());
            Assert.Equal(4, new S1ChordAngle(S1Angle.FromRadians(Math.PI)).Length2);
            Assert.Equal(Math.PI, new S1ChordAngle(S1Angle.FromRadians(Math.PI)).Radians());
            Assert.Equal(S1Angle.Infinity, new S1ChordAngle(S1Angle.Infinity).ToAngle());
            Assert.True(new S1ChordAngle(S1Angle.FromRadians(-1)).Radians() < 0);
            Assert2.Near(1.0, new S1ChordAngle(S1Angle.FromRadians(1.0)).Radians());
        }

        [Fact]
        public void Test_S1ChordAngle_Successor()
        {
            Assert.Equal(S1ChordAngle.Zero, S1ChordAngle.Negative.Successor());
            Assert.Equal(S1ChordAngle.Infinity, S1ChordAngle.Straight.Successor());
            Assert.Equal(S1ChordAngle.Infinity, S1ChordAngle.Infinity.Successor());
            S1ChordAngle x = S1ChordAngle.Negative;
            for (int i = 0; i < 10; ++i)
            {
                Assert.True(x < x.Successor());
                x = x.Successor();
            }
        }

        [Fact]
        public void Test_S1ChordAngle_Predecessor()
        {
            Assert.Equal(S1ChordAngle.Straight, S1ChordAngle.Infinity.Predecessor());
            Assert.Equal(S1ChordAngle.Negative, S1ChordAngle.Zero.Predecessor());
            Assert.Equal(S1ChordAngle.Negative, S1ChordAngle.Negative.Predecessor());
            S1ChordAngle x = S1ChordAngle.Infinity;
            for (int i = 0; i < 10; ++i)
            {
                Assert.True(x > x.Predecessor());
                x = x.Predecessor();
            }
        }

        [Fact]
        public void Test_S1ChordAngle_Arithmetic()
        {
            S1ChordAngle zero = S1ChordAngle.Zero;
            S1ChordAngle degree30 = S1ChordAngle.FromDegrees(30);
            S1ChordAngle degree60 = S1ChordAngle.FromDegrees(60);
            S1ChordAngle degree90 = S1ChordAngle.FromDegrees(90);
            S1ChordAngle degree120 = S1ChordAngle.FromDegrees(120);
            S1ChordAngle degree180 = S1ChordAngle.Straight;
            Assert.Equal(0, (zero + zero).Degrees());
            Assert.Equal(0, (zero - zero).Degrees());
            Assert.Equal(0, (degree60 - degree60).Degrees());
            Assert.Equal(0, (degree180 - degree180).Degrees());
            Assert.Equal(0, (zero - degree60).Degrees());
            Assert.Equal(0, (degree30 - degree90).Degrees());
            Assert2.Near(60, (degree60 + zero).Degrees());
            Assert2.Near(60, (degree60 - zero).Degrees());
            Assert2.Near(60, (zero + degree60).Degrees());
            Assert2.Near(90, (degree30 + degree60).Degrees());
            Assert2.Near(90, (degree60 + degree30).Degrees());
            Assert2.Near(60, (degree90 - degree30).Degrees());
            Assert2.Near(30, (degree90 - degree60).Degrees());
            Assert.Equal(180, (degree180 + zero).Degrees());
            Assert.Equal(180, (degree180 - zero).Degrees());
            Assert.Equal(180, (degree90 + degree90).Degrees());
            Assert.Equal(180, (degree120 + degree90).Degrees());
            Assert.Equal(180, (degree120 + degree120).Degrees());
            Assert.Equal(180, (degree30 + degree180).Degrees());
            Assert.Equal(180, (degree180 + degree180).Degrees());
        }

        [Fact]
        public void Test_S1ChordAngle_ArithmeticPrecision()
        {
            // Verifies that S1ChordAngle is capable of adding and subtracting angles
            // extremely accurately up to Pi/2 radians.  (Accuracy continues to be good
            // well beyond this value but degrades as angles approach Pi.)
            S1ChordAngle kEps = S1ChordAngle.FromRadians(1e-15);
            S1ChordAngle k90 = S1ChordAngle.Right;
            S1ChordAngle k90MinusEps = k90 - kEps;
            S1ChordAngle k90PlusEps = k90 + kEps;
            double kMaxError = 2 * S2.DoubleEpsilon;
            Assert2.Near(k90MinusEps.Radians(), S2.M_PI_2 - kEps.Radians(), kMaxError);
            Assert2.Near(k90PlusEps.Radians(), S2.M_PI_2 + kEps.Radians(), kMaxError);
            Assert2.Near((k90 - k90MinusEps).Radians(), kEps.Radians(), kMaxError);
            Assert2.Near((k90PlusEps - k90).Radians(), kEps.Radians(), kMaxError);
            Assert2.Near((k90MinusEps + kEps).Radians(), S2.M_PI_2, kMaxError);
        }

        [Fact]
        public void Test_S1ChordAngle_Trigonometry()
        {
            const int kIters = 20;
            for (int iter = 0; iter <= kIters; ++iter)
            {
                double radians = Math.PI * iter / kIters;
                S1ChordAngle angle = new(S1Angle.FromRadians(radians));
                Assert2.Near(Math.Sin(radians), angle.Sin(), S2.DoubleError);
                Assert2.Near(Math.Cos(radians), angle.Cos(), S2.DoubleError);
                // Since the tan(x) is unbounded near Pi/4, we map the result back to an
                // angle before comparing.  (The assertion is that the result is equal to
                // the tangent of a nearby angle.)
                Assert2.Near(Math.Atan(Math.Tan(radians)), Math.Atan(angle.Tan()), S2.DoubleError);
            }

            // Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
            S1ChordAngle angle90 = S1ChordAngle.FromLength2(2);
            S1ChordAngle angle180 = S1ChordAngle.FromLength2(4);
            Assert.Equal(1, angle90.Sin());
            Assert.Equal(0, angle90.Cos());
            Assert.Equal(double.PositiveInfinity, angle90.Tan());
            Assert.Equal(0, angle180.Sin());
            Assert.Equal(-1, angle180.Cos());
            Assert.Equal(0, angle180.Tan());
        }

        [Fact]
        public void Test_S1ChordAngle_PlusError()
        {
            Assert.Equal(S1ChordAngle.Negative, S1ChordAngle.Negative.PlusError(5));
            Assert.Equal(S1ChordAngle.Infinity, S1ChordAngle.Infinity.PlusError(-5));
            Assert.Equal(S1ChordAngle.Straight, S1ChordAngle.Straight.PlusError(5));
            Assert.Equal(S1ChordAngle.Zero, S1ChordAngle.Zero.PlusError(-5));
            Assert.Equal(S1ChordAngle.FromLength2(1.25),
                      S1ChordAngle.FromLength2(1).PlusError(0.25));
            Assert.Equal(S1ChordAngle.FromLength2(0.75),
                      S1ChordAngle.FromLength2(1).PlusError(-0.25));
        }

        [Fact]
        public void Test_S1ChordAngle_GetS2PointConstructorMaxError()
        {
            // Check that the error bound returned by GetS2PointConstructorMaxError() is
            // large enough.
            for (var iter = 0; iter < 100000; ++iter)
            {
                S2Testing.Random.Reset(iter);  // Easier to reproduce a specific case.
                var x = S2Testing.RandomPoint();
                var y = S2Testing.RandomPoint();
                if (S2Testing.Random.OneIn(10))
                {
                    // Occasionally test a point pair that is nearly identical or antipodal.
                    S1Angle r = S1Angle.FromRadians(S2.DoubleError * S2Testing.Random.RandDouble());
                    y = S2.InterpolateAtDistance(r, x, y);
                    if (S2Testing.Random.OneIn(2)) y = -y;
                }
                S1ChordAngle dist = new(x, y);
                double error = dist.GetS2PointConstructorMaxError();
                var er1 = S2Pred.CompareDistance(x, y, dist.PlusError(error));
                if (er1 > 0)
                {

                }
                Assert.True(er1 <= 0);
                var er2 = S2Pred.CompareDistance(x, y, dist.PlusError(-error));
                if (er1 < 0)
                {

                }
                Assert.True(er2 >= 0);
            }
        }
    }
}
