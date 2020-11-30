using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;
using Xunit.Abstractions;

namespace S2Geometry
{
    public class S2LatLngTests
    {
        private readonly ITestOutputHelper _logger;

        public S2LatLngTests(ITestOutputHelper logger) { _logger = logger; }

        [Fact]
        public void Test_S2LatLng_TestBasic()
        {
            S2LatLng ll_rad = S2LatLng.FromRadians(S2Constants.M_PI_4, S2Constants.M_PI_2);
            Assert.Equal(S2Constants.M_PI_4, ll_rad.LatRadians);
            Assert.Equal(S2Constants.M_PI_2, ll_rad.LngRadians);
            Assert.True(ll_rad.IsValid);
            S2LatLng ll_deg = S2LatLng.FromDegrees(45, 90);
            Assert.Equal(ll_rad, ll_deg);
            Assert.True(ll_deg.IsValid);
            Assert.False(S2LatLng.FromDegrees(-91, 0).IsValid);
            Assert.False(S2LatLng.FromDegrees(0, 181).IsValid);

            S2LatLng bad = S2LatLng.FromDegrees(120, 200);
            Assert.False(bad.IsValid);
            S2LatLng better = bad.Normalized;
            Assert.True(better.IsValid);
            Assert.Equal(S1Angle.FromDegrees(90), better.Lat);
            Assert2.Near(S1Angle.FromDegrees(-160).Radians, better.LngRadians);

            bad = S2LatLng.FromDegrees(-100, -360);
            Assert.False(bad.IsValid);
            better = bad.Normalized;
            Assert.True(better.IsValid);
            Assert.Equal(S1Angle.FromDegrees(-90), better.Lat);
            Assert2.Near(0.0, better.LngRadians);

            Assert.True((S2LatLng.FromDegrees(10, 20) + S2LatLng.FromDegrees(20, 30)).
                        ApproxEquals(S2LatLng.FromDegrees(30, 50)));
            Assert.True((S2LatLng.FromDegrees(10, 20) - S2LatLng.FromDegrees(20, 30)).
                        ApproxEquals(S2LatLng.FromDegrees(-10, -10)));
            Assert.True((0.5 * S2LatLng.FromDegrees(10, 20)).
                        ApproxEquals(S2LatLng.FromDegrees(5, 10)));

            // Check that Invalid() returns an invalid point.
            S2LatLng invalid = S2LatLng.Invalid;
            Assert.False(invalid.IsValid);

            // Check that the defaultructor sets latitude and longitude to 0.
            S2LatLng default_ll = S2LatLng.Center;
            Assert.True(default_ll.IsValid);
            Assert.Equal(0, default_ll.LatRadians);
            Assert.Equal(0, default_ll.LngRadians);
        }

        [Fact]
        public void Test_S2LatLng_TestConversion()
        {
            // Test special cases: poles, "date line"
            Assert2.Near(90.0, new S2LatLng(S2LatLng.FromDegrees(90.0, 65.0).ToPoint()).Lat.Degrees);
            Assert.Equal(-S2Constants.M_PI_2, new S2LatLng(S2LatLng.FromRadians(-S2Constants.M_PI_2, 1).ToPoint()).LatRadians);
            Assert2.Near(180.0, Math.Abs(new S2LatLng(S2LatLng.FromDegrees(12.2, 180.0).ToPoint()).Lng.Degrees));
            Assert.Equal(Math.PI, Math.Abs(new S2LatLng(S2LatLng.FromRadians(0.1, -Math.PI).ToPoint()).LngRadians));

            // Test a bunch of random points.
            for (int i = 0; i < 100000; ++i)
            {
                S2Point p = S2Testing.RandomPoint();
                Assert.True(S2PointUtil.ApproxEquals(p, new S2LatLng(p).ToPoint()));
            }
        }

        [Fact]
        public void Test_S2LatLng_TestDistance()
        {
            Assert.Equal(0.0, S2LatLng.FromDegrees(90, 0).GetDistance(S2LatLng.FromDegrees(90, 0)).Radians);
            Assert2.Near(77.0, S2LatLng.FromDegrees(-37, 25).GetDistance(S2LatLng.FromDegrees(-66, -155)).Degrees, 1e-13);
            Assert2.Near(115.0, S2LatLng.FromDegrees(0, 165).GetDistance(S2LatLng.FromDegrees(0, -80)).Degrees, 1e-13);
            Assert2.Near(180.0, S2LatLng.FromDegrees(47, -127).GetDistance(S2LatLng.FromDegrees(-47, 53)).Degrees, 2e-6);
        }

        [Fact]
        public void Test_S2LatLng_TestToString()
        {
            var values = new (double lat, double lng, double expected_lat, double expected_lng)[]
                {
                    (0, 0, 0, 0),
                    (1.5, 91.7, 1.5, 91.7),
                    (9.9, -0.31, 9.9, -0.31),
                    (Math.Sqrt(2), -Math.Sqrt(5), 1.414214, -2.236068),
                    (91.3, 190.4, 90, -169.6),
                    (-100, -710, -90, 10),
                };
            int i = 0;
            foreach (var v in values)
            {
                _logger.WriteLine("Iteration " + i++);
                S2LatLng p = S2LatLng.FromDegrees(v.lat, v.lng);
                string output = p.ToStringDegrees();

                var splitted = output.Split(',').Select(t => Convert.ToDouble(t)).ToList();
                Assert.Equal(2, splitted.Count);
                var lat = splitted[0];
                var lng = splitted[1];
                Assert2.Near(v.expected_lat, lat, 1e-8);
                Assert2.Near(v.expected_lng, lng, 1e-8);
            }
        }

        // Test the variant that returns a string.
        [Fact]
        public void Test_S2LatLng_TestToStringReturnsString()
        {
            var s = S2LatLng.FromDegrees(0, 1).ToStringDegrees();
            Assert.Equal(S2LatLng.FromDegrees(0, 1).ToStringDegrees(), s);
        }

        [Fact]
        public void Test_S2LatLng_TestHashCode()
        {
            var map = new Dictionary<S2LatLng, int>
                {
                    [S2LatLng.FromDegrees(0, 10)] = 1,
                    [S2LatLng.FromDegrees(2, 12)] = 2,
                    [S2LatLng.FromDegrees(5, 15)] = 3,
                    [S2LatLng.FromDegrees(7, 17)] = 4,
                    [S2LatLng.FromDegrees(11, 19)] = 5
                };
            Assert.Equal(5, map.Count);
            Assert.Equal(1, map[S2LatLng.FromDegrees(0, 10)]);
            Assert.Equal(2, map[S2LatLng.FromDegrees(2, 12)]);
            Assert.Equal(3, map[S2LatLng.FromDegrees(5, 15)]);
            Assert.Equal(4, map[S2LatLng.FromDegrees(7, 17)]);
            Assert.Equal(5, map[S2LatLng.FromDegrees(11, 19)]);
        }
    }
}
