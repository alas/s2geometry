using System;
using Xunit;

namespace S2Geometry
{
    public class S2EarthTests
    {
        [Fact]
        public void Test_S2EarthTest_TestAngleConversion()
        {
            Assert.Equal(1, S2Earth.ToAngle(S2Earth.RadiusMeters).Radians);
            Assert.Equal(1, S2Earth.ToChordAngle(S2Earth.RadiusMeters).Radians);
            Assert.Equal(2 * S2Earth.RadiusMeters, S2Earth.ToDistance(S1Angle.FromRadians(2)));
            Assert.Equal(2 * S2Earth.RadiusMeters, S2Earth.ToDistance(S1ChordAngle.FromRadians(2)));
            Assert.Equal(1, S2Earth.ToRadians(S2Earth.RadiusMeters));
            Assert.Equal(S2Earth.ToMeters(S1Angle.FromDegrees(180)), S2Earth.RadiusMeters * Math.PI);
            Assert.Equal(S2Earth.ToMeters(S1ChordAngle.FromDegrees(180)), S2Earth.RadiusMeters * Math.PI);
            Assert.Equal(S2Earth.ToKm(S1Angle.FromRadians(0.5)), 0.5 * S2Earth.RadiusKm);
            Assert.Equal(S2Earth.ToKm(S1ChordAngle.FromRadians(0.5)), 0.5 * S2Earth.RadiusKm);
            Assert.Equal(1, S2Earth.KmToRadians(S2Earth.RadiusMeters / 1000));
            Assert.Equal(0.5 * S2Earth.RadiusKm, S2Earth.RadiansToKm(0.5));
            Assert.Equal(0.3, S2Earth.MetersToRadians(S2Earth.RadiansToKm(0.3) * 1000));
            Assert.Equal(2500, S2Earth.RadiansToMeters(S2Earth.KmToRadians(2.5)));
        }

        [Fact]
        public void Test_S2EarthTest_TestSolidAngleConversion()
        {
            Assert.Equal(1, S2Earth.SquareKmToSteradians(Math.Pow(S2Earth.RadiusMeters / 1000, 2)));
            Assert.Equal(Math.Pow(0.5 * S2Earth.RadiusKm, 2), S2Earth.SteradiansToSquareKm(Math.Pow(0.5, 2)));
            Assert.Equal(Math.Pow(0.3, 2), S2Earth.SquareMetersToSteradians(Math.Pow(S2Earth.RadiansToKm(0.3) * 1000, 2)));
            Assert2.Equal(Math.Pow(2500, 2),
                S2Earth.SteradiansToSquareMeters(Math.Pow(S2Earth.KmToRadians(2.5), 2)));
        }

        [Fact]
        public void Test_S2EarthTest_TestToLongitudeRadians()
        {
            var earth_radius = S2Earth.RadiusMeters;

            // At the equator, ToLongitudeRadians behaves exactly like ToRadians.
            Assert.Equal(1, S2Earth.ToLongitudeRadians(earth_radius, 0));

            // The closer we get to the poles, the more radians we need to go the same
            // distance.
            Assert.True(S2Earth.ToLongitudeRadians(earth_radius, 0.5) > S2Earth.ToLongitudeRadians(earth_radius, 0.4));

            // At the poles, we should return 2PI radians instead of dividing by 0.
            Assert.Equal(S2Earth.ToLongitudeRadians(earth_radius, S2Constants.M_PI_2), S2Constants.M_2_PI);

            // Within epsilon of the poles, we should still return 2PI radians instead
            // of directing the caller to take thousands of radians around.
            Assert.Equal(S2Earth.ToLongitudeRadians(earth_radius, S2Constants.M_PI_2 - 1e-4), S2Constants.M_2_PI);
        }

        [Fact]
        public void Test_S2EarthTest_TestGetInitialBearing()
        {
            var test_configs = new (string description, S2LatLng a, S2LatLng b, S1Angle bearing)[]{
      ("Westward on equator", S2LatLng.FromDegrees(0, 50), S2LatLng.FromDegrees(0, 100), S1Angle.FromDegrees(90)),
      ("Eastward on equator", S2LatLng.FromDegrees(0, 50), S2LatLng.FromDegrees(0, 0), S1Angle.FromDegrees(-90)),
      ("Northward on meridian", S2LatLng.FromDegrees(16, 28), S2LatLng.FromDegrees(81, 28), S1Angle.FromDegrees(0)),
      ("Southward on meridian", S2LatLng.FromDegrees(24, 64), S2LatLng.FromDegrees(-27, 64), S1Angle.FromDegrees(180)),
      ("Towards north pole", S2LatLng.FromDegrees(12, 76), S2LatLng.FromDegrees(90, 50), S1Angle.FromDegrees(0)),
      ("Towards south pole", S2LatLng.FromDegrees(-35, 105), S2LatLng.FromDegrees(-90, -120), S1Angle.FromDegrees(180)),
      ("Spain to Japan", S2LatLng.FromDegrees(40.4379332, -3.749576), S2LatLng.FromDegrees(35.6733227, 139.6403486), S1Angle.FromDegrees(29.2)),
      ("Japan to Spain", S2LatLng.FromDegrees(35.6733227, 139.6403486), S2LatLng.FromDegrees(40.4379332, -3.749576), S1Angle.FromDegrees(-27.2)),
  };

            foreach (var config in test_configs)
            {
                S1Angle bearing = S2Earth.GetInitialBearing(config.a, config.b);
                S1Angle angle_diff = S1Angle.FromRadians((bearing - config.bearing).Normalize().Abs);
                Assert.True(angle_diff.Degrees <= 1e-2); // $"GetInitialBearing() test failed on: {config.description}. Expected {config.bearing}, got {bearing}";
            }
        }

        [Fact]
        public void Test_S2EarthTest_TestGetDistance()
        {
            S2Point north = new S2Point(0, 0, 1);
            S2Point south = new S2Point(0, 0, -1);
            S2Point west = new S2Point(0, -1, 0);

            Assert.Equal(Math.PI * S2Earth.RadiusMeters, S2Earth.GetDistance(north, south));
            Assert.Equal(0, S2Earth.GetDistanceKm(west, west));
            Assert.Equal(S2Constants.M_PI_2 * S2Earth.RadiusMeters, S2Earth.GetDistanceMeters(north, west));
            Assert.Equal(S2Earth.GetDistance(S2LatLng.FromDegrees(0, -90), S2LatLng.FromDegrees(-90, -38)), S2Earth.GetDistance(west, south));
            Assert.Equal(S2Earth.GetDistanceKm(S2LatLng.FromRadians(0, 0.6), S2LatLng.FromRadians(0, -0.4)), S2Earth.RadiusKm);
            Assert2.Equal(
                S2Earth.GetDistanceMeters(S2LatLng.FromDegrees(80, 27), S2LatLng.FromDegrees(55, -153)),
                1000 * S2Earth.RadiusKm * S2Constants.M_PI_4);
        }
    }
}
