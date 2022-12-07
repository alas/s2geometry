namespace S2Geometry;

using static S2Earth;
using Meters = Single;

public class S2EarthTests
{
    [Fact]
    internal void Test_S2EarthTest_TestAngleConversion()
    {
        // Functions that use meters:
        Assert2.DoubleEqual(MetersToAngle(RadiusMeters).Radians, 1);
        Assert2.DoubleEqual(MetersToChordAngle(RadiusMeters).Radians(), 1);
        Assert2.DoubleEqual(MetersToRadians(RadiansToKm(0.3) * 1000), 0.3);
        Assert2.DoubleEqual(ToMeters(S1Angle.FromDegrees(180)), RadiusMeters * S2.M_PI);
        Assert2.DoubleEqual(ToMeters(S1ChordAngle.FromDegrees(180)), RadiusMeters * S2.M_PI);
        Assert2.DoubleEqual(RadiansToMeters(KmToRadians(2.5)), 2500);

        // Functions that use kilometers:
        Assert2.DoubleEqual(KmToAngle(RadiusKm).Radians, 1);
        Assert2.DoubleEqual(KmToChordAngle(RadiusKm).Radians(), 1);
        Assert2.DoubleEqual(KmToRadians(RadiusMeters / 1000), 1);
        Assert2.DoubleEqual(ToKm(S1Angle.FromRadians(0.5)), 0.5 * RadiusKm);
        Assert2.DoubleEqual(ToKm(S1ChordAngle.FromRadians(0.5)), 0.5 * RadiusKm);
        Assert2.DoubleEqual(RadiansToKm(0.5), 0.5 * RadiusKm);

        // Functions that use util::units::Meters (which only has "float" precision,
        // but fortunately Radius() is exactly representable as "float"):
        Assert2.DoubleEqual(ToAngle(Radius).Radians, 1);
        Assert2.DoubleEqual(ToChordAngle(Radius).Radians(), 1);
        Assert2.DoubleEqual(ToRadians(Radius), 1);
        Assert2.FloatEqual(ToDistance(S1Angle.FromRadians(2)), 2 * RadiusMeters);
        Assert2.FloatEqual(ToDistance(S1ChordAngle.FromRadians(0.5)), 0.5 * RadiusMeters);
        Assert2.FloatEqual(RadiansToDistance(1.5), 1.5 * RadiusMeters);
    }

    [Fact]
    internal void Test_S2EarthTest_TestSolidAngleConversion()
    {
        Assert2.DoubleEqual(1, SquareKmToSteradians(Math.Pow(RadiusMeters / 1000, 2)));
        Assert2.DoubleEqual(Math.Pow(0.5 * RadiusKm, 2),
            SteradiansToSquareKm(Math.Pow(0.5, 2)));
        Assert2.DoubleEqual(Math.Pow(0.3, 2),
            SquareMetersToSteradians(Math.Pow(RadiansToKm(0.3) * 1000, 2)));
        Assert2.DoubleEqual(Math.Pow(2500, 2),
            SteradiansToSquareMeters(Math.Pow(KmToRadians(2.5), 2)));
    }

    [Fact]
    internal void Test_S2EarthTest_TestToLongitudeRadians()
    {
        var earth_radius = Radius;

        // At the equator, ToLongitudeRadians behaves exactly like ToRadians.
        Assert2.DoubleEqual(1, ToLongitudeRadians(earth_radius, 0));

        // The closer we get to the poles, the more radians we need to go the same
        // distance.
        Assert.True(
            ToLongitudeRadians(earth_radius, 0.5) >
            ToLongitudeRadians(earth_radius, 0.4));

        // At the poles, we should return 2PI radians instead of dividing by 0.
        Assert2.DoubleEqual(S2.M_2_PI,
            ToLongitudeRadians(earth_radius, S2.M_PI_2));

        // Within epsilon of the poles, we should still return 2PI radians instead
        // of directing the caller to take thousands of radians around.
        Assert2.DoubleEqual(S2.M_2_PI,
            ToLongitudeRadians(earth_radius, S2.M_PI_2 - 1e-4));

        // Check that the "meters" and "kilometer" versions are compatible.
        Assert2.Equal(ToLongitudeRadians(earth_radius, 0.5),
            MetersToLongitudeRadians(earth_radius, 0.5));
        Assert2.DoubleEqual(
            ToLongitudeRadians(earth_radius, 0.5),
            KmToLongitudeRadians(earth_radius / 1000.0, 0.5));
    }

    [Fact]
    internal void Test_S2EarthTest_TestGetInitialBearing()
    {
        List<(string description, S2LatLng a, S2LatLng b, S1Angle bearing)> test_configs = new()
        {
            ("Westward on equator", S2LatLng.FromDegrees(0, 50),
            S2LatLng.FromDegrees(0, 100), S1Angle.FromDegrees(90)),
            ("Eastward on equator", S2LatLng.FromDegrees(0, 50),
            S2LatLng.FromDegrees(0, 0), S1Angle.FromDegrees(-90)),
            ("Northward on meridian", S2LatLng.FromDegrees(16, 28),
            S2LatLng.FromDegrees(81, 28), S1Angle.FromDegrees(0)),
            ("Southward on meridian", S2LatLng.FromDegrees(24, 64),
            S2LatLng.FromDegrees(-27, 64), S1Angle.FromDegrees(180)),
            ("Towards north pole", S2LatLng.FromDegrees(12, 76),
            S2LatLng.FromDegrees(90, 50), S1Angle.FromDegrees(0)),
            ("Towards south pole", S2LatLng.FromDegrees(-35, 105),
            S2LatLng.FromDegrees(-90, -120), S1Angle.FromDegrees(180)),
            ("Spain to Japan", S2LatLng.FromDegrees(40.4379332, -3.749576),
            S2LatLng.FromDegrees(35.6733227, 139.6403486), S1Angle.FromDegrees(29.2)),
            ("Japan to Spain", S2LatLng.FromDegrees(35.6733227, 139.6403486),
            S2LatLng.FromDegrees(40.4379332, -3.749576), S1Angle.FromDegrees(-27.2)),
        };

        foreach (var config in test_configs)
        {
            var bearing = GetInitialBearing(config.a, config.b);
            var angle_diff = S1Angle.Abs((bearing - config.bearing).Normalize());
            Assert.True(angle_diff.GetDegrees() <= 1e-2);
            // $"GetInitialBearing() test failed on: {config.description}." +
            // $" Expected {config.bearing}, got {bearing}";
        }
    }

    [Fact]
    internal void Test_S2EarthTest_TestGetDistance()
    {
        S2Point north = new(0, 0, 1);
        S2Point south = new(0, 0, -1);
        S2Point west = new(0, -1, 0);

        Assert2.FloatEqual((Meters)(S2.M_PI * RadiusMeters), (Meters)GetDistance(north, south));
        Assert2.DoubleEqual(0, GetDistanceKm(west, west));
        Assert2.DoubleEqual(S2.M_PI_2 * RadiusMeters, GetDistanceMeters(north, west));

        Assert2.FloatEqual(GetDistance(west, south), GetDistance(
            S2LatLng.FromDegrees(0, -90),
            S2LatLng.FromDegrees(-90, -38)));

        Assert2.DoubleEqual(RadiusKm, GetDistanceKm(
            S2LatLng.FromRadians(0, 0.6),
            S2LatLng.FromRadians(0, -0.4)));

        Assert2.DoubleEqual(1000 * RadiusKm * S2.M_PI_4, GetDistanceMeters(
            S2LatLng.FromDegrees(80, 27),
            S2LatLng.FromDegrees(55, -153)));
    }
}
