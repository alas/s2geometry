// The earth modeled as a sphere.  There are lots of convenience functions so
// that it doesn't take 2 lines of code just to do a single conversion.

using Meters = System.Single;

namespace S2Geometry;

public static class S2Earth
{
    // These functions convert between distances on the unit sphere
    // (expressed as angles subtended from the sphere's center) and
    // distances on the Earth's surface.  This is possible only because
    // the Earth is modeled as a sphere; otherwise a given angle would
    // correspond to a range of distances depending on where the
    // corresponding line segment was located.

    // Functions for converting distances to angles:
    public static S1Angle MetersToAngle(double meters) => S1Angle.FromRadians(MetersToRadians(meters));
    public static /*inline*/ S1ChordAngle MetersToChordAngle(double meters) => new(MetersToAngle(meters));
    public static double MetersToRadians(double meters) => meters / RadiusMeters;

    // Functions for converting angles to distances:
    public static double ToMeters(S1Angle angle) => angle.Radians * RadiusMeters;
    public static /*inline*/ double ToMeters(S1ChordAngle cangle) => ToMeters(cangle.ToAngle());
    public static double RadiansToMeters(double radians) => radians * RadiusMeters;

    // Like the above, but where distances are expressed in kilometers:
    public static S1Angle KmToAngle(double km) => S1Angle.FromRadians(KmToRadians(km));
    public static /*inline*/ S1ChordAngle KmToChordAngle(double km) => new(KmToAngle(km));
    public static double KmToRadians(double km) => km / RadiusKm;
    public static double ToKm(S1Angle angle) => angle.Radians * RadiusKm;
    public static /*inline*/ double ToKm(S1ChordAngle cangle) => ToKm(cangle.ToAngle());
    public static double RadiansToKm(double radians) => radians * RadiusKm;

    // Like the above, but where distances are expressed in util::units::Meters.
    //
    // CAVEAT: These versions are not as accurate because util::units::Meters
    // uses "float" rather than "double" as the underlying representation.
    // (You may lose precision if you use these functions.)
    public static S1Angle ToAngle(Meters distance) => S1Angle.FromRadians(ToRadians(distance));
    public static /*inline*/ S1ChordAngle ToChordAngle(Meters distance) => new(ToAngle(distance));
    public static double ToRadians(Meters distance) => distance / RadiusMeters;
    public static Meters ToDistance(S1Angle angle) => (Meters)ToMeters(angle);
    public static /*inline*/ Meters ToDistance(S1ChordAngle cangle) => (Meters)ToMeters(cangle);
    public static Meters RadiansToDistance(double radians) => (Meters)RadiansToMeters(radians);

    // These functions convert between areas on the unit sphere
    // (as returned by the S2 library) and areas on the Earth's surface.
    // Note that the area of a region on the unit sphere is equal to the
    // solid angle it subtends from the sphere's center (measured in steradians).
    public static double SquareKmToSteradians(double km2) => km2 / (RadiusKm * RadiusKm);
    public static double SquareMetersToSteradians(double m2) => m2 / (RadiusMeters * RadiusMeters);
    public static double SteradiansToSquareKm(double steradians) => steradians * RadiusKm * RadiusKm;
    public static double SteradiansToSquareMeters(double steradians) => steradians * RadiusMeters * RadiusMeters;

    // Convenience functions for the frequent case where you need to call
    // ToRadians in order to convert an east-west distance on the globe to
    // radians. The output is a function of how close to the poles you are
    // (i.e. at the bulge at the equator, one unit of longitude represents a
    // much farther distance). The function will never return more than 2*PI
    // radians, even if you're trying to go 100 million miles west at the north
    // pole.
    public static double MetersToLongitudeRadians(double meters, double latitude_radians)
    {
        var scalar = Math.Cos(latitude_radians);
        if (scalar == 0) return S2.M_2_PI;
        return Math.Min(MetersToRadians(meters) / scalar, S2.M_2_PI);
    }
    public static double KmToLongitudeRadians(double km, double latitude_radians) => MetersToLongitudeRadians(1000 * km, latitude_radians);

    // CAVEAT: This version is not as accurate because util::units::Meters
    // uses "float" rather than "double" as the underlying representation.
    public static double ToLongitudeRadians(Meters distance, double latitude_radians) => MetersToLongitudeRadians(distance, latitude_radians);

    // Computes the initial bearing from a to b. This is the bearing an observer
    // at point a has when facing point b. A bearing of 0 degrees is north, and it
    // increases clockwise (90 degrees is east, etc).
    //
    // If a == b, a == -b, or a is one of the Earths' poles, the return value is
    // undefined.
    //
    // Sourced from http://www.movable-type.co.uk/scripts/latlong.html.
    public static S1Angle GetInitialBearing(S2LatLng a, S2LatLng b)
    {
        var lat1 = a.Lat().Radians;
        var cosLat2 = Math.Cos(b.Lat().Radians);
        var lat_diff = b.Lat().Radians - a.Lat().Radians;
        var lng_diff = b.Lng().Radians - a.Lng().Radians;

        var x = Math.Sin(lat_diff) + Math.Sin(lat1) * cosLat2 * 2 * Haversine(lng_diff);
        var y = Math.Sin(lng_diff) * cosLat2;
        return S1Angle.FromRadians(Math.Atan2(y, x));
    }

    // Returns the distance between two points.
    public static /*inline*/ double GetDistanceMeters(S2Point a, S2Point b) => RadiansToMeters(a.Angle(b));
    public static /*inline*/ double GetDistanceMeters(S2LatLng a, S2LatLng b) => ToMeters(a.GetDistance(b));
    public static /*inline*/ double GetDistanceKm(S2Point a, S2Point b) => RadiansToKm(a.Angle(b));
    public static /*inline*/ double GetDistanceKm(S2LatLng a, S2LatLng b) => ToKm(a.GetDistance(b));

    // CAVEAT: These versions are not as accurate because util::units::Meters
    // uses "float" rather than "double" as the underlying representation.
    public static /*inline*/ Meters GetDistance(S2Point a, S2Point b) => RadiansToDistance(a.Angle(b));
    public static /*inline*/ Meters GetDistance(S2LatLng a, S2LatLng b) => ToDistance(a.GetDistance(b));

    // Returns the Earth's mean radius, which is the radius of the equivalent
    // sphere with the same surface area.  According to NASA, this value is
    // 6371.01 +/- 0.02 km.  The equatorial radius is 6378.136 km, and the polar
    // radius is 6356.752 km.  They differ by one part in 298.257.
    //
    // Reference: http://ssd.jpl.nasa.gov/phys_props_earth.html, which quotes
    // Yoder, C.F. 1995. "Astrometric and Geodetic Properties of Earth and the
    // Solar System" in Global Earth Physics, A Handbook of Physical Constants,
    // AGU Reference Shelf 1, American Geophysical Union, Table 2.
    public const double RadiusMeters = 6371010.0;
    public const double RadiusKm = 0.001 * RadiusMeters;

    // CAVEAT: This version is not as accurate because util::units::Meters
    // uses "float" rather than "double" as the underlying representation.
    public const Meters Radius = (Meters)RadiusMeters;

    // Convenience functions.

    // Returns the altitude of the lowest known point on Earth. The lowest known
    // point on Earth is the Challenger Deep with an altitude of -10898 meters
    // above the surface of the spherical earth.
    public const double LowestAltitudeMeters = -10898;
    public const double LowestAltitudeKm = 0.001 * LowestAltitudeMeters;
    public const Meters LowestAltitude = (Meters)LowestAltitudeMeters;

    // Returns the altitude of the highest known point on Earth. The highest
    // known point on Earth is Mount Everest with an altitude of 8846 meters
    // above the surface of the spherical earth.
    public const double HighestAltitudeMeters = 8846;
    public const double HighestAltitudeKm = 0.001 * HighestAltitudeMeters;
    public const Meters HighestAltitude = (Meters)HighestAltitudeMeters;

    // http://en.wikipedia.org/wiki/Haversine_formula
    // Haversine(x) has very good numerical stability around zero.
    // Haversine(x) == (1-cos(x))/2 == sin(x/2)^2; must be implemented with the
    // second form to reap the numerical benefits.
    private static double Haversine(double radians)
    {
        var sinHalf = Math.Sin(radians / 2);
        return sinHalf * sinHalf;
    }
}
