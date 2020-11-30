using System;

namespace S2Geometry
{
    // The earth modeled as a sphere.  There are lots of convenience
    // functions so that it doesn't take 2 lines of code just to do
    // a single conversion.
    public static class S2Earth
    {
        /// <summary>
        /// These functions convert between distances on the unit sphere
        /// (expressed as angles subtended from the sphere's center) and
        /// distances on the Earth's surface.  This is possible only because
        /// the Earth is modeled as a sphere; otherwise a given angle would
        /// correspond to a range of distances depending on where the
        /// corresponding line segment was located.
        ///
        /// Note that you will lose precision if you use the ToDistance() method,
        /// since Meters is a single-precision type.  If you need more precision,
        /// use one of the direct conversion methods below.
        /// </summary>
        /// <param name="distance">in meters</param>
        public static S1Angle ToAngle(double distance) => S1Angle.FromRadians(ToRadians(distance));

        /// <param name="distance">in meters</param>
        public static S1ChordAngle ToChordAngle(double distance) => new S1ChordAngle(ToAngle(distance));
        /// <returns>distance in meters</returns>
        public static double ToDistance(S1Angle angle) => ToMeters(angle);
        public static double ToDistance(S1ChordAngle cangle) => ToMeters(cangle);

        /// <summary>
        /// Convenience functions.  These methods also return a double-precision
        /// result, unlike the generic ToDistance() method.
        /// </summary>
        /// <param name="distance">in meters</param>
        /// <returns></returns>
        public static double ToRadians(double distance) => distance / RadiusMeters;
        public static double ToMeters(S1Angle angle) => angle.Radians * RadiusMeters;
        public static double ToMeters(S1ChordAngle cangle) => ToMeters(cangle.ToAngle());
        public static double ToKm(S1Angle angle) => angle.Radians * RadiusKm;
        public static double ToKm(S1ChordAngle cangle) => ToKm(cangle.ToAngle());
        public static double KmToRadians(double km) => km / RadiusKm;
        public static double RadiansToKm(double radians) => radians * RadiusKm;
        public static double MetersToRadians(double meters) => meters / RadiusMeters;
        public static double RadiansToMeters(double radians) => radians * RadiusMeters;

        // These functions convert between areas on the unit sphere
        // (as returned by the S2 library) and areas on the Earth's surface.
        // Note that the area of a region on the unit sphere is equal to the
        // solid angle it subtends from the sphere's center (measured in steradians).
        public static double SquareKmToSteradians(double km2) => km2 / (RadiusKm * RadiusKm);
        public static double SquareMetersToSteradians(double m2) => m2 / (RadiusMeters * RadiusMeters);
        public static double SteradiansToSquareKm(double steradians) => steradians * RadiusKm * RadiusKm;
        public static double SteradiansToSquareMeters(double steradians) => steradians * RadiusMeters * RadiusMeters;

        /// <summary>
        /// Convenience function for the frequent case where you need to call
        /// ToRadians in order to convert an east-west distance on the globe to
        /// radians. The output is a function of how close to the poles you are
        /// (i.e. at the bulge at the equator, one unit of longitude represents a
        /// much farther distance). The function will never return more than 2*PI
        /// radians, even if you're trying to go 100 million miles west at the north
        /// pole.
        /// </summary>
        /// <param name="distance">in meters</param>
        /// <param name="latitude_radians"></param>
        /// <returns></returns>
        public static double ToLongitudeRadians(double distance, double latitude_radians)
        {
            var scalar = Math.Cos(latitude_radians);
            if (scalar == 0) return S2Constants.M_2_PI;
            return Math.Min(ToRadians(distance) / scalar, S2Constants.M_2_PI);
        }

        // Computes the initial bearing from a to b. This is the bearing an observer
        // at point a has when facing point b. A bearing of 0 degrees is north, and it
        // increases clockwise (90 degrees is east, etc).
        // If a == b, a == -b, or a is one of the Earths' poles, the return value is
        // undefined.
        // Sourced from http://www.movable-type.co.uk/scripts/latlong.html.
        public static S1Angle GetInitialBearing(S2LatLng a, S2LatLng b)
        {
            var lat1 = a.Lat.Radians;
            var cosLat2 = Math.Cos(b.Lat.Radians);
            var lat_diff = b.Lat.Radians - a.Lat.Radians;
            var lng_diff = b.Lng.Radians - a.Lng.Radians;

            var x = Math.Sin(lat_diff) + Math.Sin(lat1) * cosLat2 * 2 * Haversine(lng_diff);
            var y = Math.Sin(lng_diff) * cosLat2;
            return S1Angle.FromRadians(Math.Atan2(y, x));
        }

        // Returns the distance between two points.  Example:
        // double miles = Miles(geostore.S2Earth.GetDistance(a, b)).value();
        //
        // Note that these methods only have single-precision accuracy, since
        // Meters is a single-precision type.  If you ned more precision, use one
        // of the methods below.
        public static double GetDistance(S2Point a, S2Point b) => ToDistance(new S1Angle(a, b));
        public static double GetDistance(S2LatLng a, S2LatLng b) => ToDistance(a.GetDistance(b));

        // Convenience functions.  These methods also return a double-precision
        // result, unlike the generic GetDistance() method.
        public static double GetDistanceKm(S2Point a, S2Point b) => RadiansToKm(a.Angle(b));
        public static double GetDistanceKm(S2LatLng a, S2LatLng b) => ToKm(a.GetDistance(b));
        public static double GetDistanceMeters(S2Point a, S2Point b) => RadiansToMeters(a.Angle(b));
        public static double GetDistanceMeters(S2LatLng a, S2LatLng b) => ToMeters(a.GetDistance(b));

        // Earth's mean radius, which is the radius of the equivalent
        // sphere with the same surface area.  According to NASA, this value is
        // 6371.01 +/- 0.02 km.  The equatorial radius is 6378.136 km, and the polar
        // radius is 6356.752 km.  They differ by one part in 298.257.
        //
        // Reference: http://ssd.jpl.nasa.gov/phys_props_earth.html, which quotes
        // Yoder, C.F. 1995. "Astrometric and Geodetic Properties of Earth and the
        // Solar System" in Global Earth Physics, A Handbook of Physical Constants,
        // AGU Reference Shelf 1, American Geophysical Union, Table 2.
        public const double RadiusKm = 0.001 * RadiusMeters;
        public const double RadiusMeters = 6371010.0;

        // altitude of the lowest known point on Earth. The lowest known
        // point on Earth is the Challenger Deep with an altitude of -10898 meters
        // above the surface of the spherical earth.
        public const double LowestAltitudeKm = 0.001 * LowestAltitudeMeters;
        public const double LowestAltitudeMeters = -10898;

        // altitude of the highest known point on Earth. The highest
        // known point on Earth is Mount Everest with an altitude of 8846 meters
        // above the surface of the spherical earth.
        public const double HighestAltitudeKm = 0.001 * HighestAltitudeMeters;
        public const double HighestAltitudeMeters = 8846;

        // http://en.wikipedia.org/wiki/Haversine_formula
        // Haversine(x) has very good numerical stability around zero.
        // Haversine(x) == (1-cos(x))/2 == sin(x/2)^2; must be implemented with the
        // second form to reap the numerical benefits.
        public static double Haversine(double radians)
        {
            var sinHalf = Math.Sin(radians / 2);
            return sinHalf * sinHalf;
        }
    }
}