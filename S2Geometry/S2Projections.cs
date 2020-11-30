// Defines a few simple map projections.  (Clients that need more complex
// projections should use a third-party library such as GeographicLib to
// implement their own projection subtypes.)

using System;

namespace S2Geometry
{
    // For the purposes of the S2 library, a projection is a function that maps
    // between S2Points and R2Vectors.  It can also define the coordinate wrapping
    // behavior along each axis.
    public abstract class Projection
    {
        // Converts a point on the sphere to a projected 2D point.
        public abstract R2Point Project(S2Point p);

        // Converts a projected 2D point to a point on the sphere.
        //
        // If wrapping is defined for a given axis (see below), then this method
        // should accept any real number for the corresponding coordinate.
        public abstract S2Point Unproject(R2Point p);

        // Convenience function equivalent to Project(ll.ToPoint()), but the
        // implementation may be more efficient.
        public abstract R2Point FromLatLng(S2LatLng ll);

        // Convenience function equivalent to S2LatLng(Unproject(p)), but the
        // implementation may be more efficient.
        public abstract S2LatLng ToLatLng(R2Point p);

        // Returns the point obtained by interpolating the given fraction of the
        // distance along the line from A to B.  Almost all projections should
        // use the default implementation of this method, which simply interpolates
        // linearly in R2 space.  Fractions < 0 or > 1 result in extrapolation
        // instead.
        //
        // The only reason to override this method is if you want edges to be
        // defined as something other than straight lines in the 2D projected
        // coordinate system.  For example, using a third-party library such as
        // GeographicLib you could define edges as geodesics over an ellipsoid model
        // of the Earth.  (Note that very few data sets define edges this way.)
        //
        // Also note that there is no reason to define a projection where edges are
        // geodesics over the sphere, because this is the native S2 interpretation.
        public static R2Point Interpolate(double f, R2Point a, R2Point b)
        {
            // Default implementation, suitable for any projection where edges are defined
            // as straight lines in the 2D projected space.
            return (1.0 - f) * a + f * b;
        }

        // Defines the coordinate wrapping distance along each axis.  If this value
        // is non-zero for a given axis, the coordinates are assumed to "wrap" with
        // the given period.  For example, if wrap_distance.Y == 360 then (x, y)
        // and (x, y + 360) should map to the same S2Point.
        //
        // This information is used to ensure that edges takes the shortest path
        // between two given points.  For example, if coordinates represent
        // (latitude, longitude) pairs in degrees and wrap_distance().Y == 360,
        // then the edge (5:179, 5:-179) would be interpreted as spanning 2 degrees
        // of longitude rather than 358 degrees.
        //
        // If a given axis does not wrap, its wrap_distance should be set to zero.
        public abstract R2Point WrapDistance();

        // Helper function that wraps the coordinates of B if necessary in order to
        // obtain the shortest edge AB.  For example, suppose that A = [170, 20],
        // B = [-170, 20], and the projection wraps so that [x, y] == [x + 360, y].
        // Then this function would return [190, 20] for point B (reducing the edge
        // length in the "x" direction from 340 to 20).
        public R2Point WrapDestination(R2Point a, R2Point b)
        {
            R2Point wrap = WrapDistance();
            double x = b.X, y = b.Y;
            // The code below ensures that "b" is unmodified unless wrapping is required.
            if (wrap.X > 0 && Math.Abs(x - a.X) > 0.5 * wrap.X)
            {
                x = a.X + Math.IEEERemainder(x - a.X, wrap.X);
            }
            if (wrap.Y > 0 && Math.Abs(y - a.Y) > 0.5 * wrap.Y)
            {
                y = a.Y + Math.IEEERemainder(y - a.Y, wrap.Y);
            }
            return new R2Point(x, y);
        }
    }

    // PlateCarreeProjection defines the "plate carree" (square plate) projection,
    // which converts points on the sphere to (longitude, latitude) pairs.
    // Coordinates can be scaled so that they represent radians, degrees, etc, but
    // the projection is always centered around (latitude=0, longitude=0).
    //
    // Note that (x, y) coordinates are backwards compared to the usual (latitude,
    // longitude) ordering, in order to match the usual convention for graphs in
    // which "x" is horizontal and "y" is vertical.
    public sealed class PlateCarreeProjection : Projection
    {
        // Constructs the plate carree projection where the x coordinates
        // (longitude) span [-x_scale, x_scale] and the y coordinates (latitude)
        // span [-x_scale/2, x_scale/2].  For example if x_scale==180 then the x
        // range is [-180, 180] and the y range is [-90, 90].
        //
        // By default coordinates are expressed in radians, i.e. the x range is
        // [-Pi, Pi] and the y range is [-Pi/2, Pi/2].
        public PlateCarreeProjection(double x_scale = Math.PI)
        {
            x_wrap_ = 2 * x_scale;
            to_radians_ = Math.PI / x_scale;
            from_radians_ = x_scale / Math.PI;
        }

        public override R2Point Project(S2Point p)
        {
            return FromLatLng(new S2LatLng(p));
        }

        public override S2Point Unproject(R2Point p)
        {
            return ToLatLng(p).ToPoint();
        }

        public override R2Point FromLatLng(S2LatLng ll)
        {
            return new R2Point(from_radians_ * ll.LngRadians,
                                from_radians_ * ll.LatRadians);
        }

        public override S2LatLng ToLatLng(R2Point p)
        {
            var rem = Math.IEEERemainder(p.X, x_wrap_);
            return S2LatLng.FromRadians(to_radians_ * p.Y, to_radians_ * rem);
        }

        public override R2Point WrapDistance()
        {
            return new R2Point(x_wrap_, 0);
        }

        private readonly double x_wrap_;
        private readonly double to_radians_;    // Multiplier to convert coordinates to radians.
        private readonly double from_radians_;  // Multiplier to convert coordinates from radians.
    };

    // MercatorProjection defines the spherical Mercator projection.  Google Maps
    // uses this projection together with WGS84 coordinates, in which case it is
    // known as the "Web Mercator" projection (see Wikipedia).  This class makes
    // no assumptions regarding the coordinate system of its input points, but
    // simply applies the spherical Mercator projection to them.
    //
    // The Mercator projection is finite in width (x) but infinite in height (y).
    // "x" corresponds to longitude, and spans a finite range such as [-180, 180]
    // (with coordinate wrapping), while "y" is a function of latitude and spans
    // an infinite range.  (As "y" coordinates get larger, points get closer to
    // the north pole but never quite reach it.)  The north and south poles have
    // infinite "y" values.  (Note that this will cause problems if you tessellate
    // a Mercator edge where one endpoint is a pole.  If you need to do this, clip
    // the edge first so that the "y" coordinate is no more than about 5 * max_x.)
    public sealed class MercatorProjection : Projection 
    {
        // Constructs a Mercator projection where "x" corresponds to longitude in
        // the range [-max_x, max_x] , and "y" corresponds to latitude and can be
        // any real number.  The horizontal and vertical scales are equal locally.
        public MercatorProjection(double max_x)
        {
            x_wrap_ = 2 * max_x;
            to_radians_ = Math.PI / max_x;
            from_radians_ = max_x / Math.PI;
        }

        public override R2Point Project(S2Point p)
        {
            return FromLatLng(new S2LatLng(p));
        }

        public override S2Point Unproject(R2Point p)
        {
            return ToLatLng(p).ToPoint();
        }

        public override R2Point FromLatLng(S2LatLng ll)
        {
            // This formula is more accurate near zero than the log(tan()) version.
            // Note that latitudes of +/- 90 degrees yield "y" values of +/- infinity.
            double sin_phi = Math.Sin(ll.LatRadians);
            double y = 0.5 * Math.Log((1 + sin_phi) / (1 - sin_phi));
            return new R2Point(from_radians_ * ll.LngRadians, from_radians_ * y);
        }

        public override S2LatLng ToLatLng(R2Point p)
        {
            // This formula is more accurate near zero than the atan(exp()) version.
            double x = to_radians_ * Math.IEEERemainder(p.X, x_wrap_);
            double k = Math.Exp(2 * to_radians_ * p.Y);
            double y = double.IsInfinity(k) ? S2Constants.M_PI_2 : Math.Asin((k - 1) / (k + 1));
            return S2LatLng.FromRadians(y, x);
        }

        public override R2Point WrapDistance()
        {
            return new R2Point(x_wrap_, 0);
        }

        private readonly double x_wrap_;
        private readonly double to_radians_;    // Multiplier to convert coordinates to radians.
        private readonly double from_radians_;  // Multiplier to convert coordinates from radians.
    }
}

