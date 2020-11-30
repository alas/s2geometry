using System;
using Xunit;

namespace S2Geometry
{
    public class S2LatLngRectBounderTests
    {
        private readonly S2LatLng kRectError = S2LatLngRectBounder.MaxErrorForTests();

        [Fact]
        public void Test_RectBounder_MaxLatitudeSimple() {
            // Check cases where the min/max latitude is attained at a vertex.
            var kCubeLat = Math.Asin(1 / Math.Sqrt(3));  // 35.26 degrees
            Assert.True(GetEdgeBound(1, 1, 1, 1, -1, -1).ApproxEquals(  // NOLINT
                new S2LatLngRect(new R1Interval(-kCubeLat, kCubeLat),
                             new S1Interval(-S2Constants.M_PI_4, S2Constants.M_PI_4)), kRectError));
            Assert.True(GetEdgeBound(1, -1, 1, 1, 1, -1).ApproxEquals(  // NOLINT
                new S2LatLngRect(new R1Interval(-kCubeLat, kCubeLat),
                             new S1Interval(-S2Constants.M_PI_4, S2Constants.M_PI_4)), kRectError));

            // Check cases where the min/max latitude occurs in the edge interior.
            // These tests expect the result to be pretty close to the middle of the
            // allowable error range (i.e., by adding 0.5 * kRectError).

            // Max latitude, CW edge
            Assert2.Near(S2Constants.M_PI_4 + 0.5 * kRectError.LatRadians,
                             GetEdgeBound(1, 1, 1, 1, -1, 1).Lat.Hi);
            // Max latitude, CCW edge
            Assert2.Near(S2Constants.M_PI_4 + 0.5 * kRectError.LatRadians,
                             GetEdgeBound(1, -1, 1, 1, 1, 1).Lat.Hi);  // NOLINT
                                                                           // Min latitude, CW edge
            Assert2.Near(-S2Constants.M_PI_4 - 0.5 * kRectError.LatRadians,
                             GetEdgeBound(1, -1, -1, -1, -1, -1).Lat.Lo);  // NOLINT
                                                                               // Min latitude, CCW edge
            Assert2.Near(-S2Constants.M_PI_4 - 0.5 * kRectError.LatRadians,
                             GetEdgeBound(-1, 1, -1, -1, -1, -1).Lat.Lo);  // NOLINT

            // Check cases where the edge passes through one of the poles.
            Assert.Equal(S2Constants.M_PI_2, GetEdgeBound(.3, .4, 1, -.3, -.4, 1).Lat.Hi);  // NOLINT
            Assert.Equal(-S2Constants.M_PI_2, GetEdgeBound(.3, .4, -1, -.3, -.4, -1).Lat.Lo);  // NOLINT
        }

        [Fact]
        public void Test_RectBounder_MaxLatitudeRandom() {
            // Check that the maximum latitude of edges is computed accurately to within
            // 3 * S2Constants.DoubleEpsilon (the expected maximum error).  We concentrate on maximum
            // latitudes near the equator and north pole since these are the extremes.

            const int kIters = 100;
            for (int iter = 0; iter < kIters; ++iter) {
                // Construct a right-handed coordinate frame (U,V,W) such that U points
                // slightly above the equator, V points at the equator, and W is slightly
                // offset from the north pole.
                S2Point u = S2Testing.RandomPoint();
                u = new S2Point(u.X, u.Y, S2Constants.DoubleEpsilon * 1e-6 * Math.Pow(1e12, S2Testing.Random.RandDouble()));  // log is uniform
                u = u.Normalized;
                S2Point v = S2PointUtil.RobustCrossProd(new S2Point(0, 0, 1), u).Normalized;
                S2Point w = S2PointUtil.RobustCrossProd(u, v).Normalized;

                // Construct a line segment AB that passes through U, and check that the
                // maximum latitude of this segment matches the latitude of U.
                S2Point a = (u - S2Testing.Random.RandDouble() * v).Normalized;
                S2Point b = (u + S2Testing.Random.RandDouble() * v).Normalized;
                S2LatLngRect ab_bound = GetEdgeBound(a, b);
                Assert2.Near(S2LatLng.Latitude(u).Radians,
                            ab_bound.Lat.Hi, kRectError.LatRadians);

                // Construct a line segment CD that passes through W, and check that the
                // maximum latitude of this segment matches the latitude of W.
                S2Point c = (w - S2Testing.Random.RandDouble() * v).Normalized;
                S2Point d = (w + S2Testing.Random.RandDouble() * v).Normalized;
                S2LatLngRect cd_bound = GetEdgeBound(c, d);
                Assert2.Near(S2LatLng.Latitude(w).Radians,
                            cd_bound.Lat.Hi, kRectError.LatRadians);
            }
        }

        [Fact]
        public void Test_RectBounder_NearlyIdenticalOrAntipodalPoints()
        {
            Assert.True(false); //TODO

            // Test pairs of points that are either:
            //  - identical
            //  - nearly or exactly proportional, e.g. (1,0,0) vs. (1+2e-16, 0, 0)
            //  - very close to each other
            // Furthermore we want to test cases where the two points are:
            //  - on a nearly-polar great circle
            //  - on a nearly-equatorial great circle
            //  - near the poles, but on any great circle
            //  - near the equator, but on any great circle
            //  - positioned arbitrarily
            // Also test the corresponding situations for antipodal points, i.e. by
            // negating one of the points so that they are almost 180 degrees apart.

            int kIters = 10000;
            for (int iter = 0; iter < kIters; ++iter) {
                S2Point a, b;
                switch (S2Testing.Random.Uniform(5)) {
                    case 0:
                        // Two nearby points on a nearly-polar great circle.
                        a = S2Testing.RandomPoint();
                        b = PerturbATowardsB(a, PointNearPole());
                        break;
                    case 1:
                        // Two nearby points on a nearly-equatorial great circle.
                        a = PointNearEquator();
                        b = PerturbATowardsB(a, PointNearEquator());
                        break;
                    case 2:
                        // Two nearby points near a pole, but on any great circle.
                        a = PointNearPole();
                        b = PerturbATowardsB(a, S2Testing.RandomPoint());
                        break;
                    case 3:
                        // Two nearby points near the equator, but on any great circle.
                        a = PointNearEquator();
                        b = PerturbATowardsB(a, S2Testing.RandomPoint());
                        break;
                    case 4:
                        // Two nearby points anywhere on the sphere.
                        a = S2Testing.RandomPoint();
                        b = PerturbATowardsB(a, S2Testing.RandomPoint());
                        break;
                    default: throw new ArgumentOutOfRangeException("Random");
                }

                //TODO, when a[0] == -90d this test fails, "case 2:"

                // The two points are chosen to be so close to each other that the min/max
                // latitudes are nearly always achieved at the edge endpoints.  The only
                // thing we need to watch out for is that the latitude error bound is
                // slightly larger if the min/max latitude occurs in the edge interior.
                S2LatLngRect expected_bound = S2LatLngRect.FromPointPair(new S2LatLng(a), new S2LatLng(b));
                S2LatLngRect bound = GetEdgeBound(a, b);
                Assert.True(bound.Contains(expected_bound));
                Assert.True(expected_bound.Expanded(kRectError).PolarClosure().Contains(bound));

                // If the two points are close enough and one point is negated (antipodal
                // points), the bound should be the entire sphere.
                if ((a - b).CrossProd(a + b).Norm <= 6.110 * S2Constants.DoubleEpsilon) {
                    Assert.Equal(S2LatLngRect.Full, GetEdgeBound(a, -b));
                }
            }
        }

        [Fact]
        public void Test_RectBounder_ExpandForSubregions() {
            // First we check the various situations where the bound contains
            // nearly-antipodal points.  The tests are organized into pairs where the
            // two bounds are similar except that the first bound meets the
            // nearly-antipodal criteria while the second does not.

            // Cases where the bound does not straddle the equator (but almost does),
            // and spans nearly 180 degrees in longitude.
            Assert.True(GetSubregionBound(3e-16, 0, 1e-14, Math.PI).IsFull);
            Assert.False(GetSubregionBound(9e-16, 0, 1e-14, Math.PI).IsFull);
            Assert.True(GetSubregionBound(1e-16, 7e-16, 1e-14, Math.PI).IsFull);
            Assert.False(GetSubregionBound(3e-16, 14e-16, 1e-14, Math.PI).IsFull);
            Assert.True(GetSubregionBound(1e-100, 14e-16, 1e-14, Math.PI).IsFull);
            Assert.False(GetSubregionBound(1e-100, 22e-16, 1e-14, Math.PI).IsFull);

            // Cases where the bound spans at most 90 degrees in longitude, and almost
            // 180 degrees in latitude.  Note that S2Constants.DoubleEpsilon is about 2.22e-16, which
            // implies that the double-precision value just below Pi/2 can be written as
            // (S2Constants.M_PI_2 - 2e-16).
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2, -1e-15, S2Constants.M_PI_2 - 7e-16, 0).
                      IsFull);
            Assert.False(GetSubregionBound(-S2Constants.M_PI_2, -1e-15, S2Constants.M_PI_2 - 30e-16, 0).
                        IsFull);
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2 + 4e-16, 0, S2Constants.M_PI_2 - 2e-16, 1e-7).
                        IsFull);
            Assert.False(GetSubregionBound(-S2Constants.M_PI_2 + 30e-16, 0, S2Constants.M_PI_2, 1e-7).
                        IsFull);
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2 + 4e-16, 0, S2Constants.M_PI_2 - 4e-16, S2Constants.M_PI_2).
                        IsFull);
            Assert.False(GetSubregionBound(-S2Constants.M_PI_2, 0, S2Constants.M_PI_2 - 30e-16, S2Constants.M_PI_2).
                        IsFull);

            // Cases where the bound straddles the equator and spans more than 90
            // degrees in longitude.  These are the cases where the critical distance is
            // between a corner of the bound and the opposite longitudinal edge.  Unlike
            // the cases above, here the bound may contain nearly-antipodal points (to
            // within 3.055 * S2Constants.DoubleEpsilon) even though the latitude and longitude ranges
            // are both significantly less than (Pi - 3.055 * S2Constants.DoubleEpsilon).
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2, 0, S2Constants.M_PI_2 - 1e-8, Math.PI - 1e-7).
                      IsFull);
            Assert.False(GetSubregionBound(-S2Constants.M_PI_2, 0, S2Constants.M_PI_2 - 1e-7, Math.PI - 1e-7).
                         IsFull);
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2 + 1e-12, -Math.PI + 1e-4, S2Constants.M_PI_2, 0).
                        IsFull);
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2 + 1e-11, -Math.PI + 1e-4, S2Constants.M_PI_2, 0).
                         IsFull);

            // Now we test cases where the bound does not contain nearly-antipodal
            // points, but it does contain points that are approximately 180 degrees
            // apart in latitude.
            Assert.True(GetSubregionBound(1.5, -S2Constants.M_PI_2, 1.5, S2Constants.M_PI_2 - 2e-16).
                        ApproxEquals(new S2LatLngRect(new R1Interval(1.5, 1.5),
                                                  S1Interval.Full), kRectError));
            Assert.True(GetSubregionBound(1.5, -S2Constants.M_PI_2, 1.5, S2Constants.M_PI_2 - 7e-16).
                        ApproxEquals(new S2LatLngRect(new R1Interval(1.5, 1.5),
                                                  new S1Interval(-S2Constants.M_PI_2, S2Constants.M_PI_2 - 7e-16)),
                                     kRectError));

            // Test the full and empty bounds.
            Assert.True(S2LatLngRectBounder.ExpandForSubregions(
                S2LatLngRect.Full).IsFull);
            Assert.True(S2LatLngRectBounder.ExpandForSubregions(
                S2LatLngRect.Empty).IsEmpty);

            // Check for cases where the bound is expanded to include one of the poles.
            Assert.True(GetSubregionBound(-S2Constants.M_PI_2 + 1e-15, 0, -S2Constants.M_PI_2 + 1e-15, 0).
                        ApproxEquals(new S2LatLngRect(new R1Interval(-S2Constants.M_PI_2, -S2Constants.M_PI_2 + 1e-15),
                                                  S1Interval.Full), kRectError));
            Assert.True(GetSubregionBound(S2Constants.M_PI_2 - 1e-15, 0, S2Constants.M_PI_2 - 1e-15, 0).
                        ApproxEquals(new S2LatLngRect(new R1Interval(S2Constants.M_PI_2 - 1e-15, S2Constants.M_PI_2),
                                                  S1Interval.Full), kRectError));
        }

        private static S2Point PerturbATowardsB(S2Point a, S2Point b) {
            var choice = S2Testing.Random.RandDouble();
            if (choice < 0.1) {
                return a;
            }
            if (choice < 0.3) {
                // Return a point that is exactly proportional to A and that still
                // satisfies S2.IsUnitLength().
                for (; ; ) {
                    var b2 = (2 - a.Norm + 5 * (S2Testing.Random.RandDouble() - 0.5) * S2Constants.DoubleEpsilon) * a;
                    if (b2 != a && b2.IsUnitLength)
                        return b2;
                }
            }
            if (choice < 0.5) {
                // Return a point such that the distance squared to A will underflow.
                return S2EdgeDistances.InterpolateAtDistance(S1Angle.FromRadians(1e-300), a, b);
            }
            // Otherwise return a point whose distance from A is near S2Constants.DoubleEpsilon such
            // that the log of the pdf is uniformly distributed.
            double distance = S2Constants.DoubleEpsilon * 1e-5 * Math.Pow(1e6, S2Testing.Random.RandDouble());
            return S2EdgeDistances.InterpolateAtDistance(S1Angle.FromRadians(distance), a, b);
        }

        private static S2Point RandomPole() {
            return new S2Point(0, 0, S2Testing.Random.OneIn(2) ? 1 : -1);
        }

        private static S2Point PointNearPole() {
            return PerturbATowardsB(RandomPole(), S2Testing.RandomPoint());
        }

        private static S2Point PointNearEquator() {
            return PerturbATowardsB(
                new S2Point(S2Testing.Random.RandDouble(),
                S2Testing.Random.RandDouble(), 0).Normalized,
                RandomPole());
        }

        private static S2LatLngRect GetSubregionBound(double x_lat, double x_lng, double y_lat, double y_lng) {
            S2LatLngRect in_ = S2LatLngRect.FromPointPair(
                S2LatLng.FromRadians(x_lat, x_lng),
                S2LatLng.FromRadians(y_lat, y_lng));
            S2LatLngRect out_ = S2LatLngRectBounder.ExpandForSubregions(in_);

            // Test that the bound is actually expanded.
            Assert.True(out_.Contains(in_));
            if (in_.Lat == S2LatLngRect.FullLat) {
                Assert.False(in_.Lat.Contains(out_.Lat));
            }
            return out_;
        }
        
        private static S2LatLngRect GetEdgeBound(S2Point a, S2Point b)
        {
            var bounder = new S2LatLngRectBounder();
            bounder.AddPoint(a);
            bounder.AddPoint(b);
            return bounder.GetBound();
        }

        private static S2LatLngRect GetEdgeBound(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            return GetEdgeBound(new S2Point(x1, y1, z1).Normalized,
                                new S2Point(x2, y2, z2).Normalized);
        }
    }
}
