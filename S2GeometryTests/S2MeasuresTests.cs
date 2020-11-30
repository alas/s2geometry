using System;
using Xunit;

namespace S2Geometry
{
    public class S2MeasuresTests
    {
        [Fact]
        public void Test_S2_AngleMethods()
        {
            S2Point pz = new S2Point(0, 0, 1);
            S2Point p000 = new S2Point(1, 0, 0);
            S2Point p045 = new S2Point(1, 1, 0).Normalized;
            S2Point p090 = new S2Point(0, 1, 0);
            S2Point p180 = new S2Point(-1, 0, 0);

            Assert2.Near(S2Measures.Angle(p000, pz, p045), S2Constants.M_PI_4);
            Assert2.Near(S2Measures.TurnAngle(p000, pz, p045), -3 * S2Constants.M_PI_4);

            Assert2.Near(S2Measures.Angle(p045, pz, p180), 3 * S2Constants.M_PI_4);
            Assert2.Near(S2Measures.TurnAngle(p045, pz, p180), -S2Constants.M_PI_4);

            Assert2.Near(S2Measures.Angle(p000, pz, p180), Math.PI);
            Assert2.Near(S2Measures.TurnAngle(p000, pz, p180), 0);

            Assert2.Near(S2Measures.Angle(pz, p000, p045), S2Constants.M_PI_2);
            Assert2.Near(S2Measures.TurnAngle(pz, p000, p045), S2Constants.M_PI_2);

            Assert2.Near(S2Measures.Angle(pz, p000, pz), 0);
            Assert2.Near(Math.Abs(S2Measures.TurnAngle(pz, p000, pz)), Math.PI);
        }

        [Fact]
        public void Test_S2_AreaMethods()
        {
            S2Point pz = new S2Point(0, 0, 1);
            S2Point p000 = new S2Point(1, 0, 0);
            S2Point p045 = new S2Point(1, 1, 0).Normalized;
            S2Point p090 = new S2Point(0, 1, 0);
            S2Point p180 = new S2Point(-1, 0, 0);

            Assert2.Near(S2Measures.Area(p000, p090, pz), S2Constants.M_PI_2);
            Assert2.Near(S2Measures.Area(p045, pz, p180), 3 * S2Constants.M_PI_4);

            // Make sure that Area() has good *relative* accuracy even for
            // very small areas.
            const double eps = 1e-10;
            S2Point pepsx = new S2Point(eps, 0, 1).Normalized;
            S2Point pepsy = new S2Point(0, eps, 1).Normalized;
            double expected1 = 0.5 * eps * eps;
            Assert2.Near(S2Measures.Area(pepsx, pepsy, pz), expected1, 1e-14 * expected1);

            // Make sure that it can handle degenerate triangles.
            S2Point pr = new S2Point(0.257, -0.5723, 0.112).Normalized;
            S2Point pq = new S2Point(-0.747, 0.401, 0.2235).Normalized;
            Assert.Equal(0, S2Measures.Area(pr, pr, pr));
            // The following test is not exact due to rounding error.
            Assert2.Near(S2Measures.Area(pr, pq, pr), 0, S2Constants.DoubleError);
            Assert.Equal(0, S2Measures.Area(p000, p045, p090));

            double max_girard = 0;
            for (int i = 0; i < 10000; ++i)
            {
                S2Point p0 = S2Testing.RandomPoint();
                S2Point d1 = S2Testing.RandomPoint();
                S2Point d2 = S2Testing.RandomPoint();
                S2Point p1 = (p0 + S2Constants.DoubleError * d1).Normalized;
                S2Point p2 = (p0 + S2Constants.DoubleError * d2).Normalized;
                // The actual displacement can be as much as 1.2e-15 due to roundoff.
                // This yields a maximum triangle area of about 0.7e-30.
                Assert.True(S2Measures.Area(p0, p1, p2) <= 0.7e-30);
                max_girard = Math.Max(max_girard, S2Measures.GirardArea(p0, p1, p2));
            }
            // This check only passes if GirardArea() uses RobustCrossProd().
            Assert.True(max_girard <= 1e-14);

            // Try a very long and skinny triangle.
            S2Point p045eps = new S2Point(1, 1, eps).Normalized;
            double expected2 = 5.8578643762690495119753e-11;  // Mathematica.
            Assert2.Near(S2Measures.Area(p000, p045eps, p090), expected2, 1e-9 * expected2);

            // Triangles with near-180 degree edges that sum to a quarter-sphere.
            const double eps2 = 1e-14;
            S2Point p000eps2 = new S2Point(1, 0.1 * eps2, eps2).Normalized;
            double quarter_area1 = S2Measures.Area(p000eps2, p000, p045) +
                                   S2Measures.Area(p000eps2, p045, p180) +
                                   S2Measures.Area(p000eps2, p180, pz) +
                                   S2Measures.Area(p000eps2, pz, p000);
            Assert2.Near(quarter_area1, Math.PI);

            // Four other triangles that sum to a quarter-sphere.
            S2Point p045eps2 = new S2Point(1, 1, eps2).Normalized;
            double quarter_area2 = S2Measures.Area(p045eps2, p000, p045) +
                                   S2Measures.Area(p045eps2, p045, p180) +
                                   S2Measures.Area(p045eps2, p180, pz) +
                                   S2Measures.Area(p045eps2, pz, p000);
            Assert2.Near(quarter_area2, Math.PI);

            // Compute the area of a hemisphere using four triangles with one near-180
            // degree edge and one near-degenerate edge.
            for (int i = 0; i < 100; ++i)
            {
                double lng = S2Constants.M_2_PI * S2Testing.Random.RandDouble();
                S2Point p0 = S2LatLng.FromRadians(1e-20, lng).Normalized.ToPoint();
                S2Point p1 = S2LatLng.FromRadians(0, lng).Normalized.ToPoint();
                double p2_lng = lng + S2Testing.Random.RandDouble();
                S2Point p2 = S2LatLng.FromRadians(0, p2_lng).Normalized.ToPoint();
                S2Point p3 = S2LatLng.FromRadians(0, lng + Math.PI).Normalized.ToPoint();
                S2Point p4 = S2LatLng.FromRadians(0, lng + 5.0).Normalized.ToPoint();
                double area = (S2Measures.Area(p0, p1, p2) + S2Measures.Area(p0, p2, p3) +
                               S2Measures.Area(p0, p3, p4) + S2Measures.Area(p0, p4, p1));
                Assert2.Near(area, S2Constants.M_2_PI, 2e-15);
            }

            // This tests a case where the triangle has zero area, but S2Measures.Area()
            // computes (dmin > 0) due to rounding errors.
            Assert.Equal(0.0, S2Measures.Area(S2LatLng.FromDegrees(-45, -170).ToPoint(),
                                    S2LatLng.FromDegrees(45, -170).ToPoint(),
                                    S2LatLng.FromDegrees(0, -170).ToPoint()));
        }
    }
}
