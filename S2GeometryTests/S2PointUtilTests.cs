using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;
using Xunit.Abstractions;

namespace S2Geometry
{
    public class S2PointUtilTests
    {
        [Fact]
        public void Test_S2_Frames() {
            S2Point z = new S2Point(0.2, 0.5, -3.3).Normalized;
            var m = S2PointUtil.GetFrame(z);
            Assert.True(S2PointUtil.ApproxEquals(m.Col(2), z));
            Assert.True(m.Col(0).IsUnitLength);
            Assert.True(m.Col(1).IsUnitLength);
            Assert2.Near(m.Det(), 1);

            Assert.True(S2PointUtil.ApproxEquals(S2PointUtil.ToFrame(m, m.Col(0)), new S2Point(1, 0, 0)));
            Assert.True(S2PointUtil.ApproxEquals(S2PointUtil.ToFrame(m, m.Col(1)), new S2Point(0, 1, 0)));
            Assert.True(S2PointUtil.ApproxEquals(S2PointUtil.ToFrame(m, m.Col(2)), new S2Point(0, 0, 1)));

            Assert.True(S2PointUtil.ApproxEquals(S2PointUtil.FromFrame(m, new S2Point(1, 0, 0)), m.Col(0)));
            Assert.True(S2PointUtil.ApproxEquals(S2PointUtil.FromFrame(m, new S2Point(0, 1, 0)), m.Col(1)));
            Assert.True(S2PointUtil.ApproxEquals(S2PointUtil.FromFrame(m, new S2Point(0, 0, 1)), m.Col(2)));
        }

        [Fact]
        public void Test_S2_Rotate() {
            for (int iter = 0; iter < 1000; ++iter) {
                S2Point axis = S2Testing.RandomPoint();
                S2Point target = S2Testing.RandomPoint();
                // Choose a distance whose logarithm is uniformly distributed.
                double distance = Math.PI * Math.Pow(1e-15, S2Testing.Random.RandDouble());
                // Sometimes choose points near the far side of the axis.
                if (S2Testing.Random.OneIn(5)) distance = Math.PI - distance;
                S2Point p = S2EdgeDistances.InterpolateAtDistance(S1Angle.FromRadians(distance),
                                                              axis, target);
                // Choose the rotation angle.
                double angle = S2Constants.M_2_PI * Math.Pow(1e-15, S2Testing.Random.RandDouble());
                if (S2Testing.Random.OneIn(3)) angle = -angle;
                if (S2Testing.Random.OneIn(10)) angle = 0;
                TestRotate(p, axis, S1Angle.FromRadians(angle));
            }
        }

        [Fact]
        public void Test_S2_OriginTest() {
            // To minimize the number of expensive Sign() calculations,
            // S2PointUtil.Origin should not be nearly collinear with any commonly used edges.
            // Two important categories of such edges are:
            //
            //  - edges along a line of longitude (reasonably common geographically)
            //  - S2Cell edges (used extensively when computing S2Cell coverings)
            //
            // This implies that the origin:
            //
            //  - should not be too close to either pole (since all lines of longitude
            //    converge at the poles)
            //  - should not be colinear with edges of any S2Cell except for very small
            //    ones (which are used less frequently)
            //
            // The point chosen below is about 66km from the north pole towards the East
            // Siberian Sea.  The purpose of the STtoUV(2/3) calculation is to keep the
            // origin as far away as possible from the longitudinal edges of large
            // S2Cells.  (The line of longitude through the chosen point is always 1/3
            // or 2/3 of the way across any S2Cell with longitudinal edges that it
            // passes through.)

            Assert.Equal(new S2Point(-0.01, 0.01 * S2Coords.STtoUV(2.0 / 3), 1).Normalized, S2PointUtil.Origin);

            // Check that the origin is not too close to either pole.  (We don't use
            // S2Earth because we don't want to depend on that package.)
            var distance_km = Math.Acos(S2PointUtil.Origin.Z) * S2Earth.RadiusKm;
            Assert.True(distance_km >= 50.0);

            // Check that S2PointUtil.Origin is not collinear with the edges of any large
            // S2Cell.  We do this is two parts.  For S2Cells that belong to either
            // polar face, we simply need to check that S2PointUtil.Origin is not nearly
            // collinear with any edge of any cell that contains it (except for small
            // cells < 3 meters across).
            Assert.True(GetMinExpensiveLevel(S2PointUtil.Origin) >= 22);

            // For S2Cells that belong to the four non-polar faces, only longitudinal
            // edges can possibly be colinear with S2PointUtil.Origin.  We check these edges
            // by projecting S2PointUtil.Origin onto the equator, and then testing all S2Cells
            // that contain this point to make sure that none of their edges are nearly
            // colinear with S2PointUtil.Origin (except for small cells < 3 meters across).
            S2Point equator_point = new S2Point(S2PointUtil.Origin.X, S2PointUtil.Origin.Y, 0);
            Assert.True(GetMinExpensiveLevel(equator_point) >= 22);
        }

        private static void TestRotate(S2Point p, S2Point axis, S1Angle angle) {
            S2Point result = S2PointUtil.Rotate(p, axis, angle);

            // "result" should be unit length.
            Assert.True(result.IsUnitLength);

            // "result" and "p" should be the same distance from "axis".
            const double kMaxPositionError = S2Constants.DoubleError;
            Assert.True((new S1Angle(result, axis) - new S1Angle(p, axis)).Abs <= kMaxPositionError);

            // Check that the rotation angle is correct.  We allow a fixed error in the
            // *position* of the result, so we need to convert this into a rotation
            // angle.  The allowable error can be very large as "p" approaches "axis".
            double axis_distance = p.CrossProd(axis).Norm;
            double max_rotation_error;
            if (axis_distance < kMaxPositionError) {
                max_rotation_error = S2Constants.M_2_PI;
            } else {
                max_rotation_error = Math.Asin(kMaxPositionError / axis_distance);
            }
            double actual_rotation = S2Measures.TurnAngle(p, axis, result) + Math.PI;
            double rotation_error = Math.IEEERemainder(angle.Radians - actual_rotation,
                                              S2Constants.M_2_PI);
            Assert.True(rotation_error <= max_rotation_error);
        }

        // Given a point P, return the minimum level at which an edge of some S2Cell
        // parent of P is nearly collinear with S2PointUtil.Origin.  This is the minimum
        // level for which Sign() may need to resort to expensive calculations in
        // order to determine which side of an edge the origin lies on.
        private static int GetMinExpensiveLevel(S2Point p) {
            S2CellId id = new S2CellId(p);
            for (int level = 0; level <= S2Constants.kMaxCellLevel; ++level) {
                S2Cell cell = new S2Cell(id.Parent(level));
                for (int k = 0; k < 4; ++k) {
                    S2Point a = cell.GetVertex(k);
                    S2Point b = cell.GetVertex(k + 1);
                    if (S2Pred.TriageSign(a, b, S2PointUtil.Origin, a.CrossProd(b)) == 0) {
                        return level;
                    }
                }
            }
            return S2Constants.kMaxCellLevel + 1;
        }
    }
}
