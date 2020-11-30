using System;
using System.Collections.Generic;
using Xunit;

namespace S2Geometry
{
    public class S2CapTests
    {
        // About 9 times the double-precision roundoff relative error.
        private const double kEps = 1e-15;
        private const double kTinyRad = 1e-10;
        private const double kDegreeEps = 1e-13;
        // The distance from the center of a face to one of its vertices.
        private static readonly double kFaceRadius = Math.Atan(Math.Sqrt(2));

        [Fact]
        public void Test_S2Cap_Basic()
        {
            // Test basic properties of empty and full caps.
            S2Cap empty = S2Cap.Empty;
            S2Cap full = S2Cap.Full;
            Assert.True(empty.IsValid);
            Assert.True(empty.IsEmpty);
            Assert.True(empty.Complement().IsFull);
            Assert.True(full.IsValid);
            Assert.True(full.IsFull);
            Assert.True(full.Complement().IsEmpty);
            Assert.Equal(2, full.Height);
            Assert2.Near(180.0, full.Radius.Degrees);

            // Test the S1Angle constructor using out-of-range arguments.
            Assert.True(new S2Cap(new S2Point(1, 0, 0), S1Angle.FromRadians(-20)).IsEmpty);
            Assert.True(new S2Cap(new S2Point(1, 0, 0), S1Angle.FromRadians(5)).IsFull);
            Assert.True(new S2Cap(new S2Point(1, 0, 0), S1Angle.Infinity).IsFull);

            // Check that the default S2Cap is identical to Empty().
            var default_empty = S2Cap.Empty;
            Assert.True(default_empty.IsValid);
            Assert.True(default_empty.IsEmpty);
            Assert.Equal(empty.Center, default_empty.Center);
            Assert.Equal(empty.Height, default_empty.Height);

            // Containment and intersection of empty and full caps.
            Assert.True(empty.Contains(empty));
            Assert.True(full.Contains(empty));
            Assert.True(full.Contains(full));
            Assert.False(empty.InteriorIntersects(empty));
            Assert.True(full.InteriorIntersects(full));
            Assert.False(full.InteriorIntersects(empty));

            // Singleton cap containing the x-axis.
            S2Cap xaxis = S2Cap.FromPoint(new S2Point(1, 0, 0));
            Assert.True(xaxis.Contains(new S2Point(1, 0, 0)));
            Assert.False(xaxis.Contains(new S2Point(1, 1e-20, 0)));
            Assert.Equal(0, xaxis.Radius.Radians);

            // Singleton cap containing the y-axis.
            S2Cap yaxis = S2Cap.FromPoint(new S2Point(0, 1, 0));
            Assert.False(yaxis.Contains(xaxis.Center));
            Assert.Equal(0, xaxis.Height);

            // Check that the complement of a singleton cap is the full cap.
            S2Cap xcomp = xaxis.Complement();
            Assert.True(xcomp.IsValid);
            Assert.True(xcomp.IsFull);
            Assert.True(xcomp.Contains(xaxis.Center));

            // Check that the complement of the complement is *not* the original.
            Assert.True(xcomp.Complement().IsValid);
            Assert.True(xcomp.Complement().IsEmpty);
            Assert.False(xcomp.Complement().Contains(xaxis.Center));

            // Check that very small caps can be represented accurately.
            // Here "kTinyRad" is small enough that unit vectors perturbed by this
            // amount along a tangent do not need to be renormalized.
            S2Cap tiny = new S2Cap(new S2Point(1, 2, 3).Normalized, S1Angle.FromRadians(kTinyRad));
            var tangent = tiny.Center.CrossProd(new S2Point(3, 2, 1)).Normalized;
            Assert.True(tiny.Contains(tiny.Center + 0.99 * kTinyRad * tangent));
            Assert.False(tiny.Contains(tiny.Center + 1.01 * kTinyRad * tangent));

            // Basic tests on a hemispherical cap.
            S2Cap hemi = S2Cap.FromCenterHeight(new S2Point(1, 0, 1).Normalized, 1);
            Assert.Equal(-hemi.Center, hemi.Complement().Center);
            Assert.Equal(1, hemi.Complement().Height);
            Assert.True(hemi.Contains(new S2Point(1, 0, 0)));
            Assert.False(hemi.Complement().Contains(new S2Point(1, 0, 0)));
            Assert.True(hemi.Contains(new S2Point(1, 0, -(1 - kEps)).Normalized));
            Assert.False(hemi.InteriorContains(new S2Point(1, 0, -(1 + kEps)).Normalized));

            // A concave cap.  Note that the error bounds for point containment tests
            // increase with the cap angle, so we need to use a larger error bound
            // here.  (It would be painful to do this everywhere, but this at least
            // gives an example of how to compute the maximum error.)
            S2Point center = GetLatLngPoint(80, 10);
            S1ChordAngle radius = new S1ChordAngle(S1Angle.FromDegrees(150));
            double max_error =
                radius.GetS2PointConstructorMaxError() +
                radius.S1AngleConstructorMaxError +
                3 * S2Constants.DoubleEpsilon;  // GetLatLngPoint() error
            S2Cap concave = new S2Cap(center, radius);
            S2Cap concave_min = new S2Cap(center, radius.PlusError(-max_error));
            S2Cap concave_max = new S2Cap(center, radius.PlusError(max_error));
            Assert.True(concave_max.Contains(GetLatLngPoint(-70, 10)));
            Assert.False(concave_min.Contains(GetLatLngPoint(-70, 10)));
            Assert.True(concave_max.Contains(GetLatLngPoint(-50, -170)));
            Assert.False(concave_min.Contains(GetLatLngPoint(-50, -170)));

            // Cap containment tests.
            Assert.False(empty.Contains(xaxis));
            Assert.False(empty.InteriorIntersects(xaxis));
            Assert.True(full.Contains(xaxis));
            Assert.True(full.InteriorIntersects(xaxis));
            Assert.False(xaxis.Contains(full));
            Assert.False(xaxis.InteriorIntersects(full));
            Assert.True(xaxis.Contains(xaxis));
            Assert.False(xaxis.InteriorIntersects(xaxis));
            Assert.True(xaxis.Contains(empty));
            Assert.False(xaxis.InteriorIntersects(empty));
            Assert.True(hemi.Contains(tiny));
            Assert.True(hemi.Contains(new S2Cap(new S2Point(1, 0, 0), S1Angle.FromRadians(S2Constants.M_PI_4 - kEps))));
            Assert.False(hemi.Contains(new S2Cap(new S2Point(1, 0, 0), S1Angle.FromRadians(S2Constants.M_PI_4 + kEps))));
            Assert.True(concave.Contains(hemi));
            Assert.True(concave.InteriorIntersects(hemi.Complement()));
            Assert.False(concave.Contains(S2Cap.FromCenterHeight(-concave.Center, 0.1)));
        }

        [Fact]
        public void Test_S2Cap_AddEmptyCapToNonEmptyCap()
        {
            S2Cap non_empty_cap = new S2Cap(new S2Point(1, 0, 0), S1Angle.FromDegrees(10));
            double initial_area = non_empty_cap.Area;
            non_empty_cap = non_empty_cap.AddCap(S2Cap.Empty);
            Assert.Equal(initial_area, non_empty_cap.Area);
        }

        [Fact]
        public void Test_S2Cap_AddNonEmptyCapToEmptyCap()
        {
            S2Cap empty = S2Cap.Empty;
            S2Cap non_empty_cap = new S2Cap(new S2Point(1, 0, 0), S1Angle.FromDegrees(10));
            empty = empty.AddCap(non_empty_cap);
            Assert.Equal(non_empty_cap.Area, empty.Area);
        }

        [Fact]
        public void Test_S2Cap_GetRectBound()
        {
            // Empty and full caps.
            Assert.True(S2Cap.Empty.GetRectBound().IsEmpty);
            Assert.True(S2Cap.Full.GetRectBound().IsFull);

            // Maximum allowable error for latitudes and longitudes measured in
            // degrees.  (Assert2.Near isn't sufficient.)

            // Cap that includes the south pole.
            S2LatLngRect rect = new S2Cap(GetLatLngPoint(-45, 57),
                S1Angle.FromDegrees(50)).GetRectBound();
            Assert2.Near(rect.LatLo.Degrees, -90, kDegreeEps);
            Assert2.Near(rect.LatHi.Degrees, 5, kDegreeEps);
            Assert.True(rect.Lng.IsFull);

            // Cap that is tangent to the north pole.
            rect = new S2Cap(new S2Point(1, 0, 1).Normalized,
                S1Angle.FromRadians(S2Constants.M_PI_4 + 1e-16)).GetRectBound();
            Assert2.Near(rect.Lat.Lo, 0, kEps);
            Assert2.Near(rect.Lat.Hi, S2Constants.M_PI_2, kEps);
            Assert.True(rect.Lng.IsFull);

            var p = new S2Point(1, 0, 1).Normalized;
            var rb1 = S1Angle.FromDegrees(45 + 5e-15);
            rect = new S2Cap(p, rb1).GetRectBound();
            Assert2.Near(rect.LatLo.Degrees, 0, kDegreeEps);
            Assert2.Near(rect.LatHi.Degrees, 90, kDegreeEps);
            Assert.True(rect.Lng.IsFull);

            // The eastern hemisphere.
            var rb2 = S1Angle.FromRadians(S2Constants.M_PI_2 + 2e-16);
            rect = new S2Cap(new S2Point(0, 1, 0), rb2).GetRectBound();
            Assert2.Near(rect.LatLo.Degrees, -90, kDegreeEps);
            Assert2.Near(rect.LatHi.Degrees, 90, kDegreeEps);
            Assert.True(rect.Lng.IsFull);

            // A cap centered on the equator.
            rect = new S2Cap(GetLatLngPoint(0, 50), S1Angle.FromDegrees(20)).GetRectBound();
            Assert2.Near(rect.LatLo.Degrees, -20, kDegreeEps);
            Assert2.Near(rect.LatHi.Degrees, 20, kDegreeEps);
            Assert2.Near(rect.LngLo.Degrees, 30, kDegreeEps);
            Assert2.Near(rect.LngHi.Degrees, 70, kDegreeEps);

            // A cap centered on the north pole.
            rect = new S2Cap(GetLatLngPoint(90, 123), S1Angle.FromDegrees(10)).GetRectBound();
            Assert2.Near(rect.LatLo.Degrees, 80, kDegreeEps);
            Assert2.Near(rect.LatHi.Degrees, 90, kDegreeEps);
            Assert.True(rect.Lng.IsFull);
        }

        [Fact]
        public void Test_S2Cap_S2CellMethods()
        {
            // For each cube face, we construct some cells on
            // that face and some caps whose positions are relative to that face,
            // and then check for the expected intersection/containment results.

            for (var face = 0; face < 6; ++face) {
                // The cell consisting of the entire face.
                S2Cell root_cell = S2Cell.FromFace(face);

                // A leaf cell at the midpoint of the v=1 edge.
                S2Cell edge_cell = new S2Cell(S2Coords.FaceUVtoXYZ(face, 0, 1 - kEps));

                // A leaf cell at the u=1, v=1 corner.
                S2Cell corner_cell = new S2Cell(S2Coords.FaceUVtoXYZ(face, 1 - kEps, 1 - kEps));

                // Quick check for full and empty caps.
                Assert.True(S2Cap.Full.Contains(root_cell));
                Assert.False(S2Cap.Empty.MayIntersect(root_cell));

                // Check intersections with the bounding caps of the leaf cells that are
                // adjacent to 'corner_cell' along the Hilbert curve.  Because this corner
                // is at (u=1,v=1), the curve stays locally within the same cube face.
                S2CellId first = corner_cell.Id.Advance(-3);
                S2CellId last = corner_cell.Id.Advance(4);
                for (S2CellId id = first; id < last; id = id.Next) {
                    S2Cell cell = new S2Cell(id);
                    Assert.Equal(id == corner_cell.Id,
                              cell.GetCapBound().Contains(corner_cell));
                    Assert.Equal(id.Parent().Contains(corner_cell.Id),
                              cell.GetCapBound().MayIntersect(corner_cell));
                }

                var anti_face = (face + 3) % 6;  // Opposite face.
                for (var cap_face = 0; cap_face < 6; ++cap_face) {
                    // A cap that barely contains all of 'cap_face'.
                    S2Point center = S2Coords.GetNorm(cap_face);
                    S2Cap covering = new S2Cap(center, S1Angle.FromRadians(kFaceRadius + kEps));
                    Assert.Equal(cap_face == face, covering.Contains(root_cell));
                    Assert.Equal(cap_face != anti_face, covering.MayIntersect(root_cell));
                    Assert.Equal(center.DotProd(edge_cell.GetCenter()) > 0.1,
                              covering.Contains(edge_cell));
                    Assert.Equal(covering.MayIntersect(edge_cell), covering.Contains(edge_cell));
                    Assert.Equal(cap_face == face, covering.Contains(corner_cell));
                    Assert.Equal(center.DotProd(corner_cell.GetCenter()) > 0,
                              covering.MayIntersect(corner_cell));

                    // A cap that barely intersects the edges of 'cap_face'.
                    S2Cap bulging = new S2Cap(center, S1Angle.FromRadians(S2Constants.M_PI_4 + kEps));
                    Assert.False(bulging.Contains(root_cell));
                    Assert.Equal(cap_face != anti_face, bulging.MayIntersect(root_cell));
                    Assert.Equal(cap_face == face, bulging.Contains(edge_cell));
                    Assert.Equal(center.DotProd(edge_cell.GetCenter()) > 0.1,
                              bulging.MayIntersect(edge_cell));
                    Assert.False(bulging.Contains(corner_cell));
                    Assert.False(bulging.MayIntersect(corner_cell));

                    // A singleton cap.
                    S2Cap singleton = new S2Cap(center, S1Angle.Zero);
                    Assert.Equal(cap_face == face, singleton.MayIntersect(root_cell));
                    Assert.False(singleton.MayIntersect(edge_cell));
                    Assert.False(singleton.MayIntersect(corner_cell));
                }
            }
        }

        [Fact]
        public void Test_S2Cap_GetCellUnionBoundLevel1Radius() {
            // Check that a cap whose radius is approximately the width of a level 1
            // S2Cell can be covered by only 3 faces.
            S2Cap cap = new S2Cap(new S2Point(1, 1, 1).Normalized,
            S1Angle.FromRadians(S2Metrics.kMinWidth.GetValue(1)));
            var covering = new List<S2CellId>();
            cap.GetCellUnionBound(covering);
            Assert.Equal(3, covering.Count);
        }

        [Fact]
        public void Test_S2Cap_Expanded() {
            Assert.True(S2Cap.Empty.Expanded(S1Angle.FromRadians(2)).IsEmpty);
            Assert.True(S2Cap.Full.Expanded(S1Angle.FromRadians(2)).IsFull);
            S2Cap cap50 = new S2Cap(new S2Point(1, 0, 0), S1Angle.FromDegrees(50));
            S2Cap cap51 = new S2Cap(new S2Point(1, 0, 0), S1Angle.FromDegrees(51));
            Assert.True(cap50.Expanded(S1Angle.FromRadians(0)).ApproxEquals(cap50));
            Assert.True(cap50.Expanded(S1Angle.FromDegrees(1)).ApproxEquals(cap51));
            Assert.False(cap50.Expanded(S1Angle.FromDegrees(129.99)).IsFull);
            Assert.True(cap50.Expanded(S1Angle.FromDegrees(130.01)).IsFull);
        }

        [Fact]
        public void Test_S2Cap_GetCentroid() {
            // Empty and full caps.
            Assert.Equal(new S2Point(), S2Cap.Empty.GetCentroid());
            Assert.True(S2Cap.Full.GetCentroid().Norm <= S2Constants.DoubleError);

            // Random caps.
            for (int i = 0; i < 100; ++i) {
                S2Point center = S2Testing.RandomPoint();
                double height = S2Testing.Random.UniformDouble(0.0, 2.0);
                S2Cap cap = S2Cap.FromCenterHeight(center, height);
                S2Point centroid = cap.GetCentroid();
                S2Point expected = center * (1.0 - height / 2.0) * cap.Area;
                Assert.True((expected - centroid).Norm <= S2Constants.DoubleError);
            }
        }

        [Fact]
        public void Test_S2Cap_Union() {
            // Two caps which have the same center but one has a larger radius.
            S2Cap a = new S2Cap(GetLatLngPoint(50.0, 10.0), S1Angle.FromDegrees(0.2));
            S2Cap b = new S2Cap(GetLatLngPoint(50.0, 10.0), S1Angle.FromDegrees(0.3));
            Assert.True(b.Contains(a));
            Assert.Equal(b, a.Union(b));

            // Two caps where one is the full cap.
            Assert.True(a.Union(S2Cap.Full).IsFull);

            // Two caps where one is the empty cap.
            Assert.Equal(a, a.Union(S2Cap.Empty));

            // Two caps which have different centers, one entirely encompasses the other.
            S2Cap c = new S2Cap(GetLatLngPoint(51.0, 11.0), S1Angle.FromDegrees(1.5));
            Assert.True(c.Contains(a));
            Assert.Equal(a.Union(c).Center, c.Center);
            Assert.Equal(a.Union(c).Radius, c.Radius);

            // Two entirely disjoint caps.
            S2Cap d = new S2Cap(GetLatLngPoint(51.0, 11.0), S1Angle.FromDegrees(0.1));
            Assert.False(d.Contains(a));
            Assert.False(d.Intersects(a));
            Assert.True(a.Union(d).ApproxEquals(d.Union(a)));
            Assert2.Near(50.4588, new S2LatLng(a.Union(d).Center).Lat.Degrees, 0.001);
            Assert2.Near(10.4525, new S2LatLng(a.Union(d).Center).Lng.Degrees, 0.001);
            Assert2.Near(0.7425, a.Union(d).Radius.Degrees, 0.001);

            // Two partially overlapping caps.
            S2Cap e = new S2Cap(GetLatLngPoint(50.3, 10.3), S1Angle.FromDegrees(0.2));
            Assert.False(e.Contains(a));
            Assert.True(e.Intersects(a));
            Assert.True(a.Union(e).ApproxEquals(e.Union(a)));
            Assert2.Near(50.1500, new S2LatLng(a.Union(e).Center).Lat.Degrees, 0.001);
            Assert2.Near(10.1495, new S2LatLng(a.Union(e).Center).Lng.Degrees, 0.001);
            Assert2.Near(0.3781, a.Union(e).Radius.Degrees, 0.001);

            // Two very large caps, whose radius sums to in excess of 180 degrees, and
            // whose centers are not antipodal.
            S2Cap f = new S2Cap(new S2Point(0, 0, 1).Normalized, S1Angle.FromDegrees(150));
            S2Cap g = new S2Cap(new S2Point(0, 1, 0).Normalized, S1Angle.FromDegrees(150));
            Assert.True(f.Union(g).IsFull);

            // Two non-overlapping hemisphere caps with antipodal centers.
            S2Cap hemi = S2Cap.FromCenterHeight(new S2Point(0, 0, 1).Normalized, 1);
            Assert.True(hemi.Union(hemi.Complement()).IsFull);
        }

        [Fact]
        public void Test_S2Cap_EncodeDecode() {
            S2Cap cap = S2Cap.FromCenterHeight(new S2Point(3, 2, 1).Normalized, 1);
            Encoder encoder = new Encoder();
            cap.Encode(encoder);
            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            var (success, decoded_cap) = S2Cap.DecodeStatic(decoder);
            Assert.True(success);
            Assert.Equal(cap, decoded_cap);
        }

        public static S2Point GetLatLngPoint(double lat_degrees, double lng_degrees)
        {
            return S2LatLng.FromDegrees(lat_degrees, lng_degrees).ToPoint();
        }
    }
}
