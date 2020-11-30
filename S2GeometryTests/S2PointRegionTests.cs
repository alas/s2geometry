using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;
using Xunit.Abstractions;

namespace S2Geometry
{
    public class S2PointRegionTests
    {
        [Fact]
        public void Test_S2PointRegionTest_Basic() {
            S2Point p = new S2Point(1, 0, 0);
            S2PointRegion r0 = new S2PointRegion(p);
            Assert.Equal(r0.Point, p);
            Assert.True(r0.Contains(p));
            Assert.True(r0.Contains(r0.Point));
            Assert.False(r0.Contains(new S2Point(1, 0, 1)));
            S2PointRegion r0_clone = (S2PointRegion)r0.Clone();
            Assert.Equal(r0_clone.Point, r0.Point);
            Assert.Equal(r0.GetCapBound(), S2Cap.FromPoint(p));
            S2LatLng ll = new S2LatLng(p);
            Assert.Equal(r0.GetRectBound(), new S2LatLngRect(ll, ll));

            // The leaf cell containing a point is still much larger than the point.
            S2Cell cell = new S2Cell(p);
            Assert.False(r0.Contains(cell));
            Assert.True(r0.MayIntersect(cell));
        }

        [Fact]
        public void Test_S2PointRegionTest_EncodeAndDecode() {
            S2Point p = new S2Point(0.6, 0.8, 0);
            S2PointRegion r = new S2PointRegion(p);

            Encoder encoder = new Encoder();
            r.Encode(encoder);

            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            // S2PointRegion decoded_r = new S2PointRegion(new S2Point(1, 0, 0));
            var (_, decoded_r) = S2PointRegion.DecodeStatic(decoder);
            S2Point decoded_p = decoded_r.Point;

            Assert.Equal(0.6, decoded_p[0]);
            Assert.Equal(0.8, decoded_p[1]);
            Assert.Equal(0.0, decoded_p[2]);
        }

        [Fact]
        public void Test_S2PointRegionTest_DecodeUnitLength() {
            // Ensure that we have the right format for the next test.
            Encoder encoder = new Encoder();
            encoder.Ensure(1 + 3 * 8);

            encoder.Put8(1);  // version number
            encoder.PutDouble(1.0);
            encoder.PutDouble(0.0);
            encoder.PutDouble(0.0);

            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            // S2PointRegion r = new S2PointRegion(new S2Point(-1, 0, 0));
            var (success, _) = S2PointRegion.DecodeStatic(decoder);
            Assert.True(success);
        }

        [Fact]
        public void Test_S2PointRegionTest_DecodeNonUnitLength() {
            // Ensure that a decode of a non-unit length vector returns false,
            // rather than S2_DCHECK fails.
            Encoder encoder = new Encoder();
            encoder.Ensure(1 + 3 * 8);

            encoder.Put8(1);  // version number
            encoder.PutDouble(1.0);
            encoder.PutDouble(1.0);
            encoder.PutDouble(1.0);

            Decoder decoder = new Decoder(encoder.Buffer, 0, encoder.Length);
            // S2PointRegion r = new S2PointRegion(new S2Point(-1, 0, 0));
            var (success, _) = S2PointRegion.DecodeStatic(decoder);
            Assert.False(success);
        }
    }
}
