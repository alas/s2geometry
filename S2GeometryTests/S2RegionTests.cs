using System;
using System.Collections.Generic;
using Xunit;

namespace S2Geometry
{
    public class S2RegionTests
    {
        #region Strings

        //////////////  These values are in version 1 encoding format.  ////////////////

        // S2Cap.
        private const string kEncodedCapEmpty = "000000000000F03F00000000000000000000000000000000000000000000F0BF";
        private const string kEncodedCapFull = "000000000000F03F000000000000000000000000000000000000000000001040";
        // S2Cap from S2Point(3, 2, 1).Normalize()
        private const string kEncodedCapFromPoint = "3F36105836A8E93F2A2460E5CE1AE13F2A2460E5CE1AD13F0000000000000000";
        // S2Cap from S2Point(0, 0, 1) with height 5
        private const string kEncodedCapFromCenterHeight = "00000000000000000000000000000000000000000000F03F0000000000001040";

        // S2CellId.
        // S2CellId from Face 0.
        private const string kEncodedCellIDFace0 = "0000000000000010";
        // S2CellId from Face 5.
        private const string kEncodedCellIDFace5 = "00000000000000B0";
        // S2CellId from Face 0 in the last S2Cell at kMaxLevel.
        private const string kEncodedCellIDFace0MaxLevel = "0100000000000020";
        // S2CellId from Face 5 in the last S2Cell at kMaxLevel.
        private const string kEncodedCellIDFace5MaxLevel = "01000000000000C0";
        // S2CellId FromFacePosLevel(3, 0x12345678, S2Constants.kMaxCellLevel - 4)
        private const string kEncodedCellIDFacePosLevel = "0057341200000060";
        // S2CellId from the 0 value.
        private const string kEncodedCellIDInvalid = "0000000000000000";

        // S2Cell.
        // S2Cell from S2Point(1, 2, 3)
        private const string kEncodedCellFromPoint = "F51392E0F35DCC43";
        // S2Cell from LatLng(39.0, -120.0) - The Lake Tahoe border corner of CA/NV.
        private const string kEncodedCellFromLatLng = "6308962A95849980";
        // S2Cell from FacePosLevel(3, 0x12345678, S2Constants.kMaxCellLevel - 4)
        private const string kEncodedCellFromFacePosLevel = "0057341200000060";
        // S2Cell from Face 0.
        private const string kEncodedCellFace0 = "0000000000000010";

        // S2CellUnion.
        // An uninitialized empty S2CellUnion.
        private const string kEncodedCellUnionEmpty = "010000000000000000";
        // S2CellUnion from an S2CellId from Face 1.
        private const string kEncodedCellUnionFace1 = "0101000000000000000000000000000030";
        // S2CellUnion from the cells {0x33, 0x8e3748fab, 0x91230abcdef83427};
        private const string kEncodedCellUnionFromCells = "0103000000000000003300000000000000AB8F74E3080000002734F8DEBC0A2391";

        // S2LatLngRect
        private const string kEncodedRectEmpty = "01000000000000F03F0000000000000000182D4454FB210940182D4454FB2109C0";
        private const string kEncodedRectFull = "01182D4454FB21F9BF182D4454FB21F93F182D4454FB2109C0182D4454FB210940";
        // S2LatLngRect from Center=(80,170), Size=(40,60)
        private const string kEncodedRectCentersize = "0165732D3852C1F03F182D4454FB21F93FF75B8A41358C03408744E74A185706C0";

        // S2Loop
        private const string kEncodedLoopEmpty = "010100000000000000000000000000000000000000000000000000F03F0000000000010000"
            + "00000000F03F0000000000000000182D4454FB210940182D4454FB2109C0";
        private const string kEncodedLoopFull = "010100000000000000000000000000000000000000000000000000F0BF010000000001182D"
            + "4454FB21F9BF182D4454FB21F93F182D4454FB2109C0182D4454FB210940";
        // S2Loop from the unit test value kCross1;
        private const string kEncodedLoopCross = "0108000000D44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA1BFB4825F3C81FDEF3F"
            + "27DCF7C958DE913F1EDD892B0BDF91BFB4825F3C81FDEF3F27DCF7C958DE913F1EDD892B0B"
            + "DF913FD44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA13FD44A8442C3F9EF3F7EDA"
            + "2AB341DC91BF27DCF7C958DEA13FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91"
            + "3FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91BFD44A8442C3F9EF3F7EDA2AB3"
            + "41DC91BF27DCF7C958DEA1BF000000000001"
            + "3EFC10E8F8DFA1BF" /*+ "799D52A246DFA1BF"  // S2Loop.bound_.Lat.Lo */
            + "3EFC10E8F8DFA13F"  // S2Loop.bound_.Lat.Hi
            + "389D52A246DF91BF"  // S2Loop.bound_.Lng.Lo
            + "389D52A246DF913F"; // S2Loop.bound_.Lng.Hi


        // S2Loop encoded using EncodeCompressed from snapped points.
        // S2Point[] snapped_loop_a_vertices = {
        //       S2CellId(S2TextFormat.MakePoint("0:178")).ToPoint(),
        //       S2CellId(S2TextFormat.MakePoint("-1:180")).ToPoint(),
        //       S2CellId(S2TextFormat.MakePoint("0:-179")).ToPoint(),
        //       S2CellId(S2TextFormat.MakePoint("1:-180")).ToPoint()};
        // snapped_loop = new S2Loop>(snapped_loop_a_vertices));
        // absl.FixedArray<S2XYZFaceSiTi> points(loop.NumVertices);
        // loop.GetXYZFaceSiTiVertices(points.data());
        // loop.EncodeCompressed(encoder, points.data(), level);
        //
        private const string kEncodedLoopCompressed = "041B02222082A222A806A0C7A991DE86D905D7C3A691F2DEE40383908880A0958805000003";

        // S2PointRegion
        // The difference between an S2PointRegion and an S2Point being encoded is the
        // presence of the encoding version number as the first 8 bits. (i.e. when
        // S2Loop and others encode a stream of S2Points, they are writing triples of
        // doubles instead of encoding the points with the version byte.)
        // S2PointRegion(S2PointUtil.Origin)
        private const string kEncodedPointOrigin = "013BED86AA997A84BF88EC8B48C53C653FACD2721A90FFEF3F";
        // S2PointRegion(S2Point(12.34, 56.78, 9.1011).Normalize())
        private const string kEncodedPointTesting = "0109AD578332DBCA3FBC9FDB9BB4E4EE3FE67E7C2CA7CEC33F";

        // S2Polygon
        // S2Polygon from S2TextFormat.MakePolygon("").
        // This is encoded in compressed format v4.
        private const string kEncodedPolygonEmpty = "041E00";
        // S2Polygon from S2TextFormat.MakePolygon("full").
        // This is encoded in compressed format v4.
        private const string kEncodedPolygonFull = "040001010B000100";
        // S2Polygon from the unit test value kCross1. Encoded in lossless format.
        private const string kEncodedPolygon1Loops = "010100010000000108000000D44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA1BFB4"
    + "825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF91BFB4825F3C81FDEF3F27DCF7C958DE"
    + "913F1EDD892B0BDF913FD44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA13FD44A84"
    + "42C3F9EF3F7EDA2AB341DC91BF27DCF7C958DEA13FB4825F3C81FDEF3F27DCF7C958DE91BF"
    + "1EDD892B0BDF913FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91BFD44A8442C3"
    + "F9EF3F7EDA2AB341DC91BF27DCF7C958DEA1BF000000000001" + "799D52A246DFA1BF"/*"3EFC10E8F8DFA1BF"*/ + "3EFC10E8" // S2Loop.bound_.Lat.Lo
    + "F8DFA13F389D52A246DF91BF389D52A246DF913F01" + "799D52A246DFA1BF"/*"3EFC10E8F8DFA1BF"*/ + "3EFC10E8F8DFA13F" // S2Polygon.bound_.Lat.Lo
    + "389D52A246DF91BF389D52A246DF913F";
        // S2Polygon from the unit test value kCross1+kCrossHole.
        // This is encoded in lossless format.
        private const string kEncodedPolygon2Loops = "010101020000000108000000D44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA1BFB4"
    + "825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF91BFB4825F3C81FDEF3F27DCF7C958DE"
    + "913F1EDD892B0BDF913FD44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA13FD44A84"
    + "42C3F9EF3F7EDA2AB341DC91BF27DCF7C958DEA13FB4825F3C81FDEF3F27DCF7C958DE91BF"
    + "1EDD892B0BDF913FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91BFD44A8442C3"
    + "F9EF3F7EDA2AB341DC91BF27DCF7C958DEA1BF0000000000013EFC10E8F8DFA1BF3EFC10E8"
    + "F8DFA13F389D52A246DF91BF389D52A246DF913F0104000000C5D7FA4B60FFEF3F1EDD892B"
    + "0BDF813F214C95C437DF81BFC5D7FA4B60FFEF3F1EDD892B0BDF813F214C95C437DF813FC5"
    + "D7FA4B60FFEF3F1EDD892B0BDF81BF214C95C437DF813FC5D7FA4B60FFEF3F1EDD892B0BDF"
    + "81BF214C95C437DF81BF000100000001900C5E3B73DF81BF900C5E3B73DF813F399D52A246"
    + "DF81BF399D52A246DF813F013EFC10E8F8DFA1BF3EFC10E8F8DFA13F389D52A246DF91BF38"
    + "9D52A246DF913F";
        // TODO(user): Generate S2Polygons that use compressed encoding format for
        // testing.

        // S2Polyline
        // An S2Polyline from an empty vector.
        private const string kEncodedPolylineEmpty = "0100000000";
        // An S2Polyline from 3 S2LatLngs {(0, 0),(0, 90),(0,180)};
        private const string kEncodedPolylineSemiEquator = "0103000000000000000000F03F00000000000000000000000000000000075C143326A6913C"
            + "000000000000F03F0000000000000000000000000000F0BF075C143326A6A13C0000000000"
            + "000000";
        // An S2Polyline from MakePolyline("0:0, 0:10, 10:20, 20:30");
        private const string kEncodedPolyline3Segments = "0104000000000000000000F03F00000000000000000000000000000000171C818C8B83EF3F"
            + "89730B7E1A3AC63F000000000000000061B46C3A039DED3FE2DC829F868ED53F89730B7E1A"
            + "3AC63F1B995E6FA10AEA3F1B2D5242F611DE3FF50B8A74A8E3D53F";

        #endregion

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2Cap()
        {
            S2Cap cap_from_point = S2Cap.FromPoint(new S2Point(3, 2, 1).Normalized);
            var cap_from_center_height = S2Cap.FromCenterHeight(new S2Point(0, 0, 1).Normalized, 5);

            var cap = TestEncodeDecode(kEncodedCapEmpty, S2Cap.Empty);
            Assert.True(S2Cap.Empty.ApproxEquals(cap));
            cap = TestEncodeDecode(kEncodedCapFull, S2Cap.Full);
            Assert.True(S2Cap.Full.ApproxEquals(cap));
            cap = TestEncodeDecode(kEncodedCapFromPoint, cap_from_point);
            Assert.True(cap_from_point.ApproxEquals(cap));
            cap = TestEncodeDecode(kEncodedCapFromCenterHeight, cap_from_center_height);
            Assert.True(cap_from_center_height.ApproxEquals(cap));
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2Cell() {
            S2Cell cell_from_point = new S2Cell(new S2Point(1, 2, 3));
            S2Cell cell_from_latlng = new S2Cell(S2LatLng.FromDegrees(39.0, -120.0));
            S2Cell cell_from_face_pos_lvl = S2Cell.FromFacePosLevel(3, 0x12345678, S2Constants.kMaxCellLevel - 4);
            S2Cell cell_from_from_face = S2Cell.FromFace(0);

            var cell = TestEncodeDecode(kEncodedCellFromPoint, cell_from_point);
            Assert.Equal(cell_from_point, cell);
            cell = TestEncodeDecode(kEncodedCellFromLatLng, cell_from_latlng);
            Assert.Equal(cell_from_latlng, cell);
            cell = TestEncodeDecode(kEncodedCellFromFacePosLevel, cell_from_face_pos_lvl);
            Assert.Equal(cell_from_face_pos_lvl, cell);
            cell = TestEncodeDecode(kEncodedCellFace0, cell_from_from_face);
            Assert.Equal(cell_from_from_face, cell);
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2CellUnion() {
            S2CellUnion cu_empty = new S2CellUnion();
            S2CellUnion cu_face1 = new S2CellUnion(new List<S2CellId>{ S2CellId.FromFace(1)});
            // Cell ids taken from S2CellUnion EncodeDecode test.
            S2CellUnion cu_latlngs = S2CellUnion.FromNormalized(new List<S2CellId>
                {
                    new S2CellId(0x33),
                    new S2CellId(0x8e3748fab),
                    new S2CellId(0x91230abcdef83427),
                });

            var cu = TestEncodeDecode(kEncodedCellUnionEmpty, cu_empty);
            Assert.Equal(cu_empty, cu);
            cu = TestEncodeDecode(kEncodedCellUnionFace1, cu_face1);
            Assert.Equal(cu_face1, cu);
            cu = TestEncodeDecode(kEncodedCellUnionFromCells, cu_latlngs);
            Assert.Equal(cu_latlngs, cu);
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2LatLngRect()
        {
            var rectFromCenterSize = S2LatLngRect.FromCenterSize(
                S2LatLng.FromDegrees(80, 170), S2LatLng.FromDegrees(40, 60));

            var rect = TestEncodeDecode(kEncodedRectEmpty, S2LatLngRect.Empty);
            Assert.True(S2LatLngRect.Empty.ApproxEquals(rect));
            rect = TestEncodeDecode(kEncodedRectFull, S2LatLngRect.Full);
            Assert.True(S2LatLngRect.Full.ApproxEquals(rect));
            rect = TestEncodeDecode(kEncodedRectCentersize, rectFromCenterSize);
            Assert.True(rectFromCenterSize.ApproxEquals(rect));
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2Loop()
        {
            const string kCross1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1";
            // string kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

            var loop_cross = S2TextFormat.MakeLoopOrDie(kCross1);

            var loop = TestEncodeDecode(kEncodedLoopEmpty, S2Loop.kEmpty);
            Assert.True(S2Loop.kEmpty == loop);
            loop = TestEncodeDecode(kEncodedLoopFull, S2Loop.kFull);
            Assert.True(S2Loop.kFull == loop);
            loop = TestEncodeDecode(kEncodedLoopCross, loop_cross);
            Assert.True(loop_cross == loop);
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2PointRegion() {
            S2PointRegion point_origin = new S2PointRegion(S2PointUtil.Origin);
            S2PointRegion point_testing = new S2PointRegion(new S2Point(12.34, 56.78, 9.1011).Normalized);
            TestEncodeDecode(kEncodedPointOrigin, point_origin);
            TestEncodeDecode(kEncodedPointTesting, point_testing);
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2Polygon() {
            string kCross1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1";
            string kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

            var polygon_empty = S2TextFormat.MakePolygonOrDie("");
            var polygon_full = S2TextFormat.MakeVerbatimPolygonOrDie("full");
            var polygon_cross = S2TextFormat.MakePolygonOrDie(kCross1);
            var polygon_cross_hole = S2TextFormat.MakePolygonOrDie(kCross1 + ";" + kCrossCenterHole);

            var polygon = TestEncodeDecode(kEncodedPolygonEmpty, polygon_empty);
            Assert.True(polygon_empty == polygon);
            polygon = TestEncodeDecode(kEncodedPolygonFull, polygon_full);
            Assert.True(polygon_full == polygon);
            polygon = TestEncodeDecode(kEncodedPolygon1Loops, polygon_cross);
            Assert.True(polygon_cross == polygon);
            polygon = TestEncodeDecode(kEncodedPolygon2Loops, polygon_cross_hole);
            Assert.True(polygon_cross_hole == polygon);
        }

        [Fact]
        public void Test_S2RegionEncodeDecodeTest_S2Polyline()
        {
            var vertices = Array.Empty<S2Point>();
            var polyline_empty = new S2Polyline(vertices);
            var polyline_semi = new S2Polyline(new[]
                {
                    S2LatLng.FromDegrees(0, 0),
                    S2LatLng.FromDegrees(0, 90),
                    S2LatLng.FromDegrees(0, 180),
                });
            var polyline_3segments = S2TextFormat.MakePolylineOrDie("0:0, 0:10, 10:20, 20:30");

            var polyline = TestEncodeDecode(kEncodedPolylineEmpty, polyline_empty);
            Assert.True(polyline_empty == polyline);
            polyline = TestEncodeDecode(kEncodedPolylineSemiEquator, polyline_semi);
            Assert.True(polyline_semi == polyline);
            polyline = TestEncodeDecode(kEncodedPolyline3Segments, polyline_3segments);
            Assert.True(polyline_3segments == polyline);
        }

        // Encode/Decode not yet implemented for these types.
        // S2R2Rect
        // S2RegionIntersection
        // S2RegionUnion
        // TestEncodeDecode tests that the input encodes to match the expected
        // golden data, and then returns the decode of the data into dst.
        private static Region TestEncodeDecode<Region>(string golden, Region src) where Region : ICoder, IS2Region
        {
            Encoder encoder = new();
            src.Encode(encoder);

            var str = encoder.Buffer.ToHexa(encoder.Length);
            Assert.Equal(golden, str);

            var decoder = new Decoder(encoder.Buffer, 0, encoder.Buffer.Length);
            var funcDict = new Dictionary<Type, Func<Decoder, (bool success, ICoder result)>>
                {
                    { typeof(S2Cap), (d) => S2Cap.DecodeStatic(d) },
                    { typeof(S2Cell), (d) => S2Cell.DecodeStatic(d) },
                    { typeof(S2CellUnion), (d) => S2CellUnion.DecodeStatic(d) },
                    { typeof(S2LatLngRect), (d) => S2LatLngRect.DecodeStatic(d) },
                    { typeof(S2Loop), (d) => S2Loop.DecodeStatic(d) },
                    { typeof(S2PointRegion), (d) => S2PointRegion.DecodeStatic(d) },
                    { typeof(S2Polygon), (d) => S2Polygon.DecodeStatic(d) },
                    { typeof(S2Polyline), (d) => S2Polyline.DecodeStatic(d) },
                    { typeof(S2CellId), (d) => S2CellId.DecodeStatic(d) },
                };
            if (funcDict.ContainsKey(typeof(Region)))
            {
                var func = funcDict[typeof(Region)];
                var (success, result) = func(decoder);
                if (success) return (Region)result;
            }

            var mi = typeof(Region).GetMember("DecodeStatic").GetValue(0) as System.Reflection.MethodInfo;
            var (success_r, result_r) = ((bool, Region))mi.Invoke(null, new[] { decoder });
            if (success_r) return result_r;

            throw new NotImplementedException($"TestEncodeDecode for type: {typeof(Region).FullName}");
        }

        // TODO(user): When the remaining types implement Encode/Decode, add their
        // test cases here. i.e. S2R2Rect, S2RegionIntersection, S2RegionUnion, etc.
    }
}
