namespace S2Geometry;

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

    // S2PointRegion
    // The difference between an S2PointRegion and an S2Point being encoded is the
    // presence of the encoding version number as the first 8 bits. (i.e. when
    // S2Loop and others encode a stream of S2Points, they are writing triples of
    // doubles instead of encoding the points with the version byte.)
    // S2PointRegion(S2.Origin)
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
    internal void Test_S2RegionEncodeDecodeTest_S2Cap()
    {
        S2Cap cap_from_point = S2Cap.FromPoint(new S2Point(3, 2, 1).Normalize());
        var cap_from_center_height = S2Cap.FromCenterHeight(new S2Point(0, 0, 1).Normalize(), 5);

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
    internal void Test_S2RegionEncodeDecodeTest_S2Cell()
    {
        S2Cell cell_from_point = new(new S2Point(1, 2, 3));
        S2Cell cell_from_latlng = new(S2LatLng.FromDegrees(39.0, -120.0));
        S2Cell cell_from_face_pos_lvl = S2Cell.FromFacePosLevel(3, 0x12345678, S2.kMaxCellLevel - 4);
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
    internal void Test_S2RegionEncodeDecodeTest_S2CellUnion()
    {
        S2CellUnion cu_empty = new();
        S2CellUnion cu_face1 = new(new List<S2CellId>{ S2CellId.FromFace(1)});
        // Cell ids taken from S2CellUnion EncodeDecode test.
        S2CellUnion cu_latlngs = S2CellUnion.FromNormalized(
            [
                new S2CellId(0x33),
                new S2CellId(0x8e3748fab),
                new S2CellId(0x91230abcdef83427),
            ]);

        var cu = TestEncodeDecode(kEncodedCellUnionEmpty, cu_empty);
        Assert.Equal(cu_empty, cu);
        cu = TestEncodeDecode(kEncodedCellUnionFace1, cu_face1);
        Assert.Equal(cu_face1, cu);
        cu = TestEncodeDecode(kEncodedCellUnionFromCells, cu_latlngs);
        Assert.Equal(cu_latlngs, cu);
    }

    [Fact]
    internal void Test_S2RegionEncodeDecodeTest_S2LatLngRect()
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
    internal void Test_S2RegionEncodeDecodeTest_S2Loop()
    {
        const string kCross1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1";
        //const string kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

        S2Loop loop_empty = S2Loop.KEmpty;
        S2Loop loop_full = S2Loop.KFull;
        var loop_cross = MakeLoopOrDie(kCross1);

        var loop = TestEncodeDecode(kEncodedLoopEmpty, loop_empty);
        Assert.True(loop_empty == loop);
        loop = TestEncodeDecode(kEncodedLoopFull, loop_full);
        Assert.True(loop_full == loop);
        loop = TestEncodeDecode(kEncodedLoopCross, loop_cross);
        Assert.True(loop_cross == loop);
    }

    [Fact]
    internal void Test_S2RegionEncodeDecodeTest_S2PointRegion()
    {
        S2PointRegion point_origin = new(S2.Origin);
        S2PointRegion point_testing = new(new S2Point(12.34, 56.78, 9.1011).Normalize());
        TestEncodeDecode(kEncodedPointOrigin, point_origin);
        TestEncodeDecode(kEncodedPointTesting, point_testing);
    }

    [Fact]
    internal void Test_S2RegionEncodeDecodeTest_S2Polygon()
    {
        string kCross1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1";
        string kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

        var polygon_empty = MakePolygonOrDie("");
        var polygon_full = MakeVerbatimPolygonOrDie("full");
        var polygon_cross = MakePolygonOrDie(kCross1);
        var polygon_cross_hole = MakePolygonOrDie(kCross1 + ";" + kCrossCenterHole);

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
    internal void Test_S2RegionEncodeDecodeTest_S2Polyline()
    {
        var vertices = Array.Empty<S2Point>();
        var polyline_empty = new S2Polyline(vertices);
        var polyline_semi = new S2Polyline(
            [
                S2LatLng.FromDegrees(0, 0),
                S2LatLng.FromDegrees(0, 90),
                S2LatLng.FromDegrees(0, 180),
            ]);
        var polyline_3segments = MakePolylineOrDie("0:0, 0:10, 10:20, 20:30");

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
    private static Region TestEncodeDecode<Region>(string golden, Region src) where Region : IDecoder<Region>
    {
        Encoder encoder = new();
        src.Encode(encoder);

        Assert.Equal(golden, encoder.HexString());

        var decoder = encoder.GetDecoder();
        var (success, dst) = Region.Decode(decoder);
        Assert.True(success);
        return dst!;
    }

    // TODO(user): When the remaining types implement Encode/Decode, add their
    // test cases here. i.e. S2R2Rect, S2RegionIntersection, S2RegionUnion, etc.
}
