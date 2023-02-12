namespace S2Geometry;

public class S2PointCompressionTests
{
    #region Fields, Constants

#pragma warning disable IDE0051 // Quitar miembros privados no utilizados
    private const int s2point_compression_bm_level = 30; // Level to encode at for benchmarks.
    private const double s2point_compression_bm_radius_km = 1000.0; // Radius to use for loop for benchmarks.
#pragma warning restore IDE0051 // Quitar miembros privados no utilizados

    private readonly Encoder encoder_ = new();

    // Four vertex loop near the corner of faces 0, 1, and 2.
    private readonly S2Point[] loop_4_;
    // Four vertex loop near the corner of faces 0, 1, and 2;
    // unsnapped.
    private readonly S2Point[] loop_4_unsnapped_;
    // Four vertex loop near the corner of faces 0, 1, and 2;
    // snapped to level 14.
    private readonly S2Point[] loop_4_level_14_;
    // 100 vertex loop near the corner of faces 0, 1, and 2.
    private readonly S2Point[] loop_100_;
    // 100 vertex loop near the corner of faces 0, 1, and 2;
    // unsnapped.
    private readonly S2Point[] loop_100_unsnapped_;
    // 100 vertex loop near the corner of faces 0, 1, and 2;
    // 15 points snapped to kMakCellLevel, the others not snapped.
    private readonly S2Point[] loop_100_mixed_15_;
    // 100 vertex loop near the corner of faces 0, 1, and 2;
    // 25 points snapped to kMakCellLevel, the others not snapped.
    private readonly S2Point[] loop_100_mixed_25_;
    // 100 vertex loop near the corner of faces 0, 1, and 2;
    // snapped to level 22.
    private readonly S2Point[] loop_100_level_22_;
    // A loop with two vertices on each of three faces.
    private readonly S2Point[] loop_multi_face_;
    // A straight line of 100 vertices on face 0 that should compress well.
    private readonly S2Point[] line_;

    #endregion

    public S2PointCompressionTests()
    {
        loop_4_ = MakeRegularPoints(4, 0.1, S2.kMaxCellLevel);

        S2Point center = new S2Point(1.0, 1.0, 1.0).Normalize();
        S1Angle radius = S2Testing.KmToAngle(0.1);
        loop_4_unsnapped_ = S2Testing.MakeRegularPoints(center, radius, 4);

        // Radius is 100m, so points are about 141 meters apart.
        // Snapping to level 14 will move them by < 47m.
        loop_4_level_14_ = MakeRegularPoints(4, 0.1, 14);

        loop_100_ = MakeRegularPoints(100, 0.1, S2.kMaxCellLevel);

        loop_100_unsnapped_ = S2Testing.MakeRegularPoints(center, radius, 100);

        loop_100_mixed_15_ = S2Testing.MakeRegularPoints(center, radius, 100);
        for (int i = 0; i < 15; ++i)
        {
            loop_100_mixed_15_[3 * i] = SnapPointToLevel(loop_100_mixed_15_[3 * i],
                                                         S2.kMaxCellLevel);
        }

        loop_100_mixed_25_ = S2Testing.MakeRegularPoints(center, radius, 100);
        for (int i = 0; i < 25; ++i)
        {
            loop_100_mixed_25_[4 * i] = SnapPointToLevel(loop_100_mixed_25_[4 * i],
                                                         S2.kMaxCellLevel);
        }

        // Circumference is 628m, so points are about 6 meters apart.
        // Snapping to level 22 will move them by < 2m.
        loop_100_level_22_ = MakeRegularPoints(100, 0.1, 22);

        var multi_face_points = new S2Point[6];
        multi_face_points[0] = S2.FaceUVtoXYZ(0, -0.5, 0.5).Normalize();
        multi_face_points[1] = S2.FaceUVtoXYZ(1, -0.5, 0.5).Normalize();
        multi_face_points[2] = S2.FaceUVtoXYZ(1, 0.5, -0.5).Normalize();
        multi_face_points[3] = S2.FaceUVtoXYZ(2, -0.5, 0.5).Normalize();
        multi_face_points[4] = S2.FaceUVtoXYZ(2, 0.5, -0.5).Normalize();
        multi_face_points[5] = S2.FaceUVtoXYZ(2, 0.5, 0.5).Normalize();
        loop_multi_face_ = SnapPointsToLevel(multi_face_points, S2.kMaxCellLevel);

        var line_points = new S2Point[100];
        for (int i = 0; i < line_points.Length; ++i)
        {
            double s = 0.01 + 0.005 * i;
            double t = 0.01 + 0.009 * i;
            double u = S2.STtoUV(s);
            double v = S2.STtoUV(t);
            line_points[i] = S2.FaceUVtoXYZ(0, u, v).Normalize();
        }
        line_ = SnapPointsToLevel(line_points, S2.kMaxCellLevel);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_RoundtripsEmpty()
    {
        // Just check this doesn't crash.
        Encode(Array.Empty<S2Point>(), S2.kMaxCellLevel);
        Decode(S2.kMaxCellLevel, Array.Empty<S2Point>());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_RoundtripsFourVertexLoop()
    {
        Roundtrip(loop_4_, S2.kMaxCellLevel);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_RoundtripsFourVertexLoopUnsnapped()
    {
        Roundtrip(loop_4_unsnapped_, S2.kMaxCellLevel);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_FourVertexLoopSize()
    {
        Encode(loop_4_, S2.kMaxCellLevel);
        // It would take 32 bytes uncompressed.
        Assert.Equal(39, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_RoundtripsFourVertexLevel14Loop()
    {
        int level = 14;
        Roundtrip(loop_4_level_14_, level);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_FourVertexLevel14LoopSize()
    {
        int level = 14;
        Encode(loop_4_level_14_, level);
        // It would take 4 bytes per vertex without compression.
        Assert.Equal(23, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_Roundtrips100VertexLoop()
    {
        Roundtrip(loop_100_, S2.kMaxCellLevel);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_Roundtrips100VertexLoopUnsnapped()
    {
        Roundtrip(loop_100_unsnapped_, S2.kMaxCellLevel);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_Roundtrips100VertexLoopMixed15()
    {
        Roundtrip(loop_100_mixed_15_, S2.kMaxCellLevel);
        Assert.Equal(2381, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_Roundtrips100VertexLoopMixed25()
    {
        Roundtrip(loop_100_mixed_25_, S2.kMaxCellLevel);
        Assert.Equal(2131, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_OneHundredVertexLoopSize()
    {
        Encode(loop_100_, S2.kMaxCellLevel);
        Assert.Equal(257, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_OneHundredVertexLoopUnsnappedSize()
    {
        Encode(loop_100_unsnapped_, S2.kMaxCellLevel);
        Assert.Equal(2756, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_Roundtrips100VertexLevel22Loop()
    {
        int level = 22;
        Roundtrip(loop_100_level_22_, level);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_OneHundredVertexLoopLevel22Size()
    {
        Encode(loop_100_level_22_, 22);
        Assert.Equal(148, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_MultiFaceLoop()
    {
        Roundtrip(loop_multi_face_, S2.kMaxCellLevel);
    }

    [Fact]
    internal void Test_S2PointCompressionTest_StraightLineCompressesWell()
    {
        Roundtrip(line_, S2.kMaxCellLevel);
        // About 1 byte / vertex.
        Assert.Equal(line_.Length + 17, encoder_.Length());
    }

    [Fact]
    internal void Test_S2PointCompressionTest_FirstPointOnFaceEdge()
    {
        // This test used to trigger a bug in which EncodeFirstPointFixedLength()
        // tried to encode a pi/qi value of (2**level) in "level" bits (which did
        // not work out so well).  The fix is documented in SiTitoPiQi().
        //
        // The test data consists of two points, where the first point is exactly on
        // an S2Cell face edge (with ti == S2.kMaxSiTi), and the second point is
        // encodable at snap level 8.  This used to cause the code to try encoding
        // qi = 256 in 8 bits.
        var points = new S2PointCompression.S2XYZFaceSiTi[]
        {
            new S2PointCompression.S2XYZFaceSiTi(
                new S2Point(0.054299323861222645, -0.70606358900180299, 0.70606358900180299),
                2,                      // face
                956301312, 2147483648,  // si, ti
                -1                      // level
            ),
            new S2PointCompression.S2XYZFaceSiTi(
                new S2Point(0.056482651436986935, -0.70781701406865505, 0.70413406726388494),
                4,                    // face
                4194304, 1195376640,  // si, ti
                8                     // level
            )
        };

        Encoder encoder = new();
        S2PointCompression.S2EncodePointsCompressed(points, 8, encoder);
        var decoder = encoder.GetDecoder();
        var result = new S2Point[2];
        Assert.True(S2PointCompression.S2DecodePointsCompressed(decoder, 8, result, 0));
        Assert.True(result[0] == points[0].XYZ);
        Assert.True(result[1] == points[1].XYZ);
    }

    private static S2Point SnapPointToLevel(S2Point point, int level) =>
        new S2CellId(point).Parent(level).ToPoint();

    private static S2Point[] SnapPointsToLevel(S2Point[] points, int level)
    {
        var snapped_points = new S2Point[points.Length];
        for (int i = 0; i < points.Length; ++i)
        {
            snapped_points[i] = SnapPointToLevel(points[i], level);
        }
        return snapped_points;
    }

    // Make a regular loop around the corner of faces 0, 1, and 2 with the
    // specified radius in meters (on the earth) and number of vertices.
    private static S2Point[] MakeRegularPoints(int num_vertices, double radius_km, int level)
    {
        var center = new S2Point(1.0, 1.0, 1.0).Normalize();
        var radius_angle = S2Testing.KmToAngle(radius_km);

        var unsnapped_points = S2Testing.MakeRegularPoints(center, radius_angle, num_vertices);

        return SnapPointsToLevel(unsnapped_points.ToArray(), level);
    }

    private static S2PointCompression.S2XYZFaceSiTi[] MakeXYZFaceSiTiPoints(S2Point[] points) =>
        points.Select(t =>
        {
            var cellLevel = S2.XYZtoFaceSiTi(t, out var face, out var si, out var ti);
            return new S2PointCompression.S2XYZFaceSiTi(t, face, si, ti, cellLevel);
        }).ToArray();

    private void Encode(S2Point[] points, int level)
    {
        var pts = MakeXYZFaceSiTiPoints(points);
        S2PointCompression.S2EncodePointsCompressed(pts, level, encoder_);
    }

    private void Decode(int level, S2Point[] points)
    {
        var decoder_ = encoder_.GetDecoder();
        Assert.True(S2PointCompression.S2DecodePointsCompressed(decoder_, level, points, 0));
    }

    private void Roundtrip(S2Point[] loop, int level)
    {
        Encode(loop, level);
        var points = new S2Point[loop.Length];
        Decode(level, points);
        Assert.True(loop.SequenceEqual(points)); // "Decoded points" + points.ToDebugString() + "do not match original points" + loop.ToDebugString();
    }
}
