namespace S2Geometry;

public class EncodedS2PointVectorTests
{
    private const int kBlockSize = 16;  // Number of deltas per block in implementation.
    private readonly ITestOutputHelper _logger;

    public EncodedS2PointVectorTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_Empty()
    {
        TestEncodedS2PointVector(Array.Empty<S2Point>(), CodingHint.FAST, 1);

        // Test that an empty vector uses the UNCOMPRESSED encoding.
        TestEncodedS2PointVector(Array.Empty<S2Point>(), CodingHint.COMPACT, 1);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_OnePoint()
    {
        TestEncodedS2PointVector(new S2Point[] { new S2Point(1, 0, 0) }, CodingHint.FAST, 25);

        // Encoding: header (2 bytes), block count (1 byte), block offsets (1 byte),
        // block header (1 byte), delta (1 byte).
        TestEncodedS2PointVector(new S2Point[] { new S2Point(1, 0, 0) }, CodingHint.COMPACT, 6);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_OnePointWithExceptionsNoOverlap()
    {
        // Test encoding a block with one point when other blocks have exceptions
        // (which changes the encoding for all blocks).  The case below yields
        // delta_bits == 8 and overlap_bits == 0.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
        // Block 0: block header (1 byte), 16 deltas (16 bytes), exception (24 bytes)
        // Block 1: block header (1 byte), delta (1 byte)
        S2Point a = new(1, 0, 0);
        S2Point[] points = {
                new S2Point(1, 2, 3).Normalize(), a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                a  // Second block
            };
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 48);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_OnePointWithExceptionsWithOverlap()
    {
        // Test encoding a block with one point when other blocks have exceptions
        // (which changes the encoding for all blocks).  The case below yields
        // delta_bits == 8 and overlap_bits == 4.
        //
        // Encoding: header (2 bytes), base (2 bytes), block count (1 byte),
        //           block offsets (2 bytes)
        // Block 0: header (1 byte), offset (2 bytes), 16 deltas (16 bytes),
        //          exception (24 bytes)
        // Block 1: header (1 byte), offset (2 bytes), delta (1 byte)
        S2Point a = new S2CellId(0x946df618d0000000).ToPoint();
        S2Point b = new S2CellId(0x947209e070000000).ToPoint();
        S2Point[] points = {
                new S2Point(1, 2, 3).Normalize(), a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                b  // Second block
             };
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 54);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_CellIdWithException()
    {
        // Test one point encoded as an S2CellId with one point encoded as an
        // exception.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (1 byte),
        // block header (1 byte), two deltas (2 bytes), exception (24 bytes).
        TestEncodedS2PointVector(new S2Point[]
            { MakeCellIdOrDie("1/23").ToPoint(), new S2Point(0.1, 0.2, 0.3).Normalize()},
            CodingHint.COMPACT, 31);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_PointsAtMultipleLevels()
    {
        // Test that when points at multiple levels are present, the level with the
        // most points is chosen (preferring the smallest level in case of ties).
        // (All other points are encoded as exceptions.)

        // In this example, the two points at level 5 (on face 1) should be encoded.
        // It is possible to tell which points are encoded by the length of the
        // encoding (since different numbers of "base" bytes are encoded).
        //
        // Encoding: header (2 bytes), base (1 byte), block count (1 byte), block
        // offsets (1 byte), block header (1 byte), 5 deltas (5 bytes), S2Point
        // exceptions (72 bytes).
        TestEncodedS2PointVector(new S2Point[]
            {
                    MakeCellIdOrDie("2/11001310230102").ToPoint(),
                    MakeCellIdOrDie("1/23322").ToPoint(),
                    MakeCellIdOrDie("3/3").ToPoint(),
                    MakeCellIdOrDie("1/23323").ToPoint(),
                    MakeCellIdOrDie("2/12101023022012").ToPoint()
            },
            CodingHint.COMPACT, 83);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_NoOverlapOrExtraDeltaBitsNeeded()
    {
        // This function tests the case in GetBlockCodes() where values can be
        // encoded using the minimum number delta bits and no overlap.  From the
        // comments there:
        //
        //   Example 1: d_min = 0x72, d_max = 0x7e.  The range is 0x0c.  This can be
        //   encoded using delta_bits = 4 and overlap_bits = 0, which allows us to
        //   represent an offset of 0x70 and a maximum delta of 0x0f, so that we can
        //   encode values up to 0x7f.
        //
        // To set up this test, we need at least two blocks: one to set the global
        // minimum value, and the other to encode a specific range of deltas.  To
        // make things easier, the first block has a minimum value of zero.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
        // Block 0: header (1 byte), 8 deltas (8 bytes)
        // Block 1: header (1 byte), offset (1 byte), 4 deltas (2 bytes)
        const int level = 3;
        var points = new S2Point[kBlockSize]
            .Fill(EncodedValueToPoint(0, level));
        points = new List<S2Point>(points)
            {
                EncodedValueToPoint(0x72, level),
                EncodedValueToPoint(0x74, level),
                EncodedValueToPoint(0x75, level),
                EncodedValueToPoint(0x7e, level)
            }.ToArray();
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 10 + kBlockSize / 2);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_OverlapNeeded()
    {
        // Like the above, but tests the following case:
        //
        //   Example 2: d_min = 0x78, d_max = 0x84.  The range is 0x0c, but in this
        //   case it is not sufficient to use delta_bits = 4 and overlap_bits = 0
        //   because we can again only represent an offset of 0x70, so the maximum
        //   delta of 0x0f only lets us encode values up to 0x7f.  However if we
        //   increase the overlap to 4 bits then we can represent an offset of 0x78,
        //   which lets us encode values up to 0x78 + 0x0f = 0x87.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
        // Block 0: header (1 byte), 8 deltas (8 bytes)
        // Block 1: header (1 byte), offset (1 byte), 4 deltas (2 bytes)
        const int level = 3;
        var points = new S2Point[kBlockSize]
            .Fill(EncodedValueToPoint(0, level));
        points = new List<S2Point>(points)
            {
                EncodedValueToPoint(0x78, level),
                EncodedValueToPoint(0x7a, level),
                EncodedValueToPoint(0x7c, level),
                EncodedValueToPoint(0x84, level),
            }.ToArray();
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 10 + kBlockSize / 2);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_ExtraDeltaBitsNeeded()
    {
        // Like the above, but tests the following case:
        //
        //   Example 3: d_min = 0x08, d_max = 0x104.  The range is 0xfc, so we should
        //   be able to use 8-bit deltas.  But even with a 4-bit overlap, we can still
        //   only encode offset = 0 and a maximum value of 0xff.  (We don't allow
        //   bigger overlaps because statistically they are not worthwhile.)  Instead
        //   we increase the delta size to 12 bits, which handles this case easily.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
        // Block 0: header (1 byte), 8 deltas (8 bytes)
        // Block 1: header (1 byte), 4 deltas (6 bytes)
        const int level = 3;
        var points = new S2Point[kBlockSize]
            .Fill(EncodedValueToPoint(0, level));
        points = new List<S2Point>(points)
            {
                EncodedValueToPoint(0x08, level),
                EncodedValueToPoint(0x4e, level),
                EncodedValueToPoint(0x82, level),
                EncodedValueToPoint(0x104, level),
            }.ToArray();
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 13 + kBlockSize / 2);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_ExtraDeltaBitsAndOverlapNeeded()
    {
        // Like the above, but tests the following case:
        //
        //   Example 4: d_min = 0xf08, d_max = 0x1004.  The range is 0xfc, so we
        //   should be able to use 8-bit deltas.  With 8-bit deltas and no overlap, we
        //   have offset = 0xf00 and a maximum encodable value of 0xfff.  With 8-bit
        //   deltas and a 4-bit overlap, we still have offset = 0xf00 and a maximum
        //   encodable value of 0xfff.  Even with 12-bit deltas, we have offset = 0
        //   and we can still only represent 0xfff.  However with delta_bits = 12 and
        //   overlap_bits = 4, we can represent offset = 0xf00 and a maximum encodable
        //   value of 0xf00 + 0xfff = 0x1eff.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
        // Block 0: header (1 byte), 8 deltas (8 bytes)
        // Block 1: header (1 byte), offset (1 byte), 4 deltas (6 bytes)
        const int level = 5;
        var ini = new S2Point[kBlockSize]
            .Fill(EncodedValueToPoint(0, level));
        var points = new List<S2Point>(ini)
            {
                EncodedValueToPoint(0xf08, level),
                EncodedValueToPoint(0xf4e, level),
                EncodedValueToPoint(0xf82, level),
                EncodedValueToPoint(0x1004, level),
            }.ToArray();
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 14 + kBlockSize / 2);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_SixtyFourBitOffset()
    {
        // Tests a case where a 64-bit block offset is needed.
        //
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
        // Block 0: header (1 byte), 8 deltas (8 bytes)
        // Block 1: header (1 byte), offset (8 bytes), 2 deltas (1 byte)
        const int level = S2.kMaxCellLevel;
        var points = new S2Point[kBlockSize]
            .Fill(S2CellId.Begin(level).ToPoint());
        points = new List<S2Point>(points)
            {
                S2CellId.End(level).Prev().ToPoint(),
                S2CellId.End(level).Prev().Prev().ToPoint(),
            }.ToArray();
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 16 + kBlockSize / 2);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_AllExceptionsBlock()
    {
        // The encoding consists of two blocks; the first contains 16 encodable
        // values, while the second contains two exceptions.
        var points = new S2Point[kBlockSize].Fill(EncodedValueToPoint(0, S2.kMaxCellLevel));
        points = new List<S2Point>(points)
            {
                new S2Point(0.1, 0.2, 0.3).Normalize(),
                new S2Point(0.3, 0.2, 0.1).Normalize(),
            }.ToArray();
        // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes).
        // 1st block header (1 byte), 16 deltas (16 bytes).
        // 2nd block header (1 byte), 2 deltas (1 byte), 2 exceptions (48 bytes).
        TestEncodedS2PointVector(points, CodingHint.COMPACT, 72);

        // Encoding: header (2 bytes), 18 S2Points (432 bytes).
        TestEncodedS2PointVector(points, CodingHint.FAST, 434);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_FirstAtAllLevels()
    {
        // Test encoding the first S2CellId at each level (which also happens to have
        // the maximum face, si, and ti values).  All such S2CellIds can be encoded in
        // 6 bytes because most of the bits are zero.
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            _logger.WriteLine("Level = " + level);
            TestEncodedS2PointVector(new S2Point[] { S2CellId.Begin(level).ToPoint() },
                                     CodingHint.COMPACT, 6);
        }
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_LastAtAllLevels()
    {
        // Test encoding the last S2CellId at each level.  It turns out that such
        // S2CellIds have the largest possible face and ti values, and the minimum
        // possible si value at that level.  Such S2CellIds can be encoded in 6 to 13
        // bytes depending on the level.
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            _logger.WriteLine("Level = " + level);
            // Note that 8 bit deltas are used to encode blocks of size 1, which
            // reduces the size of "base" from ((level + 2) / 4) to (level / 4) bytes.
            int expected_size = 6 + level / 4;
            TestEncodedS2PointVector(new S2Point[] { S2CellId.End(level).Prev().ToPoint() },
                CodingHint.COMPACT, expected_size);
        }
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_MaxFaceSiTiAtAllLevels()
    {
        // Similar to the test above, but tests encoding the S2CellId at each level
        // whose face, si, and ti values are all maximal.  This turns out to be the
        // S2CellId whose human-readable form is 5/222...22 (0xb555555555555555),
        // however for clarity we construct it using S2CellId.FromFaceIJ.
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            _logger.WriteLine("Level = " + level);
            S2CellId id = S2CellId.FromFaceIJ(5, S2.kLimitIJ - 1, S2.kLimitIJ - 1)
                          .Parent(level);

            // This encoding is one byte bigger than the previous test at levels 7, 11,
            // 15, 19, 23, and 27.  This is because in the previous test, the
            // odd-numbered value bits are all zero (except for the face number), which
            // reduces the number of base bits needed by exactly 1.  The encoding size
            // at level==3 is unaffected because for singleton blocks, the lowest 8
            // value bits are encoded in the delta.
            int expected_size = (level < 4) ? 6 : 6 + (level + 1) / 4;
            TestEncodedS2PointVector(new S2Point[] { id.ToPoint() },
                                     CodingHint.COMPACT, expected_size);
        }
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_LastTwoPointsAtAllLevels()
    {
        // Test encoding the last two S2CellIds at each level.
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            _logger.WriteLine("Level = " + level);
            S2CellId id = S2CellId.End(level).Prev();
            // Notice that this costs only 4 bits more than encoding the last S2CellId
            // by itself (see LastAtAllLevels).  This is because encoding a block of
            // size 1 uses 8-bit deltas (which reduces the size of "base" by 4 bits),
            // while this test uses two 4-bit deltas.
            int expected_size = 6 + (level + 2) / 4;
            TestEncodedS2PointVector(new S2Point[] { id.ToPoint(), id.Prev().ToPoint() },
                                     CodingHint.COMPACT, expected_size);
        }
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_ManyDuplicatePointsAtAllLevels()
    {
        // Test encoding 32 copies of the last S2CellId at each level.  This uses
        // between 27 and 38 bytes depending on the level.  (Note that the encoding
        // can use less than 1 byte per point in this situation.)
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            _logger.WriteLine("Level = " + level);
            S2CellId id = S2CellId.End(level).Prev();
            // Encoding: header (2 bytes), base ((level + 2) / 4 bytes), block count
            // (1 byte), block offsets (2 bytes), block headers (2 bytes), 32 deltas
            // (16 bytes).  At level 30 the encoding size goes up by 1 byte because
            // we can't encode an 8 byte "base" value, so instead this case uses a
            // base of 7 bytes plus a one-byte offset in each of the 2 blocks.
            int expected_size = 23 + (level + 2) / 4;
            if (level == 30) expected_size += 1;
            var points = new S2Point[32].Fill(id.ToPoint());
            TestEncodedS2PointVector(points, CodingHint.COMPACT, expected_size);
        }
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_SnappedFractalLoops()
    {
        S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
#if DEBUG
        const int kMaxPoints = 3 << 10;
#else
            const int kMaxPoints = 3 << 14;
#endif

        for (int num_points = 3; num_points <= kMaxPoints; num_points *= 4)
        {
            uint s2polygon_size = 0, lax_polygon_size = 0;
            for (int i = 0; i < 10; ++i)
            {
                S2Testing.Fractal fractal = new();
                fractal.SetLevelForApproxMaxEdges(num_points);
                var frame = S2Testing.GetRandomFrame();
                var loop = fractal.MakeLoop(frame, S2Testing.KmToAngle(10));
                List<S2Point> points = new();
                for (int j = 0; j < loop.NumVertices; ++j)
                {
                    points.Add(new S2CellId(loop.Vertex(j)).ToPoint());
                }
                S2Polygon s2polygon = new(new S2Loop(points));
                Encoder encoder = new();
                s2polygon.Encode(encoder);
                s2polygon_size += (uint)encoder.Length();
                // S2LaxPolygonShape has 2 extra bytes of overhead to encode one loop.
                lax_polygon_size +=
                    (uint)TestEncodedS2PointVector(points.ToArray(), CodingHint.COMPACT, -1) + 2;
            }
            _logger.WriteLine($"n={num_points:d5}  s2={s2polygon_size:d9}  lax={lax_polygon_size:d9}");
        }
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_RoundtripEncodingFast()
    {
        TestRoundtripEncoding(CodingHint.FAST);
    }

    [Fact]
    internal void Test_EncodedS2PointVectorTest_RoundtripEncodingCompact()
    {
        TestRoundtripEncoding(CodingHint.COMPACT);
    }

    private static int TestEncodedS2PointVector(S2Point[] expected, CodingHint hint, Int64 expected_bytes)
    {
        Encoder encoder = new();
        EncodedS2PointVector.EncodeS2PointVector(expected, hint, encoder);
        if (expected_bytes >= 0)
        {
            Assert.Equal(expected_bytes, encoder.Length());
        }

        //// Allocate storage aligned to double.  Since there is a varint at the
        //// beginning of the encoding, the following S2Points will be unaligned,
        //// so we test unaligned S2Point access for UNCOMPRESSED encodings.
        //auto aligned = make_unique<double[]>(encoder.length() / sizeof(double) + 1);
        //std::memcpy(aligned.get(), encoder.base(), encoder.length());

        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedS2PointVector.Init(decoder);
        Assert.True(success);
        Assert.Equal(actual!.Decode(), expected);
        return encoder.Length();
    }

    // In order to make it easier to construct tests that encode particular
    // values, this function duplicates the part of EncodedS2PointVector that
    // converts an encoded 64-bit value back to an S2Point.
    private static S2Point EncodedValueToPoint(UInt64 value, int level)
    {
        BitsInterleave.DeinterleaveUInt32(value, out var sj, out var tj);
        int shift = S2.kMaxCellLevel - level;
        uint si = (((sj << 1) | 1) << shift) & 0x7fffffff;
        uint ti = (((tj << 1) | 1) << shift) & 0x7fffffff;
        int face = (int)(((sj << shift) >> 30) | (((tj << (shift + 1)) >> 29) & 4));
        return S2.FaceUVtoXYZ(face,
            S2.STtoUV(S2.SiTitoST(si)),
            S2.STtoUV(S2.SiTitoST(ti))).Normalize();
    }

    private static void TestRoundtripEncoding(CodingHint hint)
    {
        // Ensures that the EncodedS2PointVector can be encoded and decoded without
        // loss.
        const int level = 3;
        var pointsArray = new S2Point[kBlockSize];
        pointsArray.Fill(EncodedValueToPoint(0, level));
        var pointsList = pointsArray.ToList();
        pointsList.Add(EncodedValueToPoint(0x78, level));
        pointsList.Add(EncodedValueToPoint(0x7a, level));
        pointsList.Add(EncodedValueToPoint(0x7c, level));
        pointsList.Add(EncodedValueToPoint(0x84, level));
        pointsArray = pointsList.ToArray();

        EncodedS2PointVector a_vector = new();
        Encoder a_encoder = new();

        EncodedS2PointVector b_vector = new();
        Encoder b_encoder = new();

        // Encode and decode from a vector<S2Point>.
        {
            EncodedS2PointVector.EncodeS2PointVector(pointsArray, hint, a_encoder);
            var decoder = a_encoder.GetDecoder();
            var (success, a_vector_) = EncodedS2PointVector.Init(decoder);
            Assert.True(success);
            a_vector = a_vector_!;
        }
        Assert.Equal(pointsList, a_vector.Decode());

        // Encode and decode from an EncodedS2PointVector.
        {
            a_vector.Encode(b_encoder);
            var decoder = b_encoder.GetDecoder();
            var (success, b_vector_) = EncodedS2PointVector.Init(decoder);
            Assert.True(success);
            b_vector = b_vector_!;
        }
        Assert.Equal(pointsList, b_vector.Decode());
    }
}
