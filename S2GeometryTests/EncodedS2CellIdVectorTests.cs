namespace S2Geometry;

public class EncodedS2CellIdVectorTests
{
    [Fact]
    internal void Test_EncodedS2CellIdVector_Empty()
    {
        TestEncodedS2CellIdVector(new List<S2CellId> { }, 2);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_None()
    {
        TestEncodedS2CellIdVector(new List<S2CellId> { S2CellId.None }, 3);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_NoneNone()
    {
        TestEncodedS2CellIdVector(new List<S2CellId> { S2CellId.None, S2CellId.None }, 4);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_Sentinel()
    {
        TestEncodedS2CellIdVector(new List<S2CellId> { S2CellId.Sentinel }, 10);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_MaximumShiftCell()
    {
        // Tests the encoding of a single cell at level 2, which corresponds the
        // maximum encodable shift value (56).
        TestEncodedS2CellIdVector(new List<S2CellId> { MakeCellIdOrDie("0/00") }, 3);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_SentinelSentinel()
    {
        TestEncodedS2CellIdVector(new List<S2CellId> { S2CellId.Sentinel, S2CellId.Sentinel }, 11);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_NoneSentinelNone()
    {
        TestEncodedS2CellIdVector(new List<S2CellId>
            { S2CellId.None, S2CellId.Sentinel, S2CellId.None}, 26);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_InvalidCells()
    {
        // Tests that cells with an invalid LSB can be encoded.
        TestEncodedS2CellIdVector(new List<UInt64> { 0x6, 0xe, 0x7e }, 5);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_OneByteLeafCells()
    {
        // Tests that (1) if all cells are leaf cells, the low bit is not encoded,
        // and (2) this can be indicated using the standard 1-byte header.
        TestEncodedS2CellIdVector(new List<UInt64> { 0x3, 0x7, 0x177 }, 5);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_OneByteLevel29Cells()
    {
        // Tests that (1) if all cells are at level 29, the low bit is not encoded,
        // and (2) this can be indicated using the standard 1-byte header.
        TestEncodedS2CellIdVector(new List<UInt64> { 0xc, 0x1c, 0x47c }, 5);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_OneByteLevel28Cells()
    {
        // Tests that (1) if all cells are at level 28, the low bit is not encoded,
        // and (2) this can be indicated using the extended 2-byte header.
        TestEncodedS2CellIdVector(new List<UInt64> { 0x30, 0x70, 0x1770 }, 6);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_OneByteMixedCellLevels()
    {
        // Tests that cells at mixed levels can be encoded in one byte.
        TestEncodedS2CellIdVector(new List<UInt64> { 0x300, 0x1c00, 0x7000, 0xff00 }, 6);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_OneByteMixedCellLevelsWithPrefix()
    {
        // Tests that cells at mixed levels can be encoded in one byte even when
        // they share a multi-byte prefix.
        TestEncodedS2CellIdVector(new List<UInt64>{
                0x1234567800000300, 0x1234567800001c00,
                0x1234567800007000, 0x123456780000ff00}, 10);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_OneByteRangeWithBaseValue()
    {
        // Tests that cells can be encoded in one byte by choosing a base value
        // whose bit range overlaps the delta values.
        // 1 byte header, 3 bytes base, 1 byte size, 4 bytes deltas
        TestEncodedS2CellIdVector(new List<UInt64>{
                0x00ffff0000000000, 0x0100fc0000000000,
                0x0100500000000000, 0x0100330000000000}, 9);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_MaxShiftRange()
    {
        byte[] bytes = {
            (31 << 3)  // 31 -> add 29 to bytes[1].
            + 1,       // Number of encoded cell IDs.
            27,        // 27+29 is the maximum supported shift.
            1, 0       // Encoded cell ID. Not important.
        };
        Decoder decoder = new(bytes, 0, bytes.Length);
        var (success, _) = EncodedS2CellIdVector.Init(decoder);
        Assert.True(success);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_ShiftOutOfRange()
    {
        byte[] bytes = {
            (31 << 3)  // 31 -> add 29 to bytes[1].
            + 1,       // Number of encoded cell IDs.
            28,        // 28+29 is greater than the maximum supported shift of 56.
            1, 0       // Encoded cell ID. Not important.
        };
        Decoder decoder = new(bytes, 0, bytes.Length);
        var (success, _) = EncodedS2CellIdVector.Init(decoder);
        Assert.False(success);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_SixFaceCells()
    {
        List<S2CellId> ids = new();
        for (int face = 0; face < 6; ++face)
        {
            ids.Add(S2CellId.FromFace(face));
        }
        TestEncodedS2CellIdVector(ids, 8);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_FourLevel10Children()
    {
        List<S2CellId> ids = new();
        var parent = MakeCellIdOrDie("3/012301230");
        for (var id = parent.ChildBegin();
             id != parent.ChildEnd(); id = id.Next())
        {
            ids.Add(id);
        }
        TestEncodedS2CellIdVector(ids, 8);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_FractalS2ShapeIndexCells()
    {
        S2Testing.Fractal fractal = new();
        fractal.SetLevelForApproxMaxEdges(3 * 1024);
        S2Point center = MakePointOrDie("47.677:-122.206");
        MutableS2ShapeIndex index = new();
        index.Add(new S2Loop.Shape(
            fractal.MakeLoop(S2.GetFrame(center), S1Angle.FromDegrees(1))));
        List<S2CellId> ids = new();
        foreach (var shape in index.GetNewEnumerable())
        {
            ids.Add(shape.Item1);
        }
        Assert.Equal(966, ids.Count);
        TestEncodedS2CellIdVector(ids, 2902);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_CoveringCells()
    {
        List<UInt64> ids = new()
        {
            0x414a617f00000000, 0x414a61c000000000, 0x414a624000000000,
            0x414a63c000000000, 0x414a647000000000, 0x414a64c000000000,
            0x414a653000000000, 0x414a704000000000, 0x414a70c000000000,
            0x414a714000000000, 0x414a71b000000000, 0x414a7a7c00000000,
            0x414a7ac000000000, 0x414a8a4000000000, 0x414a8bc000000000,
            0x414a8c4000000000, 0x414a8d7000000000, 0x414a8dc000000000,
            0x414a914000000000, 0x414a91c000000000, 0x414a924000000000,
            0x414a942c00000000, 0x414a95c000000000, 0x414a96c000000000,
            0x414ab0c000000000, 0x414ab14000000000, 0x414ab34000000000,
            0x414ab3c000000000, 0x414ab44000000000, 0x414ab4c000000000,
            0x414ab6c000000000, 0x414ab74000000000, 0x414ab8c000000000,
            0x414ab94000000000, 0x414aba1000000000, 0x414aba3000000000,
            0x414abbc000000000, 0x414abe4000000000, 0x414abec000000000,
            0x414abf4000000000, 0x46b5454000000000, 0x46b545c000000000,
            0x46b5464000000000, 0x46b547c000000000, 0x46b5487000000000,
            0x46b548c000000000, 0x46b5494000000000, 0x46b54a5400000000,
            0x46b54ac000000000, 0x46b54b4000000000, 0x46b54bc000000000,
            0x46b54c7000000000, 0x46b54c8004000000, 0x46b54ec000000000,
            0x46b55ad400000000, 0x46b55b4000000000, 0x46b55bc000000000,
            0x46b55c4000000000, 0x46b55c8100000000, 0x46b55dc000000000,
            0x46b55e4000000000, 0x46b5604000000000, 0x46b560c000000000,
            0x46b561c000000000, 0x46ca424000000000, 0x46ca42c000000000,
            0x46ca43c000000000, 0x46ca444000000000, 0x46ca45c000000000,
            0x46ca467000000000, 0x46ca469000000000, 0x46ca5fc000000000,
            0x46ca604000000000, 0x46ca60c000000000, 0x46ca674000000000,
            0x46ca679000000000, 0x46ca67f000000000, 0x46ca684000000000,
            0x46ca855000000000, 0x46ca8c4000000000, 0x46ca8cc000000000,
            0x46ca8e5400000000, 0x46ca8ec000000000, 0x46ca8f0100000000,
            0x46ca8fc000000000, 0x46ca900400000000, 0x46ca98c000000000,
            0x46ca994000000000, 0x46ca99c000000000, 0x46ca9a4000000000,
            0x46ca9ac000000000, 0x46ca9bd500000000, 0x46ca9e4000000000,
            0x46ca9ec000000000, 0x46caf34000000000, 0x46caf4c000000000,
            0x46caf54000000000
        };
        Assert.Equal(97, ids.Count);
        TestEncodedS2CellIdVector(ids, 488);
    }

    [Fact]
    internal void Test_EncodedS2CellIdVector_LowerBoundLimits()
    {
        // Test seeking before the beginning and past the end of the vector.
        var first = S2CellId.Begin(S2.kMaxCellLevel);
        var last = S2CellId.End(S2.kMaxCellLevel).Prev();
        Encoder encoder = new();
        var cell_ids = MakeEncodedS2CellIdVector(
            new() { first, last }, encoder);
        Assert.Equal(0, cell_ids.LowerBound(S2CellId.None));
        Assert.Equal(0, cell_ids.LowerBound(first));
        Assert.Equal(1, cell_ids.LowerBound(first.Next()));
        Assert.Equal(1, cell_ids.LowerBound(last.Prev()));
        Assert.Equal(1, cell_ids.LowerBound(last));
        Assert.Equal(2, cell_ids.LowerBound(last.Next()));
        Assert.Equal(2, cell_ids.LowerBound(S2CellId.Sentinel));
    }

    // Encodes the given vector and returns the corresponding
    // EncodedS2CellIdVector (which points into the Encoder's data buffer).
    private static EncodedS2CellIdVector MakeEncodedS2CellIdVector(List<S2CellId> input, Encoder encoder)
    {
        EncodedS2CellIdVector.EncodeS2CellIdVector(input, encoder);
        var decoder = encoder.GetDecoder();
        var (success, cell_ids) = EncodedS2CellIdVector.Init(decoder);
        Assert.True(success);
        return cell_ids!;
    }

    // Encodes the given vector and checks that it has the expected size and
    // contents.
    private static void TestEncodedS2CellIdVector(List<S2CellId> expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedS2CellIdVector actual = MakeEncodedS2CellIdVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var dec = actual.Decode();
        Assert.True(dec.SequenceEqual(expected));
    }

    // Like the above, but accepts a UInt64[] rather than a S2CellId[].
    private static void TestEncodedS2CellIdVector(List<UInt64> raw_expected, int expected_bytes)
    {
        List<S2CellId> expected = new();
        foreach (UInt64 raw_id in raw_expected)
        {
            expected.Add(new(raw_id));
        }
        TestEncodedS2CellIdVector(expected, expected_bytes);
    }
}
