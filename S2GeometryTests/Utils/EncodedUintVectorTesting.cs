namespace S2Geometry;

internal class EncodedUIntVectorTesting
{
    // TestEncodedUIntVector
    internal static void TestEncodedUIntVector_16(ushort[] expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedUIntVector<ushort>.EncodeUIntVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector<ushort>.Init(decoder);
        Assert.True(success);
        Assert.Equal(actual!.Decode(), expected);
    }

    // TestLowerBound
    internal static void TestLowerBound_16(int bytes_per_value, int num_values)
    {
        var v = MakeSortedTestVector_16(bytes_per_value, num_values);
        Encoder encoder = new();
        var actual = MakeEncodedVector_16(v, encoder);
        foreach (var x in v)
        {
            Assert.Equal(v.GetLowerBound(x), actual.LowerBound(x));
            if (x > 0)
            {
                Assert.Equal(v.GetLowerBound(x - 1), actual.LowerBound((ushort)(x - 1)));
            }
        }
    }

    private static ushort[] MakeSortedTestVector_16(int bytes_per_value, int num_values)
    {
        Assert.True(bytes_per_value <= sizeof(ushort));
        ushort limit_value = (ushort)(~0 >> (8 * (sizeof(ushort) - bytes_per_value)));
        List<ushort> values = [];
        for (int i = 0; i + 1 < num_values; ++i)
        {
            values.Add((ushort)(limit_value * ((double)i / (num_values - 1))));
        }
        // The last value needs special handling since casting it to "double" loses
        // precision when T == UInt64.
        values.Add(limit_value);
        Assert.True(values.IsSorted());
        return [.. values];
    }

    internal static EncodedUIntVector<ushort> MakeEncodedVector_16(ushort[] values, Encoder encoder)
    {
        EncodedUIntVector<ushort>.EncodeUIntVector(values, encoder);
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector<ushort>.Init(decoder);
        Assert.True(success);
        return actual!;
    }

    // TestEncodedUIntVector
    internal static void TestEncodedUIntVector_32(uint[] expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedUIntVector<uint>.EncodeUIntVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector<uint>.Init(decoder);
        Assert.True(success);
        Assert.Equal(actual!.Decode(), expected);
    }

    // TestLowerBound
    internal static void TestLowerBound_32(int bytes_per_value, int num_values)
    {
        var v = MakeSortedTestVector_32(bytes_per_value, num_values);
        Encoder encoder = new();
        var actual = MakeEncodedVector_32(v, encoder);
        foreach (var x in v)
        {
            Assert.Equal(v.GetLowerBound(x), actual.LowerBound(x));
            if (x > 0)
            {
                Assert.Equal(v.GetLowerBound(x - 1), actual.LowerBound((uint)(x - 1)));
            }
        }
    }

    private static uint[] MakeSortedTestVector_32(int bytes_per_value, int num_values)
    {
        Assert.True(bytes_per_value <= sizeof(uint));
        uint limit_value = (uint)(~0 >> (8 * (sizeof(uint) - bytes_per_value)));
        List<uint> values = [];
        for (int i = 0; i + 1 < num_values; ++i)
        {
            values.Add((uint)(limit_value * ((double)i / (num_values - 1))));
        }
        // The last value needs special handling since casting it to "double" loses
        // precision when T == UInt64.
        values.Add(limit_value);
        Assert.True(values.IsSorted());
        return [.. values];
    }

    internal static EncodedUIntVector<uint> MakeEncodedVector_32(uint[] values, Encoder encoder)
    {
        EncodedUIntVector<uint>.EncodeUIntVector(values, encoder);
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector<uint>.Init(decoder);
        Assert.True(success);
        return actual!;
    }

    // TestEncodedUIntVector
    internal static void TestEncodedUIntVector_64(ulong[] expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedUIntVector<ulong>.EncodeUIntVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector<ulong>.Init(decoder);
        Assert.True(success);
        Assert.Equal(actual!.Decode(), expected);
    }

    // TestLowerBound
    internal static void TestLowerBound_64(int bytes_per_value, int num_values)
    {
        var v = MakeSortedTestVector_64(bytes_per_value, num_values);
        Encoder encoder = new();
        var actual = MakeEncodedVector_64(v, encoder);
        foreach (var x in v)
        {
            Assert.Equal(v.GetLowerBound(x), actual.LowerBound(x));
            if (x > 0)
            {
                Assert.Equal(v.GetLowerBound(x - 1), actual.LowerBound((ulong)(x - 1)));
            }
        }
    }

    private static ulong[] MakeSortedTestVector_64(int bytes_per_value, int num_values)
    {
        Assert.True(bytes_per_value <= sizeof(ulong));
        ulong limit_value = (ulong)(~0 >> (8 * (sizeof(ulong) - bytes_per_value)));
        List<ulong> values = [];
        for (int i = 0; i + 1 < num_values; ++i)
        {
            values.Add((ulong)(limit_value * ((double)i / (num_values - 1))));
        }
        // The last value needs special handling since casting it to "double" loses
        // precision when T == UInt64.
        values.Add(limit_value);
        Assert.True(values.IsSorted());
        return [.. values];
    }

    internal static EncodedUIntVector<ulong> MakeEncodedVector_64(ulong[] values, Encoder encoder)
    {
        EncodedUIntVector<ulong>.EncodeUIntVector(values, encoder);
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector<ulong>.Init(decoder);
        Assert.True(success);
        return actual!;
    }
}
