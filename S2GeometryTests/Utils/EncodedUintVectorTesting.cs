namespace S2Geometry;

internal class EncodedUIntVectorTesting
{
    // TestEncodedUIntVector
    internal static void TestEncodedUIntVector_16(UInt16[] expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedUIntVector2<UInt16>.EncodeUIntVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector2<UInt16>.Init(decoder);
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
                Assert.Equal(v.GetLowerBound(x - 1), actual.LowerBound((UInt16)(x - 1)));
            }
        }
    }

    private static UInt16[] MakeSortedTestVector_16(int bytes_per_value, int num_values)
    {
        Assert.True(bytes_per_value <= sizeof(UInt16));
        UInt16 limit_value = (UInt16)(~0 >> (8 * (sizeof(UInt16) - bytes_per_value)));
        List<UInt16> values = new();
        for (int i = 0; i + 1 < num_values; ++i)
        {
            values.Add((UInt16)(limit_value * ((double)i / (num_values - 1))));
        }
        // The last value needs special handling since casting it to "double" loses
        // precision when T == UInt64.
        values.Add(limit_value);
        Assert.True(values.IsSorted());
        return values.ToArray();
    }

    internal static EncodedUIntVector2<UInt16> MakeEncodedVector_16(UInt16[] values, Encoder encoder)
    {
        EncodedUIntVector2<UInt16>.EncodeUIntVector(values, encoder);
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector2<UInt16>.Init(decoder);
        Assert.True(success);
        return actual!;
    }

    // TestEncodedUIntVector
    internal static void TestEncodedUIntVector_32(UInt32[] expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedUIntVector2<UInt32>.EncodeUIntVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector2<UInt32>.Init(decoder);
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
                Assert.Equal(v.GetLowerBound(x - 1), actual.LowerBound((UInt32)(x - 1)));
            }
        }
    }

    private static UInt32[] MakeSortedTestVector_32(int bytes_per_value, int num_values)
    {
        Assert.True(bytes_per_value <= sizeof(UInt32));
        UInt32 limit_value = (UInt32)(~0 >> (8 * (sizeof(UInt32) - bytes_per_value)));
        List<UInt32> values = new();
        for (int i = 0; i + 1 < num_values; ++i)
        {
            values.Add((UInt32)(limit_value * ((double)i / (num_values - 1))));
        }
        // The last value needs special handling since casting it to "double" loses
        // precision when T == UInt64.
        values.Add(limit_value);
        Assert.True(values.IsSorted());
        return values.ToArray();
    }

    internal static EncodedUIntVector2<UInt32> MakeEncodedVector_32(UInt32[] values, Encoder encoder)
    {
        EncodedUIntVector2<UInt32>.EncodeUIntVector(values, encoder);
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector2<UInt32>.Init(decoder);
        Assert.True(success);
        return actual!;
    }

    // TestEncodedUIntVector
    internal static void TestEncodedUIntVector_64(UInt64[] expected, int expected_bytes)
    {
        Encoder encoder = new();
        EncodedUIntVector2<UInt64>.EncodeUIntVector(expected, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector2<UInt64>.Init(decoder);
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
                Assert.Equal(v.GetLowerBound(x - 1), actual.LowerBound((UInt64)(x - 1)));
            }
        }
    }

    private static UInt64[] MakeSortedTestVector_64(int bytes_per_value, int num_values)
    {
        Assert.True(bytes_per_value <= sizeof(UInt64));
        UInt64 limit_value = (UInt64)(~0 >> (8 * (sizeof(UInt64) - bytes_per_value)));
        List<UInt64> values = new();
        for (int i = 0; i + 1 < num_values; ++i)
        {
            values.Add((UInt64)(limit_value * ((double)i / (num_values - 1))));
        }
        // The last value needs special handling since casting it to "double" loses
        // precision when T == UInt64.
        values.Add(limit_value);
        Assert.True(values.IsSorted());
        return values.ToArray();
    }

    internal static EncodedUIntVector2<UInt64> MakeEncodedVector_64(UInt64[] values, Encoder encoder)
    {
        EncodedUIntVector2<UInt64>.EncodeUIntVector(values, encoder);
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedUIntVector2<UInt64>.Init(decoder);
        Assert.True(success);
        return actual!;
    }
}
