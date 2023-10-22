namespace S2Geometry;

public class EncodedStringVectorTests
{
    [Fact]
    internal void Test_EncodedStringVectorTest_Empty() =>
        TestEncodedStringVector([], 1);

    [Fact]
    internal void Test_EncodedStringVectorTest_EmptyString() =>
        TestEncodedStringVector([""], 2);

    [Fact]
    internal void Test_EncodedStringVectorTest_RepeatedEmptyStrings() =>
        TestEncodedStringVector(["", "", ""], 4);

    [Fact]
    internal void Test_EncodedStringVectorTest_OneString() =>
        TestEncodedStringVector(["apples"], 8);

    [Fact]
    internal void Test_EncodedStringVectorTest_TwoStrings() =>
        TestEncodedStringVector(["fuji", "mutsu"], 12);

    [Fact]
    internal void Test_EncodedStringVectorTest_TwoBigStrings() =>
        TestEncodedStringVector([new string('x', 10000), new string('y', 100000)], 110007);

    private static void TestEncodedStringVector(string[] input, int expected_bytes)
    {
        Encoder encoder = new();
        StringVectorEncoder.Encode(input, encoder);
        Assert.Equal(expected_bytes, encoder.Length());
        var decoder = encoder.GetDecoder();
        var (success, actual) = EncodedStringVector.Init(decoder);
        Assert.True(success);
        Assert.Equal(actual!.Decode(), input);

        // Check that `EncodedStringVector::Encode` produces the same result as
        // `StringVectorEncoder::Encode`, as documented.
        Encoder reencoder=new();
        actual.Encode(reencoder);
        Assert.True(encoder.Equals(reencoder));
    }
}
