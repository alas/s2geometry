using Xunit;

namespace S2Geometry
{
    public class EncodedStringVectorTests
    {
        [Fact]
        public void Test_EncodedStringVectorTest_Empty() {
            TestEncodedStringVector(System.Array.Empty<string>(), 1);
        }

        [Fact]
        public void Test_EncodedStringVectorTest_EmptyString() {
            TestEncodedStringVector(new[] { "" }, 2);
        }

        [Fact]
        public void Test_EncodedStringVectorTest_RepeatedEmptyStrings() {
            TestEncodedStringVector(new[] { "", "", "" }, 4);
        }

        [Fact]
        public void Test_EncodedStringVectorTest_OneString() {
            TestEncodedStringVector(new[] { "apples" }, 8);
        }

        [Fact]
        public void Test_EncodedStringVectorTest_TwoStrings() {
            TestEncodedStringVector(new[] { "fuji", "mutsu" }, 12);
        }

        [Fact]
        public void Test_EncodedStringVectorTest_TwoBigStrings() {
            TestEncodedStringVector(new[] { new string('x', 10000), new string('y', 100000) }, 110007);
        }

        private static void TestEncodedStringVector(string[] input, int expected_bytes)
        {
            Encoder encoder = new();
            StringVectorEncoder.Encode(input, encoder);
            Assert.Equal(expected_bytes, encoder.Length);
            Decoder decoder = new(encoder.Buffer, 0, encoder.Length);
            Assert.True(EncodedStringVector.Init(decoder, out var actual));
            Assert.Equal(actual.Decode(), input);
        }
    }
}
