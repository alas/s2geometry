using System;
using System.Runtime.InteropServices;
using Xunit;

namespace S2Geometry
{
    public class EncodedUintVectorTests
    {
        static EncodedUintVectorTests()
        {
            Assert.Equal(16, Marshal.SizeOf(typeof(EncodedUintVector_UInt64)));
        }

        [Fact]
        public void Test_EncodedUintVectorTest_Empty() {
            TestEncodedUintVector(Array.Empty<uint>(), 1);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_Zero() {
            TestEncodedUintVector(new UInt64[] { 0 }, 2);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_RepeatedZeros() {
            TestEncodedUintVector(new UInt16[] { 0, 0, 0}, 4);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_MaxInt() {
            TestEncodedUintVector(new UInt64[]{ ~0UL}, 9);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_OneByte() {
            TestEncodedUintVector(new UInt64[]{ 0, 255, 1, 254}, 5);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_TwoBytes() {
            TestEncodedUintVector(new UInt64[]{ 0, 255, 256, 254}, 9);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_ThreeBytes() {
            TestEncodedUintVector(new UInt64[]{ 0xffffff, 0x0102, 0, 0x050403}, 13);
        }

        [Fact]
        public void Test_EncodedUintVectorTest_EightBytes() {
            TestEncodedUintVector(new UInt64[]{ ~0UL, 0, 0x0102030405060708}, 25);
        }

        [Fact]
        public void Test_EncodedUintVector_LowerBound()
        {
            for (int bytes_per_value = 8; bytes_per_value <= 8; ++bytes_per_value)
            {
                EncodedUintVectorTesting.TestLowerBound_64(bytes_per_value, 10);
                if (bytes_per_value <= 4)
                {
                    EncodedUintVectorTesting.TestLowerBound_32(bytes_per_value, 500);
                    if (bytes_per_value <= 2)
                    {
                        EncodedUintVectorTesting.TestLowerBound_16(bytes_per_value, 100);
                    }
                }
            }
        }

        private static void TestEncodedUintVector<T>(T[] expected, int expected_bytes)
        {
            switch (expected)
            {
                case ushort[] _ushort:
                    EncodedUintVectorTesting.TestEncodedUintVector_16(_ushort, expected_bytes);
                    break;
                case uint[] _uint:
                    EncodedUintVectorTesting.TestEncodedUintVector_32(_uint, expected_bytes);
                    break;
                case ulong[] _ulong:
                    EncodedUintVectorTesting.TestEncodedUintVector_64(_ulong, expected_bytes);
                    break;
                default: throw new NotSupportedException(typeof(T).FullName);
            }
        }
    }
}