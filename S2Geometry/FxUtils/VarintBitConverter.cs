using System;

namespace S2Geometry
{
    /// <summary>
    /// from: https://github.com/topas/VarintBitConverter/blob/master/src/VarintBitConverter/VarintBitConverter.cs
    /// </summary>
    public class VarintBitConverter
    {
        /// <summary>
        /// Returns the specified byte value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">Byte value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(byte value, byte[] arr, ref int offset)
        {
            GetVarintBytes((ulong)value, arr, ref offset);
        }

        /// <summary>
        /// Returns the specified 16-bit signed value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">16-bit signed value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(short value, byte[] arr, ref int offset)
        {
            var zigzag = EncodeZigZag(value, 16);
            GetVarintBytes((ulong)zigzag, arr, ref offset);
        }

        /// <summary>
        /// Returns the specified 16-bit unsigned value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">16-bit unsigned value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(ushort value, byte[] arr, ref int offset)
        {
            GetVarintBytes((ulong)value, arr, ref offset);
        }

        /// <summary>
        /// Returns the specified 32-bit signed value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">32-bit signed value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(int value, byte[] arr, ref int offset)
        {
            var zigzag = EncodeZigZag(value, 32);
            GetVarintBytes((ulong)zigzag, arr, ref offset);
        }

        /// <summary>
        /// Returns the specified 32-bit unsigned value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">32-bit unsigned value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(uint value, byte[] arr, ref int offset)
        {
            GetVarintBytes((ulong)value, arr, ref offset);
        }

        /// <summary>
        /// Returns the specified 64-bit signed value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">64-bit signed value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(long value, byte[] arr, ref int offset)
        {
            var zigzag = EncodeZigZag(value, 64);
            GetVarintBytes((ulong)zigzag, arr, ref offset);
        }

        /// <summary>
        /// Returns the specified 64-bit unsigned value as varint encoded array of bytes.   
        /// </summary>
        /// <param name="value">64-bit unsigned value</param>
        /// <returns>Varint array of bytes.</returns>
        public static void GetVarintBytes(ulong value, byte[] arr, ref int offset)
        {
            do
            {
                var byteVal = value & 0x7f;
                value >>= 7;

                if (value != 0)
                {
                    byteVal |= 0x80;
                }

                arr[offset++] = (byte)byteVal;

            } while (value != 0);
        }

        /// <summary>
        /// Returns byte value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>Byte value</returns>
        public static (byte result, int offset) ToByte(byte[] bytes, int offset)
        {
            var (result, newoffset) = ToTarget(bytes, 8, offset);
            return ((byte)result, newoffset);
        }

        /// <summary>
        /// Returns 16-bit signed value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>16-bit signed value</returns>
        public static (short result, int offset) ToInt16(byte[] bytes, int offset)
        {
            var (zigzag, newoffset) = ToTarget(bytes, 16, offset);
            return ((short)DecodeZigZag(zigzag), newoffset);
        }

        /// <summary>
        /// Returns 16-bit usigned value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>16-bit usigned value</returns>
        public static (ushort result, int offset) ToUInt16(byte[] bytes, int offset)
        {
            var (result, newoffset) = ToTarget(bytes, 16, offset);
            return ((ushort)result, newoffset);
        }

        /// <summary>
        /// Returns 32-bit signed value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>32-bit signed value</returns>
        public static (int result, int offset) ToInt32(byte[] bytes, int offset)
        {
            var (zigzag, newoffset) = ToTarget(bytes, 32, offset);
            return ((int)DecodeZigZag(zigzag), newoffset);
        }

        /// <summary>
        /// Returns 32-bit unsigned value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>32-bit unsigned value</returns>
        public static (uint result, int offset) ToUInt32(byte[] bytes, int offset)
        {
            var (result, newoffset) = ToTarget(bytes, 32, offset);
            return ((uint)result, newoffset);
        }

        /// <summary>
        /// Returns 64-bit signed value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>64-bit signed value</returns>
        public static (long result, int offset) ToInt64(byte[] bytes, int offset)
        {
            var (zigzag, newoffset) = ToTarget(bytes, 64, offset);
            return (DecodeZigZag(zigzag), newoffset);
        }

        /// <summary>
        /// Returns 64-bit unsigned value from varint encoded array of bytes.
        /// </summary>
        /// <param name="bytes">Varint encoded array of bytes.</param>
        /// <returns>64-bit unsigned value</returns>
        public static (ulong result, int count) ToUInt64(byte[] bytes, int offset)
        {
            return ToTarget(bytes, 64, offset);
        }

        private static long EncodeZigZag(long value, int bitLength)
        {
            return (value << 1) ^ (value >> (bitLength - 1));
        }

        private static long DecodeZigZag(ulong value)
        {
            if ((value & 0x1) == 0x1)
            {
                return (-1 * ((long)(value >> 1) + 1));
            }

            return (long)(value >> 1);
        }

        private static (ulong result, int offset) ToTarget(byte[] bytes, int sizeBites, int offset)
        {
            int shift = 0;
            ulong result = 0;
            ulong byteValue;

            do
            {
                byteValue = bytes[offset++];
                ulong tmp = byteValue & 0x7f;
                result |= tmp << shift;

                if (shift > sizeBites)
                {
                    throw new ArgumentOutOfRangeException(nameof(bytes), "Byte array is too large.");
                }

                shift += 7;

            } while ((byteValue & 0x80) == 0x80);
            
            return (result, offset);
        }
    }
}
