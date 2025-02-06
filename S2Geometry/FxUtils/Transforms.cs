// Data transforms that can help code more efficiently.
//
// ZigZag Transform
//
// Good for varint coding small signed integers centered around 0.
//
//       Int32 .     UInt32
// -------------------------
//           0 .          0
//          -1 .          1
//           1 .          2
//          -2 .          3
//         ... .        ...
//  2147483647 . 4294967294
// -2147483648 . 4294967295
//
//        >> encode >>
//        << decode <<

namespace S2Geometry;

public static class Transforms
{
    public static uint ZigZagEncode(Int32 n)
    {
        // We need the cast to avoid an arithmetic shift.
        uint sign = ((uint)n) >> 31;
        return ((uint)n << 1) ^ (0u - sign);
    }

    public static Int32 ZigZagDecode(uint n)
    {
        return (Int32)((n >> 1) ^ (0u - (n & 1)));
    }

    public static ulong ZigZagEncode64(long n)
    {
        // We need the cast to avoid an arithmetic shift.
        ulong sign = ((ulong)n) >> 63;
        return ((ulong)n << 1) ^ (0u - sign);
    }

    public static long ZigZagDecode64(ulong n)
    {
        return (long)((n >> 1) ^ (0u - (n & 1)));
    }
}
