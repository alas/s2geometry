using System;

namespace S2Geometry
{
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
	public static class Transforms
	{
		public static UInt32 ZigZagEncode(Int32 n)
		{
			// We need the cast to avoid an arithmetic shift.
			UInt32 sign = ((UInt32)n) >> 31;
			return ((UInt32)n << 1) ^ (0u - sign);
		}

		public static Int32 ZigZagDecode(UInt32 n)
		{
			return (Int32)((n >> 1) ^ (0u - (n & 1)));
		}

		public static UInt64 ZigZagEncode64(Int64 n)
		{
			// We need the cast to avoid an arithmetic shift.
			UInt64 sign = ((UInt64)n) >> 63;
			return ((UInt64)n << 1) ^ (0u - sign);
		}

		public static Int64 ZigZagDecode64(UInt64 n)
		{
			return (Int64)((n >> 1) ^ (0u - (n & 1)));
		}
	}
}
