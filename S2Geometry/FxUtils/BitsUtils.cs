namespace S2Geometry;

public static class BitsUtils
{
    public static bool IsBitSet(uint num, int pos) => (num & (1 << pos)) != 0;

    public static int SetBit(int num, int pos) => num | (1 << pos);

    #region FindMSBSetNonZero

    public static int FindMSBSetNonZero64(UInt64 n) => Log2FloorNonZero64(n);

    public static int FindMSBSetNonZero(UInt32 n) => Log2FloorNonZero(n);

    private static int Log2FloorNonZero64(UInt64 n)
    {
        var topbits = (UInt32)(n >> 32);
        if (topbits == 0)
        {
            // Top bits are zero, so scan in bottom bits
            return Log2FloorNonZero((UInt32)n);
        }
        else
        {
            return 32 + Log2FloorNonZero(topbits);
        }
    }

    public static int Log2Floor64(UInt64 n)
    {
        var topbits = (UInt32)(n >> 32);
        if (topbits == 0)
        {
            // Top bits are zero, so scan in bottom bits
            return Log2Floor((UInt32)n);
        }
        else
        {
            return 32 + Log2FloorNonZero(topbits);
        }
    }

    public static int Log2FloorNonZero(UInt32 n) => Log2Floor(n);

    private static int Log2Floor(UInt32 n)
    {
        if (n == 0) return -1;

        var log = 0;
        var value = n;
        for (var i = 4; i >= 0; --i)
        {
            var shift = (1 << i);
            var x = value >> shift;
            if (x != 0)
            {
                value = x;
                log += shift;
            }
        }
        //Assert.True(value == 1);
        return log;
    }

    #endregion

    #region FindLSBSetNonZero

    public static int FindLSBSetNonZero64(UInt64 n)
    {
        var bottombits = (UInt32)n;
        if (bottombits == 0)
        {
            // Bottom bits are zero, so scan in top bits
            return 32 + FindLSBSetNonZero((UInt32)(n >> 32));
        }
        else
        {
            return FindLSBSetNonZero(bottombits);
        }
    }

    public static int FindLSBSetNonZero(uint n)
    {
        var rc = 31;
        var shift = 1 << 4;
        for (var i = 4; i >= 0; --i)
        {
            var x = n << shift;
            if (x != 0)
            {
                n = x;
                rc -= shift;
            }
            shift >>= 1;
        }
        return rc;
    }

    #endregion
}
