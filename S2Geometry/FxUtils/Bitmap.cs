namespace S2Geometry;

public class Bitmap64
{
    private const int kIntBits = 8 * sizeof(ulong);

    private int KLogIntBits { get; } = (int)Log2(kIntBits, 0);

    // A value of type `Word` with all bits set to one. If `Word` !=
    // `ArithmeticWord`, i.e. `Word` gets promoted, then this is a `Word` with all
    // bits set to one then promoted to `ArithmeticWord`. In other words, the
    // lowest N bits are one, and all other bits are zero, where N is the width of
    // `Word`.
    //
    // For example, for uint32, this is 0xFFFFFFFF. However, if `Word` is
    // `unit8_t`, and `ArithmeticWord` is equal to `uint32` and kAllOnesWord is
    // 0x000000FF.
    private const ulong kAllOnesWord = ~0UL;


    private readonly ulong[] Bits;

    private int Length { get; }

    public Bitmap64(int length, bool initialValue = false)
    {
        this.Length = length;
        this.Bits = new ulong[(length / kIntBits) + 1];

        if (initialValue)
        {
            this.SetAllBits(initialValue);
        }
    }

    private void SetAllBits(bool value)
    {
        var bucketVal = value ? ulong.MaxValue : ulong.MinValue;
        for (int index = 0; index < this.Bits.Length; index++)
        {
            this.Bits[index] = bucketVal;
        }
    }

    public bool Get(int index)
    {
        if (!this.IsValidIndex(index))
        {
            throw new ArgumentOutOfRangeException(nameof(index));
        }

        return GetInternal(index);
    }

    private bool GetInternal(int index) =>
        (this.Bits[index / kIntBits] & (0x1ul << (index % kIntBits))) > 0;

    public void Set(int index, bool value)
    {
        if (!this.IsValidIndex(index))
        {
            throw new ArgumentOutOfRangeException(nameof(index));
        }

        SetInternal(index, value);
    }

    private void SetInternal(int index, bool value)
    {
        if (value)
        {
            this.Bits[index / kIntBits] |= 0x1ul << (index % kIntBits);
        }
        else
        {
            this.Bits[index / kIntBits] &= ~(0x1ul << (index % kIntBits));
        }
    }

    public void Toggle(int index)
    {
        if (!this.IsValidIndex(index))
        {
            throw new ArgumentOutOfRangeException(nameof(index));
        }

        SetInternal(index, !GetInternal(index));
    }

    public int GetOnesCount()
    {
        int trueBitsCount = 0;

        for (int bucketIndex = 0; bucketIndex < this.Bits.Length; bucketIndex++)
        {
            var bucket = this.Bits[bucketIndex];
            int bitIndex = 0;
            while (bucket > 0 && ((bucketIndex * kIntBits) + bitIndex) < this.Length)
            {
                if ((bucket & 0x1) > 0)
                {
                    trueBitsCount++;
                }

                bucket >>= 1;
                bitIndex++;
            }
        }

        return trueBitsCount;
    }

    // If no bits in this bitmap are set, returns false. Otherwise returns true
    // and puts the index of the first set bit in this bitmap in *index. Note
    // that *index is modified even if we return false.
    public bool FindFirstSetBit(out int index) 
    {
        var indexTmp = 0;
        var result = FindNextSetBit(ref indexTmp);
        index = indexTmp;
        return result;
    }

    public bool FindNextSetBit(ref int index) 
    {
        return FindNextSetBitBeforeLimit(ref index, Length);
    }

    private bool FindNextSetBitBeforeLimit(ref int index, int limit)
    {
        if (!this.IsValidIndex((int)index))
        {
            throw new ArgumentOutOfRangeException(nameof(index));
        }

        int index_as_size_t = index;
        bool result = FindNextSetBitInVector(this.Bits, ref index_as_size_t, limit);
        if (result)
        {
            index = index_as_size_t;
        }
        return result;
    }

    // A static version of FindNextSetBitBeforeLimit that can be called
    // by other clients that have an array of words in their hands,
    // layed out in the same way as BitMap.  Scans bits in "*words"
    // starting at bit "*bit_index", looking for a set bit.  If it finds
    // a set bit before reaching bit index "bit_limit", sets
    // "*bit_index" to the bit index and returns true.  Otherwise
    // returns false.  Will not dereference "words" past
    // "words[(bit_limit+31)/kIntBits]".
    private bool FindNextSetBitInVector(ulong[] words, ref int bit_index, int bit_limit)
    {
        return FindNextBitInVector(/*complement=*/false, words, ref bit_index, bit_limit);
    }

    private bool FindNextBitInVector(bool complement, ulong[] words, ref int bit_index_inout, int limit)
    {
        var bit_index = bit_index_inout;
        if (bit_index >= limit) return false;

        // From now on limit != 0, since if it was we would have returned false.
        int int_index = bit_index >> KLogIntBits;
        var one_word = words[int_index];
        if (complement) one_word = ~one_word;

        // Simple optimization where we can immediately return true if the first
        // bit is set.  This helps for cases where many bits are set, and doesn't
        // hurt too much if not.
        var first_bit_offset = bit_index & (kIntBits - 1);
        if ((one_word & (ulong)(1 << first_bit_offset)) != 0) return true;

        // First word is special - we need to mask off leading bits
        one_word &= kAllOnesWord << first_bit_offset;

        // Loop through all but the last word.  Note that 'limit' is one
        // past the last bit we want to check, and we don't want to read
        // past the end of "words".  E.g. if size_ == kIntBits only words[0] is
        // valid, so we want to avoid reading words[1] when limit == kIntBits.
        int last_int_index = (limit - 1) >> KLogIntBits;
        while (int_index<last_int_index)
        {
            if (one_word != 0)
            {
                bit_index_inout =
                    (int_index << KLogIntBits) + BitsUtils.FindLSBSetNonZero64(one_word);
                return true;
            }
            one_word = words[++int_index];
            if (complement) one_word = ~one_word;
        }

        // Last word is special - we may need to mask off trailing bits.  Note that
        // 'limit' is one past the last bit we want to check, and if limit is a
        // multiple of kIntBits we want to check all bits in this word.
        one_word &= ~((kAllOnesWord - 1) << ((limit - 1) & (kIntBits - 1)));
        if (one_word != 0)
        {
            bit_index_inout =
                (int_index << KLogIntBits) + BitsUtils.FindLSBSetNonZero64(one_word);
            return true;
        }
        return false;
    }

    private static UInt32 Log2(UInt32 n, UInt32 p = 0)
    {
        return (n <= 1) ? p : Log2(n / 2, p + 1);
    }

    private bool IsValidIndex(int index)
    {
        return index >= 0 && index < this.Length;
    }
}
