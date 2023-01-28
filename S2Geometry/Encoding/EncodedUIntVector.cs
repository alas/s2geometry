// This class represents an encoded vector of unsigned integers of type T.
// Values are decoded only when they are accessed.  This allows for very fast
// initialization and no additional memory use beyond the encoded data.
// The encoded data is not owned by this class; typically it points into a
// large contiguous buffer that contains other encoded data as well.
//
// This is one of several helper classes that allow complex data structures to
// be initialized from an encoded format in constant time and then decoded on
// demand.  This can be a big performance advantage when only a small part of
// the data structure is actually used.
//
// Values are encoded using a fixed number of bytes per value, where the
// number of bytes depends on the largest value present.
//
// REQUIRES: T is an unsigned integer type.
// REQUIRES: 2 <= sizeof(T) <= 8

namespace S2Geometry;

using System;
using System.Numerics;

public class EncodedUIntVector<T> where T : IUnsignedNumber<T>, IBinaryInteger<T>, IMinMaxValue<T>, IBitwiseOperators<T, T, T>//, IShiftOperators<T, T, T>
{
    // Returns the size of the original vector.
    public int Count { get; }

    private byte[] Data { get; }

    private int Offset { get; }

    private byte Len { get; }

    static EncodedUIntVector()
    {
        MyDebug.Assert(typeof(T) is IUnsignedNumber<T>, "Unsupported signed integer");
        MyDebug.Assert((typeof(T).SizeOf() & 0xe) != 0, "Unsupported integer length");
    }

    /// <summary>
    /// Constructs an uninitialized object; requires Init() to be called.
    /// 
    /// Note(Alas): added the parameters to this constructor, Init does not need to be called
    /// </summary>
    public EncodedUIntVector(int count, byte[] data, int offset, byte len)
    {
        Count = count;
        Data = data;
        Offset = offset;
        Len = len;
    }

    /// <summary>
    /// Encodes a vector of unsigned integers in a format that can later be
    /// decoded as an EncodedUintVector.
    ///
    /// REQUIRES: T is an unsigned integer type.
    /// REQUIRES: 2 <= sizeof(T) <= 8
    /// REQUIRES: "encoder" uses the default constructor, so that its buffer
    ///           can be enlarged as necessary by calling Ensure(int).
    /// </summary>
    public static void EncodeUIntVector(ReadOnlySpan<T> v, Encoder encoder)
    {
        // The encoding is as follows:
        //
        //   varint64: (v.size() * sizeof(T)) | (len - 1)
        //   array of v.size() elements ["len" bytes each]
        //
        // Note that we don't allow (len == 0) since this would require an extra bit
        // to encode the length.

        T one_bits = T.One;  // Ensures len >= 1.
        foreach (var x in v) one_bits |= x;
        var one_bits_ulong = (ulong)Convert.ChangeType(one_bits, typeof(ulong));
        var len = (BitsUtils.FindMSBSetNonZero64(one_bits_ulong) >> 3) + 1;
        MyDebug.Assert(len >= 1 && len <= 8);

        // Note that the multiplication is optimized into a bit shift.
        encoder.Ensure(Encoder.kVarintMax64 + v.Length * len);
        var size_len = (v.Length * typeof(T).SizeOf()) | (len - 1);
        encoder.PutVarInt64(size_len);
        foreach (var x in v)
        {
            EncodeUIntWithLength(x, len, encoder);
        }
    }

    /// <summary>
    /// Initializes the EncodedUintVector.  Returns false on errors, leaving the
    /// vector in an unspecified state.
    ///
    /// REQUIRES: The Decoder data buffer must outlive this object.
    /// </summary>
    public static (bool success, EncodedUIntVector<T>? vector) Init(Decoder decoder)
    {

        if (!decoder.TryGetVarUInt64(out var size_len)) return (false, null);
        var size = typeof(T).SizeOf();
        var count = (int)(size_len / (ulong)size);  // Optimized into bit shift.
        var len = (byte)((size_len & ((ulong)(size - 1))) + 1);
        if (count > int.MaxValue / size) return (false, null);
        int bytes = count * len;
        if (decoder.Avail() < bytes) return (false, null);
        var data = decoder.Buffer;
        var offset = decoder.Offset;
        decoder.Skip(bytes);
        return (true, new(count, data, offset, len));
    }

    // Resets the vector to be empty.
    // 
    // Note(Alas): removed to make fields readonly, don't reuse, create a new instance
    //public void Clear()
    //{
    //    Count = 0;
    //    data_ = null;
    //}

    /// <summary>
    /// Encodes an unsigned integer in little-endian format using "length" bytes.
    /// (The client must ensure that the encoder's buffer is large enough.)
    ///
    /// REQUIRES: T is an unsigned integer type.
    /// REQUIRES: 2 <= sizeof(T) <= 8
    /// REQUIRES: 0 <= length <= sizeof(T)
    /// REQUIRES: value < 256 ** length
    /// REQUIRES: encoder->avail() >= length
    /// </summary>
    public static void EncodeUIntWithLength(T value, int length, Encoder encoder)
    {
        //MyDebug.Assert(value is IUnsignedNumber<T>, "Unsupported signed integer");
        //MyDebug.Assert((typeof(T).SizeOf() & 0xe) != 0, "Unsupported integer length");
        MyDebug.Assert(length >= 0 && length <= typeof(T).SizeOf());
        MyDebug.Assert(encoder.Avail() >= length);

        var bytes = new byte[length];
        value.WriteLittleEndian(bytes, 0);
        foreach (var b in bytes)
        {
            encoder.Put8(b);
        }
    }

    /// <summary>
    /// Decodes and consumes a variable-length integer consisting of "length" bytes
    /// in little-endian format.  Returns false if not enough bytes are available.
    ///
    /// REQUIRES: T is an unsigned integer type.
    /// REQUIRES: 2 <= sizeof(T) <= 8
    /// REQUIRES: 0 <= length <= sizeof(T)
    /// </summary>
    public static bool DecodeUIntWithLength(int length, Decoder decoder, out T result)
    {
        result = T.Zero;
        if (decoder.Avail() < length) return false;
        result = GetUIntWithLength(decoder.Buffer, decoder.Offset, length);
        decoder.Skip(length);
        return true;
    }

    /// <summary>
    /// Decodes a variable-length integer consisting of "length" bytes starting at
    /// "ptr" in little-endian format.
    ///
    /// REQUIRES: T is an unsigned integer type.
    /// REQUIRES: 2 <= sizeof(T) <= 8
    /// REQUIRES: 0 <= length <= sizeof(T)
    /// </summary>
    public static T GetUIntWithLength(byte[] ptr, int index, int length)
    {
        var size = typeof(T).SizeOf();
        //MyDebug.Assert(value is IUnsignedNumber<T>, "Unsupported signed integer");
        //MyDebug.Assert((size & 0xe) != 0, "Unsupported integer length");
        MyDebug.Assert(length >= 0 && length <= size);

        // Note that the following code is faster than any of the following:
        //
        //  - A loop that repeatedly loads and shifts one byte.
        //  - memcpying "length" bytes to a local variable of type T.
        //  - A switch statement that handles each length optimally.
        //
        // The following code is slightly faster:
        //
        //   T mask = (length == 0) ? 0 : ~T{0} >> 8 * (sizeof(T) - length);
        //   return *reinterpret_cast<const T*>(ptr) & mask;
        //
        // However this technique is unsafe because in extremely rare cases it might
        // access out-of-bounds heap memory.  (This can only happen if "ptr" is
        // within (sizeof(T) - length) bytes of the end of a memory page and the
        // following page in the address space is unmapped.)

        if ((length & size) != 0)
        {
            if (size == 8) return (T)Convert.ChangeType(BitConverter.ToUInt64(ptr, index), typeof(T));
            if (size == 4) return (T)Convert.ChangeType(BitConverter.ToUInt32(ptr, index), typeof(T));
            if (size == 2) return (T)Convert.ChangeType(BitConverter.ToUInt16(ptr, index), typeof(T));
            MyDebug.Assert(size == 1);
            return (T)Convert.ChangeType(ptr[index], typeof(T));
        }
        T x = T.Zero;
        index += length;
        if (size > 4 && (length & 4) != 0)
        {
            index -= sizeof(UInt32);
            x = (T)Convert.ChangeType(BitConverter.ToUInt32(ptr, index), typeof(T));
        }
        if (size > 2 && (length & 2) != 0)
        {
            index -= sizeof(UInt16);
            x = (x << 16) + (T)Convert.ChangeType(BitConverter.ToUInt16(ptr, index), typeof(T));
        }
        if (size > 1 && (length & 1) != 0)
        {
            x = (x << 8) + (T)Convert.ChangeType(ptr[--index], typeof(T));
        }
        return x;
    }

    /// <summary>
    /// Returns the element at the given index.  
    /// </summary>
    public T this[int i]
    {
        get
        {
            MyDebug.Assert(i >= 0 && i < Count);
            return GetUIntWithLength(Data, Offset + i * Len, Len);
        }
    }

    /// <summary>
    /// Returns the index of the first element x such that (x >= target), or
    /// size() if no such element exists.
    ///
    /// REQUIRES: The vector elements are sorted in non-decreasing order.
    /// </summary>
    public int LowerBound(T target)
    {
        MyDebug.Assert((typeof(T).SizeOf() & 0xe) != 0, "Unsupported integer length");
        MyDebug.Assert(Len >= 1 && Len <= typeof(T).SizeOf());

        // TODO(ericv): Consider using the unused 28 bits of "len_" to store the
        // last result of lower_bound() to be used as a hint.  This should help in
        // common situation where the same element is looked up repeatedly.  This
        // would require declaring the new field (length_lower_bound_hint_) as
        // mutable std::atomic<uint32> (accessed using std::memory_order_relaxed)
        // with a custom copy constructor that resets the hint component to zero.
        return Len switch
        {
            1 => LowerBound(target, 1),
            2 => LowerBound(target, 2),
            3 => LowerBound(target, 3),
            4 => LowerBound(target, 4),
            5 => LowerBound(target, 5),
            6 => LowerBound(target, 6),
            7 => LowerBound(target, 7),
            _ => LowerBound(target, 8),
        };
    }

    private int LowerBound(T target, int length)
    {
        var lo = 0;
        var hi = Count;
        while (lo < hi)
        {
            var mid = (lo + hi) >> 1;
            T value = GetUIntWithLength(Data, Offset + mid * length, length);
            if (value.CompareTo(target) < 0)
            {
                lo = mid + 1;
            }
            else
            {
                hi = mid;
            }
        }
        return lo;
    }

    /// <summary>
    /// Decodes and returns the entire original vector.
    /// </summary>
    public T[] Decode()
    {
        var result = new T[Count];
        for (int i = 0; i < Count; i++)
        {
            result[i] = this[i];
        }
        return result;
    }

    /// <summary>
    /// The encoding must be identical to StringVectorEncoder::Encode().
    /// </summary>
    public void Encode(Encoder encoder)
    {
        var size_len = (Count * typeof(T).SizeOf()) | (Len - 1);

        encoder.Ensure(Encoder.kVarintMax64 + size_len);
        encoder.PutVarInt64(size_len);
        encoder.PutN(Data, Offset, Count * Len);
    }
}
