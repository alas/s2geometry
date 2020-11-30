using System;
namespace S2Geometry
{
    /// <summary>
    /// This class represents an encoded vector of unsigned integers of type T.
    /// Values are decoded only when they are accessed.  This allows for very fast
    /// initialization and no additional memory use beyond the encoded data.
    /// The encoded data is not owned by this class; typically it points into a
    /// large contiguous buffer that contains other encoded data as well.
    ///
    /// This is one of several helper classes that allow complex data structures to
    /// be initialized from an encoded format in constant time and then decoded on
    /// demand.  This can be a big performance advantage when only a small part of
    /// the data structure is actually used.
    ///
    /// Values are encoded using a fixed number of bytes per value, where the
    /// number of bytes depends on the largest value present.
    ///
    /// REQUIRES: T is an unsigned integer type.
    /// REQUIRES: 2 <= sizeof(T) <= 8
    /// </summary>
    public class EncodedUintVector
    {
        static EncodedUintVector()
        {
#if BIGENDIAN || ARM
            throw new NotImplementedException("Not implemented on big-endian architectures");
#endif

            // UInt16 checks
            Assert.True(UInt16.MinValue == 0); // Unsupported signed integer
            Assert.True((sizeof(UInt16) & 0xE) != 0); // Unsupported integer length
            Assert.True(2 <= sizeof(UInt16)); // Unsupported integer length
            Assert.True(sizeof(UInt16) <= 8); // Unsupported integer length

            // UInt32 checks
            Assert.True(UInt32.MinValue == 0); // Unsupported signed integer
            Assert.True((sizeof(UInt32) & 0xE) != 0); // Unsupported integer length
            Assert.True(2 <= sizeof(UInt32)); // Unsupported integer length
            Assert.True(sizeof(UInt32) <= 8); // Unsupported integer length

            // UInt64 checks
            Assert.True(UInt64.MinValue == 0); // Unsupported signed integer
            Assert.True((sizeof(UInt64) & 0xE) != 0); // Unsupported integer length
            Assert.True(2 <= sizeof(UInt64)); // Unsupported integer length
            Assert.True(sizeof(UInt64) <= 8); // Unsupported integer length
        }

        // Encodes a vector of unsigned integers in a format that can later be
        // decoded as an EncodedUintVector.
        //
        // REQUIRES: T is an unsigned integer type.
        // REQUIRES: 2 <= sizeof(T) <= 8
        // REQUIRES: "encoder" uses the default constructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public static void EncodeUintVector(UInt16[] v, Encoder encoder)
        {
            // The encoding is as follows:
            //
            //   varint64: (v.size() * sizeof(T)) | (len - 1)
            //   array of v.size() elements ["len" bytes each]
            //
            // Note that we don't allow (len == 0) since this would require an extra bit
            // to encode the length.

            UInt16 one_bits = 1;  // Ensures len >= 1.
            foreach (var x in v)
                one_bits |= x;
            int len = (BitsUtils.FindMSBSetNonZero64(one_bits) >> 3) + 1;
            Assert.True(len >= 1 && len <= 8);

            // Note that the multiplication is optimized into a bit shift.
            encoder.Ensure(Encoder.kVarintMax64 + v.Length * len);
            UInt64 size_len = ((UInt64)(UInt32)v.Length * sizeof(UInt16)) | (UInt32)(len - 1);
            encoder.PutVarInt64(size_len);
            foreach (var x in v)
            {
                EncodeUintWithLength(x, len, encoder);
            }
        }

        /// <summary>
        /// Encodes an unsigned integer in little-endian format using "length" bytes.
        /// (The client must ensure that the encoder's buffer is large enough.)
        ///
        /// REQUIRES: T is an unsigned integer type.
        /// REQUIRES: 2 <= sizeof(T) <= 8
        /// REQUIRES: 0 <= length <= sizeof(T)
        /// REQUIRES: value < 256 ** length
        /// REQUIRES: encoder.Avail() >= length
        /// </summary>
        public static void EncodeUintWithLength(UInt16 value, int length, Encoder encoder)
        {
            Assert.True(length >= 0 && length <= sizeof(UInt16));
            Assert.True(encoder.Avail() >= length);

            while (--length >= 0)
            {
                encoder.Put8((byte)value);
                value = (UInt16)(value >> 8);
            }
            Assert.True(value == 0);
        }

        // Decodes and consumes a variable-length integer consisting of "length" bytes
        // in little-endian format.  Returns false if not enough bytes are available.
        //
        // REQUIRES: T is an unsigned integer type.
        // REQUIRES: 2 <= sizeof(T) <= 8
        // REQUIRES: 0 <= length <= sizeof(T)
        public static bool DecodeUintWithLength(int length, Decoder decoder, out UInt16 result)
        {
            result = default;
            if (decoder.Avail() < length) return false;
            result = EncodedUintVector_UInt16.GetUintWithLength(decoder.buf_, decoder.offset_, length);
            decoder.Skip(length);
            return true;
        }

        // Encodes a vector of unsigned integers in a format that can later be
        // decoded as an EncodedUintVector.
        //
        // REQUIRES: T is an unsigned integer type.
        // REQUIRES: 2 <= sizeof(T) <= 8
        // REQUIRES: "encoder" uses the default constructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public static void EncodeUintVector(UInt32[] v, Encoder encoder)
        {
            // The encoding is as follows:
            //
            //   varint64: (v.size() * sizeof(T)) | (len - 1)
            //   array of v.size() elements ["len" bytes each]
            //
            // Note that we don't allow (len == 0) since this would require an extra bit
            // to encode the length.

            UInt32 one_bits = 1;  // Ensures len >= 1.
            foreach (var x in v)
                one_bits |= x;
            int len = (BitsUtils.FindMSBSetNonZero64(one_bits) >> 3) + 1;
            Assert.True(len >= 1 && len <= 8);

            // Note that the multiplication is optimized into a bit shift.
            encoder.Ensure(Encoder.kVarintMax64 + v.Length * len);
            UInt64 size_len = ((UInt64)(UInt32)v.Length * sizeof(UInt32)) | (UInt32)(len - 1);
            encoder.PutVarInt64(size_len);
            foreach (var x in v)
            {
                EncodeUintWithLength(x, len, encoder);
            }
        }

        /// <summary>
        /// Encodes an unsigned integer in little-endian format using "length" bytes.
        /// (The client must ensure that the encoder's buffer is large enough.)
        ///
        /// REQUIRES: T is an unsigned integer type.
        /// REQUIRES: 2 <= sizeof(T) <= 8
        /// REQUIRES: 0 <= length <= sizeof(T)
        /// REQUIRES: value < 256 ** length
        /// REQUIRES: encoder.Avail() >= length
        /// </summary>
        public static void EncodeUintWithLength(UInt32 value, int length, Encoder encoder)
        {
            Assert.True(length >= 0 && length <= sizeof(UInt32));
            Assert.True(encoder.Avail() >= length);

            while (--length >= 0)
            {
                encoder.Put8((byte)value);
                value = (UInt32)(value >> 8);
            }
            Assert.True(value == 0);
        }

        // Decodes and consumes a variable-length integer consisting of "length" bytes
        // in little-endian format.  Returns false if not enough bytes are available.
        //
        // REQUIRES: T is an unsigned integer type.
        // REQUIRES: 2 <= sizeof(T) <= 8
        // REQUIRES: 0 <= length <= sizeof(T)
        public static bool DecodeUintWithLength(int length, Decoder decoder, out UInt32 result)
        {
            result = default;
            if (decoder.Avail() < length) return false;
            result = EncodedUintVector_UInt32.GetUintWithLength(decoder.buf_, decoder.offset_, length);
            decoder.Skip(length);
            return true;
        }

        // Encodes a vector of unsigned integers in a format that can later be
        // decoded as an EncodedUintVector.
        //
        // REQUIRES: T is an unsigned integer type.
        // REQUIRES: 2 <= sizeof(T) <= 8
        // REQUIRES: "encoder" uses the default constructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public static void EncodeUintVector(UInt64[] v, Encoder encoder)
        {
            // The encoding is as follows:
            //
            //   varint64: (v.size() * sizeof(T)) | (len - 1)
            //   array of v.size() elements ["len" bytes each]
            //
            // Note that we don't allow (len == 0) since this would require an extra bit
            // to encode the length.

            UInt64 one_bits = 1;  // Ensures len >= 1.
            foreach (var x in v)
                one_bits |= x;
            int len = (BitsUtils.FindMSBSetNonZero64(one_bits) >> 3) + 1;
            Assert.True(len >= 1 && len <= 8);

            // Note that the multiplication is optimized into a bit shift.
            encoder.Ensure(Encoder.kVarintMax64 + v.Length * len);
            UInt64 size_len = ((UInt64)(UInt32)v.Length * sizeof(UInt64)) | (UInt32)(len - 1);
            encoder.PutVarInt64(size_len);
            foreach (var x in v)
            {
                EncodeUintWithLength(x, len, encoder);
            }
        }

        /// <summary>
        /// Encodes an unsigned integer in little-endian format using "length" bytes.
        /// (The client must ensure that the encoder's buffer is large enough.)
        ///
        /// REQUIRES: T is an unsigned integer type.
        /// REQUIRES: 2 <= sizeof(T) <= 8
        /// REQUIRES: 0 <= length <= sizeof(T)
        /// REQUIRES: value < 256 ** length
        /// REQUIRES: encoder.Avail() >= length
        /// </summary>
        public static void EncodeUintWithLength(UInt64 value, int length, Encoder encoder)
        {
            Assert.True(length >= 0 && length <= sizeof(UInt64));
            Assert.True(encoder.Avail() >= length);

            while (--length >= 0)
            {
                encoder.Put8((byte)value);
                value = (UInt64)(value >> 8);
            }
            Assert.True(value == 0);
        }

        // Decodes and consumes a variable-length integer consisting of "length" bytes
        // in little-endian format.  Returns false if not enough bytes are available.
        //
        // REQUIRES: T is an unsigned integer type.
        // REQUIRES: 2 <= sizeof(T) <= 8
        // REQUIRES: 0 <= length <= sizeof(T)
        public static bool DecodeUintWithLength(int length, Decoder decoder, out UInt64 result)
        {
            result = default;
            if (decoder.Avail() < length) return false;
            result = EncodedUintVector_UInt64.GetUintWithLength(decoder.buf_, decoder.offset_, length);
            decoder.Skip(length);
            return true;
        }
    }

    public class EncodedUintVector_UInt16
    {
        // Returns the size of the original vector.
        public int Count { get; private set; }

        /// <summary>
        /// Decodes a variable-length integer consisting of "length" bytes starting at
        /// "ptr" in little-endian format.
        ///
        /// REQUIRES: 0 <= length <= sizeof(T)
        /// </summary>
        public static UInt16 GetUintWithLength(byte[] ptr, int index, int length)
        {
            Assert.True(length >= 0 && length <= sizeof(UInt16));

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

            if ((length & sizeof(UInt16)) != 0)
            {
                return BitConverter.ToUInt16(ptr, index);
            }
            UInt16 x = 0;
            index += length;
             if ((length & 1) != 0)
            {
                x = (UInt16)(((UInt16)(x << 8)) + ptr[--index]);
            }
            return x;
        }

        // Initializes the EncodedUintVector.  Returns false on errors, leaving the
        // vector in an unspecified state.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder)
        {
            if (!decoder.TryGetVar64(out var size_len)) return false;
            Count = (int)(size_len / sizeof(UInt16));  // Optimized into bit shift.
            len_ = (Byte)((size_len & ((UInt16)(sizeof(UInt16) - 1))) + 1);
            if (Count > int.MaxValue / sizeof(UInt16)) return false;
            int bytes = Count * len_;
            if (decoder.Avail() < bytes) return false;
            data_ = decoder.buf_;
            decoder.Skip(bytes);
            return true;
        }

        // Resets the vector to be empty.
        public void Clear()
        {
            Count = 0;
            data_ = null;
        }

        // Returns the element at the given index.  
        public UInt16 this[int i]
        {
            get
            {
                Assert.True(i >= 0 && i < Count);
                return GetUintWithLength(data_, i * len_, len_);
            }
        }

        // Returns the index of the first element x such that (x >= target), or
        // size() if no such element exists.
        //
        // REQUIRES: The vector elements are sorted in non-decreasing order.
        public int LowerBound(UInt16 target)
        {
            Assert.True(len_ >= 1 && len_ <= sizeof(UInt16));

            // TODO(ericv): Consider using the unused 28 bits of "len_" to store the
            // last result of GetLowerBound() to be used as a hint.  This should help in
            // common situation where the same element is looked up repeatedly.  This
            // would require declaring the new field (length_GetLowerBound_hint_) as
            // mutable std.atomic<uint32> (accessed using std.memory_order_relaxed)
            // with a custom copy constructor that resets the hint component to zero.
            return len_ switch
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

        private int LowerBound(UInt16 target, int length)
        {
            var lo = 0;
            var hi = Count;
            while (lo < hi)
            {
                var mid = (lo + hi) >> 1;
                UInt16 value = GetUintWithLength(data_, mid * length, length);
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

        // Decodes and returns the entire original vector.
        public UInt16[] Decode()
        {
            var result = new UInt16[Count];
            for (int i = 0; i < Count; i++)
            {
                result[i] = this[i];
            }
            return result;
        }

        private byte[] data_;
        private byte len_;
    }

    public class EncodedUintVector_UInt32
    {
        // Returns the size of the original vector.
        public int Count { get; private set; }

        /// <summary>
        /// Decodes a variable-length integer consisting of "length" bytes starting at
        /// "ptr" in little-endian format.
        ///
        /// REQUIRES: 0 <= length <= sizeof(T)
        /// </summary>
        public static UInt32 GetUintWithLength(byte[] ptr, int index, int length)
        {
            Assert.True(length >= 0 && length <= sizeof(UInt32));

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

            if ((length & sizeof(UInt32)) != 0)
            {
                return BitConverter.ToUInt32(ptr, index);
            }
            UInt32 x = 0;
            index += length;
             if ((length & 2) != 0)
            {
                index -= sizeof(UInt16);
                x = (x << 16) + BitConverter.ToUInt16(ptr, index);
            }
            if ((length & 1) != 0)
            {
                x = (UInt32)(((UInt32)(x << 8)) + ptr[--index]);
            }
            return x;
        }

        // Initializes the EncodedUintVector.  Returns false on errors, leaving the
        // vector in an unspecified state.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder)
        {
            if (!decoder.TryGetVar64(out var size_len)) return false;
            Count = (int)(size_len / sizeof(UInt32));  // Optimized into bit shift.
            len_ = (Byte)((size_len & ((UInt32)(sizeof(UInt32) - 1))) + 1);
            if (Count > int.MaxValue / sizeof(UInt32)) return false;
            int bytes = Count * len_;
            if (decoder.Avail() < bytes) return false;
            data_ = decoder.buf_;
            decoder.Skip(bytes);
            return true;
        }

        // Resets the vector to be empty.
        public void Clear()
        {
            Count = 0;
            data_ = null;
        }

        // Returns the element at the given index.  
        public UInt32 this[int i]
        {
            get
            {
                Assert.True(i >= 0 && i < Count);
                return GetUintWithLength(data_, i * len_, len_);
            }
        }

        // Returns the index of the first element x such that (x >= target), or
        // size() if no such element exists.
        //
        // REQUIRES: The vector elements are sorted in non-decreasing order.
        public int LowerBound(UInt32 target)
        {
            Assert.True(len_ >= 1 && len_ <= sizeof(UInt32));

            // TODO(ericv): Consider using the unused 28 bits of "len_" to store the
            // last result of GetLowerBound() to be used as a hint.  This should help in
            // common situation where the same element is looked up repeatedly.  This
            // would require declaring the new field (length_GetLowerBound_hint_) as
            // mutable std.atomic<uint32> (accessed using std.memory_order_relaxed)
            // with a custom copy constructor that resets the hint component to zero.
            return len_ switch
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

        private int LowerBound(UInt32 target, int length)
        {
            var lo = 0;
            var hi = Count;
            while (lo < hi)
            {
                var mid = (lo + hi) >> 1;
                UInt32 value = GetUintWithLength(data_, mid * length, length);
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

        // Decodes and returns the entire original vector.
        public UInt32[] Decode()
        {
            var result = new UInt32[Count];
            for (int i = 0; i < Count; i++)
            {
                result[i] = this[i];
            }
            return result;
        }

        private byte[] data_;
        private byte len_;
    }

    public class EncodedUintVector_UInt64
    {
        // Returns the size of the original vector.
        public int Count { get; private set; }

        /// <summary>
        /// Decodes a variable-length integer consisting of "length" bytes starting at
        /// "ptr" in little-endian format.
        ///
        /// REQUIRES: 0 <= length <= sizeof(T)
        /// </summary>
        public static UInt64 GetUintWithLength(byte[] ptr, int index, int length)
        {
            Assert.True(length >= 0 && length <= sizeof(UInt64));

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

            if ((length & sizeof(UInt64)) != 0)
            {
                return BitConverter.ToUInt64(ptr, index);
            }
            UInt64 x = 0;
            index += length;
             if ((length & 4) != 0)
            {
                index -= sizeof(UInt32);
                x = BitConverter.ToUInt32(ptr, index);
            }
            if ((length & 2) != 0)
            {
                index -= sizeof(UInt16);
                x = (x << 16) + BitConverter.ToUInt16(ptr, index);
            }
            if ((length & 1) != 0)
            {
                x = (UInt64)(((UInt64)(x << 8)) + ptr[--index]);
            }
            return x;
        }

        // Initializes the EncodedUintVector.  Returns false on errors, leaving the
        // vector in an unspecified state.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder)
        {
            if (!decoder.TryGetVar64(out var size_len)) return false;
            Count = (int)(size_len / sizeof(UInt64));  // Optimized into bit shift.
            len_ = (Byte)((size_len & ((UInt64)(sizeof(UInt64) - 1))) + 1);
            if (Count > int.MaxValue / sizeof(UInt64)) return false;
            int bytes = Count * len_;
            if (decoder.Avail() < bytes) return false;
            data_ = decoder.buf_;
            decoder.Skip(bytes);
            return true;
        }

        // Resets the vector to be empty.
        public void Clear()
        {
            Count = 0;
            data_ = null;
        }

        // Returns the element at the given index.  
        public UInt64 this[int i]
        {
            get
            {
                Assert.True(i >= 0 && i < Count);
                return GetUintWithLength(data_, i * len_, len_);
            }
        }

        // Returns the index of the first element x such that (x >= target), or
        // size() if no such element exists.
        //
        // REQUIRES: The vector elements are sorted in non-decreasing order.
        public int LowerBound(UInt64 target)
        {
            Assert.True(len_ >= 1 && len_ <= sizeof(UInt64));

            // TODO(ericv): Consider using the unused 28 bits of "len_" to store the
            // last result of GetLowerBound() to be used as a hint.  This should help in
            // common situation where the same element is looked up repeatedly.  This
            // would require declaring the new field (length_GetLowerBound_hint_) as
            // mutable std.atomic<uint32> (accessed using std.memory_order_relaxed)
            // with a custom copy constructor that resets the hint component to zero.
            return len_ switch
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

        private int LowerBound(UInt64 target, int length)
        {
            var lo = 0;
            var hi = Count;
            while (lo < hi)
            {
                var mid = (lo + hi) >> 1;
                UInt64 value = GetUintWithLength(data_, mid * length, length);
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

        // Decodes and returns the entire original vector.
        public UInt64[] Decode()
        {
            var result = new UInt64[Count];
            for (int i = 0; i < Count; i++)
            {
                result[i] = this[i];
            }
            return result;
        }

        private byte[] data_;
        private byte len_;
    }
}
