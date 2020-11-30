using System;
using System.Collections.Generic;
namespace S2Geometry
{
    // This class represents an encoded vector of S2CellIds.  Values are decoded
    // only when they are accessed.  This allows for very fast initialization and
    // no additional memory use beyond the encoded data.  The encoded data is not
    // owned by this class; typically it points into a large contiguous buffer
    // that contains other encoded data as well.
    //
    // This is one of several helper classes that allow complex data structures to
    // be initialized from an encoded format inant time and then decoded on
    // demand.  This can be a big performance advantage when only a small part of
    // the data structure is actually used.
    //
    // The S2CellIds do not need to be sorted or at the same level.  The
    // implementation is biased towards minimizing decoding times rather than
    // space usage.
    //
    // NOTE: If your S2CellIds represent S2Points that have been snapped to
    // S2CellId centers, then EncodedS2PointVector is both faster and smaller.
    public class EncodedS2CellIdVector
    {
        // Constructs an uninitialized object; requires Init() to be called.
        public EncodedS2CellIdVector() { }

        // Encodes a vector of S2CellIds in a format that can later be decoded as an
        // EncodedS2CellIdVector.  The S2CellIds do not need to be valid.
        //
        // REQUIRES: "encoder" uses the defaultructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public static void EncodeS2CellIdVector(List<S2CellId> v, Encoder encoder)
        {
            // v[i] is encoded as (base + (deltas[i] << shift)).
            //
            // "base" consists of 0-7 bytes, and is always shifted so that its bytes are
            // the most-significant bytes of a UInt64.
            //
            // "deltas" is an EncodedUintVector<UInt64>, which means that all deltas
            // have a fixed-length encoding determined by the largest delta.
            //
            // "shift" is in the range 0..56.  The shift value is odd only if all
            // S2CellIds are at the same level, in which case the bit at position
            // (shift - 1) is automatically set to 1 in "base".
            //
            // "base" (3 bits) and "shift" (6 bits) are encoded in either one or two
            // bytes as follows:
            //
            //   - if (shift <= 4 or shift is even), then 1 byte
            //   - otherwise 2 bytes
            //
            // Note that (shift == 1) means that all S2CellIds are leaf cells, and
            // (shift == 2) means that all S2CellIds are at level 29.
            //
            // The full encoded format is as follows:
            //
            //  Byte 0, bits 0-2: base length (0-7 bytes)
            //  Byte 0, bits 3-7: encoded shift (see below)
            //  Byte 1: extended shift code (only needed for odd shifts >= 5)
            //  Followed by 0-7 bytes of "base"
            //  Followed by an EncodedUintVector of deltas.

            UInt64 v_or = 0;
            UInt64 v_and = ~0UL;
            UInt64 v_min = ~0UL;
            UInt64 v_max = 0;

            foreach (var cellid in v)
            {
                v_or |= cellid.Id;
                v_and &= cellid.Id;
                v_min = Math.Min(v_min, cellid.Id);
                v_max = Math.Max(v_max, cellid.Id);
            }

            // These variables represent the values that will used during encoding.
            UInt64 e_base = 0;        // Base value.
            int e_base_len = 0;       // Number of bytes to represent "base".
            int e_shift = 0;          // Delta shift.
            int e_max_delta_msb = 0;  // Bit position of the MSB of the largest delta.
            if (v_or > 0)
            {
                // We only allow even shifts, unless all values have the same low bit (in
                // which case the shift is odd and the preceding bit is implicitly on).
                // There is no point in allowing shifts > 56 since deltas are encoded in
                // at least 1 byte each.
                e_shift = Math.Min(56, BitsUtils.FindLSBSetNonZero64(v_or) & ~1);
                if ((v_and & (1UL << e_shift)) != 0) ++e_shift;  // All S2CellIds same level.

                // "base" consists of the "base_len" most significant bytes of the minimum
                // S2CellId.  We consider all possible values of "base_len" (0..7) and
                // choose the one that minimizes the total encoding size.
                UInt64 e_bytes = ~0UL;  // Best encoding size so far.
                for (int len = 0; len < 8; ++len)
                {
                    // "t_base" is the base value being tested (first "len" bytes of v_min).
                    // "t_max_delta_msb" is the most-significant bit position of the largest
                    // delta (or zero if there are no deltas, i.e. if v.size() == 0).
                    // "t_bytes" is the total size of the variable portion of the encoding.
                    UInt64 t_base = v_min & ~(~0UL >> (8 * len));
                    int t_max_delta_msb = Math.Max(0, BitsUtils.Log2Floor64((v_max - t_base) >> e_shift));
                    UInt64 t_bytes = (ulong)(len + v.Count * ((t_max_delta_msb >> 3) + 1));
                    if (t_bytes < e_bytes)
                    {
                        e_base = t_base;
                        e_base_len = len;
                        e_max_delta_msb = t_max_delta_msb;
                        e_bytes = t_bytes;
                    }
                }

                // It takes one extra byte to encode odd shifts (i.e., the case where all
                // S2CellIds are at the same level), so check whether we can get the same
                // encoding size per delta using an even shift.
                if (((e_shift & 1) != 0) && (e_max_delta_msb & 7) != 7) --e_shift;
            }

            Assert.True(e_base_len <= 7);
            Assert.True(e_shift <= 56);
            encoder.Ensure(2 + e_base_len);

            // As described above, "shift" and "base_len" are encoded in 1 or 2 bytes.
            // "shift_code" is 5 bits:
            //   values <= 28 represent even shifts in the range 0..56
            //   values 29, 30 represent odd shifts 1 and 3
            //   value 31 indicates that the shift is odd and encoded in the next byte
            int shift_code = e_shift >> 1;
            if ((e_shift & 1) != 0) shift_code = Math.Min(31, shift_code + 29);
            encoder.Put8((byte)((shift_code << 3) | e_base_len));
            if (shift_code == 31) {
                encoder.Put8((byte)(e_shift >> 1));  // Shift is always odd, so 3 bits unused.
            }
            // Encode the "base_len" most-significant bytes of "base".
            UInt64 base_bytes = e_base >> (64 - 8 * Math.Max(1, e_base_len));
            EncodedUintVector.EncodeUintWithLength(base_bytes, e_base_len, encoder);

            // Finally, encode the vector of deltas.
            var deltas = new UInt64[v.Count];
            for (var i = 0; i < v.Count; i++)
            {
                var cellid = v[i];
                deltas[i] = (cellid.Id - e_base) >> e_shift;
            }
            EncodedUintVector.EncodeUintVector(deltas, encoder);
        }

        // Initializes the EncodedS2CellIdVector.
        //
        // REQUIRES: The Decoder data buffer must outlive this object.
        public bool Init(Decoder decoder)
        {
            // All encodings have at least 2 bytes (one for our header and one for the
            // EncodedUintVector header), so this is safe.
            if (decoder.Avail() < 2) return false;

            // Invert the encoding of (shift_code, base_len) described above.
            int code_plus_len = decoder.Get8();
            int shift_code = code_plus_len >> 3;
            if (shift_code == 31)
            {
                shift_code = 29 + decoder.Get8();
            }
            // Decode the "base_len" most-significant bytes of "base".
            int base_len = code_plus_len & 7;
            if (!EncodedUintVector.DecodeUintWithLength(base_len, decoder, out base_)) return false;
            base_ <<= 64 - 8 * Math.Max(1, base_len);

            // Invert the encoding of "shift_code" described above.
            if (shift_code >= 29)
            {
                shift_ = (byte)(2 * (shift_code - 29) + 1);
                base_ |= 1UL << (shift_ - 1);
            }
            else
            {
                shift_ = (byte)(2 * shift_code);
            }
            return deltas_.Init(decoder);
        }
        // Returns the size of the original vector.
        public int Count => deltas_.Count;

        // Returns the element at the given index.
        public S2CellId this[int i] => new S2CellId((deltas_[i] << shift_) + base_);

        // Returns the index of the first element x such that (x >= target), or
        // size() if no such element exists.
        //
        // REQUIRES: The vector elements are sorted in non-decreasing order.
        public int LowerBound(S2CellId target)
        {
            // We optimize the search by converting "target" to the corresponding delta
            // value and then searching directly in the deltas_ vector.
            //
            // To invert operator[], we essentially compute ((target - base_) >> shift_)
            // except that we also need to round up when shifting.  The first two cases
            // ensure that "target" doesn't wrap around past zero when we do this.
            if (target.Id <= base_) return 0;
            if (target >= S2CellId.End(S2Constants.kMaxCellLevel)) return Count;
            return deltas_.LowerBound((target.Id - base_ + (1UL << shift_) - 1) >> shift_);
        }

        // Decodes and returns the entire original vector.
        public S2CellId[] Decode()
        {
            var result = new S2CellId[Count];
            for (int i = 0; i < Count; ++i)
            {
                result[i] = this[i];
            }
            return result;
        }

        // Values are decoded as (base_ + (deltas_[i] << shift_)).
        private readonly EncodedUintVector_UInt64 deltas_ = new EncodedUintVector_UInt64();
        private UInt64 base_;
        private byte shift_;
    }
}
