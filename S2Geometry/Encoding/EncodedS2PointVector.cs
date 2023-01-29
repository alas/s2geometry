// This class represents an encoded vector of S2Points.  Values are decoded
// only when they are accessed.  This allows for very fast initialization and
// no additional memory use beyond the encoded data.  The encoded data is not
// owned by this class; typically it points into a large contiguous buffer
// that contains other encoded data as well.
//
// This is one of several helper classes that allow complex data structures to
// be initialized from an encoded format inant time and then decoded on
// demand.  This can be a big performance advantage when only a small part of
// the data structure is actually used.

namespace S2Geometry;

public class EncodedS2PointVector
{
    #region Fields and Constants

    // To save space (especially for vectors of length 0, 1, and 2), the encoding
    // format is encoded in the low-order 3 bits of the vector size.  Up to 7
    // encoding formats are supported (only 2 are currently defined).  Additional
    // formats could be supported by using "7" as an overflow indicator and
    // encoding the actual format separately, but it seems unlikely we will ever
    // need to do that.
    public const int kEncodingFormatBits = 3;
    public const byte kEncodingFormatMask = (1 << kEncodingFormatBits) - 1;

    // S2CellIds are represented in a special 64-bit format and are encoded in
    // fixed-size blocks.  kBlockSize represents the number of values per block.
    // Block sizes of 4, 8, 16, and 32 were tested and kBlockSize == 16 seems to
    // offer the best compression.  (Note that kBlockSize == 32 requires some code
    // modifications which have since been removed.)
    public const int kBlockShift = 4;
    public const int kBlockSize = 1 << kBlockShift;

    // Used to indicate that a point must be encoded as an exception (a 24-byte
    // S2Point) rather than as an S2CellId.
    public const UInt64 kException = ~0UL;

    public required Format Format_ { private get; init; }
    private UInt32 size_;
    private S2Point[]? points;
    private EncodedStringVector? blocks;
    private UInt64 base_;
    private byte level;
    private bool have_exceptions;

    #endregion

    // Initializes the EncodedS2PointVector.
    //
    // REQUIRES: The Decoder data buffer must outlive this object.
    public static (bool Success, EncodedS2PointVector? Shape) Init(Decoder decoder)
    {
        if (decoder.Avail() < 1) return (false, null);

        // Peek at the format but don't advance the decoder; the format-specific
        // Init functions will do that.
        var res = new EncodedS2PointVector
        {
            Format_ = (Format)(decoder.Peek8() & kEncodingFormatMask)
        };
        var success = res.Format_ switch
        {
            Format.UNCOMPRESSED => res.InitUncompressedFormat(decoder),
            Format.CELL_IDS => res.InitCellIdsFormat(decoder),
            _ => false,
        };
        if (!success) return (false, null);

        return (true, res);
    }

    // Returns the size of the original vector.
    public int Count() => (int)size_;

    // Returns the element at the given index.
    public S2Point this[int i]
    {
        get
        {
            return Format_ switch
            {
                Format.UNCOMPRESSED => points[i],
                Format.CELL_IDS => DecodeCellIdsFormat(i),
                _ => throw new ApplicationException("Unrecognized format")
            };
        }
    }

    // Decodes and returns the entire original vector.
    public S2Point[] Decode()
    {
        var points = new S2Point[size_];
        for (int i = 0; i < size_; ++i)
        {
            points[i] = this[i];
        }
        return points;
    }

    public void Encode(Encoder encoder)
    {
        // The encoding must be identical to EncodeS2PointVector().
        switch (Format_)
        {
            case Format.UNCOMPRESSED:
                EncodeS2PointVectorFast(points, encoder);
                break;
            case Format.CELL_IDS:
                // This is a full decode/encode dance, and not at all efficient.
                EncodeS2PointVectorCompact(Decode(), encoder);
                break;
            default:
                throw new ApplicationException("Unknown Format: " + (int)Format_);
        }
    }

    // TODO(ericv): Consider adding a method that returns an adjacent pair of
    // points.  This would save some decoding overhead.


    // Encodes a vector of S2Points in a format that can later be decoded as an
    // EncodedS2PointVector.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static void EncodeS2PointVector(S2Point[] points, CodingHint hint, Encoder encoder)
    {
        switch (hint)
        {
            case CodingHint.FAST:
                EncodeS2PointVectorFast(points, encoder);
                break;
            case CodingHint.COMPACT:
                EncodeS2PointVectorCompact(points, encoder);
                break;
            default:
                throw new ApplicationException("Unknown CodingHint: " + (int)hint);
        }
    }

    // Encodes a vector of points, optimizing for (encoding and decoding) speed.
    private static void EncodeS2PointVectorFast(S2Point[] points, Encoder encoder)
    {
#if BIGENDIAN || ARM
            throw new NotImplementedException("Not implemented on big-endian architectures");
#endif

        // This function always uses the UNCOMPRESSED encoding.  The header consists
        // of a varint64 in the following format:
        //
        //   bits 0-2:  encoding format (UNCOMPRESSED)
        //   bits 3-63: vector size
        //
        // This is followed by an array of S2Points in little-endian order.
        encoder.Ensure(Encoder.kVarintMax64 + points.Length * SizeHelper.SizeOf(typeof(S2Point)));
        UInt64 size_format = (ulong)(points.Length << kEncodingFormatBits | (byte)Format.UNCOMPRESSED);
        encoder.Put64(size_format);
        encoder.PutPoints(points);
    }

    // Encodes a vector of points, optimizing for space.
    public static void EncodeS2PointVectorCompact(S2Point[] points, Encoder encoder)
    {
        // OVERVIEW
        // --------
        //
        // We attempt to represent each S2Point as the center of an S2CellId.  All
        // S2CellIds must be at the same level.  Any points that cannot be encoded
        // exactly as S2CellId centers are stored as exceptions using 24 bytes each.
        // If there are so many exceptions that the Format.CELL_IDS encoding does not save
        // significant space, we give up and use the uncompressed encoding.
        //
        // The first step is to choose the best S2CellId level.  This requires
        // converting each point to (face, si, ti) coordinates and checking whether
        // the point can be represented exactly as an S2CellId center at some level.
        // We then build a histogram of S2CellId levels (just like the similar code
        // in S2Polygon.Encode) and choose the best level (or give up, if there are
        // not enough S2CellId-encodable points).
        //
        // The simplest approach would then be to take all the S2CellIds and
        // right-shift them to remove all theant bits at the chosen level.
        // This would give the best spatial locality and hence the smallest deltas.
        // However instead we give up some spatial locality and use the similar but
        // faster transformation described below.
        //
        // Each encodable point is first converted to the (sj, tj) representation
        // defined below:
        //
        //   sj = (((face & 3) << 30) | (si >> 1)) >> (30 - level);
        //   tj = (((face & 4) << 29) | ti) >> (31 - level);
        //
        // These two values encode the (face, si, ti) tuple using (2 * level + 3)
        // bits.  To see this, recall that "si" and "ti" are 31-bit values that all
        // share a common suffix consisting of a "1" bit followed by (30 - level)
        // "0" bits.  The code above right-shifts these values to remove the
        // constant bits and then prepends the bits for "face", yielding a total of
        // (level + 2) bits for "sj" and (level + 1) bits for "tj".
        //
        // We then combine (sj, tj) into one 64-bit value by interleaving bit pairs:
        //
        //   v = InterleaveBitPairs(sj, tj);
        //
        // (We could also interleave individual bits, but it is faster this way.)
        // The result is similar to right-shifting an S2CellId by (61 - 2 * level),
        // except that it is faster to decode and the spatial locality is not quite
        // as good.
        //
        // The 64-bit values are divided into blocks of size kBlockSize, and then
        // each value is encoded as the sum of a base value, a per-block offset, and
        // a per-value delta within that block:
        //
        //   v[i,j] = base + offset[i] + delta[i, j]
        //
        // where "i" represents a block and "j" represents an entry in that block.
        //
        // The deltas for each block are encoded using a fixed number of 4-bit nibbles
        // (1-16 nibbles per delta).  This allows any delta to be accessed inant
        // time.
        //
        // The "offset" for each block is a 64-bit value encoded in 0-8 bytes.  The
        // offset is left-shifted such that it overlaps the deltas by a configurable
        // number of bits (either 0 or 4), called the "overlap".  The overlap and
        // offset length (0-8 bytes) are specified per block.  The reason for the
        // overlap is that it allows fewer delta bits to be used in some cases.  For
        // example if base == 0 and the range within a block is 0xf0 to 0x110, then
        // rather than using 12-bits deltas with an offset of 0, the overlap lets us
        // use 8-bits deltas with an offset of 0xf0 (saving 7 bytes per block).
        //
        // The global minimum value "base" is encoded using 0-7 bytes starting with
        // the most-significant non-zero bit possible for the chosen level.  For
        // example, if (level == 7) then the encoded values have at most 17 bits, so
        // if "base" is encoded in 1 byte then it is shifted to occupy bits 9-16.
        //
        // Example: at level == 15, there are at most 33 non-zero value bits.  The
        // following shows the bit positions covered by "base", "offset", and "delta"
        // assuming that "base" and "offset" are encoded in 2 bytes each, deltas are
        // encoded in 2 nibbles (1 byte) each, and "overlap" is 4 bits:
        //
        //   Base:             1111111100000000-----------------
        //   Offset:           -------------1111111100000000----
        //   Delta:            -------------------------00000000
        //   Overlap:                                   ^^^^
        //
        // The numbers (0 or 1) in this diagram denote the byte number of the encoded
        // value.  Notice that "base" is shifted so that it starts at the leftmost
        // possible bit, "delta" always starts at the rightmost possible bit (bit 0),
        // and "offset" is shifted so that it overlaps "delta" by the chosen "overlap"
        // (either 0 or 4 bits).  Also note that all of these values are summed, and
        // therefore each value can affect higher-order bits due to carries.
        //
        // NOTE(ericv): Encoding deltas in 4-bit rather than 8-bit length increments
        // reduces encoded sizes by about 7%.  Allowing a 4-bit overlap between the
        // offset and deltas reduces encoded sizes by about 1%.  Both optimizations
        // make the code more complex but don't affect running times significantly.
        //
        // ENCODING DETAILS
        // ----------------
        //
        // Now we can move on to the actual encodings.  First, there is a 2 byte
        // header encoded as follows:
        //
        //  Byte 0, bits 0-2: encoding_format (Format.CELL_IDS)
        //  Byte 0, bit  3:   have_exceptions
        //  Byte 0, bits 4-7: (last_block_size - 1)
        //  Byte 1, bits 0-2: base_bytes
        //  Byte 1, bits 3-7: level (0-30)
        //
        // This is followed by an EncodedStringVector containing the encoded blocks.
        // Each block contains kBlockSize (8) values.  The total size of the
        // EncodeS2PointVector is not stored explicity, but instead is calculated as
        //
        //     num_values == kBlockSize * (num_blocks - 1) + last_block_size .
        //
        // (An empty vector has num_blocks == 0 and last_block_size == kBlockSize.)
        //
        // Each block starts with a 1 byte header containing the following:
        //
        //  Byte 0, bits 0-2: (offset_bytes - overlap_nibbles)
        //  Byte 0, bit  3:   overlap_nibbles
        //  Byte 0, bits 4-7: (delta_nibbles - 1)
        //
        // "overlap_nibbles" is either 0 or 1 (indicating an overlap of 0 or 4 bits),
        // while "offset_bytes" is in the range 0-8 (indicating the number of bytes
        // used to encode the offset for this block).  Note that some combinations
        // cannot be encoded: in particular, offset_bytes == 0 can only be encoded
        // with an overlap of 0 bits, and offset_bytes == 8 can only be encoded with
        // an overlap of 4 bits.  This allows us to encode offset lengths of 0-8
        // rather than just 0-7 without using an extra bit.  (Note that the
        // combinations that can't be encoded are not useful anyway.)
        //
        // The header is followed by "offset_bytes" bytes for the offset, and then
        // (4 * delta_nibbles) bytes for the deltas.
        //
        // If there are any points that could not be represented as S2CellIds, then
        // "have_exceptions" in the header is true.  In that case the delta values
        // within each block are encoded as (delta + kBlockSize), and values
        // 0...kBlockSize-1 are used to represent exceptions.  If a block has
        // exceptions, they are encoded immediately following the array of deltas,
        // and are referenced by encoding the corresponding exception index
        // 0...kBlockSize-1 as the delta.
        //
        // TODO(ericv): A vector containing a single leaf cell is currently encoded as
        // 13 bytes (2 byte header, 7 byte base, 1 byte block count, 1 byte block
        // length, 1 byte block header, 1 byte delta).  However if this case occurs
        // often, a better solution would be implement a separate format that encodes
        // the leading k bytes of an S2CellId.  It would have a one-byte header
        // consisting of the encoding format (3 bits) and the number of bytes encoded
        // (3 bits), followed by the S2CellId bytes.  The extra 2 header bits could be
        // used to store single points using other encodings, e.g. E7.
        //
        // If we had used 8-value blocks, we could have used the extra bit in the
        // first byte of the header to indicate that there is only one value, and
        // then skip the 2nd byte of header and the EncodedStringVector.  But this
        // would be messy because it also requires special cases while decoding.
        // Essentially this would be a sub-format within the Format.CELL_IDS format.

        // 1. Compute (level, face, si, ti) for each point, build a histogram of
        // levels, and determine the optimal level to use for encoding (if any).
        int level = ChooseBestLevel(points, out CellPoint[] cell_points);
        if (level < 0)
        {
            EncodeS2PointVectorFast(points, encoder);
            return;
        }

        // 2. Convert the points into encodable 64-bit values.  We don't use the
        // S2CellId itself because it requires a somewhat more complicated bit
        // interleaving operation.
        //
        // TODO(ericv): Benchmark using shifted S2CellIds instead.
        UInt64[] values = ConvertCellsToValues(cell_points, level, out bool have_exceptions);

        // 3. Choose the global encoding parameter "base" (consisting of the bit
        // prefix shared by all values to be encoded).
        UInt64 base_ = ChooseBase(values, level, have_exceptions, out int base_bits);

        // Now encode the output, starting with the 2-byte header (see above).
        int num_blocks = (values.Length + kBlockSize - 1) >> kBlockShift;
        int base_bytes = base_bits >> 3;
        encoder.Ensure(2 + base_bytes);
        int last_block_count = values.Length - kBlockSize * (num_blocks - 1);
        MyDebug.Assert(last_block_count >= 0);
        MyDebug.Assert(last_block_count <= kBlockSize);
        MyDebug.Assert(base_bytes <= 7);
        MyDebug.Assert(level <= 30);
        encoder.Put8((byte)(
            (byte)Format.CELL_IDS |
            ((have_exceptions ? 1 : 0) << 3) |
            ((last_block_count - 1) << 4)));
        encoder.Put8((byte)(base_bytes | (level << 3)));

        // Next we encode 0-7 bytes of "base".
        int base_shift = BaseShift(level, base_bits);
        EncodedUIntVector<byte>.EncodeUIntWithLength((byte)(base_ >> base_shift), base_bytes, encoder);

        // Now we encode the contents of each block.
        var blocks = new StringVectorEncoder();
        var exceptions = new List<S2Point>();
        UInt64 offset_bytes_sum = 0UL;
        UInt64 delta_nibbles_sum = 0UL;
        UInt64 exceptions_sum = 0UL;
        for (int i = 0; i < values.Length; i += kBlockSize)
        {
            int block_size = Math.Min(kBlockSize, values.Length - i);
            var code = GetBlockCode(new ArraySegment<UInt64>(values, i, block_size), base_, have_exceptions);

            // Encode the one-byte block header (see above).
            Encoder block = blocks.AddViaEncoder();
            int offset_bytes = code.OffsetBits >> 3;
            int delta_nibbles = code.DeltaBits >> 2;
            int overlap_nibbles = code.OverlapBits >> 2;
            block.Ensure(1 + offset_bytes + (kBlockSize / 2) * delta_nibbles);
            MyDebug.Assert(offset_bytes - overlap_nibbles <= 7);
            MyDebug.Assert(overlap_nibbles <= 1);
            MyDebug.Assert(delta_nibbles <= 16);
            block.Put8((byte)((offset_bytes - overlap_nibbles) | (overlap_nibbles << 3) | (delta_nibbles - 1) << 4));

            // Determine the offset for this block, and whether there are exceptions.
            UInt64 offset = ~0UL;
            int num_exceptions = 0;
            for (int j = 0; j < block_size; ++j)
            {
                if (values[i + j] == kException)
                {
                    num_exceptions += 1;
                }
                else
                {
                    MyDebug.Assert(values[i + j] >= base_);
                    offset = Math.Min(offset, values[i + j] - base_);
                }
            }
            if (num_exceptions == block_size) offset = 0;

            // Encode the offset.
            int offset_shift = code.DeltaBits - code.OverlapBits;
            offset &= ~BitMask(offset_shift);
            MyDebug.Assert((offset == 0) == (offset_bytes == 0));
            if (offset > 0)
            {
                EncodedUIntVector<ulong>.EncodeUIntWithLength(offset >> offset_shift, offset_bytes, block);
            }

            // Encode the deltas, and also gather any exceptions present.
            int delta_bytes = (delta_nibbles + 1) >> 1;
            exceptions.Clear();
            for (int j = 0; j < block_size; ++j)
            {
                UInt64 delta;
                if (values[i + j] == kException)
                {
                    delta = (ulong)exceptions.Count;
                    exceptions.Add(points[i + j]);
                }
                else
                {
                    MyDebug.Assert(values[i + j] >= offset + base_);
                    delta = values[i + j] - (offset + base_);
                    if (have_exceptions)
                    {
                        MyDebug.Assert(delta <= (~0UL - kBlockSize));
                        delta += kBlockSize;
                    }
                }
                MyDebug.Assert(delta <= BitMask(code.DeltaBits));
                if ((delta_nibbles & 1) != 0 && (j & 1) != 0)
                {
                    // Combine this delta with the high-order 4 bits of the previous delta.
                    byte last_byte = block.Last();
                    block.RemoveLast(1);
                    delta = (delta << 4) | (uint)(last_byte & 0xf);
                }
                EncodedUIntVector<ulong>.EncodeUIntWithLength(delta, delta_bytes, block);
            }
            // Append any exceptions to the end of the block.
            if (num_exceptions > 0)
            {
                int exceptions_bytes = exceptions.Count * SizeHelper.SizeOf(typeof(S2Point));
                block.Ensure(exceptions_bytes);
                block.PutPoints(exceptions);
            }
            offset_bytes_sum += (UInt64)offset_bytes;
            delta_nibbles_sum += (UInt64)delta_nibbles;
            exceptions_sum += (UInt64)num_exceptions;
        }
        blocks.Encode(encoder);
    }

    private bool InitUncompressedFormat(Decoder decoder)
    {
#if BIGENDIAN || ARM
            // TODO(b/231674214): Make this work on platforms that don't support
            // unaligned 64-bit little-endian reads, e.g. by falling back to
            //
            //   bit_cast<double>(little_endian.Load64()).
            //
            // Maybe the compiler is smart enough that we can do this all the time,
            // but more likely we will need two cases using the #ifdef above.
            // (Note that even ARMv7 does not support unaligned 64-bit loads.)
            throw new NotImplementedException("Needs architecture with 64-bit little-endian unaligned loads");
#endif
        if (!decoder.TryGetVarUInt64(out var size)) return false;
        size >>= kEncodingFormatBits;

        // Note that the encoding format supports up to 2**59 vertices, but we
        // currently only support decoding up to 2**32 vertices.
        if (size > UInt32.MaxValue) return false;
        size_ = (uint)size;

        int bytes = (int)(size_ * SizeHelper.SizeOf(typeof(S2Point)));
        if (decoder.Avail() < bytes) return false;

        var pointBuffer = new S2Point[size_];
        decoder.GetPoints(pointBuffer, 0, (int)size_);
        points = pointBuffer;
        //decoder.Skip(bytes); //note(Alas): in the original code GetPoints was a direct buffer access without advancing the pointer so here there is no need to skip
        return true;
    }
    private bool InitCellIdsFormat(Decoder decoder)
    {
        // This function inverts the encodings documented above.
        // First we decode the two-byte header.
        if (decoder.Avail() < 2) return false;
        byte header1 = decoder.Get8();
        byte header2 = decoder.Get8();
        MyDebug.Assert((header1 & 7) == (byte)Format.CELL_IDS);
        int last_block_count, base_bytes;
        have_exceptions = (header1 & 8) != 0;
        last_block_count = (header1 >> 4) + 1;
        base_bytes = header2 & 7;
        level = (byte)(header2 >> 3);

        // Decode the base value (if any).
        if (!EncodedUIntVector<ulong>.DecodeUIntWithLength(base_bytes, decoder, out ulong base_)) return false;
        this.base_ = base_ << BaseShift(level, base_bytes << 3);

        // Initialize the vector of encoded blocks.
        var (success, shape) = EncodedStringVector.Init(decoder);
        if (!success) return false;

        blocks = shape!;
        size_ = (UInt32)(kBlockSize * (blocks.Count() - 1) + last_block_count);
        return true;
    }
    private S2Point DecodeCellIdsFormat(int i)
    {
        // This function inverts the encodings documented above.

        // First we decode the block header.
        var (bytearr, offsetbytearr) = blocks.GetStart(i >> kBlockShift);
        byte header = bytearr[offsetbytearr++];
        int overlap_nibbles = (header >> 3) & 1;
        int offset_bytes = (header & 7) + overlap_nibbles;
        int delta_nibbles = (header >> 4) + 1;

        // Decode the offset for this block.
        int offset_shift = (delta_nibbles - overlap_nibbles) << 2;
        UInt64 offset = EncodedUIntVector<ulong>.GetUIntWithLength(bytearr, offsetbytearr, offset_bytes) << offset_shift;
        offsetbytearr += offset_bytes;

        // Decode the delta for the requested value.
        int delta_nibble_offset = (i & (kBlockSize - 1)) * delta_nibbles;
        int delta_bytes = (delta_nibbles + 1) >> 1;
        var delta_offset = offsetbytearr + (delta_nibble_offset >> 1);
        UInt64 delta = EncodedUIntVector<ulong>.GetUIntWithLength(bytearr, delta_offset, delta_bytes);
        delta >>= (delta_nibble_offset & 1) << 2;
        delta &= BitMask(delta_nibbles << 2);

        // Test whether this point is encoded as an exception.
        if (have_exceptions)
        {
            if (delta < kBlockSize)
            {
                int block_size = Math.Min(kBlockSize, (int)(size_ - (uint)(i & ~(kBlockSize - 1))));
                offsetbytearr += (block_size * delta_nibbles + 1) >> 1;
                offsetbytearr += (int)delta * SizeHelper.SizeOf(typeof(S2Point));
                var buff = new double[3];
                Buffer.BlockCopy(bytearr, offsetbytearr, buff, 0, SizeHelper.SizeOf(typeof(S2Point)));
                return new S2Point(buff);
            }
            delta -= kBlockSize;
        }

        // Otherwise convert the 64-bit value back to an S2Point.
        UInt64 value = base_ + offset + delta;
        int shift = S2.kMaxCellLevel - level;

        // The S2CellId version of the following code is:
        //   return S2CellId(((value << 1) | 1) << (2 * shift)).ToPoint();
        DeinterleaveUInt32BitPairs(value, out uint sj, out uint tj);
        int si = (int)((((sj << 1) | 1) << shift) & 0x7fffffff);
        int ti = (int)((((tj << 1) | 1) << shift) & 0x7fffffff);
        int face = (int)(((sj << shift) >> 30) | (((tj << (shift + 1)) >> 29) & 4));
        return S2.FaceUVtoXYZ(face, S2.STtoUV(S2.SiTitoST((uint)si)),
                               S2.STtoUV(S2.SiTitoST((uint)ti))).Normalize();
    }

    // Like util_bits.InterleaveUInt32, but interleaves bit pairs rather than
    // individual bits.  This format is faster to decode than the fully interleaved
    // format, and produces the same results for our use case.
    private static UInt64 InterleaveUInt32BitPairs(UInt32 val0, UInt32 val1)
    {
        UInt64 v0 = val0, v1 = val1;
        v0 = (v0 | (v0 << 16)) & 0x0000ffff0000ffff;
        v1 = (v1 | (v1 << 16)) & 0x0000ffff0000ffff;
        v0 = (v0 | (v0 << 8)) & 0x00ff00ff00ff00ff;
        v1 = (v1 | (v1 << 8)) & 0x00ff00ff00ff00ff;
        v0 = (v0 | (v0 << 4)) & 0x0f0f0f0f0f0f0f0f;
        v1 = (v1 | (v1 << 4)) & 0x0f0f0f0f0f0f0f0f;
        v0 = (v0 | (v0 << 2)) & 0x3333333333333333;
        v1 = (v1 | (v1 << 2)) & 0x3333333333333333;
        return v0 | (v1 << 2);
    }

    // This code is about 50% faster than util_bits.DeinterleaveUInt32, which
    // uses a lookup table.  The speed advantage is expected to be even larger in
    // code that mixes bit interleaving with other significant operations since it
    // doesn't require keeping a 256-byte lookup table in the L1 data cache.
    private static void DeinterleaveUInt32BitPairs(UInt64 code, out UInt32 val0, out UInt32 val1)
    {
        UInt64 v0 = code, v1 = code >> 2;
        v0 &= 0x3333333333333333;
        v0 |= v0 >> 2;
        v1 &= 0x3333333333333333;
        v1 |= v1 >> 2;
        v0 &= 0x0f0f0f0f0f0f0f0f;
        v0 |= v0 >> 4;
        v1 &= 0x0f0f0f0f0f0f0f0f;
        v1 |= v1 >> 4;
        v0 &= 0x00ff00ff00ff00ff;
        v0 |= v0 >> 8;
        v1 &= 0x00ff00ff00ff00ff;
        v1 |= v1 >> 8;
        v0 &= 0x0000ffff0000ffff;
        v0 |= v0 >> 16;
        v1 &= 0x0000ffff0000ffff;
        v1 |= v1 >> 16;
        val0 = (UInt32)v0;
        val1 = (UInt32)v1;
    }

    // Returns a bit mask with "n" low-order 1 bits, for 0 <= n <= 64.
    private static UInt64 BitMask(int n) => (n == 0) ? 0 : (~0UL >> (64 - n));

    // Returns the maximum number of bits per value at the given S2CellId level.
    private static int MaxBitsForLevel(int level) => 2 * level + 3;

    // Returns the number of bits that "base" should be right-shifted in order to
    // encode only its leading "base_bits" bits, assuming that all points are
    // encoded at the given S2CellId level.
    private static int BaseShift(int level, int base_bits) => Math.Max(0, MaxBitsForLevel(level) - base_bits);

    // Returns the S2CellId level for which the greatest number of the given points
    // can be represented as the center of an S2CellId.  Initializes "cell_points"
    // to contain the S2CellId representation of each point (if any).  Returns -1
    // if there is no S2CellId that would result in significant space savings.
    private static int ChooseBestLevel(S2Point[] points, out CellPoint[] cell_points)
    {
        cell_points = new CellPoint[points.Length];

        // Count the number of points at each level.
        var level_counts = new int[S2.kMaxCellLevel + 1].Fill(0);
        for (var i = 0; i < points.Length; i++)
        {
            int level = S2.XYZtoFaceSiTi(points[i], out int face, out uint si, out uint ti);
            cell_points[i] = new CellPoint(level, face, si, ti);
            if (level >= 0) ++level_counts[level];
        }
        // Choose the level for which the most points can be encoded.
        int best_level = 0;
        for (int level = 1; level <= S2.kMaxCellLevel; ++level)
        {
            if (level_counts[level] > level_counts[best_level])
            {
                best_level = level;
            }
        }
        // The uncompressed encoding is smaller *and* faster when very few of the
        // points are encodable as S2CellIds.  The Format.CELL_IDS encoding uses about 1
        // extra byte per point in this case, consisting of up to a 3 byte
        // EncodedStringVector offset for each block, a 1 byte block header, and 4
        // bits per delta (encoding an exception number from 0-7), for a total of 8
        // bytes per block.  This represents a space overhead of about 4%, so we
        // require that at least 5% of the input points should be encodable as
        // S2CellIds in order for the Format.CELL_IDS format to be worthwhile.
        const double kMinEncodableFraction = 0.05;
        if (level_counts[best_level] <= kMinEncodableFraction * points.Length)
        {
            return -1;
        }
        return best_level;
    }

    // Given a vector of points in CellPoint format and an S2CellId level that has
    // been chosen for encoding, returns a vector of 64-bit values that should be
    // encoded in order to represent these points.  Points that cannot be
    // represented losslessly as the center of an S2CellId at the chosen level are
    // indicated by the value "kException".  "have_exceptions" is set to indicate
    // whether any exceptions were present.
    private static UInt64[] ConvertCellsToValues(CellPoint[] cell_points, int level, out bool have_exceptions)
    {
        var values = new UInt64[cell_points.Length];
        have_exceptions = false;
        int shift = S2.kMaxCellLevel - level;
        for (var i = 0; i < cell_points.Length; i++)
        {
            var cp = cell_points[i];
            if (cp.Level != level)
            {
                values[i] = kException;
                have_exceptions = true;
            }
            else
            {
                // Note that bit 31 of tj is always zero, and that bits are interleaved in
                // such a way that bit 63 of the result is always zero.
                //
                // The S2CellId version of the following code is:
                // UInt64 v = S2CellId.FromFaceIJ(cp.face, cp.si >> 1, cp.ti >> 1).
                //            parent(level).id() >> (2 * shift + 1);
                UInt32 sj = (((UInt32)((cp.Face & 3) << 30)) | (cp.Si >> 1)) >> shift;
                UInt32 tj = (((UInt32)((cp.Face & 4) << 29)) | cp.Ti) >> (shift + 1);
                UInt64 v = InterleaveUInt32BitPairs(sj, tj);
                MyDebug.Assert(v <= BitMask(MaxBitsForLevel(level)));
                values[i] = v;
            }
        }
        return values;
    }

    private static UInt64 ChooseBase(UInt64[] values, int level, bool have_exceptions, out int base_bits)
    {
        // Find the minimum and maximum non-exception values to be represented.
        UInt64 v_min = kException, v_max = 0;
        foreach (var v in values)
        {
            if (v != kException)
            {
                v_min = Math.Min(v_min, v);
                v_max = Math.Max(v_max, v);
            }
        }
        if (v_min == kException)
        {
            base_bits = 0;
            return 0;
        }

        // Generally "base" is chosen as the bit prefix shared by v_min and v_max.
        // However there are a few adjustments we need to make.
        //
        // 1. Encodings are usually smaller if the bits represented by "base" and
        // "delta" do not overlap.  Usually the shared prefix rule does this
        // automatically, but if v_min == v_max or there are special circumstances
        // that increase delta_bits (such as values.size() == 1) then we need to
        // make an adjustment.
        //
        // 2. The format only allows us to represent up to 7 bytes (56 bits) of
        // "base", so we need to ensure that "base" conforms to this requirement.
        int min_delta_bits = (have_exceptions || values.Length == 1) ? 8 : 4;
        int excluded_bits = Math.Max(BitsUtils.GetBitWidth(v_min ^ v_max),
                                Math.Max(min_delta_bits, BaseShift(level, 56)));
        UInt64 base_ = v_min & ~BitMask(excluded_bits);

        // Determine how many bytes are needed to represent this prefix.
        if (base_ == 0)
        {
            base_bits = 0;
        }
        else
        {
            int low_bit = BitsUtils.FindLSBSetNonZero64(base_);
            base_bits = (MaxBitsForLevel(level) - low_bit + 7) & ~7;
        }

        // Since base_bits has been rounded up to a multiple of 8, we may now be
        // able to represent additional bits of v_min.  In general this reduces the
        // final encoded size.
        //
        // NOTE(ericv): A different strategy for choosing "base" is to encode all
        // blocks under the assumption that "base" equals v_min exactly, and then
        // set base equal to the minimum-length prefix of "v_min" that allows these
        // encodings to be used.  This strategy reduces the encoded sizes by
        // about 0.2% relative to the strategy here, but is more complicated.
        return v_min & ~BitMask(BaseShift(level, base_bits));
    }

    // Returns true if the range of values [d_min, d_max] can be encoded using the
    // specified parameters (delta_bits, overlap_bits, and have_exceptions).
    private static bool CanEncode(UInt64 d_min, UInt64 d_max, int delta_bits, int overlap_bits, bool have_exceptions)
    {
        // "offset" can't represent the lowest (delta_bits - overlap_bits) of d_min.
        d_min &= ~BitMask(delta_bits - overlap_bits);

        // The maximum delta is reduced by kBlockSize if any exceptions exist, since
        // deltas 0..kBlockSize-1 are used to indicate exceptions.
        UInt64 max_delta = BitMask(delta_bits);
        if (have_exceptions)
        {
            if (max_delta < kBlockSize) return false;
            max_delta -= kBlockSize;
        }
        // The first test below is necessary to avoid 64-bit overflow.
        return (d_min > ~max_delta) || (d_min + max_delta >= d_max);
    }

    // Given a vector of 64-bit values to be encoded and an S2CellId level, returns
    // the optimal encoding parameters that should be used to encode each block.
    // Also returns the global minimum value "base_" and the number of bits that
    // should be used to encode it ("base_bits").
    private static BlockCode GetBlockCode(ArraySegment<UInt64> values, UInt64 base_, bool have_exceptions)
    {
        // "b_min" and "b_max"n are the minimum and maximum values within this block.
        UInt64 b_min = kException, b_max = 0;
        foreach (UInt64 v in values)
        {
            if (v != kException)
            {
                b_min = Math.Min(b_min, v);
                b_max = Math.Max(b_max, v);
            }
        }
        if (b_min == kException)
        {
            // All values in this block are exceptions.
            return new BlockCode(4, 0, 0);
        }

        // Adjust the min/max values so that they are relative to "base_".
        b_min -= base_;
        b_max -= base_;

        // Determine the minimum possible delta length and overlap that can be used
        // to encode this block.  The block will usually be encodable using the
        // number of bits in (b_max - b_min) rounded up to a multiple of 4.  If this
        // is not possible, the preferred solution is to shift "offset" so that the
        // delta and offset values overlap by 4 bits (since this only costs an
        // average of 4 extra bits per block).  Otherwise we increase the delta size
        // by 4 bits.  Certain cases require that both of these techniques are used.
        //
        // Example 1: b_min = 0x72, b_max = 0x7e.  The range is 0x0c.  This can be
        // encoded using delta_bits = 4 and overlap_bits = 0, which allows us to
        // represent an offset of 0x70 and a maximum delta of 0x0f, so that we can
        // encode values up to 0x7f.
        //
        // Example 2: b_min = 0x78, b_max = 0x84.  The range is 0x0c, but in this
        // case it is not sufficient to use delta_bits = 4 and overlap_bits = 0
        // because we can again only represent an offset of 0x70, so the maximum
        // delta of 0x0f only lets us encode values up to 0x7f.  However if we
        // increase the overlap to 4 bits then we can represent an offset of 0x78,
        // which lets us encode values up to 0x78 + 0x0f = 0x87.
        //
        // Example 3: b_min = 0x08, b_max = 0x104.  The range is 0xfc, so we should
        // be able to use 8-bit deltas.  But even with a 4-bit overlap, we can still
        // only encode offset = 0 and a maximum value of 0xff.  (We don't allow
        // bigger overlaps because statistically they are not worthwhile.)  Instead
        // we increase the delta size to 12 bits, which handles this case easily.
        //
        // Example 4: b_min = 0xf08, b_max = 0x1004.  The range is 0xfc, so we
        // should be able to use 8-bit deltas.  With 8-bit deltas and no overlap, we
        // have offset = 0xf00 and a maximum encodable value of 0xfff.  With 8-bit
        // deltas and a 4-bit overlap, we still have offset = 0xf00 and a maximum
        // encodable value of 0xfff.  Even with 12-bit deltas, we have offset = 0
        // and we can still only represent 0xfff.  However with delta_bits = 12 and
        // overlap_bits = 4, we can represent offset = 0xf00 and a maximum encodable
        // value of 0xf00 + 0xfff = 0x1eff.
        //
        // It is possible to show that this last example is the worst case, i.e.  we
        // do not need to consider increasing delta_bits or overlap_bits further.
        int delta_bits = (Math.Max(1, BitsUtils.GetBitWidth(b_max - b_min) - 1) + 3) & ~3;
        int overlap_bits = 0;
        if (!CanEncode(b_min, b_max, delta_bits, 0, have_exceptions))
        {
            if (CanEncode(b_min, b_max, delta_bits, 4, have_exceptions))
            {
                overlap_bits = 4;
            }
            else
            {
                MyDebug.Assert(delta_bits <= 60);
                delta_bits += 4;
                if (!CanEncode(b_min, b_max, delta_bits, 0, have_exceptions))
                {
                    MyDebug.Assert(CanEncode(b_min, b_max, delta_bits, 4, have_exceptions));
                    overlap_bits = 4;
                }
            }
        }

        // When the block size is 1 and no exceptions exist, we have delta_bits == 4
        // and overlap_bits == 0 which wastes 4 bits.  We fix this below, which
        // among other things reduces the encoding size for single leaf cells by one
        // byte.  (Note that when exceptions exist, delta_bits == 8 and overlap_bits
        // may be 0 or 4.  These cases are covered by the unit tests.)
        if (values.Count == 1 && !have_exceptions)
        {
            MyDebug.Assert(delta_bits == 4 && overlap_bits == 0);
            delta_bits = 8;
        }

        // Now determine the number of bytes needed to encode "offset", given the
        // chosen delta length.
        UInt64 max_delta = BitMask(delta_bits) - (UInt64)(have_exceptions ? kBlockSize : 0);
        int offset_bits = 0;
        if (b_max > max_delta)
        {
            // At least one byte of offset is required.  Round up the minimum offset
            // to the next encodable value, and determine how many bits it has.
            int offset_shift = delta_bits - overlap_bits;
            UInt64 mask = BitMask(offset_shift);
            UInt64 min_offset = (b_max - max_delta + mask) & ~mask;
            MyDebug.Assert(min_offset > 0);
            offset_bits = (BitsUtils.FindMSBSetNonZero64(min_offset) + 1 - offset_shift + 7) & ~7;
            // A 64-bit offset can only be encoded with an overlap of 4 bits.
            if (offset_bits == 64) overlap_bits = 4;
        }
        return new BlockCode(delta_bits, offset_bits, overlap_bits);
    }

    // TODO(ericv): Use atomic_flag to cache the last point decoded in
    // a thread-safe way.  This reduces benchmark times for actual polygon
    // operations (e.g. S2ClosestEdgeQuery) by about 15%.

    // We use a tagged union to represent multiple formats, as opposed to an
    // abstract base class or templating.  This represents the best compromise
    // between performance, space, and convenience.  Note that the overhead of
    // checking the tag is trivial and will typically be branch-predicted
    // perfectly.
    //
    // TODO(ericv): Once additional formats have been implemented, consider
    // using variant<> instead.  It's unclear whether this would have
    // better or worse performance than the current approach.
    public enum Format : byte
    {
        UNCOMPRESSED = 0,
        CELL_IDS = 1,
    }

    // Represents a point that can be encoded as an S2CellId center.
    // (If such an encoding is not possible then level < 0.)
    private readonly struct CellPoint
    {
        public sbyte Level { get; }
        public sbyte Face { get; }
        public UInt32 Si { get; }
        public UInt32 Ti { get; }

        // Constructor necessary in order to narrow "int" arguments to "sbyte".
        public CellPoint(int level, int face, UInt32 si, UInt32 ti)
        {
            Level = (sbyte)level;
            Face = (sbyte)face;
            Si = si;
            Ti = ti;
        }
    }

    // Represents the encoding parameters to be used for a given block (consisting
    // of kBlockSize encodable 64-bit values).  See below.
    private readonly record struct BlockCode(
        int DeltaBits,     // Delta length in bits (multiple of 4)
        int OffsetBits,    // Offset length in bits (multiple of 8)
        int OverlapBits);  // {Delta, Offset} overlap in bits (0 or 4)
}

