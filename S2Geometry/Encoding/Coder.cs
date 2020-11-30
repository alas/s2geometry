using System;
using System.Collections.Generic;

namespace S2Geometry
{
    // Class for encoding data into a memory buffer
    public class Encoder
    {
#pragma warning disable IDE1006 // Estilos de nombres
        // Support for variable length encoding with 7 bits per byte
        // (these are just simple wrappers around the Varint module)
        public const int kVarintMax32 = 5;
        public const int kVarintMax64 = 10;
#pragma warning restore IDE1006 // Estilos de nombres

        static Encoder()
        {
            Assert.True(sizeof(Byte) == 1);
            Assert.True(sizeof(UInt16) == 2);
            Assert.True(sizeof(UInt32) == 4);
            Assert.True(sizeof(UInt64) == 8);
        }

        // Initialize encoder to encode into "buf"
        public Encoder()
        {
            Reset(Array.Empty<byte>(), 0);
        }

        public Encoder(byte[] buf, int maxn)
        {
            Reset(buf, maxn);
        }

        public void Reset(byte[] buf, int maxn)
        {
            Buffer = buf;
            limit_ = maxn;
            offset_ = 0;
        }
        public void Clear()
        {
            offset_ = 0;
        }

        // Encoding routines.  Note that these do not check bounds
        public void Put8(sbyte v) { Put8((byte)v); }
        public void Put8(byte v)
        {
            Assert.True(Avail()>= sizeof(byte));
            Buffer[offset_] = v;
            offset_ += sizeof(byte);
        }
        public void Put16(Int16 v) { Put16((UInt16)v); }
        public void Put16(UInt16 v)
        {
            Assert.True(Avail() >= sizeof(UInt16));
            var bytes = BitConverter.GetBytes(v);
            Buffer[offset_] = bytes[0];
            Buffer[++offset_] = bytes[1];
        }
        public void Put32(Int32 v) { Put32((UInt32)v); }
        public void Put32(UInt32 v)
        {
            Assert.True(Avail() >= sizeof(UInt32));
            var bytes = BitConverter.GetBytes(v);
            Buffer[offset_] = bytes[0];
            Buffer[offset_ + 1] = bytes[1];
            Buffer[offset_ + 2] = bytes[2];
            Buffer[offset_ + 3] = bytes[3];
            offset_ += 4;
        }
        public void Put64(Int64 v) { Put64((UInt64)v); }
        public void Put64(UInt64 v)
        {
            Assert.True(Avail() >= sizeof(UInt64));
            var bytes = BitConverter.GetBytes(v);
            Buffer[offset_] = bytes[0];
            Buffer[offset_ + 1] = bytes[1];
            Buffer[offset_ + 2] = bytes[2];
            Buffer[offset_ + 3] = bytes[3];
            Buffer[offset_ + 4] = bytes[4];
            Buffer[offset_ + 5] = bytes[5];
            Buffer[offset_ + 6] = bytes[6];
            Buffer[offset_ + 7] = bytes[7];
            offset_ += 8;
        }
        public void PutN(byte[] mem)
        {
            PutN(mem, mem.Length);
        }
        public void PutN(byte[] mem, int n)
        {
            System.Buffer.BlockCopy(mem, 0, Buffer, offset_, n);
            offset_ += n;
        }
        public void PutPoints(IEnumerable<S2Point> points)
        {
            foreach (var v in points)
            {
                PutDouble(v[0]);
                PutDouble(v[1]);
                PutDouble(v[2]);
            }
        }
        public void PutFloat(float f)
        {
            Assert.True(Avail() >= sizeof(float));
            var bytes = BitConverter.GetBytes(f);
            Buffer[offset_] = bytes[0];
            Buffer[offset_ + 1] = bytes[1];
            Buffer[offset_ + 2] = bytes[2];
            Buffer[offset_ + 3] = bytes[3];
            offset_ += 4;
        }
        public void PutDouble(double d)
        {
            Assert.True(Avail() >= sizeof(double));
            var bytes = BitConverter.GetBytes(d);
            Buffer[offset_] = bytes[0];
            Buffer[offset_ + 1] = bytes[1];
            Buffer[offset_ + 2] = bytes[2];
            Buffer[offset_ + 3] = bytes[3];
            Buffer[offset_ + 4] = bytes[4];
            Buffer[offset_ + 5] = bytes[5];
            Buffer[offset_ + 6] = bytes[6];
            Buffer[offset_ + 7] = bytes[7];
            offset_ += 8;
        }

        public void PutVarInt32(int v)
        {
            VarintBitConverter.GetVarintBytes(v, Buffer, ref offset_);
        }
        public void PutVarUInt32(uint v)
        {
            VarintBitConverter.GetVarintBytes(v, Buffer, ref offset_);
        }
        public void PutVarInt64(long v)
        {
            VarintBitConverter.GetVarintBytes(v, Buffer, ref offset_);
        }
        public void PutVarUInt64(ulong v)
        {
            VarintBitConverter.GetVarintBytes(v, Buffer, ref offset_);
        }

        // Return number of bytes encoded so far
        public int Length
        {
            get
            {
                Assert.True(offset_ >= 0);
                Assert.True(offset_ <= limit_);  // Catch the buffer overflow.
                return offset_;
            }
        }

        // Return number of bytes of space remaining in buffer
        public int Avail()
        {
            Assert.True(limit_ >= offset_);
            return limit_ - offset_;
        }

        // REQUIRES: Encoder was created with the 0-argumentructor interface.
        //
        // This interface ensures that at least "N" more bytes are available
        // in the underlying buffer by resizing the buffer (if necessary).
        //
        // Note that no bounds checking is done on any of the put routines,
        // so it is the client's responsibility to call Ensure() at
        // appropriate intervals to ensure that enough space is available
        // for the data being added.
        public void Ensure(int N)
        {
            if (Avail() < N)
            {
                EnsureSlowPath(N);
            }
        }

        // Advances the write pointer by "N" bytes.
        public void Skip(int N)
        {
            offset_ += N;
        }

        // REQUIRES: length() >= N
        // Removes the last N bytes out of the encoded buffer
        public void RemoveLast(int N)
        {
            Assert.True(Length>= N);
            offset_ -= N;
        }

        public byte Last() => Buffer[offset_ - 1];

        private void EnsureSlowPath(int N)
        {
            Assert.True(Avail() < N);

            // Double buffer size, but make sure we always have at least N extra bytes
            int current_len = Length;
            int new_capacity = Math.Max(current_len + N, 2 * current_len);

            var new_buffer = new byte[new_capacity];
            if (Buffer != null)
            {
                System.Buffer.BlockCopy(Buffer, 0, new_buffer, 0, current_len);
            }
            Buffer = new_buffer;

            limit_ = new_capacity;
            Assert.True(Avail() >= N);
        }

        public byte[] Buffer { get; private set; }
        private int limit_;
        private int offset_;
    }

    public interface ICoder
    {
        void Encode(Encoder encoder);
    }

    public interface IEncodeInit
    {
        void Encode(Encoder encoder);
        bool Init(Decoder decoder);
    }

    // Class for decoding data from a memory buffer
    public class Decoder
    {
        static Decoder()
        {
            Assert.True(sizeof(Byte) == 1);
            Assert.True(sizeof(UInt16) == 2);
            Assert.True(sizeof(UInt32) == 4);
            Assert.True(sizeof(UInt64) == 8);
        }

        // Initialize decoder to decode from "buf"
        public Decoder(byte[] buf, int offset, int maxn)
        {
            Reset(buf, offset, maxn);
        }
        public void Reset(byte[] buf, int offset, int maxn)
        {
            Buffer = buf;
            Offset = offset;
            limit_ = maxn;
        }

        // Decoding routines.  Note that these do not check bounds
        public byte Get8()
        {
            byte v = Buffer[Offset];
            Offset += sizeof(byte);
            return v;
        }
        public UInt16 Get16()
        {
            UInt16 v = BitConverter.ToUInt16(Buffer, Offset);
            Offset += sizeof(UInt16);
            return v;
        }
        public UInt32 Get32()
        {
            UInt32 v = BitConverter.ToUInt32(Buffer, Offset);
            Offset += sizeof(UInt32);
            return v;
        }
        public UInt64 Get64()
        {
            UInt64 v = BitConverter.ToUInt64(Buffer, Offset);
            Offset += sizeof(UInt64);
            return v;
        }
        public float GetFloat()
        {
            float v = BitConverter.ToSingle(Buffer, Offset);
            Offset += sizeof(float);
            return v;
        }
        public double GetDouble()
        {
            double v = BitConverter.ToDouble(Buffer, Offset);
            Offset += sizeof(double);
            return v;
        }
        public void GetN<T>(T[] mem, int offset, int n)
        {
            System.Buffer.BlockCopy(Buffer, Offset, mem, offset, n);
            Offset += n;
        }
        public void GetPoints(S2Point[] points, int offset, int count)
        {
            var limit = offset + count;
            for (var i = offset; i < limit; i++)
            {
                points[i] = new S2Point(
                    GetDouble(),
                    GetDouble(),
                    GetDouble());
            }
        }
        public bool TryGetVarUInt32(out uint v)
        {
            if (Avail() < sizeof(int))
            {
                v = 0;
                return false;
            }
            var (result, newoffset) = VarintBitConverter.ToUInt32(Buffer, Offset);
            v = result;
            Offset = newoffset;
            return true;
        }
        public bool TryGetVarInt32(out int v)
        {
            if (Avail() < sizeof(int))
            {
                v = 0;
                return false;
            }
            var (result, newoffset) = VarintBitConverter.ToInt32(Buffer, Offset);
            v = result;
            Offset = newoffset;
            return true;
        }
        public bool TryGetVarUInt64(out ulong v)
        {
            var (result, newoffset) = VarintBitConverter.ToUInt64(Buffer, Offset);
            v = result;
            Offset = newoffset;
            return true;
        }
        public bool TryGetVarInt64(out long v)
        {
            var (result, newoffset) = VarintBitConverter.ToInt64(Buffer, Offset);
            v = result;
            Offset = newoffset;
            return true;
        }

        public void Skip(int n)
        {
            Offset += n;
        }
        // Return number of available bytes to read
        public int Avail()
        {
            Assert.True(limit_>= Offset);
            return limit_ - Offset;
        }

        public byte[] Buffer { get; private set; }
        public int Offset { get; private set; }
        private int limit_;
    }
}
