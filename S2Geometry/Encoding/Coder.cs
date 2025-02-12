namespace S2Geometry;

// Class for encoding data into a memory buffer
public class Encoder(byte[] buffer, int Limit) : IEquatable<Encoder>
{
    // Support for variable length encoding with 7 bits per byte
    // (these are just simple wrappers around the Varint module)
    public const int kVarintMax32 = 5;
    public const int kVarintMax64 = 10;

    private byte[] Buffer = buffer;
    private int Offset;

    static Encoder()
    {
        MyDebug.Assert(sizeof(byte) == 1);
        MyDebug.Assert(sizeof(ushort) == 2);
        MyDebug.Assert(sizeof(uint) == 4);
        MyDebug.Assert(sizeof(ulong) == 8);
    }

    // Initialize encoder to encode into "buf"
    public Encoder() : this([], 0) { }

    #region Put

    // Encoding routines.  Note that these do not check bounds
    public void Put8(sbyte v) => Put8((byte)v);
    public void Put8(byte v)
    {
        MyDebug.Assert(Avail() >= sizeof(byte));
        Buffer[Offset] = v;
        Offset += sizeof(byte);
    }
    public void Put16(short v) => Put16((ushort)v);
    public void Put16(ushort v)
    {
        MyDebug.Assert(Avail() >= sizeof(ushort));
        var bytes = BitConverter.GetBytes(v);
        Buffer[Offset] = bytes[0];
        Buffer[++Offset] = bytes[1];
    }
    public void Put32(Int32 v) => Put32((uint)v);
    public void Put32(uint v)
    {
        MyDebug.Assert(Avail() >= sizeof(uint));
        var bytes = BitConverter.GetBytes(v);
        Buffer[Offset] = bytes[0];
        Buffer[Offset + 1] = bytes[1];
        Buffer[Offset + 2] = bytes[2];
        Buffer[Offset + 3] = bytes[3];
        Offset += 4;
    }
    public void Put64(long v) => Put64((ulong)v);
    public void Put64(ulong v)
    {
        MyDebug.Assert(Avail() >= sizeof(ulong));
        var bytes = BitConverter.GetBytes(v);
        Buffer[Offset] = bytes[0];
        Buffer[Offset + 1] = bytes[1];
        Buffer[Offset + 2] = bytes[2];
        Buffer[Offset + 3] = bytes[3];
        Buffer[Offset + 4] = bytes[4];
        Buffer[Offset + 5] = bytes[5];
        Buffer[Offset + 6] = bytes[6];
        Buffer[Offset + 7] = bytes[7];
        Offset += 8;
    }
    public void PutN(byte[] mem) => PutN(mem, mem.Length);
    public void PutN(byte[] mem, int n)
    {
        System.Buffer.BlockCopy(mem, 0, Buffer, Offset, n);
        Offset += n;
    }
    public void PutN(byte[] mem, int offset, int n)
    {
        System.Buffer.BlockCopy(mem, offset, Buffer, Offset, n);
        Offset += n;
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
        MyDebug.Assert(Avail() >= sizeof(float));
        var bytes = BitConverter.GetBytes(f);
        Buffer[Offset] = bytes[0];
        Buffer[Offset + 1] = bytes[1];
        Buffer[Offset + 2] = bytes[2];
        Buffer[Offset + 3] = bytes[3];
        Offset += 4;
    }
    public void PutDouble(double d)
    {
        MyDebug.Assert(Avail() >= sizeof(double));
        var bytes = BitConverter.GetBytes(d);
        Buffer[Offset] = bytes[0];
        Buffer[Offset + 1] = bytes[1];
        Buffer[Offset + 2] = bytes[2];
        Buffer[Offset + 3] = bytes[3];
        Buffer[Offset + 4] = bytes[4];
        Buffer[Offset + 5] = bytes[5];
        Buffer[Offset + 6] = bytes[6];
        Buffer[Offset + 7] = bytes[7];
        Offset += 8;
    }
    public void PutEncoder(Encoder encoder)
    {
        var n = encoder.Length();
        System.Buffer.BlockCopy(encoder.Buffer, 0, Buffer, Offset, n);
        Offset += n;
    }
    public void PutVarInt32(int v) => VarintBitConverter.GetVarintBytes(v, Buffer, ref Offset);
    public void PutVarUInt32(uint v) => VarintBitConverter.GetVarintBytes(v, Buffer, ref Offset);
    public void PutVarInt64(long v) => VarintBitConverter.GetVarintBytes(v, Buffer, ref Offset);
    public void PutVarUInt64(ulong v) => VarintBitConverter.GetVarintBytes(v, Buffer, ref Offset);

    #endregion

    public void Clear() => Offset = 0;

    // Return number of bytes encoded so far
    public int Length()
    {
        MyDebug.Assert(Offset >= 0);
        MyDebug.Assert(Offset <= Limit, "Catch the buffer overflow.");
        return Offset;
    }

    // Return number of bytes of space remaining in buffer
    public int Avail()
    {
        MyDebug.Assert(Limit >= Offset);
        return Limit - Offset;
    }

    // Return capacity of buffer.
    public int Capacity() => Limit - Offset;

    // REQUIRES: Encoder was created with the 0-argument constructor or 0-argument
    // reset().
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

    // Advances the write pointer by "N" bytes. It returns the position of the
    // pointer before the skip (in other words start of the skipped bytes).
    public void Skip(int N) => Offset += N;

    // REQUIRES: length() >= N
    // Removes the last N bytes out of the encoded buffer
    public void RemoveLast(int N)
    {
        MyDebug.Assert(Length() >= N);
        Offset -= N;
    }

    public byte Last() => Buffer[Offset - 1];

    private void EnsureSlowPath(int N)
    {
        MyDebug.Assert(Avail() < N);

        // Double buffer size, but make sure we always have at least N extra bytes
        int current_len = Length();
        int new_capacity = Math.Max(current_len + N, 2 * current_len);

        var new_buffer = new byte[new_capacity];
        if (Buffer is not null)
        {
            System.Buffer.BlockCopy(Buffer, 0, new_buffer, 0, current_len);
        }
        Buffer = new_buffer;

        Limit = new_capacity;
        MyDebug.Assert(Avail() >= N);
    }

    public string HexString() => Convert.ToHexString(Buffer, 0, Length());

    public Decoder GetDecoder(int truncate = 0) => new (Buffer, 0, Length() - truncate);

    #region IEquatable

    public override bool Equals(object? obj)
    {
        return obj is Encoder other && Equals(other);
    }

    public bool Equals(Encoder? other)
    {
        if (other is null) return false;

        return Buffer.SequenceEqual(other.Buffer);
    }

    public static bool operator ==(Encoder x, Encoder y) => x.Equals(y);
    public static bool operator !=(Encoder x, Encoder y) => !x.Equals(y);

    public override int GetHashCode()
    {
        return Buffer.GetHashCode();
    }

    #endregion
}

// Class for decoding data from a memory buffer
public class Decoder(byte[] buffer, int offset, int limit)
{
    public byte[] Buffer { get; private set; } = buffer;
    public int Offset { get; private set; } = offset;
    public int Limit { get; private set; } = limit;

    static Decoder()
    {
        MyDebug.Assert(sizeof(byte) == 1);
        MyDebug.Assert(sizeof(ushort) == 2);
        MyDebug.Assert(sizeof(uint) == 4);
        MyDebug.Assert(sizeof(ulong) == 8);
    }

    #region Get

    // Decoding routines.  Note that these do not check bounds
    public byte Peek8()
    {
        byte v = Buffer[Offset];
        return v;
    }
    public byte Get8()
    {
        byte v = Buffer[Offset];
        Offset += sizeof(byte);
        return v;
    }
    public ushort Get16()
    {
        ushort v = BitConverter.ToUInt16(Buffer, Offset);
        Offset += sizeof(ushort);
        return v;
    }
    public uint Get32()
    {
        uint v = BitConverter.ToUInt32(Buffer, Offset);
        Offset += sizeof(uint);
        return v;
    }
    public ulong Get64()
    {
        ulong v = BitConverter.ToUInt64(Buffer, Offset);
        Offset += sizeof(ulong);
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

    #endregion

    public void Skip(int n)
    {
        Offset += n;
    }

    // Return number of available bytes to read
    public int Avail()
    {
        MyDebug.Assert(Limit >= Offset);
        return Limit - Offset;
    }
}
