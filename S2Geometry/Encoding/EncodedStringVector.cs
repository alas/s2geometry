// This class allows an EncodedStringVector to be created by adding strings
// incrementally.  It also supports adding strings that are the output of
// another Encoder.  For example, to create a vector of encoded S2Polygons,
// you can do this:
//
// void EncodePolygons(S2Polygon[] polygons, Encoder encoder) {
//   StringVectorEncoder encoded_polygons;
//   foreach (var polygon : polygons) {
//     polygon.Encode(encoded_polygons.AddViaEncoder());
//   }
//   encoded_polygons.Encode(encoder);
// }

namespace S2Geometry;

public class StringVectorEncoder
{
    public StringVectorEncoder() => data_ = new Encoder();

    // Adds a string to the encoded vector.
    public void Add(string str)
    {
        offsets_.Add((ulong)data_.Length());
        data_.Ensure(str.Length);
        data_.PutN(System.Text.Encoding.ASCII.GetBytes(str));
    }

    // Adds a string to the encoded vector by means of the given Encoder.  The
    // string consists of all output added to the encoder before the next call
    // to any method of this class (after which the encoder is no longer valid).
    public Encoder AddViaEncoder()
    {
        offsets_.Add((ulong)data_.Length());
        return data_;
    }

    // Appends the EncodedStringVector representation to the given Encoder.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder)
    {
        offsets_.Add((ulong)data_.Length());
        // We don't encode the first element of "offsets_", which is always zero.
        var newarr = offsets_.Skip(1).Take(offsets_.Count - 1).ToArray();
        EncodedUintVector.EncodeUintVector(newarr, encoder);
        encoder.Ensure(data_.Length());
        encoder.PutEncoder(data_);
    }

    // Encodes a vector of strings in a format that can later be decoded as an
    // EncodedStringVector.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static void Encode(string[] v, Encoder encoder)
    {
        StringVectorEncoder string_vector = new();
        foreach (var str in v) string_vector.Add(str);
        string_vector.Encode(encoder);
    }

    // A vector consisting of the starting offset of each string in the
    // encoder's data buffer, plus a final entry pointing just past the end of
    // the last string.
    private readonly List<UInt64> offsets_ = new();
    private readonly Encoder data_;
}

// This class represents an encoded vector of strings.  Values are decoded
// only when they are accessed.  This allows for very fast initialization and
// no additional memory use beyond the encoded data.  The encoded data is not
// owned by this class; typically it points into a large contiguous buffer
// that contains other encoded data as well.
//
// This is one of several helper classes that allow complex data structures to
// be initialized from an encoded format inant time and then decoded on
// demand.  This can be a big performance advantage when only a small part of
// the data structure is actually used.
public class EncodedStringVector : IInitEncoder<EncodedStringVector>
{
    // Constructs an uninitialized object; requires Init() to be called.
    public EncodedStringVector() { }

    // Initializes the EncodedStringVector.  Returns false on errors, leaving
    // the vector in an unspecified state.
    //
    // REQUIRES: The Decoder data buffer must outlive this object.
    public static (bool, EncodedStringVector?) Init(Decoder decoder)
    {
        EncodedStringVector shape = new();
        if (!shape.offsets_.Init(decoder)) return (false, null);

        shape.data_ = decoder.Buffer;
        if (shape.offsets_.Count > 0)
        {
            var length = (int)shape.offsets_[^1];
            if (decoder.Avail() < length) return (false, null);
            decoder.Skip(length);
        }

        return (true, shape);
    }

    // Resets the vector to be empty.
    public void Clear()
    {
        offsets_.Clear();
        data_ = null;
    }

    // Returns the size of the original vector.
    public int Size()
    {
        return offsets_.Count;
    }

    // Returns the string at the given index.
    public string this[int i]
    {
        get
        {
            var start = (i == 0) ? 0 : (int)offsets_[i - 1];
            var limit = (int)offsets_[i];
            var countBytes = limit - start;
            var count = countBytes / sizeof(char);
            var buff = new char[count];
            Buffer.BlockCopy(data_, start, buff, 0, countBytes);
            return new string(buff);
        }
    }

    // Returns a Decoder initialized with the string at the given index.
    public Decoder GetDecoder(int i)
    {
        var start = (i == 0) ? 0 : (int)offsets_[i - 1];
        var limit = (int)offsets_[i];
        return new Decoder(data_, start, limit - start);
    }

    // Returns a pointer to the start of the string at the given index.  This is
    // faster than operator[] but returns an unbounded string.
    public (byte[] buffer, int offset) GetStart(int i)
    {
        UInt64 start = (i == 0) ? 0 : offsets_[i - 1];
        return (data_, (int)start);
    }

    // Returns the entire vector of original strings.  Requires that the
    // data buffer passed to the constructor persists until the result vector is
    // no longer needed.
    public string[] Decode()
    {
        int n = Size();
        var result = new string[n];
        for (int i = 0; i < n; ++i)
        {
            result[i] = this[i];
        }
        return result;
    }

    // The encoding must be identical to StringVectorEncoder.Encode().
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        offsets_.Encode(encoder);

        if (offsets_.Count > 0)
        {
            var length = (int)offsets_[^1];
            encoder.Ensure(length);
            encoder.PutN(data_, length);
        }
    }

    private readonly EncodedUintVector_UInt64 offsets_ = new();
    private byte[]? data_;
}
