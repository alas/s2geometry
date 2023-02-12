/*namespace S2Geometry;

internal class S2Coder_Testing<T>
{
    // Encodes a type T into a new Encoder instance.
    internal Encoder EncodeToEncoder(S2Coder<T> coder, T shape)
    {
        Encoder encoder = new();
        coder.Encode(encoder, shape);
        return encoder;
    }

    // Decodes a type T encoded into an Encoder instance.
    internal bool DecodeFromEncoded(S2Coder<T> coder, Encoder encoder, T value,
                       out S2Error error)
    {
        Decoder decoder = encoder.GetDecoder();
        return coder.Decode(decoder, value, error);
    }

    // Encodes then decodes a type T via S2Coder<T>, returning decoded object.
    internal T RoundTrip(S2Coder<T> coder, T shape, out S2Error error)
    {
        Encoder encoder = EncodeToEncoder(coder, shape);
        T val;
        Assert.True(DecodeFromEncoded(coder, encoder, val, out error));
        return val;
    }
}*/
