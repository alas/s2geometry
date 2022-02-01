// Helper functions for encoding/decoding S2Shapes and S2Shape vectors.
// Design goals:
//
//  - Allow control over encoding tradeoffs (i.e., speed vs encoding size).
//
//  - Allow control over decoding tradeoffs (e.g., whether to decode data
//    immediately or lazily).
//
//  - Support client-defined S2Shape types.
//
//  - Support client-defined encodings of standard S2Shape types.
//
//  - Support custom encodings of shape vectors; e.g., if all shapes are of a
//    known type, then there is no need to tag them individually.

namespace S2Geometry;

public static class S2ShapeUtilCoding
{
    // A function that appends a serialized representation of the given shape to
    // the given Encoder.  The encoding should *not* include any type information
    // (e.g., shape.type_tag()); the caller is responsible for encoding this
    // separately if necessary.
    //
    // The purpose of ShapeEncoder is to allow implementing custom encodings
    // without subclassing the affected S2Shape types.  For example, you can
    // define an encoder that wraps one of the standard functions and adds
    // exceptions:
    //
    // void MyShapeEncoder(S2Shape shape, Encoder encoder) {
    //   if (shape.type_tag() == S2LaxPolygonShape::kTypeTag) {
    //     MyEncodeLaxPolygon(shape, encoder);
    //     return true;
    //   } else {
    //     return CompactEncodeShape(shape, encoder);
    //   }
    // }
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public delegate bool ShapeEncoder(S2Shape shape, Encoder encoder);

    // A ShapeEncoder that can encode all standard S2Shape types, preferring fast
    // (but larger) encodings.  For example, points are typically represented as
    // uncompressed S2Point values (24 bytes each).
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static bool FastEncodeShape(S2Shape shape, Encoder encoder)
    {
        var tag = shape.GetTypeTag();
        if (tag == S2Shape.TypeTag.kNoTypeTag)
        {
            System.Diagnostics.Debug.WriteLine("(DFATAL) Unsupported S2Shape type: " + tag);
            return false;
        }
        // Update the following constant when adding new S2Shape encodings.
        Assert.True(tag < S2Shape.TypeTag.kNextAvailableTypeTag);
        shape.Encode(encoder, CodingHint.FAST);
        return true;
    }

    // A ShapeEncoder that encode all standard S2Shape types, preferring
    // compact (but slower) encodings.  For example, points that have been snapped
    // to S2CellId centers will be encoded using at most 8 bytes.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static bool CompactEncodeShape(S2Shape shape, Encoder encoder)
    {
        var tag = shape.GetTypeTag();
        if (tag == S2Shape.TypeTag.kNoTypeTag)
        {
            System.Diagnostics.Debug.WriteLine("(DFATAL) Unsupported S2Shape type: " + tag);
            return false;
        }
        // Update the following constant when adding new S2Shape encodings.
        Assert.True(tag < S2Shape.TypeTag.kNextAvailableTypeTag);
        shape.Encode(encoder, CodingHint.COMPACT);
        return true;
    }

    // A function that decodes an S2Shape of the given type, consuming data from
    // the given Decoder.  Returns null on errors.
    public delegate S2Shape? ShapeDecoder(S2Shape.TypeTag tag, Decoder decoder);

    // A ShapeDecoder that fully decodes an S2Shape of the given type.  After this
    // function returns, the underlying Decoder data is no longer needed.
    public static S2Shape? FullDecodeShape(S2Shape.TypeTag tag, Decoder decoder)
    {
        Dictionary<S2Shape.TypeTag, Func<Decoder, (bool, S2Shape?)>> dic = new()
        {
            { S2Shape.TypeTag.S2Polygon, d => S2Polygon.OwningShape.Init(d) },
            { S2Shape.TypeTag.S2Polyline, d => S2Polyline.OwningShape.Init(d) },
            { S2Shape.TypeTag.S2PointVectorShape, d => S2PointVectorShape.Init(d) },
            { S2Shape.TypeTag.S2LaxPolylineShape, d => S2LaxPolylineShape.Init(d) },
            { S2Shape.TypeTag.S2LaxPolygonShape, d => S2LaxPolygonShape.Init(d) },
        };
        var speificInit = dic.GetValueOrDefault(tag);
        if (speificInit is not null)
        {
            var (success, shape) = speificInit(decoder);
            if (success)
                return shape;
            else
                return null;
        }

        // "Unsupported S2Shape type: " + tag;
        return null;
    }

    // A ShapeDecoder that prefers to decode the given S2Shape lazily (as data is
    // accessed).  This is only possible when the given shape type (e.g.,
    // LaxPolygonShape) has an alternate implementation that can work directly
    // with encoded data (e.g., EncodedLaxPolygonShape).  All other shape types
    // are handled by decoding them fully (e.g., S2Polygon.Shape).
    public static S2Shape? LazyDecodeShape(S2Shape.TypeTag tag, Decoder decoder)
    {
        Dictionary<S2Shape.TypeTag, Func<Decoder, (bool, S2Shape?)>> dic = new()
        {
            { S2Shape.TypeTag.S2PointVectorShape, d => EncodedS2PointVectorShape.Init(d) },
            { S2Shape.TypeTag.S2LaxPolylineShape, d => EncodedS2LaxPolylineShape.Init(d) },
            { S2Shape.TypeTag.S2LaxPolygonShape, d => EncodedS2LaxPolygonShape.Init(d) },
        };
        var speificInit = dic.GetValueOrDefault(tag);
        if (speificInit is not null)
        {
            var (success, shape) = speificInit(decoder);
            if (success)
                return shape;
            else
                return null;
        }
        
        return FullDecodeShape(tag, decoder);
    }

    // Encodes the shapes in the given S2ShapeIndex.  Each shape is encoded with a
    // type tag allows it to be decoded into an S2Shape of the appropriate type.
    // "shape_encoder" allows control over the encoding strategy.  Note that when
    // an S2ShapeIndex is also being encoded, it should be encoded *after* the
    // shape vector, like this:
    //
    //   S2ShapeUtil.CompactEncodeTaggedShapes(index, encoder);
    //   index.Encode(encoder);
    //
    // This is because when the index is decoded, the shape vector is required as
    // a parameter.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static bool EncodeTaggedShapes(S2ShapeIndex index, ShapeEncoder shape_encoder, Encoder encoder)
    {
        var shape_vector = new StringVectorEncoder();
        foreach (var shape in index)
        {
            var sub_encoder = shape_vector.AddViaEncoder();
            if (shape == null) continue;  // Encode as zero bytes.

            sub_encoder.Ensure(Encoder.kVarintMax32);
            sub_encoder.PutVarInt32((int)shape.GetTypeTag());
            if (!shape_encoder(shape, sub_encoder)) return false;
        }
        shape_vector.Encode(encoder);
        return true;
    }

    // Convenience function that calls EncodeTaggedShapes using FastEncodeShape as
    // the ShapeEncoder.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static bool FastEncodeTaggedShapes(S2ShapeIndex index, Encoder encoder)
    {
        return EncodeTaggedShapes(index, FastEncodeShape, encoder);
    }

    // Convenience function that calls EncodeTaggedShapes using CompactEncodeShape
    // as the ShapeEncoder.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static bool CompactEncodeTaggedShapes(S2ShapeIndex index, Encoder encoder)
    {
        return EncodeTaggedShapes(index, CompactEncodeShape, encoder);
    }

    // A ShapeFactory that decodes a vector generated by EncodeTaggedShapes()
    // above.  Example usage:
    //
    //   index.Init(decoder, S2ShapeUtil.FullDecodeShapeFactory(decoder));
    //
    // Note that the S2Shape vector must be encoded *before* the S2ShapeIndex
    // (like the example code for EncodeTaggedShapes), since the shapes need to be
    // decoded before the index.
    //
    // REQUIRES: The Decoder data buffer must outlive all calls to the given
    //           ShapeFactory (not including its destructor).
    public class TaggedShapeFactory : S2ShapeIndex.ShapeFactory
    {
        // Returns an empty vector and/or null S2Shapes on decoding errors.
        public TaggedShapeFactory(ShapeDecoder shape_decoder, Decoder decoder)
        {
            decoder_ = decoder;
            shape_decoder_ = shape_decoder;
            var (success, shape) = EncodedStringVector.Init(decoder);
            if (success)
            {
                encoded_shapes_ = shape!;
                encoded_shapes_.Clear();
            }
        }

        public override int Count => encoded_shapes_?.Size() ?? 0;
        public override S2Shape? this[int shape_id]
        {
            get
            {
                var decoder = encoded_shapes_?.GetDecoder(shape_id);
                return decoder is null || !decoder.TryGetVarUInt32(out var tag) ? null
                    : shape_decoder_((S2Shape.TypeTag)tag, decoder);
            }
        }

        public override object CustomClone() => new TaggedShapeFactory(shape_decoder_, decoder_);

        private readonly ShapeDecoder shape_decoder_;
        private readonly EncodedStringVector? encoded_shapes_;
        private readonly Decoder decoder_;
    }

    // Convenience function that calls TaggedShapeFactory using FullDecodeShape
    // as the ShapeDecoder.
    public static TaggedShapeFactory FullDecodeShapeFactory(Decoder decoder)
    {
        return new TaggedShapeFactory(FullDecodeShape, decoder);
    }

    // Convenience function that calls TaggedShapeFactory using LazyDecodeShape
    // as the ShapeDecoder.
    public static TaggedShapeFactory LazyDecodeShapeFactory(Decoder decoder)
    {
        return new TaggedShapeFactory(LazyDecodeShape, decoder);
    }

    // A ShapeFactory that simply returns shapes from the given vector.
    //
    // REQUIRES: Each shape is requested at most once.  (This implies that when
    // the ShapeFactory is passed to an S2ShapeIndex, S2ShapeIndex.Minimize must
    // not be called.)  Additional requests for the same shape return null.
    public class VectorShapeFactory : S2ShapeIndex.ShapeFactory
    {
        public VectorShapeFactory(List<S2Shape> shapes)
        {
            shared_shapes_ = shapes;
        }

        public override int Count => shared_shapes_.Count;
        public override S2Shape this[int shape_id]
        {
            get
            {
                return shared_shapes_[shape_id];
            }
        }

        public override object CustomClone() => new VectorShapeFactory(shared_shapes_);

        // Since this class is copyable, we need to access the shape vector through
        // a shared pointer.
        private readonly List<S2Shape> shared_shapes_;
    }

    // A ShapeFactory that returns the single given S2Shape.  Useful for testing.
    public static VectorShapeFactory SingletonShapeFactory(S2Shape shape)
    {
        return new VectorShapeFactory(new List<S2Shape> { shape });
    }

    // A ShapeFactory that wraps the shapes from the given index.  Used for testing.
    public class WrappedShapeFactory : S2ShapeIndex.ShapeFactory
    {
        // REQUIRES: The given index must persist for the lifetime of this object.
        public WrappedShapeFactory(S2ShapeIndex index) { index_ = index; }

        public override int Count => index_.NumShapeIds();
        public override S2Shape this[int shape_id]
        {
            get
            {
                var shape = index_.Shape(shape_id);
                if (shape == null) return null;
                return new S2WrappedShape(shape);
            }
        }

        public override object CustomClone() => new WrappedShapeFactory(index_);

        private readonly S2ShapeIndex index_;
    }

    // Encodes the shapes in the given index, which must all have the same type.
    // The given Shape type does *not* need to have a "type tag" assigned.
    // This is useful for encoding experimental or locally defined types, or when
    // the S2Shape type in a given index is known in advance.
    //
    // REQUIRES: "Shape" must have an Encode(Encoder*, s2coding::CodingHint) method.
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public static void EncodeHomogeneousShapes(S2ShapeIndex index, Encoder encoder,
        CodingHint hint = CodingHint.COMPACT)
    {
        var shape_vector = new StringVectorEncoder();
        foreach (var shape in index)
        {
            Assert.True(shape != null);
            (shape as IEncoder)!.Encode(shape_vector.AddViaEncoder(), hint);
        }
        shape_vector.Encode(encoder);
    }

    // A ShapeFactory that decodes shapes of a given fixed type (e.g.,
    // EncodedS2LaxPolylineShape).  Example usage:
    //
    // REQUIRES: The Shape type must have an Init(Decoder) method.
    // REQUIRES: The Decoder data buffer must outlive the returned ShapeFactory
    //           and all shapes returned by that factory.
    public class HomogeneousShapeFactory<TShape> : S2ShapeIndex.ShapeFactory
        where TShape : S2Shape, IInitEncoder<TShape>, new()
    {
        // Returns an empty vector and/or null S2Shapes on decoding errors.
        public HomogeneousShapeFactory(Decoder decoder)
        {
            decoder_ = decoder;
            var (success, shape) = EncodedStringVector.Init(decoder);
            if (!success)
            {
                encoded_shapes_ = shape!;
                encoded_shapes_.Clear();
            }
        }

        public override int Count => encoded_shapes_?.Size() ?? 0;
        public override S2Shape? this[int shape_id]
        {
            get
            {
                var decoder = encoded_shapes_?.GetDecoder(shape_id);
                if (decoder is not null)
                {
                    var (success, shape) = TShape.Init(decoder);
                    if (success)
                        return shape;
                }
                return null;
            }
        }

        public override object CustomClone() => new HomogeneousShapeFactory<TShape>(decoder_);

        private readonly EncodedStringVector? encoded_shapes_;
        private readonly Decoder decoder_;
    }
}

