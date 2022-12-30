// Note that the "real" testing of these methods is in s2loop_measures_test
// and s2polyline_measures_test.  This file only checks the handling of shapes
// of different dimensions and shapes with multiple edge chains.

using static S2Geometry.S2ShapeUtilCoding;

namespace S2Geometry;

public class S2ShapeUtilCodingTests
{
    [Fact]
    internal void Test_FastEncodeShape_S2Polygon()
    {
        var polygon = MakePolygonOrDie("0:0, 0:1, 1:0");
        var shape = new S2Polygon.Shape(polygon);
        Encoder encoder = new();
        Assert.True(S2ShapeUtilCoding.FastEncodeShape(shape, encoder));
        var decoder = encoder.Decoder();
        var shape2 = S2ShapeUtilCoding.FullDecodeShape(S2Shape.TypeTag.S2Polygon, decoder);
        Assert.True(shape.Polygon == ((S2Polygon.Shape)shape2).Polygon);
    }

    [Fact]
    internal void Test_FastEncodeTaggedShapes_MixedShapes()
    {
        // Tests encoding/decoding a collection of points, polylines, and polygons.
        var index = MakeIndexOrDie(
            "0:0 | 0:1 # 1:1, 1:2, 1:3 # 2:2; 2:3, 2:4, 3:3");
        Encoder encoder = new();
        Assert.True(S2ShapeUtilCoding.FastEncodeTaggedShapes(index, encoder));
        index.Encode(encoder);
        var decoder = encoder.Decoder();
        MutableS2ShapeIndex decoded_index = new();
        Assert.True(decoded_index.Init(decoder, S2ShapeUtilCoding.FullDecodeShapeFactory(decoder)));
        Assert.Equal(index.ToDebugString(), decoded_index.ToDebugString());
    }

    [Fact]
    internal void Test_DecodeTaggedShapes_DecodeFromByteString()
    {
        var index = MakeIndexOrDie(
            "0:0 | 0:1 # 1:1, 1:2, 1:3 # 2:2; 2:3, 2:4, 3:3");
        var polygon = MakePolygonOrDie("0:0, 0:4, 4:4, 4:0");
        var polyline = MakePolylineOrDie("1:1, 1:2, 1:3");
        index.Add(new S2LaxPolylineShape(polyline));
        index.Add(new S2LaxPolygonShape(polygon));
        var bytes = Convert.FromHexString(
            "2932007C00E4002E0192010310000000000000F03F000000000000000000000000000000" +      
            "008AAFF597C0FEEF3F1EDD892B0BDF913F00000000000000000418B4825F3C81FDEF3F27" +      
            "DCF7C958DE913F1EDD892B0BDF913FD44A8442C3F9EF3FCE5B5A6FA6DDA13F1EDD892B0B" +      
            "DF913FAE0218F586F3EF3F3C3F66D2BBCAAA3F1EDD892B0BDF913F05010220FC7FB8B805" +      
            "F6EF3F28516A6D8FDBA13F27DCF7C958DEA13F96E20626CAEFEF3F4BF8A48399C7AA3F27" +      
            "DCF7C958DEA13F96B6DB0611E7EF3FC0221C80C6D8B13F27DCF7C958DEA13FE2337CCA8F"+      
            "E9EF3F6C573C9B60C2AA3F0EC9EF48C7CBAA3F0C0001040418B4825F3C81FDEF3F27DCF7"+      
            "C958DE913F1EDD892B0BDF913FD44A8442C3F9EF3FCE5B5A6FA6DDA13F1EDD892B0BDF91"+      
            "3FAE0218F586F3EF3F3C3F66D2BBCAAA3F1EDD892B0BDF913F05010120000000000000F0"+      
            "3F00000000000000000000000000000000F6FF70710BECEF3F28516A6D8FDBB13F000000"+      
            "00000000003C4A985423D8EF3F199E8D966CD0B13F28516A6D8FDBB13FF6FF70710BECEF"+      
            "3F000000000000000028516A6D8FDBB13F28C83900010003010403040504070400073807"+      
            "0E1B24292B3213000009030002130000110300092B00010001000000010D000002230410"+      
            "04020400020113082106110A4113000111030101");
        Decoder decoder=new(bytes, 0, bytes.Length);
        MutableS2ShapeIndex decoded_index=new();
        Assert.True(decoded_index.Init(decoder, S2ShapeUtilCoding.FullDecodeShapeFactory(decoder)));
        Assert.Equal(S2TextFormat.ToDebugString(index),
            S2TextFormat.ToDebugString(decoded_index));
    }

    [Fact]
    internal void Test_DecodeTaggedShapes_DecodeFromEncoded()
    {
        Encoder encoder=new();

        // Make an encoded shape.
        Assert.True(S2ShapeUtilCoding.FastEncodeShape(
            new S2PointVectorShape(ParsePointsOrDie("0:0, 0:1").ToArray()),
            encoder));
        var decoder = encoder.Decoder();
        var (success, encoded_shape) = EncodedS2PointVectorShape.Init(decoder);

        // Encode the encoded form.
        Encoder reencoder = new();
        Assert.True(S2ShapeUtilCoding.FastEncodeShape(encoded_shape!, reencoder));
        var encoded_decoder = reencoder.Decoder();

        // We can decode the shape in either full or lazy form from the same bytes.
        var full_shape = S2ShapeUtilCoding.FullDecodeShape(S2PointVectorShape.kTypeTag, encoded_decoder);
        Assert.Equal(full_shape!.GetTypeTag(), S2PointVectorShape.kTypeTag);

        encoded_decoder = reencoder.Decoder();
        var lazy_shape = S2ShapeUtilCoding.LazyDecodeShape(S2PointVectorShape.kTypeTag, encoded_decoder);
        Assert.Equal(lazy_shape!.GetTypeTag(), S2PointVectorShape.kTypeTag);
    }

    [Fact]
    internal void Test_SingletonShapeFactory_S2Polygon()
    {
        var polygon = MakePolygonOrDie("0:0, 0:1, 1:0");
        var shape = new S2Polygon.Shape(polygon);
        var shape_factory = SingletonShapeFactory(shape);

        // Returns the shape the first time.
        var shape2 = shape_factory[0];
        Assert.True(polygon.Equals(((S2Polygon.Shape)shape2).Polygon));

        // And nullptr after that.
        var shape3 = shape_factory[0];
        Assert.True(shape3 is null);
    }
}
