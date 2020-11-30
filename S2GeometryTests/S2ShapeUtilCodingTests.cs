using Xunit;

namespace S2Geometry
{
    // Note that the "real" testing of these methods is in s2loop_measures_test
    // and s2polyline_measures_test.  This file only checks the handling of shapes
    // of different dimensions and shapes with multiple edge chains.
    public class S2ShapeUtilCodingTests
    {
        [Fact]
        public void Test_FastEncodeShape_S2Polygon()
        {
            var polygon = S2TextFormat.MakePolygonOrDie("0:0, 0:1, 1:0");
            var shape = new S2Polygon.Shape(polygon);
            Encoder encoder = new();
            S2ShapeUtilCoding.FastEncodeShape(shape, encoder);
            Decoder decoder = new(encoder.Buffer, 0, encoder.Length);
            var shape2 = S2ShapeUtilCoding.FullDecodeShape(S2Shape.TypeTag.S2Polygon, decoder);
            Assert.True(shape.Polygon == ((S2Polygon.Shape)shape2).Polygon);
        }

        [Fact]
        public void Test_FastEncodeTaggedShapes_MixedShapes()
        {
            // Tests encoding/decoding a collection of points, polylines, and polygons.
            var index = S2TextFormat.MakeIndexOrDie(
                "0:0 | 0:1 # 1:1, 1:2, 1:3 # 2:2; 2:3, 2:4, 3:3");
            Encoder encoder = new();
            S2ShapeUtilCoding.FastEncodeTaggedShapes(index, encoder);
            index.Encode(encoder);
            Decoder decoder = new(encoder.Buffer, 0, encoder.Length);
            MutableS2ShapeIndex decoded_index = new();
            decoded_index.Init(decoder, S2ShapeUtilCoding.FullDecodeShapeFactory(decoder));
            Assert.Equal(index.ToDebugString(), decoded_index.ToDebugString());
        }
    }
}
