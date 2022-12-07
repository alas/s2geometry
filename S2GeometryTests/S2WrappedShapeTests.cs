namespace S2Geometry;

public class S2WrappedShapeTests
{
    [Fact]
    internal void Test_S2WrappedShape_Coverage()
    {
        // Tests that all the S2Shape methods are implemented.

        var shape = MakeLaxPolygonOrDie("0:0; 1:1, 1:2, 2:1");
        S2WrappedShape wrapped_shape = new(shape);
        S2ShapeUtil_Testing.ExpectEqual(wrapped_shape, shape);
    }
}
