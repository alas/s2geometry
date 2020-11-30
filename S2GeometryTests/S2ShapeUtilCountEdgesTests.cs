using S2Geometry.S2ShapeUtil;
using Xunit;

namespace S2Geometry
{
    public class S2ShapeUtilCountEdgesTests
    {
        [Fact]
        public void Test_CountEdgesUpTo_StopsEarly()
        {
            var index = S2TextFormat.MakeIndexOrDie(
                "0:0 | 0:1 | 0:2 | 0:3 | 0:4 # 1:0, 1:1 | 1:2, 1:3 | 1:4, 1:5, 1:6 #");
            // Verify the test parameters.
            Assert.Equal(4, index.NumShapeIds());
            Assert.Equal(5, index.Shape(0).NumEdges);
            Assert.Equal(1, index.Shape(1).NumEdges);
            Assert.Equal(1, index.Shape(2).NumEdges);
            Assert.Equal(2, index.Shape(3).NumEdges);

            Assert.Equal(9, index.GetCountEdges());
            Assert.Equal(5, index.GetCountEdgesUpTo(1));
            Assert.Equal(5, index.GetCountEdgesUpTo(5));
            Assert.Equal(6, index.GetCountEdgesUpTo(6));
            Assert.Equal(9, index.GetCountEdgesUpTo(8));
        }
    }
}
