using Xunit;
using static S2Geometry.S2TextFormat;

namespace S2Geometry
{
    public class S2ContainsVertexQueryTests
    {
        [Fact]
        public void Test_S2ContainsVertexQuery_Undetermined()
        {
            S2ContainsVertexQuery q = new(MakePointOrDie("1:2"));
            q.AddEdge(MakePointOrDie("3:4"), 1);
            q.AddEdge(MakePointOrDie("3:4"), -1);
            Assert.Equal(0, q.ContainsSign());
        }

        [Fact]
        public void Test_S2ContainsVertexQuery_ContainedWithDuplicates()
        {
            // The S2::RefDir reference direction points approximately due west.
            // Containment is determined by the unmatched edge immediately clockwise.
            S2ContainsVertexQuery q = new(MakePointOrDie("0:0"));
            q.AddEdge(MakePointOrDie("3:-3"), -1);
            q.AddEdge(MakePointOrDie("1:-5"), 1);
            q.AddEdge(MakePointOrDie("2:-4"), 1);
            q.AddEdge(MakePointOrDie("1:-5"), -1);
            Assert.Equal(1, q.ContainsSign());
        }

        [Fact]
        public void Test_S2ContainsVertexQuery_NotContainedWithDuplicates()
        {
            // The S2::RefDir reference direction points approximately due west.
            // Containment is determined by the unmatched edge immediately clockwise.
            S2ContainsVertexQuery q = new(MakePointOrDie("1:1"));
            q.AddEdge(MakePointOrDie("1:-5"), 1);
            q.AddEdge(MakePointOrDie("2:-4"), -1);
            q.AddEdge(MakePointOrDie("3:-3"), 1);
            q.AddEdge(MakePointOrDie("1:-5"), -1);
            Assert.Equal(-1, q.ContainsSign());
        }

        // Tests that S2ContainsVertexQuery is compatible with S2::AngleContainsVertex.
        [Fact]
        public void Test_S2ContainsVertexQuery_CompatibleWithAngleContainsVertex()
        {
            var points = S2Testing.MakeRegularPoints(MakePointOrDie("89:1"),
                S1Angle.FromDegrees(5), 10);
            S2PointLoopSpan loop=new(points);
            for (int i = 0; i < loop.Count; ++i)
            {
                S2Point a = loop[i];
                S2Point b = loop[i + 1];
                S2Point c = loop[i + 2];
                S2ContainsVertexQuery q=new(b);
                q.AddEdge(a, -1);
                q.AddEdge(c, 1);
                Assert.Equal(q.ContainsSign() > 0, S2.AngleContainsVertex(a, b, c));
            }
        }

        // Tests compatibility with S2::AngleContainsVertex() for a degenerate edge.
        [Fact]
        public void Test_S2ContainsVertexQuery_CompatibleWithAngleContainsVertexDegenerate()
        {
            S2Point a=new(1, 0, 0), b=new(0, 1, 0);
            S2ContainsVertexQuery q=new(b);
            q.AddEdge(a, -1);
            q.AddEdge(a, 1);
            Assert.Equal(q.ContainsSign() > 0, S2.AngleContainsVertex(a, b, a));
        }
    }
}
