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
            // The S2.Ortho reference direction points approximately due west.
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
            // The S2.Ortho reference direction points approximately due west.
            // Containment is determined by the unmatched edge immediately clockwise.
            S2ContainsVertexQuery q = new(MakePointOrDie("1:1"));
            q.AddEdge(MakePointOrDie("1:-5"), 1);
            q.AddEdge(MakePointOrDie("2:-4"), -1);
            q.AddEdge(MakePointOrDie("3:-3"), 1);
            q.AddEdge(MakePointOrDie("1:-5"), -1);
            Assert.Equal(-1, q.ContainsSign());
        }

        [Fact]
        public void Test_S2ContainsVertexQuery_MatchesLoopContainment()
        {
            // Check that the containment function defined is compatible with S2Loop
            // (which at least currently does not use this class).
            var loop = S2Loop.MakeRegularLoop(MakePointOrDie("89:-179"),
                                                S1Angle.FromDegrees(10), 1000);
            for (int i = 1; i <= loop.NumVertices; ++i)
            {
                S2ContainsVertexQuery q = new(loop.Vertex(i));
                q.AddEdge(loop.Vertex(i - 1), -1);
                q.AddEdge(loop.Vertex(i + 1), 1);
                Assert.Equal(q.ContainsSign() > 0, loop.Contains(loop.Vertex(i)));
            }
        }
    }
}
