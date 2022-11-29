namespace S2Geometry;

using S2BuilderUtil;
using static S2Builder;
using static S2Builder.GraphOptions;
public class S2BuilderUtil_FindPolygonDegeneraciesTests
{
    [Fact]
    public void Test_FindPolygonDegeneracies_EmptyPolygon() {
        ExpectDegeneracies("", Array.Empty<TestDegeneracy>());
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_NoDegeneracies() {
        ExpectDegeneracies("0:0, 0:1, 1:0", Array.Empty<TestDegeneracy>());
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_PointShell() {
        ExpectDegeneracies("0:0", new[] { new TestDegeneracy("0:0, 0:0", false) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_SiblingPairShells() {
        ExpectDegeneracies("0:0, 0:1, 1:0; 1:0, 0:1, 0:0", new[] {
            new TestDegeneracy("0:0, 0:1", false), new TestDegeneracy("0:1, 0:0", false),
            new TestDegeneracy("0:1, 1:0", false), new TestDegeneracy("1:0, 0:1", false),
            new TestDegeneracy("0:0, 1:0", false), new TestDegeneracy("1:0, 0:0", false),
        });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_AttachedSiblingPairShells() {
        ExpectDegeneracies("0:0, 0:1, 1:0; 1:0, 2:0", new[] {
            new TestDegeneracy("1:0, 2:0", false), new TestDegeneracy("2:0, 1:0", false) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_AttachedSiblingPairHoles() {
        ExpectDegeneracies("0:0, 0:3, 3:0; 0:0, 1:1", new[] {
                 new TestDegeneracy("0:0, 1:1", true), new TestDegeneracy("1:1, 0:0", true) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_AttachedSiblingPairShellsAndHoles() {
        ExpectDegeneracies("0:0, 0:3, 3:0; 3:0, 1:1; 3:0, 5:5", new[] {
            new TestDegeneracy("3:0, 1:1", true), new TestDegeneracy("1:1, 3:0", true),
            new TestDegeneracy("3:0, 5:5", false), new TestDegeneracy("5:5, 3:0", false) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_DegenerateShellsOutsideLoop() {
        ExpectDegeneracies("0:0, 0:3, 3:3, 3:0; 4:4, 5:5; 6:6", new[] {
            new TestDegeneracy("4:4, 5:5", false), new TestDegeneracy("5:5, 4:4", false),
            new TestDegeneracy("6:6, 6:6", false) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_DegenerateHolesWithinLoop() {
        ExpectDegeneracies("0:0, 0:5, 5:5, 5:0; 1:1, 2:2; 3:3", new[]{
            new TestDegeneracy("1:1, 2:2", true), new TestDegeneracy("2:2, 1:1", true),
            new TestDegeneracy("3:3, 3:3", true) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_PointHoleWithinFull() {
        ExpectDegeneracies("full; 0:0", new[] { new TestDegeneracy("0:0, 0:0", true) });
    }

    [Fact]
    public void Test_FindPolygonDegeneracies_SiblingPairHolesWithinFull() {
        ExpectDegeneracies("full; 0:0, 0:1, 1:0; 1:0, 0:1, 0:0", new []{
            new TestDegeneracy("0:0, 0:1", true), new TestDegeneracy("0:1, 0:0", true),
            new TestDegeneracy("0:1, 1:0", true), new TestDegeneracy("1:0, 0:1", true),
            new TestDegeneracy("0:0, 1:0", true), new TestDegeneracy("1:0, 0:0", true)
        });
    }

    private struct TestDegeneracy : IEquatable<TestDegeneracy>, IComparable<TestDegeneracy>
    {
        public string EdgeStr;
        public bool IsHole;

        public TestDegeneracy(string _edge_str, bool _is_hole)
        { EdgeStr = _edge_str; IsHole = _is_hole; }

        public bool Equals(TestDegeneracy other) => EdgeStr == other.EdgeStr && IsHole == other.IsHole;
        public override bool Equals(object? obj) => obj is TestDegeneracy td && Equals(td);
        public override int GetHashCode() => HashCode.Combine(EdgeStr, IsHole);

        public static bool operator ==(TestDegeneracy x, TestDegeneracy y) => Equals(x, y);
        public static bool operator !=(TestDegeneracy left, TestDegeneracy right) => !Equals(left, right);

        public int CompareTo(TestDegeneracy other)
        {
            if (EdgeStr.CompareTo(other.EdgeStr) != 0)
                return EdgeStr.CompareTo(other.EdgeStr);

            return IsHole.CompareTo(other.IsHole);
        }
        public static bool operator <(TestDegeneracy x, TestDegeneracy y)
        {
            return x.CompareTo(y) < 0;
        }
        public static bool operator >(TestDegeneracy x, TestDegeneracy y)
        {
            return x.CompareTo(y) > 0;
        }

        public override string ToString() => $"{(IsHole ? "Hole(" : "Shell(")}{EdgeStr}) ";
    }

    private class DegeneracyCheckingLayer : Layer
    {
        public DegeneracyCheckingLayer(TestDegeneracy[] expected)
        { expected_ = expected; }
        public override GraphOptions GraphOptions_() {
            return new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
                                DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
        }
        public override void Build(Graph g, out S2Error error)
        {
            var degeneracies = PolygonDegeneracy.FindPolygonDegeneracies(g, out error);
            // Convert the output into a human-readable format.
            List<TestDegeneracy> actual = new();
            foreach (var degeneracy in degeneracies)
            {
                var edge = g.GetEdge((int)degeneracy.EdgeId);
                var points = new S2Point[] { g.Vertex(edge.ShapeId), g.Vertex(edge.EdgeId) };
                actual.Add(new TestDegeneracy(S2TextFormat.ToDebugString(points), degeneracy.IsHole));
            }
            actual = new SortedSet<TestDegeneracy>(actual).ToList();
            Assert.True(expected_.SequenceEqual(actual));
            Assert.Equal(PolygonDegeneracy.IsFullyDegenerate(g), degeneracies.Count == g.NumEdges);
        }

        private readonly TestDegeneracy[] expected_;
    }

    private static void ExpectDegeneracies(string polygon_str, TestDegeneracy[] expected)
    {
        S2Builder builder = new(new Options());
        builder.StartLayer(new DegeneracyCheckingLayer(expected));
        var polygon = MakeLaxPolygonOrDie(polygon_str);
        builder.AddIsFullPolygonPredicate((IsFullPolygonPredicate)((Graph graph, out S2Error error) => {
            error = S2Error.OK;
            return (bool)polygon.GetReferencePoint().Contained;
        }));
        builder.AddShape(polygon);
        Assert.True(builder.Build(out var error));
    }
}
