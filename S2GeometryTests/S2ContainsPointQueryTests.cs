namespace S2Geometry;

public class S2ContainsPointQueryTests
{
    [Fact]
    internal void Test_S2ContainsPointQuery_VertexModelOpen()
    {
        var index = MakeIndexOrDie("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6");
        S2ContainsPointQueryOptions options = new(S2VertexModel.OPEN);
        var q = index.MakeS2ContainsPointQuery(options);
        Assert.False(q.Contains(MakePointOrDie("0:0")));
        Assert.False(q.Contains(MakePointOrDie("-1:1")));
        Assert.False(q.Contains(MakePointOrDie("1:1")));
        Assert.False(q.Contains(MakePointOrDie("0:2")));
        Assert.False(q.Contains(MakePointOrDie("0:3")));
        Assert.False(q.Contains(MakePointOrDie("0:5")));
        Assert.False(q.Contains(MakePointOrDie("0:7")));
        Assert.False(q.Contains(MakePointOrDie("2:6")));
        Assert.True(q.Contains(MakePointOrDie("1:6")));
        Assert.False(q.Contains(MakePointOrDie("10:10")));

        // Test the last few cases using the Init() method instead.
        S2ContainsPointQuery<MutableS2ShapeIndex> q2 = new(index, options);
        Assert.False(q2.ShapeContains(index.Shape(1)!, MakePointOrDie("1:6")));
        Assert.True(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("1:6")));
        Assert.False(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("0:5")));
        Assert.False(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("0:7")));
    }

    [Fact]
    internal void Test_S2ContainsPointQuery_VertexModelSemiOpen()
    {
        var index = MakeIndexOrDie("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6");
        S2ContainsPointQueryOptions options = new(S2VertexModel.SEMI_OPEN);
        var q = index.MakeS2ContainsPointQuery(options);
        Assert.False(q.Contains(MakePointOrDie("0:0")));
        Assert.False(q.Contains(MakePointOrDie("-1:1")));
        Assert.False(q.Contains(MakePointOrDie("1:1")));
        Assert.False(q.Contains(MakePointOrDie("0:2")));
        Assert.False(q.Contains(MakePointOrDie("0:5")));
        Assert.True(q.Contains(MakePointOrDie("0:7")));  // Contained vertex.
        Assert.False(q.Contains(MakePointOrDie("2:6")));
        Assert.True(q.Contains(MakePointOrDie("1:6")));
        Assert.False(q.Contains(MakePointOrDie("10:10")));

        // Test the last few cases using the Init() method instead.
        S2ContainsPointQuery<MutableS2ShapeIndex> q2 = new(index, options);
        Assert.False(q2.ShapeContains(index.Shape(1)!, MakePointOrDie("1:6")));
        Assert.True(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("1:6")));
        Assert.False(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("0:5")));
        Assert.True(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("0:7")));
    }

    [Fact]
    internal void Test_S2ContainsPointQuery_VertexModelClosed()
    {
        var index = MakeIndexOrDie("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6");
        S2ContainsPointQueryOptions options = new(S2VertexModel.CLOSED);
        var q = index.MakeS2ContainsPointQuery(options);
        Assert.True(q.Contains(MakePointOrDie("0:0")));
        Assert.True(q.Contains(MakePointOrDie("-1:1")));
        Assert.True(q.Contains(MakePointOrDie("1:1")));
        Assert.False(q.Contains(MakePointOrDie("0:2")));
        Assert.True(q.Contains(MakePointOrDie("0:5")));
        Assert.True(q.Contains(MakePointOrDie("0:7")));
        Assert.True(q.Contains(MakePointOrDie("2:6")));
        Assert.True(q.Contains(MakePointOrDie("1:6")));
        Assert.False(q.Contains(MakePointOrDie("10:10")));

        // Test the last few cases using the Init() method instead.
        S2ContainsPointQuery<MutableS2ShapeIndex> q2 = new(index, options);
        Assert.False(q2.ShapeContains(index.Shape(1)!, MakePointOrDie("1:6")));
        Assert.True(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("1:6")));
        Assert.True(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("0:5")));
        Assert.True(q2.ShapeContains(index.Shape(2)!, MakePointOrDie("0:7")));
    }

    [Fact]
    internal void Test_S2ContainsPointQuery_GetContainingShapes()
    {
        // Also tests ShapeContains().
        int kNumVerticesPerLoop = 10;
        S1Angle kMaxLoopRadius = S2Testing.KmToAngle(10);
        S2Cap center_cap = new(S2Testing.RandomPoint(), kMaxLoopRadius);
        MutableS2ShapeIndex index = [];
        for (int i = 0; i < 100; ++i)
        {
            var loop = S2Loop.MakeRegularLoop(
                S2Testing.SamplePoint(center_cap),
                S2Testing.Random.RandDouble() * kMaxLoopRadius, kNumVerticesPerLoop);
            index.Add(new S2Loop.Shape(loop));
        }
        var query = index.MakeS2ContainsPointQuery();
        for (int i = 0; i < 100; ++i)
        {
            S2Point p = S2Testing.SamplePoint(center_cap);
            List<S2Shape> expected = [];
            foreach (var shape in index)
            {
                var loop = ((S2Loop.Shape)shape).Loop;
                if (loop.Contains(p))
                {
                    Assert.True(query.ShapeContains(shape, p));
                    expected.Add(shape);
                }
                else
                {
                    Assert.False(query.ShapeContains(shape, p));
                }
            }
            var actual = query.GetContainingShapes(p);
            Assert.Equal(expected, actual);
        }
    }

    [Fact]
    internal void Test_S2ContainsPointQuery_VisitIncidentEdges()
    {
        var index = MakeIndexOrDie("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2");
        ExpectIncidentEdgeIds([new(0, 0)], index, MakePointOrDie("0:0"));
        ExpectIncidentEdgeIds([new(0, 1), new(1, 0)], index, MakePointOrDie("1:1"));
        ExpectIncidentEdgeIds([new(1, 0), new(2, 0), new(2, 2)], index, MakePointOrDie("1:2"));
        ExpectIncidentEdgeIds([new(2, 0), new(2, 1)], index, MakePointOrDie("1:3"));
        ExpectIncidentEdgeIds([new(2, 1), new(2, 2)], index, MakePointOrDie("2:2"));
    }

    private static void ExpectIncidentEdgeIds(EdgeVector expected,
        MutableS2ShapeIndex index, S2Point p)
    {
        EdgeVector actual = [];
        var q = index.MakeS2ContainsPointQuery();
        Assert.True(
            q.VisitIncidentEdges(p, e =>
            {
                actual.Add(e.Id);
                return true;
            }));
        Assert.Equal(expected, actual);
    }
}
