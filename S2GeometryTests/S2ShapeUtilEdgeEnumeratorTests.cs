namespace S2Geometry;

using static S2ShapeUtil;

public class S2ShapeUtilEdgeEnumeratorTests
{
    [Fact]
    public void Test_S2ShapeUtilEdgeEnumeratorTest_Empty() =>
        Verify(MakeIndexOrDie("##"));

    [Fact]
    public void Test_S2ShapeutilEdgeEnumeratorTest_Points() =>
        Verify(MakeIndexOrDie("0:0|1:1##"));

    [Fact]
    public void Test_S2ShapeutilEdgeEnumeratorTest_Lines() =>
        Verify(MakeIndexOrDie("#0:0,10:10|5:5,5:10|1:2,2:1#"));

    [Fact]
    public void Test_S2ShapeutilEdgeEnumeratorTest_Polygons() =>
        Verify(MakeIndexOrDie("##10:10,10:0,0:0|-10:-10,-10:0,0:0,0:-10"));


    [Fact]
    public void Test_S2ShapeutilEdgeEnumeratorTest_Collection() =>
        Verify(MakeIndexOrDie("1:1|7:2#1:1,2:2,3:3|2:2,1:7#" +
            "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0"));

    [Fact]
    public void Test_S2ShapeutilEdgeEnumeratorTest_Remove()
    {
        var index = MakeIndexOrDie("1:1|7:2#1:1,2:2,3:3|2:2,1:7#" +
            "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
        index.Release(0);
        Verify(index);
    }


    [Fact]
    public void Test_S2ShapeutilEdgeEnumeratorTest_AssignmentAndEquality()
    {
        var index1 = MakeIndexOrDie("1:1|7:2#1:1,2:2,3:3|2:2,1:7#" +
            "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
        var index2 = MakeIndexOrDie("1:1|7:2#1:1,2:2,3:3|2:2,1:7#" +
            "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
        var it1 = new EdgeEnumerator(index1);
        var it2 = new EdgeEnumerator(index2);

        // Different indices.
        Assert.True(it1 != it2);

        it1 = (EdgeEnumerator)it2.CustomClone();
        Assert.Equal(it1, it2);

        it1.MoveNext();
        Assert.True(it1 != it2);

        it2.MoveNext();
        Assert.Equal(it1, it2);
    }

    // Returns the full list of edges in g.
    // The edges are collected from points, lines, and polygons in that order.
    private static List<S2Shape.Edge> GetEdges(S2ShapeIndex index)
    {
        var result = new List<S2Shape.Edge>();
        foreach (var shape in index)
        {
            if (shape == null) continue;

            for (var j = 0; j < shape.NumEdges(); ++j)
            {
                result.Add(shape.GetEdge(j));
            }
        }
        return result;
    }

    // Verifies that the edges produced by an EdgeEnumerator matches GetEdges.
    private static void Verify(S2ShapeIndex index)
    {
        var expected = GetEdges(index);

        int i = 0;
        var it = new EdgeEnumerator(index);
        while (it.MoveNext())
        {
            Assert.True(i < expected.Count);
            Assert.Equal(expected[i], it.Current);
            i++;
        }
    }
}
