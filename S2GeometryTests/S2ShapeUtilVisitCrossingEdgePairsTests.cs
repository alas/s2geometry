namespace S2Geometry;

using static S2ShapeUtil;

public class S2ShapeUtilVisitCrossingEdgePairsTests
{
    public readonly record struct EdgePair(ShapeEdgeId Item1, ShapeEdgeId Item2)
    {
        public EdgePair(Int32Int32 id1, Int32Int32 id2)
            : this(new ShapeEdgeId(id1.Item1, id1.Item2), new ShapeEdgeId(id2.Item1, id2.Item2))
        {
        }

        public override string ToString() => $"({Item1},{Item2})";
    }

    private readonly ITestOutputHelper _logger;

    public S2ShapeUtilVisitCrossingEdgePairsTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    public void Test_GetCrossingEdgePairs_NoIntersections()
    {
        var index = new MutableS2ShapeIndex();
        TestGetCrossingEdgePairs(index, CrossingType.ALL);
        TestGetCrossingEdgePairs(index, CrossingType.INTERIOR);
    }

    [Fact]
    public void Test_GetCrossingEdgePairs_EdgeGrid()
    {
        var kGridSize = 10;  // (kGridSize + 1) * (kGridSize + 1) crossings
        var index = new MutableS2ShapeIndex();
        var shape = new S2EdgeVectorShape();
        for (int i = 0; i <= kGridSize; ++i)
        {
            shape.Add(S2LatLng.FromDegrees(0, i).ToPoint(),
                      S2LatLng.FromDegrees(kGridSize, i).ToPoint());
            shape.Add(S2LatLng.FromDegrees(i, 0).ToPoint(),
                      S2LatLng.FromDegrees(i, kGridSize).ToPoint());
        }
        index.Add(shape);
        TestGetCrossingEdgePairs(index, CrossingType.ALL);
        TestGetCrossingEdgePairs(index, CrossingType.INTERIOR);
    }

    [Fact]
    public void Test_FindSelfIntersection_Basic()
    {
        // Coordinates are (lat,lng), which can be visualized as (y,x).
        TestHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false);
        TestHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true);  // duplicate vertex
        TestHasCrossing("0:0, 0:1, 1:0, 1:1", true);  // edge crossing
        TestHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true);  // duplicate edge
        TestHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true);  // reversed edge
        TestHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0", true);  // vertex crossing
    }

    // Return true if any loop crosses any other loop (including vertex crossings
    // and duplicate edges), or any loop has a self-intersection (including
    // duplicate vertices).
    private static bool HasSelfIntersection(MutableS2ShapeIndex index)
    {
        if (EdgePairs.FindSelfIntersection(index, out _))
        {
            return true;
        }
        return false;
    }

    // This function recursively verifies that HasCrossing returns the given
    // result for all possible cyclic permutations of the loop vertices for the
    // given set of loops.
    private void TestHasCrossingPermutations(ref List<S2Loop> loops, int i, bool has_crossing)
    {
        if (i == loops.Count)
        {
            MutableS2ShapeIndex index = new();
            S2Polygon polygon = new(loops);
            index.Add(new S2Polygon.Shape(polygon));
            Assert.Equal(has_crossing, HasSelfIntersection(index));
            loops = polygon.Release();
        }
        else
        {
            var orig_loop = loops[i];
            for (int j = 0; j < orig_loop.NumVertices; ++j)
            {
                var vertices = new List<S2Point>();
                for (int k = 0; k < orig_loop.NumVertices; ++k)
                {
                    vertices.Add(orig_loop.Vertex(j + k));
                }
                loops[i] = new S2Loop(vertices, S2Debug.DISABLE);
                TestHasCrossingPermutations(ref loops, i + 1, has_crossing);
            }
            loops[i] = orig_loop;
        }
    }

    // Given a string reprsenting a polygon, and a boolean indicating whether this
    // polygon has any self-intersections or loop crossings, verify that all
    // HasSelfIntersection returns the expected result for all possible cyclic
    // permutations of the loop vertices.
    private void TestHasCrossing(string polygon_str, bool has_crossing)
    {
        // Set S2Debug.DISABLE to allow invalid polygons.
        var polygon = MakePolygonOrDie(polygon_str);
        var loops = polygon.Release();
        TestHasCrossingPermutations(ref loops, 0, has_crossing);
    }

    // A set of edge pairs within an S2ShapeIndex.

    private static List<EdgePair> GetCrossings(S2ShapeIndex index, CrossingType type)
    {
        List<EdgePair> edge_pairs = new();
        EdgePairs.VisitCrossingEdgePairs(index, type, (ShapeEdge a, ShapeEdge b, bool bo) =>
        {
            edge_pairs.Add(new(a.Id, b.Id));
            return true;  // Continue visiting.
        });
        if (edge_pairs.Count > 1)
        {
            edge_pairs = new SortedSet<EdgePair>(edge_pairs).ToList();
        }
        return edge_pairs;
    }

    private static List<EdgePair> GetCrossingEdgePairsBruteForce(S2ShapeIndex index, CrossingType type)
    {
        List<EdgePair> result = new();
        var min_sign = (type == CrossingType.ALL) ? 0 : 1;
        var a_iter = new EdgeEnumerator(index);
        while (a_iter.MoveNext())
        {
            var a = a_iter.Current;
            var b_iter = (EdgeEnumerator)a_iter.CustomClone();
            while (b_iter.MoveNext())
            {
                var b = b_iter.Current;
                if (S2.CrossingSign(a.V0, a.V1, b.V0, b.V1) >= min_sign)
                {
                    result.Add(new(a_iter.GetShapeEdgeId(), b_iter.GetShapeEdgeId()));
                }
            }
        }
        return result;
    }

    private void TestGetCrossingEdgePairs(S2ShapeIndex index, CrossingType type)
    {
        var expected = GetCrossingEdgePairsBruteForce(index, type);
        var actual = GetCrossings(index, type);
        if (!actual.SequenceEqual(expected))
        {
            _logger.WriteLine($@"Unexpected edge pairs; see details below.
Expected number of edge pairs: {expected.Count}
Actual number of edge pairs: {actual.Count}");
            foreach (var edge_pair in expected)
            {
                if (actual.Count(t => t.Equals(edge_pair)) != 1)
                {
                    _logger.WriteLine("Missing value: " + edge_pair);
                }
            }
            foreach (var edge_pair in actual)
            {
                if (expected.Count(t => t.Equals(edge_pair)) != 1)
                {
                    _logger.WriteLine("Extra value: " + edge_pair);
                }
            }
            Assert.True(false);
        }
    }
}
