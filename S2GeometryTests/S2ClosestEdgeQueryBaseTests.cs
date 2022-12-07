namespace S2Geometry;

// This is a proof-of-concept prototype of a possible S2FurthestEdgeQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.
using FurthestEdgeQuery = S2ClosestEdgeQueryBase<S2MaxDistance>;

// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_edge_query_test.cc.
public class S2ClosestEdgeQueryBaseTests
{
    [Fact]
    internal void Test_S2ClosestEdgeQueryBase_MaxDistance()
    {
        var index = MakeIndexOrDie("0:0 | 1:0 | 2:0 | 3:0 # #");
        FurthestEdgeQuery query = new(index);
        FurthestEdgeQuery.Options options = new();
        options.MaxResults = (1);
        FurthestPointTarget target = new(MakePointOrDie("4:0"));
        var results = query.FindClosestEdges(target, options);
        Assert.Single(results);
        Assert.Equal(0, results[0].ShapeId);
        Assert.Equal(0, results[0].EdgeId);
        Assert2.Near(4, results[0].Distance.ToS1ChordAngle().ToAngle().GetDegrees(), 1e-13);
    }

    private sealed class FurthestPointTarget : S2MaxDistancePointTarget
    {
        internal FurthestPointTarget(S2Point point) : base(point) { }

        public override int MaxBruteForceIndexSize => 10;  // Arbitrary.
    }
}
