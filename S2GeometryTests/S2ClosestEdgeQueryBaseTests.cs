using Xunit;

// This is a proof-of-concept prototype of a possible S2FurthestEdgeQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.
using FurthestEdgeQuery = S2Geometry.S2ClosestEdgeQueryBase<S2Geometry.S2MaxDistance>;

namespace S2Geometry
{
    // This file contains some basic tests of the templating support.  Testing of
    // the actual algorithms is in s2closest_edge_query_test.cc.
    public class S2ClosestEdgeQueryBaseTests
    {
        [Fact]
        public void Test_S2ClosestEdgeQueryBase_MaxDistance()
        {
            var index = S2TextFormat.MakeIndexOrDie("0:0 | 1:0 | 2:0 | 3:0 # #");
            FurthestEdgeQuery query = new(index);
            FurthestEdgeQuery.Options options = new();
            options.MaxResults = (1);
            FurthestPointTarget target = new(S2TextFormat.MakePointOrDie("4:0"));
            var results = query.FindClosestEdges(target, options);
            Assert.Single(results);
            Assert.Equal(0, results[0].ShapeId);
            Assert.Equal(0, results[0].EdgeId);
            Assert2.Near(4, results[0].Distance.ToS1ChordAngle().ToAngle().Degrees, 1e-13);
        }

        private sealed class FurthestPointTarget : S2MaxDistancePointTarget
        {
            public FurthestPointTarget(S2Point point) : base(point) { }

            public override int MaxBruteForceIndexSize => 10;  // Arbitrary.
        }
    }
}
