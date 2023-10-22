namespace S2Geometry;

// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_point_query_test.cc.
public class S2ClosestPointQueryBaseTests
{
    [Fact]
    internal void Test_S2ClosestPointQueryBase_MaxDistance()
    {
        S2PointIndex<int> index = [];
        var points = ParsePointsOrDie("0:0, 1:0, 2:0, 3:0");
        for (int i = 0; i < points.Count; ++i)
        {
            index.Add(points[i], i);
        }
        FurthestPointQuery<int> query = new(index);
        FurthestPointQuery<int>.Options options = new()
        {
            MaxResults = 1
        };
        FurthestPointTarget target = new(MakePointOrDie("4:0"));
        var results = query.FindClosestPoints(target, options);
        Assert.Single(results);
        Assert.Equal(points[0], results[0].Point);
        Assert.Equal(0, results[0].Data);
        Assert2.Near(4, results[0].Distance.ToS1ChordAngle().ToAngle().GetDegrees(), 1e-13);
    }

    internal sealed class FurthestPointTarget : S2MaxDistancePointTarget
    {
        internal FurthestPointTarget(S2Point point) : base(point) { }
        public override int MaxBruteForceIndexSize => 10;  // Arbitrary.
    }

    // This is a proof-of-concept prototype of a possible S2FurthestPointQuery
    // class.  The purpose of this test is just to make sure that the code
    // compiles and does something reasonable.
    private class FurthestPointQuery<Data>
        : S2ClosestPointQueryBase<S2MaxDistance, Data>
        where Data : IComparable<Data>
    {
        internal FurthestPointQuery(S2PointIndex<Data> index) : base(index) { }
    }
}
