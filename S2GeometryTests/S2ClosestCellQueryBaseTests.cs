// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_cell_query_test.cc.

namespace S2Geometry;

// This is a proof-of-concept prototype of a possible S2FurthestCellQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.
using FurthestCellQuery = S2ClosestCellQueryBase<S2MaxDistance>;

public class S2ClosestCellQueryBaseTests
{
    [Fact]
    internal void Test_S2ClosestCellQueryBase_MaxDistance()
    {
        S2CellIndex index = new();
        index.Add(MakeCellUnionOrDie("0/123, 0/22, 0/3"), 1 /*label*/);
        index.Build();
        FurthestCellQuery query = new(index);
        FurthestCellQuery.Options options = new();
        options.MaxResults = (1);
        FurthestPointTarget target = new(MakeCellIdOrDie("3/123").ToPoint());
        var results = query.FindClosestCells(target, options);
        Assert.Single(results);
        Assert.Equal("0/123", results[0].CellId.ToString());
        Assert.Equal(1, results[0].Label);
        Assert.Equal(4.0, results[0].Distance.ToS1ChordAngle().Length2);
    }

    private sealed class FurthestPointTarget : S2MaxDistancePointTarget
    {
        internal FurthestPointTarget(S2Point point) : base(point) { }

        public override int MaxBruteForceIndexSize => 10;  // Arbitrary.
    }
}
