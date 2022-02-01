namespace S2Geometry;

using static S2ShapeUtil;

public class S2ShapeUtilRangeEnumeratorTests
{
    [Fact]
    public void Test_RangeEnumerator_Next()
    {
        // Create an index with one point each on S2CellId faces 0, 1, and 2.
        var index = MakeIndexOrDie("0:0 | 0:90 | 90:0 # #");
        RangeEnumerator it = new(index);
        Assert.Equal(0UL, it.Id.Face());
        it.MoveNext();
        Assert.Equal(1UL, it.Id.Face());
        it.MoveNext();
        Assert.Equal(2UL, it.Id.Face());
        it.MoveNext();
        Assert.Equal(S2CellId.Sentinel, it.Id);
        Assert.True(it.Done());
    }

    [Fact]
    public void Test_RangeEnumerator_EmptyIndex()
    {
        var empty = MakeIndexOrDie("# #");
        var non_empty = MakeIndexOrDie("0:0 # #");
        RangeEnumerator empty_it = new(empty);
        RangeEnumerator non_empty_it = new(non_empty);
        Assert.False(non_empty_it.Done());
        Assert.True(empty_it.Done());

        empty_it.SeekTo(non_empty_it);
        Assert.True(empty_it.Done());

        empty_it.SeekBeyond(non_empty_it);
        Assert.True(empty_it.Done());

        empty_it.SeekTo(empty_it);
        Assert.True(empty_it.Done());

        empty_it.SeekBeyond(empty_it);
        Assert.True(empty_it.Done());
    }
}
