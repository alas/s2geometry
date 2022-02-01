namespace S2Geometry;

using Index = S2PointIndex<int>;
using PointData = TreeNode<int>;

public class S2PointIndexTests
{
    private readonly Index index_ = new();
    private readonly List<PointData> contents_ = new();

    [Fact]
    public void Test_S2PointIndexTest_NoPoints() => Verify();

    [Fact]
    public void Test_S2PointIndexTest_DuplicatePoints()
    {
        for (int i = 0; i < 10; ++i)
        {
            Add(new S2Point(1, 0, 0), 123);  // All points have same Data argument.
        }
        Verify();
        // Now remove half of the points.
        for (int i = 0; i < 5; ++i)
        {
            Remove(new S2Point(1, 0, 0), 123);
        }
        Verify();
    }

    [Fact]
    public void Test_S2PointIndexTest_RandomPoints()
    {
        for (int i = 0; i < 100; ++i)
        {
            Add(S2Testing.RandomPoint(), S2Testing.Random.Uniform(100));
        }
        Verify();
        // Now remove some of the points.
        for (int i = 0; i < 10; ++i)
        {
            int pos;
            do
            {
                pos = index_.Seek(S2Testing.GetRandomCellId(S2.kMaxCellLevel));
            } while (pos >= index_.NumPoints());
            var val = index_.Get(pos);
            Remove(val!.Value.Point, val.Value.Data);
            Verify();
        }
    }

    [Fact]
    public void Test_S2PointIndex_EmptyData()
    {
        // Verify that points can be added and removed with an empty Data class.
        var index = new Index
            {
                new S2Point(1, 0, 0)
            };
        index.Remove(new S2Point(1, 0, 0));
        Assert.Equal(0, index.NumPoints());
    }

    private void Add(S2Point point, int data)
    {
        index_.Add(point, data);
        contents_.Add(new PointData(new S2CellId(point), point, data));
    }

    private void Remove(S2Point point, int data)
    {
        // If there are multiple copies, remove only one.
        contents_.Remove(new PointData(new S2CellId(point), point, data));
        index_.Remove(point, data);  // Invalidates "point".
    }

    private void Verify()
    {
        VerifyContents();
        VerifyEnumeratorMethods();
    }

    private void VerifyContents()
    {
        var remaining = contents_;
        var it = index_.GetEnumerator();
        while (it.MoveNext())
        {
            var pos = remaining.IndexOf(it.Current);
            Assert.True(pos != -1);
            remaining.RemoveAt(pos);
        }
        Assert.False(remaining.Any());
    }

    private void VerifyEnumeratorMethods()
    {
        // Iterate through all the cells in the index.
        var prev_cellid = S2CellId.None;
        var min_cellid = S2CellId.Begin(S2.kMaxCellLevel);
        foreach (var it in index_)
        {
            var cellid = it.Id;
            Assert.Equal(cellid, new S2CellId(it.Point));
            Assert.True(cellid >= prev_cellid);

            int pos2;
            if (cellid == prev_cellid)
            {
                pos2 = index_.Seek(cellid);
            }

            // Generate a cellunion that covers the range of empty leaf cells between
            // the last cell and this one.  Then make sure that seeking to any of
            // those cells takes us to the immediately following cell.
            if (cellid > prev_cellid)
            {
                foreach (var skipped in S2CellUnion.FromBeginEnd(min_cellid, cellid))
                {
                    pos2 = index_.Seek(skipped);
                    var val = index_.Get(pos2);
                    Assert.Equal(cellid, val?.Id);
                }
            }
            // Test Prev(), Next(), and Seek().
            if (prev_cellid.IsValid())
            {
                pos2 = index_.Seek(cellid);
                Assert.True(--pos2 >= 0);
                var val = index_.Get(pos2);
                Assert.Equal(prev_cellid, val?.Id);
                pos2++;
                val = index_.Get(pos2);
                Assert.Equal(cellid, val?.Id);
                pos2 = index_.Seek(prev_cellid);
                val = index_.Get(pos2);
                Assert.Equal(prev_cellid, val?.Id);
            }
            prev_cellid = cellid;
            min_cellid = cellid.Next();
        }
    }
}
