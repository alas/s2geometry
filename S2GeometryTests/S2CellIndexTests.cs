namespace S2Geometry;

using LabelledCell = S2CellIndex.LabelledCell;

public class S2CellIndexTests
{
    private readonly S2CellIndex _index = new();
    private readonly List<LabelledCell> contents_ = new();
    
    [Fact]
    internal void Test_S2CellIndexTest_Empty()
    {
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_OneFaceCell()
    {
        Add("0/", 0);
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_OneLeafCell()
    {
        Add("1/012301230123012301230123012301", 12);
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_DuplicateValues()
    {
        Add("0/", 0);
        Add("0/", 0);
        Add("0/", 1);
        Add("0/", 17);
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_DisjointCells()
    {
        Add("0/", 0);
        Add("3/", 0);
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_NestedCells()
    {
        // Tests nested cells, including cases where several cells have the same
        // range_min() or range_max() and with randomly ordered labels.
        Add("1/", 3);
        Add("1/0", 15);
        Add("1/000", 9);
        Add("1/00000", 11);
        Add("1/012", 6);
        Add("1/01212", 5);
        Add("1/312", 17);
        Add("1/31200", 4);
        Add("1/3120000", 10);
        Add("1/333", 20);
        Add("1/333333", 18);
        Add("5/", 3);
        Add("5/3", 31);
        Add("5/3333", 27);
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_RandomCellUnions()
    {
        // Construct cell unions from random S2CellIds at random levels.  Note that
        // because the cell level is chosen uniformly, there is a very high
        // likelihood that the cell unions will overlap.
        for (int i = 0; i < 100; ++i) {
            Add(GetRandomCellUnion(), i);
        }
        QuadraticValidate();
    }

    [Fact]
    internal void Test_S2CellIndexTest_ContentsEnumeratorSuppressesDuplicates()
    {
        // Checks that ContentsIterator stops reporting values once it reaches a
        // node of the cell tree that was visited by the previous call to Begin().
        Add("2/1", 1);
        Add("2/1", 2);
        Add("2/10", 3);
        Add("2/100", 4);
        Add("2/102", 5);
        Add("2/1023", 6);
        Add("2/31", 7);
        Add("2/313", 8);
        Add("2/3132", 9);
        Add("3/1", 10);
        Add("3/12", 11);
        Add("3/13", 12);
        QuadraticValidate();

        var contents = new S2CellIndex.ContentsEnumerator(_index);
        ExpectContents("1/123", contents, Array.Empty<(string, int)>());
        ExpectContents("2/100123", contents, new[] { ("2/1", 1), ("2/1", 2), ("2/10", 3), ("2/100", 4) });
        // Check that a second call with the same key yields no additional results.
        ExpectContents("2/100123", contents, Array.Empty<(string, int)>());
        // Check that seeking to a different branch yields only the new values.
        ExpectContents("2/10232", contents, new[] { ("2/102", 5), ("2/1023", 6) });
        // Seek to a node with a different root.
        ExpectContents("2/313", contents, new[] { ("2/31", 7), ("2/313", 8) });
        // Seek to a descendant of the previous node.
        ExpectContents("2/3132333", contents, new[] { ("2/3132", 9) });
        // Seek to an ancestor of the previous node.
        ExpectContents("2/213", contents, Array.Empty<(string, int)>());
        // A few more tests of incremental reporting.
        ExpectContents("3/1232", contents, new (string, int)[] { ("3/1", 10), ("3/12", 11) });
        ExpectContents("3/133210", contents, new[] { ("3/13", 12) });
        ExpectContents("3/133210", contents, Array.Empty<(string, int)>());
        ExpectContents("5/0", contents, Array.Empty<(string, int)>());

        // Now try moving backwards, which is expected to yield values that were
        // already reported above.
        ExpectContents("3/13221", contents, new[] { ("3/1", 10), ("3/13", 12) });
        ExpectContents("2/31112", contents, new[] { ("2/31", 7) });
    }

    [Fact]
    internal void Test_S2CellIndexTest_IntersectionOptimization()
    {
        // Tests various corner cases for the binary search optimization in
        // VisitIntersectingCells.

        Add("1/001", 1);
        Add("1/333", 2);
        Add("2/00", 3);
        Add("2/0232", 4);
        Build();
        TestIntersection(MakeCellUnion(new[] { "1/010", "1/3"}));
        TestIntersection(MakeCellUnion(new[] { "2/010", "2/011", "2/02"}));
    }

    [Fact]
    internal void Test_S2CellIndexTest_IntersectionRandomUnions()
    {
        // Construct cell unions from random S2CellIds at random levels.  Note that
        // because the cell level is chosen uniformly, there is a very high
        // likelihood that the cell unions will overlap.
        for (int i = 0; i < 100; ++i) {
            Add(GetRandomCellUnion(), i);
        }
        Build();
        // Now repeatedly query a cell union constructed in the same way.
        for (int i = 0; i < 200; ++i) {
            TestIntersection(GetRandomCellUnion());
        }
    }

    [Fact]
    internal void Test_S2CellIndexTest_IntersectionSemiRandomUnions() {
        // This test also uses random S2CellUnions, but the unions are specially
        // constructed so that interesting cases are more likely to arise.
        for (int iter = 0; iter < 200; ++iter) {
            Clear();
            var id = S2CellId.FromDebugString("1/0123012301230123");
            var target = new List<S2CellId>();
            for (int i = 0; i < 100; ++i) {
                if (S2Testing.Random.OneIn(10)) Add(id, i);
                if (S2Testing.Random.OneIn(4)) target.Add(id);
                if (S2Testing.Random.OneIn(2)) id = id.NextWrap();
                if (S2Testing.Random.OneIn(6) && !id.IsFace()) id = id.Parent();
                if (S2Testing.Random.OneIn(6) && !id.IsLeaf()) id = id.ChildBegin();
            }
            Build();
            TestIntersection(new S2CellUnion(target));
        }
    }

    // Creates a cell union from a small number of random cells at random levels.
    private static S2CellUnion GetRandomCellUnion()
    {
        var ids = new List<S2CellId>();
        for (int j = 0; j < 10; ++j)
        {
            ids.Add(S2Testing.GetRandomCellId());
        }
        return new S2CellUnion(ids);
    }

    private static S2CellUnion MakeCellUnion(string[] strs)
    {
        var ids = new List<S2CellId>();
        foreach (var str in strs)
        {
            ids.Add(S2CellId.FromDebugString(str));
        }
        return new S2CellUnion(ids);
    }

    // Adds the (cell_id, label) pair to index_ and also contents_ (which is
    // used for independent validation).
    private void Add(S2CellId cell_id, int label)
    {
        _index.Add(cell_id, label);
        contents_.Add(new(cell_id, label));
    }

    private void Add(string cell_str, int label)
    {
        Add(S2CellId.FromDebugString(cell_str), label);
    }

    private void Add(S2CellUnion cell_union, int label)
    {
        _index.Add(cell_union, label);
        foreach (var cell_id in cell_union)
        {
            contents_.Add(new(cell_id, label));
        }
    }

    private void Build() { _index.Build(); }
    private void Clear() { _index.Clear(); }

    // Verifies that the index computes the correct set of (cell_id, label) pairs
    // for every possible leaf cell.  The running time of this function is
    // quadratic in the size of the index.
    private void QuadraticValidate()
    {
        Build();
        VerifyCellEnumerator();
        VerifyIndexContents();
        VerifyRangeEnumerators();
    }

    // Verifies that S2CellIndex.CellIterator visits each (cell_id, label) pair
    // exactly once.
    private void VerifyCellEnumerator()
    {
        var actual = new List<LabelledCell>();
        foreach (var it in _index.GetCellEnumerable())
        {
            actual.Add(new(it.CellId, it.Label));
        }
        ExpectEqual(contents_, actual);
    }
    // Verifies that "expected" and "actual" have the same contents.  Note that
    // duplicate values are allowed.
    private static void ExpectEqual(List<LabelledCell> expected, List<LabelledCell> actual)
    {
        expected.Sort();
        actual.Sort();
        Assert.Equal(expected, actual);
    }

    private void VerifyRangeEnumerators()
    {
        // Test Finish(), which is not otherwise tested below.
        // And also for non-empty ranges.
        var non_empty_finish = _index.GetNERNEnum();
        non_empty_finish.Finish();
        Assert.True(non_empty_finish.Done());

        // Iterate through all the ranges in the index.  We simultaneously iterate
        // through the non-empty ranges and check that the correct ranges are found.
        S2CellId prev_start = S2CellId.None;
        S2CellId non_empty_prev_start = S2CellId.None;
        var it = _index.GetRangeNodeEnumerator();
        var non_empty = _index.GetNERNEnum();
        non_empty.SetPosition(0);
        while (it.MoveNext())
        {
            // Check that seeking in the current range takes us to this range.
            var it2 = _index.GetRangeNodeEnumerator();
            S2CellId start = it.Current.StartId;
            it2.Seek(it.Current.StartId);
            Assert.Equal(start, it2.Current.StartId);
            it2.Seek(it.GetLimitId().Prev());
            Assert.Equal(start, it2.Current.StartId);

            // And also for non-empty ranges.
            var non_empty2 = _index.GetNERNEnum();
            S2CellId non_empty_start = non_empty.Current.StartId;
            non_empty2.Seek(it.Current.StartId);
            Assert.Equal(non_empty_start, non_empty2.Current.StartId);
            non_empty2.Seek(it.GetLimitId().Prev());
            Assert.Equal(non_empty_start, non_empty2.Current.StartId);

            // Test Prev() and Next().
            if (it2.MovePrevious())
            {
                Assert.Equal(prev_start, it2.Current.StartId);
                it2.MoveNext();
                Assert.Equal(start, it2.Current.StartId);
            }
            else
            {
                Assert.Equal(start, it2.Current.StartId);
                Assert.Equal(S2CellId.None, prev_start);
            }

            // And also for non-empty ranges.
            if (non_empty2.MovePrevious())
            {
                Assert.Equal(non_empty_prev_start, non_empty2.Current.StartId);
                non_empty2.MoveNext();
                Assert.Equal(non_empty_start, non_empty2.Current.StartId);
            }
            else
            {
                Assert.Equal(non_empty_start, non_empty2.Current.StartId);
                Assert.Equal(S2CellId.None, non_empty_prev_start);
            }

            // Keep the non-empty iterator synchronized with the regular one.
            if (!it.Current.IsEmpty)
            {
                Assert.Equal(it.Current.StartId, non_empty.Current.StartId);
                Assert.Equal(it.GetLimitId(), non_empty.GetLimitId());
                Assert.False(non_empty.Done());
                non_empty_prev_start = non_empty_start;
                non_empty.MoveNext();
            }
            prev_start = start;
        }
        // Verify that the NonEmptyRangeIterator is also finished.
        Assert.True(non_empty.Done());
    }

    // Verifies that RangeIterator and ContentsIterator can be used to determine
    // the exact set of (s2cell_id, label) pairs that contain any leaf cell.
    private void VerifyIndexContents()
    {
        // "min_cellid" is the first S2CellId that has not been validated yet.
        S2CellId min_cell_id = S2CellId.Begin(S2.kMaxCellLevel);
        var range = _index.GetRangeNodeEnumerator();
        while (range.MoveNext())
        {
            var rn = range.Current;
            Assert.Equal(min_cell_id, rn.StartId);
            Assert.True(min_cell_id < range.GetLimitId());
            Assert.True(range.GetLimitId().IsLeaf());
            min_cell_id = range.GetLimitId();

            // Build a list of expected (cell_id, label) pairs for this range.
            List<LabelledCell> expected = new();
            foreach (var x in contents_)
            {
                if (x.CellId.RangeMin() <= rn.StartId &&
                    x.CellId.RangeMax().Next() >= range.GetLimitId())
                {
                    // The cell contains the entire range.
                    expected.Add(x);
                }
                else
                {
                    // Verify that the cell does not intersect the range.
                    Assert.False(x.CellId.RangeMin() <= range.GetLimitId().Prev() &&
                                    x.CellId.RangeMax() >= rn.StartId);
                }
            }
            List<LabelledCell> actual = new();
            var contents = new S2CellIndex.ContentsEnumerator(_index);
            contents.StartUnion(rn);
            while (contents.MoveNext())
            {
                actual.Add(new(contents.Current.CellId, contents.Current.Label));
            }
            ExpectEqual(expected, actual);
        }
        Assert.Equal(S2CellId.End(S2.kMaxCellLevel), min_cell_id);
    }

    // Tests that VisitIntersectingCells() and GetIntersectingLabels() return
    // correct results for the given target.
    private void TestIntersection(S2CellUnion target)
    {
        List<LabelledCell> expected = new();
        List<LabelledCell> actual = new();
        LabelSet expected_labels = new();
        foreach (var it in _index.GetCellEnumerable())
        {
            if (target.Intersects(it.CellId))
            {
                expected.Add(new(it.CellId, it.Label));
                expected_labels.Add(it.Label);
            }
        }
        _index.VisitIntersectingCells(
            target, (S2CellId cell_id, int label) =>
            {
                actual.Add(new(cell_id, label));
                return true;
            });
        ExpectEqual(expected, actual);
        var actual_labels = _index.GetIntersectingLabels(target);
        Assert.Equal(expected_labels, actual_labels.ToList());
    }

    // Given an S2CellId "target_str" in human-readable form, expects that the
    // first leaf cell contained by this target will intersect the exact set of
    // (cell_id, label) pairs given by "expected_strs".
    private void ExpectContents(string target_str, S2CellIndex.ContentsEnumerator contents, (string, int)[] expected_strs)
    {
        var expected = expected_strs.Select(t => new LabelledCell(S2CellId.FromDebugString(t.Item1), t.Item2))
            .ToList();

        var range = _index.GetRangeNodeEnumerator();
        range.Seek(S2CellId.FromDebugString(target_str).RangeMin());
        contents.StartUnion(range.Current);

        var actual = new List<LabelledCell>();
        while (contents.MoveNext())
        {
            actual.Add(new(contents.Current.CellId, contents.Current.Label));
        }
        ExpectEqual(expected, actual);
    }
}
