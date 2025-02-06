namespace S2Geometry;

public class S2CellUnionTests(ITestOutputHelper logger)
{
    private readonly ITestOutputHelper _logger = logger;

    [Fact]
    internal void Test_S2CellUnion_DefaultConstructor()
    {
        List<S2CellId> ids = [];
        S2CellUnion empty = new(ids);
        Assert.True(empty.IsEmpty());
    }

    [Fact]
    internal void Test_S2CellUnion_S2CellIdConstructor()
    {
        var face1_id = S2CellId.FromFace(1);
        S2CellUnion face1_union = new(new List<S2CellId> { face1_id });
        Assert.Equal(1, face1_union.Size());
        Assert.Equal(face1_id, face1_union.CellId(0));
    }

    [Fact]
    internal void Test_S2CellUnion_WholeSphere()
    {
        var whole_sphere = S2CellUnion.WholeSphere();
        Assert.Equal(whole_sphere.LeafCellsCovered(), 6 * (1UL << 60));
        whole_sphere.Expand(0);
        Assert.Equal(whole_sphere, S2CellUnion.WholeSphere());
    }

    [Fact]
    internal void Test_S2CellUnion_DuplicateCellsNotValid()
    {
        S2CellId id = new(new S2Point(1, 0, 0));
        var cell_union = FromVerbatimNoChecks([id, id]);
        Assert.False(cell_union.IsValid());
    }

    [Fact]
    internal void Test_S2CellUnion_UnsortedCellsNotValid()
    {
        var id = new S2CellId(new S2Point(1, 0, 0)).Parent(10);
        var cell_union = FromVerbatimNoChecks([id, id.Prev()]);
        Assert.False(cell_union.IsValid());
    }

    [Fact]
    internal void Test_S2CellUnion_InvalidCellIdNotValid()
    {
        Assert.False(S2CellId.None.IsValid());
        var cell_union = FromVerbatimNoChecks([S2CellId.None]);
        Assert.False(cell_union.IsValid());
    }

    [Fact]
    internal void Test_S2CellUnion_InvalidCellIdNotValidWithDebugFlag()
    {
        // Manually save and restore flag, to preserve test state in opensource
        // without gflags.
        Assert.False(S2CellId.None.IsValid());
        var cell_union = S2CellUnion.FromVerbatimNoCheck([S2CellId.None]);
        Assert.False(cell_union.IsValid());
    }

    [Fact]
    internal void Test_S2CellUnion_IsNormalized()
    {
        var id = new S2CellId(new S2Point(1, 0, 0)).Parent(10);
        var cell_union = S2CellUnion.FromVerbatim(
            [id.Child(0), id.Child(1), id.Child(2), id.Child(3)]);
        Assert.True(cell_union.IsValid());
        Assert.False(cell_union.IsNormalized());
    }

    [Fact]
    internal void Test_S2CellUnion_Normalize()
    {
        // Try a bunch of random test cases, and keep track of average
        // statistics for normalization (to see if they agree with the
        // analysis above).
        double in_sum = 0, out_sum = 0;
        const int kIters = 2000;
        for (int i = 0; i < kIters; ++i)
        {
            List<S2CellId> input = [];
            List<S2CellId> expected = [];
            AddCells(S2CellId.None, false, input, expected);
            in_sum += input.Count;
            out_sum += expected.Count;
            S2CellUnion cellunion = new(input);
            Assert.Equal(expected.Count, cellunion.Size());
            for (var j = 0; j < expected.Count; ++j)
            {
                Assert.Equal(expected[j], cellunion.CellId(j));
            }

            // Test GetCapBound().
            var cap = cellunion.GetCapBound();
            foreach (var id in cellunion)
            {
                Assert.True(cap.Contains(new S2Cell(id)));
            }

            // Test Contains(S2CellId) and Intersects(S2CellId).
            foreach (var input_id in input)
            {
                Assert.True(cellunion.Contains(input_id));
                Assert.True(cellunion.Contains(input_id.ToPoint()));
                Assert.True(cellunion.Intersects(input_id));
                if (!input_id.IsFace())
                {
                    Assert.True(cellunion.Intersects(input_id.Parent()));
                    if (input_id.Level() > 1)
                    {
                        Assert.True(cellunion.Intersects(input_id.Parent().Parent()));
                        Assert.True(cellunion.Intersects(input_id.Parent(0)));
                    }
                }
                if (!input_id.IsLeaf())
                {
                    Assert.True(cellunion.Contains(input_id.ChildBegin()));
                    Assert.True(cellunion.Intersects(input_id.ChildBegin()));
                    Assert.True(cellunion.Contains(input_id.ChildEnd().Prev()));
                    Assert.True(cellunion.Intersects(input_id.ChildEnd().Prev()));
                    Assert.True(cellunion.Contains(input_id.ChildBegin(S2.kMaxCellLevel)));
                    Assert.True(cellunion.Intersects(input_id.ChildBegin(S2.kMaxCellLevel)));
                }
            }
            foreach (var expected_id in expected)
            {
                if (!expected_id.IsFace())
                {
                    Assert.True(!cellunion.Contains(expected_id.Parent()));
                    Assert.True(!cellunion.Contains(expected_id.Parent(0)));
                }
            }

            // Test Contains(S2CellUnion), Intersects(S2CellUnion), Union(),
            // Intersection(), and Difference().
            List<S2CellId> x = [];
            List<S2CellId> y = [];
            List<S2CellId> x_or_y = [];
            List<S2CellId> x_and_y = [];
            foreach (var input_id in input)
            {
                var in_x = S2Testing.Random.OneIn(2);
                var in_y = S2Testing.Random.OneIn(2);
                if (in_x) x.Add(input_id);
                if (in_y) y.Add(input_id);
                if (in_x || in_y) x_or_y.Add(input_id);
            }
            S2CellUnion xcells = new(x);
            S2CellUnion ycells = new(y);
            S2CellUnion x_or_y_expected = new(x_or_y);
            var x_or_y_cells = xcells.Union(ycells);
            Assert.True(x_or_y_cells == x_or_y_expected);

            // Compute the intersection of "x" with each cell of "y",
            // check that this intersection is correct, and append the
            // results to x_and_y_expected.
            foreach (var yid in ycells)
            {
                var ucells = xcells.Intersection(yid);
                foreach (var xid in xcells)
                {
                    if (xid.Contains(yid))
                    {
                        Assert.True(ucells.Size() == 1 && ucells.CellId(0) == yid);
                    }
                    else if (yid.Contains(xid))
                    {
                        Assert.True(ucells.Contains(xid));
                    }
                }
                foreach (var uid in ucells)
                {
                    Assert.True(xcells.Contains(uid));
                    Assert.True(yid.Contains(uid));
                }
                x_and_y.AddRange(ucells);
            }
            S2CellUnion x_and_y_expected = new(x_and_y);
            var x_and_y_cells = xcells.Intersection(ycells);
            Assert.True(x_and_y_cells == x_and_y_expected);

            var x_minus_y_cells = xcells.Difference(ycells);
            var y_minus_x_cells = ycells.Difference(xcells);
            Assert.True(xcells.Contains(x_minus_y_cells));
            Assert.True(!x_minus_y_cells.Intersects(ycells));
            Assert.True(ycells.Contains(y_minus_x_cells));
            Assert.True(!y_minus_x_cells.Intersects(xcells));
            Assert.True(!x_minus_y_cells.Intersects(y_minus_x_cells));

            var diff_intersection_union = x_minus_y_cells.Union(y_minus_x_cells).Union(x_and_y_cells);
            Assert.True(diff_intersection_union == x_or_y_cells);

            List<S2CellId> test = [];
            List<S2CellId> dummy = [];
            AddCells(S2CellId.None, false, test, dummy);
            foreach (var test_id in test)
            {
                var contains = false;
                var intersects = false;
                foreach (var expected_id in expected)
                {
                    if (expected_id.Contains(test_id)) contains = true;
                    if (expected_id.Intersects(test_id)) intersects = true;
                }
                Assert.Equal(contains, cellunion.Contains(test_id));
                Assert.Equal(intersects, cellunion.Intersects(test_id));
            }
        }
        _logger.WriteLine($"avg in {in_sum / kIters:2f}, avg out {out_sum / kIters:2f}");
    }

    [Fact]
    internal void Test_S2CellUnion_Expand()
    {
        // This test generates coverings for caps of random sizes, expands
        // the coverings by a random radius, and then make sure that the new
        // covering covers the expanded cap.  It also makes sure that the
        // new covering is not too much larger than expected.

        S2RegionCoverer coverer = new();
        for (var i = 0; i < 1000; ++i)
        {
            _logger.WriteLine($"Iteration {i}");
            var cap = S2Testing.GetRandomCap(
                S2Cell.AverageArea(S2.kMaxCellLevel), S2.M_4_PI);

            // Expand the cap area by a random factor whose log is uniformly
            // distributed between 0 and log(1e2).
            var expanded_cap = S2Cap.FromCenterHeight(
                cap.Center, Math.Min(2.0, Math.Pow(1e2, S2Testing.Random.RandDouble()) * cap.Height()));

            var radius = (expanded_cap.Radius - cap.Radius).Radians();
            var max_level_diff = S2Testing.Random.Uniform(8);

            // Generate a covering for the original cap, and measure the maximum
            // distance from the cap center to any point in the covering.
            coverer.Options_.MaxCells = 1 + S2Testing.Random.Skewed(10);
            var covering = coverer.GetCovering(cap);
            S2Testing.CheckCovering(cap, covering, true);
            var covering_radius = GetRadius(covering, cap.Center);

            // This code duplicates the logic in Expand(min_radius, max_level_diff)
            // that figures out an appropriate cell level to use for the expansion.
            int min_level = S2.kMaxCellLevel;
            foreach (var id in covering)
            {
                min_level = Math.Min(min_level, id.Level());
            }
            var expand_level = Math.Min(min_level + max_level_diff,
                S2.kMinWidth.GetLevelForMinValue(radius));

            // Generate a covering for the expanded cap, and measure the new maximum
            // distance from the cap center to any point in the covering.
            covering = covering.Expand(S1Angle.FromRadians(radius), max_level_diff);
            S2Testing.CheckCovering(expanded_cap, covering, false);
            double expanded_covering_radius = GetRadius(covering, cap.Center);

            // If the covering includes a tiny cell along the boundary, in theory the
            // maximum angle of the covering from the cap center can increase by up to
            // twice the maximum length of a cell diagonal.
            Assert.True(expanded_covering_radius - covering_radius <=
                      2 * S2.kMaxDiag.GetValue(expand_level));
        }
    }

    [Fact]
    internal void Test_S2CellUnion_EncodeDecode()
    {
        var cell_ids = new List<S2CellId>{new(0x33),
                           new(0x8e3748fab),
                           new(0x91230abcdef83427)};
        var cell_union = S2CellUnion.FromVerbatim(cell_ids);

        Encoder encoder = new();
        cell_union.Encode(encoder);
        var decoder = encoder.GetDecoder();
        var (success, decoded_cell_union) = S2CellUnion.Decode(decoder);
        Assert.True(success);
        Assert.Equal(cell_union, decoded_cell_union);
    }

    [Fact]
    internal void Test_S2CellUnion_EncodeDecodeEmpty()
    {
        S2CellUnion empty_cell_union = new();

        Encoder encoder = new();
        empty_cell_union.Encode(encoder);
        var decoder = encoder.GetDecoder();
        var (success, decoded_cell_union) = S2CellUnion.Decode(decoder);
        Assert.True(success);
        Assert.Equal(empty_cell_union, decoded_cell_union);
    }

    [Fact]
    internal void Test_S2CellUnion_FromMinMax()
    {
        // Check the very first leaf cell and face cell.
        S2CellId face1_id = S2CellId.FromFace(0);
        TestFromMinMax(face1_id.RangeMin(), face1_id.RangeMin());
        TestFromMinMax(face1_id.RangeMin(), face1_id.RangeMax());

        // Check the very last leaf cell and face cell.
        S2CellId face5_id = S2CellId.FromFace(5);
        TestFromMinMax(face5_id.RangeMin(), face5_id.RangeMax());
        TestFromMinMax(face5_id.RangeMax(), face5_id.RangeMax());

        // Check random ranges of leaf cells.
        for (int iter = 0; iter < 100; ++iter)
        {
            S2CellId x = S2Testing.GetRandomCellId(S2.kMaxCellLevel);
            S2CellId y = S2Testing.GetRandomCellId(S2.kMaxCellLevel);
            if (x > y)
            {
                (y, x) = (x, y);
            }
            TestFromMinMax(x, y);
        }
    }

    [Fact]
    internal void Test_S2CellUnion_FromBeginEnd()
    {
        // Since FromMinMax() is implemented in terms of FromBeginEnd(), we
        // focus on test cases that generate an empty range.
        S2CellId initial_id = S2CellId.FromFace(3);

        // Test an empty range before the minimum S2CellId.
        S2CellUnion cell_union = new(new List<S2CellId> { initial_id });
        S2CellId id_begin = S2CellId.Begin(S2.kMaxCellLevel);
        cell_union.InitFromBeginEnd(id_begin, id_begin);
        Assert.True(cell_union.IsEmpty());

        // Test an empty range after the maximum S2CellId.
        cell_union = new S2CellUnion(new List<S2CellId> { initial_id });
        S2CellId id_end = S2CellId.End(S2.kMaxCellLevel);
        cell_union.InitFromBeginEnd(id_end, id_end);
        Assert.True(cell_union.IsEmpty());

        // Test the full sphere.
        cell_union = S2CellUnion.FromBeginEnd(id_begin, id_end);
        Assert.Equal(6, cell_union.Size());
        foreach (S2CellId id in cell_union)
        {
            Assert.True(id.IsFace());
        }
    }

    [Fact]
    internal void Test_S2CellUnion_EmptyMutableOps()
    {
        S2CellUnion empty_cell_union = new();

        // Normalize()
        empty_cell_union.Normalize();
        Assert.True(empty_cell_union.IsEmpty());

        // Denormalize(...)
        List<S2CellId> output=[];
        empty_cell_union.Denormalize(0, 2, output);
        Assert.True(empty_cell_union.IsEmpty());

        // Pack(...)
        empty_cell_union.Pack();
        Assert.True(empty_cell_union.IsEmpty());

        // Expand(...)
        empty_cell_union.Expand(S1Angle.FromRadians(1), 20);
        Assert.True(empty_cell_union.IsEmpty());
        empty_cell_union.Expand(10);
        Assert.True(empty_cell_union.IsEmpty());
    }

    [Fact]
    internal void Test_S2CellUnion_EmptyAndNonEmptyBooleanOps()
    {
        S2CellUnion empty_cell_union = new();
        S2CellId face1_id = S2CellId.FromFace(1);
        S2CellUnion non_empty_cell_union = new(new List<S2CellId> { face1_id });

        // Contains(...)
        Assert.False(empty_cell_union.Contains(face1_id));
        Assert.True(non_empty_cell_union.Contains(face1_id));
        Assert.True(empty_cell_union.Contains(empty_cell_union));
        Assert.True(non_empty_cell_union.Contains(empty_cell_union));
        Assert.False(empty_cell_union.Contains(non_empty_cell_union));
        Assert.True(non_empty_cell_union.Contains(non_empty_cell_union));

        // Intersects(...)
        Assert.False(empty_cell_union.Intersects(face1_id));
        Assert.True(non_empty_cell_union.Intersects(face1_id));
        Assert.False(empty_cell_union.Intersects(empty_cell_union));
        Assert.False(non_empty_cell_union.Intersects(empty_cell_union));
        Assert.False(empty_cell_union.Intersects(non_empty_cell_union));
        Assert.True(non_empty_cell_union.Intersects(non_empty_cell_union));

        // Union(...)
        Assert.Equal(empty_cell_union, empty_cell_union.Union(empty_cell_union));
        Assert.Equal(non_empty_cell_union, non_empty_cell_union.Union(empty_cell_union));
        Assert.Equal(non_empty_cell_union, empty_cell_union.Union(non_empty_cell_union));
        Assert.Equal(non_empty_cell_union,
                  non_empty_cell_union.Union(non_empty_cell_union));

        // Intersection(...)
        Assert.Equal(empty_cell_union, empty_cell_union.Intersection(face1_id));
        Assert.Equal(non_empty_cell_union, non_empty_cell_union.Intersection(face1_id));
        Assert.Equal(empty_cell_union, empty_cell_union.Intersection(empty_cell_union));
        Assert.Equal(empty_cell_union,
                  non_empty_cell_union.Intersection(empty_cell_union));
        Assert.Equal(empty_cell_union,
                  empty_cell_union.Intersection(non_empty_cell_union));
        Assert.Equal(non_empty_cell_union,
                  non_empty_cell_union.Intersection(non_empty_cell_union));

        // Difference(...)
        Assert.Equal(empty_cell_union, empty_cell_union.Difference(empty_cell_union));
        Assert.Equal(non_empty_cell_union,
                  non_empty_cell_union.Difference(empty_cell_union));
        Assert.Equal(empty_cell_union,
                  empty_cell_union.Difference(non_empty_cell_union));
        Assert.Equal(new S2CellUnion(),
                  non_empty_cell_union.Difference(non_empty_cell_union));
    }

    [Fact]
    internal void Test_S2CellUnion_Clear()
    {
        S2CellId face1_id = S2CellId.FromFace(1);
        S2CellUnion face1_union = new(new List<S2CellId> { face1_id });

        Assert.Equal(1, face1_union.Size());
        Assert.True(1 == face1_union.CellIds.Count);
        Assert.True(1 <= face1_union.CellIds.Capacity);

        face1_union.Clear();
        Assert.Equal(0, face1_union.Size());
        Assert.True(0 == face1_union.CellIds.Count);
        Assert.Equal(0, face1_union.CellIds.Capacity);
    }

    [Fact]
    internal void Test_S2CellUnion_RefuseToDecode()
    {
        List<S2CellId> cellids = [];
        S2CellId id = S2CellId.Begin(S2.kMaxCellLevel);
        for (int i = 0; i <= S2CellUnion.Union_decode_max_num_cells; ++i)
        {
            cellids.Add(id);
            id = id.Next();
        }
        S2CellUnion cell_union = S2CellUnion.FromVerbatim(cellids);
        Encoder encoder = new();
        cell_union.Encode(encoder);
        var decoder = encoder.GetDecoder();
        var (success, _) = S2CellUnion.Decode(decoder);
        Assert.False(success);
    }

    [Fact]
    internal void Test_S2CellUnion_Release()
    {
        S2CellId face1_id = S2CellId.FromFace(1);
        S2CellUnion face1_union = new(new List<S2CellId> { face1_id });
        Assert.Equal(1, face1_union.Size());
        Assert.Equal(face1_id, face1_union.CellId(0));

        var released = face1_union.Release();
        Assert.True(1 == released.Count);
        Assert.Equal(face1_id, released[0]);
        Assert.Equal(0, face1_union.Size());
    }

    [Fact]
    internal void Test_S2CellUnion_LeafCellsCovered()
    {
        S2CellUnion cell_union = new();
        Assert.Equal(0UL, cell_union.LeafCellsCovered());

        var ids = new List<S2CellId>
        {
            // One leaf cell on face 0.
            S2CellId.FromFace(0).ChildBegin(S2.kMaxCellLevel)
        };
        cell_union = new S2CellUnion(ids);
        Assert.Equal(1UL, cell_union.LeafCellsCovered());

        // Face 0 itself (which includes the previous leaf cell).
        ids.Add(S2CellId.FromFace(0));
        cell_union = new S2CellUnion(ids);
        Assert.Equal(1UL << 60, cell_union.LeafCellsCovered());
        // Five faces.
        cell_union = cell_union.Expand(0);
        Assert.Equal(5UL << 60, cell_union.LeafCellsCovered());
        // Whole world.
        cell_union = cell_union.Expand(0);
        Assert.Equal(6UL << 60, cell_union.LeafCellsCovered());

        // Add some disjoint cells.
        ids.Add(S2CellId.FromFace(1).ChildBegin(1));
        ids.Add(S2CellId.FromFace(2).ChildBegin(2));
        ids.Add(S2CellId.FromFace(2).ChildEnd(2).Prev());
        ids.Add(S2CellId.FromFace(3).ChildBegin(14));
        ids.Add(S2CellId.FromFace(4).ChildBegin(27));
        ids.Add(S2CellId.FromFace(4).ChildEnd(15).Prev());
        ids.Add(S2CellId.FromFace(5).ChildBegin(30));
        cell_union = new S2CellUnion(ids);
        ulong expected = 1UL + (1UL << 6) + (1UL << 30) + (1UL << 32) +
            (2UL << 56) + (1UL << 58) + (1UL << 60);
        Assert.Equal(expected, cell_union.LeafCellsCovered());
    }

    [Fact]
    internal void Test_S2CellUnion_WorksInContainers()
    {
        var ids = new List<S2CellId> { S2CellId.FromFace(1) };
        var union_vector = new List<S2CellUnion> { new(ids) };
        Assert.Equal(ids, union_vector.Last().CellIds);
    }

    [Fact]
    internal void Test_S2CellUnion_ToStringEmpty()
    {
        Assert.Equal(new S2CellUnion().ToString(), "Size:0 S2CellIds:");
    }

    [Fact]
    internal void Test_S2CellUnion_ToStringOneCell()
    {
        Assert.Equal(new S2CellUnion(new List<S2CellId> { S2CellId.FromFace(1)}).ToString(),
            "Size:1 S2CellIds:3");
    }

    [Fact]
    internal void Test_S2CellUnion_ToStringTwoCells()
    {
        Assert.Equal(
            new S2CellUnion(new List<S2CellId> { S2CellId.FromFace(1), S2CellId.FromFace(2)}).ToString(),
            "Size:2 S2CellIds:3,5");
    }

    [Fact]
    internal void Test_S2CellUnion_ToStringOver500Cells()
    {
        List<S2CellId> ids=[];
        new S2CellUnion(new List<S2CellId> { S2CellId.FromFace(1)}).Denormalize(6, 1, ids);  // 4096 cells
        var result = S2CellUnion.FromVerbatim(ids).ToString();
        Assert.Equal(result.Count(t => t == ','), 500);
        Assert.Equal(result[^4..], ",...");
    }

    [Fact]
    internal void Test_S2CellUnion_IntersectionOneInputNormalized()
    {
        S2CellId id = S2CellId.FromFace(3);  // arbitrary
        S2CellUnion parent=new(new List<S2CellId> { id });
        S2CellUnion children = S2CellUnion.FromVerbatim(
          [id.Child(0), id.Child(1), id.Child(2), id.Child(3)]);
        S2CellUnion intersection = parent.Intersection(children);
        Assert.Equal(intersection, children);
    }

    // Return the maximum geodesic distance from "axis" to any point of
    // "covering".
    private static double GetRadius(S2CellUnion covering, S2Point axis)
    {
        double max_dist = 0;
        foreach (S2CellId id in covering)
        {
            S2Cell cell = new(id);
            for (int j = 0; j < 4; ++j)
            {
                S2Point a = cell.Vertex(j);
                S2Point b = cell.Vertex(j + 1);
                double dist;
                // The maximum distance is not always attained at a cell vertex: if at
                // least one vertex is in the opposite hemisphere from "axis" then the
                // maximum may be attained along an edge.  We solve this by computing
                // the minimum distance from the edge to (-axis) instead.  We can't
                // simply do this all the time because S2EdgeDistances.GetDistance() has
                // poor accuracy when the result is close to Pi.
                //
                // TODO(ericv): Improve S2EdgeDistances.GetDistance() accuracy near Pi.
                if (a.Angle(axis) > S2.M_PI_2 || b.Angle(axis) > S2.M_PI_2)
                {
                    dist = Math.PI - S2.GetDistance(-axis, a, b).Radians;
                }
                else
                {
                    dist = a.Angle(axis);
                }
                max_dist = Math.Max(max_dist, dist);
            }
        }
        return max_dist;
    }

    private void AddCells(S2CellId id, bool selected, List<S2CellId> input, List<S2CellId> expected)
    {
        // Decides whether to add "id" and/or some of its descendants to the
        // test case.  If "selected" is true, then the region covered by "id"
        // *must* be added to the test case (either by adding "id" itself, or
        // some combination of its descendants, or both).  If cell ids are to
        // the test case "input", then the corresponding expected result after
        // simplification is added to "expected".

        if (id == S2CellId.None)
        {
            // Initial call: decide whether to add cell(s) from each face.
            for (int face = 0; face < 6; ++face)
            {
                AddCells(S2CellId.FromFace(face), false, input, expected);
            }
            return;
        }
        if (id.IsLeaf())
        {
            // The rnd.OneIn() call below ensures that the parent of a leaf cell
            // will always be selected (if we make it that far down the hierarchy).
            Assert.True(selected);
            input.Add(id);
            return;
        }
        // The following code ensures that the probability of selecting a cell
        // at each level is approximately the same, i.e. we test normalization
        // of cells at all levels.
        if (!selected && S2Testing.Random.OneIn(S2.kMaxCellLevel - id.Level()))
        {
            // Once a cell has been selected, the expected output is predetermined.
            // We then make sure that cells are selected that will normalize to
            // the desired output.
            expected.Add(id);
            selected = true;
        }

        // With the rnd.OneIn()ants below, this function adds an average
        // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
        // level at which the cell was first selected (level 15 on average).
        // Therefore the average number of input cells in a test case is about
        // (5/6 * 15 * 6) = 75.  The average number of output cells is about 6.

        // If a cell is selected, we add it to "input" with probability 5/6.
        bool added = false;
        if (selected && !S2Testing.Random.OneIn(6))
        {
            input.Add(id);
            added = true;
        }
        int num_children = 0;
        S2CellId child = id.ChildBegin();
        for (int pos = 0; pos < 4; ++pos, child = child.Next())
        {
            // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
            // This intentionally may result in a cell and some of its children
            // being included in the test case.
            //
            // If the cell is not selected, on average we recurse on one child.
            // We also make sure that we do not recurse on all 4 children, since
            // then we might include all 4 children in the input case by accident
            // (in which case the expected output would not be correct).
            if (S2Testing.Random.OneIn(selected ? 12 : 4) && num_children < 3)
            {
                AddCells(child, selected, input, expected);
                ++num_children;
            }
            // If this cell was selected but the cell itself was not added, we
            // must ensure that all 4 children (or some combination of their
            // descendants) are added.
            if (selected && !added) AddCells(child, selected, input, expected);
        }
    }

    // Creates a possibly invalid S2CellUnion without any checks.
    private static S2CellUnion FromVerbatimNoChecks(List<S2CellId> cell_ids)
    {
        return new S2CellUnion(cell_ids, false);
    }

    private static void TestFromMinMax(S2CellId min_id, S2CellId max_id)
    {
        var cell_union = S2CellUnion.FromMinMax(min_id, max_id);
        var cell_ids = cell_union.CellIds;

        Assert.True(cell_ids.Count > 0);
        Assert.Equal(min_id, cell_ids.First().RangeMin());
        Assert.Equal(max_id, cell_ids.Last().RangeMax());
        for (int i = 1; i < cell_ids.Count; ++i)
        {
            Assert.Equal(cell_ids[i].RangeMin(), cell_ids[i - 1].RangeMax().Next());
        }
        Assert.True(cell_union.IsNormalized());
    }
}
