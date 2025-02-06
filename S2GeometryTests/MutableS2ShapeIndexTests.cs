namespace S2Geometry;

using static MutableS2ShapeIndex;

public class MutableS2ShapeIndexTests
{
    // This test harness owns a MutableS2ShapeIndex for convenience.
    private readonly MutableS2ShapeIndex index_ = [];

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_SpaceUsed()
    {
        index_.Add(new S2EdgeVectorShape(new(1, 0, 0), new(0, 1, 0)));
        Assert.False(index_.IsFresh());
        int size_before = index_.SpaceUsed();
        Assert.False(index_.IsFresh());

        QuadraticValidate(index_);
        int size_after = index_.SpaceUsed();

        Assert.True(index_.IsFresh());

        Assert.True(size_after > size_before);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_NoEdges()
    {
        TestEnumeratorMethods(index_);
        TestEncodeDecode(index_);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_OneEdge()
    {
        Assert.Equal(0, index_.Add(new S2EdgeVectorShape(
            new(1, 0, 0), new(0, 1, 0))));
        QuadraticValidate(index_);
        TestEnumeratorMethods(index_);
        TestEncodeDecode(index_);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_ShrinkToFitOptimization()
    {
        // This used to trigger a bug in the ShrinkToFit optimization.  The loop
        // below contains almost all of face 0 except for a small region in the
        // 0/00000 subcell.  That subcell is the only one that contains any edges.
        // This caused the index to be built only in that subcell.  However, all the
        // other cells on that face should also have index entries, in order to
        // indicate that they are contained by the loop.
        var loop = S2Loop.MakeRegularLoop(
            new S2Point(1, 0.5, 0.5).Normalize(), S1Angle.FromDegrees(89), 100);
        index_.Add(new S2Loop.Shape(loop));
        QuadraticValidate(index_);
        TestEncodeDecode(index_);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_LoopsSpanningThreeFaces()
    {
        int kNumEdges = 100;
        // Validation is quadratic
        // Construct two loops consisting of kNumEdges vertices each, centered
        // around the cube vertex at the start of the Hilbert curve.
        S2Testing.ConcentricLoopsPolygon(
            new S2Point(1, -1, -1).Normalize(), 2, kNumEdges, out var polygon);
        var loops = polygon.Release();
        foreach (var loop in loops)
        {
            index_.Add(new S2Loop.Shape(loop));
        }
        QuadraticValidate(index_);
        TestEnumeratorMethods(index_);
        TestEncodeDecode(index_);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_ManyIdenticalEdges()
    {
        int kNumEdges = 100;  // Validation is quadratic
        var a = new S2Point(0.99, 0.99, 1).Normalize();
        var b = new S2Point(-0.99, -0.99, 1).Normalize();
        for (int i = 0; i < kNumEdges; ++i)
        {
            Assert.Equal(i, index_.Add(new S2EdgeVectorShape(a, b)));
        }
        QuadraticValidate(index_);
        TestEnumeratorMethods(index_);
        TestEncodeDecode(index_);
        // Since all edges span the diagonal of a face, no subdivision should
        // have occurred (with the default index options).
        for (var it = new MutableS2ShapeIndex.Enumerator(index_, S2ShapeIndex.InitialPosition.BEGIN);
            !it.Done(); it.MoveNext())
        {
            Assert.Equal(0, it.Id.Level());
        }
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_DegenerateEdge()
    {
        // This test verifies that degenerate edges are supported.  The following
        // point is a cube face vertex, and so it should be indexed in 3 cells.
        var a = new S2Point(1, 1, 1).Normalize();
        S2EdgeVectorShape shape = new();
        shape.Add(a, a);
        index_.Add(shape);
        QuadraticValidate(index_);
        TestEncodeDecode(index_);
        // Check that exactly 3 index cells contain the degenerate edge.
        int count = 0;
        for (MutableS2ShapeIndex.Enumerator it = new(index_, S2ShapeIndex.InitialPosition.BEGIN);
            !it.Done(); it.MoveNext(), ++count)
        {
            Assert.True(it.Id.IsLeaf());
            Assert.Equal(1, it.Cell.NumClipped());
            Assert.Equal(1, it.Cell.Clipped(0).NumEdges);
        }
        Assert.Equal(3, count);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_ManyTinyEdges()
    {
        // This test adds many edges to a single leaf cell, to check that
        // subdivision stops when no further subdivision is possible.
        int kNumEdges = 100;  // Validation is quadratic
                              // Construct two points in the same leaf cell.
        var a = new S2CellId(new S2Point(1, 0, 0)).ToPoint();
        var b = (a + new S2Point(0, 1e-12, 0)).Normalize();
        S2EdgeVectorShape shape = new();
        for (int i = 0; i < kNumEdges; ++i)
        {
            shape.Add(a, b);
        }
        index_.Add(shape);
        QuadraticValidate(index_);
        TestEncodeDecode(index_);
        // Check that there is exactly one index cell and that it is a leaf cell.
        MutableS2ShapeIndex.Enumerator it = new(index_, S2ShapeIndex.InitialPosition.BEGIN);
        Assert.True(it.MoveNext());
        Assert.True(it.Id.IsLeaf());
        Assert.False(it.MoveNext());
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_SimpleUpdates()
    {
        // Add 5 loops one at a time, then release them one at a time,
        // validating the index at each step.
        S2Testing.ConcentricLoopsPolygon(new(1, 0, 0), 5, 20, out var polygon);
        for (int i = 0; i < polygon.NumLoops(); ++i)
        {
            index_.Add(new S2Loop.Shape(polygon.Loop(i)));
            QuadraticValidate(index_);
        }
        for (int id = 0; id < polygon.NumLoops(); ++id)
        {
            index_.Release(id);
            Assert.Null(index_.Shape(id));
            QuadraticValidate(index_);
            TestEncodeDecode(index_);
        }
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RandomUpdates()
    {
        // Set the temporary memory budget such that at least one shape needs to be
        // split into multiple update batches (namely, the "5 concentric rings"
        // polygon below which needs ~25KB of temporary space).
        //absl.FlagSaver fs;
        //absl.SetFlag(&FLAGS_s2shape_index_tmp_memory_budget, 10000);

        // Allow the seed to be varied from the command line.
        S2Testing.Random.Reset(S2Testing.Random.RandomSeed);

        // A few polylines.
        index_.Add(new S2Polyline.OwningShape(
            MakePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")));
        index_.Add(new S2Polyline.OwningShape(
            MakePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")));
        index_.Add(new S2Polyline.OwningShape(
            MakePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")));

        // A loop that used to trigger an indexing bug.
        index_.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
            new S2Point(1, 0.5, 0.5).Normalize(), S1Angle.FromDegrees(89), 20)));

        // Five concentric loops.
        S2Testing.ConcentricLoopsPolygon(new S2Point(1, -1, -1).Normalize(), 5, 20, out var polygon5);
        for (int i = 0; i < polygon5.NumLoops(); ++i)
        {
            index_.Add(new S2Loop.Shape(polygon5.Loop(i)));
        }

        // Two clockwise loops around S2Cell cube vertices.
        index_.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
            new S2Point(-1, 1, 1).Normalize(), S1Angle.FromRadians(Math.PI - 0.001), 10)));
        index_.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
            new S2Point(-1, -1, -1).Normalize(), S1Angle.FromRadians(Math.PI - 0.001), 10)));

        // A shape with no edges and no interior.
        index_.Add(new S2Loop.Shape(S2Loop.KEmpty));

        // A shape with no edges that covers the entire sphere.
        index_.Add(new S2Loop.Shape(S2Loop.KFull));

        List<S2Shape> released = [];
        List<int> added = [];
        added.Iota(0, index_.NumShapeIds());
        QuadraticValidate(index_);
        TestEncodeDecode(index_);
        for (int iter = 0; iter < 100; ++iter)
        {
            // Choose some shapes to add and release.
            int num_updates = 1 + S2Testing.Random.Skewed(5);
            for (int n = 0; n < num_updates; ++n)
            {
                if (S2Testing.Random.OneIn(2) && added.Count!=0)
                {
                    int i = S2Testing.Random.Uniform(added.Count);
                    released.Add(index_.Release(added[i]));
                    added.RemoveAt(i);
                }
                else if (released.Count!=0)
                {
                    int i = S2Testing.Random.Uniform(released.Count);
                    var shape = released[i];
                    index_.Add(released[i]);  // Changes shape->id().
                    released.RemoveAt(i);
                    added.Add(shape.Id);
                    //S2_VLOG(1) << "  Added shape " << shape->id()
                            //<< " (" << shape << ")";
                }
            }
            QuadraticValidate(index_);
            TestEncodeDecode(index_);
        }
    }

    [Fact]
    internal void Test_MutableS2ShapeIndex_ConstMethodsThreadSafe()
    {
        // Ensure that lazy updates are thread-safe.  In other words, make sure that
        // nothing bad happens when multiple threads call "const" methods that
        // cause pending updates to be applied.
        LazyUpdatesTest test = new();

        // The number of readers should be large enough so that it is likely that
        // several readers will be running at once (with a multiple-core CPU).
        const int kNumReaders = 8;
        const int kIters = 100;
        test.Run(kNumReaders, kIters);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndex_MixedGeometry()
    {
        // This test used to trigger a bug where the presence of a shape with an
        // interior could cause shapes that don't have an interior to suddenly
        // acquire one.  This would cause extra S2ShapeIndex cells to be created
        // that are outside the bounds of the given geometry.
        List<S2Polyline> polylines =
        [
            MakePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6"),
            MakePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6"),
            MakePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6"),
        ];
        MutableS2ShapeIndex index = [];
        foreach (var polyline in polylines)
        {
            index.Add(new S2Polyline.OwningShape(polyline));
        }
        S2Loop loop = new(new S2Cell(S2CellId.Begin(S2.kMaxCellLevel)));
        index.Add(new S2Loop.Shape(loop));
        // No geometry intersects face 1, so there should be no index cells there.
        Assert.Equal(S2CellRelation.DISJOINT, index.LocateCell(S2CellId.FromFace(1)).cellRelation);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_LinearSpace()
    {
        // Build an index that requires FLAGS_s2shape_index_min_short_edge_fraction
        // to be non-zero in order to use a non-quadratic amount of space.

        // Uncomment the following line to check whether this test works properly.
        // FLAGS_s2shape_index_min_short_edge_fraction = 0;

        // Set the maximum number of "short" edges per cell to 1 so that we can
        // implement this test using a smaller index.
        Options options = new()
        {
            MaxEdgesPerCell = 1
        };
        index_.Options_ = options;
        index_.Init();

        // The idea is to create O(n) copies of a single long edge, along with O(n)
        // clusters of (M + 1) points equally spaced along the long edge, where "M"
        // is the max_edges_per_cell() parameter.  The edges are divided such that
        // there are equal numbers of long and short edges; this maximizes the index
        // size when FLAGS_s2shape_index_min_short_edge_fraction is set to zero.
        const int kNumEdges = 100;  // Validation is quadratic
        int edges_per_cluster = options.MaxEdgesPerCell + 1;
        int num_clusters = kNumEdges / 2 / edges_per_cluster;

        // Create the long edges.
        S2Point a = new(1, 0, 0), b = new(0, 1, 0);
        for (int i = 0; i < kNumEdges / 2; ++i)
        {
            index_.Add(new S2EdgeVectorShape(a, b));
        }
        // Create the clusters of short edges.
        for (int k = 0; k < num_clusters; ++k)
        {
            var p = S2.Interpolate(a, b, k / (num_clusters - 1.0));
            var points = new S2Point[edges_per_cluster];
            points.Fill(p);
            index_.Add(new S2PointVectorShape(points));
        }
        QuadraticValidate(index_);

        // The number of index cells should not exceed the number of clusters.
        int cell_count = 0;
        foreach (var it in index_)
        {
            ++cell_count;
        }
        Assert.True(cell_count <= num_clusters);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_LongIndexEntriesBound()
    {
        // This test demonstrates that the c2 = 366 upper bound (using default
        // parameter values) mentioned in the .cc file is achievable.

        // Set the maximum number of "short" edges per cell to 1 so that we can test
        // using a smaller index.
        Options options = new()
        {
            MaxEdgesPerCell = 1
        };
        index_.Options_ = options;
        index_.Init();

        // This is a worst-case edge AB that touches as many cells as possible at
        // level 30 while still being considered "short" at level 29.  We create an
        // index consisting of two copies of this edge plus a full polygon.
        var a = S2.FaceSiTitoXYZ(0, 0, (1 << 30) + 0).Normalize();
        var b = S2.FaceSiTitoXYZ(0, 0, (1 << 30) + 6).Normalize();
        for (int i = 0; i < 2; ++i)
        {
            index_.Add(new S2EdgeVectorShape(a, b));
        }
        index_.Add(new S2LaxPolygonShape([]));

        // Count the number of index cells at each level.
        var counts = new int[S2.kMaxCellLevel + 1];
        for (MutableS2ShapeIndex.Enumerator it = new(index_, S2ShapeIndex.InitialPosition.BEGIN);
             !it.Done(); it.MoveNext())
        {
            ++counts[it.Id.Level()];
        }
        int sum = 0;
        for (int i = 0; i < counts.Length; ++i)
        {
            //S2_LOG(INFO) << i << ": " << counts[i];
            sum += counts[i];
        }
        Assert.Equal(366, sum);
    }

    // Test move-construct and move-assign functionality of `S2Shape`.  It has an id
    // value which is set when it's added to an index.  So we can create two
    // `S2LaxPolygonShape`s, add them to an index, then
    [Fact]
    internal void Test_MutableS2ShapeIndex_ShapeIdSwaps()
    {
        MutableS2ShapeIndex index =
        [
            MakeLaxPolylineOrDie("1:1, 2:2"),
            MakeLaxPolylineOrDie("3:3, 4:4"),
            MakeLaxPolylineOrDie("5:5, 6:6"),
        ];

        var a = (S2LaxPolylineShape)index.Shape(1)!;
        var b = (S2LaxPolylineShape)index.Shape(2)!;
        Assert.Equal(a.Id, 1);
        Assert.Equal(b.Id, 2);

        // Verify move construction moves the id value.
        S2LaxPolylineShape c=new(a);
        Assert.Equal(c.Id, 1);
        Assert.Equal(c, MakeLaxPolylineOrDie("3:3, 4:4"));

        // Verify move assignment moves the id value.
        S2LaxPolylineShape d;
        d = b;
        Assert.Equal(d.Id, 2);
        Assert.Equal(d, MakeLaxPolylineOrDie("5:5, 6:6"));
    }

    // NOTE(ericv): The tests below are all somewhat fragile since they depend on
    // the internal BatchGenerator heuristics; if these heuristics change
    // (including constants) then the tests below may need to change as well.

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RemoveFullPolygonBatch()
    {
        TestBatchGenerator(0, [], 100 /*bytes*/, 7,
            [new(new(7, 0), new(7, 0), 0)]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_AddFullPolygonBatch()
    {
        TestBatchGenerator(0, [0], 100 /*bytes*/, 7,
            [new(new(7, 0), new(8, 0), 0)]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RemoveManyEdgesInOneBatch()
    {
        // Test removing more edges than would normally fit in a batch.  For good
        // measure we also add two full polygons in the same batch.
        TestBatchGenerator(1000, [0, 0], 100 /*bytes*/, 7,
            [new(new(7, 0), new(9, 0), 1000)]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RemoveAndAddEdgesInOneBatch()
    {
        // Test removing and adding edges in one batch.
        TestBatchGenerator(3, [4, 5], 10000 /*bytes*/, 7,
            [new(new(7, 0), new(9, 0), 12)]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RemoveAndAddEdgesInTwoBatches()
    {
        // Test removing many edges and then adding a few.
        TestBatchGenerator(1000, [3], 1000 /*bytes*/, 7,
            [
                new(new(7, 0), new(7, 0), 1000),
                new(new(7, 0), new(8, 0), 3)
            ]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RemoveAndAddEdgesAndFullPolygonsInTwoBatches()
    {
        // Like the above, but also add two full polygons such that one polygon is
        // processed in each batch.
        TestBatchGenerator(1000, [0, 3, 0], 1000 /*bytes*/, 7,
            [
                new(new(7, 0), new(8, 0), 1000),
                new(new(8, 0), new(10, 0), 3)
            ]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_SeveralShapesInOneBatch()
    {
        // Test adding several shapes in one batch.
        TestBatchGenerator(0, [3, 4, 5], 10000 /*bytes*/, 7,
            [new(new(7, 0), new(10, 0), 12)]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_GroupSmallShapesIntoBatches()
    {
        // Test adding several small shapes that must be split into batches.
        // 10000 bytes ~= temporary space to process 48 edges.
        TestBatchGenerator(0, [20, 20, 20, 20, 20], 10000 /*bytes*/, 7,
            [
                new(new(7, 0), new(9, 0), 40),
                new(new(9, 0), new(11, 0), 40),
                new(new(11, 0), new(12, 0), 20)
            ]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_AvoidPartialShapeInBatch()
    {
        // Test adding a small shape followed by a large shape that won't fit in the
        // same batch as the small shape, but will fit in its own separate batch.
        // 10000 bytes ~= temporary space to process 48 edges.
        TestBatchGenerator(0, [20, 40, 20], 10000 /*bytes*/, 7,
            [
                new(new(7, 0), new(8, 0), 20),
                new(new(8, 0), new(9, 0), 40),
                new(new(9, 0), new(10, 0), 20)
        ]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_SplitShapeIntoTwoBatches()
    {
        // Test adding a few small shapes, then a large shape that can be split
        // across the remainder of the first batch plus the next batch.  The first
        // two batches should have the same amount of remaining space relative to
        // their maximum size.  (For 10000 bytes of temporary space, the ideal batch
        // sizes are 48, 46, 45.)
        //
        // Note that we need a separate batch for the full polygon at the end, even
        // though it has no edges, because partial shapes must always be the last
        // shape in their batch.
        TestBatchGenerator(0, [20, 60, 0], 10000 /*bytes*/, 7,
            [
                new(new(7, 0), new(8, 21), 41),
                new(new(8, 21), new(9, 0), 39),
                new(new(9, 0), new(10, 0), 0)
            ]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_RemoveEdgesAndAddPartialShapeInSameBatch()
    {
        // Test a batch that consists of removing some edges and then adding a
        // partial shape.  We also check that the small shape at the end is put into
        // its own batch, since partial shapes must be the last shape in their batch.
        TestBatchGenerator(20, [60, 5], 10000 /*bytes*/, 7,
            [
                new(new(7, 0), new(7, 21), 41),
                new(new(7, 21), new(8, 0), 39),
                new(new(8, 0), new(9, 0), 5)
            ]);
    }

    [Fact]
    internal void Test_MutableS2ShapeIndexTest_SplitShapeIntoManyBatches()
    {
        // Like the above except that the shape is split into 10 batches.  With
        // 10000 bytes of temporary space, the ideal batch sizes are 63, 61, 59, 57,
        // 55, 53, 51, 49, 48, 46.  The first 8 batches are as full as possible,
        // while the last two batches have the same amount of remaining space
        // relative to their ideal size.  There is also a small batch at the end.
        TestBatchGenerator(0, [20, 500, 5], 10000 /*bytes*/, 7,
            [
                new(new(7, 0), new(8, 43), 63),
                new(new(8, 43), new(8, 104), 61),
                new(new(8, 104), new(8, 163), 59),
                new(new(8, 163), new(8, 220), 57),
                new(new(8, 220), new(8, 275), 55),
                new(new(8, 275), new(8, 328), 53),
                new(new(8, 328), new(8, 379), 51),
                new(new(8, 379), new(8, 428), 49),
                new(new(8, 428), new(8, 465), 37),
                new(new(8, 465), new(9, 0), 35),
                new(new(9, 0), new(10, 0), 5)
            ]);
    }

    // Verifies that that every cell of the index contains the correct edges, and
    // that no cells are missing from the index.  The running time of this
    // function is quadratic in the number of edges.
    private static void QuadraticValidate(MutableS2ShapeIndex index)
    {
        // Iterate through a sequence of nonoverlapping cell ids that cover the
        // sphere and include as a subset all the cell ids used in the index.  For
        // each cell id, verify that the expected set of edges is present.

        // "min_cellid" is the first S2CellId that has not been validated yet.
        var min_cellid = S2CellId.Begin(S2.kMaxCellLevel);
        var it = index.GetNewEnumerator(S2ShapeIndex.InitialPosition.BEGIN);
        var done = it.MoveNext();
        for (;;)
        {
            // Generate a list of S2CellIds ("skipped cells") that cover the gap
            // between the last cell we validated and the next cell in the index.
            S2CellUnion skipped = new();
            if (!done)
            {
                S2CellId cellid = it.Id;
                Assert.True(cellid >= min_cellid);
                skipped.InitFromBeginEnd(min_cellid, cellid.RangeMin());
                min_cellid = cellid.RangeMax().Next();
            }
            else
            {
                // Validate the empty cells beyond the last cell in the index.
                skipped.InitFromBeginEnd(min_cellid, S2CellId.End(S2.kMaxCellLevel));
            }
            // Iterate through all the shapes, simultaneously validating the current
            // index cell and all the skipped cells.
            int num_edges = 0;              // all edges in the cell
            int num_short_edges = 0;        // "short" edges
            int num_containing_shapes = 0;  // shapes containing cell's entry vertex
            for (int id = 0; id < index.NumShapeIds(); ++id)
            {
                S2Shape? shape = index.Shape(id);
                S2ClippedShape? clipped = null;
                if (!done) clipped = it.Current.FindClipped(id);

                // First check that contains_center() is set correctly.
                foreach (var skipped_id in skipped)
                {
                    ValidateInterior(shape, skipped_id, false);
                }
                if (!done)
                {
                    bool contains_center = clipped is not null && clipped.ContainsCenter;
                    ValidateInterior(shape, it.Id, contains_center);
                    S2PaddedCell pcell = new(it.Id, kCellPadding);
                    if (shape is not null)
                    {
                        num_containing_shapes +=
                            shape.ContainsBruteForce(pcell.GetEntryVertex()) ? 1 : 0;
                    }
                }
                // If this shape has been released, it should not be present at all.
                if (shape is null)
                {
                    Assert.Null(clipped);
                    continue;
                }
                // Otherwise check that the appropriate edges are present.
                for (int e = 0; e < shape.NumEdges(); ++e)
                {
                    var edge = shape.GetEdge(e);
                    for (int j = 0; j < skipped.Size(); ++j)
                    {
                        ValidateEdge(edge.V0, edge.V1, skipped.CellId(j), false);
                    }
                    if (!done)
                    {
                        bool has_edge = clipped is not null && clipped.ContainsEdge(e);
                        ValidateEdge(edge.V0, edge.V1, it.Id, has_edge);
                        int max_level = GetEdgeMaxLevel(edge);
                        if (has_edge)
                        {
                            ++num_edges;
                            if (it.Id.Level() < max_level) ++num_short_edges;
                        }
                    }
                }
            }
            // This mirrors the calculation in MutableS2ShapeIndex.MakeIndexCell().
            // It is designed to ensure that the index size is always linear in the
            // number of indexed edges.
            int max_short_edges = Math.Max(
                index.Options_.MaxEdgesPerCell,
                (int)(
                    s2shape_index_min_short_edge_fraction *
                    (num_edges + num_containing_shapes)));
            Assert.True(num_short_edges < max_short_edges);
            if (done) break;

            done = it.MoveNext();
        }
    }

    // Given an edge and a cell id, determines whether or not the edge should be
    // present in that cell and verify that this matches "index_has_edge".
    //
    // Verify that "index_has_edge" is true if and only if the edge AB intersects
    // the given cell id.
    private static void ValidateEdge(S2Point a, S2Point b, S2CellId id, bool index_has_edge)
    {
        // Expand or shrink the padding slightly to account for errors in the
        // function we use to test for intersection (IntersectsRect).
        var padding = kCellPadding;
        padding += (index_has_edge ? 1 : -1) * S2EdgeClipping.kIntersectsRectErrorUVDist;
        var bound = id.BoundUV().Expanded(padding);
        Assert.Equal(S2EdgeClipping.ClipToPaddedFace(a, b, id.Face(), padding, out var a_uv, out var b_uv)
                  && S2EdgeClipping.IntersectsRect(a_uv, b_uv, bound),
                  index_has_edge);
    }

    // Given a shape and a cell id, determines whether or not the shape contains
    // the cell center and verify that this matches "index_contains_center".
    private static void ValidateInterior(S2Shape? shape, S2CellId id, bool index_contains_center)
    {
        if (shape is null)
        {
            Assert.False(index_contains_center);
        }
        else
        {
            Assert.Equal(shape.ContainsBruteForce(id.ToPoint()),
                index_contains_center);
        }
    }

    // Verifies that the index can be encoded and decoded without change.
    private static void TestEncodeDecode(MutableS2ShapeIndex index)
    {
        Encoder encoder = new();
        index.Encode(encoder);
        var decoder = encoder.GetDecoder();
        MutableS2ShapeIndex index2 = [];
        Assert.True(index2.Init(decoder, new S2ShapeUtilCoding.WrappedShapeFactory(index)));
        S2ShapeUtil_Testing.ExpectEqual(index, index2);
    }

    // Converts the given vector of batches to a human-readable form.
    private static string BatchDescriptorsToString(List<BatchDescriptor> batches)
    {
        return string.Join(", ", batches);
    }

    // Verifies that removing and adding the given combination of shapes with
    // the given memory budget yields the expected vector of batches.
    private static void TestBatchGenerator(
        int num_edges_removed, int[] shape_edges_added,
        long tmp_memory_budget, int shape_id_begin,
        List<BatchDescriptor> expected_batches)
    {
        //absl.FlagSaver fs;
        //absl.SetFlag(s2shape_index_tmp_memory_budget, tmp_memory_budget);

        int num_edges_added = 0;
        foreach (var n in shape_edges_added) num_edges_added += n;
        BatchGenerator bgen = new(num_edges_removed, num_edges_added,
                                                 shape_id_begin);
        for (int i = 0; i < shape_edges_added.Length; ++i)
        {
            bgen.AddShape(shape_id_begin + i, shape_edges_added[i]);
        }
        var actual_batches = bgen.Finish();
        Assert.Equal(BatchDescriptorsToString(actual_batches), BatchDescriptorsToString(expected_batches));
    }

    private static void TestEnumeratorMethods(MutableS2ShapeIndex index)
    {
        MutableS2ShapeIndex.Enumerator it = new(index, S2ShapeIndex.InitialPosition.BEGIN);
        Assert.False(it.MovePrevious());
        it.Finish();
        Assert.True(it.Done());
        List<S2CellId> ids = [];
        MutableS2ShapeIndex.Enumerator it2 = new(index);
        S2CellId min_cellid = S2CellId.Begin(S2.kMaxCellLevel);
        for (it.Reset(); !it.Done(); it.MoveNext())
        {
            var cellid = it.Id;
            var skipped = S2CellUnion.FromBeginEnd(min_cellid, cellid.RangeMin());
            foreach (var skipped_id in skipped)
            {
                Assert.False(it2.Locate(skipped_id.ToPoint()));
                Assert.Equal(S2CellRelation.DISJOINT, it2.Locate(skipped_id));
                it2.Reset();
                it2.Seek(skipped_id);
                Assert.Equal(cellid, it2.Id);
            }
            if (ids.Count!=0)
            {
                it2 = it;
                Assert.True(it2.MovePrevious());
                Assert.Equal(ids.Last(), it2.Id);
                it2.MoveNext();
                Assert.Equal(cellid, it2.Id);
                it2.Seek(ids.Last());
                Assert.Equal(ids.Last(), it2.Id);
            }
            it2.Reset();
            Assert.Equal(cellid.ToPoint(), it.Center);
            Assert.True(it2.Locate(it.Center));
            Assert.Equal(cellid, it2.Id);
            it2.Reset();
            Assert.Equal(S2CellRelation.INDEXED, it2.Locate(cellid));
            Assert.Equal(cellid, it2.Id);
            if (!cellid.IsFace())
            {
                it2.Reset();
                Assert.Equal(S2CellRelation.SUBDIVIDED, it2.Locate(cellid.Parent()));
                Assert.True(it2.Id <= cellid);
                Assert.True(it2.Id >= cellid.Parent().RangeMin());
            }
            if (!cellid.IsLeaf())
            {
                for (int i = 0; i < 4; ++i)
                {
                    it2.Reset();
                    Assert.Equal(S2CellRelation.INDEXED, it2.Locate(cellid.Child(i)));
                    Assert.Equal(cellid, it2.Id);
                }
            }
            ids.Add(cellid);
            min_cellid = cellid.RangeMax().Next();
        }
    }

    // A test that repeatedly updates "index_" in one thread and attempts to
    // concurrently read the index_ from several other threads.  When all threads
    // have finished reading, the first thread makes another update.
    //
    // Note that we only test concurrent read access, since MutableS2ShapeIndex
    // requires all updates to be single-threaded and not concurrent with any
    // reads.
    internal class LazyUpdatesTest : ReaderWriterTest
    {
        internal override void WriteOp()
        {
            index_.Clear();
            int num_vertices = 4 * S2Testing.Random.Skewed(10);  // Up to 4K vertices
            S2Loop loop = S2Loop.MakeRegularLoop(
                S2Testing.RandomPoint(), S2Testing.KmToAngle(5), num_vertices);
            index_.Add(new S2Loop.Shape(loop));
        }

        internal override void ReadOp()
        {
            foreach (var it in index_)
            {
                continue;  // NOLINT
            }
        }

        protected MutableS2ShapeIndex index_ = [];
    }
}
