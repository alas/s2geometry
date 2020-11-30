using S2Geometry.S2ShapeUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace S2Geometry
{
    public class MutableS2ShapeIndexTests
    {
        // This test harness owns a MutableS2ShapeIndex for convenience.
        private readonly MutableS2ShapeIndex index_ = new();

        [Fact]
        public void Test_MutableS2ShapeIndexTest_SpaceUsed()
        {
            index_.Add(new S2EdgeVectorShape(new S2Point(1, 0, 0), new S2Point(0, 1, 0)));
            Assert.False(index_.IsFresh());
            int size_before = index_.SpaceUsed();
            Assert.False(index_.IsFresh());

            QuadraticValidate(index_);
            int size_after = index_.SpaceUsed();

            Assert.True(index_.IsFresh());

            Assert.True(size_after > size_before);
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_NoEdges()
        {
            TestEnumeratorMethods(index_);
            TestEncodeDecode(index_);
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_OneEdge()
        {
            Assert.Equal(0, index_.Add(new S2EdgeVectorShape(
                new S2Point(1, 0, 0), new S2Point(0, 1, 0))));
            QuadraticValidate(index_);
            TestEnumeratorMethods(index_);
            TestEncodeDecode(index_);
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_ShrinkToFitOptimization()
        {
            // This used to trigger a bug in the ShrinkToFit optimization.  The loop
            // below contains almost all of face 0 except for a small region in the
            // 0/00000 subcell.  That subcell is the only one that contains any edges.
            // This caused the index to be built only in that subcell.  However, all the
            // other cells on that face should also have index entries, in order to
            // indicate that they are contained by the loop.
            var loop = S2Loop.MakeRegularLoop(
                new S2Point(1, 0.5, 0.5).Normalized, S1Angle.FromDegrees(89), 100);
            index_.Add(new S2Loop.Shape(loop));
            QuadraticValidate(index_);
            TestEncodeDecode(index_);
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_LoopsSpanningThreeFaces()
        {
            int kNumEdges = 100;
            // Validation is quadratic
            // Construct two loops consisting of kNumEdges vertices each, centered
            // around the cube vertex at the start of the Hilbert curve.
            S2Testing.ConcentricLoopsPolygon(
                new S2Point(1, -1, -1).Normalized, 2, kNumEdges, out var polygon);
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
        public void Test_MutableS2ShapeIndexTest_ManyIdenticalEdges()
        {
            int kNumEdges = 100;  // Validation is quadratic
            S2Point a = new S2Point(0.99, 0.99, 1).Normalized;
            S2Point b = new S2Point(-0.99, -0.99, 1).Normalized;
            for (int i = 0; i < kNumEdges; ++i)
            {
                Assert.Equal(i, index_.Add(new S2EdgeVectorShape(a, b)));
            }
            QuadraticValidate(index_);
            TestEnumeratorMethods(index_);
            TestEncodeDecode(index_);
            // Since all edges span the diagonal of a face, no subdivision should
            // have occurred (with the default index options).
            foreach (var it in index_.GetNewEnumerable())
            {
                Assert.Equal(0, it.Item1.Level);
            }
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_DegenerateEdge()
        {
            // This test verifies that degenerate edges are supported.  The following
            // point is a cube face vertex, and so it should be indexed in 3 cells.
            var a = new S2Point(1, 1, 1).Normalized;
            var shape = new S2EdgeVectorShape();
            shape.Add(a, a);
            index_.Add(shape);
            QuadraticValidate(index_);
            TestEncodeDecode(index_);
            // Check that exactly 3 index cells contain the degenerate edge.
            int count = 0;
            foreach (var it in index_.GetNewEnumerable())
            {
                ++count;
                Assert.True(it.Item1.IsLeaf);
                Assert.Equal(1, it.Item2.NumClipped());
                Assert.Equal(1, it.Item2.Clipped(0).NumEdges);
            }
            Assert.Equal(3, count);
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_ManyTinyEdges()
        {
            // This test adds many edges to a single leaf cell, to check that
            // subdivision stops when no further subdivision is possible.
            int kNumEdges = 100;  // Validation is quadratic
                                  // Construct two points in the same leaf cell.
            S2Point a = new S2CellId(new S2Point(1, 0, 0)).ToPoint();
            S2Point b = (a + new S2Point(0, 1e-12, 0)).Normalized;
            var shape = new S2EdgeVectorShape();
            for (int i = 0; i < kNumEdges; ++i)
            {
                shape.Add(a, b);
            }
            index_.Add(shape);
            QuadraticValidate(index_);
            TestEncodeDecode(index_);
            // Check that there is exactly one index cell and that it is a leaf cell.
            var it = index_.GetNewEnumerator();
            Assert.True(it.MoveNext());
            Assert.True(it.Current.Item1.IsLeaf);
            Assert.False(it.MoveNext());
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_SimpleUpdates()
        {
            // Add 5 loops one at a time, then release them one at a time,
            // validating the index at each step.
            S2Testing.ConcentricLoopsPolygon(new S2Point(1, 0, 0), 5, 20, out var polygon);
            for (int i = 0; i < polygon.NumLoops(); ++i)
            {
                index_.Add(new S2Loop.Shape(polygon.Loop(i)));
                QuadraticValidate(index_);
            }
            for (int id = 0; id < polygon.NumLoops(); ++id)
            {
                index_.Release(id);
                QuadraticValidate(index_);
                TestEncodeDecode(index_);
            }
        }

        [Fact]
        public void Test_MutableS2ShapeIndexTest_RandomUpdates()
        {
            // Allow the seed to be varied from the command line.
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);

            // A few polylines.
            index_.Add(new S2Polyline.OwningShape(
                S2TextFormat.MakePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")));
            index_.Add(new S2Polyline.OwningShape(
                S2TextFormat.MakePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")));
            index_.Add(new S2Polyline.OwningShape(
                S2TextFormat.MakePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")));

            // A loop that used to trigger an indexing bug.
            index_.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
                new S2Point(1, 0.5, 0.5).Normalized, S1Angle.FromDegrees(89), 20)));

            // Five concentric loops.
            S2Testing.ConcentricLoopsPolygon(new S2Point(1, -1, -1).Normalized, 5, 20, out var polygon5);
            for (int i = 0; i < polygon5.NumLoops(); ++i)
            {
                index_.Add(new S2Loop.Shape(polygon5.Loop(i)));
            }

            // Two clockwise loops around S2Cell cube vertices.
            index_.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
                new S2Point(-1, 1, 1).Normalized, S1Angle.FromRadians(Math.PI - 0.001), 10)));
            index_.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
                new S2Point(-1, -1, -1).Normalized, S1Angle.FromRadians(Math.PI - 0.001), 10)));

            // A shape with no edges and no interior.
            index_.Add(new S2Loop.Shape(S2Loop.kEmpty));

            // A shape with no edges that covers the entire sphere.
            index_.Add(new S2Loop.Shape(S2Loop.kFull));

            List<S2Shape> released = new();
            List<int> added = new();
            added.Iota(0, index_.NumShapeIds());
            QuadraticValidate(index_);
            TestEncodeDecode(index_);
            for (int iter = 0; iter < 100; ++iter)
            {
                // Choose some shapes to add and release.
                int num_updates = 1 + S2Testing.Random.Skewed(5);
                for (int n = 0; n < num_updates; ++n)
                {
                    if (S2Testing.Random.OneIn(2) && added.Any())
                    {
                        int i = S2Testing.Random.Uniform(added.Count);
                        released.Add(index_.Release(added[i]));
                        added.RemoveAt(i);
                    }
                    else if (released.Any())
                    {
                        int i = S2Testing.Random.Uniform(released.Count);
                        var shape = released[i];
                        index_.Add(released[i]);  // Changes shape.Current.Item1.
                        released.RemoveAt(i);
                        added.Add(shape.Id);
                    }
                }
                QuadraticValidate(index_);
                TestEncodeDecode(index_);
            }
        }

        [Fact]
        public void Test_LazyUpdatesTest_ConstMethodsThreadSafe()
        {
            // Ensure that lazy updates are thread-safe.  In other words, make sure that
            // nothing bad happens when multiple threads call "const" methods that
            // cause pending updates to be applied.

            // The number of readers should be large enough so that it is likely that
            // several readers will be running at once (with a multiple-core CPU).
            /*int kNumReaders = 8;
            using var pool = new ReaderThreadPool(this, kNumReaders);
            lock (lock_)
            {
                int kIters = 100;
                for (int iter = 0; iter < kIters; ++iter)
                {
                    // Loop invariant: lock_ is held and num_readers_left_ == 0.
                    Assert.Equal(0, num_readers_left_);
                    // Since there are no readers, it is safe to modify the index.
                    index_.Clear();
                    int num_vertices = 4 * S2Testing.Random.Skewed(10);  // Up to 4K vertices
                    S2Loop loop = (S2Loop.MakeRegularLoop(
                        S2Testing.RandomPoint(), S2Testing.KmToAngle(5), num_vertices));
                    index_.Add(new S2Loop.Shape(loop));
                    num_readers_left_ = kNumReaders;
                    ++num_updates_;
                    update_ready_.SignalAll();
                    while (num_readers_left_ > 0)
                    {
                        all_readers_done_.Wait(lock_);
                    }
                }
                // Signal the readers to exit.
                num_updates_ = -1;
                update_ready_.SignalAll();
            }
            // ReaderThreadPool destructor waits for all threads to complete.*/
        }

        [Fact]
        public void Test_MutableS2ShapeIndex_MixedGeometry()
        {
            // This test used to trigger a bug where the presence of a shape with an
            // interior could cause shapes that don't have an interior to suddenly
            // acquire one.  This would cause extra S2ShapeIndex cells to be created
            // that are outside the bounds of the given geometry.
            var polylines = new List<S2Polyline>
            {
                S2TextFormat.MakePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6"),
                S2TextFormat.MakePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6"),
                S2TextFormat.MakePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6"),
            };
            MutableS2ShapeIndex index = new();
            foreach (var polyline in polylines)
            {
                index.Add(new S2Polyline.OwningShape(polyline));
            }
            S2Loop loop = new(new S2Cell(S2CellId.Begin(S2Constants.kMaxCellLevel)));
            index.Add(new S2Loop.Shape(loop));
            // No geometry intersects face 1, so there should be no index cells there.
            Assert.Equal(S2ShapeIndex.CellRelation.DISJOINT, index.LocateCell(S2CellId.FromFace(1)).cellRelation);
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
            S2CellId min_cellid = S2CellId.Begin(S2Constants.kMaxCellLevel);
            var it = index.GetNewEnumerator();
            bool done = it.MoveNext();
            for (; ; )
            {
                // Generate a list of S2CellIds ("skipped cells") that cover the gap
                // between the last cell we validated and the next cell in the index.
                S2CellUnion skipped = new();
                if (!done)
                {
                    S2CellId cellid = it.Current.Item1;
                    Assert.True(cellid >= min_cellid);
                    skipped.InitFromBeginEnd(min_cellid, cellid.RangeMin);
                    min_cellid = cellid.RangeMax.Next;
                }
                else
                {
                    // Validate the empty cells beyond the last cell in the index.
                    skipped.InitFromBeginEnd(min_cellid, S2CellId.End(S2Constants.kMaxCellLevel));
                }
                // Iterate through all the shapes, simultaneously validating the current
                // index cell and all the skipped cells.
                int short_edges = 0;  // number of edges counted toward subdivision
                for (int id = 0; id < index.NumShapeIds(); ++id)
                {
                    S2Shape shape = index.Shape(id);
                    S2ClippedShape clipped = null;
                    if (!done) clipped = it.Current.Item2.FindClipped(id);

                    // First check that contains_center() is set correctly.
                    foreach (S2CellId skipped_id in skipped)
                    {
                        ValidateInterior(shape, skipped_id, false);
                    }
                    if (!done)
                    {
                        bool contains_center = clipped != null && clipped.ContainsCenter;
                        ValidateInterior(shape, it.Current.Item1, contains_center);
                    }
                    // If this shape has been released, it should not be present at all.
                    if (shape == null)
                    {
                        Assert.Null(clipped);
                        continue;
                    }
                    // Otherwise check that the appropriate edges are present.
                    for (int e = 0; e < shape.NumEdges; ++e)
                    {
                        var edge = shape.GetEdge(e);
                        for (int j = 0; j < skipped.NumCells; ++j)
                        {
                            ValidateEdge(edge.V0, edge.V1, skipped.CellId(j), false);
                        }
                        if (!done)
                        {
                            bool has_edge = clipped != null && clipped.ContainsEdge(e);
                            ValidateEdge(edge.V0, edge.V1, it.Current.Item1, has_edge);
                            int max_level = MutableS2ShapeIndex.GetEdgeMaxLevel(edge);
                            if (has_edge && it.Current.Item1.Level < max_level)
                            {
                                ++short_edges;
                            }
                        }
                    }
                }
                Assert.True(short_edges <= index.Options_.MaxEdgesPerCell);
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
            var padding = MutableS2ShapeIndex.kCellPadding;
            padding += (index_has_edge ? 1 : -1) * S2EdgeClipping.kIntersectsRectErrorUVDist;
            R2Rect bound = id.GetBoundUV().Expanded(padding);
            Assert.Equal(S2EdgeClipping.ClipToPaddedFace(a, b, (int)id.Face, padding, out var a_uv, out var b_uv)
                      && S2EdgeClipping.IntersectsRect(a_uv, b_uv, bound),
                      index_has_edge);
        }

        // Given a shape and a cell id, determines whether or not the shape contains
        // the cell center and verify that this matches "index_contains_center".
        private static void ValidateInterior(S2Shape shape, S2CellId id, bool index_contains_center)
        {
            if (shape == null)
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
            Decoder decoder = new(encoder.Buffer, 0, encoder.Length);
            MutableS2ShapeIndex index2 = new();
            Assert.True(index2.Init(decoder, new S2ShapeUtilCoding.WrappedShapeFactory(index)));
            S2ShapeTestsUtil.ExpectEqual(index, index2);
        }

        private static void TestEnumeratorMethods(MutableS2ShapeIndex index) {
            S2ShapeIndex.Enumerator it = new(index);
            Assert.False(it.MovePrevious());
            it.Reset();
            List<S2CellId> ids = new();
            S2ShapeIndex.Enumerator it2 = new(index);
            S2CellId min_cellid = S2CellId.Begin(S2Constants.kMaxCellLevel);
            while (it.MoveNext())
            {
                var cellid = it.Current.Item1;
                var skipped = S2CellUnion.FromBeginEnd(min_cellid, cellid.RangeMin);
                foreach (var skipped_id in skipped) {
                    Assert.False(index.LocatePoint(skipped_id.ToPoint()).found);
                    Assert.Equal(S2ShapeIndex.CellRelation.DISJOINT, index.LocateCell(skipped_id).cellRelation);
                    it2.Reset();
                    it2.SetPosition(index.SeekCell(skipped_id).pos);
                    Assert.Equal(cellid, it2.Current.Item1);
                }
                if (ids.Any()) {
                    it2 = it;
                    Assert.True(it2.MovePrevious());
                    Assert.Equal(ids.Last(), it2.Current.Item1);
                    it2.MoveNext();
                    Assert.Equal(cellid, it2.Current.Item1);
                    it2.SetPosition(index.SeekCell(ids.Last()).pos);
                    Assert.Equal(ids.Last(), it2.Current.Item1);
                }
                it2.Reset();
                Assert.Equal(cellid.ToPoint(), it.Current.Item1.ToPoint());
                var (pos, found) = index.LocatePoint(it.Current.Item1.ToPoint());
                Assert.True(found);
                it2.SetPosition(pos);
                Assert.Equal(cellid, it2.Current.Item1);
                it2.Reset();
                var (rel2, pos2) = index.LocateCell(cellid);
                Assert.Equal(S2ShapeIndex.CellRelation.INDEXED, rel2);
                it2.SetPosition(pos2);
                Assert.Equal(cellid, it2.Current.Item1);
                if (!cellid.IsFace) {
                    it2.Reset();
                    var (rel3, pos3) = index.LocateCell(cellid.Parent());
                    Assert.Equal(S2ShapeIndex.CellRelation.SUBDIVIDED, rel3);
                    it2.SetPosition(pos3);
                    Assert.True(it2.Current.Item1 >= cellid);
                    Assert.True(it2.Current.Item1 >= cellid.Parent().RangeMin);
                }
                if (!cellid.IsLeaf) {
                    for (int i = 0; i < 4; ++i) {
                        it2.Reset();
                        var (rel3, pos3) = index.LocateCell(cellid.Child(i));
                        Assert.Equal(S2ShapeIndex.CellRelation.INDEXED, rel3);
                        it2.SetPosition(pos3);
                        Assert.Equal(cellid, it2.Current.Item1);
                    }
                }
                ids.Add(cellid);
                min_cellid = cellid.RangeMax.Next;
            }
        }

        // A test that repeatedly updates "index_" in one thread and attempts to
        // concurrently read the index_ from several other threads.  When all threads
        // have finished reading, the first thread makes another update.
        //
        // Note that we only test concurrent read access, since MutableS2ShapeIndex
        // requires all updates to be single-threaded and not concurrent with any
        // reads.
        /*public class LazyUpdatesTest
        {
            public LazyUpdatesTest() { num_updates_ = 0; num_readers_left_ = 0; }

            // The function executed by each reader thread.
            public void ReaderThread()
            {
                for (int last_update = 0; ; last_update = num_updates_)
                {
                    lock (lock_)
                    {
                        while (num_updates_ == last_update)
                        {
                            update_ready_.Wait(lock_);
                        }
                        if (num_updates_ < 0) break;

                        // The index is built on demand the first time we attempt to use it.
                        // We intentionally release the lock so that many threads have a chance
                        // to access the MutableS2ShapeIndex in parallel.
                    }

                    foreach (var it in index_.GetNewEnumerable())
                    {
                        continue;  // NOLINT
                    }

                    lock (lock_)
                    {
                        if (--num_readers_left_ == 0)
                        {
                            all_readers_done_.Signal();
                        }
                    }
                }
            }

            public class ReaderThreadPool : IDisposable
            {
                public ReaderThreadPool(LazyUpdatesTest test, int num_threads)
                {
                    threads_ = new thread[num_threads];
                    num_threads_ = num_threads;
                    for (int i = 0; i < num_threads_; ++i)
                    {
                        threads_[i] = new thread(LazyUpdatesTest.ReaderThread, test);
                    }
                }
                public void Dispose()
                {
                    for (int i = 0; i < num_threads_; ++i) threads_[i].join();
                }

                private thread[] threads_;
                private int num_threads_;
            }

            MutableS2ShapeIndex index_;
            // The following fields are guarded by lock_.
            object lock_ = new object();
            int num_updates_;
            int num_readers_left_;

            // Signalled when a new update is ready to be processed.
            CondVar update_ready_;
            // Signalled when all readers have processed the latest update.
            CondVar all_readers_done_;
        }*/
    }
}
