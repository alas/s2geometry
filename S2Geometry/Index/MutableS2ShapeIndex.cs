using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

namespace S2Geometry
{
    using S2ShapeUtil;
    using S2ShapeIndexIdCell = KeyData<S2CellId, S2ShapeIndexCell>;

    // MutableS2ShapeIndex is a class for in-memory indexing of polygonal geometry.
    // The objects in the index are known as "shapes", and may consist of points,
    // polylines, and/or polygons, possibly overlapping.  The index makes it very
    // fast to answer queries such as finding nearby shapes, measuring distances,
    // testing for intersection and containment, etc.
    //
    // MutableS2ShapeIndex allows not only building an index, but also updating it
    // incrementally by adding or removing shapes (hence its name).  It is one of
    // several implementations of the S2ShapeIndex interface.  MutableS2ShapeIndex
    // is designed to be compact; usually it is smaller than the underlying
    // geometry being indexed.  It is capable of indexing up to hundreds of
    // millions of edges.  The index is also fast to construct.
    //
    // There are a number of built-in classes that work with S2ShapeIndex objects.
    // Generally these classes accept any collection of geometry that can be
    // represented by an S2ShapeIndex, i.e. any combination of points, polylines,
    // and polygons.  Such classes include:
    //
    // - S2ContainsPointQuery: returns the shape(s) that contain a given point.
    //
    // - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
    //                       S2CellId, or S2ShapeIndex.
    //
    // - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
    //
    // - S2BooleanOperation: computes boolean operations such as union,
    //                       and boolean predicates such as containment.
    //
    // - S2ShapeIndexRegion: computes approximations for a collection of geometry.
    //
    // - S2ShapeIndexBufferedRegion: computes approximations that have been
    //                               expanded by a given radius.
    //
    // Here is an example showing how to build an index for a set of polygons, and
    // then then determine which polygon(s) contain each of a set of query points:
    //
    //   void TestContainment(S2Point[] points,
    //                        S2Polygon[] polygons) {
    //     MutableS2ShapeIndex index;
    //     foreach (var polygon in polygons) {
    //       index.Add(new S2Polygon.Shape(polygon));
    //     }
    //     var query = index.MakeS2ContainsPointQuery();
    //     foreach (var& point in points) {
    //       foreach (var shape in query.GetContainingShapes(point)) {
    //         S2Polygon* polygon = polygons[shape.id()];
    //         ... do something with (point, polygon) ...
    //       }
    //     }
    //   }
    //
    // This example uses S2Polygon.Shape, which is one example of an S2Shape
    // object.  S2Polyline and S2Loop also have nested Shape classes, and there are
    // additional S2Shape types defined in *_shape.h.
    //
    // Internally, MutableS2ShapeIndex is essentially a map from S2CellIds to the
    // set of shapes that intersect each S2CellId.  It is adaptively refined to
    // ensure that no cell contains more than a small number of edges.
    //
    // For efficiency, updates are batched together and applied lazily on the
    // first subsequent query.  Locking is used to ensure that MutableS2ShapeIndex
    // has the same thread-safety properties as "vector": methods are
    // thread-safe, while non-methods are not thread-safe.  This means that
    // if one thread updates the index, you must ensure that no other thread is
    // reading or updating the index at the same time.
    //
    // TODO(ericv): MutableS2ShapeIndex has an Encode() method that allows the
    // index to be serialized.  An encoded S2ShapeIndex can be decoded either into
    // its original form (MutableS2ShapeIndex) or into an EncodedS2ShapeIndex.
    // The key property of EncodedS2ShapeIndex is that it can be constructed
    // instantaneously, since the index is kept in its original encoded form.
    // Data is decoded only when an operation needs it.  For example, to determine
    // which shapes(s) contain a given query point only requires decoding the data
    // in the S2ShapeIndexCell that contains that point.
    public sealed class MutableS2ShapeIndex : S2ShapeIndex, IDisposable
    {
        // The default maximum number of edges per cell (not counting "long" edges).
        // If a cell has more than this many edges, and it is not a leaf cell, then it
        // is subdivided.  This flag can be overridden via MutableS2ShapeIndex.Options.
        // Reasonable values range from 10 to about 50 or so.  Small values makes
        // queries faster, while large values make construction faster and use less memory.
        private const int s2shape_index_default_max_edges_per_cell = 10;

        // Attempt to limit the amount of temporary memory allocated while building or
        // updating a MutableS2ShapeIndex to at most this value.  This is achieved by
        // splitting the updates into multiple batches when necessary.  (The memory
        // required is proportional to the number of edges being updated at once.)
        //
        // Note that this limit is not a hard guarantee, for several reasons:
        //  (1) the memory estimates are only approximations;
        //  (2) all edges in a given shape are added or removed at once, so shapes
        //      with huge numbers of edges may exceed the budget;
        //  (3) shapes being removed are always processed in a single batch.  (This
        //      could be fixed, but it seems better to keep the code simpler for now.)
        //
        // If more memory than this is needed, updates will automatically be split
        // into batches internally.
        private const int s2shape_index_tmp_memory_budget_mb = 100;

        // The cell size relative to the length of an edge at which it is first
        // considered to be "long".  Long edges do not contribute toward the decision
        // to subdivide a cell further.  For example, a value of 2.0 means that the
        // cell must be at least twice the size of the edge in order for that edge to
        // be counted.  There are two reasons for not counting long edges: (1) such
        // edges typically need to be propagated to several children, which increases
        // time and memory costs without much benefit, and (2) in pathological cases,
        // many long edges close together could force subdivision to continue all the
        // way to the leaf cell level.
        //
        // The size and speed of the index are typically not very sensitive to this
        // parameter.  Reasonable values range from 0.1 to 10, with smaller values
        // causing more aggressive subdivision of long edges grouped closely together.
        private const double s2shape_index_cell_into_long_edge_ratio = 1.0;

        // The amount by which cells are "padded" to compensate for numerical errors
        // when clipping line segments to cell boundaries.
        //
        // The total error when clipping an edge comes from two sources:
        // (1) Clipping the original spherical edge to a cube face (the "face edge").
        //     The maximum error in this step is S2EdgeClipping.kFaceClipErrorUVCoord.
        // (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
        //     The maximum error in this step is S2.kEdgeClipErrorUVCoord.
        // Finally, since we encounter the same errors when clipping query edges, we
        // double the total error so that we only need to pad edges during indexing
        // and not at query time.
#pragma warning disable IDE1006 // Estilos de nombres
        public static readonly double kCellPadding = 2 * (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kEdgeClipErrorUVCoord);
#pragma warning restore IDE1006 // Estilos de nombres

        // When adding a new encoding, be aware that old binaries will not be able
        // to decode it.
#pragma warning disable IDE1006 // Estilos de nombres
        public const byte kCurrentEncodingVersionNumber = 0;
#pragma warning restore IDE1006 // Estilos de nombres

        // Options that affect construction of the MutableS2ShapeIndex.
        public class Options
        {
            public Options()
            {
                MaxEdgesPerCell = s2shape_index_default_max_edges_per_cell;
            }

            // The maximum number of edges per cell.  If a cell has more than this
            // many edges that are not considered "long" relative to the cell size,
            // then it is subdivided.  (Whether an edge is considered "long" is
            // controlled by --s2shape_index_cell_into_long_edge_ratio flag.)
            //
            // Values between 10 and 50 represent a reasonable balance between memory
            // usage, construction time, and query time.  Small values make queries
            // faster, while large values make construction faster and use less memory.
            // Values higher than 50 do not save significant additional memory, and
            // query times can increase substantially, especially for algorithms that
            // visit all pairs of potentially intersecting edges (such as polygon
            // validation), since this is quadratic in the number of edges per cell.
            //
            // Note that the *average* number of edges per cell is generally slightly
            // less than half of the maximum value defined here.
            //
            // Defaults to value given by --s2shape_index_default_max_edges_per_cell.
            public int MaxEdgesPerCell { get; set; }
        }

        // Creates a MutableS2ShapeIndex that uses the default option settings.
        // Option values may be changed by calling Init().
        public MutableS2ShapeIndex() { Options_ = new Options(); InitChannel(); }

        // Create a MutableS2ShapeIndex with the given options.
        public MutableS2ShapeIndex(Options options) { Options_ = options; InitChannel(); }

        public void Dispose() { Clear(); }

        // Initialize a MutableS2ShapeIndex with the given options.  This method may
        // only be called when the index is empty (i.e. newly created or Reset() has
        // just been called).
        public void Init(Options options)
        {
            Assert.True(!shapes_.Any());
            Options_ = options;
        }

        // The options supplied for this index.
        public Options Options_ { get; private set; }

        // The number of distinct shape ids that have been assigned.  This equals
        // the number of shapes in the index provided that no shapes have ever been
        // removed.  (Shape ids are not reused.)
        public override int NumShapeIds() => shapes_.Count;

        // Returns a pointer to the shape with the given id, or null if the shape
        // has been removed from the index.
        public override S2Shape Shape(int id) => shapes_[id];

        // Minimizes memory usage by requesting that any data structures that can be
        // rebuilt should be discarded.  This method invalidates all iterators.
        //
        // Like all non-methods, this method is not thread-safe.
        public override void Minimize()
        {
            // TODO(ericv): Implement.  In theory we should be able to discard the
            // entire index and rebuild it the next time it is needed.
        }

        // Appends an encoded representation of the S2ShapeIndex to "encoder".
        //
        // This method does not encode the S2Shapes in the index; it is the client's
        // responsibility to encode them separately.  For example:
        //
        //   S2ShapeUtil.CompactEncodeTaggedShapes(index, encoder);
        //   index.Encode(encoder);
        //
        // REQUIRES: "encoder" uses the default constructor, so that its buffer
        //           can be enlarged as necessary by calling Ensure(int).
        public void Encode(Encoder encoder)
        {
            // The version number is encoded in 2 bits, under the assumption that by the
            // time we need 5 versions the first version can be permanently retired.
            // This only saves 1 byte, but that's significant for very small indexes.
            encoder.Ensure(Encoder.kVarintMax64);
            var max_edges = (UInt64)Options_.MaxEdgesPerCell;
            encoder.PutVarUInt64(max_edges << 2 | kCurrentEncodingVersionNumber);

            var cell_ids = new List<S2CellId>(cell_map_.Count);
            var encoded_cells = new StringVectorEncoder();

            foreach (var it in GetNewEnumerable())
            {
                cell_ids.Add(it.Item1);
                it.Item2.Encode(NumShapeIds(), encoded_cells.AddViaEncoder());
            }
            EncodedS2CellIdVector.EncodeS2CellIdVector(cell_ids, encoder);
            encoded_cells.Encode(encoder);
        }

        // Decodes an S2ShapeIndex, returning true on success.
        //
        // This method does not decode the S2Shape objects in the index; this is
        // the responsibility of the client-provided function "shape_factory"
        // (see s2shapeutil_coding.h).  Example usage:
        //
        //   index.Init(decoder, S2ShapeUtil.LazyDecodeShapeFactory(decoder));
        //
        // Note that the S2Shape vector must be encoded *before* the S2ShapeIndex in
        // this example.
        public bool Init(Decoder decoder, ShapeFactory shape_factory)
        {
            Clear();
            if (!decoder.TryGetVarUInt64(out var max_edges_version)) return false;
            var version = (int)(max_edges_version & 3);
            if (version != kCurrentEncodingVersionNumber) return false;
            Options_.MaxEdgesPerCell = (int)(max_edges_version >> 2);
            var num_shapes = (UInt32)shape_factory.Count;
            shapes_.Capacity = (int)num_shapes;
            for (var shape_id = 0; shape_id < num_shapes; ++shape_id)
            {
                var shape = shape_factory[shape_id];
                if (shape != null) shape.SetId(shape_id);
                shapes_.Add(shape);
            }

            var cell_ids = new EncodedS2CellIdVector();
            if (!cell_ids.Init(decoder)) return false;

            if (!EncodedStringVector.Init(decoder, out var encoded_cells)) return false;

            for (var i = 0; i < cell_ids.Count; ++i)
            {
                var id = cell_ids[i];
                var cell = new S2ShapeIndexCell();
                var decoder2 = encoded_cells.GetDecoder(i);
                if (!cell.Decode((int)num_shapes, decoder2)) return false;
                // todo: maybe addsorted
                cell_map_.Add(new S2ShapeIndexIdCell(id, cell));
            }
            return true;
        }

        #region IEnumerator

        public override IReversableEnumerator<S2ShapeIndexIdCell> GetNewEnumerator()
        {
            MaybeApplyUpdates();

            return new ReversableEnumerator<S2ShapeIndexIdCell>(cell_map_);
        }

        public override IEnumerable<S2ShapeIndexIdCell> GetNewEnumerable()
        {
            MaybeApplyUpdates();

            return cell_map_;
        }
        public override int GetEnumerableCount()
        {
            MaybeApplyUpdates();

            return cell_map_.Count;
        }

        public override (int pos, bool found) SeekCell(S2CellId target)
        {
            MaybeApplyUpdates();

            var tup = new S2ShapeIndexIdCell(target, null);
            var pos = cell_map_.BinarySearch(tup);
            if (pos < 0) return (~pos, false);

            while (pos > 0 && tup.CompareTo(cell_map_[pos]) == 0) pos--;

            return (pos + 1, true);
        }

        public override S2CellId GetCellId(int pos)
        {
            // MaybeApplyUpdates has already been called

            if (pos >= cell_map_.Count || pos < 0) return null;

            return cell_map_[pos].Item1;
        }
        public override S2ShapeIndexIdCell? GetIndexCell(int pos)
        {
            if (pos >= cell_map_.Count || pos < 0) return null;

            return cell_map_[pos];
        }

        #endregion

        // Takes ownership of the given shape and adds it to the index.  Also
        // assigns a unique id to the shape (shape.id()) and returns that id.
        // Shape ids are assigned sequentially starting from 0 in the order shapes
        // are added.  Invalidates all iterators and their associated data.
        public int Add(S2Shape shape)
        {
            // Additions are processed lazily by ApplyUpdates().
            var id = NumShapeIds();
            shape.SetId(id);
            shapes_.Add(shape);
            IncChannel();
            return id;
        }

        // Removes the given shape from the index and return ownership to the caller.
        // Invalidates all iterators and their associated data.
        public S2Shape Release(int shape_id)
        {
            // This class updates itself lazily, because it is much more efficient to
            // process additions and removals in batches.  However this means that when
            // a shape is removed, we need to make a copy of all its edges, since the
            // client is free to delete "shape" once this call is finished.

            Assert.True(shapes_[shape_id] != null);
            var shape = shapes_[shape_id];
            if (shape_id < pending_additions_begin_)
            {
                var num_edges = shape.NumEdges;
                var edges = new List<S2Shape.Edge>(num_edges);
                for (var e = 0; e < num_edges; ++e)
                {
                    edges.Add(shape.GetEdge(e));
                }

                if (pending_removals_ == null)
                {
                    pending_removals_ = new List<RemovedShape>();
                }
                var removed = new RemovedShape(shape.Id, shape.Dimension() == 2,
                    shape.ContainsBruteForce(K_InteriorTrackerOrigin()), edges);
                pending_removals_.Add(removed);
            }
            /*else
            {
                // We are removing a shape that has not yet been added to the index,
                // so there is nothing else to do.
            }*/
            IncChannel();
            return shape;
        }

        // Defines the initial focus point of MutableS2ShapeIndex.InteriorTracker
        // (the start of the S2CellId space-filling curve).
        //
        // TODO(ericv): Move InteriorTracker here to avoid the need for this method.
        private static S2Point K_InteriorTrackerOrigin()
        {
            return S2Coords.FaceUVtoXYZ(0, -1, -1).Normalized;
        }

        // Resets the index to its original state and returns ownership of all
        // shapes to the caller.  This method is much more efficient than removing
        // all shapes one at a time.
        public void ReleaseAll()
        {
            cell_map_.Clear();
            pending_additions_begin_ = 0;
            if (pending_removals_ != null) pending_removals_.Clear();
            ResetChannel();
            shapes_.Clear();
        }

        // Resets the index to its original state and deletes all shapes.  Any
        // options specified via Init() are preserved.
        public void Clear()
        {
            ReleaseAll();
        }

        // Returns the number of bytes currently occupied by the index (including any
        // unused space at the end of vectors, etc). It has the same thread safety
        // as the other "const" methods (see introduction).
        public override int SpaceUsed()
        {
            var size = Marshal.SizeOf(this);
            size += NumShapeIds() * Marshal.SizeOf(typeof(S2Shape));
            // cell_map_ itself is already included in sizeof(*this).
            size += cell_map_.Count - Marshal.SizeOf(cell_map_);
            size += cell_map_.Count * Marshal.SizeOf(typeof(S2ShapeIndexCell));
            foreach (var it in GetNewEnumerable())
            {
                var cell = it.Item2;
                size += cell.NumClipped() * Marshal.SizeOf(typeof(S2ClippedShape));
                for (var s = 0; s < cell.NumClipped(); ++s)
                {
                    var clipped = cell.Clipped(s);
                    size += clipped.NumEdges* sizeof(Int32);
                }
            }
            if (pending_removals_ != null)
            {
                size += pending_removals_.Count * Marshal.SizeOf(typeof(RemovedShape));
            }

            return size;
        }

        // Calls to Add() and Release() are normally queued and processed on the
        // first subsequent query (in a thread-safe way).  This has many advantages,
        // the most important of which is that sometimes there *is* no subsequent
        // query, which lets us avoid building the index completely.
        //
        // This method forces any pending updates to be applied immediately.
        // Calling this method is rarely a good idea.  (One valid reason is to
        // exclude the cost of building the index from benchmark results.)
        // 
        // It is the client's responsibility to ensure correct thread synchronization.
        public void ForceBuild()
        {
            ApplyUpdatesInternal();
        }

        // Returns true if there are no pending updates that need to be applied.
        // This can be useful to avoid building the index unnecessarily, or for
        // choosing between two different algorithms depending on whether the index
        // is available.
        //
        // The returned index status may be slightly out of date if the index was
        // built in a different thread.  This is fine for the intended use (as an
        // efficiency hint), but it should not be used by internal methods  (see
        // MaybeApplyUpdates).
        public bool IsFresh()
        {
            lock (_channelLock)
            {
                return _channelPendingOperations == 0;
            }
        }
        public void IncChannel()
        {
            lock (_channelLock)
            {
                _channelPendingOperations++;
            }
        }
        public void ResetChannel()
        {
            lock (_channelLock)
            {
                _channelPendingOperations = 0;
            }
        }

        // Return true if this is the first update to the index.
        private bool IsFirstUpdate()
        {
            // Note that it is not sufficient to check whether cell_map_ is empty, since
            // entries are added during the update process.
            return pending_additions_begin_ == 0;
        }
        // Given that the given shape is being updated, return true if it is being
        // removed (as opposed to being added).
        private bool IsShapeBeingRemoved(int shape_id)
        {
            // All shape ids being removed are less than all shape ids being added.
            return shape_id < pending_additions_begin_;
        }
        
        // Ensure that any pending updates have been applied.  This method must be
        // called before accessing the cell_map_ field, even if the index_status_
        // appears to be FRESH, because a memory barrier is required in order to
        // ensure that all the index updates are visible if the updates were done in
        // another thread.
        private void MaybeApplyUpdates()
        {
            if (!IsFresh())
            {
                ApplyUpdatesThreadSafe();
            }
        }

        // Apply any pending updates in a thread-safe way.
        private void ApplyUpdatesThreadSafe()
        {
            lock (_channelLock)
            {
                ApplyUpdatesInternal();
            }
        }

        // This method updates the index by applying all pending additions and
        // removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
        private void ApplyUpdatesInternal()
        {
            // Check whether we have so many edges to process that we should process
            // them in multiple batches to save memory.  Building the index can use up
            // to 20x as much memory (per edge) as the final index size.
            var batches = GetUpdateBatches();
            foreach (var batch in batches)
            {
                var all_edges = new Array6<List<FaceEdge>>(() => new List<FaceEdge>());

                ReserveSpace(batch, all_edges);
                var tracker = new InteriorTracker();
                if (pending_removals_ != null)
                {
                    // The first batch implicitly includes all shapes being removed.
                    foreach (var pending_removal in pending_removals_)
                    {
                        RemoveShape(pending_removal, all_edges, tracker);
                    }
                    pending_removals_ = null;
                }
                for (var id = pending_additions_begin_; id < batch.AdditionsEnd; ++id)
                {
                    AddShape(id, all_edges, tracker);
                }
                for (var face = 0; face < 6; ++face)
                {
                    UpdateFaceEdges(face, all_edges[face], tracker);
                    // Save memory by clearing vectors after we are done with them.
                    all_edges[face].Clear();
                }
                pending_additions_begin_ = batch.AdditionsEnd;
            }
            ResetChannel();
        }
        // Count the number of edges being updated, and break them into several
        // batches if necessary to reduce the amount of memory needed.  (See the
        // documentation for FLAGS_s2shape_index_tmp_memory_budget_mb.)
        private List<BatchDescriptor> GetUpdateBatches()
        {
            var res = new List<BatchDescriptor>();
            // Count the edges being removed and added.
            var num_edges_removed = 0;
            if (pending_removals_ != null)
            {
                foreach (var pending_removal in pending_removals_)
                {
                    num_edges_removed += pending_removal.Edges.Count;
                }
            }
            var num_edges_added = 0;
            var shapesCount = NumShapeIds();
            for (var id = pending_additions_begin_; id < shapesCount; ++id)
            {
                var shape = Shape(id);
                if (shape == null) continue;
                num_edges_added += shape.NumEdges;
            }
            var num_edges = num_edges_removed + num_edges_added;

            // The following memory estimates are based on heap profiling.
            //
            // The final size of a MutableS2ShapeIndex depends mainly on how finely the
            // index is subdivided, as controlled by Options.max_edges_per_cell() and
            // --s2shape_index_default_max_edges_per_cell. For realistic values of
            // max_edges_per_cell() and shapes with moderate numbers of edges, it is
            // difficult to get much below 8 bytes per edge.  [The minimum possible size
            // is 4 bytes per edge (to store a 32-bit edge id in an S2ClippedShape) plus
            // 24 bytes per shape (for the S2ClippedShape itself plus a pointer in the
            // shapes_ vector.]
            //
            // The temporary memory consists mainly of the FaceEdge and ClippedEdge
            // structures plus a ClippedEdge pointer for every level of recursive
            // subdivision.  For very large indexes this can be 200 bytes per edge.
            const int kFinalBytesPerEdge = 8;
            const int kTmpBytesPerEdge = 200;
            const int kTmpMemoryBudgetBytes = s2shape_index_tmp_memory_budget_mb << 20;

            // We arbitrarily limit the number of batches just as a safety measure.
            // With the current default memory budget of 100 MB, this limit is not
            // reached even when building an index of 350 million edges.
            const int kMaxUpdateBatches = 100;

            if (num_edges * kTmpBytesPerEdge <= kTmpMemoryBudgetBytes)
            {
                // We can update all edges at once without exceeding kTmpMemoryBudgetBytes.
                res.Add(new BatchDescriptor { AdditionsEnd = shapesCount, NumEdges = num_edges });
                return res;
            }
            // Otherwise, break the updates into up to several batches, where the size
            // of each batch is chosen so that all batches use approximately the same
            // high-water memory.  GetBatchSizes() returns the recommended number of
            // edges in each batch.
            var batch_sizes = new List<int>();
            GetBatchSizes(num_edges, kMaxUpdateBatches, kFinalBytesPerEdge,
                          kTmpBytesPerEdge, kTmpMemoryBudgetBytes, batch_sizes);

            // We always process removed edges in a single batch, since (1) they already
            // take up a lot of memory because we have copied all their edges, and (2)
            // AbsorbIndexCell() uses (shapes_[id] == null) to detect when a shape is
            // being removed, so in order to split the removals into batches we would
            // need a different approach (e.g., temporarily add fake entries to shapes_
            // and restore them back to null as shapes are actually removed).
            num_edges = 0;
            if (pending_removals_ != null)
            {
                num_edges += num_edges_removed;
                if (num_edges >= batch_sizes[0])
                {
                    res.Add(new BatchDescriptor { AdditionsEnd = pending_additions_begin_, NumEdges = num_edges });
                    num_edges = 0;
                }
            }
            // Keep adding shapes to each batch until the recommended number of edges
            // for that batch is reached, then move on to the next batch.
            for (var id = pending_additions_begin_; id < shapesCount; ++id)
            {
                var shape = Shape(id);
                if (shape == null) continue;
                num_edges += shape.NumEdges;
                if (num_edges >= batch_sizes[res.Count])
                {
                    res.Add(new BatchDescriptor { AdditionsEnd = id + 1, NumEdges = num_edges });
                    num_edges = 0;
                }
            }
            // Some shapes have no edges.  If a shape with no edges is the last shape to
            // be added or removed, then the final batch may not include it, so we fix
            // that problem here.
            res[^1] = new BatchDescriptor { AdditionsEnd = shapesCount, NumEdges = res[^1].NumEdges };
            Assert.True(res.Count <= kMaxUpdateBatches);
            return res;
        }

        // Given "num_items" items, each of which uses "high_water_bytes_per_item" while it
        // is being updated but only "final_bytes_per_item" in the end, divide the
        // items into batches that have approximately the same *total* memory usage
        // consisting of the temporary memory needed for the items in the current
        // batch plus the final size of all the items that have already been
        // processed.  Use the fewest number of batches (but never more than
        // "max_batches") such that the total memory usage does not exceed the
        // combined final size of all the items plus "preferred_max_bytes_per_batch".
        private static void GetBatchSizes(int num_items, int max_batches, double final_bytes_per_item, double high_water_bytes_per_item, double preferred_max_bytes_per_batch, List<int> batch_sizes)
        {
            // This code tries to fit all the data into the same memory space
            // ("total_budget_bytes") at every iteration.  The data consists of some
            // number of processed items (at "final_bytes_per_item" each), plus some
            // number being updated (at "high_water_bytes_per_item" each).  The space occupied
            // by the items being updated is the "free space".  At each iteration, the
            // free space is multiplied by (1 - final_bytes_per_item/high_water_bytes_per_item)
            // as the items are converted into their final form.
            var final_bytes = num_items * final_bytes_per_item;
            var final_bytes_ratio = final_bytes_per_item / high_water_bytes_per_item;
            var free_space_multiplier = 1 - final_bytes_ratio;

            // The total memory budget is the greater of the final size plus the allowed
            // temporary memory, or the minimum amount of memory required to limit the
            // number of batches to "max_batches".
            var total_budget_bytes = Math.Max(
                final_bytes + preferred_max_bytes_per_batch,
                final_bytes / (1 - Math.Pow(free_space_multiplier, max_batches)));

            // "max_batch_items" is the number of items in the current batch.
            var max_batch_items = total_budget_bytes / high_water_bytes_per_item;
            batch_sizes.Clear();
            for (var i = 0; i + 1 < max_batches && num_items > 0; ++i)
            {
                var batch_items =
                    Math.Min(num_items, (int)(max_batch_items + 1));
                batch_sizes.Add(batch_items);
                num_items -= batch_items;
                max_batch_items *= free_space_multiplier;
            }
            Assert.True(batch_sizes.Count <= max_batches);
        }
        // Reserve an appropriate amount of space for the top-level face edges in the
        // current batch.  This data structure uses about half of the temporary memory
        // needed during index construction.  Furthermore, if the arrays are grown via
        // push_back() then up to 10% of the total run time consists of copying data
        // as these arrays grow, so it is worthwhile to preallocate space for them.
        private void ReserveSpace(BatchDescriptor batch, Array6<List<FaceEdge>> all_edges)
        {
            // If the number of edges is relatively small, then the fastest approach is
            // to simply reserve space on every face for the maximum possible number of
            // edges.  We use a different threshold for this calculation than for
            // deciding when to break updates into batches, because the cost/benefit
            // ratio is different.  (Here the only extra expense is that we need to
            // sample the edges to estimate how many edges per face there are.)
            const int kMaxCheapBytes = 30 << 20;  // 30 MB
            var kMaxCheapEdges = kMaxCheapBytes / (6 * Marshal.SizeOf(typeof(FaceEdge)));
            if (batch.NumEdges <= kMaxCheapEdges)
            {
                for (var face = 0; face < 6; ++face)
                {
                    all_edges[face] = new List<FaceEdge> { Capacity = batch.NumEdges };
                }
                return;
            }
            // Otherwise we estimate the number of edges on each face by taking a random
            // sample.  The goal is to come up with an estimate that is fast and
            // accurate for non-pathological geometry.  If our estimates happen to be
            // wrong, the vector will still grow automatically - the main side effects
            // are that memory usage will be larger (by up to a factor of 3), and
            // constructing the index will be about 10% slower.
            //
            // Given a desired sample size, we choose equally spaced edges from
            // throughout the entire data set.  We use a Bresenham-type algorithm to
            // choose the samples.
            var kDesiredSampleSize = 10000;
            var sample_interval = Math.Max(1, batch.NumEdges / kDesiredSampleSize);

            // Initialize "edge_id" to be midway through the first sample interval.
            // Because samples are equally spaced the actual sample size may differ
            // slightly from the desired sample size.
            var edge_id = sample_interval / 2;
            var actual_sample_size = (batch.NumEdges + edge_id) / sample_interval;
            var face_count = new int[6] { 0, 0, 0, 0, 0, 0 };
            if (pending_removals_ != null)
            {
                foreach (var removed in pending_removals_)
                {
                    edge_id += removed.Edges.Count;
                    while (edge_id >= sample_interval)
                    {
                        edge_id -= sample_interval;
                        face_count[S2Coords.GetFace(removed.Edges[edge_id].V0)] += 1;
                    }
                }
            }
            for (int id = pending_additions_begin_; id < batch.AdditionsEnd; ++id)
            {
                var shape = Shape(id);
                if (shape == null) continue;
                edge_id += shape.NumEdges;
                while (edge_id >= sample_interval)
                {
                    edge_id -= sample_interval;
                    // For speed, we only count the face containing one endpoint of the
                    // edge.  In general the edge could span all 6 faces (with padding), but
                    // it's not worth the expense to compute this more accurately.
                    face_count[S2Coords.GetFace(shape.GetEdge(edge_id).V0)] += 1;
                }
            }
            // Now given the raw face counts, compute a confidence interval such that we
            // will be unlikely to allocate too little space.  Computing accurate
            // binomial confidence intervals is expensive and not really necessary.
            // Instead we use a simple approximation:
            //  - For any face with at least 1 sample, we use at least a 4-sigma
            //    confidence interval.  (The chosen width is adequate for the worst case
            //    accuracy, which occurs when the face contains approximately 50% of the
            //    edges.)  Assuming that our sample is representative, the probability of
            //    reserving too little space is approximately 1 in 30,000.
            //  - For faces with no samples at all, we don't bother reserving space.
            //    It is quite likely that such faces are truly empty, so we save time
            //    and memory this way.  If the face does contain some edges, there will
            //    only be a few so it is fine to let the vector grow automatically.
            // On average, we reserve 2% extra space for each face that has geometry.

            // kMaxSemiWidth is the maximum semi-width over all probabilities p of a
            // 4-sigma binomial confidence interval with a sample size of 10,000.
            const double kMaxSemiWidth = 0.02;
            var sample_ratio = 1.0 / actual_sample_size;
            for (var face = 0; face < 6; ++face)
            {
                if (face_count[face] == 0) continue;
                var fraction = sample_ratio * face_count[face] + kMaxSemiWidth;
                all_edges[face].Capacity = 1 + (int)(fraction * batch.NumEdges);
            }
        }
        // Clip all edges of the given shape to the six cube faces, add the clipped
        // edges to "all_edges", and start tracking its interior if necessary.
        private void AddShape(int id, Array6<List<FaceEdge>> all_edges, InteriorTracker tracker)
        {
            var shape = Shape(id);
            if (shape == null)
            {
                return;  // This shape has already been removed.
            }
            var has_interior = shape.Dimension() == 2;
            if (has_interior)
            {
                tracker.AddShape(id, shape.ContainsBruteForce(tracker.Focus()));
            }
            var num_edges = shape.NumEdges;
            for (var e = 0; e < num_edges; ++e)
            {
                var edge = shape.GetEdge(e);
                var max_level = GetEdgeMaxLevel(edge);
                AddFaceEdge(all_edges, id, e, max_level, has_interior, edge);
            }
        }
        private static void RemoveShape(RemovedShape removed, Array6<List<FaceEdge>> all_edges, InteriorTracker tracker)
        {
            if (removed.HasInterior)
            {
                tracker.AddShape(removed.ShapeId, removed.ContainsTrackerOrigin);
            }
            foreach (var edge in removed.Edges)
            {
                var max_level = GetEdgeMaxLevel(edge);
                AddFaceEdge(all_edges, removed.ShapeId, -1, max_level, removed.HasInterior, edge);
            }
        }
        private static void AddFaceEdge(Array6<List<FaceEdge>> all_edges, Int32 shape_id, Int32 edge_id, Int32 max_level, bool has_interior, S2Shape.Edge edge)
        {
            R2Point edge_a, edge_b;
            // Fast path: both endpoints are on the same face, and are far enough from
            // the edge of the face that don't intersect any (padded) adjacent face.
            var a_face = S2Coords.GetFace(edge.V0);
            if (a_face == S2Coords.GetFace(edge.V1))
            {
                S2Coords.ValidFaceXYZtoUV(a_face, edge.V0, out edge_a);
                S2Coords.ValidFaceXYZtoUV(a_face, edge.V1, out edge_b);
                var kMaxUV = 1 - kCellPadding;
                if (Math.Abs(edge_a[0]) <= kMaxUV && Math.Abs(edge_a[1]) <= kMaxUV &&
                    Math.Abs(edge_b[0]) <= kMaxUV && Math.Abs(edge_b[1]) <= kMaxUV)
                {
                    var fe = new FaceEdge(shape_id, edge_id, max_level, has_interior, edge, edge_a, edge_b);
                    all_edges[a_face].Add(fe);
                    return;
                }
            }
            // Otherwise we simply clip the edge to all six faces.
            for (int face = 0; face < 6; ++face)
            {
                var padded = S2EdgeClipping.ClipToPaddedFace(edge.V0, edge.V1, face, kCellPadding, out edge_a, out edge_b);
                if (padded)
                {
                    var fe = new FaceEdge(shape_id, edge_id, max_level, has_interior, edge, edge_a, edge_b);
                    all_edges[face].Add(fe);
                }
            }
        }
        // Given a face and a vector of edges that intersect that face, add or remove
        // all the edges from the index.  (An edge is added if shapes_[id] is not
        // null, and removed otherwise.)
        private void UpdateFaceEdges(int face, List<FaceEdge> face_edges, InteriorTracker tracker)
        {
            var num_edges = face_edges.Count;
            if (num_edges == 0 && !tracker.ShapeIds().Any()) return;

            // Create the initial ClippedEdge for each FaceEdge.  Additional clipped
            // edges are created when edges are split between child cells.  We create
            // two arrays, one containing the edge data and another containing pointers
            // to those edges, so that during the recursion we only need to copy
            // pointers in order to propagate an edge to the correct child.
            var clipped_edge_storage = new List<ClippedEdge>(num_edges);
            var clipped_edges = new List<ClippedEdge>(num_edges);
            var bound = R2Rect.Empty;
            for (var e = 0; e < num_edges; ++e)
            {
                var clipped = new ClippedEdge(face_edges[e],
                    R2Rect.FromPointPair(face_edges[e].A, face_edges[e].B));
                clipped_edge_storage.Add(clipped);
                clipped_edges.Add(clipped_edge_storage.Last());
                bound = bound.AddRect(clipped.Bound);
            }
            // Construct the initial face cell containing all the edges, and then update
            // all the edges in the index recursively.
            var alloc = new EdgeAllocator();
            var face_id = S2CellId.FromFace(face);
            var pcell = new S2PaddedCell(face_id, kCellPadding);

            // "disjoint_from_index" means that the current cell being processed (and
            // all its descendants) are not already present in the index.
            bool disjoint_from_index = IsFirstUpdate();
            if (num_edges > 0)
            {
                var shrunk_id = ShrinkToFit(pcell, bound);
                if (shrunk_id != pcell.Id)
                {
                    // All the edges are contained by some descendant of the face cell.  We
                    // can save a lot of work by starting directly with that cell, but if we
                    // are in the interior of at least one shape then we need to create
                    // index entries for the cells we are skipping over.
                    SkipCellRange(face_id.RangeMin, shrunk_id.RangeMin,
                                  tracker, alloc, disjoint_from_index);
                    pcell = new S2PaddedCell(shrunk_id, kCellPadding);
                    UpdateEdges(pcell, clipped_edges, tracker, alloc, disjoint_from_index);
                    SkipCellRange(shrunk_id.RangeMax.Next, face_id.RangeMax.Next,
                                  tracker, alloc, disjoint_from_index);
                    return;
                }
            }
            // Otherwise (no edges, or no shrinking is possible), subdivide normally.
            UpdateEdges(pcell, clipped_edges, tracker, alloc, disjoint_from_index);
        }
        private S2CellId ShrinkToFit(S2PaddedCell pcell, R2Rect bound)
        {
            var shrunk_id = pcell.ShrinkToFit(bound);
            if (!IsFirstUpdate() && shrunk_id != pcell.Id)
            {
                // Don't shrink any smaller than the existing index cells, since we need
                // to combine the new edges with those cells.
                // Use InitStale() to avoid applying updated recursively.
                var (r, pos) = LocateCell(shrunk_id);
                if (r == CellRelation.INDEXED) { shrunk_id = GetCellId(pos); }
            }
            return shrunk_id;
        }
        // Skip over the cells in the given range, creating index cells if we are
        // currently in the interior of at least one shape.
        private void SkipCellRange(S2CellId begin, S2CellId end, InteriorTracker tracker, EdgeAllocator alloc, bool disjoint_from_index)
        {
            // If we aren't in the interior of a shape, then skipping over cells is easy.
            if (!tracker.ShapeIds().Any()) return;

            // Otherwise generate the list of cell ids that we need to visit, and create
            // an index entry for each one.
            foreach (var skipped_id in S2CellUnion.FromBeginEnd(begin, end))
            {
                var d = new List<ClippedEdge>();
                UpdateEdges(new S2PaddedCell(skipped_id, kCellPadding), d, tracker, alloc, disjoint_from_index);
            }
        }

        // Given a cell and a set of ClippedEdges whose bounding boxes intersect that
        // cell, add or remove all the edges from the index.  Temporary space for
        // edges that need to be subdivided is allocated from the given EdgeAllocator.
        // "disjoint_from_index" is an optimization hint indicating that cell_map_
        // does not contain any entries that overlap the given cell.
        private void UpdateEdges(S2PaddedCell pcell, List<ClippedEdge> edges, InteriorTracker tracker, EdgeAllocator alloc, bool disjoint_from_index)
        {
            // Cases where an index cell is not needed should be detected before this.
            Assert.True(edges.Any() || tracker.ShapeIds().Any());

            // This function is recursive with a maximum recursion depth of 30
            // (S2Constants.kMaxCellLevel).  Note that using an explicit stack does not seem
            // to be any faster based on profiling.

            // Incremental updates are handled as follows.  All edges being added or
            // removed are combined together in "edges", and all shapes with interiors
            // are tracked using "tracker".  We subdivide recursively as usual until we
            // encounter an existing index cell.  At this point we "absorb" the index
            // cell as follows:
            //
            //   - Edges and shapes that are being removed are deleted from "edges" and
            //     "tracker".
            //   - All remaining edges and shapes from the index cell are added to
            //     "edges" and "tracker".
            //   - Continue subdividing recursively, creating new index cells as needed.
            //   - When the recursion gets back to the cell that was absorbed, we
            //     restore "edges" and "tracker" to their previous state.
            //
            // Note that the only reason that we include removed shapes in the recursive
            // subdivision process is so that we can find all of the index cells that
            // contain those shapes efficiently, without maintaining an explicit list of
            // index cells for each shape (which would be expensive in terms of memory).
            bool index_cell_absorbed = false;
            if (!disjoint_from_index)
            {
                // There may be existing index cells contained inside "pcell".  If we
                // encounter such a cell, we need to combine the edges being updated with
                // the existing cell contents by "absorbing" the cell.
                // Use InitStale() to avoid applying updated recursively.
                var (r, pos) = LocateCell(pcell.Id);
                if (r == S2ShapeIndex.CellRelation.DISJOINT)
                {
                    disjoint_from_index = true;
                }
                else if (r == S2ShapeIndex.CellRelation.INDEXED)
                {
                    // Absorb the index cell by transferring its contents to "edges" and
                    // deleting it.  We also start tracking the interior of any new shapes.
                    AbsorbIndexCell(pcell, GetIndexCell(pos).Value, edges, tracker, alloc);
                    index_cell_absorbed = true;
                    disjoint_from_index = true;
                }
                else
                {
                    Assert.True(S2ShapeIndex.CellRelation.SUBDIVIDED == r);
                }
            }

            // If there are existing index cells below us, then we need to keep
            // subdividing so that we can merge with those cells.  Otherwise,
            // MakeIndexCell checks if the number of edges is small enough, and creates
            // an index cell if possible (returning true when it does so).
            if (!disjoint_from_index || !MakeIndexCell(pcell, edges, tracker))
            {
                // Reserve space for the edges that will be passed to each child.  This is
                // important since otherwise the running time is dominated by the time
                // required to grow the vectors.  The amount of memory involved is
                // relatively small, so we simply reserve the maximum space for every child.
                int num_edges = edges.Count;
                var child_edges = new List<ClippedEdge>[][]   // [i][j]
                {
                   new [] { new List<ClippedEdge>(num_edges), new List<ClippedEdge>(num_edges) },
                   new [] { new List<ClippedEdge>(num_edges), new List<ClippedEdge>(num_edges) },
                };

                // Remember the current size of the EdgeAllocator so that we can free any
                // edges that are allocated during edge splitting.
                var alloc_size = alloc.Count();

                // Compute the middle of the padded cell, defined as the rectangle in
                // (u,v)-space that belongs to all four (padded) children.  By comparing
                // against the four boundaries of "middle" we can determine which children
                // each edge needs to be propagated to.
                var middle = pcell.Middle;

                // Build up a vector edges to be passed to each child cell.  The (i,j)
                // directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
                // Note that the vast majority of edges are propagated to a single child.
                // This case is very fast, consisting of between 2 and 4 floating-point
                // comparisons and copying one pointer.  (ClipVAxis is inline.)
                for (var e = 0; e < num_edges; ++e)
                {
                    var edge = edges[e];
                    if (edge.Bound[0].Hi <= middle[0].Lo)
                    {
                        // Edge is entirely contained in the two left children.
                        ClipVAxis(edge, middle[1], child_edges[0], alloc);
                    }
                    else if (edge.Bound[0].Lo >= middle[0].Hi)
                    {
                        // Edge is entirely contained in the two right children.
                        ClipVAxis(edge, middle[1], child_edges[1], alloc);
                    }
                    else if (edge.Bound[1].Hi <= middle[1].Lo)
                    {
                        // Edge is entirely contained in the two lower children.
                        child_edges[0][0].Add(ClipUBound(edge, 1, middle[0].Hi, alloc));
                        child_edges[1][0].Add(ClipUBound(edge, 0, middle[0].Lo, alloc));
                    }
                    else if (edge.Bound[1].Lo >= middle[1].Hi)
                    {
                        // Edge is entirely contained in the two upper children.
                        child_edges[0][1].Add(ClipUBound(edge, 1, middle[0].Hi, alloc));
                        child_edges[1][1].Add(ClipUBound(edge, 0, middle[0].Lo, alloc));
                    }
                    else
                    {
                        // The edge bound spans all four children.  The edge itself intersects
                        // either three or four (padded) children.
                        var left = ClipUBound(edge, 1, middle[0].Hi, alloc);
                        ClipVAxis(left, middle[1], child_edges[0], alloc);
                        var right = ClipUBound(edge, 0, middle[0].Lo, alloc);
                        ClipVAxis(right, middle[1], child_edges[1], alloc);
                    }
                }
                // Now recursively update the edges in each child.  We call the children in
                // increasing order of S2CellId so that when the index is first constructed,
                // all insertions into cell_map_ are at the end (which is much faster).
                for (var pos = 0; pos < 4; ++pos)
                {
                    pcell.GetChildIJ(pos, out int i, out int j);
                    if (child_edges[i][j].Any() || tracker.ShapeIds().Any())
                    {
                        UpdateEdges(new S2PaddedCell(pcell, i, j), child_edges[i][j],
                                    tracker, alloc, disjoint_from_index);
                    }
                }
                // Free any temporary edges that were allocated during clipping.
                alloc.Reset(alloc_size);
            }
            if (index_cell_absorbed)
            {
                // Restore the state for any edges being removed that we are tracking.
                tracker.RestoreStateBefore(pending_additions_begin_);
            }
        }
        // Absorb an index cell by transferring its contents to "edges" and/or
        // "tracker", and then delete this cell from the index.  If "edges" includes
        // any edges that are being removed, this method also updates their
        // InteriorTracker state to correspond to the exit vertex of this cell, and
        // saves the InteriorTracker state by calling SaveAndClearStateBefore().  It
        // is the caller's responsibility to restore this state by calling
        // RestoreStateBefore() when processing of this cell is finished.
        private void AbsorbIndexCell(S2PaddedCell pcell, S2ShapeIndexIdCell item, List<ClippedEdge> edges, InteriorTracker tracker, EdgeAllocator alloc)
        {
            Assert.True(pcell.Id== item.Item1);

            // When we absorb a cell, we erase all the edges that are being removed.
            // However when we are finished with this cell, we want to restore the state
            // of those edges (since that is how we find all the index cells that need
            // to be updated).  The edges themselves are restored automatically when
            // UpdateEdges returns from its recursive call, but the InteriorTracker
            // state needs to be restored explicitly.
            //
            // Here we first update the InteriorTracker state for removed edges to
            // correspond to the exit vertex of this cell, and then save the
            // InteriorTracker state.  This state will be restored by UpdateEdges when
            // it is finished processing the contents of this cell.
            if (tracker.IsActive&& edges.Any() &&
                IsShapeBeingRemoved(edges[0].FaceEdge.ShapeId))
            {
                // We probably need to update the InteriorTracker.  ("Probably" because
                // it's possible that all shapes being removed do not have interiors.)
                if (!tracker.AtCellId(pcell.Id))
                {
                    tracker.MoveTo(pcell.GetEntryVertex());
                }
                tracker.DrawTo(pcell.GetExitVertex());
                tracker.NextCellId = pcell.Id.Next;
                foreach (var edge in edges)
                {
                    var face_edge = edge.FaceEdge;
                    if (!IsShapeBeingRemoved(face_edge.ShapeId))
                    {
                        break;  // All shapes being removed come first.
                    }
                    if (face_edge.HasInterior)
                    {
                        tracker.TestEdge(face_edge.ShapeId, face_edge.Edge);
                    }
                }
            }
            // Save the state of the edges being removed, so that it can be restored
            // when we are finished processing this cell and its children.  We don't
            // need to save the state of the edges being added because they aren't being
            // removed from "edges" and will therefore be updated normally as we visit
            // this cell and its children.
            tracker.SaveAndClearStateBefore(pending_additions_begin_);

            // Create a FaceEdge for each edge in this cell that isn't being removed.
            var face_edges = alloc.FaceEdges;
            face_edges.Clear();
            bool tracker_moved = false;
            var cell = item.Item2;
            for (var s = 0; s < cell.NumClipped(); ++s)
            {
                var clipped = cell.Clipped(s);
                int shape_id = clipped.ShapeId;
                var shape = Shape(shape_id);
                if (shape == null) continue;  // This shape is being removed.
                int num_edges = clipped.NumEdges;

                // If this shape has an interior, start tracking whether we are inside the
                // shape.  UpdateEdges() wants to know whether the entry vertex of this
                // cell is inside the shape, but we only know whether the center of the
                // cell is inside the shape, so we need to test all the edges against the
                // line segment from the cell center to the entry vertex.
                var edge_shape_id = shape.Id;
                var edge_has_interior = shape.Dimension() == 2;
                if (edge_has_interior)
                {
                    tracker.AddShape(shape_id, clipped.ContainsCenter);
                    // There might not be any edges in this entire cell (i.e., it might be
                    // in the interior of all shapes), so we delay updating the tracker
                    // until we see the first edge.
                    if (!tracker_moved && num_edges > 0)
                    {
                        tracker.MoveTo(pcell.GetCenter());
                        tracker.DrawTo(pcell.GetEntryVertex());
                        tracker.NextCellId = pcell.Id;
                        tracker_moved = true;
                    }
                }

                for (int i = 0; i < num_edges; ++i)
                {
                    int e = clipped.Edge(i);
                    var edge_edge_id = e;
                    var edge_edge = shape.GetEdge(e);
                    var edge_max_level = GetEdgeMaxLevel(edge_edge);
                    if (edge_has_interior) tracker.TestEdge(shape_id, edge_edge);
                    var face = (int)pcell.Id.Face;
                    if (!S2EdgeClipping.ClipToPaddedFace(edge_edge.V0, edge_edge.V1, face, kCellPadding, out R2Point a, out R2Point b))
                    {
                        throw new ApplicationException("Invariant failure in MutableS2ShapeIndex");
                    }
                    face_edges.Add(new FaceEdge(edge_shape_id, edge_edge_id, edge_max_level, edge_has_interior, edge_edge, a, b));
                }
            }
            // Now create a ClippedEdge for each FaceEdge, and put them in "new_edges".
            var new_edges = new List<ClippedEdge>();
            foreach (var face_edge in face_edges)
            {
                var bound = S2EdgeClipping.GetClippedEdgeBound(face_edge.A, face_edge.B, pcell.Bound);
                var clipped = new ClippedEdge(face_edge, bound);
                alloc.AddClippedEdge(clipped);
                new_edges.Add(clipped);
            }
            // Discard any edges from "edges" that are being removed, and append the
            // remainder to "new_edges".  (This keeps the edges sorted by shape id.)
            for (int i = 0; i < edges.Count; ++i)
            {
                var clipped = edges[i];
                if (!IsShapeBeingRemoved(clipped.FaceEdge.ShapeId))
                {
                    new_edges.AddRange(edges.Skip(i));
                    break;
                }
            }
            // Update the edge list and delete this cell from the index.
            edges.Clear();
            edges.AddRange(new_edges);
            // todo: if the list is sorted get the index with binarysearch and use RemoveAt
            cell_map_.RemoveAll(t => t.Item1 == pcell.Id);
        }

        // Return the first level at which the edge will *not* contribute towards
        // the decision to subdivide.
        public static int GetEdgeMaxLevel(S2Shape.Edge edge)
        {
            // Compute the maximum cell size for which this edge is considered "long".
            // The calculation does not need to be perfectly accurate, so we use Norm()
            // rather than Angle() for speed.
            double cell_size = (edge.V0 - edge.V1).Norm * s2shape_index_cell_into_long_edge_ratio;
            // Now return the first level encountered during subdivision where the
            // average cell size is at most "cell_size".
            return S2Metrics.kAvgEdge.GetLevelForMaxValue(cell_size);
        }
        // Return the number of distinct shapes that are either associated with the
        // given edges, or that are currently stored in the InteriorTracker.
        private static int CountShapes(List<ClippedEdge> edges, List<int> cshape_ids)
        {
            int count = 0;
            int last_shape_id = -1;
            var cnext = 0;  // Next shape
            foreach (var edge in edges)
            {
                if (edge.FaceEdge.ShapeId != last_shape_id)
                {
                    ++count;
                    last_shape_id = edge.FaceEdge.ShapeId;
                    // Skip over any containing shapes up to and including this one,
                    // updating "count" appropriately.
                    for (; cnext != cshape_ids.Count; ++cnext)
                    {
                        if (cshape_ids[cnext] > last_shape_id) break;
                        if (cshape_ids[cnext] < last_shape_id) ++count;
                    }
                }
            }
            // Count any remaining containing shapes.
            count += cshape_ids.Count - cnext;
            return count;
        }
        // Attempt to build an index cell containing the given edges, and return true
        // if successful.  (Otherwise the edges should be subdivided further.)
        private bool MakeIndexCell(S2PaddedCell pcell, List<ClippedEdge> edges, InteriorTracker tracker)
        {
            if (!edges.Any() && !tracker.ShapeIds().Any())
            {
                // No index cell is needed.  (In most cases this situation is detected
                // before we get to this point, but this can happen when all shapes in a
                // cell are removed.)
                return true;
            }

            // Count the number of edges that have not reached their maximum level yet.
            // Return false if there are too many such edges.
            int count = 0;
            foreach (ClippedEdge edge in edges)
            {
                count += (pcell.Level< edge.FaceEdge.MaxLevel) ? 1 : 0;
                if (count > Options_.MaxEdgesPerCell)
                    return false;
            }

            // Possible optimization: Continue subdividing as long as exactly one child
            // of "pcell" intersects the given edges.  This can be done by finding the
            // bounding box of all the edges and calling ShrinkToFit():
            //
            // S2CellId cellid = pcell.ShrinkToFit(GetRectBound(edges));
            //
            // Currently this is not beneficial; it slows down construction by 4-25%
            // (mainly computing the union of the bounding rectangles) and also slows
            // down queries (since more recursive clipping is required to get down to
            // the level of a spatial index cell).  But it may be worth trying again
            // once "contains_center" is computed and all algorithms are modified to
            // take advantage of it.

            // We update the InteriorTracker as follows.  For every S2Cell in the index
            // we construct two edges: one edge from entry vertex of the cell to its
            // center, and one from the cell center to its exit vertex.  Here "entry"
            // and "exit" refer the S2CellId ordering, i.e. the order in which points
            // are encountered along the S2 space-filling curve.  The exit vertex then
            // becomes the entry vertex for the next cell in the index, unless there are
            // one or more empty intervening cells, in which case the InteriorTracker
            // state is unchanged because the intervening cells have no edges.

            // Shift the InteriorTracker focus point to the center of the current cell.
            if (tracker.IsActive&& edges.Any())
            {
                if (!tracker.AtCellId(pcell.Id))
                {
                    tracker.MoveTo(pcell.GetEntryVertex());
                }
                tracker.DrawTo(pcell.GetCenter());
                TestAllEdges(edges, tracker);
            }
            // Allocate and fill a new index cell.  To get the total number of shapes we
            // need to merge the shapes associated with the intersecting edges together
            // with the shapes that happen to contain the cell center.
            var cshape_ids = tracker.ShapeIds();
            int num_shapes = CountShapes(edges, cshape_ids);
            var cell = new S2ShapeIndexCell();
            cell.AddCapacity(num_shapes);

            // To fill the index cell we merge the two sources of shapes: "edge shapes"
            // (those that have at least one edge that intersects this cell), and
            // "containing shapes" (those that contain the cell center).  We keep track
            // of the index of the next intersecting edge and the next containing shape
            // as we go along.  Both sets of shape ids are already sorted.
            int enext = 0;
            using var cnext = cshape_ids.GetEnumerator();
            var cnextHasNext = cnext.MoveNext();
            for (int i = 0; i < num_shapes; ++i)
            {
                var clipped = new S2ClippedShape();
                cell.AddShape(clipped);
                int eshape_id = NumShapeIds(), cshape_id = eshape_id;  // Sentinels
                if (enext != edges.Count)
                {
                    eshape_id = edges[enext].FaceEdge.ShapeId;
                }
                if (cnextHasNext)
                {
                    cshape_id = cnext.Current;
                }
                int ebegin = enext;
                if (cshape_id < eshape_id)
                {
                    // The entire cell is in the shape interior.
                    clipped.Init(cshape_id, 0);
                    clipped.ContainsCenter = (true);
                    cnextHasNext = cnext.MoveNext();
                }
                else
                {
                    // Count the number of edges for this shape and allocate space for them.
                    while (enext < edges.Count &&
                           edges[enext].FaceEdge.ShapeId == eshape_id)
                    {
                        ++enext;
                    }
                    clipped.Init(eshape_id, enext - ebegin);
                    for (int e = ebegin; e < enext; ++e)
                    {
                        clipped.SetEdge(e - ebegin, edges[e].FaceEdge.EdgeId);
                    }
                    if (cshape_id == eshape_id)
                    {
                        clipped.ContainsCenter = (true);
                        cnextHasNext = cnext.MoveNext();
                    }
                }
            }
            // UpdateEdges() visits cells in increasing order of S2CellId, so during
            // initial construction of the index all insertions happen at the end.  It
            // is much faster to give an insertion hint in this case.  Otherwise the
            // hint doesn't do much harm.  With more effort we could provide a hint even
            // during incremental updates, but this is probably not worth the effort.
            cell_map_.Add(new S2ShapeIndexIdCell(pcell.Id, cell));
            // todo: maybe add sorted

            // Shift the InteriorTracker focus point to the exit vertex of this cell.
            if (tracker.IsActive&& edges.Any())
            {
                tracker.DrawTo(pcell.GetExitVertex());
                TestAllEdges(edges, tracker);
                tracker.NextCellId = pcell.Id.Next;
            }
            return true;
        }

        // Call tracker.TestEdge() on all edges from shapes that have interiors.
        private static void TestAllEdges(List<ClippedEdge> edges, InteriorTracker tracker)
        {
            foreach (var edge in edges)
            {
                var face_edge = edge.FaceEdge;
                if (face_edge.HasInterior)
                {
                    tracker.TestEdge(face_edge.ShapeId, face_edge.Edge);
                }
            }
        }

        // Given an edge and two bound endpoints that need to be updated, allocate and
        // return a new edge with the updated bound.
        private static ClippedEdge UpdateBound(ClippedEdge edge, int u_end, double u, int v_end, double v, EdgeAllocator alloc)
        {
            var d = new double[2];
            d[u_end] = u;
            d[1 - u_end] = edge.Bound[0][1 - u_end];
            var x = new R1Interval(d);
            
            d[v_end] = v;
            d[1 - v_end] = edge.Bound[1][1 - v_end];
            var y = new R1Interval(d);
            
            var res = new ClippedEdge(edge.FaceEdge, new R2Rect(x, y));
            Assert.True(!res.Bound.IsEmpty);
            Assert.True(edge.Bound.Contains(res.Bound));
            alloc.AddClippedEdge(res);
            return res;
        }

        // Given an edge, clip the given endpoint (lo=0, hi=1) of the u-axis so that
        // it does not extend past the given value.
        private static ClippedEdge ClipUBound(ClippedEdge edge, int u_end, double u, EdgeAllocator alloc)
        {
            // First check whether the edge actually requires any clipping.  (Sometimes
            // this method is called when clipping is not necessary, e.g. when one edge
            // endpoint is in the overlap area between two padded child cells.)
            if (u_end == 0)
            {
                if (edge.Bound[0].Lo >= u) return edge;
            }
            else
            {
                if (edge.Bound[0].Hi <= u) return edge;
            }
            // We interpolate the new v-value from the endpoints of the original edge.
            // This has two advantages: (1) we don't need to store the clipped endpoints
            // at all, just their bounding box; and (2) it avoids the accumulation of
            // roundoff errors due to repeated interpolations.  The result needs to be
            // clamped to ensure that it is in the appropriate range.
            var e = edge.FaceEdge;
            double v = edge.Bound[1].Project(
                S2EdgeClipping.InterpolateDouble(u, e.A[0], e.B[0], e.A[1], e.B[1]));

            // Determine which endpoint of the v-axis bound to update.  If the edge
            // slope is positive we update the same endpoint, otherwise we update the
            // opposite endpoint.
            int v_end = u_end ^ (((e.A[0] > e.B[0]) != (e.A[1] > e.B[1])) ? 1 : 0);
            return UpdateBound(edge, u_end, u, v_end, v, alloc);
        }
        // Given an edge, clip the given endpoint (lo=0, hi=1) of the v-axis so that
        // it does not extend past the given value.
        private static ClippedEdge ClipVBound(ClippedEdge edge, int v_end, double v, EdgeAllocator alloc)
        {
            // See comments in ClipUBound.
            if (v_end == 0)
            {
                if (edge.Bound[1].Lo >= v) return edge;
            }
            else
            {
                if (edge.Bound[1].Hi <= v) return edge;
            }
            var e = edge.FaceEdge;
            double u = edge.Bound[0].Project(
                S2EdgeClipping.InterpolateDouble(v, e.A[1], e.B[1], e.A[0], e.B[0]));
            int u_end = v_end ^ (((e.A[0] > e.B[0]) != (e.A[1] > e.B[1])) ? 1 : 0);
            return UpdateBound(edge, u_end, u, v_end, v, alloc);
        }

        // Given an edge and an interval "middle" along the v-axis, clip the edge
        // against the boundaries of "middle" and add the edge to the corresponding
        // children.
        private static void ClipVAxis(ClippedEdge edge, R1Interval middle, List<ClippedEdge>[]/*[2]*/ child_edges, EdgeAllocator alloc)
        {
            if (edge.Bound[1].Hi <= middle.Lo)
            {
                // Edge is entirely contained in the lower child.
                child_edges[0].Add(edge);
            }
            else if (edge.Bound[1].Lo >= middle.Hi)
            {
                // Edge is entirely contained in the upper child.
                child_edges[1].Add(edge);
            }
            else
            {
                // The edge bound spans both children.
                child_edges[0].Add(ClipVBound(edge, 1, middle.Hi, alloc));
                child_edges[1].Add(ClipVBound(edge, 0, middle.Lo, alloc));
            }
        }

        // The shapes in the index, accessed by their shape id.  Removed shapes are
        // replaced by null pointers.
        private readonly List<S2Shape> shapes_ = new();

        // A map from S2CellId to the set of clipped shapes that intersect that
        // cell.  The cell ids cover a set of non-overlapping regions on the
        // sphere.  Note that this field is updated lazily (see below).  Const
        // methods *must* call MaybeApplyUpdates() before accessing this field.
        // (The easiest way to achieve this is simply to use an Iterator.)
        private readonly List<S2ShapeIndexIdCell> cell_map_ = new(); // gtl.btree_map

        // The id of the first shape that has been queued for addition but not
        // processed yet.
        private int pending_additions_begin_ = 0;

        // The representation of an edge that has been queued for removal.
        private record RemovedShape(
            int ShapeId,
            bool HasInterior, //Belongs to a shape of dimension 2.
            bool ContainsTrackerOrigin,
            List<S2Shape.Edge> Edges);

        // The set of shapes that have been queued for removal but not processed
        // yet.  Note that we need to copy the edge data since the caller is free to
        // destroy the shape once Release() has been called.  This field is present
        // only when there are removed shapes to process (to save memory).
        private List<RemovedShape> pending_removals_;

        // Additions and removals are queued and processed on the first subsequent
        // query.  There are several reasons to do this:
        //
        //  - It is significantly more efficient to process updates in batches.
        //  - Often the index will never be queried, in which case we can save both
        //    the time and memory required to build it.  Examples:
        //     + S2Loops that are created simply to pass to an S2Polygon.  (We don't
        //       need the S2Loop index, because S2Polygon builds its own index.)
        //     + Applications that load a database of geometry and then query only
        //       a small fraction of it.
        //     + Applications that only read and write geometry (Decode/Encode).
        //
        // The main drawback is that we need to go to some extra work to ensure that
        // "const" methods are still thread-safe.  Note that the goal is *not* to
        // make this class thread-safe in general, but simply to hide the fact that
        // we defer some of the indexing work until query time.
        private int _channelPendingOperations;
        private object _channelLock;

        private void InitChannel()
        {
            _channelPendingOperations = 0;
            _channelLock = new object();
        }

        // FaceEdge and ClippedEdge store temporary edge data while the index is being
        // updated.  FaceEdge represents an edge that has been projected onto a given
        // face, while ClippedEdge represents the portion of that edge that has been
        // clipped to a given S2Cell.
        //
        // While it would be possible to combine all the edge information into one
        // structure, there are two good reasons for separating it:
        //
        //  - Memory usage.  Separating the two classes means that we only need to
        //    store one copy of the per-face data no matter how many times an edge is
        //    subdivided, and it also lets us delay computing bounding boxes until
        //    they are needed for processing each face (when the dataset spans
        //    multiple faces).
        //
        //  - Performance.  UpdateEdges is significantly faster on large polygons when
        //    the data is separated, because it often only needs to access the data in
        //    ClippedEdge and this data is cached more successfully.
        [StructLayout(LayoutKind.Sequential)]
        public record FaceEdge
        (
            Int32 ShapeId,      // The shape that this edge belongs to
            Int32 EdgeId,       // Edge id within that shape
            Int32 MaxLevel,     // Not desirable to subdivide this edge beyond this level
            bool HasInterior,   // Belongs to a shape of dimension 2.
            S2Shape.Edge Edge,   // The edge endpoints
            R2Point A, R2Point B // The edge endpoints, clipped to a given face
        );

        public record ClippedEdge
        (
            FaceEdge FaceEdge,  // The original unclipped edge
            R2Rect Bound        // Bounding box for the clipped portion
        );

        // Given a set of shapes, InteriorTracker keeps track of which shapes contain
        // a particular point (the "focus").  It provides an efficient way to move the
        // focus from one point to another and incrementally update the set of shapes
        // which contain it.  We use this to compute which shapes contain the center
        // of every S2CellId in the index, by advancing the focus from one cell center
        // to the next.
        //
        // Initially the focus is at the start of the S2CellId space-filling curve.
        // We then visit all the cells that are being added to the MutableS2ShapeIndex
        // in increasing order of S2CellId.  For each cell, we draw two edges: one
        // from the entry vertex to the center, and another from the center to the
        // exit vertex (where "entry" and "exit" refer to the points where the
        // space-filling curve enters and exits the cell).  By counting edge crossings
        // we can incrementally compute which shapes contain the cell center.  Note
        // that the same set of shapes will always contain the exit point of one cell
        // and the entry point of the next cell in the index, because either (a) these
        // two points are actually the same, or (b) the intervening cells in S2CellId
        // order are all empty, and therefore there are no edge crossings if we follow
        // this path from one cell to the other.
        public class InteriorTracker
        {
            // Constructs the InteriorTracker.  You must call AddShape() for each shape
            // that will be tracked before calling MoveTo() or DrawTo().
            //
            // As shapes are added, we compute which ones contain the start of the
            // S2CellId space-filling curve by drawing an edge from S2PointUtil.Origin to this
            // point and counting how many shape edges cross this edge.
            public InteriorTracker()
            {
                IsActive = false; b_ = Origin();
                next_cellid_ = S2CellId.Begin(S2Constants.kMaxCellLevel);
            }

            // Returns the initial focus point when the InteriorTracker is constructed
            // (corresponding to the start of the S2CellId space-filling curve).
            public static S2Point Origin()
            {
                // The start of the S2CellId space-filling curve.
                return S2Coords.FaceUVtoXYZ(0, -1, -1).Normalized;
            }

            // Returns the current focus point (see above).
            public S2Point Focus() { return b_; }

            // Returns true if any shapes are being tracked.
            public bool IsActive { get; private set; }
            // Adds a shape whose interior should be tracked.  "is_inside" indicates
            // whether the current focus point is inside the shape.  Alternatively, if
            // the focus point is in the process of being moved (via MoveTo/DrawTo), you
            // can also specify "is_inside" at the old focus point and call TestEdge()
            // for every edge of the shape that might cross the current DrawTo() line.
            // This updates the state to correspond to the new focus point.
            //
            // REQUIRES: shape.dimension() == 2
            public void AddShape(Int32 shape_id, bool is_inside)
            {
                IsActive = true;
                if (is_inside)
                {
                    ToggleShape(shape_id);
                }
            }

            // Moves the focus to the given point.  This method should only be used when
            // it is known that there are no edge crossings between the old and new
            // focus locations; otherwise use DrawTo().
            public void MoveTo(S2Point b) { b_ = b; }

            // Moves the focus to the given point.  After this method is called,
            // TestEdge() should be called with all edges that may cross the line
            // segment between the old and new focus locations.
            public void DrawTo(S2Point b)
            {
                a_ = b_;
                b_ = b;
                crosser_.Init(a_, b_);
            }

            // Indicates that the given edge of the given shape may cross the line
            // segment between the old and new focus locations (see DrawTo).
            // REQUIRES: shape.dimension() == 2
            public void TestEdge(Int32 shape_id, S2Shape.Edge edge)
            {
                if (crosser_.EdgeOrVertexCrossing(edge.V0, edge.V1))
                {
                    ToggleShape(shape_id);
                }
            }

            // The set of shape ids that contain the current focus.
            public List<int> ShapeIds() { return shape_ids_; }

            // Indicates that the last argument to MoveTo() or DrawTo() was the entry
            // vertex of the given S2CellId, i.e. the tracker is positioned at the start
            // of this cell.  By using this method together with at_cellid(), the caller
            // can avoid calling MoveTo() in cases where the exit vertex of the previous
            // cell is the same as the entry vertex of the current cell.
            public S2CellId NextCellId
            {
                get => next_cellid_;
                set => next_cellid_ = value.RangeMin;
            }
            private S2CellId next_cellid_;

            // Returns true if the focus is already at the entry vertex of the given
            // S2CellId (provided that the caller calls set_next_cellid() as each cell
            // is processed).
            public bool AtCellId(S2CellId cellid) => cellid.RangeMin == next_cellid_;

            // Makes an internal copy of the state for shape ids below the given limit,
            // and then clear the state for those shapes.  This is used during
            // incremental updates to track the state of added and removed shapes
            // separately.
            public void SaveAndClearStateBefore(Int32 limit_shape_id)
            {
                Assert.True(!saved_ids_.Any());
                var limit = LowerBound(limit_shape_id);
                saved_ids_.AddRange(shape_ids_.Take(limit));
                shape_ids_.RemoveRange(0, limit);
            }

            // Restores the state previously saved by SaveAndClearStateBefore().  This
            // only affects the state for shape_ids below "limit_shape_id".
            public void RestoreStateBefore(Int32 limit_shape_id)
            {
                shape_ids_.RemoveRange(0, LowerBound(limit_shape_id));
                shape_ids_.InsertRange(0, saved_ids_);
                saved_ids_.Clear();
            }

            // Removes "shape_id" from shape_ids_ if it exists, otherwise insert it.
            private void ToggleShape(int shape_id)
            {
                // Since shape_ids_.size() is typically *very* small (0, 1, or 2), it turns
                // out to be significantly faster to maintain a sorted array rather than
                // using an STL set or btree_set.
                if (!shape_ids_.Any())
                {
                    shape_ids_.Add(shape_id);
                }
                else if (shape_ids_[0] == shape_id)
                {
                    shape_ids_.RemoveAt(0);
                }
                else
                {
                    int pos = 0;
                    while (pos < shape_id)
                    {
                        if (++pos == shape_ids_.Count)
                        {
                            shape_ids_.Add(shape_id);
                            return;
                        }
                    }
                    if (pos == shape_id)
                    {
                        shape_ids_.RemoveAt(pos);
                    }
                    else
                    {
                        shape_ids_.Insert(pos, shape_id);
                    }
                }
            }

            // Returns a pointer to the first entry "x" where x >= shape_id.
            //
            // Like lower_bound(shape_ids_.begin(), shape_ids_.end(), shape_id), but
            // implemented with linear rather than binary search because the number of
            // shapes being tracked is typically very small.
            private int LowerBound(Int32 shape_id)
            {
                var pos = 0;
                while (pos != shape_ids_.Count && shape_ids_[pos] < shape_id) { ++pos; }
                return pos;
            }

            private S2Point a_, b_;
            private readonly S2EdgeCrosser crosser_ = new S2EdgeCrosser();
            private readonly List<int> shape_ids_ = new List<int>();

            // Shape ids saved by SaveAndClearStateBefore().  The state is never saved
            // recursively so we don't need to worry about maintaining a stack.
            private readonly List<int> saved_ids_ = new List<int>();
        }

        // A BatchDescriptor represents a set of pending updates that will be applied
        // at the same time.  The batch consists of all updates with shape ids between
        // the current value of "ShapeIndex.pending_additions_begin_" (inclusive) and
        // "additions_end" (exclusive).  The first batch to be processed also
        // implicitly includes all shapes being removed.  "num_edges" is the total
        // number of edges that will be added or removed in this batch.
        private readonly struct BatchDescriptor
        {
            public readonly int AdditionsEnd { get; init; }
            public readonly int NumEdges { get; init; }
        }

        // EdgeAllocator provides temporary storage for new ClippedEdges that are
        // created during indexing.  It is essentially a stack model, where edges are
        // allocated as the recursion goes down and freed as it comes back up.
        //
        // It also provides a mutable vector of FaceEdges that is used when
        // incrementally updating the index (see AbsorbIndexCell).
        public class EdgeAllocator
        {
            public EdgeAllocator() { size_ = 0; }

            // Return a pointer to a newly allocated edge.  The EdgeAllocator
            // retains ownership.
            public void AddClippedEdge(ClippedEdge edge)
            {
                if (size_ >= clipped_edges_.Count)
                {
                    clipped_edges_.Add(edge);
                }
                size_++;
            }
            // Return the number of allocated edges.
            public int Count() { return size_; }

            // Reset the allocator to only contain the first "size" allocated edges.
            public void Reset(int size) { size_ = size; }

            // We can't use ClippedEdge[] because edges are not allowed to move
            // once they have been allocated.  Instead we keep a pool of allocated edges
            // that are all deleted together at the end.
            private int size_;
            private readonly List<ClippedEdge> clipped_edges_ = new List<ClippedEdge>();

            // On the other hand, we can use FaceEdge[] because they are allocated
            // only at one level during the recursion (namely, the level at which we
            // absorb an existing index cell).
            public List<FaceEdge> FaceEdges = new List<FaceEdge>();
        }
    }
}
