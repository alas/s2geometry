// MutableS2ShapeIndex is a class for in-memory indexing of polygonal geometry.
// The objects in the index are known as "shapes", and may consist of points,
// polylines, and/or polygons, possibly overlapping.  The index makes it very
// fast to answer queries such as finding nearby shapes, measuring distances,
// testing for intersection and containment, etc.  It is one of several
// implementations of the S2ShapeIndex interface (see EncodedS2ShapeIndex).
//
// MutableS2ShapeIndex allows not only building an index, but also updating it
// incrementally by adding or removing shapes (hence its name).  It is designed
// to be compact; usually the index is smaller than the underlying geometry.
// It is capable of indexing up to hundreds of millions of edges.  The index is
// also fast to construct.  The index size and construction time are guaranteed
// to be linear in the number of input edges.
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
// - S2ShapeIndexRegion: can be used together with S2RegionCoverer to
//                       approximate geometry as a set of S2CellIds.
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
// MutableS2ShapeIndex has an Encode() method that allows the index to be
// serialized.  An encoded S2ShapeIndex can be decoded either into its
// original form (MutableS2ShapeIndex) or into an EncodedS2ShapeIndex.  The
// key property of EncodedS2ShapeIndex is that it can be constructed
// instantaneously, since the index is kept in its original encoded form.
// Data is decoded only when an operation needs it.  For example, to determine
// which shapes(s) contain a given query point only requires decoding the data
// in the S2ShapeIndexCell that contains that point.

namespace S2Geometry;

using System;
using System.Runtime.InteropServices;
using System.Threading;
using ShapeEdgeId = S2ShapeUtil.ShapeEdgeId;

public sealed class MutableS2ShapeIndex : S2ShapeIndex, IDisposable
{
    #region Fields and Properties

    // The default maximum number of edges per cell (not counting 'long' edges).
    // If a cell has more than this many edges, and it is not a leaf cell, then it
    // is subdivided.  This flag can be overridden via MutableS2ShapeIndex.Options.
    // Reasonable values range from 10 to about 50 or so.  Small values makes
    // queries faster, while large values make construction faster and use less memory.
    private const int s2shape_index_default_max_edges_per_cell = 10;

    // FLAGS_s2shape_index_tmp_memory_budget
    //
    // Attempt to limit the amount of temporary memory allocated while building or
    // updating a MutableS2ShapeIndex to at most this number of bytes.  This is 
    // achieved by splitting the updates into multiple batches when necessary.
    // (The memory required is proportional to the number of edges being updated
    // at once.)
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
    private const long s2shape_index_tmp_memory_budget = 100L << 20 /*100 MB*/;

    // FLAGS_s2shape_index_cell_size_to_long_edge_ratio
    //
    // The maximum cell size, relative to an edge's length, for which that edge is
    // considered 'long'.  Cell size is defined as the average edge length of all
    // cells at a given level.  For example, a value of 2.0 means that an edge E
    // is long at cell level k iff the average edge length at level k is at most
    // twice the length of E.  Long edges edges are not counted towards the
    // max_edges_per_cell() limit because such edges typically need to be
    // propagated to several children, which increases time and memory costs
    // without commensurate benefits.
    //
    // The maximum cell size, relative to an edge's length, for which that
    // edge is considered 'long'.  Long edges are not counted towards the
    // max_edges_per_cell() limit.  The size and speed of the index are
    // typically not very sensitive to this parameter.  Reasonable values range
    // from 0.1 to 10, with smaller values causing more aggressive subdivision
    // of long edges grouped closely together.
    private const double s2shape_index_cell_size_to_long_edge_ratio = 1.0;

    // FLAGS_s2shape_index_min_short_edge_fraction
    //
    // The minimum fraction of 'short' edges that must be present in a cell in
    // order for it to be subdivided.  If this parameter is non-zero then the
    // total index size and construction time are guaranteed to be linear in the
    // number of input edges; this prevents the worst-case quadratic space and
    // time usage that can otherwise occur with certain input configurations.
    // Specifically, the maximum index size is
    //
    //     O((c1 + c2 * (1 - f) / f) * n)
    //
    // where n is the number of input edges, f is this parameter value, and
    // constant c2 is roughly 20 times larger than constant c1.  (The exact values
    // of c1 and c2 depend on the cell_size_to_long_edge_ratio and
    // max_edges_per_cell parameters and certain properties of the input geometry
    // such as whether it consists of O(1) shapes, whether it includes polygons,
    // and whether the polygon interiors are disjoint.)
    //
    // Reasonable parameter values range from 0.1 up to perhaps 0.95.  The main
    // factors to consider when choosing this parameter are:
    //
    //  - For pathological geometry, larger values result in indexes that are
    //    smaller and faster to construct but have worse query performance (due to
    //    having more edges per cell).  However note that even a setting of 0.1
    //    reduces the worst case by 100x compared with a setting of 0.001.
    //
    //  - For normal geometry, values up to about 0.8 result in indexes that are
    //    virtually unchanged except for a slight increase in index construction
    //    time (proportional to the parameter value f) for very large inputs.
    //    With millions of edges, indexing time increases by about (15% * f),
    //    e.g. a parameter value of 0.5 slows down indexing for very large inputs
    //    by about 7.5%.  (Indexing time for small inputs is not affected.)
    //
    //  - Values larger than about 0.8 start to affect index construction even for
    //    normal geometry, resulting in smaller indexes and faster construction
    //    times but gradually worse query performance.
    //
    // Essentially this parameter provides control over a space-time tradeoff that
    // largely affects only pathological geometry.  The default value of 0.2 was
    // chosen to make index construction as fast as possible while still
    // protecting against possible quadratic space usage.
    //
    // The minimum fraction of 'short' edges that must be present in a cell in
    // order for it to be subdivided.If this parameter is non-zero then the
    // total index size and construction time are guaranteed to be linear in the
    // number of input edges, where the constant of proportionality has the
    // form (c1 + c2* (1 - f) / f).  Reasonable values range from 0.1 to
    // perhaps 0.95.  Values up to about 0.8 have almost no effect on 'normal'
    // geometry except for a small increase in index construction time
    // (proportional to f) for very large inputs.For worst-case geometry,
    // larger parameter values result in indexes that are smaller and faster
    // to construct but have worse query performance(due to having more edges
    // per cell).  Essentially this parameter provides control over a space-time
    // tradeoff that largely affects only pathological geometry.
    public const double s2shape_index_min_short_edge_fraction = 0.2;

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
    public static readonly double kCellPadding = 2 * (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kEdgeClipErrorUVCoord);

    // When adding a new encoding, be aware that old binaries will not be able
    // to decode it.
    public const byte kCurrentEncodingVersionNumber = 0;

    // The following memory estimates are based on heap profiling.

    // The batch sizes during a given update gradually decrease as the space
    // occupied by the index itself grows.  In order to do this, we need a
    // conserative lower bound on how much the index grows per edge.
    //
    // The final size of a MutableS2ShapeIndex depends mainly on how finely the
    // index is subdivided, as controlled by Options.max_edges_per_cell() and
    // --s2shape_index_default_max_edges_per_cell. For realistic values of
    // max_edges_per_cell() and shapes with moderate numbers of edges, it is
    // difficult to get much below 8 bytes per edge.  (The minimum possible size
    // is 4 bytes per edge (to store a 32-bit edge id in an S2ClippedShape) plus
    // 24 bytes per shape (for the S2ClippedShape itself plus a pointer in the
    // shapes_ vector.)  Note that this value is a lower bound; a typical final
    // index size is closer to 24 bytes per edge.
    private const int kFinalBytesPerEdge = 8;

    // The temporary memory consists mainly of the FaceEdge and ClippedEdge
    // structures plus a ClippedEdge pointer for every level of recursive
    // subdivision.  For very large indexes this can be 200 bytes per edge.
    // subdivision.  This can be more than 220 bytes per edge even for typical
    // geometry.  (The pathological worst case is higher, but we don't use this to
    // determine the batch sizes.)
    private const int kTmpBytesPerEdge = 226;

    // We arbitrarily limit the number of batches as a safety measure.  With the
    // current default memory budget of 100 MB, this limit is not reached even
    // when building an index of 350 million edges.
    private const int kMaxBatches = 100;

    // The options supplied for this index.
    public Options Options_ { get; set; }

    // The shapes in the index, accessed by their shape id.  Removed shapes are
    // replaced by null pointers.
    private readonly List<S2Shape?> shapes_ = [];

    // A map from S2CellId to the set of clipped shapes that intersect that
    // cell.  The cell ids cover a set of non-overlapping regions on the
    // sphere.  Note that this field is updated lazily (see below).  Const
    // methods *must* call MaybeApplyUpdates() before accessing this field.
    // (The easiest way to achieve this is simply to use an Iterator.)
    private readonly List<S2ShapeIndexIdCell> cell_map_ = []; // gtl.btree_map

    // The id of the first shape that has been queued for addition but not
    // processed yet.
    private int pending_additions_begin_ = 0;

    // Reads and writes to this field are guarded by "lock_".
    private IndexStatus index_status_ = IndexStatus.FRESH;

    // The set of shapes that have been queued for removal but not processed
    // yet.  Note that we need to copy the edge data since the caller is free to
    // destroy the shape once Release() has been called.  This field is present
    // only when there are removed shapes to process (to save memory).
    private List<RemovedShape>? pending_removals_;

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
    private int _channelPendingOperations = 0;
    private readonly Lock _channelLock = new();

    private readonly S2MemoryTracker.Client mem_tracker_ = new();

    #endregion

    #region Constructor

    // Create a MutableS2ShapeIndex with the given options.
    // Option values may be changed by calling Init().
    public MutableS2ShapeIndex(Options? options = null)
    {
        Options_ = options ?? new Options();
        Init();
    }

    #endregion

    public void Dispose() => Clear();

    // Initialize a MutableS2ShapeIndex with the given options.  This method may
    // only be called when the index is empty (i.e. newly created or Clear() has
    // just been called).  May be called before or after set_memory_tracker().
    public void Init()
    {
        MyDebug.Assert(shapes_.Count==0);
        // Memory tracking is not affected by this method.
    }

    // The number of distinct shape ids that have been assigned.  This equals
    // the number of shapes in the index provided that no shapes have ever been
    // removed.  (Shape ids are not reused.)
    public override int NumShapeIds() => shapes_.Count;

    // Returns a pointer to the shape with the given id, or null if the shape
    // has been removed from the index.
    public override S2Shape? Shape(int id) => shapes_[id];

    // Minimizes memory usage by requesting that any data structures that can be
    // rebuilt should be discarded.  This method invalidates all iterators.
    //
    // Like all non-methods, this method is not thread-safe.
    public override void Minimize()
    {
        mem_tracker_.Tally(-mem_tracker_.ClientUsageBytes);
        cell_map_.Clear();
        pending_additions_begin_ = 0;
        pending_removals_?.Clear();
        ResetChannel();
        MarkIndexStale();
        if (mem_tracker_.IsActive()) mem_tracker_.Tally(SpaceUsed());
    }

    // Appends an encoded representation of the S2ShapeIndex to "encoder".
    //
    // This method does not encode the S2Shapes in the index; it is the client's
    // responsibility to encode them separately.  For example:
    //
    //   S2ShapeUtil.CompactEncodeTaggedShapes(index, encoder);
    //   index.Encode(encoder);
    //
    // The encoded size is typically much smaller than the in-memory size.
    // Here are a few examples:
    //
    //  Number of edges     In-memory space used     Encoded size  (%)
    //  --------------------------------------------------------------
    //                8                      192                8   4%
    //              768                   18,264            2,021  11%
    //        3,784,212               80,978,992       17,039,020  21%
    //
    // The encoded form also has the advantage of being a contiguous block of
    // memory.
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder)
    {
        // The version number is encoded in 2 bits, under the assumption that by the
        // time we need 5 versions the first version can be permanently retired.
        // This only saves 1 byte, but that's significant for very small indexes.
        encoder.Ensure(Encoder.kVarintMax64);
        var max_edges = (ulong)Options_.MaxEdgesPerCell;
        encoder.PutVarUInt64(max_edges << 2 | kCurrentEncodingVersionNumber);

        // The index will be built anyway when we iterate through it, but building
        // it in advance lets us size the cell_ids vector correctly.
        ForceBuild();

        List<S2CellId> cell_ids = new(cell_map_.Count);
        StringVectorEncoder encoded_cells = new();

        for (var it = new Enumerator(this, S2ShapeIndex.InitialPosition.BEGIN);
            !it.Done(); it.MoveNext())
        {
            cell_ids.Add(it.Id);
            it.Cell!.Encode(NumShapeIds(), encoded_cells.AddViaEncoder());
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
        var num_shapes = (uint)shape_factory.Count;
        shapes_.Capacity = (int)num_shapes;
        for (var shape_id = 0; shape_id < num_shapes; ++shape_id)
        {
            var shape = shape_factory[shape_id];
            shape?.SetId(shape_id);
            shapes_.Add(shape);
        }

        var (success, cell_ids) = EncodedS2CellIdVector.Init(decoder);
        if (!success) return false;

        var (success2, encoded_cells) = EncodedStringVector.Init(decoder);
        if (!success2) return false;

        for (var i = 0; i < cell_ids!.Count; ++i)
        {
            var id = cell_ids[i];
            S2ShapeIndexCell cell = new();
            var decoder2 = encoded_cells!.GetDecoder(i);
            if (!cell.Decode((int)num_shapes, decoder2)) return false;
            // todo: maybe addsorted
            cell_map_.Add(new S2ShapeIndexIdCell(id, cell));
        }
        return true;
    }

    #region IEnumerator

    public override EnumeratorBase<S2ShapeIndexCell> GetNewEnumerator(InitialPosition pos)
    {
        MaybeApplyUpdates();

        return new Enumerator(this, pos);
    }

    public override (int pos, bool found) SeekCell(S2CellId target)
    {
        MaybeApplyUpdates();

        var tup = new S2ShapeIndexIdCell(target, new());
        var pos = cell_map_.BinarySearch(tup);
        if (pos < 0) return (~pos, false);

        while (pos > 0 && tup.CompareTo(cell_map_[pos]) == 0) pos--;

        return (pos + 1, true);
    }

    public override S2CellId? GetCellId(int pos)
    {
        // MaybeApplyUpdates has already been called

        if (pos >= cell_map_.Count || pos < 0) return null;

        return cell_map_[pos].Item1;
    }
    public override S2ShapeIndexCell? GetCell(int pos)
    {
        if (pos >= cell_map_.Count || pos < 0) return null;

        return cell_map_[pos].Item2;
    }

    #endregion

    // Specifies that memory usage should be tracked and/or limited by the given
    // S2MemoryTracker.  For example:
    //
    //   S2MemoryTracker tracker;
    //   tracker.set_limit(500 << 20);  // 500 MB memory limit
    //   MutableS2ShapeIndex index;
    //   index.set_memory_tracker(&tracker);
    //
    // If the memory limit is exceeded, an appropriate status is returned in
    // memory_tracker()->error() and any partially built index is discarded
    // (equivalent to calling Minimize()).
    //
    // This method may be called multiple times in order to switch from one
    // memory tracker to another or stop memory tracking altogether (by passing
    // nullptr) in which case the memory usage due to this index is subtracted.
    //
    // REQUIRES: The lifetime of "tracker" must exceed the lifetime of the index
    //           unless set_memory_tracker(nullptr) is called to stop memory
    //           tracking before the index destructor is called.
    //
    //           This implies that the S2MemoryTracker must be declared *before*
    //           the MutableS2ShapeIndex in the example above.
    //
    // CAVEATS:
    //
    //  - This method is not const and is therefore not thread-safe.
    //
    //  - Does not track memory used by the S2Shapes in the index.
    //
    //  - While the index representation itself is tracked very accurately,
    //    the temporary data needed for index construction is tracked using
    //    heuristics and may be underestimated or overestimated.
    //
    //  - Temporary memory usage is typically 10x larger than the final index
    //    size, however it can be reduced by specifying a suitable value for
    //    FLAGS_s2shape_index_tmp_memory_budget (the default is 100 MB).  If
    //    more temporary memory than this is needed during construction, index
    //    updates will be split into multiple batches in order to keep the
    //    estimated temporary memory usage below this limit.
    //
    //  - S2MemoryTracker.limit() has no effect on how much temporary memory
    //    MutableS2ShapeIndex will attempt to use during index construction; it
    //    simply causes an error to be returned when the limit would otherwise
    //    be exceeded.  If you set a memory limit smaller than 100MB and want to
    //    reduce memory usage rather than simply generating an error then you
    //    should also set FLAGS_s2shape_index_tmp_memory_budget appropriately.
    public S2MemoryTracker? MemoryTracker
    {
        get { return mem_tracker_.Tracker; }
        set
        {
            mem_tracker_.Tally(-mem_tracker_.ClientUsageBytes);
            mem_tracker_.Init(value);
            if (mem_tracker_.IsActive()) mem_tracker_.Tally(SpaceUsed());
        }
    }

    // Called to set the index status when the index needs to be rebuilt.
    private void MarkIndexStale()
    {
        // The UPDATING status can only be changed in ApplyUpdatesThreadSafe().
        if (index_status_ == IndexStatus.UPDATING) return;

        // If a memory tracking error has occurred we set the index status to FRESH
        // in order to prevent us from attempting to rebuild it.
        var status = (shapes_.Count==0 || !mem_tracker_.Ok()) ? IndexStatus.FRESH : IndexStatus.STALE;
        index_status_ = status;
    }

    // Takes ownership of the given shape and adds it to the index.  Also
    // assigns a unique id to the shape (shape.id()) and returns that id.
    // Shape ids are assigned sequentially starting from 0 in the order shapes
    // are added.  Invalidates all iterators and their associated data.
    //
    // Note that this method is not affected by S2MemoryTracker, i.e. shapes can
    // continue to be added even once the specified limit has been reached.
    public int Add(S2Shape shape)
    {
        // Additions are processed lazily by ApplyUpdates().  Note that in order to
        // avoid unexpected client behavior, this method continues to add shapes
        // even once the specified S2MemoryTracker limit has been exceeded.
        var id = NumShapeIds();
        shape.SetId(id);
        mem_tracker_.AddSpace(shapes_, 1);
        shapes_.Add(shape);
        IncChannel();
        MarkIndexStale();
        return id;
    }

    // Removes the given shape from the index and return ownership to the caller.
    // Invalidates all iterators and their associated data.
    public S2Shape Release(int shape_id)
    {
        // This class updates itself lazily, because it is much more efficient to
        // process additions and removals in batches.  However this means that when
        // a shape is removed we need to make a copy of all its edges, since the
        // client is free to delete "shape" once this call is finished.

        var shape = shapes_[shape_id];
        MyDebug.Assert(shape is not null);
        if (shape_id < pending_additions_begin_)
        {
            var num_edges = shape.NumEdges();
            var edges = new List<S2Shape.Edge>(num_edges);
            for (var e = 0; e < num_edges; ++e)
            {
                edges.Add(shape.GetEdge(e));
            }

            if (pending_removals_ is null)
            {
                if (!mem_tracker_.Tally(Marshal.SizeOf<List<RemovedShape>>()))
                {
                    Minimize();
                    return shape;
                }
                pending_removals_ = [];
            }
            var removed = new RemovedShape(shape.Id, shape.Dimension() == 2,
                shape.ContainsBruteForce(InteriorTracker.Origin()), edges);

            if (!mem_tracker_.AddSpace(removed.Edges, num_edges) ||
                !mem_tracker_.AddSpace(pending_removals_, 1))
            {
                Minimize();
                return shape;
            }
        }
        /*else
        {
            // We are removing a shape that has not yet been added to the index,
            // so there is nothing else to do.
        }*/
        MarkIndexStale();
        IncChannel();
        return shape;
    }

    // Resets the index to its original state and returns ownership of all
    // shapes to the caller.  This method is much more efficient than removing
    // all shapes one at a time.
    public void ReleaseAll()
    {
        shapes_.Clear();
        Minimize();
    }

    // Resets the index to its original state and deletes all shapes.  Any
    // options specified via Init() are preserved.
    public void Clear() => ReleaseAll();

    // Returns the number of bytes currently occupied by the index (including any
    // unused space at the end of vectors, etc). It has the same thread safety
    // as the other "const" methods (see introduction).
    public override int SpaceUsed()
    {
        var size = SizeHelper.SizeOf(this);
        size += NumShapeIds() * SizeHelper.SizeOf<S2Shape>();
        // cell_map_ itself is already included in sizeof(*this).
        size += cell_map_.Count - SizeHelper.SizeOf(cell_map_);
        size += cell_map_.Count * SizeHelper.SizeOf<S2ShapeIndexCell>();
        MutableS2ShapeIndex.Enumerator it_;
        for (it_ = new(this, InitialPosition.BEGIN); !it_.Done(); it_.MoveNext())
        {
            var cell = it_.Cell;
            size += cell!.NumClipped() * SizeHelper.SizeOf<S2ClippedShape>();
            for (var s = 0; s < cell.NumClipped(); ++s)
            {
                var clipped = cell.Clipped(s);
                size += clipped.NumEdges * sizeof(Int32);
            }
        }
        if (pending_removals_ is not null)
        {
            size += Marshal.SizeOf(pending_removals_);
            size += pending_removals_.Count * Marshal.SizeOf<RemovedShape>();

            foreach (var removed in pending_removals_)
            {
                size += removed.Edges.Capacity * Marshal.SizeOf<S2Shape.Edge>();
            }
        }

        return size;
    }

    // Calls to Add() and Release() are normally queued and processed on the
    // first subsequent query (in a thread-safe way).  Building the index lazily
    // in this way has several advantages, the most important of which is that
    // sometimes there *is* no subsequent query and the index doesn't need to be
    // built at all.
    //
    // In contrast, ForceBuild() causes any pending updates to be applied
    // immediately.  It is thread-safe and may be called simultaneously with
    // other "const" methods (see notes on thread safety above).  Similarly this
    // method is "const" since it does not modify the visible index contents.
    //
    // ForceBuild() should not normally be called since it prevents lazy index
    // construction (which is usually benficial).  Some reasons to use it
    // include:
    //
    //  - To exclude the cost of building the index from benchmark results.
    //  - To ensure that the first subsequent query is as fast as possible.
    //  - To ensure that the index can be built successfully without exceeding a
    //    specified S2MemoryTracker limit (see the constructor for details).
    //
    // Note that this method is thread-safe.
    public void ForceBuild() => MaybeApplyUpdates();

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
            ResetChannelInternal();
        }
    }
    public void ResetChannelInternal() => _channelPendingOperations = 0;

    // Given that the given shape is being updated, return true if it is being
    // removed (as opposed to being added).
    private bool IsShapeBeingRemoved(int shape_id) =>
        // All shape ids being removed are less than all shape ids being added.
        shape_id < pending_additions_begin_;

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
        if (index_status_ == IndexStatus.FRESH)
        {
        }
        else if (index_status_ == IndexStatus.UPDATING)
        {
            // Wait until the updating thread is finished.  We do this by attempting
            // to lock a mutex that is held by the updating thread.  When this mutex
            // is unlocked the index_status_ is guaranteed to be FRESH.
#if DEBUG
#pragma warning disable IDE0031 // Use null propagation
            if (ApplyUpdatesTask is not null)
#pragma warning restore IDE0031 // Use null propagation
#endif
                ApplyUpdatesTask.Wait();
        }
        else
        {
            lock (_channelLock)
            {
                index_status_ = IndexStatus.UPDATING;
                ApplyUpdatesTask = ApplyUpdatesInternalAsync();
                index_status_ = IndexStatus.FRESH;
            }
        }
    }
    private Task? ApplyUpdatesTask;

    // This method updates the index by applying all pending additions and
    // removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
    private async Task ApplyUpdatesInternalAsync()
    {
        await Task.Run(() =>
        {
            // Check whether we have so many edges to process that we should process
            // them in multiple batches to save memory.  Building the index can use up
            // to 20x as much memory (per edge) as the final index size.
            var batches = GetUpdateBatches();
            for (int i = 0; i < batches.Count; i++)
            {
                var batch = batches[i];
                if (mem_tracker_.IsActive())
                {
                    MyDebug.Assert(mem_tracker_.ClientUsageBytes == SpaceUsed(), "Invariant.");
                }
                Array6<List<FaceEdge>> all_edges = new(() => []);

                ReserveSpace(ref batch, all_edges);
                if (!mem_tracker_.Ok())
                {
                    batches[i] = batch;
                    Minimize();
                    return;
                }

                InteriorTracker tracker = new();
                if (pending_removals_ is not null)
                {
                    // The first batch implicitly includes all shapes being removed.
                    foreach (var pending_removal in pending_removals_)
                    {
                        RemoveShape(pending_removal, all_edges, tracker);
                    }
                    pending_removals_ = null;
                }

                // A batch consists of zero or more full shapes followed by zero or one
                // partial shapes.  The loop below handles all such cases.
                var begin = batch.Begin;
                var beginShapeId = begin.ShapeId;
                var beginEdgeId = begin.EdgeId;
                for (; new ShapeEdgeId(beginShapeId, beginEdgeId) < batch.End;
                    ++beginShapeId, beginEdgeId = 0)
                {
                    var shape = shapes_[beginShapeId];
                    if (shape is null) continue;  // Already removed.
                    int edges_end = beginShapeId == batch.End.ShapeId
                        ? batch.End.EdgeId
                        : shape.NumEdges();
                    AddShape(shape, beginEdgeId, edges_end, all_edges, tracker);
                }
                batches[i] = new(new(beginShapeId, beginEdgeId), batch.End, batch.NumEdges);

                for (var face = 0; face < 6; ++face)
                {
                    UpdateFaceEdges(face, all_edges[face], tracker);
                    // Save memory by clearing vectors after we are done with them.
                    all_edges[face].Clear();
                }

                pending_additions_begin_ = batch.End.ShapeId;
                if (beginEdgeId > 0 && batch.End.EdgeId == 0)
                {
                    // We have just finished adding the edges of shape that was split over
                    // multiple batches.  Now we need to mark the interior of the shape, if
                    // any, by setting contains_center() on the appropriate index cells.
                    FinishPartialShape(tracker.PartialShapeId);
                }
                if (mem_tracker_.IsActive())
                {
                    mem_tracker_.Tally(-mem_tracker_.ClientUsageBytes);
                    if (!mem_tracker_.Tally(SpaceUsed()))
                    {
                        Minimize();
                        return;
                    }
                }
            }
            ResetChannelInternal();
        });
    }

    // Count the number of edges being updated, and break them into several
    // batches if necessary to reduce the amount of memory needed.  (See the
    // documentation for FLAGS_s2shape_index_tmp_memory_budget_mb.)
    private List<BatchDescriptor> GetUpdateBatches()
    {
        // Count the edges being removed and added.
        var num_edges_removed = 0;
        if (pending_removals_ is not null)
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
            var shape = shapes_[id];
            if (shape != null)
            {
                num_edges_added += shape.NumEdges();
            }
        }

        BatchGenerator batch_gen = new(num_edges_removed, num_edges_added,
                                 pending_additions_begin_);
        for (int id = pending_additions_begin_; id < shapes_.Count; ++id)
        {
            var shape = shapes_[id];
            if (shape is not null) batch_gen.AddShape(id, shape.NumEdges());
        }
        return batch_gen.Finish();
    }

    // Reserve an appropriate amount of space for the top-level face edges in the
    // current batch.  This data structure uses about half of the temporary memory
    // needed during index construction.  Furthermore, if the arrays are grown via
    // push_back() then up to 10% of the total run time consists of copying data
    // as these arrays grow, so it is worthwhile to preallocate space for them.
    private void ReserveSpace(ref BatchDescriptor batch, Array6<List<FaceEdge>> all_edges)
    {
        // The following accounts for the temporary space needed for everything
        // except the FaceEdge vectors (which are allocated separately below).
        long other_usage = batch.NumEdges * (kTmpBytesPerEdge - Marshal.SizeOf<FaceEdge>());

        // If the number of edges is relatively small, then the fastest approach is
        // to simply reserve space on every face for the maximum possible number of
        // edges.  (We use a different threshold for this calculation than for
        // deciding when to break updates into batches because the cost/benefit
        // ratio is different.  Here the only extra expense is that we need to
        // sample the edges to estimate how many edges per face there are, and
        // therefore we generally use a lower threshold.)
        long kMaxCheapBytes =
            Math.Min(s2shape_index_tmp_memory_budget / 2,
                30L << 20 /*30 MB*/);
        long face_edge_usage = batch.NumEdges * 6 * Marshal.SizeOf<FaceEdge>();
        if (face_edge_usage <= kMaxCheapBytes)
        {
            if (!mem_tracker_.TallyTemp(face_edge_usage + other_usage))
            {
                return;
            }
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
        if (pending_removals_ is not null)
        {
            foreach (var removed in pending_removals_)
            {
                edge_id += removed.Edges.Count;
                while (edge_id >= sample_interval)
                {
                    edge_id -= sample_interval;
                    face_count[S2.GetFace(removed.Edges[edge_id].V0)] += 1;
                }
            }
        }
        var begin = batch.Begin;
        var beginShapeId = begin.ShapeId;
        var beginEdgeId = begin.EdgeId;
        for (; new ShapeEdgeId(beginShapeId, beginEdgeId) < batch.End;
             ++beginShapeId, beginEdgeId = 0)
        {
            var shape = shapes_[begin.ShapeId];
            if (shape is null) continue;  // Already removed.
            int edges_end = begin.ShapeId == batch.End.ShapeId ? batch.End.EdgeId
                : shape.NumEdges();
            edge_id += edges_end - begin.EdgeId;
            while (edge_id >= sample_interval)
            {
                edge_id -= sample_interval;
                // For speed, we only count the face containing one endpoint of the
                // edge.  In general the edge could span all 6 faces (with padding), but
                // it's not worth the expense to compute this more accurately.
                face_count[S2.GetFace(shape.GetEdge(edge_id + begin.EdgeId).V0)] += 1;
            }
        }
        batch = new(new(beginShapeId, beginEdgeId), batch.End, batch.NumEdges);
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
        // On average, we reserve 2% extra space for each face that has geometry
        // (which could be up to 12% extra space overall, but typically 2%).

        // kMaxSemiWidth is the maximum semi-width over all probabilities p of a
        // 4-sigma binomial confidence interval with a sample size of 10,000.
        const double kMaxSemiWidth = 0.02;

        // First estimate the total amount of memory we are about to allocate.
        double multiplier = 1.0;
        for (int face = 0; face < 6; ++face)
        {
            if (face_count[face] != 0) multiplier += kMaxSemiWidth;
        }
        face_edge_usage = Convert.ToInt32(multiplier * batch.NumEdges * Marshal.SizeOf<FaceEdge>());
        if (!mem_tracker_.TallyTemp(face_edge_usage + other_usage))
        {
            return;
        }
        var sample_ratio = 1.0 / actual_sample_size;
        for (var face = 0; face < 6; ++face)
        {
            if (face_count[face] == 0) continue;
            var fraction = sample_ratio * face_count[face] + kMaxSemiWidth;
            all_edges[face].Capacity = 1 + (int)(fraction * batch.NumEdges);
        }
    }

    // Clips the edges of the given shape to the six cube faces, add the clipped
    // edges to "all_edges", and start tracking its interior if necessary.
    private static void AddShape(S2Shape shape, int edges_begin, int edges_end, Array6<List<FaceEdge>> all_edges, InteriorTracker tracker)
    {
        var shapeId = shape.Id;
        var has_interior = false;
        if (shape.Dimension() == 2)
        {
            // To add a single shape with an interior over multiple batches, we first
            // add all the edges without tracking the interior.  After all edges have
            // been added, the interior is updated in a separate step by setting the
            // contains_center() flags appropriately.
            if (edges_begin > 0 || edges_end < shape.NumEdges())
            {
                tracker.PartialShapeId = shapeId;
            }
            else
            {
                has_interior = true;
                tracker.AddShape(
                    shapeId,
                    shape.ContainsBruteForce(tracker.Focus()));
            }
        }

        for (var e = edges_begin; e < edges_end; ++e)
        {
            var edge = shape.GetEdge(e);
            var max_level = GetEdgeMaxLevel(edge);
            AddFaceEdge(all_edges, shapeId, e, max_level, has_interior, edge);
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

    private void FinishPartialShape(int shape_id)
    {
        if (shape_id < 0) return;  // The partial shape did not have an interior.
        var shape = shapes_[shape_id];

        // Filling in the interior of a partial shape can grow the cell_map_
        // significantly, however the new cells have just one shape and no edges.
        // The following is a rough estimate of how much extra memory is needed
        // based on experiments.  It assumes that one new cell is required for every
        // 10 shape edges, and that the cell map uses 50% more space than necessary
        // for the new entries because they are inserted between existing entries
        // (which means that the btree nodes are not full).
        if (mem_tracker_.IsActive())
        {
            var new_usage = Convert.ToInt64(
                SpaceUsed() - mem_tracker_.ClientUsageBytes +
                0.1 * shape!.NumEdges() * (1.5 * Marshal.SizeOf<S2ShapeIndexCell>() + //cell_map_.value_type
                                            Marshal.SizeOf<S2ShapeIndexCell>() +
                                            Marshal.SizeOf<S2ClippedShape>()));
            if (!mem_tracker_.TallyTemp(new_usage)) return;
        }

        // All the edges of the partial shape have already been indexed, now we just
        // need to set the contains_center() flags appropriately.  We use a fresh
        // InteriorTracker for this purpose since we don't want to continue tracking
        // the interior state of any other shapes in this batch.
        //
        // We have implemented this below in the simplest way possible, namely by
        // scanning through the entire index.  In theory it would be more efficient
        // to keep track of the set of index cells that were modified when the
        // partial shape's edges were added, and then visit only those cells.
        // However in practice any shape that is added over multiple batches is
        // likely to occupy most or all of the index anyway, so it is faster and
        // simpler to just iterate through the entire index.
        //
        // "tmp_edges" below speeds up large polygon index construction by 3-12%.
        List<S2Shape.Edge> tmp_edges = [];  // Temporary storage.
        InteriorTracker tracker=new();
        tracker.AddShape(shape_id,
                         shape!.ContainsBruteForce(tracker.Focus()));
        S2CellId begin = S2CellId.Begin(S2.kMaxCellLevel);
        for (var index_it = 0; ; ++index_it)
        {
            if (tracker.ShapeIds().Count!=0)
            {
                // Check whether we need to add new cells that are entirely contained by
                // the partial shape.
                S2CellId fill_end = (index_it != cell_map_.Count)
                    ? cell_map_[index_it].Item1.RangeMin()
                    : S2CellId.End(S2.kMaxCellLevel);
                if (begin != fill_end)
                {
                    foreach (var cellid_u in S2CellUnion.FromBeginEnd(begin, fill_end))
                    {
                        S2ShapeIndexCell icell = new();
                        S2ClippedShape clipped = new(shape_id, 0, true);
                        icell.AddShape(clipped);
                        index_it = cell_map_.AddSorted(new(cellid_u, icell));
                        ++index_it;
                    }
                }
            }
            if (index_it == cell_map_.Count) break;

            // Now check whether the current index cell needs to be updated.
            var cellid = cell_map_[index_it].Item1;
            var cell = cell_map_[index_it].Item2;
            int n = cell.NumClipped();
            if (n > 0 && cell.Clipped(n - 1).ShapeId == shape_id)
            {
                // This cell contains edges of the partial shape.  If the partial shape
                // contains the center of this cell, we must update the index.
                S2PaddedCell pcell = new(cellid, kCellPadding);
                if (!tracker.AtCellId(cellid))
                {
                    tracker.MoveTo(pcell.GetEntryVertex());
                }
                tracker.DrawTo(pcell.GetCenter());
                var clipped = cell.Clipped(n - 1);
                int num_edges = clipped.NumEdges;
                MyDebug.Assert(num_edges > 0);
                for (int i = 0; i < num_edges; ++i)
                {
                    tmp_edges.Add(shape!.GetEdge(clipped.Edges[i]));
                }
                foreach (var edge in tmp_edges)
                {
                    tracker.TestEdge(shape_id, edge);
                }
                if (tracker.ShapeIds().Count!=0)
                {
                    // The partial shape contains the center of this index cell.
                    clipped.ContainsCenter = true;
                }
                tracker.DrawTo(pcell.GetExitVertex());
                foreach (var edge in tmp_edges)
                {
                    tracker.TestEdge(shape_id, edge);
                }
                tracker.NextCellId = cellid.Next();
                tmp_edges.Clear();

            }
            else if (tracker.ShapeIds().Count!=0)
            {
                // The partial shape contains the center of an existing index cell that
                // does not intersect any of its edges.
                S2ClippedShape clipped = new(shape_id, 0, true);
                cell.AddShape(clipped);
            }
            begin = cellid.RangeMax().Next();
        }
    }

    private static void AddFaceEdge(Array6<List<FaceEdge>> all_edges, Int32 shape_id, Int32 edge_id, Int32 max_level, bool has_interior, S2Shape.Edge edge)
    {
        R2Point edge_a, edge_b;
        // Fast path: both endpoints are on the same face, and are far enough from
        // the edge of the face that don't intersect any (padded) adjacent face.
        var a_face = S2.GetFace(edge.V0);
        if (a_face == S2.GetFace(edge.V1))
        {
            S2.ValidFaceXYZtoUV(a_face, edge.V0, out edge_a);
            S2.ValidFaceXYZtoUV(a_face, edge.V1, out edge_b);
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
        if (num_edges == 0 && tracker.ShapeIds().Count==0) return;

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
        // all its descendants) are not already present in the index.  It is set to
        // true during the recursion whenever we detect that the current cell is
        // disjoint from the index.  We could save a tiny bit of work by setting
        // this flag to true here on the very first update, however currently there
        // is no easy way to check that.  (It's not sufficient to test whether
        // cell_map_.empty() or pending_additions_begin_ == 0.)
        bool disjoint_from_index = false;
        if (num_edges > 0)
        {
            var shrunk_id = ShrinkToFit(pcell, bound);
            if (shrunk_id != pcell.Id)
            {
                // All the edges are contained by some descendant of the face cell.  We
                // can save a lot of work by starting directly with that cell, but if we
                // are in the interior of at least one shape then we need to create
                // index entries for the cells we are skipping over.
                SkipCellRange(face_id.RangeMin(), shrunk_id.RangeMin(),
                              tracker, alloc, disjoint_from_index);
                pcell = new S2PaddedCell(shrunk_id, kCellPadding);
                UpdateEdges(pcell, clipped_edges, tracker, alloc, disjoint_from_index);
                SkipCellRange(shrunk_id.RangeMax().Next(), face_id.RangeMax().Next(),
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
        if (shrunk_id != pcell.Id)
        {
            // Don't shrink any smaller than the existing index cells, since we need
            // to combine the new edges with those cells.  Use InitStale() to avoid
            // applying updates recursively.
            var (r, pos) = LocateCell(shrunk_id);
            if (r == S2CellRelation.INDEXED) { shrunk_id = GetCellId(pos)!.Value; }
        }
        return shrunk_id;
    }

    // Skip over the cells in the given range, creating index cells if we are
    // currently in the interior of at least one shape.
    private void SkipCellRange(S2CellId begin, S2CellId end, InteriorTracker tracker, EdgeAllocator alloc, bool disjoint_from_index)
    {
        // If we aren't in the interior of a shape, then skipping over cells is easy.
        if (tracker.ShapeIds().Count==0) return;

        // Otherwise generate the list of cell ids that we need to visit, and create
        // an index entry for each one.
        foreach (var skipped_id in S2CellUnion.FromBeginEnd(begin, end))
        {
            var d = new List<ClippedEdge>();
            UpdateEdges(new S2PaddedCell(skipped_id, kCellPadding), d, tracker, alloc, disjoint_from_index);
        }
    }

    // Given an edge and an interval "middle" along the v-axis, clip the edge
    // against the boundaries of "middle" and add the edge to the corresponding
    // children.
    //ABSL_ATTRIBUTE_ALWAYS_INLINE  // ~8% faster
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

    // Given a cell and a set of ClippedEdges whose bounding boxes intersect that
    // cell, add or remove all the edges from the index.  Temporary space for
    // edges that need to be subdivided is allocated from the given EdgeAllocator.
    // "disjoint_from_index" is an optimization hint indicating that cell_map_
    // does not contain any entries that overlap the given cell.
    private void UpdateEdges(S2PaddedCell pcell, List<ClippedEdge> edges, InteriorTracker tracker, EdgeAllocator alloc, bool disjoint_from_index)
    {
        // Cases where an index cell is not needed should be detected before this.
        MyDebug.Assert(edges.Count!=0 || tracker.ShapeIds().Count!=0);

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
            // the existing cell contents by "absorbing" the cell.  We use InitStale()
            // to avoid applying updates recursively.
            Enumerator iter = new(this);
            var r = iter.Locate(pcell.Id);
            if (r == S2CellRelation.DISJOINT)
            {
                disjoint_from_index = true;
            }
            else if (r == S2CellRelation.INDEXED)
            {
                // Absorb the index cell by transferring its contents to "edges" and
                // deleting it.  We also start tracking the interior of any new shapes.
                AbsorbIndexCell(pcell, iter, edges, tracker, alloc);
                index_cell_absorbed = true;
                disjoint_from_index = true;
            }
            else
            {
                MyDebug.Assert(S2CellRelation.SUBDIVIDED == r);
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
                   [new List<ClippedEdge>(num_edges), new List<ClippedEdge>(num_edges)],
                   [new List<ClippedEdge>(num_edges), new List<ClippedEdge>(num_edges)],
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
                if (child_edges[i][j].Count!=0 || tracker.ShapeIds().Count!=0)
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
    private void AbsorbIndexCell(S2PaddedCell pcell, Enumerator iter, List<ClippedEdge> edges, InteriorTracker tracker, EdgeAllocator alloc)
    {
        MyDebug.Assert(pcell.Id == iter.Id);

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
        // it is finished processing the contents of this cell.  (Note in the test
        // below that removed edges are always sorted before added edges.)
        if (tracker.IsActive && edges.Count!=0 &&
            IsShapeBeingRemoved(edges[0].FaceEdge.ShapeId))
        {
            // We probably need to update the InteriorTracker.  ("Probably" because
            // it's possible that all shapes being removed do not have interiors.)
            if (!tracker.AtCellId(pcell.Id))
            {
                tracker.MoveTo(pcell.GetEntryVertex());
            }
            tracker.DrawTo(pcell.GetExitVertex());
            tracker.NextCellId = pcell.Id.Next();
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
        // Save the state of the edges being removed so that it can be restored when
        // we are finished processing this cell and its children.  Below we not only
        // remove those edges but also add new edges whose state only needs to be
        // tracked within this subtree.  We don't need to save the state of the
        // edges being added because they aren't being removed from "edges" and will
        // therefore be updated normally as we visit this cell and its children.
        tracker.SaveAndClearStateBefore(pending_additions_begin_);

        // Create a FaceEdge for each edge in this cell that isn't being removed.
        var face_edges = alloc.FaceEdges;
        face_edges.Clear();
        bool tracker_moved = false;
        var cell = iter.Cell;
        for (var s = 0; s < cell.NumClipped(); ++s)
        {
            var clipped = cell.Clipped(s);
            int shape_id = clipped.ShapeId;
            var shape = shapes_[shape_id];
            if (shape is null) continue;  // This shape is being removed.
            int num_edges = clipped.NumEdges;

            // If this shape has an interior, start tracking whether we are inside the
            // shape.  UpdateEdges() wants to know whether the entry vertex of this
            // cell is inside the shape, but we only know whether the center of the
            // cell is inside the shape, so we need to test all the edges against the
            // line segment from the cell center to the entry vertex.
            var edge_shape_id = shape_id;
            var edge_has_interior = shape.Dimension() == 2 && shape_id != tracker.PartialShapeId;
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
                int e = clipped.Edges[i];
                var edge_edge_id = e;
                var edge_edge = shape.GetEdge(e);
                var edge_max_level = GetEdgeMaxLevel(edge_edge);
                if (edge_has_interior) tracker.TestEdge(shape_id, edge_edge);
                var face = pcell.Id.Face();
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

    // Returns the first level for which the given edge will be considered "long",
    // i.e. it will not count towards the max_edges_per_cell() limit.
    public static int GetEdgeMaxLevel(S2Shape.Edge edge)
    {
        // Compute the maximum cell edge length for which this edge is considered
        // "long".  The calculation does not need to be perfectly accurate, so we
        // use Norm() rather than Angle() for speed.
        double max_cell_edge = (edge.V0 - edge.V1).Norm() * s2shape_index_cell_size_to_long_edge_ratio;
        // Now return the first level encountered during subdivision where the
        // average cell edge length at that level is at most "max_cell_edge".
        return S2.kAvgEdge.GetLevelForMaxValue(max_cell_edge);
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
        if (edges.Count==0 && tracker.ShapeIds().Count==0)
        {
            // No index cell is needed.  (In most cases this situation is detected
            // before we get to this point, but this can happen when all shapes in a
            // cell are removed.)
            return true;
        }

        // We can show using amortized analysis that the total index size is
        //
        //     O(c1 * n + c2 * (1 - f) / f * n)
        //
        // where n is the number of input edges (and where we also count an "edge"
        // for each shape with an interior but no edges), f is the value of
        // FLAGS_s2shape_index_min_short_edge_fraction, and c1 and c2 are constants
        // where c2 is about 20 times larger than c1.
        //
        // First observe that the space used by a MutableS2ShapeIndex is
        // proportional to the space used by all of its index cells, and the space
        // used by an S2ShapeIndexCell is proportional to the number of edges that
        // intersect that cell plus the number of shapes that contain the entire
        // cell ("containing shapes").  Define an "index entry" as an intersecting
        // edge or containing shape stored by an index cell.  Our goal is then to
        // bound the number of index entries.
        //
        // We divide the index entries into two groups.  An index entry is "short"
        // if it represents an edge that was considered short in that index cell's
        // parent, and "long" otherwise.  (Note that the long index entries also
        // include the containing shapes mentioned above.)  We then bound the
        // maximum number of both types of index entries by associating them with
        // edges that were considered short in those index cells' parents.
        //
        // First consider the short index entries for a given edge E.  Let S be the
        // set of index cells that intersect E and where E was considered short in
        // those index cells' parents.  Since E was short in each parent cell, the
        // width of those parent cells is at least some fraction "g" of E's length
        // (as controlled by FLAGS_s2shape_index_cell_size_to_long_edge_ratio).
        // Therefore the minimum width of each cell in S is also at least some
        // fraction of E's length (i.e., g / 2).  This implies that there are at most
        // a constant number c1 of such cells, since they all intersect E and do not
        // overlap, which means that there are at most (c1 * n) short entries in
        // total.
        //
        // With index_cell_size_to_long_edge_ratio = 1.0 (the default value), it can
        // be shown that c1 = 10.  In other words, it is not possible for a given
        // edge to intersect more than 10 index cells where it was considered short
        // in those cells' parents.  The value of c1 can be reduced as low c1 = 4 by
        // increasing index_cell_size_to_long_edge_ratio to about 3.1.  (The reason
        // the minimum value is 3.1 rather than 2.0 is that this ratio is defined in
        // terms of the average edge length of cells at a given level, rather than
        // their minimum width, and 2 * (S2.kAvgEdge / S2.kMinWidth) ~= 3.1.)
        //
        // Next we consider the long index entries.  Let c2 be the maximum number of
        // index cells where a given edge E was considered short in those cells'
        // parents.  (Unlike the case above, we do not require that these cells
        // intersect E.)  Because the minimum width of each parent cell is at least
        // some fraction of E's length and the parent cells at a given level do not
        // overlap, there can be at most a small constant number of index cells at
        // each level where E is considered short in those cells' parents.  For
        // example, consider a very short edge E that intersects the midpoint of a
        // cell edge at level 0.  There are 16 cells at level 30 where E was
        // considered short in the parent cell, 12 cells at each of levels 29..2, and
        // 4 cells at levels 1 and 0 (pretending that all 6 face cells share a common
        // "parent").  This yields a total of c2 = 360 index cells.  This is actually
        // the worst case for index_cell_size_to_long_edge_ratio >= 3.1; with the
        // default value of 1.0 it is possible to have a few more index cells at
        // levels 29 and 30, for a maximum of c2 = 366 index cells.
        //
        // The code below subdivides a given cell only if
        //
        //     s > f * (s + l)
        //
        // where "f" is the min_short_edge_fraction parameter, "s" is the number of
        // short edges that intersect the cell, and "l" is the number of long edges
        // that intersect the cell plus an upper bound on the number of shapes that
        // contain the entire cell.  (It is an upper bound rather than an exact count
        // because we use the number of shapes that contain an arbitrary vertex of
        // the cell.)  Note that the number of long index entries in each child of
        // this cell is at most "l" because no child intersects more edges than its
        // parent or is entirely contained by more shapes than its parent.
        //
        // The inequality above can be rearranged to give
        //
        //    l < s * (1 - f) / f
        //
        // This says that each long index entry in a child cell can be associated
        // with at most (1 - f) / f edges that were considered short when the parent
        // cell was subdivided.  Furthermore we know that there are at most c2 index
        // cells where a given edge was considered short in the parent cell.  Since
        // there are only n edges in total, this means that the maximum number of
        // long index entries is at most
        //
        //    c2 * (1 - f) / f * n
        //
        // and putting this together with the result for short index entries gives
        // the desired bound.
        //
        // There are a variety of ways to make this bound tighter, e.g. when "n" is
        // relatively small.  For example when the indexed geometry satisfies the
        // requirements of S2BooleanOperation (i.e., shape interiors are disjoint)
        // and the min_short_edge_fraction parameter is not too large, then the
        // constant c2 above is only about half as big (i.e., c2 ~= 180).  This is
        // because the worst case under these circumstances requires having many
        // shapes whose interiors overlap.

        // Continue subdividing if the proposed index cell would contain too many
        // edges that are "short" relative to its size (as controlled by the
        // FLAGS_s2shape_index_cell_size_to_long_edge_ratio parameter).  Usually "too
        // many" means more than options_.max_edges_per_cell(), but this value might
        // be increased if the cell has a lot of long edges and/or containing shapes.
        // This strategy ensures that the total index size is linear (see above).
        if (edges.Count > Options_.MaxEdgesPerCell)
        {
            int max_short_edges =
                Math.Max(Options_.MaxEdgesPerCell,
                    (int)s2shape_index_min_short_edge_fraction *
                        (edges.Count + tracker.ShapeIds().Count));
            int count = 0;
            foreach (var edge in edges)
            {
                count += (pcell.Level < edge.FaceEdge.MaxLevel) ? 1 : 0;
                if (count > max_short_edges) return false;
            }
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
        if (tracker.IsActive && edges.Count!=0)
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
                S2ClippedShape clipped = new(cshape_id, 0, true);
                cnextHasNext = cnext.MoveNext();
                cell.AddShape(clipped);
            }
            else
            {
                // Count the number of edges for this shape and allocate space for them.
                while (enext < edges.Count &&
                       edges[enext].FaceEdge.ShapeId == eshape_id)
                {
                    ++enext;
                }
                S2ClippedShape clipped = new(eshape_id, enext - ebegin);
                for (int e = ebegin; e < enext; ++e)
                {
                    clipped.Edges[e - ebegin] = edges[e].FaceEdge.EdgeId;
                }
                if (cshape_id == eshape_id)
                {
                    clipped.ContainsCenter = true;
                    cnextHasNext = cnext.MoveNext();
                }
                cell.AddShape(clipped);
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
        if (tracker.IsActive && edges.Count!=0)
        {
            tracker.DrawTo(pcell.GetExitVertex());
            TestAllEdges(edges, tracker);
            tracker.NextCellId = pcell.Id.Next();
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
        MyDebug.Assert(!res.Bound.IsEmpty());
        MyDebug.Assert(edge.Bound.Contains(res.Bound));
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

    private enum IndexStatus
    {
        STALE,     // There are pending updates.
        UPDATING,  // Updates are currently being applied.
        FRESH,     // There are no pending updates.
    }

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

    /// <summary>
    /// The representation of an edge that has been queued for removal.
    /// </summary>
    private record RemovedShape(
        int ShapeId,
        bool HasInterior, //Belongs to a shape of dimension 2.
        bool ContainsTrackerOrigin,
        List<S2Shape.Edge> Edges);

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
    public record FaceEdge(
        Int32 ShapeId,         // The shape that this edge belongs to
        Int32 EdgeId,          // Edge id within that shape
        Int32 MaxLevel,        // Not desirable to subdivide this edge beyond this level
        bool HasInterior,      // Belongs to a shape of dimension 2.
        S2Shape.Edge Edge,     // The edge endpoints
        R2Point A, R2Point B); // The edge endpoints, clipped to a given face

    public record ClippedEdge(
        FaceEdge FaceEdge,  // The original unclipped edge
        R2Rect Bound);      // Bounding box for the clipped portion

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
        // S2CellId space-filling curve by drawing an edge from S2.Origin to this
        // point and counting how many shape edges cross this edge.
        public InteriorTracker()
        {
            b_ = Origin();
            next_cellid_ = S2CellId.Begin(S2.kMaxCellLevel);
        }

        // Returns the initial focus point when the InteriorTracker is constructed
        // (corresponding to the start of the S2CellId space-filling curve).
        public static S2Point Origin()
        {
            // The start of the S2CellId space-filling curve.
            return S2.FaceUVtoXYZ(0, -1, -1).Normalize();
        }

        // Returns the current focus point (see above).
        public S2Point Focus() { return b_; }

        // Returns true if any shapes are being tracked.
        public bool IsActive { get; private set; } = false;

        // Adds a shape whose interior should be tracked.  "contains_focus" indicates
        // whether the current focus point is inside the shape.  Alternatively, if the
        // focus point is in the process of being moved (via MoveTo/DrawTo), you can
        // also specify "contains_focus" at the old focus point and call TestEdge()
        // for every edge of the shape that might cross the current DrawTo() line.
        // This updates the state to correspond to the new focus point.
        //
        // REQUIRES: shape.dimension() == 2
        public void AddShape(Int32 shape_id, bool contains_focus)
        {
            IsActive = true;
            if (contains_focus)
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
            set => next_cellid_ = value.RangeMin();
        }
        private S2CellId next_cellid_;

        // Returns true if the focus is already at the entry vertex of the given
        // S2CellId (provided that the caller calls set_next_cellid() as each cell
        // is processed).
        public bool AtCellId(S2CellId cellid) => cellid.RangeMin() == next_cellid_;

        // Makes an internal copy of the state for shape ids below the given limit,
        // and then clear the state for those shapes.  This is used during
        // incremental updates to track the state of added and removed shapes
        // separately.
        public void SaveAndClearStateBefore(Int32 limit_shape_id)
        {
            MyDebug.Assert(saved_ids_.Count==0);
            var limit = LowerBound(limit_shape_id);
            saved_ids_.AddRange(shape_ids_.Take(limit));
            shape_ids_.RemoveRange(0, limit);
            saved_is_active_ = IsActive;
        }

        // Restores the state previously saved by SaveAndClearStateBefore().  This
        // only affects the state for shape_ids below "limit_shape_id".
        public void RestoreStateBefore(Int32 limit_shape_id)
        {
            shape_ids_.RemoveRange(0, LowerBound(limit_shape_id));
            shape_ids_.InsertRange(0, saved_ids_);
            saved_ids_.Clear();
            IsActive = saved_is_active_;
        }

        // Removes "shape_id" from shape_ids_ if it exists, otherwise insert it.
        private void ToggleShape(int shape_id)
        {
            // Since shape_ids_.size() is typically *very* small (0, 1, or 2), it turns
            // out to be significantly faster to maintain a sorted array rather than
            // using an STL set or btree_set.
            if (shape_ids_.Count==0)
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
        private readonly S2EdgeCrosser crosser_ = new();
        private readonly List<int> shape_ids_ = [];

        // Shape ids saved by SaveAndClearStateBefore().  The state is never saved
        // recursively so we don't need to worry about maintaining a stack.
        private readonly List<int> saved_ids_ = [];

        // As an optimization, we also save is_active_ so that RestoreStateBefore()
        // can deactivate the tracker again in the case where the shapes being added
        // and removed do not have an interior, but some existing shapes do.
        private bool saved_is_active_;

        // If non-negative, indicates that only some edges of the given shape are
        // being added and therefore its interior should not be tracked yet.
        public int PartialShapeId { get; set; } = -1;
    }

    // A BatchDescriptor represents a set of pending updates that will be applied
    // at the same time.  The batch consists of all edges in (shape id, edge id)
    // order from "begin" (inclusive) to "end" (exclusive).  Note that the last
    // shape in a batch may have only some of its edges added.  The first batch
    // also implicitly includes all shapes being removed.  "num_edges" is the
    // total number of edges that will be added or removed in this batch.
    //
    // REQUIRES: If end.edge_id != 0, it must refer to a valid edge.
    public readonly record struct BatchDescriptor(ShapeEdgeId Begin, ShapeEdgeId End, int NumEdges)
    {
        public override string ToString() => $"({Begin}, {End}, {NumEdges})";
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
        private readonly List<ClippedEdge> clipped_edges_ = [];

        // On the other hand, we can use FaceEdge[] because they are allocated
        // only at one level during the recursion (namely, the level at which we
        // absorb an existing index cell).
        public List<FaceEdge> FaceEdges = [];
    }

    // The purpose of BatchGenerator is to divide large updates into batches such
    // that all batches use approximately the same amount of high-water memory.
    // This class is defined here so that it can be tested independently.
    public record class BatchGenerator
    {
        // Given the total number of edges that will be removed and added, prepares
        // to divide the edges into batches.  "shape_id_begin" identifies the first
        // shape whose edges will be added.
        public BatchGenerator(int num_edges_removed, int num_edges_added,
            int shape_id_begin)
        {
            max_batch_sizes_ = GetMaxBatchSizes(num_edges_removed, num_edges_added);
            batch_begin_ = new(shape_id_begin, 0);
            shape_id_end_ = shape_id_begin;

            if (max_batch_sizes_.Count > 1)
            {
                /*S2_VLOG(1) << "Removing " << num_edges_removed << ", adding "
                        << num_edges_added << " edges in " << max_batch_sizes_.size()
                        << " batches";*/

                // Duplicate the last entry to simplify next_max_batch_size().
                max_batch_sizes_.Add(max_batch_sizes_.Last());

                // We process edge removals before additions, and edges are always removed
                // in a single batch.  The reasons for this include: (1) removed edges use
                // quite a bit of memory (about 50 bytes each) and this space can be freed
                // immediately when we process them in one batch; (2) removed shapes are
                // expected to be small fraction of the index size in typical use cases
                // (e.g. incremental updates of large indexes), and (3) AbsorbIndexCell()
                // uses (shape(id) == nullptr) to detect when a shape is being removed, so
                // in order to split the removed shapes into multiple batches we would need
                // a different approach (e.g., temporarily adding fake entries to shapes_
                // and restoring them back to nullptr as shapes are removed).  Removing
                // individual shapes over multiple batches would be even more work.
                batch_size_ = num_edges_removed;
            }
        }

        // Indicates that the given shape will be added to the index.  Shapes with
        // few edges will be grouped together into a single batch, while shapes with
        // many edges will be split over several batches if necessary.
        public void AddShape(int shape_id, int num_edges)
        {
            int batch_remaining = MaxBatchSize() - batch_size_;
            if (num_edges <= batch_remaining)
            {
                ExtendBatch(num_edges);
            }
            else if (num_edges <= NextMaxBatchSize())
            {
                // Avoid splitting shapes across batches unnecessarily.
                FinishBatch(0, new ShapeEdgeId(shape_id, 0));
                ExtendBatch(num_edges);
            }
            else
            {
                // This shape must be split across at least two batches.  We simply fill
                // each batch until the remaining edges will fit in two batches, and then
                // divide those edges such that both batches have the same amount of
                // remaining space relative to their maximum size.
                int e_begin = 0;
                while (batch_remaining + NextMaxBatchSize() < num_edges)
                {
                    e_begin += batch_remaining;
                    FinishBatch(batch_remaining, new ShapeEdgeId(shape_id, e_begin));
                    num_edges -= batch_remaining;
                    batch_remaining = MaxBatchSize();
                }
                // Figure out how many edges to add to the current batch so that it will
                // have the same amount of remaining space as the next batch.
                int n = (num_edges + batch_remaining - NextMaxBatchSize()) / 2;
                FinishBatch(n, new ShapeEdgeId(shape_id, e_begin + n));
                FinishBatch(num_edges - n, new ShapeEdgeId(shape_id + 1, 0));
            }
            shape_id_end_ = shape_id + 1;
        }

        // Returns a vector describing each batch.  This method should be called
        // once all shapes have been added.
        public List<BatchDescriptor> Finish()
        {
            // We must generate at least one batch even when num_edges_removed ==
            // num_edges_added == 0, because some shapes have an interior but no edges.
            // (Specifically, the full polygon has this property.)
            if (batches_.Count==0 || shape_id_end_ != batch_begin_.ShapeId)
            {
                FinishBatch(0, new ShapeEdgeId(shape_id_end_, 0));
            }
            return batches_;
        }

        // Returns a vector indicating the maximum number of edges in each batch.
        // (The actual batch sizes are adjusted later in order to avoid splitting
        // shapes between batches unnecessarily.)
        //
        // Divides "num_edges" edges into batches where each batch needs about the
        // same total amount of memory.  (The total memory needed by a batch consists
        // of the temporary memory needed to process the edges in that batch plus the
        // final representations of the edges that have already been indexed.)  It
        // uses the fewest number of batches (up to kMaxBatches) such that the total
        // memory usage does not exceed the combined final size of all the edges plus
        // FLAGS_s2shape_index_tmp_memory_budget.  Returns a vector of sizes
        // indicating the desired number of edges in each batch.
        private static List<int> GetMaxBatchSizes(int num_edges_removed, int num_edges_added)
        {
            // Check whether we can update all the edges at once.
            int num_edges_total = num_edges_removed + num_edges_added;
            const double tmp_memory_budget_bytes =
                s2shape_index_tmp_memory_budget;
            if (num_edges_total * kTmpBytesPerEdge <= tmp_memory_budget_bytes)
            {
                return [num_edges_total];
            }

            // Each batch is allowed to use up to "total_budget_bytes".  The memory
            // usage consists of some number of edges already added by previous batches
            // (at kFinalBytesPerEdge each), plus some number being updated in the
            // current batch (at kTmpBytesPerEdge each).  The available free space is
            // multiplied by (1 - kFinalBytesPerEdge / kTmpBytesPerEdge) after each
            // batch is processed as edges are converted into their final form.
            double final_bytes = num_edges_added * kFinalBytesPerEdge;
            const double kFinalBytesRatio = 1.0 * kFinalBytesPerEdge /
                                                kTmpBytesPerEdge;
            const double kTmpSpaceMultiplier = 1 - kFinalBytesRatio;

            // The total memory budget is the greater of the final size plus the allowed
            // temporary memory, or the minimum amount of memory required to limit the
            // number of batches to "kMaxBatches".
            var total_budget_bytes = Math.Max(
                final_bytes + tmp_memory_budget_bytes,
                final_bytes / (1 - MathUtils.IPow(kTmpSpaceMultiplier, kMaxBatches - 1)));

            // "ideal_batch_size" is the number of edges in the current batch before
            // rounding to an integer.
            double ideal_batch_size = total_budget_bytes / kTmpBytesPerEdge;

            // Removed edges are always processed in the first batch, even if this might
            // use more memory than requested (see the BatchGenerator constructor).
            var batch_sizes = new List<int>();
            int num_edges_left = num_edges_added;
            if (num_edges_removed > ideal_batch_size)
            {
                batch_sizes.Add(num_edges_removed);
            }
            else
            {
                num_edges_left += num_edges_removed;
            }
            for (int i = 0; num_edges_left > 0; ++i)
            {
                int batch_size = (int)(ideal_batch_size + 1);
                batch_sizes.Add(batch_size);
                num_edges_left -= batch_size;
                ideal_batch_size *= kTmpSpaceMultiplier;
            }
            MyDebug.Assert(batch_sizes.Count <= kMaxBatches);
            return [.. batch_sizes];
        }

        // Returns the maximum number of edges in the current batch.
        private int MaxBatchSize() { return max_batch_sizes_[batch_index_]; }

        // Returns the maximum number of edges in the next batch.
        private int NextMaxBatchSize() { return max_batch_sizes_[batch_index_ + 1]; }

        // Adds the given number of edges to the current batch.
        private void ExtendBatch(int num_edges)
        {
            batch_size_ += num_edges;
        }

        // Adds the given number of edges to the current batch, ending with the edge
        // just before "batch_end", and then starts a new batch.
        private void FinishBatch(int num_edges, ShapeEdgeId batch_end)
        {
            ExtendBatch(num_edges);
            batches_.Add(new BatchDescriptor(batch_begin_, batch_end, batch_size_));
            batch_begin_ = batch_end;
            batch_index_edges_left_ -= batch_size_;
            while (batch_index_edges_left_ < 0)
            {
                batch_index_edges_left_ += MaxBatchSize();
                batch_index_ += 1;
                batch_size_ = 0;
            }
        }

        // A vector representing the ideal number of edges in each batch; the batch
        // sizes gradually decrease to ensure that each batch uses approximately the
        // same total amount of memory as the index grows.  The actual batch sizes
        // are then adjusted based on how many edges each shape has in order to
        // avoid splitting shapes between batches unnecessarily.
        private readonly List<int> max_batch_sizes_;

        // The maximum size of the current batch is determined by how many edges
        // have been added to the index so far.  For example if GetBatchSizes()
        // returned {100, 70, 50, 30} and we have added 0 edges, the current batch
        // size is 100.  But if we have already added 90 edges then the current
        // batch size would be 70, and if have added 150 edges the batch size would
        // be 50.  We keep track of (1) the current index into batch_sizes and (2)
        // the number of edges remaining before we increment the batch index.
        private int batch_index_ = 0;
        private int batch_index_edges_left_ = 0;

        private ShapeEdgeId batch_begin_;  // The start of the current batch.
        private int shape_id_end_;         // One beyond the last shape to be added.
        private int batch_size_ = 0;       // The number of edges in the current batch.
        private readonly List<BatchDescriptor> batches_ = [];  // The completed batches so far.
    }

    public new sealed class Enumerator : EnumeratorBase<S2ShapeIndexCell>
    {
        private MutableS2ShapeIndex Index { get; }
        private int iter_;
        private readonly int end_;

        // Constructs an iterator positioned as specified.  By default iterators
        // are unpositioned, since this avoids an extra seek in this situation
        // where one of the seek methods (such as Locate) is immediately called.
        //
        // If you want to position the iterator at the beginning, e.g. in order to
        // loop through the entire index, do this instead:
        //
        //   for (MutableS2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN);
        //        !it.done(); it.Next()) { ... }
        //
        // Initializes an iterator for the given MutableS2ShapeIndex.  This method
        // may also be called in order to restore an iterator to a valid state
        // after the underlying index has been updated (although it is usually
        // easier just to declare a new iterator whenever required, since iterator
        // construction is cheap).
        //
        // Initialize an iterator for the given MutableS2ShapeIndex without
        // applying any pending updates.  This can be used to observe the actual
        // current state of the index without modifying it in any way.
        public Enumerator(MutableS2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED)
        {
            index.MaybeApplyUpdates();
            Index = index;
            end_ = Index.cell_map_.Count;
            if (pos == InitialPosition.BEGIN)
            {
                iter_ = 0;
            }
            else
            {
                iter_ = end_;
            }
            Refresh();
        }

        // Inherited non-virtual methods:
        //   S2CellId id() const;
        //   bool done() const;
        //   S2Point center() const;
        public override S2ShapeIndexCell? Cell =>
            // Since MutableS2ShapeIndex always sets the "cell_" field, we can skip the
            // logic in the base class that conditionally calls GetCell().
            CellInternal;

        #region S2CellIterator API

        /// <summary>
        /// Begin
        /// </summary>
        public override void Reset()
        {
            // Make sure that the index has not been modified since Init() was called.
            MyDebug.Assert(Index.IsFresh());
            iter_ = 0;
            Refresh();
        }

        public override void Finish()
        {
            iter_ = end_;
            Refresh();
        }

        public override bool MoveNext()
        {
            MyDebug.Assert(!Done());
            ++iter_;
            Refresh();
            return !Done();
        }
        public override bool MovePrevious()
        {
            if (iter_ == 0) return false;
            --iter_;
            Refresh();
            return true;
        }
        public override void Seek(S2CellId target)
        {
            iter_ = Index.cell_map_.GetLowerBound(new(target, default));
            Refresh();
        }

        public override bool Locate(S2Point target) => LocateImpl(this, target);

        public override S2CellRelation Locate(S2CellId target) => LocateImpl(this, target);

        public override void SetPosition(int position)
        {
            iter_ = position;
            Refresh();
        }

        public override void Dispose() { }

        #endregion

        protected override S2ShapeIndexCell GetCell() => throw new NotImplementedException();

        // Updates the IteratorBase fields.
        private void Refresh()
        {
            if (iter_ == end_)
            {
                SetFinished();
            }
            else
            {
                SetState(Index.GetCellId(iter_)!.Value, Index.GetCell(iter_));
            }
        }
    }
}
