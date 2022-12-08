namespace S2Geometry;

public static partial class S2ShapeIndexX
{
    // Returns an S2ShapeIndexRegion that wraps the given S2ShapeIndex.  Note that
    // it is efficient to return S2ShapeIndexRegion objects by value.
    public static S2ShapeIndexRegion<T> MakeS2ShapeIndexRegion<T>(this T index)
        where T : S2ShapeIndex => new(index);
}

// This class wraps an S2ShapeIndex object with the additional methods needed
// to implement the S2Region API, in order to allow S2RegionCoverer to compute
// S2CellId coverings of arbitrary collections of geometry.
//
// It also contains a method VisitIntersectingShapes() that may be used to
// efficiently visit all shapes that intersect an arbitrary S2CellId (not
// limited to cells in the index).
//
// These methods could conceivably be made part of S2ShapeIndex itself, but
// there are several advantages to having a separate class:
//
//  - The class can be templated in order to avoid virtual calls and memory
//    allocation (for iterators) when the concrete S2ShapeIndex type is known.
//
//  - Implementing these methods efficiently requires an S2ShapeIndex iterator,
//    and this design allows a single iterator to be allocated and reused.
//
//  - S2Region.Clone() is not a good fit for the S2ShapeIndex API because
//    it can't be implemented for some subtypes (e.g., EncodedS2ShapeIndex).
//
// Example usage:
//
// S2CellUnion GetCovering(S2ShapeIndex& index) {
//   S2RegionCoverer coverer;
//   coverer.Options().set_max_cells(20);
//   S2CellUnion covering;
//   coverer.GetCovering(index.MakeS2ShapeIndexRegion(), &covering);
//   return covering;
// }
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
public sealed class S2ShapeIndexRegion<TIndex> : IS2Region<S2ShapeIndexRegion<TIndex>> where TIndex : S2ShapeIndex
{
    #region Fields, Constants

    // This class is not thread-safe!
    private readonly S2ContainsPointQuery<TIndex> contains_query_;

    #endregion

    #region Constructors

    // Rather than calling this constructor, which requires specifying the
    // S2ShapeIndex type explicitly, the preferred idiom is to call
    // MakeS2ShapeIndexRegion() instead.  For example:
    //
    //   coverer.GetCovering(index.MakeS2ShapeIndexRegion(), &covering);
    public S2ShapeIndexRegion(TIndex index)
    {
        contains_query_ = new S2ContainsPointQuery<TIndex>(index);
    }

    #endregion

    #region S2ShapeIndexRegion

    public TIndex Index()
    {
        return contains_query_.Index;
    }

    // Computes the smallest S2Cell that covers the S2Cell range (first, last) and
    // adds this cell to "cell_ids".
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private static void CoverRange(S2CellId first, S2CellId last, List<S2CellId> cell_ids)
    {
        if (first == last)
        {
            // The range consists of a single index cell.
            cell_ids.Add(first);
        }
        else
        {
            // Add the lowest common ancestor of the given range.
            int level = first.CommonAncestorLevel(last);
            System.Diagnostics.Debug.Assert(level >= 0);
            cell_ids.Add(first.Parent(level));
        }
    }
    /*// Returns true if the indexed shape "clipped" in the indexed cell "id"
    // contains the point "p".
    //
    // REQUIRES: id.contains(S2CellId(p))
    bool Contains(S2CellId id, const S2ClippedShape& clipped,
                const S2Point& p) const;*/
    // REQUIRES: iter_.id() contains "p".
    private bool Contains(S2CellId id, S2ClippedShape clipped, S2Point p)
    {
        return contains_query_.ShapeContains(id, clipped, p);
    }

    // Returns true if any edge of the indexed shape "clipped" intersects the
    // cell "target".  It may also return true if an edge is very close to
    // "target"; the maximum error is less than 10 * S2Constants.DoubleEpsilon radians (about
    // 15 nanometers).
    private bool AnyEdgeIntersects(S2ClippedShape clipped, S2Cell target)
    {
        var bound = target.BoundUV.Expanded(kMaxError);
        var face = target.Face;
        var shape = Index().Shape(clipped.ShapeId);
        int num_edges = clipped.NumEdges;
        for (int i = 0; i < num_edges; ++i)
        {
            var edge = shape.GetEdge(clipped.Edge(i));
            if (S2EdgeClipping.ClipToPaddedFace(edge.V0, edge.V1, face, kMaxError, out var p0, out var p1) &&
                S2EdgeClipping.IntersectsRect(p0, p1, bound))
            {
                return true;
            }
        }
        return false;
    }
    private static readonly double kMaxError = S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kIntersectsRectErrorUVDist;

    #endregion

    #region ICustomCloneable

    /// <returns>
    /// returns a *shallow* copy; it does not make a copy of the
    /// underlying S2ShapeIndex.
    /// </returns>
    public object CustomClone() => new S2ShapeIndexRegion<TIndex>(Index());

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        var covering = new List<S2CellId>();
        GetCellUnionBound(covering);
        return new S2CellUnion(covering).GetCapBound();
    }
    public S2LatLngRect GetRectBound()
    {
        var covering = new List<S2CellId>();
        GetCellUnionBound(covering);
        return new S2CellUnion(covering).GetRectBound();
    }

    // This method currently returns at most 4 cells, unless the index spans
    // multiple faces in which case it may return up to 6 cells.
    public void GetCellUnionBound(List<S2CellId> cell_ids)
    {
        // We find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  There are
        // two cases:
        //
        //  - If the index intersects two or more faces, then for each intersected
        //    face we add one cell to the covering.  Rather than adding the entire
        //    face, instead we add the smallest S2Cell that covers the S2ShapeIndex
        //    cells within that face.
        //
        //  - If the index intersects only one face, then we first find the smallest
        //    cell S that contains the index cells (just like the case above).
        //    However rather than using the cell S itself, instead we repeat this
        //    process for each of its child cells.  In other words, for each
        //    child cell C we add the smallest S2Cell C' that covers the index cells
        //    within C.  This extra step is relatively cheap and produces much
        //    tighter coverings when the S2ShapeIndex consists of a small region
        //    near the center of a large S2Cell.
        //
        // The following code uses only a single Iterator object because creating an
        // Iterator may be relatively expensive for some S2ShapeIndex types (e.g.,
        // it may involve memory allocation).
        cell_ids = new List<S2CellId>(6);

        // Find the last S2CellId in the index.
        var count = Index().GetEnumerableCount();
        if (count < 1) return;  // Empty index.

        var pos = count - 1;
        var last_index_id = Index().GetCellId(pos)!.Value;
        pos = 0;
        var cellId = Index().GetCellId(pos)!.Value;
        if (cellId != last_index_id)
        {
            // The index has at least two cells.  Choose an S2CellId level such that
            // the entire index can be spanned with at most 6 cells (if the index
            // spans multiple faces) or 4 cells (it the index spans a single face).
            int level = cellId.CommonAncestorLevel(last_index_id) + 1;

            // For each cell C at the chosen level, we compute the smallest S2Cell
            // that covers the S2ShapeIndex cells within C.
            var last_id = last_index_id.Parent(level);
            for (var id = cellId.Parent(level); id != last_id; id = id.Next())
            {
                // If the cell C does not contain any index cells, then skip it.
                if (id.RangeMax() < cellId) continue;

                // Find the range of index cells contained by C and then shrink C so
                // that it just covers those cells.
                var first = cellId;
                var (pos2, _) = Index().SeekCell(id.RangeMax().Next());
                pos2--;
                cellId = Index().GetCellId(pos2)!.Value;
                CoverRange(first, cellId, cell_ids);
                cellId = Index().GetCellId(pos2 + 1)!.Value;
            }
        }
        CoverRange(cellId, last_index_id, cell_ids);
    }

    // Returns true if "target" is contained by any single shape.  If the cell
    // is covered by a union of different shapes then it may return false.
    //
    // The implementation is conservative but not exact; if a shape just barely
    // contains the given cell then it may return false.  The maximum error is
    // less than 10 * S2Constants.DoubleEpsilon radians (or about 15 nanometers).
    public bool Contains(S2Cell target)
    {
        var (relation, pos) = Index().LocateCell(target.Id);

        // If the relation is DISJOINT, then "target" is not contained.  Similarly if
        // the relation is SUBDIVIDED then "target" is not contained, since index
        // cells are subdivided only if they (nearly) intersect too many edges.
        if (relation != S2ShapeIndex.CellRelation.INDEXED) return false;

        // Otherwise, the iterator points to an index cell containing "target".
        // If any shape contains the target cell, we return true.
        var icell = Index().GetIndexCell(pos);
        System.Diagnostics.Debug.Assert(icell.Value.Item1.Contains(target.Id));
        var cell = icell.Value.Item2;
        for (int s = 0; s < cell.NumClipped(); ++s)
        {
            var clipped = cell.Clipped(s);
            // The shape contains the target cell iff the shape contains the cell
            // center and none of its edges intersects the (padded) cell interior.
            if (icell.Value.Item1 == target.Id)
            {
                if (clipped.NumEdges == 0 && clipped.ContainsCenter) return true;
            }
            else
            {
                // It is faster to call AnyEdgeIntersects() before Contains().
                if (Index().Shape(clipped.ShapeId).Dimension() == 2 &&
                    !AnyEdgeIntersects(clipped, target) &&
                    Contains(icell.Value.Item1, clipped, target.Center()))
                {
                    return true;
                }
            }
        }
        return false;
    }

    // Returns true if any shape intersects "target".
    //
    // The implementation is conservative but not exact; if a shape is just
    // barely disjoint from the given cell then it may return true.  The maximum
    // error is less than 10 * S2Constants.DoubleEpsilon radians (or about 15 nanometers).
    public bool MayIntersect(S2Cell target)
    {
        var (relation, pos) = Index().LocateCell(target.Id);

        // If "target" does not overlap any index cell, there is no intersection.
        if (relation == S2ShapeIndex.CellRelation.DISJOINT) return false;

        // If "target" is subdivided into one or more index cells, then there is an
        // intersection to within the S2ShapeIndex error bound.
        if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) return true;

        // Otherwise, the iterator points to an index cell containing "target".
        //
        // If "target" is an index cell itself, there is an intersection because index
        // cells are created only if they have at least one edge or they are
        // entirely contained by the loop.
        var icell = Index().GetIndexCell(pos);
        System.Diagnostics.Debug.Assert(icell.Value.Item1.Contains(target.Id));
        if (icell.Value.Item1 == target.Id) return true;

        // Test whether any shape intersects the target cell or contains its center.
        var cell = icell.Value.Item2;
        for (int s = 0; s < cell.NumClipped(); ++s)
        {
            var clipped = cell.Clipped(s);
            if (AnyEdgeIntersects(clipped, target)) return true;
            if (Contains(icell.Value.Item1, clipped, target.Center()))
            {
                return true;
            }
        }
        return false;
    }

    // A function that is called with shapes that intersect a target S2Cell.
    // "contains_target" means that the shape fully contains the target S2Cell.
    // The function should return true to continue visiting intersecting shapes,
    // or false to terminate the algorithm early.
    //
    // Note that the API allows non-const access to the visited shapes.
    //
    // ENSURES: shape != nullptr
    public delegate bool ShapeVisitor(S2Shape shape, bool contains_target);

    // Visits all shapes that intersect "target", terminating early if the
    // "visitor" return false (in which case VisitIntersectingShapes returns
    // false as well).  Each shape is visited at most once.
    //
    // This method can also be used to visit all shapes that fully contain
    // "target" (VisitContainingShapes) by simply having the ShapeVisitor
    // function immediately return true when "contains_target" is false.
    public bool VisitIntersectingShapes(S2Cell target, ShapeVisitor visitor)
    {
        var (cellRelation, pos) = Index().LocateCell(target.Id);
        S2ShapeIndex.Enumerator iter_ = new(Index());
        iter_.SetPosition(pos);
        switch (cellRelation)
        {
            case S2ShapeIndex.CellRelation.DISJOINT:
                return true;

            case S2ShapeIndex.CellRelation.SUBDIVIDED:
                {
                    // A shape contains the target cell iff it appears in at least one cell,
                    // it contains the center of all cells, and it has no edges in any cell.
                    // It is easier to keep track of whether a shape does *not* contain the
                    // target cell because boolean values default to false.
                    Dictionary<int, bool> shape_not_contains = new();
                    for (var max = target.Id.RangeMax();
                        !iter_.Done() && iter_.Id <= max; iter_.MoveNext())
                    {
                        var cell = iter_.Cell;
                        for (int s = 0; s < cell.NumClipped(); ++s)
                        {
                            var clipped = cell.Clipped(s);
                            shape_not_contains[clipped.ShapeId] |=
                                clipped.NumEdges > 0 || !clipped.ContainsCenter;
                        }
                    }
                    // TODO(user,b/210097200): Use structured bindings when we require
                    // C++17 in opensource.
                    foreach (var (shape_id, not_contains) in shape_not_contains)
                    {
                        if (!visitor(Index().Shape(shape_id), !not_contains)) return false;
                    }
                    return true;
                }

            case S2ShapeIndex.CellRelation.INDEXED:
                {
                    var cell = iter_.Cell;
                    for (int s = 0; s < cell.NumClipped(); ++s)
                    {
                        // The shape contains the target cell iff the shape contains the cell
                        // center and none of its edges intersects the (padded) cell interior.
                        var clipped = cell.Clipped(s);
                        bool contains = false;
                        if (iter_.Id == target.Id)
                        {
                            contains = clipped.NumEdges == 0 && clipped.ContainsCenter;
                        }
                        else
                        {
                            if (!AnyEdgeIntersects(clipped, target))
                            {
                                if (!Contains(iter_.Id, clipped, target.Center()))
                                {
                                    continue;  // Disjoint.
                                }
                                contains = true;
                            }
                        }
                        if (!visitor(Index().Shape(clipped.ShapeId), contains)) return false;
                    }
                    return true;
                }
        }

        throw new Exception("(FATAL) Unhandled S2ShapeIndex::CellRelation");
    }

    // Returns true if the given point is contained by any two-dimensional shape
    // (i.e., polygon).  Boundaries are treated as being semi-open (i.e., the
    // same rules as S2Polygon).  Zero and one-dimensional shapes are ignored by
    // this method (if you need more flexibility, see S2BooleanOperation).
    public bool Contains(S2Point p)
    {
        var (pos, found) = Index().LocatePoint(p);
        if (found)
        {
            var icell = Index().GetIndexCell(pos);
            var cell = icell.Value.Item2;
            for (int s = 0; s < cell.NumClipped(); ++s)
            {
                if (Contains(icell.Value.Item1, cell.Clipped(s), p))
                {
                    return true;
                }
            }
        }
        return false;
    }

    #endregion

    #region IEquatable

    public override bool Equals(object? other)
    {
        throw new NotImplementedException();
    }

    public bool Equals(S2ShapeIndexRegion<TIndex>? other)
    {
        throw new NotImplementedException();
    }

    public override int GetHashCode()
    {
        throw new NotImplementedException();
    }

    #endregion
}
