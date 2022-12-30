// S2CrossingEdgeQuery is used to find edges or shapes that are crossed by
// an edge.  Here is an example showing how to index a set of polylines,
// and then find the polylines that are crossed by a given edge AB:
//
// void Test(S2Polyline[] polylines,`
//           S2Point a0, S2Point &a1) {
//   MutableS2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(absl.make_unique<S2Polyline.Shape>(polyline));
//   }
//   S2CrossingEdgeQuery query(out index);
//   foreach (var& edge in query.GetCrossingEdges(a, b, CrossingType.ALL)) {
//     System.Diagnostics.Debug.Assert(S2EdgeCrossings.CrossingSign(a0, a1, edge.v0(), edge.v1()) >= 0);
//   }
// }
//
// Note that if you need to query many edges, it is more efficient to declare
// a single S2CrossingEdgeQuery object and reuse it so that temporary storage
// does not need to be reallocated each time.
//
// If you want to find *all* pairs of crossing edges, use
// S2ShapeUtil.VisitCrossingEdgePairs() instead.

namespace S2Geometry;

using ShapeEdge = S2ShapeUtil.ShapeEdge;

public class S2CrossingEdgeQuery
{
    // For small loops it is faster to use brute force.  The threshold below was
    // determined using the benchmarks in the unit test.
    private const int kMaxBruteForceEdges = 27;

    // REQUIRES: "index" is not modified after this method is called.
    public S2CrossingEdgeQuery(S2ShapeIndex index)
    {
        Index = index;
    }

    // Default constructor; requires Init() to be called.
    public S2CrossingEdgeQuery() { }

    public S2ShapeIndex Index { get; }
    // Returns all edges that intersect the given query edge (a0,a1) and that
    // have the given CrossingType (ALL or INTERIOR).  Edges are sorted and
    // unique.
    public List<ShapeEdge> GetCrossingEdges(S2Point a0, S2Point a1, CrossingType type)
    {
        var edges = new List<ShapeEdge>();
        GetCrossingEdges(a0, a1, type, edges);
        return edges;
    }

    // A specialized version of GetCrossingEdges() that only returns the edges
    // that belong to a particular S2Shape.
    public List<ShapeEdge> GetCrossingEdges(S2Point a0, S2Point a1, S2Shape shape, CrossingType type)
    {
        var edges = new List<ShapeEdge>();
        GetCrossingEdges(a0, a1, shape, type, edges);
        return edges;
    }

    // These versions can be more efficient when they are called many times,
    // since they do not require allocating a new vector on each call.
    public void GetCrossingEdges(S2Point a0, S2Point a1, CrossingType type, List<ShapeEdge> edges)
    {
        edges.Clear();
        GetCandidates(a0, a1, tmp_candidates_);
        int min_sign = (type == CrossingType.ALL) ? 0 : 1;
        var crosser = new S2CopyingEdgeCrosser(a0, a1);
        int shape_id = -1;
        S2Shape? shape = null;
        foreach (var (shapeId, edgeId) in tmp_candidates_)
        {
            if (shapeId != shape_id)
            {
                shape_id = shapeId;
                shape = Index.Shape(shape_id);
            }
            int edge_id = edgeId;
            S2Shape.Edge b = shape.GetEdge(edge_id);
            if (crosser.CrossingSign(b.V0, b.V1) >= min_sign)
            {
                edges.Add(new(shape_id, edge_id, b));
            }
        }
    }

    public void GetCrossingEdges(S2Point a0, S2Point a1, S2Shape shape, CrossingType type, List<ShapeEdge> edges)
    {
        edges.Clear();
        GetCandidates(a0, a1, shape, tmp_candidates_);
        int min_sign = (type == CrossingType.ALL) ? 0 : 1;
        var crosser = new S2CopyingEdgeCrosser(a0, a1);
        foreach (var (_, edgeId) in tmp_candidates_)
        {
            int edge_id = edgeId;
            S2Shape.Edge b = shape.GetEdge(edge_id);
            if (crosser.CrossingSign(b.V0, b.V1) >= min_sign)
            {
                edges.Add(new(shape.Id, edge_id, b));
            }
        }
    }


    /////////////////////////// Low-Level Methods ////////////////////////////
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use, since they require the client to do
    // all the actual crossing tests.

    // Returns a superset of the edges that intersect a query edge (a0, a1).
    // This method is useful for clients that want to test intersections in some
    // other way, e.g. using S2.EdgeOrVertexCrossing().
    public List<Edge> GetCandidates(S2Point a0, S2Point a1)
    {
        var edges = new List<Edge>();
        GetCandidates(a0, a1, edges);
        return edges;
    }

    // A specialized version of GetCandidates() that only returns the edges that
    // belong to a particular S2Shape.
    public List<Edge> GetCandidates(S2Point a0, S2Point a1, S2Shape shape)
    {
        var edges = new List<Edge>();
        GetCandidates(a0, a1, shape, edges);
        return edges;
    }

    // These versions can be more efficient when they are called many times,
    // since they do not require allocating a new vector on each call.
    public void GetCandidates(S2Point a0, S2Point a1, List<Edge> edges)
    {
        edges.Clear();
        var sorted = new SortedSet<Edge>();
        int num_edges = Index.GetCountEdgesUpTo(kMaxBruteForceEdges + 1);
        VisitRawCandidates(a0, a1, (Edge id) =>
        {
            sorted.Add(id);
            return true;
        });
        edges.AddRange(sorted);
    }

    public void GetCandidates(S2Point a0, S2Point a1, S2Shape shape, List<Edge> edges)
    {
        edges.Clear();
        var sorted = new SortedSet<Edge>();
        VisitRawCandidates(a0, a1, shape, (Edge id) =>
        {
            sorted.Add(id);
            return true;
        });
        edges.AddRange(sorted);
    }

    // A function that is called with each candidate intersecting edge.  The
    // function may return false in order to request that the algorithm should
    // be terminated, i.e. no further crossings are needed.
    public delegate bool ShapeEdgeIdVisitor(Edge id);

    // Visits a superset of the edges that intersect the query edge (a0, a1),
    // terminating early if the given ShapeEdgeIdVisitor returns false (in which
    // case this function returns false as well).
    //
    // CAVEAT: Edges may be visited more than once.
    public bool VisitRawCandidates(S2Point a0, S2Point a1, ShapeEdgeIdVisitor visitor)
    {
        int num_edges = Index.GetCountEdgesUpTo(kMaxBruteForceEdges + 1);
        if (num_edges <= kMaxBruteForceEdges)
        {
            int num_shape_ids = Index.NumShapeIds();
            for (int s = 0; s < num_shape_ids; ++s)
            {
                var shape = Index.Shape(s);
                if (shape is null) continue;
                int num_shape_edges = shape.NumEdges();
                for (int e = 0; e < num_shape_edges; ++e)
                {
                    if (!visitor(new(s, e))) return false;
                }
            }
            return true;
        }
        return VisitCells(a0, a1, (S2ShapeIndexCell cell) =>
        {
            for (int s = 0; s < cell.NumClipped(); ++s)
            {
                var clipped = cell.Clipped(s);
                for (int j = 0; j < clipped.NumEdges; ++j)
                {
                    if (!visitor(new(clipped.ShapeId, clipped.Edge(j))))
                    {
                        return false;
                    }
                }
            }
            return true;
        });
    }

    public bool VisitRawCandidates(S2Point a0, S2Point a1, S2Shape shape, ShapeEdgeIdVisitor visitor)
    {
        int num_edges = shape.NumEdges();
        if (num_edges <= kMaxBruteForceEdges)
        {
            for (int e = 0; e < num_edges; ++e)
            {
                if (!visitor(new(shape.Id, e))) return false;
            }
            return true;
        }
        return VisitCells(a0, a1, (S2ShapeIndexCell cell) =>
        {
            var clipped = cell.FindClipped(shape.Id);
            if (clipped is null) return true;
            for (int j = 0; j < clipped.NumEdges; ++j)
            {
                if (!visitor(new(shape.Id, clipped.Edge(j)))) return false;
            }
            return true;
        });
    }

    // A function that is called with each S2ShapeIndexCell that might contain
    // edges intersecting the given query edge.  The function may return false
    // in order to request that the algorithm should be terminated, i.e. no
    // further crossings are needed.
    public delegate bool CellVisitor(S2ShapeIndexCell cell);

    // Visits all S2ShapeIndexCells that might contain edges intersecting the
    // given query edge (a0, a1), terminating early if the given CellVisitor
    // returns false (in which case this function returns false as well).
    //
    // NOTE: Each candidate cell is visited exactly once.
    public bool VisitCells(S2Point a0, S2Point a1, CellVisitor visitor)
    {
        visitor_ = visitor;
        var segments = new S2EdgeClipping.FaceSegmentVector();
        S2EdgeClipping.GetFaceSegments(a0, a1, segments);
        foreach (var segment in segments)
        {
            a0_ = segment.a;
            a1_ = segment.b;

            // Optimization: rather than always starting the recursive subdivision at
            // the top level face cell, instead we start at the smallest S2CellId that
            // contains the edge (the "edge root cell").  This typically lets us skip
            // quite a few levels of recursion since most edges are short.
            var edge_bound = R2Rect.FromPointPair(a0_, a1_);
            var pcell = new S2PaddedCell(S2CellId.FromFace(segment.face), 0);
            var edge_root = pcell.ShrinkToFit(edge_bound);

            // Now we need to determine how the edge root cell is related to the cells
            // in the spatial index (cell_map_).  There are three cases:
            //
            //  1. edge_root is an index cell or is contained within an index cell.
            //     In this case we only need to look at the contents of that cell.
            //  2. edge_root is subdivided into one or more index cells.  In this case
            //     we recursively subdivide to find the cells intersected by a0a1.
            //  3. edge_root does not intersect any index cells.  In this case there
            //     is nothing to do.
            var (relation, pos) = Index.LocateCell(edge_root);
            if (relation == S2ShapeIndex.CellRelation.INDEXED)
            {
                var icell = Index.GetIndexCell(pos);
                // edge_root is an index cell or is contained by an index cell (case 1).
                System.Diagnostics.Debug.Assert(icell.Value.Item1.Contains(edge_root));
                if (!visitor(icell.Value.Item2)) return false;
            }
            else if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED)
            {
                // edge_root is subdivided into one or more index cells (case 2).  We
                // find the cells intersected by a0a1 using recursive subdivision.
                if (!edge_root.IsFace()) pcell = new S2PaddedCell(edge_root, 0);
                if (!VisitCells(pcell, edge_bound)) return false;
            }
        }
        return true;
    }

    // Visits all S2ShapeIndexCells within "root" that might contain edges
    // intersecting the given query edge (a0, a1), terminating early if the
    // given CellVisitor returns false (in which case this function returns
    // false as well).
    //
    // NOTE: Each candidate cell is visited exactly once.
    //
    // REQUIRES: root.padding() == 0
    //   [This low-level method does not support padding; the argument is supplied
    //    as an S2PaddedCell in order to avoid constructing it repeatedly when
    //    this method is called using different query edges with the same root.]
    public bool VisitCells(S2Point a0, S2Point a1, S2PaddedCell root, CellVisitor visitor)
    {
        System.Diagnostics.Debug.Assert(root.Padding == 0);
        visitor_ = visitor;
        // We use padding when clipping to ensure that the result is non-empty
        // whenever the edge (a0, a1) intersects the given root cell.
        if (S2EdgeClipping.ClipToPaddedFace(a0, a1, (int)root.Id.Face(),
            S2EdgeClipping.kFaceClipErrorUVCoord, out a0_, out a1_))
        {
            R2Rect edge_bound = R2Rect.FromPointPair(a0_, a1_);
            if (root.Bound.Intersects(edge_bound))
            {
                return VisitCells(root, edge_bound);
            }
        }
        return true;
    }

    // Given a query edge AB and a cell "root", returns all S2ShapeIndex cells
    // within "root" that might contain edges intersecting AB.
    //
    // REQUIRES: root.padding() == 0 (see above)
    public void GetCells(S2Point a0, S2Point a1, S2PaddedCell root, List<S2ShapeIndexCell> cells)
    {
        cells.Clear();
        VisitCells(a0, a1, root, (S2ShapeIndexCell cell) =>
        {
            cells.Add(cell);
            return true;
        });
    }

    // Computes the index cells intersected by the current edge that are
    // descendants of "pcell" and calls visitor_ for each one.
    //
    // WARNING: This function is recursive with a maximum depth of 30.  The frame
    // size is about 2K in versions of GCC prior to 4.7 due to poor overlapping
    // of storage for temporaries.  This is fixed in GCC 4.7, reducing the frame
    // size to about 350 bytes (i.e., worst-case total stack usage of about 10K).
    private bool VisitCells(S2PaddedCell pcell, R2Rect edge_bound)
    {
        // This code uses S2PaddedCell because it has the methods we need for
        // efficient splitting, however the actual padding is required to be zero.
        System.Diagnostics.Debug.Assert(pcell.Padding == 0);

        var (pos, found) = Index.SeekCell(pcell.Id.RangeMin());
        KeyData<S2CellId, S2ShapeIndexCell>? icell = null;
        if (found)
        {
            icell = Index.GetIndexCell(pos);
        }

        if (!found || !icell.HasValue || icell.Value.Item1 > pcell.Id.RangeMax())
        {
            // The index does not contain "pcell" or any of its descendants.
            return true;
        }
        if (icell.Value.Item1 == pcell.Id)
        {
            return visitor_(icell.Value.Item2);
        }

        // Otherwise, split the edge among the four children of "pcell".
        R2Point center = pcell.Middle.Lo();
        if (edge_bound[0].Hi < center[0])
        {
            // Edge is entirely contained in the two left children.
            return ClipVAxis(edge_bound, center[1], 0, pcell);
        }
        else if (edge_bound[0].Lo >= center[0])
        {
            // Edge is entirely contained in the two right children.
            return ClipVAxis(edge_bound, center[1], 1, pcell);
        }
        else
        {
            var child_bounds = new R2Rect[2];
            SplitUBound(edge_bound, center[0], child_bounds);
            if (edge_bound[1].Hi < center[1])
            {
                // Edge is entirely contained in the two lower children.
                return (VisitCells(new S2PaddedCell(pcell, 0, 0), child_bounds[0]) &&
                        VisitCells(new S2PaddedCell(pcell, 1, 0), child_bounds[1]));
            }
            else if (edge_bound[1].Lo >= center[1])
            {
                // Edge is entirely contained in the two upper children.
                return (VisitCells(new S2PaddedCell(pcell, 0, 1), child_bounds[0]) &&
                        VisitCells(new S2PaddedCell(pcell, 1, 1), child_bounds[1]));
            }
            else
            {
                // The edge bound spans all four children.  The edge itself intersects
                // at most three children (since no padding is being used).
                return (ClipVAxis(child_bounds[0], center[1], 0, pcell) &&
                        ClipVAxis(child_bounds[1], center[1], 1, pcell));
            }
        }
    }

    // Given either the left (i=0) or right (i=1) side of a padded cell "pcell",
    // determine whether the current edge intersects the lower child, upper child,
    // or both children, and call VisitCells() recursively on those children.
    // "center" is the v-coordinate at the center of "pcell".
    private bool ClipVAxis(R2Rect edge_bound, double center, int i, S2PaddedCell pcell)
    {
        if (edge_bound[1].Hi < center)
        {
            // Edge is entirely contained in the lower child.
            return VisitCells(new S2PaddedCell(pcell, i, 0), edge_bound);
        }
        else if (edge_bound[1].Lo >= center)
        {
            // Edge is entirely contained in the upper child.
            return VisitCells(new S2PaddedCell(pcell, i, 1), edge_bound);
        }
        else
        {
            // The edge intersects both children.
            var child_bounds = new R2Rect[2];
            SplitVBound(edge_bound, center, child_bounds);
            return (VisitCells(new S2PaddedCell(pcell, i, 0), child_bounds[0]) &&
                    VisitCells(new S2PaddedCell(pcell, i, 1), child_bounds[1]));
        }
    }

    // Split the current edge into two child edges at the given u-value "u" and
    // return the bound for each child.
    private void SplitUBound(R2Rect edge_bound, double u, R2Rect[] child_bounds)
    {
        // See comments in MutableS2ShapeIndex.ClipUBound.
        double v = edge_bound[1].Project(
            S2EdgeClipping.InterpolateDouble(u, a0_[0], a1_[0], a0_[1], a1_[1]));

        // "diag_" indicates which diagonal of the bounding box is spanned by a0a1:
        // it is 0 if a0a1 has positive slope, and 1 if a0a1 has negative slope.
        int diag = ((a0_[0] > a1_[0]) != (a0_[1] > a1_[1])) ? 1 : 0;
        SplitBound(edge_bound, 0, u, diag, v, child_bounds);
    }

    // Split the current edge into two child edges at the given v-value "v" and
    // return the bound for each child.
    private void SplitVBound(R2Rect edge_bound, double v, R2Rect[] child_bounds)
    {
        double u = edge_bound[0].Project(
            S2EdgeClipping.InterpolateDouble(v, a0_[1], a1_[1], a0_[0], a1_[0]));
        int diag = ((a0_[0] > a1_[0]) != (a0_[1] > a1_[1])) ? 1 : 0;
        SplitBound(edge_bound, diag, u, 0, v, child_bounds);
    }

    // Split the current edge into two child edges at the given point (u,v) and
    // return the bound for each child.  "u_end" and "v_end" indicate which bound
    // endpoints of child 1 will be updated.
    private static void SplitBound(R2Rect edge_bound, int u_end, double u, int v_end, double v, R2Rect[] child_bounds)
    {
        var arr1 = new double[2] { edge_bound[0][0], edge_bound[0][1] };
        var arr2 = new double[2] { edge_bound[1][0], edge_bound[1][1] };
        arr1[1 - u_end] = u;
        arr2[1 - v_end] = v;
        child_bounds[0] = new R2Rect(new R1Interval(arr1), new R1Interval(arr2));
        System.Diagnostics.Debug.Assert(!child_bounds[0].IsEmpty());
        System.Diagnostics.Debug.Assert(edge_bound.Contains(child_bounds[0]));

        arr1 = new double[2] { edge_bound[0][0], edge_bound[0][1] };
        arr2 = new double[2] { edge_bound[1][0], edge_bound[1][1] };
        arr1[u_end] = u;
        arr2[v_end] = v;
        child_bounds[1] = new R2Rect(new R1Interval(arr1), new R1Interval(arr2));
        System.Diagnostics.Debug.Assert(!child_bounds[1].IsEmpty());
        System.Diagnostics.Debug.Assert(edge_bound.Contains(child_bounds[1]));
    }

    //////////// Temporary storage used while processing a query ///////////

    private R2Point a0_, a1_;
    private CellVisitor visitor_;

    // Avoids repeated allocation when methods are called many times.
    private readonly List<Edge> tmp_candidates_ = new();
}

// A parameter that controls the reporting of edge intersections.
//
//  - CrossingType.INTERIOR reports intersections that occur at a point
//    interior to both edges (i.e., not at a vertex).
//
//  - CrossingType.ALL reports all intersections, even those where two edges
//    intersect only because they share a common vertex.
public enum CrossingType { INTERIOR, ALL }
