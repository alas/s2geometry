namespace S2Geometry;

public static partial class S2ShapeUtil
{
    public static class EdgePairs
    {
        // Ensure that we don't usually need to allocate memory when collecting the
        // edges in an S2ShapeIndex cell (which by default have about 10 edges).
        private class ShapeEdgeVector : Array16<ShapeEdge> { }

        /// <summary>
        /// Appends all edges in the given S2ShapeIndexCell to the given vector.
        /// </summary>
        private static void AppendShapeEdges(S2ShapeIndex index, S2ShapeIndexCell cell, ShapeEdgeVector shape_edges)
        {
            for (int s = 0; s < cell.NumClipped(); ++s)
            {
                var clipped = cell.Clipped(s);
                var shape = index.Shape(clipped.ShapeId);
                var num_edges = clipped.NumEdges;
                for (int i = 0; i < num_edges; i++)
                {
                    shape_edges.Add(new ShapeEdge(shape, clipped.Edge(i)));
                }
            }
        }

        /// <summary>
        /// Returns a vector containing all edges in the given S2ShapeIndexCell.
        /// (The result is returned as an output parameter so that the same storage can
        /// be reused, rather than allocating a new temporary vector each time.)
        /// </summary>
        private static void GetShapeEdges(S2ShapeIndex index, S2ShapeIndexCell cell, ShapeEdgeVector shape_edges)
        {
            shape_edges.Clear();
            AppendShapeEdges(index, cell, shape_edges);
        }

        /// <summary>
        /// Returns a vector containing all edges in the given S2ShapeIndexCell vector.
        /// (The result is returned as an output parameter so that the same storage can
        /// be reused, rather than allocating a new temporary vector each time.)
        /// </summary>
        private static void GetShapeEdges(S2ShapeIndex index, List<S2ShapeIndexCell> cells, ShapeEdgeVector shape_edges)
        {
            shape_edges.Clear();
            foreach (var cell in cells)
            {
                AppendShapeEdges(index, cell, shape_edges);
            }
        }

        // A function that is called with pairs of crossing edges.  The function may
        // return false in order to request that the algorithm should be terminated,
        // i.e. no further crossings are needed.
        //
        // "is_interior" indicates that the crossing is at a point interior to both
        // edges (i.e., not at a vertex).  (The calling function already has this
        // information and it is moderately expensive to recompute.)
        public delegate bool EdgePairVisitor(ShapeEdge a, ShapeEdge b, bool is_interior);

        /// <summary>
        /// Given a vector of edges within an S2ShapeIndexCell, visit all pairs of
        /// crossing edges (of the given CrossingType).
        /// </summary>
        private static bool VisitCrossings(ShapeEdgeVector shape_edges, CrossingType type, bool need_adjacent, EdgePairVisitor visitor)
        {
            var min_crossing_sign = type == CrossingType.INTERIOR ? 1 : 0;
            var num_edges = shape_edges.Count;
            for (int i = 0; i + 1 < num_edges; i++)
            {
                var a = shape_edges[i];
                var j = i + 1;

                // A common situation is that an edge AB is followed by an edge BC.  We
                // only need to visit such crossings if "need_adjacent" is true (even if
                // AB and BC belong to different edge chains).
                if (!need_adjacent && a.V1 == shape_edges[j].V0)
                {
                    j++;
                    if (j >= num_edges) break;
                }

                var crosser = new S2EdgeCrosser(a.V0, a.V1);
                for (; j < num_edges; j++)
                {
                    var b = shape_edges[j];
                    if (crosser.C == S2Point.Empty || crosser.C != b.V0)
                    {
                        crosser.RestartAt(b.V0);
                    }
                    var sign = crosser.CrossingSign(b.V1);
                    if (sign >= min_crossing_sign)
                    {
                        if (!visitor(a, b, sign == 1)) return false;
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
        /// early if the given EdgePairVisitor function returns false (in which case
        /// VisitCrossings returns false as well).  "type" indicates whether all
        /// crossings should be visited, or only interior crossings.
        /// 
        /// If "need_adjacent" is false, then edge pairs of the form (AB, BC) may
        /// optionally be ignored (even if the two edges belong to different edge
        /// chains).  This option exists for the benefit of FindSelfIntersection(),
        /// which does not need such edge pairs (see below).
        /// </summary>
        private static bool VisitCrossings(S2ShapeIndex index, CrossingType type, bool need_adjacent, EdgePairVisitor visitor)
        {
            // TODO(ericv): Use brute force if the total number of edges is small enough
            // (using a larger threshold if the S2ShapeIndex is not constructed yet).
            var count = index.GetEnumerableCount();
            for (var pos = 0; pos < count; pos++)
            {
                var shape_edges = new ShapeEdgeVector();
                var icell = index.GetIndexCell(pos);
                var indexCell = icell.Value.Item2;
                GetShapeEdges(index, indexCell, shape_edges);
                if (!VisitCrossings(shape_edges, type, need_adjacent, visitor))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
        /// early if the given EdgePairVisitor function returns false (in which case
        /// VisitCrossings returns false as well).  "type" indicates whether all
        /// crossings should be visited, or only interior crossings.
        /// 
        /// CAVEAT: Crossings may be visited more than once.
        /// </summary>
        public static bool VisitCrossingEdgePairs(S2ShapeIndex index, CrossingType type, EdgePairVisitor visitor)
        {
            var needAdjacent = type == CrossingType.ALL;
            return VisitCrossings(index, type, needAdjacent, visitor);
        }

        /// <summary>
        /// Like the above, but visits all pairs of crossing edges where one edge comes
        /// from each S2ShapeIndex.
        /// 
        /// CAVEAT: Crossings may be visited more than once.
        /// </summary>
        public static bool VisitCrossingEdgePairs(S2ShapeIndex a_index, S2ShapeIndex b_index, CrossingType type, EdgePairVisitor visitor)
        {
            // We look for S2CellId ranges where the indexes of A and B overlap, and
            // then test those edges for crossings.

            // TODO(ericv): Use brute force if the total number of edges is small enough
            // (using a larger threshold if the S2ShapeIndex is not constructed yet).
            var ai = new RangeEnumerator(a_index);
            var bi = new RangeEnumerator(b_index);
            var ab = new IndexCrosser(a_index, b_index, type, visitor, false);  // Tests A against B
            var ba = new IndexCrosser(b_index, a_index, type, visitor, true);   // Tests B against A
            while (!ai.Done() || !bi.Done())
            {
                if (ai.RangeMax < bi.RangeMin)
                {
                    // The A and B cells don't overlap, and A precedes B.
                    ai.SeekTo(bi);
                }
                else if (bi.RangeMax < ai.RangeMin)
                {
                    // The A and B cells don't overlap, and B precedes A.
                    bi.SeekTo(ai);
                }
                else
                {
                    // One cell contains the other.  Determine which cell is larger.
                    var ab_relation = ai.Id.LowestOnBit() - bi.Id.LowestOnBit();
                    if (ab_relation > 0)
                    {
                        // A's index cell is larger.
                        if (!ab.VisitCrossings(ai, bi)) return false;
                    }
                    else if (ab_relation < 0)
                    {
                        // B's index cell is larger.
                        if (!ba.VisitCrossings(bi, ai)) return false;
                    }
                    else
                    {
                        // The A and B cells are the same.
                        if (ai.Cell.NumEdges() > 0 && bi.Cell.NumEdges() > 0)
                        {
                            if (!ab.VisitCellCellCrossings(ai.Cell, bi.Cell)) return false;
                        }
                        ai.MoveNext();
                        bi.MoveNext();
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// IndexCrosser is a helper class for finding the edge crossings between a
        /// pair of S2ShapeIndexes.  It is instantiated twice, once for the index pair
        /// (A,B) and once for the index pair (B,A), in order to be able to test edge
        /// crossings in the most efficient order.
        /// </summary>
        private class IndexCrosser
        {
            /// <summary>
            /// If "swapped" is true, the loops A and B have been swapped.  This affects
            /// how arguments are passed to the given loop relation, since for example
            /// A.Contains(B) is not the same as B.Contains(A).
            /// </summary>
            public IndexCrosser(S2ShapeIndex a_index, S2ShapeIndex b_index, CrossingType type, EdgePairVisitor visitor, bool swapped)
            {
                a_index_ = a_index;
                b_index_ = b_index;
                visitor_ = visitor;
                min_crossing_sign_ = type == CrossingType.INTERIOR ? 1 : 0;
                swapped_ = swapped;
                b_query_ = new S2CrossingEdgeQuery(b_index_);
            }

            /// <summary>
            /// Given two iterators positioned such that ai.id().Contains(bi.id()),
            /// visits all crossings between edges of A and B that intersect a.id().
            /// Terminates early and returns false if visitor_ returns false.
            /// Advances both iterators past ai.id().
            /// </summary>
            public bool VisitCrossings(RangeEnumerator ai, RangeEnumerator bi)
            {
                System.Diagnostics.Debug.Assert(ai.Id.Contains(bi.Id));
                if (ai.Cell.NumEdges() == 0)
                {
                    // Skip over the cells of B using binary search.
                    bi.SeekBeyond(ai);
                }
                else
                {
                    // If ai.id() intersects many edges of B, then it is faster to use
                    // S2CrossingEdgeQuery to narrow down the candidates.  But if it
                    // intersects only a few edges, it is faster to check all the crossings
                    // directly.  We handle this by advancing "bi" and keeping track of how
                    // many edges we would need to test.
                    int b_edges = 0;
                    b_cells_.Clear();
                    do
                    {
                        int cell_edges = bi.Cell.NumEdges();
                        if (cell_edges > 0)
                        {
                            b_edges += cell_edges;
                            if (b_edges >= kEdgeQueryMinEdges)
                            {
                                // There are too many edges, so use an S2CrossingEdgeQuery.
                                if (!VisitSubcellCrossings(ai.Cell, ai.Id)) return false;
                                bi.SeekBeyond(ai);
                                return true;
                            }
                            b_cells_.Add(bi.Cell);
                        }
                        bi.MoveNext();
                    } while (bi.Id <= ai.RangeMax);
                    if (b_cells_.Any())
                    {
                        // Test all the edge crossings directly.
                        GetShapeEdges(a_index_, ai.Cell, a_shape_edges_);
                        GetShapeEdges(b_index_, b_cells_, b_shape_edges_);
                        if (!VisitEdgesEdgesCrossings(a_shape_edges_, b_shape_edges_))
                        {
                            return false;
                        }
                    }
                }
                ai.MoveNext();
                return true;
            }
            private const int kEdgeQueryMinEdges = 23;

            /// <summary>
            /// Given two index cells, visits all crossings between edges of those cells.
            /// Terminates early and returns false if visitor_ returns false.
            /// </summary>
            public bool VisitCellCellCrossings(S2ShapeIndexCell a_cell, S2ShapeIndexCell b_cell)
            {
                // Test all edges of "a_cell" against all edges of "b_cell".
                GetShapeEdges(a_index_, a_cell, a_shape_edges_);
                GetShapeEdges(b_index_, b_cell, b_shape_edges_);
                return VisitEdgesEdgesCrossings(a_shape_edges_, b_shape_edges_);
            }

            private bool VisitEdgePair(ShapeEdge a, ShapeEdge b, bool is_interior)
            {
                if (swapped_)
                {
                    return visitor_(b, a, is_interior);
                }
                else
                {
                    return visitor_(a, b, is_interior);
                }
            }

            /// <summary>
            /// Visits all crossings of the current edge with all edges of the given index
            /// cell of B.  Terminates early and returns false if visitor_ returns false.
            /// </summary>
            private bool VisitEdgeCellCrossings(ShapeEdge a, S2ShapeIndexCell b_cell)
            {
                // Test the current edge of A against all edges of "b_cell".

                // Note that we need to use a new S2EdgeCrosser (or call Init) whenever we
                // replace the contents of b_shape_edges_, since S2EdgeCrosser requires that
                // its S2Point arguments point to values that persist between Init() calls.
                GetShapeEdges(b_index_, b_cell, b_shape_edges_);
                var crosser = new S2EdgeCrosser(a.V0, a.V1);
                foreach (var b in b_shape_edges_)
                {
                    if (crosser.C == S2Point.Empty || crosser.C != b.V0)
                    {
                        crosser.RestartAt(b.V0);
                    }
                    int sign = crosser.CrossingSign(b.V1);
                    if (sign >= min_crossing_sign_)
                    {
                        if (!VisitEdgePair(a, b, sign == 1)) return false;
                    }
                }
                return true;
            }

            // Visits all crossings of any edge in "a_cell" with any index cell of B that
            // is a descendant of "b_id".  Terminates early and returns false if
            // visitor_ returns false.
            private bool VisitSubcellCrossings(S2ShapeIndexCell a_cell, S2CellId b_id)
            {
                // Test all edges of "a_cell" against the edges contained in B index cells
                // that are descendants of "b_id".
                GetShapeEdges(a_index_, a_cell, a_shape_edges_);
                var b_root = new S2PaddedCell(b_id, 0);
                foreach (var a in a_shape_edges_)
                {
                    // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
                    // of B that might contain crossing edges.
                    var result = b_query_.VisitCells(a.V0, a.V1, b_root, (S2ShapeIndexCell cell) =>
                        VisitEdgeCellCrossings(a, cell));
                    if (!result) return false;
                }
                return true;
            }


            // Visits all crossings of any edge in "a_edges" with any edge in "b_edges".
            private bool VisitEdgesEdgesCrossings(ShapeEdgeVector a_edges, ShapeEdgeVector b_edges)
            {
                // Test all edges of "a_edges" against all edges of "b_edges".
                foreach (var a in a_edges)
                {
                    var crosser = new S2EdgeCrosser(a.V0, a.V1);
                    foreach (var b in b_edges)
                    {
                        if (crosser.C == S2Point.Empty || crosser.C != b.V0)
                        {
                            crosser.RestartAt(b.V0);
                        }
                        int sign = crosser.CrossingSign(b.V1);
                        if (sign >= min_crossing_sign_)
                        {
                            if (!VisitEdgePair(a, b, sign == 1)) return false;
                        }
                    }
                }
                return true;
            }

            private readonly S2ShapeIndex a_index_;
            private readonly S2ShapeIndex b_index_;
            private readonly EdgePairVisitor visitor_;
            private readonly int min_crossing_sign_;
            private readonly bool swapped_;

            // Temporary data declared here to avoid repeated memory allocations.
            private readonly S2CrossingEdgeQuery b_query_;
            private readonly List<S2ShapeIndexCell> b_cells_ = new();
            private readonly ShapeEdgeVector a_shape_edges_ = new();
            private readonly ShapeEdgeVector b_shape_edges_ = new();
        }

        //////////////////////////////////////////////////////////////////////

        // Helper function that formats a loop error message.  If the loop belongs to
        // a multi-loop polygon, adds a prefix indicating which loop is affected.
        private static void InitLoopError(S2ErrorCode code, string str, S2Shape.ChainPosition ap, bool is_polygon, out S2Error error)
        {
            error = new(code, str);
            if (is_polygon)
            {
                error = new(code, $"Loop {ap.ChainId}: {error.Text}");
            }
        }

        // Given two loop edges that cross (including at a shared vertex), return true
        // if there is a crossing error and set "error" to a human-readable message.
        private static bool FindCrossingError(S2Shape shape, ShapeEdge a, ShapeEdge b, bool is_interior, out S2Error error)
        {
            var is_polygon = shape.NumChains() > 1;
            var ap = shape.GetChainPosition(a.Id.EdgeId);
            var bp = shape.GetChainPosition(b.Id.EdgeId);
            if (is_interior)
            {
                if (ap.ChainId != bp.ChainId)
                {
                    error = new(S2ErrorCode.POLYGON_LOOPS_CROSS,
                        $"Loop {ap.ChainId} edge {ap.Offset} crosses loop {bp.ChainId} edge bp.offset");
                }
                else
                {
                    InitLoopError(S2ErrorCode.LOOP_SELF_INTERSECTION,
                        $"Edge {ap.Offset} crosses edge {bp.Offset}", ap, is_polygon, out error);
                }
                return true;
            }
            // Loops are not allowed to have duplicate vertices, and separate loops
            // are not allowed to share edges or cross at vertices.  We only need to
            // check a given vertex once, so we also require that the two edges have
            // the same end vertex.
            if (a.V1 != b.V1)
            {
                error = S2Error.OK;
                return false;
            }
            if (ap.ChainId == bp.ChainId)
            {
                InitLoopError(S2ErrorCode.DUPLICATE_VERTICES,
                    $"Edge {ap.Offset} has duplicate vertex with edge {bp.Offset}",
                    ap, is_polygon, out error);
                return true;
            }
            int a_len = shape.GetChain(ap.ChainId).Length;
            int b_len = shape.GetChain(bp.ChainId).Length;
            int a_next = (ap.Offset + 1 == a_len) ? 0 : ap.Offset + 1;
            int b_next = (bp.Offset + 1 == b_len) ? 0 : bp.Offset + 1;
            S2Point a2 = shape.ChainEdge(ap.ChainId, a_next).V1;
            S2Point b2 = shape.ChainEdge(bp.ChainId, b_next).V1;
            if (a.V0 == b.V0 || a.V0 == b2)
            {
                // The second edge index is sometimes off by one, hence "near".
                error = new(S2ErrorCode.POLYGON_LOOPS_SHARE_EDGE,
                    $"Loop {ap.ChainId} edge {ap.Offset} has duplicate near loop {bp.ChainId} edge {bp.Offset}");
                return true;
            }
            // Since S2ShapeIndex loops are oriented such that the polygon interior is
            // always on the left, we need to handle the case where one wedge contains
            // the complement of the other wedge.  This is not specifically detected by
            // GetWedgeRelation, so there are two cases to check for.
            //
            // Note that we don't need to maintain any state regarding loop crossings
            // because duplicate edges are detected and rejected above.
            var overlaps = S2WedgeRelations.GetWedgeRelation(a.V0, a.V1, a2, b.V0, b2) == WedgeRelation.WEDGE_PROPERLY_OVERLAPS &&
                S2WedgeRelations.GetWedgeRelation(a.V0, a.V1, a2, b2, b.V0) == WedgeRelation.WEDGE_PROPERLY_OVERLAPS;
            if (overlaps)
            {
                error = new(S2ErrorCode.POLYGON_LOOPS_CROSS,
                    $"Loop {ap.ChainId} edge {ap.Offset} crosses loop {bp.ChainId} edge {bp.Offset}");
                return true;
            }
            error = S2Error.OK;
            return false;
        }

        /// <summary>
        /// Given an S2ShapeIndex containing a single polygonal shape (e.g., an
        /// S2Polygon or S2Loop), return true if any loop has a self-intersection
        /// (including duplicate vertices) or crosses any other loop (including vertex
        /// crossings and duplicate edges) and set "error" to a human-readable error
        /// message.  Otherwise return false and leave "error" unchanged.
        /// 
        /// This method is used to implement the FindValidationError methods of S2Loop
        /// and S2Polygon.
        /// 
        /// TODO(ericv): Add an option to support S2LaxPolygonShape rules (i.e.,
        /// duplicate vertices and edges are allowed, but loop crossings are not).
        /// </summary>
        public static bool FindSelfIntersection(S2ShapeIndex index, out S2Error error)
        {
            if (index.NumShapeIds() == 0)
            {
                error = S2Error.OK;
                return false;
            }

            System.Diagnostics.Debug.Assert(index.NumShapeIds() == 1);
            var shape = index.Shape(0);

            // Visit all crossing pairs except possibly for ones of the form (AB, BC),
            // since such pairs are very common and FindCrossingError() only needs pairs
            // of the form (AB, AC).
            S2Error tmpErr = S2Error.OK;
            var res = !VisitCrossings(index, CrossingType.ALL, false /*need_adjacent*/, (a, b, is_interior) =>
                !FindCrossingError(shape, a, b, is_interior, out tmpErr));
            error = tmpErr;
            return res;
        }
    }
}
