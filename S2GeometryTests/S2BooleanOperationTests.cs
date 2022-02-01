namespace S2Geometry
{
    using static S2BooleanOperation;
    using static S2Builder.GraphOptions;
    using Loops = List<List<S2Point>>;

    public class S2BooleanOperationTests
    {
        private readonly ITestOutputHelper _logger;
        private readonly DegeneracyCoverageTest dct;

        public S2BooleanOperationTests(ITestOutputHelper logger)
        {
            _logger = logger;
            dct = new DegeneracyCoverageTest(logger);
        }

        // The polygon used in the polyline/polygon vertex tests below.
        private const string kVertexTestPolygonStr = "0:0, 0:1, 0:2, 0:3, 0:4, 0:5, 5:5, 5:4, 5:3, 5:2, 5:1, 5:0";

        // TODO(ericv): Clean up or remove these notes.
        //
        // Options to test:
        //   polygon_model:                   OPEN, SEMI_OPEN, CLOSED
        //   polyline_model:                  OPEN, SEMI_OPEN, CLOSED
        //   polyline_loops_have_boundaries:  true, false
        //   conservative:                    true, false
        //
        // Geometry combinations to test:
        //
        // Point/point:
        //  - disjoint, coincident
        // Point/polyline:
        //  - Start vertex, end vertex, interior vertex, degenerate polyline
        //  - With polyline_loops_have_boundaries: start/end vertex, degenerate polyline
        // Point/polygon:
        //  - Polygon interior, exterior, vertex
        //  - Vertex of degenerate sibling pair shell, hole
        //  - Vertex of degenerate single point shell, hole
        // Polyline/polyline:
        //  - Vertex intersection:
        //    - Start, end, interior, degenerate, loop start/end, degenerate loop
        //    - Test cases where vertex is not emitted because an incident edge is.
        //  - Edge/edge: interior crossing, duplicate, reversed, degenerate
        //  - Test that degenerate edges are ignored unless polyline has a single edge.
        //    (For example, AA has one edge but AAA has no edges.)
        // Polyline/polygon:
        //  - Vertex intersection: polyline vertex cases already covered, but test
        //    polygon normal vertex, sibling pair shell/hole, single vertex shell/hole
        //    - Also test cases where vertex is not emitted because an edge is.
        //  - Edge/edge: interior crossing, duplicate, reversed
        //  - Edge/interior: polyline edge in polygon interior, exterior
        // Polygon/polygon:
        //  - Vertex intersection:
        //    - normal vertex, sibling pair shell/hole, single vertex shell/hole
        //    - Also test cases where vertex is not emitted because an edge is.
        //    - Test that polygons take priority when there is a polygon vertex and
        //      also isolated polyline vertices.  (There should not be any points.)
        //  - Edge/edge: interior crossing, duplicate, reversed
        //  - Interior/interior: polygons in interior/exterior of other polygons

        [Fact]
        public void Test_S2BooleanOperation_DegeneratePolylines()
        {
            // Verify that degenerate polylines are preserved under all boundary models.
            Options options = new();
            var a = "# 0:0, 0:0 #";
            var b = "# #";
            options.PolylineModel_ = PolylineModel.OPEN;
            ExpectResult(OpType.UNION, options, a, b, a);
            options.PolylineModel_ = (PolylineModel.SEMI_OPEN);
            ExpectResult(OpType.UNION, options, a, b, a);
            options.PolylineModel_ = (PolylineModel.CLOSED);
            ExpectResult(OpType.UNION, options, a, b, a);
        }

        [Fact]
        public void Test_S2BooleanOperation_DegeneratePolygons()
        {
            // Verify that degenerate polygon features (single-vertex and sibling pair
            // shells and holes) are preserved under all boundary models.
            Options options = new();
            var a = "# # 0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 3:3; 6:6; 7:7, 8:8";
            var b = "# #";
            options.PolygonModel_ = (PolygonModel.OPEN);
            ExpectResult(OpType.UNION, options, a, b, a);
            options.PolygonModel_ = (PolygonModel.SEMI_OPEN);
            ExpectResult(OpType.UNION, options, a, b, a);
            options.PolygonModel_ = (PolygonModel.CLOSED);
            ExpectResult(OpType.UNION, options, a, b, a);
        }

        [Fact]
        public void Test_S2BooleanOperation_PointPoint()
        {
            Options options = new();
            var a = "0:0 | 1:0 # #";
            var b = "0:0 | 2:0 # #";
            // Note that these results have duplicates, which is correct.  Clients can
            // eliminated the duplicates with the appropriate GraphOptions.
            ExpectResult(OpType.UNION, options, a, b,
                         "0:0 | 0:0 | 1:0 | 2:0 # #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "0:0 | 0:0 # #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "1:0 # #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "1:0 | 2:0 # #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PointOpenPolyline()
        {
            // Tests operations between an open polyline and its vertices.
            //
            // The polyline "3:0, 3:0" consists of a single degenerate edge and contains
            // no points (since polyline_model() is OPEN).  Since S2BooleanOperation
            // preserves degeneracies, this means that the union includes *both* the
            // point 3:0 and the degenerate polyline 3:0, since they do not intersect.
            //
            // This test uses Options.polyline_loops_have_boundaries() == true, which
            // means that the loop "4:0, 5:0, 4:0" does not contain the vertex "4:0".
            Options options = new();
            options.PolylineModel_ = (PolylineModel.OPEN);
            var a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
            var b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "0:0 | 2:0 | 3:0 | 4:0 " +
                         "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "1:0 | 5:0 # #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "0:0 | 2:0 | 3:0 | 4:0 # #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "0:0 | 2:0 | 3:0 | 4:0" +
                         "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PointOpenPolylineLoopBoundariesFalse()
        {
            // With Options.polyline_loops_have_boundaries() == false, the loop
            // "4:0, 5:0, 4:0" has two vertices, both of which are contained.
            Options options = new();
            options.PolylineModel_ = (PolylineModel.OPEN);
            options.PolylineLoopsHaveBoundaries = (false);
            var a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
            var b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "0:0 | 2:0 | 3:0 " +
                         "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "1:0 | 4:0 | 5:0 # #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "0:0 | 2:0 | 3:0 # #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "0:0 | 2:0 | 3:0 " +
                         "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PointSemiOpenPolyline()
        {
            // Degenerate polylines are defined not contain any points under the
            // SEMI_OPEN model either, so again the point 3:0 and the degenerate
            // polyline "3:0, 3:0" do not intersect.
            //
            // The result does not depend on Options.polyline_loops_have_boundaries().
            Options options = new();
            options.PolylineModel_ = (PolylineModel.SEMI_OPEN);
            foreach (bool bool_value in new[] { false, true })
            {
                options.PolylineLoopsHaveBoundaries = (bool_value);
                var a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
                var b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
                ExpectResult(OpType.UNION, options, a, b,
                             "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
                ExpectResult(OpType.INTERSECTION, options, a, b,
                             "0:0 | 1:0 | 4:0 | 5:0 # #");
                ExpectResult(OpType.DIFFERENCE, options, a, b,
                             "2:0 | 3:0 # #");
                ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                             "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
            }
        }

        [Fact]
        public void Test_S2BooleanOperation_PointClosedPolyline()
        {
            // Under the CLOSED model, the degenerate polyline 3:0 does contain its
            // vertex.  Since polylines take precedence over points, the union of the
            // point 3:0 and the polyline 3:0 is the polyline only.  Similarly, since
            // subtracting a point from a polyline has no effect, the symmetric
            // difference includes only the polyline objects.
            //
            // The result does not depend on Options.polyline_loops_have_boundaries().
            Options options = new();
            options.PolylineModel_ = (PolylineModel.CLOSED);
            foreach (bool bool_value in new[] { false, true })
            {
                options.PolylineLoopsHaveBoundaries = (bool_value);
                var a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
                var b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
                ExpectResult(OpType.UNION, options, a, b,
                             "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
                ExpectResult(OpType.INTERSECTION, options, a, b,
                             "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #");
                ExpectResult(OpType.DIFFERENCE, options, a, b,
                             "# #");
                ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                             "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
            }
        }

        [Fact]
        public void Test_S2BooleanOperation_PointPolygonInterior()
        {
            Options options = new();  // PolygonModel is irrelevant.
                                                         // One interior point and one exterior point.
            var a = "1:1 | 4:4 # #";
            var b = "# # 0:0, 0:3, 3:0";
            ExpectResult(OpType.UNION, options, a, b,
                         "4:4 # # 0:0, 0:3, 3:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "1:1 # #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "4:4 # #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "4:4 # # 0:0, 0:3, 3:0");
        }

        [Fact]
        public void Test_S2BooleanOperation_PointOpenPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.OPEN);
            // See notes about the two vertices below.
            var a = "0:1 | 1:0 # #";
            var b = "# # 0:0, 0:1, 1:0";
            ExpectResult(OpType.UNION, options, a, b,
                         "0:1 | 1:0 # # 0:0, 0:1, 1:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "0:1 | 1:0 # #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "0:1 | 1:0 # # 0:0, 0:1, 1:0");
        }

        [Fact]
        public void Test_S2BooleanOperation_PointSemiOpenPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.SEMI_OPEN);
            // The two vertices are chosen such that the polygon contains one vertex but
            // not the other under PolygonModel.SEMI_OPEN.  (The same vertices are used
            // for all three PolygonModel options.)
            var polygon = MakePolygonOrDie("0:0, 0:1, 1:0");
            Assert.True(polygon.Contains(MakePointOrDie("0:1")));
            Assert.False(polygon.Contains(MakePointOrDie("1:0")));
            var a = "0:1 | 1:0 # #";
            var b = "# # 0:0, 0:1, 1:0";
            ExpectResult(OpType.UNION, options, a, b,
                         "1:0 # # 0:0, 0:1, 1:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "0:1 # #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "1:0 # #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "1:0 # # 0:0, 0:1, 1:0");
        }

        [Fact]
        public void Test_S2BooleanOperation_PointClosedPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.CLOSED);
            // See notes about the two vertices above.
            var a = "0:1 | 1:0 # #";
            var b = "# # 0:0, 0:1, 1:0";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:1, 1:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "0:1 | 1:0 # #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:1, 1:0");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexOpenPolylineVertex()
        {
            // Test first, last, and middle vertices of both polylines.  Also test
            // first/last and middle vertices of two polyline loops.
            //
            // Degenerate polylines are tested in PolylineEdgePolylineEdgeOverlap below.
            Options options = new();
            options.PolylineModel_ = (PolylineModel.OPEN);
            var a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
            var b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                     "| 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

            // The output consists of the portion of each input polyline that intersects
            // the opposite region, so the intersection vertex is present twice.  This
            // allows reassembling the individual polylins that intersect, if desired.
            // (Otherwise duplicates can be removed using DuplicateEdges.MERGE.)
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:1, 0:1 | 0:1, 0:1 #");

            // Note that all operations are defined such that subtracting a
            // lower-dimensional subset of an object has no effect.  In this case,
            // subtracting the middle vertex of a polyline has no effect.
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexOpenPolylineVertexLoopBoundariesFalse()
        {
            // With Options.polyline_loops_have_boundaries() == false, the 3 polyline
            // loops each have two vertices, both of which are contained.
            Options options = new();
            options.PolylineModel_ = (PolylineModel.OPEN);
            options.PolylineLoopsHaveBoundaries = (false);
            var a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
            var b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                     "| 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

            // Note that the polyline "0:3, 0:4, 0:3" only has two vertices, not three.
            // This means that 0:3 is emitted only once for that polyline, plus once for
            // the other polyline, for a total of twice.
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:1, 0:1 | 0:1, 0:1 " +
                         "| 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #");

            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexSemiOpenPolylineVertex()
        {
            // The result does not depend on Options.polyline_loops_have_boundaries().
            Options options = new();
            options.PolylineModel_ = (PolylineModel.SEMI_OPEN);
            foreach (bool bool_value in new[] { false, true })
            {
                options.PolylineLoopsHaveBoundaries = (bool_value);
                var a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
                var b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
                ExpectResult(OpType.UNION, options, a, b,
                             "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                             "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
                ExpectResult(OpType.INTERSECTION, options, a, b,
                             "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 " +
                             "| 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #");
                ExpectResult(OpType.DIFFERENCE, options, a, b,
                             "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
                ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                             "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                             "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
            }
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexClosedPolylineVertex()
        {
            Options options = new();
            options.PolylineModel_ = (PolylineModel.CLOSED);
            var a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
            var b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                     "| 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

            // Since Options.polyline_loops_have_boundaries() == true, the polyline
            // "0:3, 0:4, 0:3" has three vertices.  Therefore 0:3 is emitted twice for
            // that polyline, plus once for the other polyline, for a total of thrice.
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 " +
                         "| 0:2, 0:2 | 0:2, 0:2 " +
                         "| 0:3, 0:3 | 0:3, 0:3 | 0:3, 0:3 " +
                         "| 0:4, 0:4 | 0:4, 0:4 | 0:4, 0:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexClosedPolylineVertexLoopBoundariesFalse()
        {
            Options options = new();
            options.PolylineModel_ = (PolylineModel.CLOSED);
            options.PolylineLoopsHaveBoundaries = (false);
            var a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
            var b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                     "| 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

            // Since Options.polyline_loops_have_boundaries() == false, the polyline
            // "0:3, 0:4, 0:3" has two vertices.  Therefore 0:3 is emitted once for
            // that polyline, plus once for the other polyline, for a total of twice.
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 " +
                         "| 0:2, 0:2 | 0:2, 0:2 " +
                         "| 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 " +
                         "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_TestSemiOpenPolygonVerticesContained()
        {
            // Verify whether certain vertices of the test polygon are contained under
            // the semi-open boundary model (for use in the tests below).
            var polygon = MakePolygonOrDie(kVertexTestPolygonStr);
            Assert.True(polygon.Contains(MakePointOrDie("0:1")));
            Assert.True(polygon.Contains(MakePointOrDie("0:2")));
            Assert.True(polygon.Contains(MakePointOrDie("0:3")));
            Assert.True(polygon.Contains(MakePointOrDie("0:4")));
            Assert.False(polygon.Contains(MakePointOrDie("5:1")));
            Assert.False(polygon.Contains(MakePointOrDie("5:2")));
            Assert.False(polygon.Contains(MakePointOrDie("5:3")));
            Assert.False(polygon.Contains(MakePointOrDie("5:4")));
        }

        // Don't bother testing every PolylineModel with every PolygonModel for vertex
        // intersection, since we have already tested the PolylineModels individually
        // above.  It is sufficient to use PolylineModel.CLOSED with the various
        // PolygonModel options.
        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexOpenPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.OPEN);

            // Define some constants to reduce code duplication.
            // Test all combinations of polylines that start or end on a polygon vertex,
            // where the polygon vertex is open or closed using semi-open boundaries,
            // and where the incident edge is inside or outside the polygon.
            var a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 " +
                    "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #");
            var b = "# # " + kVertexTestPolygonStr;

            string kDifferenceResult =
                "# 0:1, 0:1 | 0:2, 0:2 | -1:3, 0:3 | 0:4, -1:4" +
                "| 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #";
            ExpectResult(OpType.UNION, options, a, b,
                         kDifferenceResult + kVertexTestPolygonStr);
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 0:1 | 0:2, 1:2 | 4:3, 5:3 | 5:4, 4:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult);
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         kDifferenceResult + kVertexTestPolygonStr);
        }

        // Like the test above, except that every polygon vertex is also incident to a
        // closed polyline vertex.  This tests that when an open vertex and a closed
        // vertex coincide with each other, the result is considered closed.
        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexOpenPolygonClosedPolylineVertex()
        {
            var kTestGeometrySuffix =
                "-2:0, 0:1 | -2:1, 0:2 | -2:2, 0:3 | -2:3, 0:4 | " +
                "7:0, 5:1 | 7:1, 5:2 | 7:2, 5:3 | 7:3, 5:4 # " +
                kVertexTestPolygonStr;

            Options options = new();
            options.PolygonModel_ = (PolygonModel.OPEN);
            var a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 " +
                      "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #");
            var b = ("# " + kTestGeometrySuffix);

            string kDifferencePrefix =
                "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2";
            ExpectResult(OpType.UNION, options, a, b,
                         kDifferencePrefix +
                         " | 0:1, 0:1 | 0:2, 0:2 | 5:3, 5:3 | 5:4, 5:4 | " +
                         kTestGeometrySuffix);
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4" +
                         "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4" +
                         "| 0:1, 0:1 | 0:2, 0:2 | 0:3, 0:3 | 0:4, 0:4" +
                         "| 5:1, 5:1 | 5:2, 5:2 | 5:3, 5:3 | 5:4, 5:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b, kDifferencePrefix + " #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         kDifferencePrefix + " | " + kTestGeometrySuffix);
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexSemiOpenPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.SEMI_OPEN);
            // Test all combinations of polylines that start or end on a polygon vertex,
            // where the polygon vertex is open or closed using semi-open boundaries,
            // and where the incident edge is inside or outside the polygon.
            //
            // The vertices at latitude 0 used below are all closed while the vertices
            // at latitude 5 are all open (see TestSemiOpenPolygonVerticesContained).
            var a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 " +
                      "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #");
            var b = "# # " + kVertexTestPolygonStr;
            string kDifferenceResult =
                "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #";
            ExpectResult(OpType.UNION, options, a, b,
                         kDifferenceResult + kVertexTestPolygonStr);
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4 " +
                         "| 4:3, 5:3 | 5:4, 4:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult);
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         kDifferenceResult + kVertexTestPolygonStr);
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineVertexClosedPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.CLOSED);
            // Test all combinations of polylines that start or end on a polygon vertex,
            // where the polygon vertex is open or closed using semi-open boundaries,
            // and where the incident edge is inside or outside the polygon.
            var a = ("# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 " +
                      "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #");
            var b = "# # " + kVertexTestPolygonStr;
            string kDifferenceResult =
                "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 #";
            ExpectResult(OpType.UNION, options, a, b,
                         kDifferenceResult + kVertexTestPolygonStr);
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4" +
                         "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult);
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         kDifferenceResult + kVertexTestPolygonStr);
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgePolylineEdgeCrossing()
        {
            // Two polyline edges that cross at a point interior to both edges.
            Options options = RoundToE(1);
            var a = "# 0:0, 2:2 #";
            var b = "# 2:0, 0:2 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 1:1 | 1:1, 1:1 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:0, 1:1, 2:2 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgePolylineEdgeOverlap()
        {
            // The PolylineModel does not affect this calculation.  In particular the
            // intersection of a degenerate polyline edge with itself is non-empty, even
            // though the edge contains no points in the OPEN and SEMI_OPEN models.
            Options options = new();
            options.PolygonModel_ = (PolygonModel.OPEN);
            // Test edges in the same and reverse directions, and degenerate edges.
            var a = "# 0:0, 1:0, 2:0, 2:5 | 3:0, 3:0 | 6:0, 5:0, 4:0 #";
            var b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0 #";
            // As usual, the expected output includes the relevant portions of *both*
            // input polylines.  Duplicates can be removed using GraphOptions.
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 1:0, 2:0, 2:5 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 " +
                         "| 6:0, 5:0, 4:0 | 4:0, 5:0 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:0, 1:0, 2:0 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 " +
                         "| 5:0, 4:0 | 4:0, 5:0 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 2:0, 2:5 | 6:0, 5:0 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 2:0, 2:5 | 6:0, 5:0 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineLoopMultipleOpenPolylineEdge()
        {
            // Here we test a polyline loop ABCA with the pairs {AA, AB} and {AA, AC}.
            // This tests not only what happens when degenerate polylines intersect loop
            // endpoints, but also what happens when polylines intersect a degenerate
            // and non-degenerate edge that overlap each other.
            Options options = new();
            options.PolylineModel_ = PolylineModel.OPEN;
            var a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
            var b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 " +
                         "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:0, 0:1 | 0:0, 0:1 | 2:2, 3:2 | 3:2, 2:2 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:1, 1:0, 0:0 | 0:0, 0:0 | 2:2, 2:3, 3:2 | 2:2, 2:2 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineLoopMultipleSemiOpenPolylineEdge()
        {
            // Like the test above but with SEMI_OPEN boundaries.  In this case ABCA
            // intersected with {AA, AB} is {AA, AB, AB} but ABCA intersected with {AA,
            // AC} is {AA, AA, AC, CA} since the chain ABCA contains its start vertex
            // but not its end vertex.
            Options options = new();
            options.PolylineModel_ = PolylineModel.SEMI_OPEN;
            var a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
            var b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 " +
                         "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:0, 0:0 | 0:0, 0:1 | 0:0, 0:1 " +
                         "| 2:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 | 3:2, 2:2 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineLoopMultipleClosedPolylineEdge()
        {
            // Like the test above but with CLOSED boundaries.  In this case ABCA
            // intersected with {AA, AB} is {AA, AA, AB, AB} since the chain ABCA
            // contains both its start vertex and end vertex.
            Options options = new();
            options.PolylineModel_ = (PolylineModel.CLOSED);
            var a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
            var b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 " +
                         "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 0:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 | 0:0, 0:1 " +
                         "| 2:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 | 3:2, 2:2 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineLoopMultiplePolylineEdgeLoopBoundariesFalse()
        {
            // Like the tests above but with polyline_loops_have_boundaries() == false.
            // In this case the result does not depend on the polyline model.  The
            // polyline AA intersects ABCA exactly once, and the intersection of ABCA
            // with {AA, AB} is {AA, AB, AB}.
            var polyline_models = new[]
            {
                PolylineModel.OPEN, PolylineModel.SEMI_OPEN, PolylineModel.CLOSED
            };
            foreach (var polyline_model in polyline_models)
            {
                Options options = new();
                options.PolylineModel_ = polyline_model;
                options.PolylineLoopsHaveBoundaries = false;
                var a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
                var b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
                ExpectResult(OpType.UNION, options, a, b,
                             "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 " +
                             "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
                ExpectResult(OpType.INTERSECTION, options, a, b,
                             "# 0:0, 0:0 | 0:0, 0:1 | 0:0, 0:1 " +
                             "| 2:2, 2:2 | 2:2, 3:2 | 3:2, 2:2 #");
                ExpectResult(OpType.DIFFERENCE, options, a, b,
                             "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
                ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                             "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
            }
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgeOpenPolygonEdgeOverlap()
        {
            Options options = new();
            options.PolygonModel_ = PolygonModel.OPEN;
            // A polygon and two polyline edges that coincide with the polygon boundary,
            // one in the same direction and one in the reverse direction.
            var a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
            var b = "# # 1:1, 1:3, 3:3, 3:1";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 1:1, 1:3, 3:3 | 3:3, 1:3 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgeSemiOpenPolygonEdgeOverlap()
        {
            var polygon = MakePolygonOrDie("1:1, 1:3, 3:3, 3:1");
            Assert.False(polygon.Contains(MakePointOrDie("1:1")));
            Assert.True(polygon.Contains(MakePointOrDie("1:3")));
            Assert.False(polygon.Contains(MakePointOrDie("3:3")));
            Assert.False(polygon.Contains(MakePointOrDie("3:1")));
            Options options = new();
            options.PolygonModel_ = (PolygonModel.SEMI_OPEN);
            var a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
            var b = "# # 1:1, 1:3, 3:3, 3:1";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:3, 1:3 | 1:1, 1:3, 3:3 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgeClosedPolygonEdgeOverlap()
        {
            Options options = new();
            options.PolygonModel_ = PolygonModel.CLOSED;
            var a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
            var b = "# # 1:1, 1:3, 3:3, 3:1";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 1:1, 1:3, 3:3, 3:1");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 1:3, 3:3 | 3:3, 1:3 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 1:1, 1:3, 3:3, 3:1");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonVertexMatching()
        {
            // This test shows that CrossingProcessor.ProcessEdgeCrossings() must set
            // a0_matches_polygon and a1_matches_polygon correctly even when (a0, a1)
            // itself is a polygon edge (or its sibling).  (It requires degenerate
            // polygon geometry to demonstrate this.)
            Options options = new();
            options.PolylineModel_ = (PolylineModel.CLOSED);
            options.PolygonModel_ = (PolygonModel.CLOSED);
            var a = "# 0:0, 1:1 # ";
            var b = "# # 0:0, 1:1";
            ExpectResult(OpType.UNION, options, a, b, "# # 0:0, 1:1");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgePolygonInterior()
        {
            Options options = new();  // PolygonModel is irrelevant.
            // One normal and one degenerate polyline edge in the polygon interior, and
            // similarly for the polygon exterior.
            var a = "# 1:1, 2:2 | 3:3, 3:3 | 6:6, 7:7 | 8:8, 8:8 # ";
            var b = "# # 0:0, 0:5, 5:5, 5:0";
            ExpectResult(OpType.UNION, options, a, b,
                         "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# 1:1, 2:2 | 3:3, 3:3 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 6:6, 7:7 | 8:8, 8:8 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEdgeIsolatedStartVertexPlusInteriorCrossing()
        {
            // Tests a polyline XYZ that when intersected with a polygon results in an
            // isolated vertex X plus a clipped portion UYV.  This case is unusual
            // because the isolated vertex is handled by creating a separate degenerate
            // S2Builder input edge XX which is added before the actual edge XY, and the
            // crossing edge information needs to be associated with XY rather than XX
            // in order for GraphEdgeClipper to be able to do its work properly.  The
            // test is constructed such that if the crossings are incorrectly associated
            // with the degenerate edge XX then not only will the output be incorrect,
            // it will also trigger an internal S2_DCHECK.
            var options = RoundToE(1);
            var a = "# 0:0, 0:10, 0:4 # ";
            var b = "# # 0:0, -5:5, 5:5";
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# 0:0, 0:0 | 0:5, 0:10, 0:5 #");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonEdgeIsolatedStartVertexPlusInteriorCrossing()
        {
            // Similar to the case above, but tests a polygon XYZ rather than a
            // polyline.  This requires using the CLOSED polygon model and computing the
            // intersection with a clockwise loop rather than subtracting a CCW loop.
            // The test is constructed such that if the crossings for the edge 0:0, 0:8
            // are incorrectly associated with the degenerate edge 0:0, then not only
            // will the output be incorrect, it will also trigger an internal S2_DCHECK.
            var options = RoundToE(1);
            options.PolygonModel_ = PolygonModel.CLOSED;
            var a = "# # 0:0, 5:5, -5:5";
            var b = "# # 1:4, 0:0, 0:8";
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 0:0; 0:5, 0:8, 0.8:5");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonVertexOpenPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = PolygonModel.OPEN;
            var a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
            var b = "# # 0:0, 5:3, 5:2";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonVertexSemiOpenPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = PolygonModel.SEMI_OPEN;
            var a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
            var b = "# # 0:0, 5:3, 5:2";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonVertexClosedPolygonVertex()
        {
            Options options = new();
            options.PolygonModel_ = PolygonModel.CLOSED;
            var a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
            var b = "# # 0:0, 5:3, 5:2";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 0:0");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
            ExpectResult(OpType.DIFFERENCE, options, b, a,
                         "# # 0:0, 5:3, 5:2");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonEdgePolygonEdgeCrossing()
        {
            // Two polygons whose edges cross at points interior to both edges.
            Options options = RoundToE(2);
            var a = "# # 0:0, 0:2, 2:2, 2:0";
            var b = "# # 1:1, 1:3, 3:3, 3:1";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:2, 1:2, 1:3, 3:3, 3:1, 2:1, 2:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 1:1, 1:2, 2:2, 2:1");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0; " +
                         "1:2, 1:3, 3:3, 3:1, 2:1, 2:2");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonEdgeOpenPolygonEdgeOverlap()
        {
            Options options = new();
            // One shape is a rectangle, the other consists of one triangle inside the
            // rectangle and one triangle outside the rectangle, where each triangle
            // shares one edge with the rectangle.  This implies that the edges are in
            // the same direction in one case and opposite directions in the other case.
            options.PolygonModel_ = PolygonModel.OPEN;
            var a = "# # 0:0, 0:4, 2:4, 2:0";
            var b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0; 0:4, 1:5, 2:4");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 0:0, 1:1, 2:0");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0, 1:1");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonEdgeSemiOpenPolygonEdgeOverlap()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.SEMI_OPEN);
            var a = "# # 0:0, 0:4, 2:4, 2:0";
            var b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:4, 1:5, 2:4, 2:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 0:0, 1:1, 2:0");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0, 1:1");
            // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are
            // normalized, i.e. the output could contain siblings pairs (which can be
            // discarded using S2Builder.GraphOptions).
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonEdgeClosedPolygonEdgeOverlap()
        {
            Options options = new();
            options.PolygonModel_ = (PolygonModel.CLOSED);
            var a = "# # 0:0, 0:4, 2:4, 2:0";
            var b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:4, 1:5, 2:4, 2:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 0:0, 1:1, 2:0; 0:4, 2:4");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0, 1:1");
            // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are
            // normalized, i.e. the output could contain siblings pairs.
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonPolygonInterior()
        {
            Options options = new();  // PolygonModel is irrelevant.
            // One loop in the interior of another polygon and one loop in the exterior.
            var a = "# # 0:0, 0:4, 4:4, 4:0";
            var b = "# # 1:1, 1:2, 2:2, 2:1; 5:5, 5:6, 6:6, 6:5";
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:0, 0:4, 4:4, 4:0; 5:5, 5:6, 6:6, 6:5");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 1:1, 1:2, 2:2, 2:1");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1; " +
                         "5:5, 5:6, 6:6, 6:5");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolygonEdgesDegenerateAfterSnapping()
        {
            Options options = RoundToE(0);
            // Two narrow rectangles forming a plus sign.
            var a = "# # 0:-1, 0:1, 0.1:1, 0.1:-1";
            var b = "# # -1:0.1, 1:0.1, 1:0, -1:0";
            // When snapping causes an output edge to become degenerate, it is still
            // emitted (since otherwise loops that contract to a single point would be
            // lost).  If the output layer doesn't want such edges, they can be removed
            // via DegenerateEdges.DISCARD.
            ExpectResult(OpType.UNION, options, a, b,
                         "# # 0:-1, 0:0, 0:1, 0:0 | " +
                         "-1:0, 0:0, 1:0, 0:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                         "# # 0:0");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                         "# # 0:-1, 0:0, 0:1, 0:0");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                         "# # 0:-1, 0:0, 0:1, 0:0 | " +
                         "-1:0, 0:0, 1:0, 0:0");
        }

        // The following class comprehensively tests the handling of degenerate
        // geometry.  It is used to implement over 4,000 individual test cases encoded
        // as a series of textual tables.
        private class DegeneracyCoverageTest //: testing.Test
        {
            private readonly ITestOutputHelper _logger;

            public DegeneracyCoverageTest(ITestOutputHelper logger) { _logger = logger; }

            // Verifies that the S2BooleanOperation results for the given OpType and
            // PolygonModel match the given set of rules (encoded as described below).
            public void Run(OpType op_type, PolygonModel polygon_model,
                List<string> rules)
            {
                Assert.Equal(rules.Count, kInputChars.Length);

                // For the symmetric operators (i.e., all except difference) we only need to
                // define half of the expected results.
                var symmetric = op_type != OpType.DIFFERENCE;

                // Possible values for Options.polyline_model().
                List<PolylineModel> kPolylineModels = new()
                {
                    PolylineModel.OPEN,
                    PolylineModel.SEMI_OPEN,
                    PolylineModel.CLOSED
                };

                // Possible values for Options.polyline_loops_have_boundaries().
                List<bool> kPolylineLoopOptions = new() { true, false };

                // The set of characters representing polyline inputs.
                const string kLineChars = "Pud";

                // The following nested loops iterate over all combinations of:
                //  - Input character 0
                //  - Input character 1
                //  - polyline model (if either input is a polyline)
                //  - whether a closed polyline loop is considered to have a boundary
                //    (if either input is a degenerate polyline)
                Options options = new();
                options.PolygonModel_ = polygon_model;
                for (int i = 0; i < kInputChars.Length; ++i)
                {
                    char ch0 = kInputChars[i];
                    var row = rules[i].Split(' ', StringSplitOptions.RemoveEmptyEntries).ToList();
                    // Verify and remove the row label.
                    Assert.Equal(ch0.ToString(), row[0]);
                    Assert.Equal("|", row[1]);
                    row.RemoveRange(0, 2);
                    var limit = symmetric ? (i + 1) : kInputChars.Length;
                    Assert.Equal(row.Count, limit);
                    for (int j = 0; j < limit; ++j)
                    {
                        char ch1 = kInputChars[j];
                        // Only iterate over polyline models if at least one input is a polyline.
                        int num_line_models = (kLineChars.Contains(ch0) ||
                            kLineChars.Contains(ch1)) ? 3 : 1;
                        for (int k = 0; k < num_line_models; ++k)
                        {
                            options.PolylineModel_ = kPolylineModels[k];
                            // Only iterate over polyline loop boundary options if at least one
                            // input is a degenerate polyline.
                            int num_polyline_loop_options = (ch0 == 'P' || ch1 == 'P') ? 2 : 1;
                            for (int m = 0; m < num_polyline_loop_options; ++m)
                            {
                                options.PolylineLoopsHaveBoundaries = kPolylineLoopOptions[m];
                                var code = row[j];
                                // Process any '<' or '>' operators in the result.
                                var choices = code.Split("<>", StringSplitOptions.RemoveEmptyEntries);
                                var result = choices[0];
                                if (choices.Length > 1)
                                {
                                    Assert.Equal(2, choices.Length);
                                    // Test whether each input contains the point A.  Note that we
                                    // can't use S2ContainsPointQuery because the containment test
                                    // must be done using the given S2BooleanOperation options.
                                    bool in0 = Contains(MakeIndex(ch0.ToString()),
                                                        MakeIndex("p"), options);
                                    bool in1 = Contains(MakeIndex(ch1.ToString()),
                                        MakeIndex("p"), options);
                                    // If the point containment conditions specified by the operators
                                    // '<' and '>' are satisfied then the result is choices[0],
                                    // otherwise it is choices[1].
                                    if ((code.Contains('<') && !in0) ||
                                        (code.Contains('>') && !in1))
                                    {
                                        result = choices[1];
                                    }
                                }
                                // Next process any '|' operators in the result.
                                choices = result.Split('|');
                                if (choices.Length > 1)
                                {
                                    Assert.Equal(3, num_line_models); //<< "No polylines present";
                                    Assert.Equal(3, choices.Length);
                                    result = choices[k];
                                }
                                TestResult(op_type, options, ch0, ch1, result);
                                if (symmetric && j != i)
                                {
                                    TestResult(op_type, options, ch1, ch0, result);
                                }
                            }
                        }
                    }
                }
            }

            // The inputs to the test cases are intended to span all possible types of
            // degenerate and non-degenerate edge configurations in the vicinity of an
            // individual input edge.  Essentially we are trying to cover all possible
            // behaviors of S2BooleanOperation.Impl.CrossingProcessor.ProcessEdge(),
            // which is responsible for determining which parts of each input edge
            // should be emitted to the output.
            //
            // For this purpose it is sufficient to consider two points A and B and the
            // possible types of degeneracies that involve just these two points.  The
            // actual locations of the points are immaterial, although for descriptive
            // purposes we suppose that B is located above A so that the edge AB is
            // considered "up" while the edge BA is considered "down".  Only one point
            // needs to be considered for point degeneracies, so we use A for this
            // purpose.  This means that all inputs consist of some subset of AA, AB,
            // and BA represented as edges of dimensions 0, 1, or 2.  For example, a
            // polyline edge intersected with a point located at one of its vertices
            // would be represented as the edge AB of dimension 1 in region 0 and the
            // edge AA of dimension 0 in region 1.
            //
            // Each possibility is represented as a single letter as follows (with the
            // corresponding edges in braces):
            //
            // Special:        . = empty {}
            //                 * = full sphere {}
            //
            // Dimension 0:    p = point {AA}
            //
            // Dimension 1:    P = point polyline {AA}
            //                 u = up edge {AB}
            //                 d = down edge {BA}
            //
            // Dimension 2:    s = point shell {AA}
            //                 S = sibling pair shell {AB, BA}
            //                 U = up edge {AB}
            //                 D = down edge {BA}
            //                 H = sibling pair hole {AB, BA}
            //                 h = point hole {AA}
            //
            // Using this encoding, the test case described above would be represented
            // as 'u' for region 0 and 'p' for region 1.  Note that while the test cases
            // focus only on what happens to the edges AB, BA, and AA, additional edges
            // are needed in order to construct valid inputs.  For example, the test
            // case U (the "up" polygon edge AB) is represented as a triangle ABC, while
            // the test case D (the "down" polygon edge BA) is represented as the
            // triangle AB(-C).
            //
            // The expected result of a given operation is then represented as a
            // sequence of the characters above.  In some cases additional symbols are
            // needed to define the expected results.  For example, in the closed
            // polygon model the union of U and D is the quadrilateral A(-C)BC which
            // does not appear in the list of symbols above.  So the results use the
            // following additional symbols:
            //
            //    Q   : the union of the "U" and "D" shapes (a quadrilateral)
            //    B   : a point polyline consisting only of the vertex B {BB}
            //
            // Finally, in some cases the expected result depends on the polyline model,
            // or on whether one of the two regions contains the vertex A.  For example,
            // the intersection of 'p' and 'U' in the semi-open model depends on whether
            // the representation of U (the triangle ABC mentioned above) contains its
            // vertex A.  Conditional results of this sort are encoded using the
            // following operators:
            //
            //    ~X  : the complement of X (where X must be either U or D)
            //   X<Y  : the result is X if region 0 contains point A, otherwise Y.
            //   X>Y  : the result is X if region 1 contains point A, otherwise Y.
            //  X|Y|Z : the result is X, Y, or Z depending on whether the polyline
            //          model is OPEN, SEMI_OPEN, or CLOSED.
            //
            // The operators above may be combined, e.g. X<>Y means the result is X if
            // both regions contain point A and Y otherwise.

            // These are the characters representing possible input regions, in the
            // order they are used in the rules listed further below.
            private const string kInputChars = ".pPudsSUDHh*";

            // Returns an S2ShapeIndex corresponding to the given string of characters.
            private static MutableS2ShapeIndex MakeIndex(string chars)
            {
                // The locations of A, B, C are arbitrary, however some tests are sensitive
                // as to whether certain polygons contain the points {A, B} or not.  If
                // these points are moved then some test results may need to be changed.
                S2Point a = new(1, 0, 0), b = new(0, 0, 1), c = new(0, 1, 0);
                var index = new MutableS2ShapeIndex();
                for (int i = 0; i < chars.Length; ++i)
                {
                    char ch = chars[i];
                    if (ch == '.')
                    {         // Empty
                    }
                    else if (ch == 'p')
                    {  // Point
                        index.Add(new S2PointVectorShape(new S2Point[] { a }));
                    }
                    else if (ch == 'P')
                    {  // Polyline consisting only of the point A
                        index.Add(new S2LaxPolylineShape(new S2Point[] { a, a }));
                    }
                    else if (ch == 'B')
                    {  // Polyline consisting only of the point B
                        index.Add(new S2LaxPolylineShape(new S2Point[] { b, b }));
                    }
                    else if (ch == 'u')
                    {  // Upwards polyline edge
                        index.Add(new S2LaxPolylineShape(new S2Point[] { a, b }));
                    }
                    else if (ch == 'd')
                    {  // Downwards polyline edge
                        index.Add(new S2LaxPolylineShape(new S2Point[] { b, a }));
                    }
                    else if (ch == 's')
                    {  // Point shell
                        index.Add(new S2LaxPolygonShape(new Loops() { new() { a } }));
                    }
                    else if (ch == 'S')
                    {  // Sibling pair shell
                        index.Add(new S2LaxPolygonShape(new Loops() { new() { a, b } }));
                    }
                    else if (ch == 'U')
                    {  // Upwards polygon edge
                        int i2 = index.Add(new S2LaxPolygonShape(new Loops() { new() { a, b, -c } }));
                        // Some test results require that the U polygon contains A but not B.
                        Assert.True(index.Shape(i2).ContainsBruteForce(a));
                        Assert.True(!index.Shape(i2).ContainsBruteForce(b));
                    }
                    else if (ch == 'D')
                    {  // Downwards polygon edge
                        int i2 = index.Add(new S2LaxPolygonShape(new Loops() { new() { b, a, c } }));
                        // Some test cases require that the D polygon excludes both A and B.
                        Assert.True(!index.Shape(i2).ContainsBruteForce(a));
                        Assert.True(!index.Shape(i2).ContainsBruteForce(b));
                    }
                    else if (ch == '~')
                    {  // Complement of following region (U or D)
                        ch = chars[++i];
                        if (ch == 'U')
                        {
                            index.Add(new S2LaxPolygonShape(new Loops() { new() { -c, b, a } }));
                        }
                        else if (ch == 'D')
                        {
                            index.Add(new S2LaxPolygonShape(new Loops { new() { c, a, b } }));
                        }
                        else
                        {
                            throw new Exception($"Unsupported character for ~ operator: {ch}");
                        }
                    }
                    else if (ch == 'Q')
                    {  // Union of 'U' and 'D' shapes
                        index.Add(new S2LaxPolygonShape(new Loops() { new() { a, c, b, -c } }));
                    }
                    else if (ch == 'H')
                    {  // Sibling pair hole
                        index.Add(new S2LaxPolygonShape(new Loops() { new() { a, b }, new() { } }));
                    }
                    else if (ch == 'h')
                    {  // Point hole
                        index.Add(new S2LaxPolygonShape(new Loops() { new() { a }, new() { } }));
                    }
                    else if (ch == '*')
                    {  // Full sphere
                        index.Add(new S2LaxPolygonShape(new Loops() { new() { } }));
                    }
                    else
                    {
                        throw new Exception($"Unknown degeneracy coverage symbol: {ch}");
                    }
                }
                return index;
            }

            // Verifies that S2BooleanOperation returns the given result for the inputs
            // represented by the characters 'ch0' and 'ch1'.
            private void TestResult(OpType op_type, Options options,
                       char ch0, char ch1, string result)
            {
                var index0 = MakeIndex(ch0.ToString());
                var index1 = MakeIndex(ch1.ToString());
                var expected = MakeIndex(result);
                _logger.WriteLine(@$"
  Operation: {op_type}
  PolygonModel: {options.PolygonModel_}
  PolylineModel: {options.PolylineModel_}
  polyline_loops_have_boundaries: {(options.PolylineLoopsHaveBoundaries ? "true" : "false")}
  Inputs: {ch0}, {ch1}
  Expected: {result}");
                ExpectResult(op_type, options, index0, index1, expected);
            }
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_OpenIntersection()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  .    pp",
                "P |  .    p<.   PP",
                "u |  .    p<.   PP<.  uu",
                "d |  .    p<.   PP<.  ud    dd",
                "s |  .     .     .     .     .     s",
                "S |  .     .     .     .     .     .     S",
                "U |  .     .     .     .     .     .     .     U",
                "D |  .     .     .     .     .     .     .     .     D",
                "H |  .     .     .     .     .     .     .     U     D     H",
                "h |  .     .     .     u     d     .     S     U     D     H     h",
                "* |  .     p     P     u     d     s     S     U     D     H     h     *",
            };
            dct.Run(OpType.INTERSECTION, PolygonModel.OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_SemiOpenIntersection()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  .    pp",
                "P |  .    p<.   PP",
                "u |  .    p<.   PP<.  uu",
                "d |  .    p<.   PP<.  ud    dd",
                "s |  .     .     .     .     .     s",
                "S |  .    p<.   P<.    .     .    s<.    S",
                "U |  .    p<.   P<.    u    P<>.  s<.    .     U",
                "D |  .    p<.   P<.   P<>.   d    s<.    .     .     D",
                "H |  .    p<.   P<.    u     d    s<.    .     U     D     H",
                "h |  .     p     P     u     d     .     S     U     D     H     h",
                "* |  .     p     P     u     d     s     S     U     D     H     h     *",
            };
            dct.Run(OpType.INTERSECTION, PolygonModel.SEMI_OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_ClosedIntersection()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  .    pp",
                "P |  .    p<.   PP",
                "u |  .    p<.   PP<.  uu",
                "d |  .    p<.   PP<.  ud    dd",
                "s |  .     p     P    P>.   P>.    s",
                "S |  .     p     P     u     d     s     S",
                "U |  .     p     P     u     d     s     S     U",
                "D |  .     p     P     u     d     s     S     S     D",
                "H |  .     p     P     u     d     s     S     U     D     H",
                "h |  .     p     P     u     d     s     S     U     D     H     h",
                "* |  .     p     P     u     d     s     S     U     D     H     h     *",
            };
            dct.Run(OpType.INTERSECTION, PolygonModel.CLOSED, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_OpenUnion()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  p    pp",
                "P |  P   P<Pp   PP",
                "u |  u   u<up   Pu    uu",
                "d |  d   d<dp   Pd    ud    dd",
                "s |  s    ps    Ps    us    ds     s",
                "S |  S    pS    PS    uS    dS     S     S",
                "U |  U    pU    PU    uU    dU     U     U     U",
                "D |  D    pD    PD    uD    dD     D     D    UD     D",
                "H |  H    pH    PH    uH    dH     H     H     H     H     H",
                "h |  h    ph    Ph   Ph>h  Ph>h    h     h     h     h     h     h",
                "* |  *     *     *     *     *     *     *     *     *     *     *     *",
            };
            dct.Run(OpType.UNION, PolygonModel.OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_SemiOpenUnion()
        {
            // CAVEAT: The results for (U,u) and (D,d) require the U polygon to contain
            // vertex A but not vertex B, and the D polygon to contain neither vertex.
            // This differs from most of the other tests, which encode the results
            // conditionally using the '<' and '>' operators.  That was not practical in
            // this case because (1) no conditional operators are defined for the 'B'
            // vertex and (2) encoding the full set of possibilites for all 12 cases
            // (i.e., the 3 polyline models and whether U contains A and/or B) would be
            // unwieldy.
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  p    pp",
                "P |  P   P<Pp   PP",
                "u |  u   u<up   Pu    uu",
                "d |  d   d<dp   Pd    ud    dd",
                "s |  s    ps    Ps    us    ds     s",
                "S |  S   S<pS  S<PS   uS    dS     S     S",
                "U |  U   U<pU  U<PU U|U|BU  dU     U     U     U",
                "D |  D   D<pD  D<PD   uD  D|BD|PBD D     D     Q     D",
                "H |  H   H<pH  H<PH    H     H     H     *     *     *     H",
                "h |  h     h     h     h     h     *    *>h   *>h   *>h   *>h    h",
                "* |  *     *     *     *     *     *     *     *     *     *     *     *",
            };
            dct.Run(OpType.UNION, PolygonModel.SEMI_OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_ClosedUnion()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  p    pp",
                "P |  P   P<Pp   PP",
                "u |  u   u<up   Pu    uu",
                "d |  d   d<dp   Pd    ud    dd",
                "s |  s     s     s    us    ds     s",
                "S |  S     S     S     S     S     S     S",
                "U |  U     U     U     U     U     U     U     U",
                "D |  D     D     D     D     D     D     D     Q     D",
                "H |  H     H     H     H     H     H     *     *     *     H",
                "h |  h     h     h     h     h     *     *     *     *     *     h",
                "* |  *     *     *     *     *     *     *     *     *     *     *     *",
            };
            dct.Run(OpType.UNION, PolygonModel.CLOSED, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_OpenDifference()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .     .     .     .     .     .     .     .     .     .     .     .",
                "p |  p     .    .>p   .>p   .>p    p     p     p     p     p     p     .",
                "P |  P     P     .    .>P   .>P    P     P     P     P     P     P     .",
                "u |  u     u     u     .   .|P|.   u     u     u     u     u    P<.    .",
                "d |  d     d     d   .|B|.   .     d     d     d     d     d    P<.    .",
                "s |  s     s     s     s     s     .     s     s     s     s     s     .",
                "S |  S     S     S     S     S     S     .     S     S     S     .     .",
                "U |  U     U     U     U     U     U     U     .     U     .     .     .",
                "D |  D     D     D     D     D     D     D     D     .     .     .     .",
                "H |  H     H     H     H     H     H     H    ~U    ~D     .     .     .",
                "h |  h     h     h     h     h     h     H    ~U    ~D     S     .     .",
                "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .",
            };
            dct.Run(OpType.DIFFERENCE, PolygonModel.OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_SemiOpenDifference()
        {
            // See SemiOpenUnion notes regarding (u,U) and (d,D).
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .     .     .     .     .     .     .     .     .     .     .     .",
                "p |  p     .    .>p   .>p   .>p    p     p    .>p   .>p    .     .     .",
                "P |  P     P     .    .>P   .>P    P     P    .>P   .>P    .     .     .",
                "u |  u     u     u     .   .|P|.   u     u   .|.|B   u     .     .     .",
                "d |  d     d     d   .|B|.   .     d     d     d   .|B|PB  .     .     .",
                "s |  s     s     s     s     s     .    .>s   .>s   .>s   .>s    s     .",
                "S |  S     S     S     S     S     S     .     .     .     S    s<.    .",
                "U |  U     U     U     U     U     U     U     .     U     .    s<.    .",
                "D |  D     D     D     D     D     D     D     D     .     .    s<.    .",
                "H |  H     H     H     H     H     H     H    ~U    ~D     .    s<.    .",
                "h |  h     h     h     h     h     h     H    ~U    ~D     S     .     .",
                "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .",
            };
            dct.Run(OpType.DIFFERENCE, PolygonModel.SEMI_OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_ClosedDifference()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .     .     .     .     .     .     .     .     .     .     .     .",
                "p |  p     .    .>p   .>p   .>p    .     .     .     .     .     .     .",
                "P |  P     P     .    .>P   .>P    .     .     .     .     .     .     .",
                "u |  u     u     u     .   .|P|.   u     .     .     .     .     .     .",
                "d |  d     d     d   .|B|.   .     d     .     .     .     .     .     .",
                "s |  s     s     s     s     s     .     .     .     .     .     s     .",
                "S |  S     S     S     S     S     S     .     .     .     S     .     .",
                "U |  U     U     U     U     U     U     U     .     U     .     .     .",
                "D |  D     D     D     D     D     D     D     D     .     .     .     .",
                "H |  H     H     H     H     H     H     H    ~U    ~D     .     .     .",
                "h |  h     h     h     h     h     h     H    ~U    ~D     S     .     .",
                "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .",
            };
            dct.Run(OpType.DIFFERENCE, PolygonModel.CLOSED, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_OpenSymmetricDifference()
        {
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  p     .",
                "P |  P   P<Pp    .",
                "u |  u   u<up  u<uP    .",
                "d |  d   d<dp  d<dP .|PB|.   .",
                "s |  s    sp    sP    su    sd     .",
                "S |  S    Sp    SP    Su    Sd     S     .",
                "U |  U    Up    UP    Uu    Ud     U     U     .",
                "D |  D    Dp    DP    Du    Dd     D     D    UD     .",
                "H |  H    Hp    HP    Hu    Hd     H     H    ~U    ~D     .",
                "h |  h    hp    hP   hP>h  hP>h    h     H    ~U    ~D     S     .",
                "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .",
            };
            dct.Run(OpType.SYMMETRIC_DIFFERENCE, PolygonModel.OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_SemiOpenSymmetricDifference()
        {
            // See SemiOpenUnion notes regarding (U,u) and (D,d).
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  p     .",
                "P |  P   P<Pp    .",
                "u |  u   u<up  u<uP    .",
                "d |  d   d<dp  d<dP .|PB|.   .",
                "s |  s    sp    sP    su    sd     .",
                "S |  S    Sp    SP    Su    Sd     S     .",
                "U |  U   U<Up  U<UP U|U|UB  Ud     U     U     .",
                "D |  D   D<Dp  D<DP   Du  D|BD|PBD D     D    UD     .",
                "H |  H     H     H     H     H     H     H    ~U    ~D     .",
                "h |  h     h     h     h     h     h     H    ~U    ~D     S     .",
                "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .",
            };
            dct.Run(OpType.SYMMETRIC_DIFFERENCE, PolygonModel.SEMI_OPEN, rules);
        }

        [Fact]
        public void Test_DegeneracyCoverageTest_ClosedSymmetricDifference()
        {
            // Note that (H,S)->H, (h,s)->h and (U,D)->UD.  In all three cases the
            // shared boundary is present on both sides and therefore these edges should
            // not be contained by the result, however this is not possible under the
            // CLOSED model.  The indicated results are the best approximation.
            List<string> rules = new()
            {
                //    .     p     P     u     d     s     S     U     D     H     h     *
                // |-----------------------------------------------------------------------
                ". |  .",
                "p |  p     .",
                "P |  P   P<Pp    .",
                "u |  u   u<up  u<uP    .",
                "d |  d   d<dp  d<dP .|PB|.   .",
                "s |  s     s     s    su    sd     .",
                "S |  S     S     S     S     S     S     .",
                "U |  U     U     U     U     U     U     U     .",
                "D |  D     D     D     D     D     D     D    UD     .",
                "H |  H     H     H     H     H     H     H    ~U    ~D     .",
                "h |  h     h     h     h     h     h     H    ~U    ~D     S     .",
                "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .",
            };
            dct.Run(OpType.SYMMETRIC_DIFFERENCE, PolygonModel.CLOSED, rules);
        }

        ///////////////////////////////////////////////////////////////////////////
        // The remaining tests are intended to cover combinations of features or
        // interesting special cases.

        [Fact]
        public void Test_S2BooleanOperation_ThreeOverlappingBars()
        {
            // Two vertical bars and a horizontal bar that overlaps both of the other
            // bars and connects them.

            // Round intersection points to E2 precision because the expected results
            // were computed in lat/lng space rather than using geodesics.
            var options = RoundToE(2);
            var a = "# # 0:0, 0:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3";
            var b = "# # 1:1, 1:4, 2:4, 2:1";
            ExpectResult(OpType.UNION, options, a, b,
                "# # 0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# # 1:1, 1:2, 2:2, 2:1; 1:3, 1:4, 2:4, 2:3");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; " +
                "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; " +
                "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3; " +
                "1:2, 1:3, 2:3, 2:2");
        }

        [Fact]
        public void Test_S2BooleanOperation_FourOverlappingBars()
        {
            // Two vertical bars and two horizontal bars.

            // Round intersection points to E2 precision because the expected results
            // were computed in lat/lng space rather than using geodesics.
            var options = RoundToE(2);
            var a = "# # 1:88, 1:93, 2:93, 2:88; -1:88, -1:93, 0:93, 0:88";
            var b = "# # -2:89, -2:90, 3:90, 3:89; -2:91, -2:92, 3:92, 3:91";
            ExpectResult(OpType.UNION, options, a, b,
                "# # -1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, " +
                "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, " +
                "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88; " +
                "0:90, 1:90, 1:91, 0:91" /*CW*/ );
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# # 1:89, 1:90, 2:90, 2:89; 1:91, 1:92, 2:92, 2:91; " +
                "-1:89, -1:90, 0:90, 0:89; -1:91, -1:92, 0:92, 0:91");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# # 1:88, 1:89, 2:89, 2:88; 1:90, 1:91, 2:91, 2:90; " +
                "1:92, 1:93, 2:93, 2:92; -1:88, -1:89, 0:89, 0:88; " +
                "-1:90, -1:91, 0:91, 0:90; -1:92, -1:93, 0:93, 0:92");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# # 1:88, 1:89, 2:89, 2:88; -1:88, -1:89, 0:89, 0:88; " +
                "1:90, 1:91, 2:91, 2:90; -1:90, -1:91, 0:91, 0:90; " +
                "1:92, 1:93, 2:93, 2:92; -1:92, -1:93, 0:93, 0:92; " +
                "-2:89, -2:90, -1:90, -1:89; -2:91, -2:92, -1:92, -1:91; " +
                "0:89, 0:90, 1:90, 1:89; 0:91, 0:92, 1:92, 1:91; " +
                "2:89, 2:90, 3:90, 3:89; 2:91, 2:92, 3:92, 3:91");
        }

        [Fact]
        public void Test_S2BooleanOperation_OverlappingDoughnuts()
        {
            // Two overlapping square doughnuts whose holes do not overlap.
            // This means that the union polygon has only two holes rather than three.

            // Round intersection points to E2 precision because the expected results
            // were computed in lat/lng space rather than using geodesics.
            var options = RoundToE(1);
            var a = "# # -1:-93, -1:-89, 3:-89, 3:-93; " +
                                "0:-92, 2:-92, 2:-90, 0:-90" /*CW*/;
            var b = "# # -3:-91, -3:-87, 1:-87, 1:-91; " +
                                "-2:-90, 0:-90, 0:-88, -2:-88" /*CW*/;
            ExpectResult(OpType.UNION, options, a, b,
                "# # -1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93; " +
                "0:-92, 2:-92, 2:-90, 1:-90, 1:-91, 0:-91; " + /*CW */
                "-2:-90, -1:-90, -1:-89, 0:-89, 0:-88, -2:-88" /* CW */ );
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# # -1:-91, -1:-90, 0:-90, 0:-91; " +
                "0:-90, 0:-89, 1:-89, 1:-90");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, " +
                "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; " +
                "-1:-90, -1:-89, 0:-89, 0:-90");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, " +
                "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; " +
                "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88,-2:-88,-2:-90,-1:-90,-1:-91; " +
                "-1:-90, -1:-89, 0:-89, 0:-90; " +
                "1:-91, 0:-91, 0:-90, 1:-90");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineEnteringRectangle()
        {
            // A polyline that enters a rectangle very close to one of its vertices.
            var options = RoundToE(1);
            var a = "# 0:0, 2:2 #";
            var b = "# # 1:1, 1:3, 3:3, 3:1";
            ExpectResult(OpType.UNION, options, a, b,
                "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# 1:1, 2:2 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# 0:0, 1:1 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1");
        }

        [Fact]
        public void Test_S2BooleanOperation_PolylineCrossingRectangleTwice()
        {
            // A polyline that crosses a rectangle in one direction, then moves to a
            // different side and crosses the rectangle in the other direction.  Note
            // that the input polyline has a self-intersection and that an extra vertex
            // is *not* added at the intersection point.  This is an important feature
            // since it allows polylines with many self-intersections (such as GPS
            // tracks) to be manipulated without possible quadratic size increases.
            var options = RoundToE(1);
            var a = "# 0:-5, 0:5, 5:0, -5:0 #";
            var b = "# # 1:1, 1:-1, -1:-1, -1:1";
            ExpectResult(OpType.UNION, options, a, b,
                "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 " +
                "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# 0:-1, 0:1 | 1:0, -1:0 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 " +
                "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1");
        }

        // This test demonstrated that S2 geometry can easily be transformed such that
        // no edge crosses the 180 degree meridian, as required by formats such as
        // GeoJSON, by simply subtracting a degenerate loop that follows the 180 degree
        // meridian.  This not only splits polylines along the meridian, it also inserts
        // the necessary extra vertices at the north/south poles.  (The only extra step
        // is that the vertices along the 180 degree meridian or at the poles may need
        // to be "doubled" into two vertices, one at longitude 180 and one at longitude
        // -180, in order to match the longitudes of the adjacent vertices.)
        [Fact]
        public void Test_S2BooleanOperation_MeridianSplitting()
        {
            // A line along the equator crossing the 180 degree meridian.
            TestMeridianSplitting("# 0:-160, 0:170 #", "# 0:-160, 0:180, 0:170 #");

            // The northern hemisphere.
            TestMeridianSplitting("# # 0:0, 0:120, 0:-120",
                        "# # 90:0, 0:180, 0:-120, 0:0, 0:120, 0:180");

            // A small square that crosses the 180th meridian.  Notice that one input
            // loop is split into two output loops.
            TestMeridianSplitting(
                "# # 9:179, 9:-179, 10:-179, 10:179",
                "# # 9.00134850712993:180, 9:-179, 10:-179, 10.0014925269841:180; " +
                "10.0014925269841:180, 10:179, 9:179, 9.00134850712993:180");

            // An annulus that crosses the 180th meridian.  This turns into two shells.
            TestMeridianSplitting(
                "# # 8:178, 8:-178, 11:-178, 11:178; 9:179, 10:179, 10:-179, 9:-179",
                "# # 10.0014925269841:180, 10:-179, 9:-179, 9.00134850712993:180, " +
                "8.00481316618607:180, 8:-178, 11:-178, 11.00654129428:180; " +
                "9.00134850712993:180, 9:179, 10:179, 10.0014925269841:180, " +
                "11.00654129428:180, 11:178, 8:178, 8.00481316618607:180");

            // An annulus that crosses the 180th meridian.  This turns into two shells.
            TestMeridianSplitting(
                "# # 8:178, 8:-178, 11:-178, 11:178; 9:179, 10:179, 10:-179, 9:-179",
                "# # 10.0014925269841:180, 10:-179, 9:-179, 9.00134850712993:180, " +
                "8.00481316618607:180, 8:-178, 11:-178, 11.00654129428:180; " +
                "9.00134850712993:180, 9:179, 10:179, 10.0014925269841:180, " +
                "11.00654129428:180, 11:178, 8:178, 8.00481316618607:180");

            // The whole world except for a small square that crosses the 180th meridian.
            // This is a single loop that visits both poles.  The result is correct
            // except that (1) +180 or -180 needs to be chosen consistently with the
            // adjacent points, and (2) each pole needs to be duplicated (once with
            // longitude -180 and once with longitude 180).
            TestMeridianSplitting(
                "# # 9:-179, 9:179, 10:179, 10:-179",
                "# # 0:180, 9.00134850712993:180, 9:179, 10:179, 10.0014925269841:180, " +
                "90:0, 10.0014925269841:180, 10:-179, 9:-179, 9.00134850712993:180, " +
                "0:180, -90:0");
        }


        private static void ComputeTestUnion(Loops a_loops,
                              Loops b_loops,
                              S1Angle snap_radius, S2LaxPolygonShape result)
        {
            MutableS2ShapeIndex a = new(), b = new();
            a.Add(new S2LaxPolygonShape(a_loops));
            b.Add(new S2LaxPolygonShape(b_loops));
            S2BooleanOperation op = new(OpType.UNION,
                                  new LaxPolygonLayer(result),
                                  new Options(
                                      new IdentitySnapFunction(snap_radius)));
            Assert.True(op.Build(a, b, out _));// << error;
            Assert.False(result.IsEmpty());
            //    << "\nS2Polygon: " << s2textformat.ToString(a)
            //    << "\nS2Polygon: " << s2textformat.ToString(b);
        }

        [Fact]
        public void Test_S2BooleanOperation_GetCrossedVertexIndexBug1()
        {
            // This test exercises a rare special case in GetCrossedVertexIndex where
            // two crossing edge chains snap to a different permutation of the same
            // vertices.  In this example one input edge crosses another edge from right
            // to left, the first edge snaps to BCD and the second snaps to ABDC, and
            // triangle BCD is CCW.  Since BCD is to the right of BD, this means that
            // the first edge has not yet crossed the second at vertex B, leaving C or D
            // as the possible crossing vertices.
            Loops a_loops = new()
            {
                new()
                {
                    new(-0.38306437985388492, -0.74921955334206214, 0.54030708099846292),
                    new(-0.3830643798552798, -0.74921955334134249, 0.5403070809984718),
                    new(-0.38306437985529124, -0.74921955334136414, 0.54030708099843361),
                    new(-0.38306437985389635, -0.74921955334208379, 0.54030708099842473),
                },
            };
            Loops b_loops = new()
            {
                new()
                {
                    new(-0.38306437985390962, -0.74921955334210588, 0.54030708099838465),
                    new(-0.38306437985527797, -0.74921955334134205, 0.54030708099847369),
                    new(-0.38306437985527941, -0.74921955334134405, 0.54030708099847014),
                    new(-0.38306437985391095, -0.74921955334210777, 0.54030708099838098),
                },
            };
            S2LaxPolygonShape result = new();
            ComputeTestUnion(a_loops, b_loops, S2.kIntersectionMergeRadiusS1Angle,
                             result);
        }

        [Fact]
        public void Test_S2BooleanOperation_GetCrossedVertexIndexBug2()
        {
            // This test exercises another rare case where the crossing vertices chosen
            // by GetCrossedVertexIndex() are not ordered correctly along the edge being
            // crossed.  This is handled by adding extra edges to the output in order to
            // link up the crossings in the correct order.
            Loops a_loops = new()
            {
                new()
                {
                    new(-0.3837392878495085, -0.7477800800281974, 0.5418201831546835),
                    new(-0.38373928785696076, -0.7477800800212292, 0.54182018315902258),
                    new(-0.38373928785701278, -0.74778008002124685, 0.5418201831589613),
                    new(-0.38373928785703426, -0.7477800800212544, 0.54182018315893576),
                    new(-0.38373947205489456, -0.74778014227795497, 0.5418199667802881),
                    new(-0.38373947204434411, -0.74778014228781997, 0.54181996677414512),
                    new(-0.38373947205872994, -0.74778014228185352, 0.54181996677219124),
                    new(-0.38373947218468357, -0.74778014288930306, 0.54181996584462788),
                    new(-0.3837396702525171, -0.74778021044361542, 0.54181973233114322),
                    new(-0.38373967023137123, -0.74778021046333043, 0.54181973231891067),
                    new(-0.38373947216030285, -0.74778014290791484, 0.54181996583620895),
                    new(-0.38373947217087578, -0.74778014289805739, 0.54181996584232528),
                    new(-0.38373947215649007, -0.74778014290402395, 0.54181996584427927),
                    new(-0.3837394720305386, -0.74778014229658485, 0.5418199667718262),
                    new(-0.38373928783585998, -0.74778008004095942, 0.54182018314673686),
                    new(-0.38373928784641037, -0.7477800800310942, 0.54182018315287972),
                    new(-0.38373928783578648, -0.74778008004093421, 0.54182018314682368),
                    new(-0.383739287835765, -0.74778008004092666, 0.54182018314684921),
                },
            };
            Loops b_loops = new()
            {
                new()
                {
                    new(-0.38373923813692823, -0.7477800632164362, 0.54182024156551456),
                    new(-0.3837392878569364, -0.74778008002122087, 0.54182018315905123),
                    new(-0.38373928784640354, -0.74778008003106944, 0.54182018315291858),
                    new(-0.38373928784638789, -0.74778008003108642, 0.54182018315290648),
                    new(-0.38373928784638023, -0.74778008003109453, 0.54182018315290048),
                    new(-0.38373928783692102, -0.74778008004124585, 0.54182018314559),
                    new(-0.38373928783691913, -0.74778008004124541, 0.54182018314559188),
                    new(-0.38373928784636568, -0.74778008003110774, 0.54182018315289271),
                    new(-0.38373928784637329, -0.74778008003109953, 0.54182018315289848),
                    new(-0.38373928783583561, -0.74778008004095109, 0.5418201831467655),
                    new(-0.38373923811582744, -0.74778006323616641, 0.54182024155322883),
                    new(-0.38373857650312843, -0.74777983961840766, 0.54182101875399913),
                    new(-0.38373857652422921, -0.74777983959867744, 0.54182101876628486),
                },
            };
            S2LaxPolygonShape result = new();
            ComputeTestUnion(a_loops, b_loops, S2.kIntersectionMergeRadiusS1Angle,
                             result);
        }

        [Fact]
        public void Test_S2BooleanOperation_GetCrossedVertexIndexBug3()
        {
            // This test exercise the special case in GetCrossedVertexIndex() that
            // requires checking the orientation of a loop.  This is done by adding up the
            // turning angles at each vertex, which in turn involves computing the edge
            // normals and measuring the angles between them.  However in this test, some
            // of the edge normals returned by S2.RobustCrossProd() used to be so small
            // that there were floating-point underflows when computing the angles between
            // them.  This was fixed by implementing the long-standing TODO of making
            // S2.RobustCrossProd() actually robust.
            Loops a_loops = new()
            {
                new()
                {
                    new(1, 0, 2.4678234835261742e-72),
                    new(0.99984769515639127, 0.017452406437283512, 1.8530922845942552e-27),
                    new(0.99740259703611311, 0.069881849826437858, 0.017452406437283512),
                },
            };
            Loops b_loops = new()
            {
                new()
                {
                    new(0.99999999999999989, 2.4674476220564615e-72, 2.4678234835261742e-72),
                    new(0.99999999999999989, 2.8837981406657438e-169, 2.4678234835261742e-72),
                    new(1, 2.8837981406657432e-169, 2.4678234835261742e-72),
                },
            };
            S2LaxPolygonShape result = new();
            ComputeTestUnion(a_loops, b_loops, S1Angle.Zero, result);
        }

        [Fact]
        public void Test_S2BooleanOperation_GetCrossedVertexIndexBug4()
        {
            // This example tests the "special case" in GetCrossedVertexIndex() in
            // situations where two edges snap to the same sequence of vertices in
            // different orders.  The first two edges (a0, a1) and (b0, b1) of the
            // following polygons cross such that after snapping, the corresponding edge
            // chains are:
            //
            //   a0 a1 . a0 b0 b1 x a1
            //   b0 b1 . b0 x b1
            //
            // where "x" is the computed intersection point of (a0, a1) and (b0, b1).
            // Previously there was a bug such that the two edge chains did not choose
            // the same vertex to represent the point where the two chains cross: the
            // (a0, a1) chain chose "x" as the crossing point while the (b0, b1) chain
            // chose "b0".  This has been fixed such that both chains now choose "x".
            // (Both "x" and "b1" happen to be valid choices in this example, but it is
            // essential that both subchains make the same choice.)

            // S2LatLng coordinates are not accurate enough to reproduce this example.
            var a_loops = new Loops{new List<S2Point>{
      // 51.5131559470858:-0.130381523356724
      new S2Point(0.62233331065911901, -0.0014161759526823048, 0.78275107466533156),
      // 51.5131892038956:-0.130404244210776
      new S2Point(0.6223328557578689, -0.0014164217071954736, 0.78275143589379825),
      MakePointOrDie("51.51317:-0.1306")
    }};
            var b_loops = new Loops{new List<S2Point>{
      // 51.5131559705551:-0.13038153939079
      new S2Point(0.62233331033809591, -0.001416176126110953, 0.78275107492024998),
      // 51.5131559705551:-0.130381539390786
      new S2Point(0.62233331033809591, -0.0014161761261109063, 0.78275107492025009),
      MakePointOrDie("51.52:-0.12"),
      MakePointOrDie("51.52:-0.14")
    }};

            S2LaxPolygonShape result = new();
            ComputeTestUnion(a_loops, b_loops, S1Angle.Zero, result);
        }

        [Fact]
        public void Test_S2BooleanOperation_FullAndEmptyResults()
        {
            // The followingants are all in S2TextFormat.MakeLaxPolygon() format.
            string kEmpty = "";
            string kFull = "full";

            // Two complementary shell/hole pairs, together with alternative shells that
            // are slightly smaller or larger than the original.
            string kShell1 = "10:0, 10:10, 20:10";
            string kHole1 = "10:0, 20:10, 10:10";
            string kShell1Minus = "11:2, 11:9, 18:9";
            string kShell1Plus = "9:-2, 9:11, 22:11";
            string kShell2 = "10:20, 10:30, 20:30";
            string kHole2 = "10:20, 20:30, 10:30";

            // The northern and southern hemispheres.
            string kNorthHemi = "0:0, 0:120, 0:-120";
            string kSouthHemi = "0:0, 0:-120, 0:120";
            // These edges deviate from kSouthHemi by slightly more than 1 degree.
            string kSouthHemiPlus = "0.5:0, 0.5:-120, 0.5:120";

            // A shell and hole that cover complementary hemispheres, such that each
            // hemisphere intersects all six S2 cube faces.  There are also alternative
            // shells that are slightly smaller or larger than the original.
            string k6FaceShell1 = "0:-45, 45:0, 45:90, 0:135, -45:180, -45:-90";
            string k6FaceHole1 = "0:-45, -45:-90, -45:180, 0:135, 45:90, 45:0";
            string k6FaceShell1Minus = "-1:-45, 44:0, 44:90, -1:135, -46:180, -46:-90";
            string k6FaceShell1Plus = "1:-45, 46:0, 46:90, 1:135, -44:180, -44:-90";

            // Two complementary shell/hole pairs that are small enough so that they will
            // disappear when the snap radius chosen above is used.
            string kAlmostEmpty1 = "2:0, 2:10, 3:0";
            string kAlmostFull1 = "2:0, 3:0, 2:10";
            string kAlmostEmpty2 = "4:0, 4:10, 5:0";
            string kAlmostFull2 = "4:0, 5:0, 4:10";

            // A polygon that intersects all 6 faces such but snaps to an empty polygon.
            string k6FaceAlmostEmpty1 = k6FaceShell1Minus + "; " + k6FaceHole1;

            // Test empty UNION results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.UNION, kEmpty, kEmpty, kEmpty);
            //  - Empty due to snapping, union does not intersect all 6 cube faces.
            ExpectPolygon(OpType.UNION, kAlmostEmpty1, kAlmostEmpty2, kEmpty);
            //  - Empty due to snapping, union intersects all 6 cube faces.
            ExpectPolygon(OpType.UNION, k6FaceAlmostEmpty1, k6FaceAlmostEmpty1, kEmpty);

            // Test full UNION results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.UNION, kEmpty, kFull, kFull);
            ExpectPolygon(OpType.UNION, kEmpty, kFull, kFull);
            ExpectPolygon(OpType.UNION, kFull, kFull, kFull);
            //  - Exact result, some input edges.
            ExpectPolygon(OpType.UNION, kFull, kShell1, kFull);
            ExpectPolygon(OpType.UNION, kHole1, kHole2, kFull);
            ExpectPolygon(OpType.UNION, kHole1, kShell1, kFull);
            //  - Full due to snapping, almost complementary polygons.
            ExpectPolygon(OpType.UNION, kHole1, kShell1Minus, kFull);
            ExpectPolygon(OpType.UNION, k6FaceHole1, k6FaceShell1Minus, kFull);

            // Test empty INTERSECTION results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.INTERSECTION, kEmpty, kEmpty, kEmpty);
            ExpectPolygon(OpType.INTERSECTION, kEmpty, kFull, kEmpty);
            ExpectPolygon(OpType.INTERSECTION, kFull, kEmpty, kEmpty);
            //  - Exact result, inputs do not both intersect all 6 cube faces.
            ExpectPolygon(OpType.INTERSECTION, kEmpty, kHole1, kEmpty);
            ExpectPolygon(OpType.INTERSECTION, kShell1, kShell2, kEmpty);
            ExpectPolygon(OpType.INTERSECTION, kShell1, kHole1, kEmpty);
            //  - Exact result, inputs both intersect all 6 cube faces.
            ExpectPolygon(OpType.INTERSECTION, k6FaceShell1, k6FaceHole1, kEmpty);
            //  - Empty due to snapping, inputs do not both intersect all 6 cube faces.
            ExpectPolygon(OpType.INTERSECTION, kShell1Plus, kHole1, kEmpty);
            //  - Empty due to snapping, inputs both intersect all 6 cube faces.
            ExpectPolygon(OpType.INTERSECTION, k6FaceShell1Plus, k6FaceHole1, kEmpty);

            // Test full INTERSECTION results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.INTERSECTION, kFull, kFull, kFull);
            //  - Full due to snapping, almost full input polygons.
            ExpectPolygon(OpType.INTERSECTION, kAlmostFull1, kAlmostFull2, kFull);

            // Test empty DIFFERENCE results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.DIFFERENCE, kEmpty, kEmpty, kEmpty);
            ExpectPolygon(OpType.DIFFERENCE, kEmpty, kFull, kEmpty);
            ExpectPolygon(OpType.DIFFERENCE, kFull, kFull, kEmpty);
            //  - Exact result, first input does not intersect all 6 cube faces.
            ExpectPolygon(OpType.DIFFERENCE, kEmpty, kShell1, kEmpty);
            ExpectPolygon(OpType.DIFFERENCE, kShell1, kFull, kEmpty);
            ExpectPolygon(OpType.DIFFERENCE, kShell1, kShell1, kEmpty);
            ExpectPolygon(OpType.DIFFERENCE, kShell1, kHole2, kEmpty);
            //  - Exact result, first input intersects all 6 cube faces.
            ExpectPolygon(OpType.DIFFERENCE, k6FaceShell1, k6FaceShell1Plus, kEmpty);
            //  - Empty due to snapping, first input does not intersect all 6 cube faces.
            ExpectPolygon(OpType.DIFFERENCE, kShell1Plus, kShell1, kEmpty);
            //  - Empty due to snapping, first input intersect all 6 cube faces.
            ExpectPolygon(OpType.DIFFERENCE, k6FaceShell1Plus, k6FaceShell1, kEmpty);

            // Test full DIFFERENCE results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.DIFFERENCE, kFull, kEmpty, kFull);
            //  - Full due to snapping, almost full/empty input polygons.
            ExpectPolygon(OpType.DIFFERENCE, kAlmostFull1, kAlmostEmpty2, kFull);

            // Test empty SYMMETRIC_DIFFERENCE results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kEmpty, kEmpty, kEmpty);
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kFull, kFull, kEmpty);
            //  - Exact result, union does not intersect all 6 cube faces.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1, kShell1, kEmpty);
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kNorthHemi, kEmpty);
            //  - Exact result, union intersects all 6 cube faces.  This case is only
            //    handled correctly due to the kBiasTowardsEmpty heuristic.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1, k6FaceShell1, kEmpty);
            //  - Empty due to snapping, union does not intersect all 6 cube faces.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1Plus, kShell1, kEmpty);
            //  - Empty due to snapping, union intersects all 6 cube faces.  This case is
            //    only handled correctly due to the kBiasTowardsEmpty heuristic.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Plus, k6FaceShell1, kEmpty);
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Minus, k6FaceShell1, kEmpty);

            // Test full SYMMETRIC_DIFFERENCE results.
            //  - Exact result, no input edges.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kFull, kEmpty, kFull);
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kEmpty, kFull, kFull);
            //  - Exact result, complementary input polygons.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1, kHole1, kFull);
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kAlmostEmpty1, kAlmostFull1, kFull);
            //  - Full due to snapping, almost complementary input polygons.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1Plus, kHole1, kFull);
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kAlmostFull1, kAlmostEmpty2, kFull);
            //  - Exact result, complementary hemispheres, at least one input does not
            //    intersect all 6 cube faces.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kSouthHemi, kFull);
            //  - Exact result, almost complementary hemispheres, at least one input does
            //    not intersect all 6 cube faces.
            ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kSouthHemiPlus, kFull);

            // TODO(ericv): The following case is not currently implemented.
            //  - Full result, complementary (to within the snap radius) input polygons
            //    each with an area of approximately 2*Pi, and both polygons intersect all
            //    6 cube faces.
#if Undefined_0
  ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1, k6FaceHole1, kFull);
  ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Plus, k6FaceHole1,kFull);
  ExpectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Minus, k6FaceHole1,kFull);
#endif
        }

        // Tests S2BooleanOperation.Equals, which computes the symmetric difference
        // between two geometries and tests whether the result is empty.
        //
        // This also indirectly tests IsEmpty(), which is used to implement Contains()
        // and Intersects().
        [Fact]
        public void Test_S2BooleanOperation_Equals()
        {
            Assert.True(TestEqual("# #", "# #"));
            Assert.True(TestEqual("# # full", "# # full"));

            Assert.False(TestEqual("# #", "# # full"));
            Assert.False(TestEqual("0:0 # #", "# #"));
            Assert.False(TestEqual("0:0 # #", "# # full"));
            Assert.False(TestEqual("# 0:0, 1:1 #", "# #"));
            Assert.False(TestEqual("# 0:0, 1:1 #", "# # full"));
            Assert.False(TestEqual("# # 0:0, 0:1, 1:0 ", "# #"));
            Assert.False(TestEqual("# # 0:0, 0:1, 1:0 ", "# # full"));
        }

        [Fact]
        public void Test_S2BooleanOperation_GetCrossedVertexIndexBug5()
        {
            // Yet another bizarre situation where two crossing edges snap (correctly) to
            // a sequence of vertices in different orders.  Using the internal vertex
            // numbers assigned by S2Builder, input edges 3 and 12 snap to the following
            // vertex sequences:
            //
            //   Input edge  3:  14, 8, 4, 9, 2, 5
            //   Input edge 12:   2, 7, 8, 9
            //
            // Furthermore input edge 3 crosses input edge 12 from left to right.
            // Schematically, here is what edge 12 crossing edge 3 looks like:
            //
            //   14-->--8-->--4-->--9-->--2-->--5
            //          |\         /     /
            //          \ \--->---/     /
            //           \             /
            //            \--<--7--<--/
            //
            // And here is what edge 3 crossing edge 12 looks like:
            //
            //             14-->--\   /---4->-\
            //                     \ /         \
            //          2-->--7-->--8----->-----9
            //         / \                     /
            //  5--<--/   \---------<---------/
            //
            // In both cases, the only possible choice of crossing vertex consistent with
            // the fact that edge 3 crosses edge 12 from left to right is vertex 9.
            // Determining this requires knowing that loop (9, 2, 7, 8) is clockwise
            // (the "special case" in GetCrossedVertexIndex).  The code previously didn't
            // have quite the correct test to decide when this was necessary.
            Loops a_loops = new()
            {
                new()
                {
                    new(0.99984769515639127, 0, 0.017452406437283512),
                    new(0.99923861495548261, 0.017441774902830158, 0.034899496702500969),
                    new(0.99847743863945992, 0.052327985223313139, 0.017452406437283512),
                    new(0.99802119662406841, 0.034851668155187324, 0.052335956242943835),
                },
            };
            Loops b_loops = new()
            {
                new()
                {
                    new(0.99802119662406841, 0.034851668155187324, 0.052335956242943835),
                    new(0.99619692339885657, 0.052208468483931986, 0.069756473744125302),
                    new(0.99802098681615425, 0.034839714972148959, 0.052347914334467859),
                    new(0.99741208276778681, 0.017411821260589495, 0.069756473744125302),
                    new(0.99741219210106513, 0.017411340538768819, 0.069755030419252628),
                    new(0.99741211642315963, 0.017409893252357169, 0.069756473744125302),
                    new(0.99984769515639116, 4.9500424645560228e-16, 0.017452406437284993),
                    new(0.99984769515639127, 3.7368529835165677e-16, 0.017452406437284632),
                    new(0.99984769515639116, 3.3065924905014365e-16, 0.017452406437284504),
                    new(0.99984769515639127, 9.9060035932242025e-16, 0.017452406437284504),
                    new(0.99969541350954794, 0.017449748351250485, 0.017452406437283512),
                },
                new()
                {
                    new(0.99984769515639116, 3.3065924905014365e-16, 0.017452406437284504),
                    new(0.99984769515639116, 3.3006856770496304e-16, 0.017452406437284504),
                    new(0.99984769515639127, 0, 0.017452406437284504),
                    new(0.99984769515639127, 0, 0.017452406437283512),
                },
            };
            S2LaxPolygonShape result = new();
            ComputeTestUnion(a_loops, b_loops, S1Angle.Zero, result);
        }

        [Fact]
        public void Test_S2BooleanOperation_GetCrossedVertexIndexBug6()
        {
            // This is another test of the code in GetCrossedVertexIndex() that checks
            // whether the B subchain contains an interior vertex of the A edge.
            Loops a_loops = new()
            {
                new()
                {
                    new(0.99870488823558456, 0.026138065586168355, 0.043650289137205818),
                    new(0.99876259434149239, 0.030513215246694664, 0.0392711578586665),
                    new(0.99984769515639127, 0.017452406437283512, 0),
                    new(0.998782023517925, 0.034862286684437908, 0.034915476003791211),
                    new(0.99878202512991221, 0.034878236872062651, 0.034899496702500969),
                    new(0.9975640502598242, 0.069756473744125302, 0),
                    new(0.99877979583714305, 0.034883478425067296, 0.034958008531414335),
                    new(0.99619692339885657, 0.052208468483931986, 0.069756473744125302),
                    new(0.99847581234813876, 0.017465633646566288, 0.052354596713645812),
                    new(0.9975640502598242, 0, 0.069756473744125302),
                    new(0.99847674250410212, 0.017444393356200013, 0.052343937746706169),
                    new(0.99847743863945992, 0.017428488520812163, 0.052335956242943835),
                    new(0.99984769515639127, 0, 0.017452406437283512),
                },
                new()
                {
                    new(0.99619692339885657, 0.052208468483931986, 0.069756473744125302),
                    new(0.99802119661969568, 0.034851668280404598, 0.052335956242943835),
                    new(0.9987605225894034, 0.030527121154938986, 0.039313018084772409),
                    new(0.99870321796526884, 0.026161932439896601, 0.043674199670139441),
                },
                new()
                {
                    new(0.99619692339885657, 0.052208468483931986, 0.069756473744125302),
                    new(0.99619692339885657, 0.06966087492121549, 0.052335956242943835),
                    new(0.99513403437078507, 0.069586550480032719, 0.069756473744125302),
                },
            };
            Loops b_loops = new()
            {
                new()
                {
                    new(0.99802200429988497, 0.034828499898458924, 0.052335977377554299),
                    new(0.99862953475457383, 0, 0.052335956242943835),
                    new(0.99923793061512223, 0.017455729388178846, 0.034912111530741322),
                    new(0.99923859085845868, 0.017443155365764275, 0.034899496702500969),
                    new(0.99923793076147094, 0.017455737780810811, 0.034912103145779166),
                    new(0.9992865072388355, 0.020934110218524152, 0.0314362764933699),
                    new(1, 0, 0),
                    new(0.99929987808789411, 0.022418034384064717, 0.029953053064335624),
                    new(0.99931406232431441, 0.02616995393092059, 0.026201876881811362),
                    new(0.99984769515639127, 0.017452406437283512, 0),
                    new(0.99930573320200933, 0.029072747464899757, 0.023298646837028814),
                    new(0.99862953475457383, 0.052335956242943835, 1.700986599320836e-73),
                    new(0.99838518277004218, 0.038347188759395717, 0.041910857059723181),
                    new(0.99619692339885668, 0.052208468483931979, 0.069756473744125289),
                },
                new()
                {
                    new(0.99802119662406841, 0.052304074592470849, 0.034899496702500969),
                    new(0.99847743834686298, 0.052327990806397578, 0.017452406437283512),
                    new(0.99619645281505653, 0.052208443821680058, 0.069763212314351342),
                    new(0.99619692339885657, 0.052208468483932, 0.069756473744125316),
                    new(0.99619692339885657, 0.052208468483931986, 0.069756473744125302),
                    new(0.99619692339885679, 0.052208468483931993, 0.069756473744125316),
                    new(0.99619692339885679, 0.052208468483931986, 0.069756473744125302),
                    new(0.99619692339885668, 0.052208468483931979, 0.069756473744125289),
                },
            };
            S2LaxPolygonShape result = new();
            ComputeTestUnion(a_loops, b_loops, S1Angle.Zero, result);
        }

        // Tests Contains() on empty and full geometries.
        [Fact]
        public void Test_S2BooleanOperation_ContainsEmptyAndFull()
        {
            var empty = MakeIndexOrDie("# #");
            var full = MakeIndexOrDie("# # full");
            Assert.True(Contains(empty, empty));
            Assert.False(Contains(empty, full));
            Assert.True(Contains(full, empty));
            Assert.True(Contains(full, full));
        }

        // Tests Intersects() on empty and full geometries.
        [Fact]
        public void Test_S2BooleanOperation_IntersectsEmptyAndFull()
        {
            var empty = MakeIndexOrDie("# #");
            var full = MakeIndexOrDie("# # full");
            Assert.False(Intersects(empty, empty));
            Assert.False(Intersects(empty, full));
            Assert.False(Intersects(full, empty));
            Assert.True(Intersects(full, full));
        }

        private static void ExpectResult(OpType op_type,
                  Options options,
                  string a_str, string b_str,
                  string expected_str)
        {
            var a = MakeIndexOrDie(a_str);
            var b = MakeIndexOrDie(b_str);
            var expected = MakeIndexOrDie(expected_str);
            ExpectResult(op_type, options, a, b, expected);
        }

        private static void ExpectResult(OpType op_type,
            Options options, S2ShapeIndex a, S2ShapeIndex b,
            S2ShapeIndex expected)
        {
            List<Layer> layers = new();
            for (int dim = 0; dim < 3; ++dim)
            {
                // Since all S2Builder polygon layers require DISCARD or DISCARD_EXCESS
                // for degenerate edges, we intentionally do not require any specific
                // multiplicity for degenerate edges and sibling pairs of dimension 2.
                GraphOptions graph_options = new(
                    EdgeType.DIRECTED,
                    (dim == 2) ? DegenerateEdges.DISCARD_EXCESS : DegenerateEdges.KEEP,
                    DuplicateEdges.KEEP,
                    (dim == 2) ? SiblingPairs.DISCARD_EXCESS : SiblingPairs.KEEP);
                layers.Add(new IndexMatchingLayer(
                    graph_options, expected, dim));
            }
            S2BooleanOperation op = new(op_type, layers, options);
            Assert.True(op.Build(a, b, out _));
            // S2BooleanOperation.OpTypeToString(op_type) << " failed:\n"
            // << "Expected result: " << s2textformat.ToString(expected) << "\n"
            // << error;

            // Now try the same thing with boolean output.
            Assert.Equal(expected.NumShapeIds() == 0,
                      IsEmpty(op_type, a, b, options));
        }

        // The intersections in the "expected" data below were computed in lat-lng
        // space (i.e., the rectangular projection), while the actual intersections
        // are computed using geodesics.  We can compensate for this by rounding the
        // intersection points to a fixed precision in degrees (e.g., 2 decimals).
        private static Options RoundToE(int exp)
        {
            Options options = new();
            options.SnapFunction_ = new IntLatLngSnapFunction(exp);
            return options;
        }

        [Fact]
        public void Test_S2BooleanOperation_SelfIntersectingPolylines()
        {
            // Two polylines that intersect at the point 2:4, and that also have
            // self-intersections at the points 2:2 and 3:4 respectively.  The
            // intersection point should always be created, but the self-intersection
            // points should be created iff split_all_crossing_polyline_edges() is true.

            Options options = RoundToE(1);
            var a = "# 0:2, 4:2, 2:0, 2:5 #";
            var b = "# 0:4, 5:4, 3:6, 3:3 #";
            ExpectResult(OpType.UNION, options, a, b,
                "# 0:2, 4:2, 2:0, 2:4, 2:5 | 0:4, 2:4, 5:4, 3:6, 3:3 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# 2:4, 2:4 | 2:4, 2:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# 0:2, 4:2, 2:0, 2:4, 2:5 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# 0:2, 4:2, 2:0, 2:4, 2:5 | 0:4, 2:4, 5:4, 3:6, 3:3 #");

            options.split_all_crossing_polyline_edges_ = true;
            ExpectResult(OpType.UNION, options, a, b,
                "# 0:2, 2:2, 4:2, 2:0, 2:2, 2:4, 2:5 " +
                "| 0:4, 2:4, 3:4, 5:4, 3:6, 3:4, 3:3 #");
            ExpectResult(OpType.INTERSECTION, options, a, b,
                "# 2:4, 2:4 | 2:4, 2:4 #");
            ExpectResult(OpType.DIFFERENCE, options, a, b,
                "# 0:2, 2:2, 4:2, 2:0, 2:2, 2:4, 2:5 #");
            ExpectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
                "# 0:2, 2:2, 4:2, 2:0, 2:2, 2:4, 2:5 " +
                "| 0:4, 2:4, 3:4, 5:4, 3:6, 3:4, 3:3 #");
        }

        // Subtracts a degenerate loop along the 180 degree meridian from the given
        // input geometry, and compares the result to "expected_str".  The inputs should
        // be in the format expected by S2TextFormat.MakeIndex().
        private static void TestMeridianSplitting(string input_str, string expected_str)
        {
            var input = MakeIndexOrDie(input_str);
            MutableS2ShapeIndex meridian = new();
            var loops = new Loops {
        new List<S2Point>{
            new S2Point(0, 0, -1), new S2Point(-1, 0, 0),
                    new S2Point(0, 0, 1), new S2Point(-1, 0, 0)}
        };
            meridian.Add(new S2LaxPolygonShape(loops));
            MutableS2ShapeIndex output = new();
            var layers = new Layer[3];
            layers[0] = new IndexedS2PointVectorLayer(output);
            // TODO(ericv): Implement s2builderutil.IndexedS2LaxPolylineVectorLayer.
            layers[1] = new IndexedS2PolylineVectorLayer(output);
            layers[2] = new IndexedLaxPolygonLayer(output);
            S2BooleanOperation op = new(OpType.DIFFERENCE, layers.ToList());
            Assert.True(op.Build(input, meridian, out _));
            Assert.Equal(expected_str, output.ToDebugString());
        }

        // Performs the given operation and compares the result to "expected_str".  All
        // arguments are in S2TextFormat.MakeLaxPolygon() format.
        private static void ExpectPolygon(OpType op_type,
            string a_str, string b_str, string expected_str)
        {
            var a = MakeIndexOrDie("# # " + a_str);
            var b = MakeIndexOrDie("# # " + b_str);
            LaxPolygonLayer.Options polygon_options = new();
            polygon_options.DegenerateBoundaries_ = (LaxPolygonLayer.Options.DegenerateBoundaries.DISCARD);
            S2LaxPolygonShape output = new();
            S2BooleanOperation op = new(op_type,
              new LaxPolygonLayer(output, polygon_options),
              new Options(
                  new IdentitySnapFunction(
                S1Angle.FromDegrees(1.1))
              ));
            Assert.True(op.Build(a, b, out _));
            Assert.Equal(expected_str, output.ToString());
        }

        // Tests whether the two S2ShapeIndexes are equal according to
        // S2BooleanOperation.Equals().
        private static bool TestEqual(string a_str, string b_str)
        {
            var a = MakeIndexOrDie(a_str);
            var b = MakeIndexOrDie(b_str);
            return S2BooleanOperation.Equals(a, b);
        }
    }
}
