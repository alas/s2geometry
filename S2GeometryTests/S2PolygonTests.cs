using Xunit;
using Xunit.Abstractions;
using S2Geometry.S2BuilderUtil;

namespace S2Geometry
{
    public class S2PolygonTests
    {
        #region Fields, Constants

        // A set of nested loops around the point 0:0 (lat:lng).
        // Every vertex of kNear0 is a vertex of kNear1.
        private const string kNearPoint = "0:0";
        private const string kNear0 = "-1:0, 0:1, 1:0, 0:-1;";
        private const string kNear1 = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;";
        private const string kNear2 = "-1:-2, -2:5, 5:-2;";
        private const string kNear3 = "-2:-2, -3:6, 6:-3;";
        private const string kNearHemi = "0:-90, -90:0, 0:90, 90:0;";

        // A set of nested loops around the point 0:180 (lat:lng).
        // Every vertex of kFar0 and kFar2 belongs to kFar1, and all
        // the loops except kFar2 are non-convex.
        private const string kFar0 = "0:179, 1:180, 0:-179, 2:-180;";
        private const string kFar1 = "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;";
        private const string kFar2 = "3:-178, 3:178, -1:179, -1:-179;";
        private const string kFar3 = "-3:-178, 4:-177, 4:177, -3:178, -2:179;";
        private const string kFarHemi = "0:-90, 60:90, -60:90;";

        // A set of nested loops around the point -90:0 (lat:lng).
        private const string kSouthPoint = "-89.9999:0.001";
        private const string kSouth0a = "-90:0, -89.99:0.01, -89.99:0;";
        private const string kSouth0b = "-90:0, -89.99:0.03, -89.99:0.02;";
        private const string kSouth0c = "-90:0, -89.99:0.05, -89.99:0.04;";
        private const string kSouth1 = "-90:0, -89.9:0.1, -89.9:-0.1;";
        private const string kSouth2 = "-90:0, -89.8:0.2, -89.8:-0.2;";
        private const string kSouthHemi = "0:-180, 0:60, 0:-60;";

        // Two different loops that surround all the Near and Far loops except
        // for the hemispheres.
        private const string kNearFar1 = "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, " +
                         "1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;";
        private const string kNearFar2 = "-2:15, -2:170, -8:-175, 8:-175, " +
                         "2:170, 2:15, 8:-4, -8:-4;";

        // Loops that result from intersection of other loops.
        private const string kFarHSouthH = "0:-180, 0:90, -60:90, 0:-90;";

        // Rectangles that form a cross, with only shared vertices, no crossing edges.
        // Optional holes outside the intersecting region.
        private const string kCross1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1;";
        private const string kCross1SideHole = "-1.5:0.5, -1.2:0.5, -1.2:-0.5, -1.5:-0.5;";
        private const string kCross2 = "1:-2, 1:-1, 1:1, 1:2, -1:2, -1:1, -1:-1, -1:-2;";
        private const string kCross2SideHole = "0.5:-1.5, 0.5:-1.2, -0.5:-1.2, -0.5:-1.5;";
        private const string kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

        // Two rectangles that intersect, but no edges cross and there's always
        // local containment (rather than crossing) at each shared vertex.
        // In this ugly ASCII art, 1 is A+B, 2 is B+C:
        //      +---+---+---+
        //      | A | B | C |
        //      +---+---+---+
        private const string kOverlap1 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
        private const string kOverlap1SideHole = "0.2:0.8, 0.8:0.8, 0.8:0.2, 0.2:0.2;";
        private const string kOverlap2 = "1:1, 2:1, 3:1, 3:0, 2:0, 1:0;";
        private const string kOverlap2SideHole = "2.2:0.8, 2.8:0.8, 2.8:0.2, 2.2:0.2;";
        private const string kOverlapCenterHole = "1.2:0.8, 1.8:0.8, 1.8:0.2, 1.2:0.2;";

        // An empty polygon.
        private const string kEmpty = "";
        // By symmetry, the intersection of the two polygons has almost half the area
        // of either polygon.
        private const string kOverlap3 = "-10:10, 0:10, 0:-10, -10:-10, -10:0";
        private const string kOverlap4 = "-10:0, 10:0, 10:-10, -10:-10";

        private const int kIters = 100;

        private readonly ITestOutputHelper _logger;

        private readonly TestCase[] test_cases = {
            // Two triangles that share an edge.
            new("4:2, 3:1, 3:3;",
                "3:1, 2:2, 3:3;",
                "",  // and
                "4:2, 3:1, 2:2, 3:3;",  // or
                "4:2, 3:1, 3:3;",  // minus
                "4:2, 3:1, 2:2, 3:3;"),  // xor  
            // Two vertical bars and a horizontal bar connecting them.
            new("0:0, 0:2, 3:2, 3:0;   0:3, 0:5, 3:5, 3:3;",
                "1:1, 1:4, 2:4, 2:1;",
                "1:1, 1:2, 2:2, 2:1;   1:3, 1:4, 2:4, 2:3;",  // and
                "0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0;",  // or
                "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;   " + // minus 
                "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;",
                "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;   " + // xor
                "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;   " +
                "1:2, 1:3, 2:3, 2:2"),
            // Two vertical bars and two horizontal bars.
            new("1:88, 1:93, 2:93, 2:88;   -1:88, -1:93, 0:93, 0:88;",
                "-2:89, -2:90, 3:90, 3:89;   -2:91, -2:92, 3:92, 3:91;",
                "1:89, 1:90, 2:90, 2:89;   1:91, 1:92, 2:92, 2:91;   "+  // and
                "-1:89, -1:90, 0:90, 0:89;   -1:91, -1:92, 0:92, 0:91;",
                "-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, "+  // or
                "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, "+
                "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88;   "+
                "0:90, 0:91, 1:91, 1:90;",
                "1:88, 1:89, 2:89, 2:88;   1:90, 1:91, 2:91, 2:90;   "  +// minus
                "1:92, 1:93, 2:93, 2:92;   -1:88, -1:89, 0:89, 0:88;   "+
                "-1:90, -1:91, 0:91, 0:90;   -1:92, -1:93, 0:93, 0:92;",
                "1:88, 1:89, 2:89, 2:88;   -1:88, -1:89, 0:89, 0:88;   " + // xor
                "1:90, 1:91, 2:91, 2:90;   -1:90, -1:91, 0:91, 0:90;   "+
                "1:92, 1:93, 2:93, 2:92;   -1:92, -1:93, 0:93, 0:92;   "+
                "-2:89, -2:90, -1:90, -1:89;   -2:91, -2:92, -1:92, -1:91;   "+
                "0:89, 0:90, 1:90, 1:89;   0:91, 0:92, 1:92, 1:91;   "+
                "2:89, 2:90, 3:90, 3:89;   2:91, 2:92, 3:92, 3:91;"),
            // Two interlocking square doughnuts.
            new ("-1:-93, -1:-89, 3:-89, 3:-93;   0:-92, 0:-90, 2:-90, 2:-92;",
                "-3:-91, -3:-87, 1:-87, 1:-91;   -2:-90, -2:-88, 0:-88, 0:-90;",
                "-1:-91, -1:-90, 0:-90, 0:-91;   0:-90, 0:-89, 1:-89, 1:-90;",  // and
                "-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93;   "+  // or
                "0:-92, 0:-91, 1:-91, 1:-90, 2:-90, 2:-92;   "+
                "-2:-90, -2:-88, 0:-88, 0:-89, -1:-89, -1:-90;",
                "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, " + // minus
                "1:-90, 1:-89, 3:-89, 3:-93;   "+
                "-1:-90, -1:-89, 0:-89, 0:-90;",
                "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, " + // xor
                "1:-90, 1:-89, 3:-89, 3:-93;   "+
                "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88, "+
                "-2:-88, -2:-90, -1:-90, -1:-91;   "+
                "-1:-90, -1:-89, 0:-89, 0:-90;   "+
                "1:-91, 0:-91, 0:-90, 1:-90;"),
            // An incredibly thin triangle intersecting a square, such that the two
            // intersection points of the triangle with the square are identical.
            // This results in a degenerate loop that needs to be handled correctly.
            new("10:44, 10:46, 12:46, 12:44;",
                "11:45, 89:45.00000000000001, 90:45;",
                "",  // Empty intersection!
                // Original square with extra vertex, and triangle disappears (due to
                // default vertex_merge_radius of S2EdgeCrossings.kIntersectionMergeRadius).
                "10:44, 10:46, 12:46, 12:45.001774937, 12:44;",  // or
                "10:44, 10:46, 12:46, 12:45.001774937, 12:44;",  // minus
                "10:44, 10:46, 12:46, 12:45.001774937, 12:44;"),  // xor
        };

        // Some standard polygons to use in the tests.
        private readonly S2Polygon empty_;
        private readonly S2Polygon full_;
        private readonly S2Polygon near_0_;
        private readonly S2Polygon near_10_;
        private readonly S2Polygon near_30_;
        private readonly S2Polygon near_32_;
        private readonly S2Polygon near_3210_;
        private readonly S2Polygon near_H3210_;

        private readonly S2Polygon far_10_;
        private readonly S2Polygon far_21_;
        private readonly S2Polygon far_321_;
        private readonly S2Polygon far_H20_;
        private readonly S2Polygon far_H3210_;

        private readonly S2Polygon south_0ab_;
        private readonly S2Polygon south_2_;
        private readonly S2Polygon south_210b_;
        private readonly S2Polygon south_H21_;
        private readonly S2Polygon south_H20abc_;

        private readonly S2Polygon nf1_n10_f2_s10abc_;

        private readonly S2Polygon nf2_n2_f210_s210ab_;

        private readonly S2Polygon f32_n0_;
        private readonly S2Polygon n32_s0b_;

        private readonly S2Polygon cross1_;
        private readonly S2Polygon cross1_side_hole_;
        private readonly S2Polygon cross1_center_hole_;
        private readonly S2Polygon cross2_;
        private readonly S2Polygon cross2_side_hole_;
        private readonly S2Polygon cross2_center_hole_;

        private readonly S2Polygon overlap1_;
        private readonly S2Polygon overlap1_side_hole_;
        private readonly S2Polygon overlap1_center_hole_;
        private readonly S2Polygon overlap2_;
        private readonly S2Polygon overlap2_side_hole_;
        private readonly S2Polygon overlap2_center_hole_;

        private readonly S2Polygon far_H_;
        private readonly S2Polygon south_H_;
        private readonly S2Polygon far_H_south_H_;

        #region SetInput

        private S2Polygon simplified;
        private S2Polygon original;

        #endregion

        #region IsValidTest

        private delegate void ModifyPolygonDelegate(S2Polygon polygon);
        private ModifyPolygonDelegate modify_polygon_hook_;

        private bool init_oriented_;
        private readonly List<List<S2Point>> vloops_ = new();

        #endregion

        #if DEBUG

        private readonly S2PolygonDecodeTest decodeTest;
        
        #endif

        #endregion

        public S2PolygonTests(ITestOutputHelper logger)
        {
            _logger = logger;
            empty_ = new S2Polygon();
            full_ = MakePolygon("full");
            near_0_ = MakePolygon(kNear0);
            near_10_ = MakePolygon(kNear0 + kNear1);
            near_30_ = MakePolygon(kNear3 + kNear0);
            near_32_ = MakePolygon(kNear2 + kNear3);
            near_3210_ = MakePolygon(kNear0 + kNear2 + kNear3 + kNear1);
            near_H3210_ = MakePolygon(kNear0 + kNear2 + kNear3 + kNearHemi + kNear1);

            far_10_ = MakePolygon(kFar0 + kFar1);
            far_21_ = MakePolygon(kFar2 + kFar1);
            far_321_ = MakePolygon(kFar2 + kFar3 + kFar1);
            far_H20_ = MakePolygon(kFar2 + kFarHemi + kFar0);
            far_H3210_ = MakePolygon(kFar2 + kFarHemi + kFar0 + kFar1 + kFar3);

            south_0ab_ = MakePolygon(kSouth0a + kSouth0b);
            south_2_ = MakePolygon(kSouth2);
            south_210b_ = MakePolygon(kSouth2 + kSouth0b + kSouth1);
            south_H21_ = MakePolygon(kSouth2 + kSouthHemi + kSouth1);
            south_H20abc_ = MakePolygon(kSouth2 + kSouth0b + kSouthHemi + kSouth0a + kSouth0c);

            nf1_n10_f2_s10abc_ = MakePolygon(kSouth0c + kFar2 + kNear1 + kNearFar1 +
                                           kNear0 + kSouth1 + kSouth0b + kSouth0a);

            nf2_n2_f210_s210ab_ = MakePolygon(kFar2 + kSouth0a + kFar1 + kSouth1 + kFar0 +
                                            kSouth0b + kNearFar2 + kSouth2 + kNear2);

            f32_n0_ = MakePolygon(kFar2 + kNear0 + kFar3);
            n32_s0b_ = MakePolygon(kNear3 + kSouth0b + kNear2);

            cross1_ = MakePolygon(kCross1);
            cross1_side_hole_ = MakePolygon(kCross1 + kCross1SideHole);
            cross1_center_hole_ = MakePolygon(kCross1 + kCrossCenterHole);
            cross2_ = MakePolygon(kCross2);
            cross2_side_hole_ = MakePolygon(kCross2 + kCross2SideHole);
            cross2_center_hole_ = MakePolygon(kCross2 + kCrossCenterHole);

            overlap1_ = MakePolygon(kOverlap1);
            overlap1_side_hole_ = MakePolygon(kOverlap1 + kOverlap1SideHole);
            overlap1_center_hole_ = MakePolygon(kOverlap1 + kOverlapCenterHole);
            overlap2_ = MakePolygon(kOverlap2);
            overlap2_side_hole_ = MakePolygon(kOverlap2 + kOverlap2SideHole);
            overlap2_center_hole_ = MakePolygon(kOverlap2 + kOverlapCenterHole);

            far_H_ = MakePolygon(kFarHemi);
            south_H_ = MakePolygon(kSouthHemi);
            far_H_south_H_ = MakePolygon(kFarHSouthH);
                        
            #region IsValidTest
            
            init_oriented_ = false;
            modify_polygon_hook_ = null;
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);

            #endregion

            #if DEBUG

            decodeTest = new S2PolygonDecodeTest(); 
            
            #endif
        }

        [Fact]
        public void Test_S2Polygon_Init() {
            CheckContains(kNear1, kNear0);
            CheckContains(kNear2, kNear1);
            CheckContains(kNear3, kNear2);
            CheckContains(kNearHemi, kNear3);
            CheckContains(kFar1, kFar0);
            CheckContains(kFar2, kFar1);
            CheckContains(kFar3, kFar2);
            CheckContains(kFarHemi, kFar3);
            CheckContains(kSouth1, kSouth0a);
            CheckContains(kSouth1, kSouth0b);
            CheckContains(kSouth1, kSouth0c);
            CheckContains(kSouthHemi, kSouth2);
            CheckContains(kNearFar1, kNear3);
            CheckContains(kNearFar1, kFar3);
            CheckContains(kNearFar2, kNear3);
            CheckContains(kNearFar2, kFar3);

            CheckContainsPoint(kNear0, kNearPoint);
            CheckContainsPoint(kNear1, kNearPoint);
            CheckContainsPoint(kNear2, kNearPoint);
            CheckContainsPoint(kNear3, kNearPoint);
            CheckContainsPoint(kNearHemi, kNearPoint);
            CheckContainsPoint(kSouth0a, kSouthPoint);
            CheckContainsPoint(kSouth1, kSouthPoint);
            CheckContainsPoint(kSouth2, kSouthPoint);
            CheckContainsPoint(kSouthHemi, kSouthPoint);
        }

        [Fact]
        public void Test_S2Polygon_OverlapFractions() {
            S2Polygon a = MakePolygon(kEmpty);
            S2Polygon b = MakePolygon(kEmpty);
            var result = S2Polygon.GetOverlapFractions(a, b);
            Assert2.Near(1.0, result.Item1);
            Assert2.Near(1.0, result.Item2);

            b = MakePolygon(kOverlap3);
            result = S2Polygon.GetOverlapFractions(a, b);
            Assert2.Near(1.0, result.Item1);
            Assert2.Near(0.0, result.Item2);

            a = MakePolygon(kOverlap4);
            result = S2Polygon.GetOverlapFractions(a, b);
            Assert2.Near(0.5, result.Item1, 1e-14);
            Assert2.Near(0.5, result.Item2, 1e-14);
        }

        [Fact]
        public void Test_S2Polygon_OriginNearPole() {
            // S2Polygon operations are more efficient if S2.Origin is near a pole.
            // (Loops that contain a pole tend to have very loose bounding boxes because
            // they span the full longitude range.  S2Polygon canonicalizes all loops so
            // that they don't contain S2.Origin, thus by placing S2.Origin near a
            // pole we minimize the number of canonical loops which contain that pole.)
            Assert.True(S2LatLng.Latitude(S2.Origin).GetDegrees() >= 80);
        }

        [Fact]
        public void Test_S2Polygon_TestApproxContainsAndDisjoint() {
            // We repeatedly choose a random cell id and intersect its bounding polygon
            // "A" with the bounding polygon "B" of one its child cells.  The result may
            // not be contained by either A or B, because the vertices of B near the
            // edge midpoints of A may be slightly outside A, and even when the crossing
            // edges are intersected, the intersection point may also be slightly
            // outside A and/or B.
            //
            // We repeat the test many times and expect that some fraction of the exact
            // tests should fail, while all of the approximate test should succeed.
            int kIters = 1000;
            int exact_contains = 0, exact_disjoint = 0;
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            for (int iter = 0; iter < kIters; ++iter) {
                S2CellId id = S2Testing.GetRandomCellId(10);
                S2Polygon parent_polygon = new(new S2Cell(id));
                S2Polygon child_polygon = new(new S2Cell(id.Child(0)));

                // Get the intersection.  There is no guarantee that the intersection will
                // be contained by A or B.  Similarly, the intersection may slightly
                // overlap an adjacent disjoint polygon C.
                S2Polygon intersection = new();
                intersection.InitToIntersection(parent_polygon, child_polygon);
                if (parent_polygon.Contains(intersection)) {
                    ++exact_contains;
                }
                Assert.True(parent_polygon.ApproxContains(
                    intersection, S2.kIntersectionMergeRadiusS1Angle));

                S2Polygon adjacent_polygon = new(new S2Cell(id.Child(1)));
                if (!adjacent_polygon.Intersects(intersection)) {
                    ++exact_disjoint;
                }
                Assert.True(adjacent_polygon.ApproxDisjoint(
                    intersection, S2.kIntersectionMergeRadiusS1Angle));
            }
            // All of the approximate results are true, so we check that at least some
            // of the exact results are false in order to make sure that this test
            // actually tests something.
            //
            // There are two vertices in each child cell that have a 50% chance of being
            // outside the parent cell.  When a vertex is outside, an intersection point
            // is computed near that vertex that also has a 50% chance of being
            // outside.  Snapping used to choose one of these vertices at random, but
            // currently the vertex whose S2CellId is smaller is always chosen.  For the
            // exact containment test, it turns out that one vertex is adjacent to a
            // lower-numbered S2CellId and the other is adjacent to a higher-numbered
            // S2CellId, which means that one vertex will always be chosen outside the
            // parent if possible, and the other will always be chosen inside if
            // possible.  This works out to an expectation that 0.5 * 0.75 = 37.5% of
            // the exact containment tests will succeed.
            //
            // For the exact disjoint test, there is one shared vertex that might be
            // replaced by a computed intersection point.  The shared vertex is inside
            // the parent 50% of the time.  Otherwise there is a 50% chance that the
            // intersection point will not be chosen for snapping because it has a
            // higher S2CellId that the shared vertex, and otherwise there is still a
            // 50-75% chance that intersection point will not be inside the adjacent
            // child cell (depending on how far the shared vertex is outside the parent
            // cell).  This means that we expect 1 - 0.5 * 0.5 * (0.25 ~ 0.5) = 87.5% to
            // 93.75% of the exact disjoint tests to succeed on average.
            Assert.True(exact_contains <= 0.40 * kIters);  // about 37.5% succeed
            Assert.True(exact_disjoint <= 0.96 * kIters);  // 87.5% - 93.75% succeed
        }

        [Fact]
        public void Test_S2PolygonTestBase_Relations() {
            TestRelation(near_10_, empty_, true, false, false);
            TestRelation(near_10_, near_10_, true, true, true);
            TestRelation(full_, near_10_, true, false, true);
            TestRelation(near_10_, near_30_, false, true, true);
            TestRelation(near_10_, near_32_, false, false, false);
            TestRelation(near_10_, near_3210_, false, true, true);
            TestRelation(near_10_, near_H3210_, false, false, false);
            TestRelation(near_30_, near_32_, true, false, true);
            TestRelation(near_30_, near_3210_, true, false, true);
            TestRelation(near_30_, near_H3210_, false, false, true);
            TestRelation(near_32_, near_3210_, false, true, true);
            TestRelation(near_32_, near_H3210_, false, false, false);
            TestRelation(near_3210_, near_H3210_, false, false, false);

            TestRelation(far_10_, far_21_, false, false, false);
            TestRelation(far_10_, far_321_, false, true, true);
            TestRelation(far_10_, far_H20_, false, false, false);
            TestRelation(far_10_, far_H3210_, false, false, false);
            TestRelation(far_21_, far_321_, false, false, false);
            TestRelation(far_21_, far_H20_, false, false, false);
            TestRelation(far_21_, far_H3210_, false, true, true);
            TestRelation(far_321_, far_H20_, false, false, true);
            TestRelation(far_321_, far_H3210_, false, false, true);
            TestRelation(far_H20_, far_H3210_, false, false, true);

            TestRelation(south_0ab_, south_2_, false, true, true);
            TestRelation(south_0ab_, south_210b_, false, false, true);
            TestRelation(south_0ab_, south_H21_, false, true, true);
            TestRelation(south_0ab_, south_H20abc_, false, true, true);
            TestRelation(south_2_, south_210b_, true, false, true);
            TestRelation(south_2_, south_H21_, false, false, true);
            TestRelation(south_2_, south_H20abc_, false, false, true);
            TestRelation(south_210b_, south_H21_, false, false, true);
            TestRelation(south_210b_, south_H20abc_, false, false, true);
            TestRelation(south_H21_, south_H20abc_, true, false, true);

            TestRelation(nf1_n10_f2_s10abc_, nf2_n2_f210_s210ab_, false, false, true);
            TestRelation(nf1_n10_f2_s10abc_, near_32_, true, false, true);
            TestRelation(nf1_n10_f2_s10abc_, far_21_, false, false, false);
            TestRelation(nf1_n10_f2_s10abc_, south_0ab_, false, false, false);
            TestRelation(nf1_n10_f2_s10abc_, f32_n0_, true, false, true);

            TestRelation(nf2_n2_f210_s210ab_, near_10_, false, false, false);
            TestRelation(nf2_n2_f210_s210ab_, far_10_, true, false, true);
            TestRelation(nf2_n2_f210_s210ab_, south_210b_, true, false, true);
            TestRelation(nf2_n2_f210_s210ab_, south_0ab_, true, false, true);
            TestRelation(nf2_n2_f210_s210ab_, n32_s0b_, true, false, true);

            TestRelation(cross1_, cross2_, false, false, true);
            TestRelation(cross1_side_hole_, cross2_, false, false, true);
            TestRelation(cross1_center_hole_, cross2_, false, false, true);
            TestRelation(cross1_, cross2_side_hole_, false, false, true);
            TestRelation(cross1_, cross2_center_hole_, false, false, true);
            TestRelation(cross1_side_hole_, cross2_side_hole_, false, false, true);
            TestRelation(cross1_center_hole_, cross2_side_hole_, false, false, true);
            TestRelation(cross1_side_hole_, cross2_center_hole_, false, false, true);
            TestRelation(cross1_center_hole_, cross2_center_hole_, false, false, true);

            // These cases_, when either polygon has a hole, test a different code path
            // from the other cases.
            TestRelation(overlap1_, overlap2_, false, false, true);
            TestRelation(overlap1_side_hole_, overlap2_, false, false, true);
            TestRelation(overlap1_center_hole_, overlap2_, false, false, true);
            TestRelation(overlap1_, overlap2_side_hole_, false, false, true);
            TestRelation(overlap1_, overlap2_center_hole_, false, false, true);
            TestRelation(overlap1_side_hole_, overlap2_side_hole_, false, false, true);
            TestRelation(overlap1_center_hole_, overlap2_side_hole_, false, false, true);
            TestRelation(overlap1_side_hole_, overlap2_center_hole_, false, false, true);
            TestRelation(overlap1_center_hole_, overlap2_center_hole_, false, false, true);

            void TestRelation(S2Polygon a, S2Polygon b, bool contains, bool contained, bool intersects)
                => TestRelationWithDesc(a, b, contains, contained, intersects, "args " + a + ", " + b);
        }

        [Fact]
        public void Test_S2PolygonTestBase_EmptyAndFull() {
            Assert.True(empty_.IsEmpty());
            Assert.False(full_.IsEmpty());
            Assert.False(empty_.IsFull());
            Assert.True(full_.IsFull());

            TestNestedPair(empty_, empty_);
            TestNestedPair(full_, empty_);
            TestNestedPair(full_, full_);
        }

        [Fact]
        public void Test_S2PolygonTestBase_Operations() {
            S2Polygon far_south = new();
            far_south.InitToIntersection(far_H_, south_H_);
            CheckEqual(far_south, far_H_south_H_, S1Angle.FromRadians(1e-15));

            int i = 0;
            foreach (var test in test_cases) {
                _logger.WriteLine($"Polygon operation test case {i++}");
                S2Polygon a = MakePolygon(test.A);
                S2Polygon b = MakePolygon(test.B);
                S2Polygon expected_a_and_b = MakePolygon(test.AAndB);
                S2Polygon expected_a_or_b = MakePolygon(test.AOrB);
                S2Polygon expected_a_minus_b = MakePolygon(test.AMinusB);
                S2Polygon expected_a_xor_b = MakePolygon(test.AXorB);

                // The intersections in the "expected" data were computed in lat-lng
                // space, while the actual intersections are computed using geodesics.
                // The error due to this depends on the length and direction of the line
                // segment being intersected, and how close the intersection is to the
                // endpoints of the segment.  The worst case is for a line segment between
                // two points at the same latitude, where the intersection point is in the
                // middle of the segment.  In this case the error is approximately
                // (p * t^2) / 8, where "p" is the absolute latitude in radians, "t" is
                // the longitude difference in radians, and both "p" and "t" are small.
                // The test cases all have small latitude and longitude differences.
                // If "p" and "t" are converted to degrees, the following error bound is
                // valid as long as (p * t^2 < 150).

                var kMaxError = S1Angle.FromRadians(1e-4);

                var a_and_b = new S2Polygon();
                var a_or_b = new S2Polygon();
                var a_minus_b = new S2Polygon();
                var a_xor_b = new S2Polygon();
                a_and_b.InitToIntersection(a, b);
                CheckEqual(a_and_b, expected_a_and_b, kMaxError);
                a_or_b.InitToUnion(a, b);
                CheckEqual(a_or_b, expected_a_or_b, kMaxError);
                TestDestructiveUnion(a, b);
                a_minus_b.InitToDifference(a, b);
                CheckEqual(a_minus_b, expected_a_minus_b, kMaxError);
                a_xor_b.InitToSymmetricDifference(a, b);
                CheckEqual(a_xor_b, expected_a_xor_b, kMaxError);
            }
        }

        [Fact]
        public void Test_S2Polygon_IntersectionSnapFunction() {
            // This tests that an intersection point is rounded to the nearest allowable
            // vertex position (using E0 coordinates, i.e. integer lat/lng values).
            var a = MakePolygon("0:0, 0:10, 1:10, 1:0");
            var b = MakePolygon("0:0, 0:10, 3:0");
            var expected = MakePolygon("0:0, 0:10, 1:7, 1:0");
            S2Polygon actual = new();
            actual.InitToIntersection(a, b, new IntLatLngSnapFunction(0));  // E0 coords
            CheckEqual(expected, actual);
        }

        [Fact]
        public void Test_S2Polygon_IntersectionPreservesLoopOrder() {
            var a = MakePolygon("0:0, 0:10, 10:10, 10:0");
            var b = MakePolygon("1:1, 1:9, 9:5; 2:2, 2:8, 8:5");
            S2Polygon actual = new();
            actual.InitToIntersection(a, b);
            Assert.Equal(b.ToDebugString(), actual.ToDebugString());
        }

        // Verifies that the bounding rectangle optimization in InitToIntersection()
        // resets the result polygon to be empty.
        [Fact]
        public void Test_S2Polygon_EmptyIntersectionClearsResult()
        {
            // The bounding rectangles of these two polygons do not intersect.
            var a = MakePolygon("0:0, 0:1, 1:0");
            var b = MakePolygon("3:3, 3:4, 4:3");

            // Initialize the result polygon to be non-empty, then verify that computing
            // the intersection clears the result.
            var result = MakePolygon("0:0, 0:1, 1:0");
            result.InitToIntersection(a, b);
            Assert.True(result.IsEmpty());

            // Repeat with the version of InitToIntersection that allows error reporting.
            S2Error error;
            result = MakePolygon("0:0, 0:1, 1:0");
            Assert.True(result.InitToIntersection(a, b,
                new IdentitySnapFunction(S1Angle.Zero), out error));
            Assert.True(result.IsEmpty());
        }

        // Verifies that S2Polygon does not destroy or replace pointers to S2Loop, so
        // caller can rely on using raw pointers.
        [Fact]
        public void Test_S2Polygon_LoopPointers() {
            var loops = new List<S2Loop>
            {
                MakeLoopOrDie("4:4, 4:6, 6:6, 6:4"),
                MakeLoopOrDie("3:3, 3:7, 7:7, 7:3"),
                MakeLoopOrDie("2:2, 2:8, 8:8, 8:2"),
                MakeLoopOrDie("1:1, 1:9, 9:9, 9:1"),
                MakeLoopOrDie("10:10, 15:15, 20:10"),
                MakeLoopOrDie("-1:-1, -9:-1, -9:-9, -1:-9"),
                MakeLoopOrDie("-5:-5, -6:-5, -6:-6, -5:-6")
            };

            var loops_raw_ptrs = new List<S2Loop>();
            foreach (var loop in loops) {
                loops_raw_ptrs.Add(loop);
            }
            S2Polygon polygon = new(loops);

            // Check that loop pointers didn't change (but could've gotten reordered).
            Assert.Equal(loops_raw_ptrs.Count, polygon.NumLoops());
            for (int i = 0; i < polygon.NumLoops(); i++) {
                Assert.Equal(1, loops_raw_ptrs.Count(t => t == polygon.Loop(i)));
            }
        }

        // The "Bug" tests are regression tests from previous versions of the algorithm.
        [Fact]
        public void Test_S2Polygon_Bug1() {
            var a_vertices = new S2Point[][]{
    new []{
      new S2Point(-0.10531193335759943, -0.80522214810955617, 0.58354664670985534),
      new S2Point(-0.10531194840431297, -0.80522215192439039, 0.58354663873039425),
      new S2Point(-0.10531192794033867, -0.80522217497559767, 0.58354661061568747),
      new S2Point(-0.10531191284235047, -0.80522217121852058, 0.58354661852470402)
    },
  };
            var b_vertices = new S2Point[][]{
     new S2Point[]{
      new S2Point(-0.10531174240075937, -0.80522236320875284, 0.58354638436119843),
      new S2Point(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
      new S2Point(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508),  // A
      new S2Point(-0.10531191288915481, -0.80522217116640804, 0.5835466185881667),   // B
      new S2Point(-0.10531191288915592, -0.8052221711664066, 0.58354661858816803),   // B
      new S2Point(-0.10531192039151964, -0.80522217309710431, 0.58354661457010204),  // A
      new S2Point(-0.10531192794033779, -0.80522217497559878, 0.58354661061568636),
      new S2Point(-0.1053117575499668, -0.80522236690813498, 0.58354637652254981),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            // Given edges do not form loops (indegree != outdegree)
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug2() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10618951389689163, -0.80546461394606728, 0.58305277875939732),
      new S2Point(-0.10618904764039243, -0.8054645437464607, 0.58305296065497536),
      new S2Point(-0.10618862643748632, -0.80546451917975415, 0.58305307130470341),
      new S2Point(-0.10617606798507535, -0.80544758470051458, 0.58307875187433833),
    },
  };
            S2Point[][] b_vertices = {
    new S2Point[]{
      new S2Point(-0.10618668131028208, -0.80544613076731553, 0.58307882755616247),
      new S2Point(-0.10618910658843225, -0.80546454998744921, 0.58305294129732887),
      new S2Point(-0.10618904764039225, -0.80546454374646081, 0.58305296065497536),
      new S2Point(-0.10618898834264634, -0.80546453817003949, 0.58305297915823251),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Given edges do not form loops (indegree != outdegree)
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug3() {
           var  a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10703494861068318, -0.80542232562508131, 0.58295659972299307),
      new S2Point(-0.10703494998722708, -0.80542232255642865, 0.58295660370995028),
      new S2Point(-0.10703495367938694, -0.80542232008675829, 0.58295660644418046),
      new S2Point(-0.10703495869785147, -0.80542231887781635, 0.58295660719304865),
      new S2Point(-0.10703496369792719, -0.80542231925353791, 0.58295660575589636),
      new S2Point(-0.10703496733984781, -0.80542232111324863, 0.58295660251780734),
      new S2Point(-0.10703496864776367, -0.80542232395864055, 0.58295659834642488),
      new S2Point(-0.10703496727121976, -0.80542232702729322, 0.58295659435946767),
      new S2Point(-0.10703496357905991, -0.80542232949696357, 0.5829565916252375),
      new S2Point(-0.10703495856059538, -0.80542233070590552, 0.58295659087636931),
      new S2Point(-0.10703495356051966, -0.80542233033018396, 0.58295659231352159),
      new S2Point(- 0.10703494991859903, -0.80542232847047324, 0.58295659555161061),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10703494861068762, -0.80542232562508098, 0.58295659972299274),
      new S2Point(-0.10703494998723152, -0.80542232255642832, 0.58295660370994995),
      new S2Point(-0.10703495367939138, -0.80542232008675796, 0.58295660644418013),
      new S2Point(-0.10703495869785591, -0.80542231887781601, 0.58295660719304832),
      new S2Point(-0.10703496369793163, -0.80542231925353758, 0.58295660575589603),
      new S2Point(-0.10703496733985225, -0.8054223211132483, 0.58295660251780701),
      new S2Point(-0.10703496864776811, -0.80542232395864022, 0.58295659834642455),
      new S2Point(-0.1070349672712242, -0.80542232702729288, 0.58295659435946734),
      new S2Point(-0.10703496357906438, -0.80542232949696346, 0.58295659162523727),
      new S2Point(-0.10703495856059982, -0.80542233070590519, 0.58295659087636897),
      new S2Point(-0.1070349535605241, -0.80542233033018362, 0.58295659231352126),
      new S2Point(-0.10703494991860348, -0.8054223284704729, 0.58295659555161028),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Given edges do not form loops (indegree != outdegree)
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug4() {
            var a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10667065556339718, -0.80657502337947207, 0.58142764201754193),
      new S2Point(-0.10667064691895933, -0.80657502457251051, 0.58142764194845853),
      new S2Point(-0.10667064691930939, -0.80657502457246333, 0.58142764194845975),
      new S2Point(-0.10667065556339746, -0.80657502337947395, 0.5814276420175396),
      new S2Point(-0.10667077559567185, -0.80657589269604968, 0.58142641405029793),
      new S2Point(-0.10667077059539463, -0.80657589232162286, 0.58142641548708696),
      new S2Point(-0.10667063827452879, -0.80657502576554818, 0.58142764187937435),
      new S2Point(-0.10667063169531328, -0.80657498170361974, 0.58142770421053058),
      new S2Point(- 0.10667064898418178, -0.8065749793175444, 0.58142770434869739),
    },
    new S2Point[]{
      new S2Point(-0.10667064691897719, -0.80657502457250896, 0.58142764194845697),
      new S2Point(-0.10667063827452879, -0.80657502576554818, 0.58142764187937435),
      new S2Point(-0.10667064691861985, -0.80657502457255736, 0.58142764194845586),
    },
  };
            var b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10667064691896312, -0.80657502457251107, 0.58142764194845697),
      new S2Point(-0.10667064691896297, -0.80657502457251007, 0.58142764194845853),
      new S2Point(-0.10667064033974753, -0.80657498051058207, 0.58142770427961399),
      new S2Point(-0.10667064076268165, -0.80657498045444342, 0.58142770427989865),
      new S2Point(-0.10667051785242875, -0.80657409963649807, 0.58142894872603923),
      new S2Point(-0.1066707756642685, -0.80657588679775971, 0.58142642222003538),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Loop 1: Edge 1 crosses edge 3
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug5() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10574444273627338, -0.80816264611829447, 0.57938868667714882),
      new S2Point(-0.10574444845633162, -0.80816268110163325, 0.57938863683652475),
      new S2Point(-0.10574444825833453, -0.80816268112970524, 0.57938863683350494),
      new S2Point(-0.10574444253827629, -0.80816264614636646, 0.57938868667412902),
      new S2Point(-0.10574408792844124, -0.80816047738475361, 0.57939177648757634),
      new S2Point(-0.10574408812643833, -0.80816047735668162, 0.57939177649059592),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.1057440881264381, -0.80816047735668017, 0.57939177649059825),
      new S2Point(-0.10574408802743954, -0.80816047737071606, 0.57939177648908835),
      new S2Point(-0.10574408812649677, -0.8081604773570521, 0.57939177649006868),
      new S2Point(-0.10574408812649701, -0.80816047735705354, 0.57939177649006646),
      new S2Point(-0.10574408802703171, -0.80816047737077379, 0.57939177648908202),
      new S2Point(-0.10574408792844098, -0.80816047738475194, 0.57939177648757834),
      new S2Point(-0.10574408792838257, -0.80816047738438168, 0.5793917764881058),
      new S2Point(-0.1057440879283823, -0.80816047738438002, 0.57939177648810791),
      new S2Point(-0.10574407993470979, -0.80816042849578984, 0.57939184613891748),
      new S2Point(-0.10574408013270691, -0.80816042846771807, 0.57939184614193739),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Loop 0 edge 8 crosses loop 1 edge 0
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug6() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10618849949725141, -0.80552159562437586, 0.58297423747304822),
      new S2Point(-0.10618849959636036, -0.80552159561106063, 0.58297423747339361),
      new S2Point(-0.10618849949722192, -0.80552159562415893, 0.5829742374733532),
      new S2Point(-0.10618834540082922, -0.80552043435619214, 0.58297587011440333),
      new S2Point(-0.10618834559910612, -0.80552043432999554, 0.58297587011448437),
      new S2Point(-0.10618849969546933, -0.80552159559774539, 0.58297423747373922),
      new S2Point(-0.10618849969546955, -0.80552159559774716, 0.582974237473737),
      new S2Point(-0.10618849969549882, -0.80552159559796233, 0.58297423747343424),
      new S2Point(-0.10618849959710704, -0.80552159561096182, 0.58297423747339394),
      new S2Point(-0.10618849949725161, -0.80552159562437742, 0.58297423747304589),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10618856154870562, -0.80552206324314812, 0.58297358004005528),
      new S2Point(-0.10618849949722212, -0.80552159562416048, 0.58297423747335086),
      new S2Point(-0.10618849969549901, -0.80552159559796388, 0.58297423747343191),
      new S2Point(-0.10618856174698249, -0.8055220632169513, 0.58297358004013622),
      new S2Point(-0.10618857104277038, -0.80552213326985989, 0.58297348155149287),
      new S2Point(-0.10618857084449349, -0.80552213329605649, 0.58297348155141182),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new ();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Loop 0 edge 0 crosses loop 1 edge 4
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug7() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10651728339354898, -0.80806023027835039, 0.57938996589599123),
      new S2Point(-0.10651728368541774, -0.80806023024121265, 0.57938996589412783),
      new S2Point(-0.10651743884289547, -0.80806147782022508, 0.5793881973990701),
      new S2Point(-0.1065172793067945, -0.80806153133252501, 0.5793881520963412),
      new S2Point(-0.10651707335497011, -0.80806158532388361, 0.57938811465868356),
      new S2Point(-0.10651593657771009, -0.80806167503227055, 0.57938819853274059),
      new S2Point(-0.10651567693742285, -0.80806182530835402, 0.57938803667826444),
      new S2Point(-0.10651496089498214, -0.80806213485510237, 0.57938773659696563),
      new S2Point(-0.10651453461919227, -0.80806229235522298, 0.57938759530083062),
      new S2Point(-0.10651448583749658, -0.80806230280784852, 0.57938758969074455),
      new S2Point(-0.10651428153471061, -0.80806061225022852, 0.57938998503506256),
      new S2Point(-0.10651428161845182, -0.8080606122395747, 0.57938998503452654),
      new S2Point(-0.10651427761078044, -0.80806057978063328, 0.57939003104095654),
      new S2Point(-0.10651427761077951, -0.80806057978062562, 0.57939003104096709),
      new S2Point(-0.10651387099203104, -0.8080572864940091, 0.5793946988282096),
      new S2Point(-0.10651387099202798, -0.80805728649398445, 0.57939469882824468),
      new S2Point(-0.10651386444607201, -0.80805723347699177, 0.57939477397218053),
      new S2Point(-0.10651386444607169, -0.8080572334769891, 0.57939477397218409),
      new S2Point(-0.106513765993723, -0.80805643609199118, 0.57939590414857456),
      new S2Point(-0.10651376671438624, -0.8080564359989727, 0.57939590414581921),
      new S2Point(-0.10651368187839319, -0.80805575808078389, 0.57939686520139033),
      new S2Point(-0.10651465698432123, -0.80805552598235797, 0.57939700963750851),
      new S2Point(-0.1065149024434091, -0.80805548225095913, 0.57939702550292815),
      new S2Point(-0.10651504788182964, -0.80805555533715756, 0.5793968968362615),
      new S2Point(-0.10651511658091152, -0.80805559604710031, 0.57939682743066534),
      new S2Point(-0.10651517919248171, -0.80805562751022852, 0.57939677204023521),
      new S2Point(-0.10651528575974038, -0.80805561374213786, 0.57939677165077275),
      new S2Point(-0.10651648823358072, -0.80805539171529139, 0.57939686023850034),
      new S2Point(-0.10651666406737116, -0.80805537863686483, 0.57939684615295572),
      new S2Point(-0.10651674780673852, -0.80805605121551227, 0.57939589274577097),
      new S2Point(-0.10651674667750256, -0.80805605136137271, 0.57939589274994641),
      new S2Point(-0.10651678418140036, -0.80805634336988752, 0.57939547860450136),
      new S2Point(-0.10651680240261223, -0.80805648524178364, 0.57939527739240138),
      new S2Point(-0.10651680240261237, -0.80805648524178486, 0.57939527739239993),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10651727337444802, -0.80806023111043901, 0.57938996657744879),
      new S2Point(-0.10651727440799089, -0.80806022882029649, 0.57938996958144073),
      new S2Point(-0.10651679374955145, -0.80805648637258243, 0.57939527740611751),
      new S2Point(-0.10651677552833975, -0.80805634450068775, 0.57939547861821594),
      new S2Point(-0.10651673802444192, -0.80805605249217261, 0.57939589276366099),
      new S2Point(-0.10651674651102909, -0.80805605138312775, 0.5793958927502102),
      new S2Point(-0.10651673915225639, -0.80805605233507238, 0.57939589277542292),
      new S2Point(-0.10651665541288889, -0.80805537975642383, 0.57939684618260878),
      new S2Point(-0.10651667272185343, -0.80805537751730583, 0.57939684612330267),
      new S2Point(-0.1065167564612207, -0.8080560500959526, 0.57939589271611924),
      new S2Point(-0.1065167553320342, -0.80805605024202609, 0.57939589271998793),
      new S2Point(-0.10651679283446101, -0.80805634223908773, 0.57939547859078699),
      new S2Point(-0.10651681105567287, -0.80805648411098374, 0.57939527737868723),
      new S2Point(-0.10651680240318392, -0.80805648524170914, 0.5793952773924006),
      new S2Point(-0.10651680240261234, -0.80805648524178475, 0.57939527739239982),
      new S2Point(-0.1065168110556733, -0.80805648411098718, 0.57939527737868224),
      new S2Point(-0.10651729169518892, -0.80806022641135866, 0.57938996976297907),
      new S2Point(-0.10651729210462238, -0.80806022661896348, 0.579389969398166),
      new S2Point(-0.1065172934126499, -0.80806022944626155, 0.57938996521453356),
      new S2Point(-0.10651729203606744, -0.80806023249651726, 0.57938996121349717),
      new S2Point(-0.1065172883437291, -0.80806023495241674, 0.57938995846713126),
      new S2Point(-0.10651728332499401, -0.80806023615590394, 0.5793899577113224),
      new S2Point(-0.10651727832462815, -0.80806023578450537, 0.57938995914858893),
      new S2Point(-0.10651727468247554, -0.80806023393773707, 0.57938996239381635),
    },
    new S2Point[]{
      new S2Point(-0.10651680240204828, -0.80805648524185858, 0.57939527739240082),
      new S2Point(-0.10651679861449742, -0.80805648573682254, 0.57939527739840524),
      new S2Point(-0.10651680240261419, -0.80805648524178353, 0.57939527739240138),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Loop 0: Edge 33 crosses edge 35
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug8() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10703872198218529, -0.80846112144645677, 0.57873424566545062),
      new S2Point(-0.10703872122182066, -0.80846111957630917, 0.57873424841857957),
      new S2Point(-0.10703873813385757, -0.80846111582010538, 0.57873425053786276),
      new S2Point(-0.1070387388942222, -0.80846111769025297, 0.57873424778473381),
      new S2Point(-0.10703873050793056, -0.80846111955286837, 0.57873424673382978),
      new S2Point(-0.1070387388942227, -0.80846111769025419, 0.57873424778473193),
      new S2Point(-0.10703919382477994, -0.80846223660916783, 0.57873260056976505),
      new S2Point(-0.10703917691274406, -0.80846224036537406, 0.57873259845047831),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10703917691274355, -0.80846224036537273, 0.57873259845047997),
      new S2Point(-0.1070391853685064, -0.8084622384873289, 0.57873259951008804),
      new S2Point(-0.10703919381027188, -0.80846223657409677, 0.57873260062144094),
      new S2Point(-0.10703919381027233, -0.80846223657409788, 0.57873260062143939),
      new S2Point(-0.10703918536876245, -0.80846223848727206, 0.57873259951012024),
      new S2Point(-0.10703919382478132, -0.80846223660917116, 0.57873260056976017),
      new S2Point(-0.10703957146434441, -0.80846316542623331, 0.57873123320737097),
      new S2Point(-0.10703955455230836, -0.8084631691824391, 0.57873123108808489),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            //  Loop 1: Edge 1 crosses edge 3
        }

        [Fact]
        public void Test_S2Polygon_Bug9() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10639937100501309, -0.80810205676564995, 0.57935329437301375),
      new S2Point(-0.10639937101137514, -0.80810205688156922, 0.57935329421015713),
      new S2Point(-0.10639937101137305, -0.80810205688156944, 0.57935329421015713),
      new S2Point(-0.106399371005011, -0.80810205676565017, 0.57935329437301375),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10639937099530022, -0.8081020567669569, 0.57935329437297489),
      new S2Point(-0.10639937102108385, -0.80810205688026293, 0.5793532942101961),
      new S2Point(-0.10639937102108181, -0.80810205688026326, 0.5793532942101961),
      new S2Point(-0.10639937099529816, -0.80810205676695701, 0.57935329437297478),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Given edges do not form loops (indegree != outdegree)
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug10() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10592889932808099, -0.80701394501854917, 0.58095400922339757),
      new S2Point(-0.10592787800899696, -0.8070140771413753, 0.58095401191158469),
      new S2Point(-0.1059270044681431, -0.80701419014619669, 0.58095401421031945),
      new S2Point(-0.10592685562894633, -0.80701420940058122, 0.58095401460194696),
      new S2Point(-0.10592685502239066, -0.80701420947920588, 0.58095401460332308),
      new S2Point(-0.10592681668594067, -0.80701421444855337, 0.5809540146902914),
      new S2Point(-0.10592586497682262, -0.8070143378130904, 0.58095401684902004),
      new S2Point(-0.10592586434121586, -0.80701433789547994, 0.58095401685046155),
      new S2Point(-0.10592585898876766, -0.80701428569270217, 0.58095409034224832),
      new S2Point(-0.10592585898876755, -0.80701428569270128, 0.58095409034224987),
      new S2Point(-0.10592571912106936, -0.8070129215545373, 0.58095601078971082),
      new S2Point(-0.10592571912106795, -0.80701292155452331, 0.58095601078973025),
      new S2Point(-0.10592546626664477, -0.80701045545315664, 0.58095948256783148),
      new S2Point(-0.10592546630689463, -0.80701045544795602, 0.58095948256771723),
      new S2Point(-0.10592538513536764, -0.80700975616910509, 0.58096046873415197),
      new S2Point(-0.10592564439344856, -0.80700971612782446, 0.58096047708524956),
      new S2Point(-0.1059267844512099, -0.80700966174311928, 0.58096034476466896),
      new S2Point(-0.10592686088387009, -0.80700965393230761, 0.58096034167862642),
      new S2Point(-0.10592691331665709, -0.80700961093727019, 0.58096039184274961),
      new S2Point(-0.10592705773734933, -0.80700947507458121, 0.58096055423665138),
      new S2Point(-0.10592721940752658, -0.80700934249808198, 0.58096070892049412),
      new S2Point(-0.10592756003095027, -0.80700933299293154, 0.58096066001769275),
      new S2Point(-0.10592832507751106, -0.80700935762745474, 0.58096048630521868),
      new S2Point(-0.1059284165295875, -0.80701007424011018, 0.58095947418602778),
      new S2Point(-0.10592841614913188, -0.80701007428931704, 0.58095947418704452),
      new S2Point(-0.10592864947042728, -0.8070119434176124, 0.58095683523192998),
      new S2Point(-0.1059286884898481, -0.80701225600079662, 0.58095639390519271),
      new S2Point(-0.10592868927069989, -0.80701225581371527, 0.58095639402269295),
      new S2Point(-0.10592869427137827, -0.80701225619024619, 0.58095639258785126),
      new S2Point(-0.10592869791375134, -0.80701225804491505, 0.58095638934738025),
      new S2Point(-0.10592869922184817, -0.80701226088076483, 0.5809563851695615),
      new S2Point(-0.10592869922184843, -0.80701226088076705, 0.58095638516955805),
      new S2Point(-0.10592869784516552, -0.80701226393793402, 0.58095638117383475),
      new S2Point(-0.10592869415258396, -0.80701226639725276, 0.58095637843085768),
      new S2Point(-0.10592868991437976, -0.80701226741266929, 0.58095637779310561),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10592564460843924, -0.80700972122716552, 0.58096046996257766),
      new S2Point(-0.10592539435053176, -0.80700975987840939, 0.58096046190138972),
      new S2Point(-0.10592547496472972, -0.80701045435596641, 0.58095948250602925),
      new S2Point(-0.10592546630689462, -0.80701045544795591, 0.58095948256771723),
      new S2Point(-0.10592546630693271, -0.80701045544826022, 0.58095948256728758),
      new S2Point(-0.1059254749287661, -0.80701045440038255, 0.5809594824508878),
      new S2Point(-0.10592572778318898, -0.80701292050174633, 0.58095601067279068),
      new S2Point(-0.1059257191207934, -0.80701292155455673, 0.58095601078973391),
      new S2Point(-0.1059257194541381, -0.80701292151405679, 0.58095601078521419),
      new S2Point(-0.10592572778319062, -0.80701292050176254, 0.58095601067276803),
      new S2Point(-0.10592586765088864, -0.80701428463992497, 0.58095409022530931),
      new S2Point(-0.10592585899855227, -0.80701428569151201, 0.58095409034211776),
      new S2Point(-0.10592585898857355, -0.80701428569272593, 0.58095409034225098),
      new S2Point(-0.10592586765088888, -0.80701428463992686, 0.58095409022530675),
      new S2Point(-0.10592587247896063, -0.80701433172842685, 0.58095402393347073),
      new S2Point(-0.10592681605007616, -0.80701420941876889, 0.58095402179319922),
      new S2Point(-0.10592685438651758, -0.80701420444942229, 0.58095402170623067),
      new S2Point(-0.10592685499307326, -0.80701420437079774, 0.58095402170485466),
      new S2Point(-0.10592685562894634, -0.80701420940058122, 0.58095401460194696),
      new S2Point(-0.10592685499689927, -0.80701420437030225, 0.58095402170484534),
      new S2Point(-0.10592700383609792, -0.80701418511591771, 0.58095402131321794),
      new S2Point(-0.10592787737695626, -0.80701407211109533, 0.58095401901448296),
      new S2Point(-0.10592889869604118, -0.80701393998826909, 0.58095401632629584),
      new S2Point(-0.10592889996012077, -0.80701395004882903, 0.58095400212049919),
      new S2Point(-0.10592787864104941, -0.80701408217165349, 0.58095400480868631),
      new S2Point(-0.10592787800903029, -0.80701407714164064, 0.58095401191120999),
      new S2Point(-0.10592787864103763, -0.80701408217165482, 0.5809540048086862),
      new S2Point(-0.10592700510019466, -0.80701419517647521, 0.58095400710742118),
      new S2Point(-0.1059270044681431, -0.80701419014619669, 0.58095401421031934),
      new S2Point(-0.10592700510018833, -0.8070141951764761, 0.58095400710742118),
      new S2Point(-0.10592685626275877, -0.80701421443063182, 0.58095400749904391),
      new S2Point(-0.10592685565826369, -0.80701421450898914, 0.58095400750041526),
      new S2Point(-0.10592685502239063, -0.80701420947920566, 0.58095401460332308),
      new S2Point(-0.10592685565826078, -0.80701421450898947, 0.58095400750041526),
      new S2Point(-0.10592681732181129, -0.80701421947833718, 0.58095400758738369),
      new S2Point(-0.10592681668594069, -0.80701421444855348, 0.58095401469029151),
      new S2Point(-0.10592681732180521, -0.80701421947833796, 0.58095400758738369),
      new S2Point(-0.10592586561269894, -0.80701434284287321, 0.58095400974611222),
      new S2Point(-0.10592586497746249, -0.80701433781815202, 0.58095401684187198),
      new S2Point(-0.10592586561268771, -0.80701434284287465, 0.58095400974611222),
      new S2Point(-0.10592586497708102, -0.80701434292526464, 0.58095400974755396),
      new S2Point(-0.10592586434121586, -0.80701433789548005, 0.58095401685046166),
      new S2Point(-0.10592585567909471, -0.80701433894825569, 0.58095401696740323),
      new S2Point(-0.1059258503266465, -0.80701428674547793, 0.58095409045919011),
      new S2Point(-0.10592571045894811, -0.80701292260731206, 0.58095601090665361),
      new S2Point(-0.10592571912060067, -0.80701292155459425, 0.58095601078971715),
      new S2Point(-0.10592571878923682, -0.80701292159485349, 0.58095601079421),
      new S2Point(-0.10592571045894694, -0.80701292260730051, 0.58095601090666993),
      new S2Point(-0.10592545760452345, -0.80701045650593073, 0.58095948268477515),
      new S2Point(-0.10592545764454649, -0.80701045650106651, 0.58095948268423492),
      new S2Point(-0.10592537647753246, -0.80700975726109381, 0.58096046879584118),
      new S2Point(-0.10592538513536764, -0.80700975616910509, 0.58096046873415197),
      new S2Point(-0.10592538413784101, -0.80700975119062324, 0.58096047583161736),
      new S2Point(-0.10592564339592514, -0.80700971114934217, 0.58096048418271495),
      new S2Point(-0.10592564439344856, -0.80700971612782446, 0.58096047708524956),
      new S2Point(-0.10592564496449927, -0.80700971099098684, 0.58096048411668999),
      new S2Point(-0.10592678502227458, -0.80700965660628099, 0.58096035179610783),
      new S2Point(-0.10592678388014524, -0.80700966687995779, 0.58096033773323019),
    },
    new S2Point[]{
      new S2Point(-0.10592585898876757, -0.80701428569270128, 0.58095409034224987),
      new S2Point(-0.10592585897888845, -0.80701428569390288, 0.58095409034238166),
      new S2Point(-0.1059258503266465, -0.80701428674547793, 0.58095409045919011),
    },
    new S2Point[]{
      new S2Point(-0.10592546626664477, -0.80701045545315664, 0.58095948256783148),
      new S2Point(-0.10592546623958927, -0.8070104554564449, 0.58095948256819674),
      new S2Point(-0.10592546626662946, -0.80701045545303429, 0.580959482568004),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            // Inconsistent loop orientations detected
        }

        [Fact]
        public void Test_S2Polygon_Bug11() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10727349803435572, -0.80875763107088172, 0.57827631008375979),
      new S2Point(-0.10727349807040805, -0.80875763112192245, 0.57827631000568813),
      new S2Point(-0.10727349807040625, -0.80875763112192278, 0.57827631000568813),
    },
    new S2Point[]{
      new S2Point(-0.1072729603486537, -0.80875606054879057, 0.57827860629945249),
      new S2Point(-0.10727299870478688, -0.80875633377729705, 0.57827821705818028),
      new S2Point(-0.10727299875560981, -0.80875633413933223, 0.57827821654242495),
      new S2Point(-0.10727309272230967, -0.80875700360375646, 0.57827726282438607),
      new S2Point(-0.10727318660000487, -0.80875767243400742, 0.57827631000742785),
      new S2Point(-0.10727349802669105, -0.80875763101356435, 0.57827631016534387),
      new S2Point(-0.10727349803435525, -0.80875763107087817, 0.57827631008376468),
      new S2Point(-0.10727349803435572, -0.80875763107088172, 0.57827631008375979),
      new S2Point(-0.1072734980420204, -0.80875763112819909, 0.57827631000217561),
      new S2Point(-0.10727318657570066, -0.80875767255391384, 0.57827630984423972),
      new S2Point(-0.10727318651657966, -0.80875767256177711, 0.57827630984420975),
      new S2Point(-0.10727318650891528, -0.80875767250445951, 0.57827630992579371),
      new S2Point(-0.10727318640981781, -0.80875767251785957, 0.57827630992543622),
      new S2Point(-0.10727309252411468, -0.80875700363055636, 0.57827726282367087),
      new S2Point(-0.10727299855741491, -0.8087563341661328, 0.57827821654170874),
      new S2Point(-0.10727299850659211, -0.8087563338040985, 0.57827821705746318),
      new S2Point(-0.10727296014242577, -0.80875606051836801, 0.57827860638025652),
      new S2Point(-0.10727296024152315, -0.80875606050496729, 0.57827860638061501),
      new S2Point(-0.10727296023340849, -0.8087560604477102, 0.57827860646219797),
      new S2Point(-0.10727348576547496, -0.80875598914629976, 0.57827860869282954),
      new S2Point(-0.1072734857817042, -0.80875598926081438, 0.57827860852966395),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.1072734857735896, -0.80875598920355718, 0.5782786086112468),
      new S2Point(-0.10727348576547457, -0.80875598914629976, 0.57827860869282954),
      new S2Point(-0.10727839137361543, -0.80875532356817348, 0.57827862950694298),
      new S2Point(-0.10727839137881608, -0.80875532356471602, 0.57827862951081388),
      new S2Point(-0.10727839143632178, -0.80875532355090063, 0.5782786295194674),
      new S2Point(-0.10727839149361706, -0.80875532355509905, 0.57827862950296649),
      new S2Point(-0.1072783915353497, -0.80875532357618651, 0.57827862946573261),
      new S2Point(-0.10727839154773799, -0.80875532360290581, 0.57827862942606567),
      new S2Point(-0.10727848921795155, -0.80875531035110082, 0.57827862984032907),
      new S2Point(-0.1072784892332832, -0.80875531046514559, 0.57827862967798682),
      new S2Point(-0.10727971608197531, -0.8087551454635169, 0.57827863284376713),
      new S2Point(-0.10727986275126807, -0.80875539440654376, 0.57827825747332484),
      new S2Point(-0.10727959167812619, -0.80875599171505064, 0.57827747239052929),
      new S2Point(-0.10727974196569352, -0.80875625444235633, 0.57827707706958686),
      new S2Point(-0.10727993501555312, -0.80875677560355186, 0.57827631237878363),
      new S2Point(-0.10727870858143702, -0.80875693828645479, 0.57827631237896882),
      new S2Point(-0.1072787085493927, -0.80875693804871851, 0.5782763127174031),
      new S2Point(-0.10727615977928232, -0.80875727704955946, 0.57827631143112901),
      new S2Point(-0.10727615977915911, -0.80875727704957578, 0.57827631143112901),
      new S2Point(-0.10727349803435751, -0.80875763107088128, 0.57827631008375968),
      new S2Point(-0.10727349803435574, -0.80875763107088183, 0.57827631008375979),
      new S2Point(-0.10727318656803594, -0.80875767249659658, 0.57827630992582391),
      new S2Point(-0.10727318650891531, -0.80875767250445962, 0.57827630992579382),
      new S2Point(-0.10727309262321218, -0.80875700361715641, 0.57827726282402847),
      new S2Point(-0.10727299865651231, -0.80875633415273218, 0.57827821654206735),
      new S2Point(-0.10727299860568951, -0.80875633379069789, 0.57827821705782179),
      new S2Point(-0.10727296024152314, -0.80875606050496718, 0.57827860638061501),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Given edges do not form loops (indegree != outdegree)
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_Bug12() {
            S2Point[][] a_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10772916872905106, -0.80699542608967267, 0.58064861015531188),
      new S2Point(-0.10772916892726483, -0.80699542606300401, 0.58064861015560143),
      new S2Point(-0.10772916892726613, -0.80699542606301333, 0.58064861015558844),
      new S2Point(-0.10772916872905235, -0.806995426089682, 0.58064861015529889),
    },
  };
            S2Point[][] b_vertices = new S2Point[][]{
    new S2Point[]{
      new S2Point(-0.10772916872905348, -0.80699542608969022, 0.58064861015528724),
      new S2Point(-0.10772916892726496, -0.80699542606300489, 0.58064861015559999),
      new S2Point(-0.10772930108168739, -0.80699639165138115, 0.58064724364290399),
      new S2Point(-0.10772930088347589, -0.80699639167806647, 0.58064724364259113),
    },
  };
            S2Polygon a = new(MakeLoops(a_vertices));
            S2Polygon b = new(MakeLoops(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            // Given edges do not form loops (indegree != outdegree)
            Assert.False(c.IsEmpty());
        }

        // This tests polygon-polyline intersections.
        // It covers the same edge cases as TestOperations and also adds some
        // extra tests for shared edges.
        [Fact]
        public void Test_S2PolygonTestBase_PolylineIntersection() {
            for (int v = 0; v < 3; ++v) {
                PolylineIntersectionSharedEdgeTest(cross1_, v, 1);
                PolylineIntersectionSharedEdgeTest(cross1_, v + 1, -1);
                PolylineIntersectionSharedEdgeTest(cross1_side_hole_, v, 1);
                PolylineIntersectionSharedEdgeTest(cross1_side_hole_, v + 1, -1);
            }

            // See comments in TestOperations about the vlue of this constant.
            var kMaxError = S1Angle.FromRadians(1e-4);

            // This duplicates some of the tests in TestOperations by
            // converting the outline of polygon A to a polyline then intersecting
            // it with the polygon B. It then converts B to a polyline and intersects
            // it with A. It then feeds all of the results into a polygon builder and
            // tests that the output is equal to doing an intersection between A and B.
            int i = 0;
            foreach (TestCase test in test_cases) {
                _logger.WriteLine("Polyline intersection test case " + i++);
                var a = MakePolygon(test.A);
                var b = MakePolygon(test.B);
                var expected_a_and_b = MakePolygon(test.AAndB);

                var points = new List<S2Point>();
                var polylines = new List<S2Polyline>();
                for (int ab = 0; ab < 2; ab++) {
                    var tmp = ab != 0 ? a : b;
                    var tmp2 = ab != 0 ? b : a;
                    for (int l = 0; l < tmp.NumLoops(); l++) {
                        points.Clear();
                        if (tmp.Loop(l).IsHole()) {
                            for (int v = tmp.Loop(l).NumVertices; v >= 0; v--) {
                                points.Add(tmp.Loop(l).Vertex(v));
                            }
                        } else {
                            for (int v = 0; v <= tmp.Loop(l).NumVertices; v++) {
                                points.Add(tmp.Loop(l).Vertex(v));
                            }
                        }
                        S2Polyline polyline = new(points.ToArray());
                        polylines.AddRange(tmp2.IntersectWithPolyline(polyline));
                    }
                }

                S2Builder builder = new(new Options());
                S2Polygon a_and_b = new();
                builder.StartLayer(new S2PolygonLayer(a_and_b));
                foreach (var polyline in polylines) {
                    builder.AddPolyline(polyline);
                }

                Assert.True(builder.Build(out var error));
                CheckEqual(a_and_b, expected_a_and_b, kMaxError);
            }
        }

        [Fact]
        public void Test_S2PolygonTestBase_Splitting() {
            // It takes too long to test all the polygons in debug mode, so we just pick
            // out some of the more interesting ones.

            SplitAndAssemble(near_10_);
            SplitAndAssemble(near_H3210_);
            SplitAndAssemble(far_H3210_);
            SplitAndAssemble(south_0ab_);
            SplitAndAssemble(south_210b_);
            SplitAndAssemble(south_H20abc_);
            SplitAndAssemble(nf1_n10_f2_s10abc_);
            SplitAndAssemble(nf2_n2_f210_s210ab_);
            SplitAndAssemble(far_H_);
            SplitAndAssemble(south_H_);
            SplitAndAssemble(far_H_south_H_);
        }

        [Fact]
        public void Test_S2Polygon_InitToCellUnionBorder() {
            // Test S2Polygon.InitToCellUnionBorder().
            // The main thing to check is that adjacent cells of different sizes get
            // merged correctly.  To do this we generate two random adjacent cells,
            // convert to polygon, and make sure the polygon only has a single loop.
            for (int iter = 0; iter < 200; ++iter) {
                _logger.WriteLine("Iteration " + iter);

                // Choose a random non-leaf cell.
                S2CellId big_cell =
                    S2Testing.GetRandomCellId(S2Testing.Random.Uniform(S2.kMaxCellLevel));
                // Get all neighbors at some smaller level.
                int small_level = big_cell.Level() +
                    S2Testing.Random.Uniform(Math.Min(16, S2.kMaxCellLevel - big_cell.Level()));
                var neighbors = new List<S2CellId>();
                big_cell.AppendAllNeighbors(small_level, neighbors);
                // Pick one at random.
                S2CellId small_cell = neighbors[S2Testing.Random.Uniform(neighbors.Count)];
                // If it's diagonally adjacent, bail out.
                var edge_neighbors = new S2CellId[4];
                big_cell.EdgeNeighbors(edge_neighbors);
                bool diagonal = true;
                for (int i = 0; i < 4; ++i) {
                    if (edge_neighbors[i].Contains(small_cell)) {
                        diagonal = false;
                    }
                }
                if (diagonal) {
                    continue;
                }

                var cells = new List<S2CellId>
                {
                    big_cell,
                    small_cell
                };
                S2CellUnion cell_union = new(cells);
                Assert.Equal(2, cell_union.Size());
                S2Polygon poly = new();
                poly.InitToCellUnionBorder(cell_union);
                Assert.Equal(1, poly.NumLoops());
                // If the conversion were perfect we could test containment, but due to
                // rounding the polygon won't always exactly contain both cells.  We can
                // at least test intersection.
                Assert.True(poly.MayIntersect(new S2Cell(big_cell)));
                Assert.True(poly.MayIntersect(new S2Cell(small_cell)));
            }
        }

        [Fact]
        public void Test_S2Polygon_UnionWithAmbgiuousCrossings() {
            S2Point[] a_vertices = {
    new S2Point(0.044856812877680216, -0.80679210859571904, 0.5891301722422051),
    new S2Point(0.044851868273159699, -0.80679240802900054, 0.5891301386444033),
    new S2Point(0.044854246527738666, -0.80679240292188514, 0.58912996457145106)
  };
            S2Point[] b_vertices = {
    new S2Point(0.044849715793028468, -0.80679253837178111, 0.58913012401412856),
    new S2Point(0.044855344598821352, -0.80679219751320641, 0.589130162266992),
    new S2Point(0.044854017712818696, -0.80679210327223405, 0.58913039235179754)
  };
            S2Polygon a = new(new S2Loop(a_vertices));
            S2Polygon b = new(new S2Loop(b_vertices));
            S2Polygon c = new();
            c.InitToUnion(a, b);
            Assert.False(c.IsEmpty());
        }

        [Fact]
        public void Test_S2Polygon_InitToSloppySupportsEmptyPolygons() {
            S2Polygon empty_polygon = new();
            S2Polygon polygon = new();
            polygon.InitToSnapped(empty_polygon);
            // InitToSloppy is further tested by SnapSplitsPolygon.
        }

        [Fact]
        public void Test_S2Polygon_InitToSnappedDoesNotRotateVertices() {
            // This particular example came from MapFacts, but in fact InitToSnapped
            // used to cyclically rotate the vertices of all "hole" loops.
            var polygon = MakePolygonOrDie(
      "49.9305505:-124.8345463, 49.9307448:-124.8299657, " +
      "49.9332101:-124.8301996, 49.9331224:-124.8341368; " +
      "49.9311087:-124.8327042, 49.9318176:-124.8312621, " +
      "49.9318866:-124.8334451");
            S2Polygon polygon2 = new(), polygon3 = new();
            polygon2.InitToSnapped(polygon);

            // Check that the first vertex is the same when converted to E7.
            Assert.Equal(S2LatLng.Latitude(polygon.Loop(0).Vertex(0)).E7(),
                      S2LatLng.Latitude(polygon2.Loop(0).Vertex(0)).E7());
            Assert.Equal(S2LatLng.Longitude(polygon.Loop(0).Vertex(0)).E7(),
                      S2LatLng.Longitude(polygon2.Loop(0).Vertex(0)).E7());

            // Check that snapping twice doesn't rotate the vertices.
            polygon3.InitToSnapped(polygon2);
            Assert.True(polygon2 == polygon3);
        }

        [Fact]
        public void Test_S2Polygon_InitToSnappedWithSnapLevel() {
            var polygon = MakePolygonOrDie("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");
            for (int level = 0; level <= S2.kMaxCellLevel; ++level) {
                S2Polygon snapped_polygon = new();
                snapped_polygon.InitToSnapped(polygon, level);
                Assert.True(snapped_polygon.IsValid());
                S1Angle merge_radius = new[]{
                    S1Angle.FromRadians(S2.kMaxDiag.GetValue(level)),
                    SnapFunction.kMaxSnapRadius}.Min();
                Assert.True(snapped_polygon.ApproxContains(polygon, merge_radius));
            }
        }

        [Fact]
        public void Test_S2Polygon_InitToSnappedIsValid_A() {
            var poly = MakePolygonOrDie(
      "53.1328020478452:6.39444903453293, 53.1328019:6.394449, " +
      "53.1327091:6.3961766, 53.1313753:6.3958652, 53.1312825:6.3975924, " +
      "53.132616:6.3979042, 53.1326161348736:6.39790423150577");
            Assert.True(poly.IsValid());
            S2Polygon poly_snapped = new();
            poly_snapped.InitToSnapped(poly);
            Assert.False(poly_snapped.FindValidationError(out _));
        }

        [Fact]
        public void Test_S2Polygon_InitToSnappedIsValid_B() {
            var poly = MakePolygonOrDie(
      "51.6621651:4.9858102, 51.6620965:4.9874227, 51.662028:4.9890355, " +
      "51.6619796006122:4.99017864445347, 51.6622335420397:4.98419752545216, " +
      "51.6622334:4.9841975; 51.66189957578:4.99206198576131, " +
      "51.6618911:4.9922612, 51.6618224:4.9938741, 51.6605122:4.993639, " +
      "51.6604437:4.9952519, 51.6603751:4.9968648, 51.6603064:4.9984777, " +
      "51.6602379:5.0000907, 51.660169:5.0017037, 51.6601003:5.0033165, " +
      "51.6600318:5.0049298, 51.659963:5.0065427, 51.6598943:5.0081561, " +
      "51.6612044207178:5.00839208571886, 51.6612732068132:5.00677860122814, " +
      "51.6612732:5.0067786, 51.6613418:5.0051654, 51.6614106:5.0035525, " +
      "51.6614793:5.0019393, 51.6615479:5.0003263, " +
      "51.6615946694783:4.99923124520759, 51.6616389353165:4.99819106536521, " +
      "51.6616852:4.9971, 51.6617538:4.995487, " +
      "51.661753964726:4.99548702962593");
            Assert.True(poly.IsValid());
            S2Polygon poly_snapped = new();
            poly_snapped.InitToSnapped(poly);
            Assert.False(poly_snapped.FindValidationError(out _));
        }

        [Fact]
        public void Test_S2Polygon_InitToSnappedIsValid_C() {
            var poly = MakePolygonOrDie(
      "53.5316236236404:19.5841192796855, 53.5416584:19.5915903, " +
      "53.5416584189104:19.5915901888287; 53.5416584:19.5915903, " +
      "53.5363122:19.62299, 53.5562817:19.6378935, 53.5616342:19.606474; " +
      "53.5616342:19.606474, 53.5916039:19.6288326, 53.5912689:19.6307982, " +
      "53.5925176:19.6317308, 53.5928526:19.6297652, 53.6015949:19.6362943, " +
      "53.6015950436033:19.6362944072725, 53.6015950814439:19.6362941852262, " +
      "53.5616342380536:19.6064737764314");
            Assert.True(poly.IsValid());
            S2Polygon poly_snapped = new();
            poly_snapped.InitToSnapped(poly);
            Assert.False(poly_snapped.FindValidationError(out _));
        }

        [Fact]
        public void Test_S2Polygon_InitToSnappedIsValid_D() {
            var poly = (MakePolygonOrDie(
      "52.0909316:4.8673826, 52.0909317627574:4.86738262858533, " +
      "52.0911338452911:4.86248482549567, 52.0911337:4.8624848, " +
      "52.0910665:4.8641176, 52.090999:4.8657502"));
            Assert.True(poly.IsValid());
            S2Polygon poly_snapped = new();
            poly_snapped.InitToSnapped(poly);
            Assert.False(poly_snapped.FindValidationError(out _));
        }

        [Fact]
        public void Test_S2Polygon_MultipleInit() {
            var polygon = MakePolygonOrDie("0:0, 0:2, 2:0");
            Assert.Equal(1, polygon.NumLoops());
            Assert.Equal(3, polygon.NumVertices);
            S2LatLngRect bound1 = polygon.GetRectBound();

            var loops = new List<S2Loop>
            {
                MakeLoopOrDie("10:0, -10:-20, -10:20"),
                MakeLoopOrDie("40:30, 20:10, 20:50")
            };
            polygon = new S2Polygon(loops);
            Assert.True(polygon.IsValid());
            Assert.Equal(2, polygon.NumLoops());
            Assert.Equal(6, polygon.NumVertices);
            Assert.True(bound1 != polygon.GetRectBound());
        }

        [Fact]
        public void Test_S2Polygon_InitSingleLoop() {
            S2Polygon polygon = new(S2Loop.kEmpty);
            Assert.True(polygon.IsEmpty());
            polygon = new S2Polygon(S2Loop.kFull);
            Assert.True(polygon.IsFull());
            polygon = new S2Polygon(MakeLoopOrDie("0:0, 0:10, 10:0"));
            Assert.Equal(3, polygon.NumVertices);
        }

        [Fact]
        public void Test_S2PolygonTestBase_TestSimpleEncodeDecode() {
            Encoder encoder = new();
            cross1_.Encode(encoder);
            var decoder = encoder.Decoder();
            var (success, decoded_polygon) = S2Polygon.Decode(decoder);
            Assert.True(success);
            Assert.True(cross1_.BoundaryEquals(decoded_polygon));
            Assert.Equal(cross1_.GetRectBound(), decoded_polygon.GetRectBound());
        }

        [Fact]
        public void Test_S2Polygon_TestEncodeDecodeDefaultPolygon() {
            S2Polygon polygon = new();
            Assert.True(TestEncodeDecode(polygon));
        }

        [Fact]
        public void Test_S2Polygon_CompressedEmptyPolygonRequires3Bytes() {
            S2Polygon empty_polygon = new();
            Encoder encoder = new();

            S2Polygon snapped_empty_polygon = new();
            snapped_empty_polygon.InitToSnapped(empty_polygon);

            snapped_empty_polygon.Encode(encoder);
            // 1 byte for version, 1 for the level, 1 for the length.
            Assert.Equal(1 + 1 + 1, encoder.Length());

            Assert.True(snapped_empty_polygon.IsEmpty());
            Assert.Equal(S2LatLngRect.Empty, snapped_empty_polygon.GetRectBound());
        }

        [Fact]
        public void Test_S2Polygon_CompressedEncodedPolygonRequires69Bytes() {
            var polygon = MakePolygonOrDie("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");

            S2Polygon snapped_polygon = new();
            snapped_polygon.InitToSnapped(polygon);

            Encoder encoder = new();
            snapped_polygon.Encode(encoder);

            // 2 loops, one with 3 vertices, one with 4.
            // Polygon:
            //   1 byte for version
            //   1 byte for level
            //   1 byte for num_loops
            // Loops:
            //   5 bytes overhead
            //   8 bytes per vertex
            Assert.Equal(1 + 1 + 1 + 2 * 5 + 7 * 8, encoder.Length());
        }

        [Fact]
        public void Test_S2PolygonTestBase_CompressedEncodedPolygonDecodesApproxEqual() {
            // To compare the boundaries, etc we want to snap first.
            S2Polygon snapped = new();
            snapped.InitToSnapped(near_30_);
            Assert.Equal(2, snapped.NumLoops());
            Assert.Equal(0, snapped.Loop(0).Depth);
            Assert.Equal(1, snapped.Loop(1).Depth);

            Encoder encoder = new();
            snapped.Encode(encoder);

            var decoder = encoder.Decoder();

            var (success, decoded_polygon) = S2Polygon.Decode(decoder);
            Assert.True(success);
            Assert.True(decoded_polygon.IsValid());
            Assert.True(snapped.BoundaryEquals(decoded_polygon));
            Assert.Equal(snapped.GetRectBound(), decoded_polygon.GetRectBound());
            Assert.Equal(snapped.NumVertices, decoded_polygon.NumVertices);
            Assert.Equal(2, decoded_polygon.NumLoops());
            Assert.Equal(0, decoded_polygon.Loop(0).Depth);
            Assert.Equal(1, decoded_polygon.Loop(1).Depth);
        }

        // This test checks that S2Polygons created directly from S2Cells behave
        // identically to S2Polygons created from the vertices of those cells; this
        // previously was not the case, because S2Cells calculate their bounding
        // rectangles slightly differently, and S2Polygons created from them just
        // copied the S2Cell bounds.
        [Fact]
        public void Test_S2Polygon_TestS2CellConstructorAndContains() {
            S2LatLng latlng = new(S1Angle.FromE6(40565459), S1Angle.FromE6(-74645276));
            S2Cell cell = new(latlng);
            S2Polygon cell_as_polygon = new(cell);
            S2Polygon empty = new();
            S2Polygon polygon_copy = new();
            polygon_copy.InitToUnion(cell_as_polygon, empty);
            Assert.True(polygon_copy.Contains(cell_as_polygon));
            Assert.True(cell_as_polygon.Contains(polygon_copy));
        }

        [Fact]
        public void Test_S2PolygonTest_Project() {
            var polygon = MakePolygon(kNear0 + kNear2);
            S2Point point;
            S2Point projected;

            // The point inside the polygon should be projected into itself.
            point = MakePointOrDie("1.1:0");
            projected = polygon.Project(point);
            Assert.True(S2.ApproxEquals(point, projected));

            // The point is on the outside of the polygon.
            point = MakePointOrDie("5.1:-2");
            projected = polygon.Project(point);
            Assert.True(S2.ApproxEquals(MakePointOrDie("5:-2"), projected));

            // The point is inside the hole in the polygon.
            point = MakePointOrDie("-0.49:-0.49");
            projected = polygon.Project(point);
            Assert.True(S2.ApproxEquals(MakePointOrDie("-0.5:-0.5"),
                                         projected, S1Angle.FromRadians(1e-6)));

            point = MakePointOrDie("0:-3");
            projected = polygon.Project(point);
            Assert.True(S2.ApproxEquals(MakePointOrDie("0:-2"), projected));
        }

        [Fact]
        public void Test_S2PolygonTestBase_GetDistance() {
            // The empty and full loops don't have boundaries.
            TestDistanceMethods(empty_, new S2Point(0, 1, 0), new S2Point());
            TestDistanceMethods(full_, new S2Point(0, 1, 0), new S2Point());

            // A polygon consisting of two nested rectangles centered around
            // S2LatLng(0,0).  Note that because lines of latitude are curved on the
            // sphere, it is not straightforward to project points onto any edge except
            // along the equator.  (The equator is the only line of latitude that is
            // also a geodesic.)
            var nested = (MakePolygonOrDie(
      "3:1, 3:-1, -3:-1, -3:1; 4:2, 4:-2, -4:-2, -4:2;"));

            // All points on the boundary of the polygon should be at distance zero.
            for (int i = 0; i < nested.NumLoops(); i++) {
                var loop = nested.Loop(i);
                for (int j = 0; j < loop.NumVertices; j++) {
                    // A vertex.
                    TestDistanceMethods(nested, loop.Vertex(j), new S2Point());
                    // A point along an edge.
                    TestDistanceMethods(nested, S2.Interpolate(
                        S2Testing.Random.RandDouble(), loop.Vertex(j), loop.Vertex(j + 1)),
                                        new S2Point());
                }
            }
            // A point outside the outer shell that projects to an edge.
            TestDistanceMethods(nested, S2LatLng.FromDegrees(0, -4.7).ToPoint(),
                                S2LatLng.FromDegrees(0, -2).ToPoint());
            // A point outside the outer shell that projects to a vertex.
            TestDistanceMethods(nested, S2LatLng.FromDegrees(6, -3).ToPoint(),
                                S2LatLng.FromDegrees(4, -2).ToPoint());
            // A point inside the polygon that projects to an outer edge.
            TestDistanceMethods(nested, S2LatLng.FromDegrees(0, 1.7).ToPoint(),
                                S2LatLng.FromDegrees(0, 2).ToPoint());
            // A point inside the polygon that projects to an inner vertex.
            TestDistanceMethods(nested, S2LatLng.FromDegrees(-3.3, -1.3).ToPoint(),
                                S2LatLng.FromDegrees(-3, -1).ToPoint());
            // A point inside the inner hole.
            TestDistanceMethods(nested, S2LatLng.FromDegrees(0, 0.1).ToPoint(),
                                S2LatLng.FromDegrees(0, 1).ToPoint());
        }

        [Fact]
        public void Test_S2PolygonTestBase_Area() {
            Assert2.Near(0.0, empty_.GetArea());
            Assert2.Near(S2.M_4_PI, full_.GetArea());
            Assert2.Near(S2.M_2_PI, south_H_.GetArea());
            Assert2.Near(Math.PI, far_H_south_H_.GetArea());

            var two_shells = MakePolygon(kCross1SideHole +kCrossCenterHole);
            Assert2.Near(
                two_shells.Loop(0).Area() + two_shells.Loop(1).Area(),
                two_shells.GetArea());

            var holey_shell = MakePolygon(kCross1 +kCrossCenterHole);
            Assert2.Near(
                holey_shell.Loop(0).Area() - holey_shell.Loop(1).Area(),
                holey_shell.GetArea());
        }

        [Fact]
        public void Test_S2Polygon_UninitializedIsValid() {
            S2Polygon polygon = new();
            Assert.True(polygon.IsValid());
        }

        [Fact]
        public void Test_IsValidTest_UnitLength() {
            // This test can only be run in optimized builds because there are
            // Assert.True(IsUnitLength()) calls scattered throughout the S2 code.
#if !DEBUG
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(1 + S2Testing.Random.Uniform(6), 3 /*min_vertices*/);
                var vloop = vloops_[S2Testing.Random.Uniform(vloops_.Count)];
                S2Point p = vloop[S2Testing.Random.Uniform(vloop.Count)];
                switch (S2Testing.Random.Uniform(3)) {
                    case 0: p = S2Point.Empty; break;
                    case 1: p *= 1e-30 * Math.Pow(1e60, S2Testing.Random.RandDouble()); break;
                    case 2: p = double.NaN * S2Point.Empty; break;
                }
                vloop[S2Testing.Random.Uniform(vloop.Count)] = p;
                CheckInvalid("unit length");
            }
#endif
        }

        [Fact]
        public void Test_IsValidTest_VertexCount() {
            for (int iter = 0; iter < kIters; ++iter) {
                var vloop = new List<S2Point>();
                if (S2Testing.Random.OneIn(2)) {
                    vloop.Add(S2Testing.RandomPoint());
                    vloop.Add(S2Testing.RandomPoint());
                }
                AddLoop(vloop);
                CheckInvalid("at least 3 vertices");
            }
        }

        [Fact]
        public void Test_IsValidTest_DuplicateVertex() {
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(1, 3 /*min_vertices*/);
                var vloop = vloops_[0];
                int n = vloop.Count;
                int i = S2Testing.Random.Uniform(n);
                int j = S2Testing.Random.Uniform(n - 1);
                vloop[i] = vloop[j + ((j >= i) ? 1 : 0)];
                CheckInvalid("duplicate vertex");
            }
        }

        [Fact]
        public void Test_IsValidTest_SelfIntersection() {
            for (int iter = 0; iter < kIters; ++iter) {
                // Use multiple loops so that we can test both holes and shells.  We need
                // at least 5 vertices so that the modified edges don't intersect any
                // nested loops.
                AddConcentricLoops(1 + S2Testing.Random.Uniform(6), 5 /*min_vertices*/);
                var vloop = vloops_[S2Testing.Random.Uniform(vloops_.Count)];
                int n = vloop.Count;
                int i = S2Testing.Random.Uniform(n);
                var tmp = vloop[i]; vloop[i] = vloop[(i + 1) % n]; vloop[(i + 1) % n] = tmp;
                CheckInvalid("crosses edge");
            }
        }

        [Fact]
        public void Test_IsValidTest_EmptyLoop() {
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(S2Testing.Random.Uniform(5), 3 /*min_vertices*/);
                AddLoop(S2Loop.kEmpty.Vertices);
                CheckInvalid("empty loop");
            }
        }

        [Fact]
        public void Test_IsValidTest_FullLoop() {
            for (int iter = 0; iter < kIters; ++iter) {
                // This is only an error if there is at least one other loop.
                AddConcentricLoops(1 + S2Testing.Random.Uniform(5), 3 /*min_vertices*/);
                AddLoop(S2Loop.kFull.Vertices);
                CheckInvalid("full loop");
            }
        }

        [Fact]
        public void Test_IsValidTest_LoopsCrossing() {
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(2, 4 /*min_vertices*/);
                // Both loops have the same number of vertices, and vertices at the same
                // index position are collinear with the center point, so we can create a
                // crossing by simply exchanging two vertices at the same index position.
                int n = vloops_[0].Count;
                int i = S2Testing.Random.Uniform(n);
                var tmp = vloops_[0][i]; vloops_[0][i] = vloops_[1][i]; vloops_[1][i] = tmp;
                if (S2Testing.Random.OneIn(2)) {
                    // By copy the two adjacent vertices from one loop to the other, we can
                    // ensure that the crossings happen at vertices rather than edges.
                    vloops_[0][(i + 1) % n] = vloops_[1][(i + 1) % n];
                    vloops_[0][(i + n - 1) % n] = vloops_[1][(i + n - 1) % n];
                }
                CheckInvalid("crosses loop");
            }
        }

        [Fact]
        public void Test_IsValidTest_DuplicateEdge() {
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(2, 4 /*min_vertices*/);
                int n = vloops_[0].Count;
                if (S2Testing.Random.OneIn(2)) {
                    // Create a shared edge (same direction in both loops).
                    int i = S2Testing.Random.Uniform(n);
                    vloops_[0][i] = vloops_[1][i];
                    vloops_[0][(i + 1) % n] = vloops_[1][(i + 1) % n];
                } else {
                    // Create a reversed edge (opposite direction in each loop) by cutting
                    // loop 0 into two halves along one of its diagonals and replacing both
                    // loops with the result.
                    int split = 2 + S2Testing.Random.Uniform(n - 3);
                    vloops_[1].Clear();
                    vloops_[1].Add(vloops_[0][0]);
                    for (int i = split; i < n; ++i) {
                        vloops_[1].Add(vloops_[0][i]);
                    }
                    vloops_[0].Capacity = split + 1;
                }
                CheckInvalid("has duplicate");
            }
        }

        [Fact]
        public void Test_IsValidTest_InconsistentOrientations() {
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(2 + S2Testing.Random.Uniform(5), 3 /*min_vertices*/);
                init_oriented_ = true;
                CheckInvalid("Inconsistent loop orientations");
            }
        }

        [Fact]
        public void Test_IsValidTest_LoopDepthNegative() {
            modify_polygon_hook_ = SetInvalidLoopDepth;
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(1 + S2Testing.Random.Uniform(4), 3 /*min_vertices*/);
                CheckInvalid("invalid loop depth");
            }
        }

        [Fact]
        public void Test_IsValidTest_LoopNestingInvalid() {
            modify_polygon_hook_ = SetInvalidLoopNesting;
            for (int iter = 0; iter < kIters; ++iter) {
                AddConcentricLoops(2 + S2Testing.Random.Uniform(4), 3 /*min_vertices*/);
                // Randomly invert all the loops in order to generate cases where the
                // outer loop encompasses almost the entire sphere.  This tests different
                // code paths because bounding box checks are not as useful.
                if (S2Testing.Random.OneIn(2)) {
                    foreach (var loop in vloops_) {
                        loop.Reverse();
                    }
                }
                CheckInvalid("Invalid nesting");
            }
        }

        [Fact]
        public void Test_IsValidTest_FuzzTest() {
            // Check that the S2Loop/S2Polygon constructors and IsValid() don't crash
            // when they receive arbitrary invalid input.  (We don't test large inputs;
            // it is assumed that the client enforces their own size limits before even
            // attempting to construct geometric objects.)
#if !DEBUG
            for (int iter = 0; iter < kIters; ++iter) {
                int num_loops = 1 + S2Testing.Random.Uniform(10);
                for (int i = 0; i < num_loops; ++i) {
                    int num_vertices = S2Testing.Random.Uniform(10);
                    var vloop = new List<S2Point>();
                    while (vloop.Count < num_vertices) {
                        // Since the number of vertices is random, we automatically test empty
                        // loops, full loops, and invalid vertex counts.  Also since most
                        // vertices are random, we automatically get self-intersections and
                        // loop crossings.  That leaves zero and NaN vertices, duplicate
                        // vertices, and duplicate edges to be created explicitly.
                        if (S2Testing.Random.OneIn(10)) {
                            // Zero vertex.
                            vloop.Add(S2Point.Empty);
                        } else if (S2Testing.Random.OneIn(10)) {
                            // NaN vertex.
                            vloop.Add(double.NaN * S2Point.Empty);
                        } else if (S2Testing.Random.OneIn(10) && vloop.Any()) {
                            // Duplicate vertex.
                            vloop.Add(vloop[S2Testing.Random.Uniform(vloop.Count)]);
                        } else if (S2Testing.Random.OneIn(10) && vloop.Count + 2 <= num_vertices) {
                            // Try to copy an edge from a random loop.
                            var other = vloops_[S2Testing.Random.Uniform(vloops_.Count)];
                            int n = other.Count;
                            if (n >= 2) {
                                int k0 = S2Testing.Random.Uniform(n);
                                int k1 = (k0 + 1) % n;
                                if (S2Testing.Random.OneIn(2)) { var tmp = k0; k0 = k1; k1 = tmp; }  // Copy reversed edge.
                                vloop.Add(other[k0]);
                                vloop.Add(other[k1]);
                            }
                        } else {
                            // Random non-unit-length point.
                            S2Point p = S2Testing.RandomPoint();
                            vloop.Add(1e-30 * Math.Pow(1e60, S2Testing.Random.RandDouble()) * p);
                        }
                    }
                    AddLoop(vloop);
                }
                CheckInvalid("");  // We could get any error message.
            }
#endif
        }

        [Fact]
        public void Test_S2PolygonSimplifierTest_NoSimplification() {
            SetInput("0:0, 0:20, 20:20, 20:0", 1.0);
            Assert.Equal(4, simplified.NumVertices);

            Assert.Equal(0, MaximumDistanceInDegrees(simplified, original, 0));
            Assert.Equal(0, MaximumDistanceInDegrees(original, simplified, 0));
        }

        // Here, 10:-2 will be removed and  0:0-20:0 will intersect two edges.
        // (The resulting polygon will in fact probably have more edges.)
        [Fact]
        public void Test_S2PolygonSimplifierTest_SimplifiedLoopSelfIntersects() {
            SetInput("0:0, 0:20, 10:-0.1, 20:20, 20:0, 10:-0.2", 0.22);

            // The simplified polygon has the same number of vertices but it should now
            // consists of two loops rather than one.
            Assert.Equal(2, simplified.NumLoops());
            Assert.True(0.22 >= MaximumDistanceInDegrees(simplified, original, 0));
            Assert.True(0.22 >= MaximumDistanceInDegrees(original, simplified, 0.22));
        }

        [Fact]
        public void Test_S2PolygonSimplifierTest_NoSimplificationManyLoops() {
            SetInput("0:0,    0:1,   1:0;   0:20, 0:21, 1:20; " +
                     "20:20, 20:21, 21:20; 20:0, 20:1, 21:0", 0.01);
            Assert.Equal(0, MaximumDistanceInDegrees(simplified, original, 0));
            Assert.Equal(0, MaximumDistanceInDegrees(original, simplified, 0));
        }

        [Fact]
        public void Test_S2PolygonSimplifierTest_TinyLoopDisappears() {
            SetInput("0:0, 0:1, 1:1, 1:0", 1.1);
            Assert.True(simplified.IsEmpty());
        }

        [Fact]
        public void Test_S2PolygonSimplifierTest_StraightLinesAreSimplified() {
            SetInput("0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0," +
                     "6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1", 0.01);
            Assert.Equal(4, simplified.NumVertices);
        }

        [Fact]
        public void Test_S2PolygonSimplifierTest_EdgeSplitInManyPieces() {
            // near_square's right four-point side will be simplified to a vertical
            // line at lng=7.9, that will cut the 9 teeth of the saw (the edge will
            // therefore be broken into 19 pieces).
            string saw =
                "1:1, 1:8, 2:2, 2:8, 3:2, 3:8, 4:2, 4:8, 5:2, 5:8," +
                "6:2, 6:8, 7:2, 7:8, 8:2, 8:8, 9:2, 9:8, 10:1";
            string near_square =
                "0:0, 0:7.9, 1:8.1, 10:8.1, 11:7.9, 11:0";
            SetInput(saw + ";" + near_square, 0.21);

            Assert.True(simplified.IsValid());
            Assert.True(0.11 >= MaximumDistanceInDegrees(simplified, original, 0));
            Assert.True(0.11 >= MaximumDistanceInDegrees(original, simplified, 0));
            // The resulting polygon's 9 little teeth are very small and disappear
            // due to the vertex_merge_radius of the polygon builder.  There remains
            // nine loops.
            Assert.Equal(9, simplified.NumLoops());
        }

        [Fact]
        public void Test_S2PolygonSimplifierTest_EdgesOverlap() {
            // Two loops, One edge of the second one ([0:1 - 0:2]) is part of an
            // edge of the first one..
            SetInput("0:0, 0:3, 1:0; 0:1, -1:1, 0:2", 0.01);
            var true_poly = MakePolygonOrDie("0:3, 1:0, 0:0, 0:1, -1:1, 0:2");
            Assert.True(simplified.BoundaryApproxEquals(true_poly, S1Angle.FromRadians(1e-15)));
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_PointsOnCellBoundaryKept() {
            S2Cell cell = new(S2CellId.FromToken("89c25c"));
            var polygon = MakeCellPolygon(cell, new[]{ "0.1:0, 0.2:0, 0.2:0.5"});
            S1Angle tolerance = new S1Angle(polygon.Loop(0).Vertex(0), polygon.Loop(0).Vertex(1)) * 1.1;
            S2Polygon simplified = new();
            simplified.InitToSimplified(polygon, new IdentitySnapFunction(tolerance));
            Assert.True(simplified.IsEmpty());
            S2Polygon simplified_in_cell = new();
            simplified_in_cell.InitToSimplifiedInCell(polygon, cell, tolerance);
            Assert.True(simplified_in_cell.BoundaryEquals(polygon));
            Assert.Equal(3, simplified_in_cell.NumVertices);
            Assert.Equal(-1, simplified.GetSnapLevel());
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_PointsInsideCellSimplified() {
            S2CellId cell_id = S2CellId.FromToken("89c25c");
            S2Cell cell = new(cell_id);
            var polygon = MakeCellPolygon(
                cell, new[]{ "0.3:0, 0.4:0, 0.4:0.5, 0.4:0.8, 0.2:0.8"});
            S1Angle tolerance = new S1Angle(polygon.Loop(0).Vertex(0), polygon.Loop(0).Vertex(1)) * 1.1;
            S2Polygon simplified = new();
            simplified.InitToSimplifiedInCell(polygon, cell, tolerance);
            Assert.True(simplified.BoundaryNear(polygon, S1Angle.FromRadians(1e-15)));
            Assert.Equal(4, simplified.NumVertices);
            Assert.Equal(-1, simplified.GetSnapLevel());
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_CellCornerKept() {
            S2Cell cell = new(S2CellId.FromToken("00001"));
            var input = MakeCellPolygon(cell, new[]{ "1:0, 1:0.05, 0.99:0"});
            S1Angle tolerance = 0.02 * new S1Angle(cell.Vertex(0), cell.Vertex(1));
            S2Polygon simplified = new();
            simplified.InitToSimplifiedInCell(input, cell, tolerance);
            Assert.True(simplified.BoundaryNear(input, S1Angle.FromRadians(1e-15)));
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_NarrowStripRemoved() {
            S2Cell cell = new(S2CellId.FromToken("00001"));
            var input = MakeCellPolygon(cell, new[]{ "0.9:0, 0.91:0, 0.91:1, 0.9:1"});
            S1Angle tolerance = 0.02 * new S1Angle(cell.Vertex(0), cell.Vertex(1));
            S2Polygon simplified = new();
            simplified.InitToSimplifiedInCell(input, cell, tolerance);
            Assert.True(simplified.IsEmpty());
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_NarrowGapRemoved() {
            S2Cell cell = new(S2CellId.FromToken("00001"));
            var input = MakeCellPolygon(
                cell, new[] { "0.7:0, 0.75:0, 0.75:1, 0.7:1", "0.76:0, 0.8:0, 0.8:1, 0.76:1"});
            var expected = MakeCellPolygon(cell, new[] { "0.7:0, 0.8:0, 0.8:1, 0.7:1"});
            S1Angle tolerance = 0.02 * new S1Angle(cell.Vertex(0), cell.Vertex(1));
            S2Polygon simplified = new();
            simplified.InitToSimplifiedInCell(input, cell, tolerance);
            Assert.True(simplified.BoundaryNear(expected, S1Angle.FromRadians(1e-15)));
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_CloselySpacedEdgeVerticesKept() {
            S2Cell cell = new(S2CellId.FromToken("00001"));
            var input = MakeCellPolygon(
                cell, new[]{ "0:0.303, 0:0.302, 0:0.301, 0:0.3, 0.1:0.3, 0.1:0.4"});
            S1Angle tolerance = 0.02 * new S1Angle(cell.Vertex(0), cell.Vertex(1));
            S2Polygon simplified = new();
            simplified.InitToSimplifiedInCell(input, cell, tolerance);
            Assert.True(simplified.BoundaryApproxEquals(input, S1Angle.FromRadians(1e-15)));
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_PolylineAssemblyBug() {
            S2Cell cell = new(S2CellId.FromToken("5701"));
            var polygon = MakePolygon(
                "55.8699252:-163.9412145, " + // South-west corner of 5701
                "54.7672352:-166.7579678, " + // North-east corner of 5701
      /* Offending part: a tiny triangle near south-east corner */
                "54.7109214:-164.6376338, " + // forced vertex, on edge 4
                "54.7140193:-164.6398404, " +
                "54.7113202:-164.6374015");  // forced vertex, on edge 4
            S1Angle tolerance = S1Angle.FromRadians(2.138358e-05);  // 136.235m
            S1Angle max_dist = S1Angle.FromRadians(2.821947e-09);  // 18mm
            S2Polygon simplified_in_cell = new();
            simplified_in_cell.InitToSimplifiedInCell(polygon, cell, tolerance,
                                                      max_dist);
            Assert.False(simplified_in_cell.IsEmpty());
        }

        [Fact]
        public void Test_InitToSimplifiedInCell_InteriorEdgesSnappedToBoundary() {
            var polygon = MakePolygonOrDie(
                "37.8011672:-122.3247322, 37.8011648:-122.3247399, " +
                "37.8011647:-122.3247403, 37.8011646:-122.3247408, " +
                "37.8011645:-122.3247411, 37.8011633:-122.3247449, " +
                "37.8011621:-122.3247334");
            S2Cell cell = new(S2CellId.FromDebugString("4/001013300"));
            S1Angle snap_radius = S2Testing.MetersToAngle(1.0);
            S1Angle boundary_tolerance =
                S1Angle.FromRadians(0.5 * S2.kMaxWidth.GetValue(S2.kMaxCellLevel - 1)) +
                IntLatLngSnapFunction.MinSnapRadiusForExponent(7);
            S2Polygon simplified_polygon = new();
            simplified_polygon.InitToSimplifiedInCell(polygon, cell, snap_radius,
                                                      boundary_tolerance);
            Assert.False(simplified_polygon.FindValidationError(out _));
        }

        // Tests that a regular polygon with many points gets simplified
        // enough.
        [Fact]
        public void Test_S2PolygonSimplifierTest_LargeRegularPolygon() {
            double kRadius = 2.0;  // in degrees
            int num_initial_points = 1000;
            int num_desired_points = 250;
            double tolerance = 1.05 * kRadius * (1 - Math.Cos(Math.PI / num_desired_points));

            SetInput(MakeRegularPolygon("0:0", num_initial_points, kRadius), tolerance);

            Assert.True(tolerance >= MaximumDistanceInDegrees(simplified, original, 0));
            Assert.True(tolerance >= MaximumDistanceInDegrees(original, simplified, 0));
            Assert.True(250 >= simplified.NumVertices);
            Assert.True(200 <= simplified.NumVertices);
        }

        [Fact]
        public void Test_S2PolygonDecodeTest_FuzzUncompressedEncoding() {
            // Some parts of the S2 library S2_DCHECK on invalid data, even if we set
            // s2debug to false or use S2Polygon.set_s2debug_override. So we
            // only run this test in opt mode.
#if DEBUG
    for (int i = 0; i < 100000; ++i) {
                decodeTest.AppendFakeUncompressedEncodingData();
    decodeTest.Test();
  }
#endif
        }

        [Fact]
        public void Test_S2PolygonDecodeTest_FuzzCompressedEncoding() {
            // Some parts of the S2 library S2_DCHECK on invalid data, even if we set
            // s2debug to false or use S2Polygon.set_s2debug_override. So we
            // only run this test in opt mode.
#if DEBUG
  for (int i = 0; i < 100000; ++i) {
                decodeTest.AppendFakeCompressedEncodingData();
                decodeTest.Test();
  }
#endif
        }

        [Fact]
        public void Test_S2PolygonDecodeTest_FuzzEverything() {
            // Some parts of the S2 library S2_DCHECK on invalid data, even if we set
            // s2debug to false or use S2Polygon.set_s2debug_override. So we
            // only run this test in opt mode.
#if DEBUG
  for (int i = 0; i < 100000; ++i) {
                decodeTest.AppendRandomData();
                decodeTest.Test();
  }
#endif
        }

        [Fact]
        public void Test_S2PolygonTestBase_FullPolygonShape() {
            var shape = new S2Polygon.Shape(full_);
            Assert.Equal(0, shape.NumEdges());
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty());
            Assert.True(shape.IsFull());
            Assert.Equal(1, shape.NumChains());
            Assert.Equal(0, shape.GetChain(0).Start);
            Assert.Equal(0, shape.GetChain(0).Length);
            Assert.True(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2PolygonTestBase_EmptyPolygonShape() {
            var shape = new S2Polygon.Shape(empty_);
            Assert.Equal(0, shape.NumEdges());
            Assert.Equal(2, shape.Dimension());
            Assert.True(shape.IsEmpty());
            Assert.False(shape.IsFull());
            Assert.Equal(0, shape.NumChains());
            Assert.False(shape.GetReferencePoint().Contained);
        }

        [Fact]
        public void Test_S2PolygonTestBase_OneLoopPolygonShape() {
            TestPolygonShape(near_0_);
        }

        [Fact]
        public void Test_S2PolygonTestBase_SeveralLoopPolygonShape() {
            TestPolygonShape(near_3210_);
        }

        [Fact]
        public void Test_S2Polygon_ManyLoopPolygonShape() {
            int kNumLoops = 100;
            int kNumVerticesPerLoop = 6;
            S2Testing.ConcentricLoopsPolygon(new S2Point(1, 0, 0),
                kNumLoops, kNumVerticesPerLoop, out var polygon);
            TestPolygonShape(polygon);
        }

        [Fact]
        public void Test_S2PolygonOwningShape_Ownership() {
            // Debug mode builds will catch any memory leak below.
            var loops = new List<S2Loop>();
            var polygon = new S2Polygon(loops);
            _ = new S2Polygon.OwningShape(polygon);
        }

        [Fact]
        public void Test_S2Polygon_PointInBigLoop() {
            // This code used to demonstrate a bug in S2ShapeIndex.
            S2LatLng center = S2LatLng.FromRadians(0.3, 2);
            S1Angle radius = S1Angle.FromDegrees(80);
            S2Polygon poly = new(S2Loop.MakeRegularLoop(center.ToPoint(), radius, 10));
            Assert.True(poly.MayIntersect(new S2Cell(new S2CellId(center))));
        }

        [Fact]
        public void Test_S2PolygonTestBase_IndexContainsOnePolygonShape() {
            MutableS2ShapeIndex index = near_0_.Index;
            Assert.Equal(1, index.NumShapeIds());
            var shape = (S2Polygon.Shape)index.Shape(0);
            Assert.Equal(near_0_, shape.Polygon);
        }

        [Fact]
        public void Test_S2PolygonTestBase_PolygonPolygonDistance() {
            // Verify that the example code for S2Polygon.index() actually works.
            S2Polygon polygon1 = near_0_;
            S2Polygon polygon2 = far_10_;
            S2ClosestEdgeQuery query = new(polygon1.Index);
            var target = new S2ClosestEdgeQuery.ShapeIndexTarget(polygon2.Index);
            S1ChordAngle distance = query.GetDistance(target);
            Assert.True(distance > new S1ChordAngle(S1Angle.FromDegrees(175)));
        }

        private static bool TestEncodeDecode(S2Polygon src)
        {
            Encoder encoder = new();
            src.Encode(encoder);
            var decoder = encoder.Decoder();
            var (_, dst) = S2Polygon.Decode(decoder);
            return src == dst;
        }

        private static S2Polygon MakePolygon(string str)
        {
            var polygon = MakeVerbatimPolygonOrDie(str);

            // Check that InitToSnapped() is idempotent.
            S2Polygon snapped1 = new(), snapped2 = new();
            snapped1.InitToSnapped(polygon);
            snapped2.InitToSnapped(snapped1);
            Assert.True(snapped1 == snapped2);

            // Check that Decode(Encode(x)) is the identity function.
            Assert.True(TestEncodeDecode(polygon));
            return polygon;
        }

        private static void CheckContains(string a_str, string b_str)
        {
            var a = MakePolygon(a_str);
            var b = MakePolygon(b_str);
            Assert.True(a.Contains(b));
            Assert.True(a.ApproxContains(b, S1Angle.FromRadians(1e-15)));
            Assert.False(a.ApproxDisjoint(b, S1Angle.FromRadians(1e-15)));
        }

        private static void CheckContainsPoint(string a_str, string b_str)
        {
            var a = MakePolygonOrDie(a_str);
            Assert.True(a.Contains(MakePointOrDie(b_str)));
        }

        private static void CheckEqual(S2Polygon a, S2Polygon b) => CheckEqual(a, b, S1Angle.Zero);

        private static void CheckEqual(S2Polygon a, S2Polygon b, S1Angle max_error)
        {
            if (a.BoundaryApproxEquals(b, max_error)) return;

            S2Builder builder = new(new Options());
            S2Polygon a2 = new(), b2 = new();
            builder.StartLayer(new S2PolygonLayer(a2));
            builder.AddPolygon(a);
            Assert.True(builder.Build(out _));
            builder.StartLayer(new S2PolygonLayer(b2));
            builder.AddPolygon(b); 
            Assert.True(builder.Build(out _));
            _ = a.ToDebugString();
            _ = b.ToDebugString();
            _ = a2.ToDebugString();
            _ = b2.ToDebugString();
            Assert.True(a2.BoundaryApproxEquals(b2, max_error));
        }

        private static void CheckComplementary(S2Polygon a, S2Polygon b)
        {
            S2Polygon b1 = new();
            b1.InitToComplement(b);
            CheckEqual(a, b1);
        }

        // Given a pair of polygons where A contains B, check that various identities
        // involving union, intersection, and difference operations hold true.
        private static void TestOneNestedPair(S2Polygon a, S2Polygon b)
        {
            Assert.True(a.Contains(b));
            Assert.Equal(!b.IsEmpty(), a.Intersects(b));
            Assert.Equal(!b.IsEmpty(), b.Intersects(a));

            S2Polygon c = new(), d = new(), e = new(), f = new(), g = new();
            c.InitToUnion(a, b);
            CheckEqual(c, a);

            d.InitToIntersection(a, b);
            CheckEqual(d, b);

            e.InitToDifference(b, a);
            Assert.True(e.IsEmpty());

            f.InitToDifference(a, b);
            g.InitToSymmetricDifference(a, b);
            CheckEqual(f, g);
        }

        // Given a pair of disjoint polygons A and B, check that various identities
        // involving union, intersection, and difference operations hold true.
        private static void TestOneDisjointPair(S2Polygon a, S2Polygon b)
        {
            Assert.False(a.Intersects(b));
            Assert.False(b.Intersects(a));
            Assert.Equal(b.IsEmpty(), a.Contains(b));
            Assert.Equal(a.IsEmpty(), b.Contains(a));

            S2Polygon ab = new(), c = new(), d = new(), e = new(), f = new(), g = new();
            S2Builder builder = new(new Options());
            builder.StartLayer(new S2PolygonLayer(ab));
            builder.AddPolygon(a);
            builder.AddPolygon(b);
            Assert.True(builder.Build(out _));

            c.InitToUnion(a, b);
            CheckEqual(c, ab);

            d.InitToIntersection(a, b);
            Assert.True(d.IsEmpty());

            e.InitToDifference(a, b);
            CheckEqual(e, a);

            f.InitToDifference(b, a);
            CheckEqual(f, b);

            g.InitToSymmetricDifference(a, b);
            CheckEqual(g, ab);
        }

        // Given polygons A and B whose union covers the sphere, check that various
        // identities involving union, intersection, and difference hold true.
        private static void TestOneCoveringPair(S2Polygon a, S2Polygon b)
        {
            Assert.Equal(a.IsFull(), a.Contains(b));
            Assert.Equal(b.IsFull(), b.Contains(a));

            S2Polygon c = new();
            c.InitToUnion(a, b);
            Assert.True(c.IsFull());
        }

        // Given polygons A and B such that both A and its complement intersect both B
        // and its complement, check that various identities involving union,
        // intersection, and difference hold true.
        private static void TestOneOverlappingPair(S2Polygon a, S2Polygon b)
        {
            Assert.False(a.Contains(b));
            Assert.False(b.Contains(a));
            Assert.True(a.Intersects(b));

            S2Polygon c = new(), d = new(), e = new(), f = new(), g = new(), h = new();
            c.InitToUnion(a, b);
            Assert.False(c.IsFull());

            d.InitToIntersection(a, b);
            Assert.False(d.IsEmpty());

            e.InitToDifference(b, a);
            Assert.False(e.IsEmpty());

            f.InitToDifference(a, b);
            g.InitToUnion(e, f);
            h.InitToSymmetricDifference(a, b);
            CheckEqual(g, h);
        }

        // Given a pair of polygons where A contains B, test various identities
        // involving A, B, and their complements.
        private static void TestNestedPair(S2Polygon a, S2Polygon b)
        {
            S2Polygon a1 = new(), b1 = new();
            a1.InitToComplement(a);
            b1.InitToComplement(b);

            TestOneNestedPair(a, b);
            TestOneNestedPair(b1, a1);
            TestOneDisjointPair(a1, b);
            TestOneCoveringPair(a, b1);
        }

        // Given a pair of disjoint polygons A and B, test various identities
        // involving A, B, and their complements.
        private static void TestDisjointPair(S2Polygon a, S2Polygon b)
        {
            S2Polygon a1 = new(), b1 = new();
            a1.InitToComplement(a);
            b1.InitToComplement(b);

            TestOneDisjointPair(a, b);
            TestOneCoveringPair(a1, b1);
            TestOneNestedPair(a1, b);
            TestOneNestedPair(b1, a);
        }

        // Given polygons A and B such that both A and its complement intersect both B
        // and its complement, test various identities involving these four polygons.
        private static void TestOverlappingPair(S2Polygon a, S2Polygon b)
        {
            S2Polygon a1 = new(), b1 = new();
            a1.InitToComplement(a);
            b1.InitToComplement(b);

            TestOneOverlappingPair(a, b);
            TestOneOverlappingPair(a1, b1);
            TestOneOverlappingPair(a1, b);
            TestOneOverlappingPair(a, b1);
        }

        // "a1" is the complement of "a", and "b1" is the complement of "b".
        private static void TestOneComplementPair(S2Polygon a, S2Polygon a1, S2Polygon b, S2Polygon b1)
        {
            // Check DeMorgan's Law and that subtraction is the same as intersection
            // with the complement.  This function is called multiple times in order to
            // test the various combinations of complements.

            S2Polygon a1_or_b = new(), a_and_b1 = new(), a_minus_b = new();
            a_and_b1.InitToIntersection(a, b1);
            a1_or_b.InitToUnion(a1, b);
            a_minus_b.InitToDifference(a, b);

            CheckComplementary(a1_or_b, a_and_b1);
            CheckEqual(a_minus_b, a_and_b1);
        }

        // Test identities that should hold for any pair of polygons A, B and their
        // complements.
        private static void TestComplements(S2Polygon a, S2Polygon b)
        {
            S2Polygon a1 = new(), b1 = new();
            a1.InitToComplement(a);
            b1.InitToComplement(b);

            TestOneComplementPair(a, a1, b, b1);
            TestOneComplementPair(a1, a, b, b1);
            TestOneComplementPair(a, a1, b1, b);
            TestOneComplementPair(a1, a, b1, b);

            // There is a lot of redundancy if we do this test for each complementary
            // pair, so we just do it once instead.
            S2Polygon a_xor_b1 = new(), a1_xor_b = new();
            a_xor_b1.InitToSymmetricDifference(a, b1);
            a1_xor_b.InitToSymmetricDifference(a1, b);
            CheckEqual(a_xor_b1, a1_xor_b);
        }

        private static void TestDestructiveUnion(S2Polygon a, S2Polygon b)
        {
            S2Polygon c = new();
            c.InitToUnion(a, b);
            var polygons = new List<S2Polygon>
            {
                (S2Polygon)a.CustomClone(),
                (S2Polygon)b.CustomClone()
            };
            var c_destructive = S2Polygon.DestructiveUnion(polygons);
            CheckEqual(c, c_destructive);
        }

        private void TestRelationWithDesc(S2Polygon a, S2Polygon b, bool contains, bool contained, bool intersects, string description)
        {
            _logger.WriteLine(description);
            Assert.Equal(contains, a.Contains(b));
            Assert.Equal(contained, b.Contains(a));
            Assert.Equal(intersects, a.Intersects(b));
            if (contains) TestNestedPair(a, b);
            if (contained) TestNestedPair(b, a);
            if (!intersects) TestDisjointPair(a, b);
            if (intersects && !(contains | contained))
            {
                TestOverlappingPair(a, b);  // See TestOverlappingPair for definition
            }
            TestDestructiveUnion(a, b);
            TestComplements(a, b);
        }

        private static List<S2Loop> MakeLoops(S2Point[][] loop_vertices)
        {
            var result = new List<S2Loop>();
            foreach (var vertices in loop_vertices)
            {
                result.Add(new S2Loop(vertices));
                Assert.False(result.Last().FindValidationError(out _));
            }
            return result;
        }

        private void PolylineIntersectionSharedEdgeTest(S2Polygon p, int start_vertex, int direction)
        {
            _logger.WriteLine($@"Polyline intersection shared edge test start={start_vertex} direction=direction");
            S2Point[] points = {p.Loop(0).Vertex(start_vertex),
                            p.Loop(0).Vertex(start_vertex + direction)};
            S2Polyline polyline = new(points);
            if (direction < 0)
            {
                var polylines = p.IntersectWithPolyline(polyline);
                Assert.Empty(polylines);
                polylines = p.SubtractFromPolyline(polyline);
                Assert.Single(polylines);
                Assert.Equal(2, polylines[0].NumVertices());
                Assert.Equal(points[0], polylines[0].Vertex(0));
                Assert.Equal(points[1], polylines[0].Vertex(1));
                Assert.False(p.Intersects(polyline));
                Assert.False(p.Contains(polyline));
            }
            else
            {
                var polylines = p.IntersectWithPolyline(polyline);
                Assert.Single(polylines);
                Assert.Equal(2, polylines[0].NumVertices());
                Assert.Equal(points[0], polylines[0].Vertex(0));
                Assert.Equal(points[1], polylines[0].Vertex(1));
                polylines = p.SubtractFromPolyline(polyline);
                Assert.Empty(polylines);
                Assert.True(p.Intersects(polyline));
                Assert.True(p.Contains(polyline));
            }
        }

        private static void CheckCoveringIsConservative(S2Polygon polygon, List<S2CellId> cells)
        {
            // Check that Contains(S2Cell) and MayIntersect(S2Cell) are implemented
            // conservatively, by comparing against the Contains/Intersect result with
            // the "cell polygon" defined by the four cell vertices.  Please note that
            // the cell polygon is *not* an exact representation of the S2Cell: cell
            // vertices are rounded from their true mathematical positions, which leads
            // to tiny cracks and overlaps between the cell polygons at different cell
            // levels.  That is why Contains(S2Cell) and MayIntersect(S2Cell) cannot be
            // implemented by simply converting the cell to an S2Polygon.  But it is
            // still useful to do this as a sanity check.  In particular:
            //
            //  - If Contains(cell) is true, the polygon must contain the cell polygon.
            //  - If the polygon intersects the cell polygon, then MayIntersect(cell)
            //    must return true.
            //
            foreach (var cell_id in cells)
            {
                var cell = new S2Cell(cell_id);
                var cell_poly = new S2Polygon(cell);
                if (polygon.Contains(cell))
                {
                    Assert.True(polygon.Contains(cell_poly));
                }
                if (polygon.Intersects(cell_poly))
                {
                    Assert.True(polygon.MayIntersect(cell));
                }
            }
        }

        // Remove a random polygon from "pieces" and return it.
        private static S2Polygon ChoosePiece(List<S2Polygon> pieces)
        {
            int i = S2Testing.Random.Uniform(pieces.Count);
            S2Polygon result = pieces[i];
            pieces.RemoveAt(i);
            return result;
        }

        private static void SplitAndAssemble(S2Polygon polygon)
        {
            // Normalize the polygon's loop structure by rebuilding it with S2Builder.
            S2Builder builder = new(new Options());
            S2Polygon expected = new();
            builder.StartLayer(new S2PolygonLayer(expected));
            builder.AddPolygon(polygon);

            Assert.True(builder.Build(out _));

            var debug = false;
#if DEBUG
            debug = true;
#endif

            for (int iter = 0; iter < (debug ? 3 : 10); ++iter)
            {
                S2RegionCoverer coverer = new();
                // Compute the minimum level such that the polygon's bounding
                // cap is guaranteed to be cut.
                double diameter = 2 * polygon.GetCapBound().Radius.Radians();
                int min_level = S2.kMaxWidth.GetLevelForMaxValue(diameter);

                // Now choose a level that has up to 500 cells in the covering.
                int level = min_level + S2Testing.Random.Uniform(debug ? 4 : 6);
                coverer.Options_.MinLevel = (min_level);
                coverer.Options_.MaxLevel = (level);
                coverer.Options_.MaxCells = (500);

                coverer.GetCovering(polygon, out var cells);
                var covering = new S2CellUnion(cells);
                S2Testing.CheckCovering(polygon, covering, false);
                CheckCoveringIsConservative(polygon, cells);
                var pieces = new List<S2Polygon>();
                foreach (var cell_id in cells)
                {
                    S2Cell cell = new(cell_id);
                    S2Polygon window = new(cell);
                    var piece = new S2Polygon();
                    piece.InitToIntersection(polygon, window);
                    pieces.Add(piece);
                }

                // Now we repeatedly remove two random pieces, compute their union, and
                // insert the result as a new piece until only one piece is left.
                //
                // We don't use S2Polygon.DestructiveUnion() because it joins the pieces
                // in a mostly deterministic order.  We don't just call random_shuffle()
                // on the pieces and repeatedly join the last two pieces in the vector
                // because this always joins a single original piece to the current union
                // rather than doing the unions according to a random tree structure.
                while (pieces.Count > 1)
                {
                    S2Polygon a = ChoosePiece(pieces);
                    S2Polygon b = ChoosePiece(pieces);
                    var c = new S2Polygon();
                    c.InitToUnion(a, b);
                    pieces.Add(c);
                }
                S2Polygon result = pieces[0];
                pieces.RemoveAt(pieces.Count - 1);

                // The moment of truth!
                var rs = result.ToDebugString();
                var es = expected.ToDebugString();
                Assert.True(expected.BoundaryNear(result, S1Angle.FromRadians(2e-15)));

                // Check that ApproxEquals produces the same result.
                if (!expected.ApproxEquals(result, S2.kIntersectionMergeRadiusS1Angle))
                {
                    S2Polygon symmetric_difference = new();
                    symmetric_difference.InitToApproxSymmetricDifference(
                        expected, result, S2.kIntersectionMergeRadiusS1Angle);
                    var ss = symmetric_difference.ToDebugString();
                    Assert.True(false);
                }
            }
        }

        // Helper function for testing the distance methods.  "boundary_x" is the
        // expected result of projecting "x" onto the polygon boundary.  For
        // convenience it can be set to S2Point() to indicate that (boundary_x == x).
        private static void TestDistanceMethods(S2Polygon polygon, S2Point x, S2Point boundary_x)
        {
            // This error is not guaranteed by the implementation but is okay for tests.
            S1Angle kMaxError = S1Angle.FromRadians(1e-15);

            if (boundary_x == new S2Point()) boundary_x = x;
            Assert.True(new S1Angle(boundary_x, polygon.ProjectToBoundary(x)) <= kMaxError);

            if (polygon.IsEmpty() || polygon.IsFull())
            {
                Assert.Equal(S1Angle.Infinity, polygon.GetDistanceToBoundary(x));
            }
            else
            {
                // Debug.Assert.Near only works with doubles.
                Assert2.Near(new S1Angle(x, boundary_x).GetDegrees(),
                            polygon.GetDistanceToBoundary(x).GetDegrees(),
                            kMaxError.GetDegrees());
            }
            if (polygon.Contains(x))
            {
                Assert.Equal(S1Angle.Zero, polygon.GetDistance(x));
                Assert.Equal(x, polygon.Project(x));
            }
            else
            {
                Assert.Equal(polygon.GetDistanceToBoundary(x), polygon.GetDistance(x));
                Assert.Equal(polygon.ProjectToBoundary(x), polygon.Project(x));
            }
        }

        private static void SetInvalidLoopDepth(S2Polygon polygon)
        {
            int i = S2Testing.Random.Uniform(polygon.NumLoops());
            if (i == 0 || S2Testing.Random.OneIn(3))
            {
                polygon.Loop(i).Depth = (-1);
            }
            else
            {
                polygon.Loop(i).Depth = (polygon.Loop(i - 1).Depth + 2);
            }
        }

        private static void SetInvalidLoopNesting(S2Polygon polygon)
        {
            int i = S2Testing.Random.Uniform(polygon.NumLoops());
            polygon.Loop(i).Invert();
        }

        // Returns the diameter of a loop (maximum distance between any two
        // points in the loop).
        private static S1Angle LoopDiameter(S2Loop loop)
        {
            S1Angle diameter = S1Angle.Zero;
            for (int i = 0; i < loop.NumVertices; ++i)
            {
                S2Point test_point = loop.Vertex(i);
                for (int j = i + 1; j < loop.NumVertices; ++j)
                {
                    diameter = new[]{ diameter,
                                   S2.GetDistance(test_point, loop.Vertex(j),
                                                           loop.Vertex(j + 1))}.Max();
                }
            }
            return diameter;
        }

        // Returns the maximum distance from any vertex of poly_a to poly_b, that is,
        // the directed Haussdorf distance of the set of vertices of poly_a to the
        // boundary of poly_b.
        //
        // Doesn't consider loops from poly_a that have diameter less than min_diameter
        // in degrees.
        private static double MaximumDistanceInDegrees(S2Polygon poly_a, S2Polygon poly_b, double min_diameter_in_degrees)
        {
            double min_distance = 360;
            bool has_big_loops = false;
            for (int l = 0; l < poly_a.NumLoops(); ++l)
            {
                var a_loop = poly_a.Loop(l);
                if (LoopDiameter(a_loop).GetDegrees() <= min_diameter_in_degrees)
                {
                    continue;
                }
                has_big_loops = true;
                for (int v = 0; v < a_loop.NumVertices; ++v)
                {
                    double distance = poly_b.GetDistance(a_loop.Vertex(v)).GetDegrees();
                    if (distance < min_distance)
                    {
                        min_distance = distance;
                    }
                }
            }
            if (has_big_loops)
            {
                return min_distance;
            }
            else
            {
                return 0.0;  // As if the first polygon were empty.
            }
        }

        // Creates a polygon from loops specified as a comma separated list of u:v
        // coordinates relative to a cell. The loop "0:0, 1:0, 1:1, 0:1" is
        // counter-clockwise.
        private static S2Polygon MakeCellPolygon(S2Cell cell, string[] strs)
        {
            var loops = new List<S2Loop>();
            foreach (var str in strs)
            {
                var points = ParseLatLngsOrDie(str);
                var loop_vertices = new List<S2Point>();
                R2Rect uv = cell.BoundUV;
                foreach (var p in points)
                {
                    double u = p.Lat().GetDegrees(), v = p.Lng().GetDegrees();
                    loop_vertices.Add(S2.FaceUVtoXYZ(cell.Face,
                        uv[0][0] * (1 - u) + uv[0][1] * u,
                        uv[1][0] * (1 - v) + uv[1][1] * v).Normalize());
                }
                loops.Add(new S2Loop(loop_vertices));
            }
            return new S2Polygon(loops);
        }

        private static S2Polygon MakeRegularPolygon(string center, int num_points, double radius_in_degrees)
        {
            S1Angle radius = S1Angle.FromDegrees(radius_in_degrees);
            return new S2Polygon(S2Loop.MakeRegularLoop(MakePointOrDie(center), radius, num_points));
        }

        private record TestCase(string A, string B, string AAndB, string AOrB, string AMinusB, string AXorB);

        private class S2PolygonDecodeTest
        {
            public S2PolygonDecodeTest() {
                data_array_ = new byte[kMaxBytes];
                encoder_ = new Encoder(data_array_, kMaxBytes);
            }

            public void AppendByte(int value)
            {
                encoder_.Put8((sbyte)value);
            }

            public void AppendInt32(int value)
            {
                encoder_.Put32(value);
            }

            public void AppendRandomData(int size)
            {
                for (int i = 0; i < size && encoder_.Avail() > 0; ++i)
                {
                    AppendByte(S2Testing.Random.Uniform(256));
                }
            }

            public void AppendRandomData()
            {
                AppendRandomData(S2Testing.Random.Uniform(kMaxBytes));
            }

            public void AppendFakeUncompressedEncodingData()
            {
                AppendByte(1);                      // polygon number
                AppendByte(0);                      // unused
                AppendByte(0);                      // "has holes" flag
                AppendInt32(PickRandomCount());     // num loops
                AppendByte(1);                      // loop version
                AppendInt32(PickRandomCount());     // num vertices
                AppendRandomData();                 // junk to fill out the buffer
            }

            public void AppendFakeCompressedEncodingData()
            {
                AppendByte(4);                      // polygon number
                AppendByte(S2Testing.Random.Uniform(50));    // snap level
                AppendInt32(PickRandomCount());     // num loops
                AppendInt32(PickRandomCount());     // num vertices
                AppendRandomData();                 // junk to fill out the buffer
            }

            private static Int32 PickRandomCount()
            {
                if (S2Testing.Random.OneIn(10))
                {
                    return -1;
                }
                if (S2Testing.Random.OneIn(10))
                {
                    return 0;
                }
                if (S2Testing.Random.OneIn(10))
                {
                    return 1000000000;
                }
                if (S2Testing.Random.OneIn(2))
                {
                    return S2Testing.Random.Uniform(1000000000);
                }
                return S2Testing.Random.Uniform(1000);
            }

            public bool Test()
            {
                decoder_ = new Decoder(data_array_, 0, encoder_.Length());
                encoder_.Clear();
                var (success, _) = S2Polygon.Decode(decoder_);
                return success;
            }

            // Maximum size of the data array.
            private const int kMaxBytes = 256;

            // The data array.
            private readonly byte[] data_array_;

            // Encoder that is used to put data into the array.
            private readonly Encoder encoder_;

            // Decoder used to extract data from the array.
            private Decoder decoder_;
        };

        private static void TestPolygonShape(S2Polygon polygon)
        {
            Assert.True(!polygon.IsFull());
            var shape = new S2Polygon.Shape(polygon);
            Assert.Equal(polygon, shape.Polygon);
            Assert.Equal(polygon.NumVertices, shape.NumEdges());
            Assert.Equal(polygon.NumLoops(), shape.NumChains());
            for (int e = 0, i = 0; i < polygon.NumLoops(); ++i)
            {
                S2Loop loop_i = polygon.Loop(i);
                Assert.Equal(e, shape.GetChain(i).Start);
                Assert.Equal(loop_i.NumVertices, shape.GetChain(i).Length);
                for (int j = 0; j < loop_i.NumVertices; ++j, ++e)
                {
                    var edge = shape.GetEdge(e);
                    Assert.Equal(loop_i.OrientedVertex(j), edge.V0);
                    Assert.Equal(loop_i.OrientedVertex(j + 1), edge.V1);
                }
            }
            Assert.Equal(2, shape.Dimension());
            Assert.False(shape.IsEmpty());
            Assert.False(shape.IsFull());
            Assert.Equal(polygon.Contains(S2.Origin),
                      shape.GetReferencePoint().Contained);
        }

        #region SetInput

        private void SetInput(S2Polygon poly, double tolerance_in_degrees)
        {
            original = poly;

            simplified = new S2Polygon();
            simplified.InitToSimplified(original, new IdentitySnapFunction(S1Angle.FromDegrees(tolerance_in_degrees)));
        }

        private void SetInput(string poly, double tolerance_in_degrees)
        {
            SetInput(MakePolygonOrDie(poly), tolerance_in_degrees);
        }

        #endregion

        #region IsValidTest

        private void AddLoop(List<S2Point> points)
        {
            vloops_.Add(points);
        }

        private void AddLoop(S2Point[]points)
        {
            vloops_.Add(points.ToList());
        }

        // Create "num_loops" nested regular loops around a common center point.
        // All loops have the same number of vertices (at least "min_vertices").
        // Furthermore, the vertices at the same index position are collinear with
        // the common center point of all the loops.  The loop radii decrease
        // exponentially in order to prevent accidental loop crossings when one of
        // the loops is modified.
        private void AddConcentricLoops(int num_loops, int min_vertices)
        {
            Assert.True(num_loops <= 10);  // Because radii decrease exponentially.
            S2Point center = S2Testing.RandomPoint();
            int num_vertices = min_vertices + S2Testing.Random.Uniform(10);
            for (int i = 0; i < num_loops; ++i)
            {
                S1Angle radius = S1Angle.FromDegrees(80 * Math.Pow(0.1, i));
                AddLoop(S2Testing.MakeRegularPoints(center, radius, num_vertices));
            }
        }

        private void Reset()
        {
            vloops_.Clear();
        }

        private void CheckInvalid(string snippet)
        {
            var loops = vloops_.Select(vloop => new S2Loop(vloop, S2Debug.DISABLE)).ToList();
            // Cannot replace with shuffle (b/65670707) since this uses an
            // incompatible random source which is also used as a source of randomness
            // in the surrounding code.
            // NOLINTNEXTLINE
            loops.Shuffle(S2Testing.Random.Next);
            var polygon = new S2Polygon(loops, init_oriented_);
            modify_polygon_hook_?.Invoke(polygon);
            Assert.True(polygon.FindValidationError(out var error));
            Assert.True(error.Text.IndexOf(snippet) != -1);
            Reset();
        }
        
        #endregion
    }
}
