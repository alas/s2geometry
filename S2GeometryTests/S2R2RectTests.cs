// Most of the S2R2Rect methods have trivial implementations in terms of the
// R1Interval class, so most of the testing is done in that unit test.

using Xunit;
using Xunit.Abstractions;

namespace S2Geometry
{
    public class S2R2RectTests
    {
        private readonly ITestOutputHelper _logger;
        public S2R2RectTests(ITestOutputHelper testOutputHelper)
        {
            _logger = testOutputHelper;
        }

        [Fact]
        public void Test_S2R2Rect_EmptyRectangles() {
            // Test basic properties of empty rectangles.
            S2R2Rect empty = S2R2Rect.Empty;
            Assert.True(empty.IsValid);
            Assert.True(empty.IsEmpty);
        }

        [Fact]
        public void Test_S2R2Rect_ConstructorsAndAccessors() {
            // Check various constructors and accessor methods.
            S2R2Rect d1 = new S2R2Rect(new R2Point(0.1, 0), new R2Point(0.25, 1));
            Assert.Equal(0.1, d1.X.Lo);
            Assert.Equal(0.25, d1.X.Hi);
            Assert.Equal(0.0, d1.Y.Lo);
            Assert.Equal(1.0, d1.Y.Hi);
            Assert.Equal(new R1Interval(0.1, 0.25), d1.X);
            Assert.Equal(new R1Interval(0, 1), d1.Y);
        }

        [Fact]
        public void Test_S2R2Rect_FromCell() {
            // FromCell, FromCellId
            Assert.Equal(new S2R2Rect(new R2Point(0, 0), new R2Point(0.5, 0.5)),
                      S2R2Rect.FromCell(S2Cell.FromFacePosLevel(0, 0, 1)));
            Assert.Equal(new S2R2Rect(new R2Point(0, 0), new R2Point(1, 1)),
                      S2R2Rect.FromCellId(S2CellId.FromFacePosLevel(0, 0, 0)));
        }

        [Fact]
        public void Test_S2R2Rect_FromCenterSize() {
            // FromCenterSize()
            Assert.True(S2R2Rect.FromCenterSize(new R2Point(0.3, 0.5), new R2Point(0.2, 0.4)).
                        ApproxEquals(new S2R2Rect(new R2Point(0.2, 0.3), new R2Point(0.4, 0.7))));
            Assert.True(S2R2Rect.FromCenterSize(new R2Point(1, 0.1), new R2Point(0, 2)).
                        ApproxEquals(new S2R2Rect(new R2Point(1, -0.9), new R2Point(1, 1.1))));
        }

        [Fact]
        public void Test_S2R2Rect_FromPoint() {
            // FromPoint(), FromPointPair()
            S2R2Rect d1 = new S2R2Rect(new R2Point(0.1, 0), new R2Point(0.25, 1));
            Assert.Equal(new S2R2Rect(d1.Lo, d1.Lo), S2R2Rect.FromPoint(d1.Lo));
            Assert.Equal(new S2R2Rect(new R2Point(0.15, 0.3), new R2Point(0.35, 0.9)),
                      S2R2Rect.FromPointPair(new R2Point(0.15, 0.9), new R2Point(0.35, 0.3)));
            Assert.Equal(new S2R2Rect(new R2Point(0.12, 0), new R2Point(0.83, 0.5)),
                      S2R2Rect.FromPointPair(new R2Point(0.83, 0), new R2Point(0.12, 0.5)));
        }

        [Fact]
        public void Test_S2R2Rect_SimplePredicates() {
            // GetCenter(), GetVertex(), Contains(R2Point), InteriorContains(R2Point).
            R2Point sw1 = new R2Point(0, 0.25);
            R2Point ne1 = new R2Point(0.5, 0.75);
            S2R2Rect r1 = new S2R2Rect(sw1, ne1);

            Assert.Equal(new R2Point(0.25, 0.5), r1.GetCenter());
            Assert.Equal(new R2Point(0, 0.25), r1.GetVertex(0));
            Assert.Equal(new R2Point(0.5, 0.25), r1.GetVertex(1));
            Assert.Equal(new R2Point(0.5, 0.75), r1.GetVertex(2));
            Assert.Equal(new R2Point(0, 0.75), r1.GetVertex(3));
            Assert.True(r1.Contains(new R2Point(0.2, 0.4)));
            Assert.False(r1.Contains(new R2Point(0.2, 0.8)));
            Assert.False(r1.Contains(new R2Point(-0.1, 0.4)));
            Assert.False(r1.Contains(new R2Point(0.6, 0.1)));
            Assert.True(r1.Contains(sw1));
            Assert.True(r1.Contains(ne1));
            Assert.False(r1.InteriorContains(sw1));
            Assert.False(r1.InteriorContains(ne1));

            // Make sure that GetVertex() returns vertices in CCW order.
            for (int k = 0; k < 4; ++k) {
                _logger.WriteLine("k="+ k);
                Assert.Equal(1, S2Pred.Sign(S2R2Rect.ToS2Point(r1.GetVertex(k - 1)),
                                          S2R2Rect.ToS2Point(r1.GetVertex(k)),
                                          S2R2Rect.ToS2Point(r1.GetVertex(k + 1))));
            }
        }

        [Fact]
        public void Test_S2R2Rect_IntervalOperations() {
            // Contains(S2R2Rect), InteriorContains(S2R2Rect),
            // Intersects(), InteriorIntersects(), Union(), Intersection().
            //
            // Much more testing of these methods is done in s1interval_test
            // and r1interval_test.

            S2R2Rect empty = S2R2Rect.Empty;
            R2Point sw1 = new R2Point(0, 0.25);
            R2Point ne1 = new R2Point(0.5, 0.75);
            S2R2Rect r1 = new S2R2Rect(sw1, ne1);
            S2R2Rect r1_mid = new S2R2Rect(new R2Point(0.25, 0.5), new R2Point(0.25, 0.5));
            S2R2Rect r_sw1 = new S2R2Rect(sw1, sw1);
            S2R2Rect r_ne1 = new S2R2Rect(ne1, ne1);

            TestIntervalOps(r1, r1_mid, "TTTT", r1, r1_mid);
            TestIntervalOps(r1, r_sw1, "TFTF", r1, r_sw1);
            TestIntervalOps(r1, r_ne1, "TFTF", r1, r_ne1);

            Assert.Equal(new S2R2Rect(new R2Point(0, 0.25), new R2Point(0.5, 0.75)), r1);
            TestIntervalOps(r1, new S2R2Rect(new R2Point(0.45, 0.1), new R2Point(0.75, 0.3)), "FFTT",
                            new S2R2Rect(new R2Point(0, 0.1), new R2Point(0.75, 0.75)),
                            new S2R2Rect(new R2Point(0.45, 0.25), new R2Point(0.5, 0.3)));
            TestIntervalOps(r1, new S2R2Rect(new R2Point(0.5, 0.1), new R2Point(0.7, 0.3)), "FFTF",
                            new S2R2Rect(new R2Point(0, 0.1), new R2Point(0.7, 0.75)),
                            new S2R2Rect(new R2Point(0.5, 0.25), new R2Point(0.5, 0.3)));
            TestIntervalOps(r1, new S2R2Rect(new R2Point(0.45, 0.1), new R2Point(0.7, 0.25)), "FFTF",
                            new S2R2Rect(new R2Point(0, 0.1), new R2Point(0.7, 0.75)),
                            new S2R2Rect(new R2Point(0.45, 0.25), new R2Point(0.5, 0.25)));

            TestIntervalOps(new S2R2Rect(new R2Point(0.1, 0.2), new R2Point(0.1, 0.3)),
                            new S2R2Rect(new R2Point(0.15, 0.7), new R2Point(0.2, 0.8)), "FFFF",
                            new S2R2Rect(new R2Point(0.1, 0.2), new R2Point(0.2, 0.8)),
                            empty);

            // Check that the intersection of two rectangles that overlap in x but not y
            // is valid, and vice versa.
            TestIntervalOps(new S2R2Rect(new R2Point(0.1, 0.2), new R2Point(0.4, 0.5)),
                            new S2R2Rect(new R2Point(0, 0), new R2Point(0.2, 0.1)), "FFFF",
                            new S2R2Rect(new R2Point(0, 0), new R2Point(0.4, 0.5)), empty);
            TestIntervalOps(new S2R2Rect(new R2Point(0, 0), new R2Point(0.1, 0.3)),
                            new S2R2Rect(new R2Point(0.2, 0.1), new R2Point(0.3, 0.4)), "FFFF",
                            new S2R2Rect(new R2Point(0, 0), new R2Point(0.3, 0.4)), empty);
        }

        [Fact]
        public void Test_S2R2Rect_AddPoint() {
            // AddPoint()
            R2Point sw1 = new R2Point(0, 0.25);
            R2Point ne1 = new R2Point(0.5, 0.75);
            S2R2Rect r1 = new S2R2Rect(sw1, ne1);

            S2R2Rect r2 = S2R2Rect.Empty;
            r2.AddPoint(new R2Point(0, 0.25));
            r2.AddPoint(new R2Point(0.5, 0.25));
            r2.AddPoint(new R2Point(0, 0.75));
            r2.AddPoint(new R2Point(0.1, 0.4));
            Assert.Equal(r1, r2);
        }

        [Fact]
        public void Test_S2R2Rect_Project() {
            S2R2Rect r1 = new S2R2Rect(new R1Interval(0, 0.5), new R1Interval(0.25, 0.75));

            Assert.Equal(new R2Point(0, 0.25), r1.Project(new R2Point(-0.01, 0.24)));
            Assert.Equal(new R2Point(0, 0.48), r1.Project(new R2Point(-5.0, 0.48)));
            Assert.Equal(new R2Point(0, 0.75), r1.Project(new R2Point(-5.0, 2.48)));
            Assert.Equal(new R2Point(0.19, 0.75), r1.Project(new R2Point(0.19, 2.48)));
            Assert.Equal(new R2Point(0.5, 0.75), r1.Project(new R2Point(6.19, 2.48)));
            Assert.Equal(new R2Point(0.5, 0.53), r1.Project(new R2Point(6.19, 0.53)));
            Assert.Equal(new R2Point(0.5, 0.25), r1.Project(new R2Point(6.19, -2.53)));
            Assert.Equal(new R2Point(0.33, 0.25), r1.Project(new R2Point(0.33, -2.53)));
            Assert.Equal(new R2Point(0.33, 0.37), r1.Project(new R2Point(0.33, 0.37)));
        }

        [Fact]
        public void Test_S2R2Rect_Expanded() {
            // Expanded()
            Assert.True(S2R2Rect.Empty.Expanded(new R2Point(0.1, 0.3)).IsEmpty);
            Assert.True(S2R2Rect.Empty.Expanded(new R2Point(-0.1, -0.3)).IsEmpty);
            Assert.True(new S2R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                        Expanded(new R2Point(0.1, 0.3)).
                        ApproxEquals(new S2R2Rect(new R2Point(0.1, 0.1), new R2Point(0.4, 1.0))));
            Assert.True(new S2R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                        Expanded(new R2Point(-0.1, 0.3)).IsEmpty);
            Assert.True(new S2R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                        Expanded(new R2Point(0.1, -0.2)).IsEmpty);
            Assert.True(new S2R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                        Expanded(new R2Point(0.1, -0.1)).
                        ApproxEquals(new S2R2Rect(new R2Point(0.1, 0.5), new R2Point(0.4, 0.6))));
            Assert.True(new S2R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).Expanded(0.1).
                        ApproxEquals(new S2R2Rect(new R2Point(0.1, 0.3), new R2Point(0.4, 0.8))));
        }

        [Fact]
        public void Test_S2R2Rect_Bounds() {
            // GetCapBound(), GetRectBound()
            S2R2Rect empty = S2R2Rect.Empty;
            Assert.True(empty.GetCapBound().IsEmpty);
            Assert.True(empty.GetRectBound().IsEmpty);
            Assert.Equal(S2Cap.FromPoint(new S2Point(1, 0, 0)),
                      new S2R2Rect(new R2Point(0.5, 0.5), new R2Point(0.5, 0.5)).GetCapBound());
            Assert.Equal(S2LatLngRect.FromPoint(S2LatLng.FromDegrees(0, 0)),
                      new S2R2Rect(new R2Point(0.5, 0.5), new R2Point(0.5, 0.5)).GetRectBound());

            for (int i = 0; i < 10; ++i) {
                _logger.WriteLine("i="+ i);
                S2R2Rect rect = S2R2Rect.FromCellId(S2Testing.GetRandomCellId());
                S2Cap cap = rect.GetCapBound();
                S2LatLngRect llrect = rect.GetRectBound();
                for (int k = 0; k < 4; ++k) {
                    S2Point v = S2R2Rect.ToS2Point(rect.GetVertex(k));
                    // v2 is a point that is well outside the rectangle.
                    S2Point v2 = (cap.Center + 3 * (v - cap.Center)).Normalized;
                    Assert.True(cap.Contains(v));
                    Assert.False(cap.Contains(v2));
                    Assert.True(llrect.Contains(v));
                    Assert.False(llrect.Contains(v2));
                }
            }
        }

        [Fact]
        public void Test_S2R2Rect_CellOperations() {
            // Contains(S2Cell), MayIntersect(S2Cell)

            S2R2Rect empty = S2R2Rect.Empty;
            TestCellOps(empty, S2Cell.FromFace(3), 0);

            // This rectangle includes the first quadrant of face 0.  It's expanded
            // slightly because cell bounding rectangles are slightly conservative.
            S2R2Rect r4 = new S2R2Rect(new R2Point(0, 0), new R2Point(0.5, 0.5));
            TestCellOps(r4, S2Cell.FromFacePosLevel(0, 0, 0), 3);
            TestCellOps(r4, S2Cell.FromFacePosLevel(0, 0, 1), 4);
            TestCellOps(r4, S2Cell.FromFacePosLevel(1, 0, 1), 0);

            // This rectangle intersects the first quadrant of face 0.
            S2R2Rect r5 = new S2R2Rect(new R2Point(0, 0.45), new R2Point(0.5, 0.55));
            TestCellOps(r5, S2Cell.FromFacePosLevel(0, 0, 0), 3);
            TestCellOps(r5, S2Cell.FromFacePosLevel(0, 0, 1), 3);
            TestCellOps(r5, S2Cell.FromFacePosLevel(1, 0, 1), 0);

            // Rectangle consisting of a single point.
            TestCellOps(new S2R2Rect(new R2Point(0.51, 0.51), new R2Point(0.51, 0.51)),
                        S2Cell.FromFace(0), 3);

            // Rectangle that intersects the bounding rectangle of face 0
            // but not the face itself.
            TestCellOps(new S2R2Rect(new R2Point(0.01, 1.001), new R2Point(0.02, 1.002)),
                        S2Cell.FromFace(0), 0);

            // Rectangle that intersects one corner of face 0.
            TestCellOps(new S2R2Rect(new R2Point(0.99, -0.01), new R2Point(1.01, 0.01)),
                        S2Cell.FromFacePosLevel(0, ~0UL >> S2CellId.kFaceBits, 5),
              3);
        }

        private static void TestIntervalOps(S2R2Rect x, S2R2Rect y, string expected_rexion, S2R2Rect expected_union, S2R2Rect expected_intersection)
        {
            // Test all of the interval operations on the given pair of intervals.
            // "expected_rexion" is a sequence of "T" and "F" characters corresponding
            // to the expected results of Contains(), InteriorContains(), Intersects(),
            // and InteriorIntersects() respectively.

            Assert.Equal(expected_rexion[0] == 'T', x.Contains(y));
            Assert.Equal(expected_rexion[1] == 'T', x.InteriorContains(y));
            Assert.Equal(expected_rexion[2] == 'T', x.Intersects(y));
            Assert.Equal(expected_rexion[3] == 'T', x.InteriorIntersects(y));

            Assert.Equal(x.Union(y) == x, x.Contains(y));
            Assert.Equal(!x.Intersection(y).IsEmpty, x.Intersects(y));

            Assert.Equal(expected_union, x.Union(y));
            Assert.Equal(expected_intersection, x.Intersection(y));

            if (y.GetSize() == new R2Point(0, 0))
            {
                S2R2Rect r = x;
                r.AddPoint(y.Lo);
                Assert.Equal(expected_union, r);
            }
        }

        private static void TestCellOps(S2R2Rect r, S2Cell cell, int level)
        {
            // Test the relationship between the given rectangle and cell:
            // 0 == no intersection, 2 == Intersects,
            // 3 == Intersects and one region contains a vertex of the other,
            // 4 == Contains

            bool vertex_contained = false;
            for (int i = 0; i < 4; ++i)
            {
                // This would be easier to do byructing an S2R2Rect from the cell,
                // but that would defeat the purpose of testing this code independently.
                if (S2Coords.FaceXYZtoUV(0, cell.GetVertexRaw(i), out var u, out var v))
                {
                    if (r.Contains(new R2Point(S2Coords.UVtoST(u), S2Coords.UVtoST(v))))
                        vertex_contained = true;
                }
                if (!r.IsEmpty && cell.Contains(S2R2Rect.ToS2Point(r.GetVertex(i))))
                    vertex_contained = true;
            }
            Assert.Equal(level >= 2, r.MayIntersect(cell));
            Assert.Equal(level >= 3, vertex_contained);
            Assert.Equal(level >= 4, r.Contains(cell));
        }
    }
}
