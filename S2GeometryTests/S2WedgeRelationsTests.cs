using Xunit;

namespace S2Geometry
{
    /// <summary>
    /// For simplicity, all of these tests use an origin of (0, 0, 1).
    /// This shouldn't matter as long as the lower-level primitives are
    /// implemented correctly.
    /// </summary>
    public class S2WedgeRelationsTests
    {
        [Fact]
        public void Test_IntersectionInOneWedge() =>
            TestWedge(new S2Point(-1, 0, 10), new S2Point(0, 0, 1), new S2Point(1, 2, 10),
                      new S2Point(0, 1, 10), new S2Point(1, -2, 10),
                      false, true, WedgeRelation.WEDGE_PROPERLY_OVERLAPS);
        [Fact]
        public void Test_IntersectionInTwoWedges() =>
            TestWedge(new S2Point(-1, -1, 10), new S2Point(0, 0, 1), new S2Point(1, -1, 10),
                      new S2Point(1, 0, 10), new S2Point(-1, 1, 10),
                      false, true, WedgeRelation.WEDGE_PROPERLY_OVERLAPS);
        [Fact]
        public void Test_NormalContainment() =>
            TestWedge(new S2Point(-1, -1, 10), new S2Point(0, 0, 1), new S2Point(1, -1, 10),
                      new S2Point(-1, 0, 10), new S2Point(1, 0, 10),
                      true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
        [Fact]
        public void Test_ContainmentWithEqualityOnOneSide() =>
            TestWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(-1, -1, 10),
                      new S2Point(2, 1, 10), new S2Point(1, -5, 10),
                      true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
        [Fact]
        public void Test_ContainmentWithEqualityOnTheOtherSide() =>
            TestWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(-1, -1, 10),
                      new S2Point(1, -2, 10), new S2Point(-1, -1, 10),
                      true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
        [Fact]
        public void Test_ContainmentWithEqualityOnBothSides() =>
            TestWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(4, -5, 10),
                      new S2Point(-2, 3, 10), new S2Point(4, -5, 10),
                      true, true, WedgeRelation.WEDGE_EQUALS);
        [Fact]
        public void Test_DisjointWithEqualityOnOneSide() =>
            TestWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(4, -5, 10),
                      new S2Point(4, -5, 10), new S2Point(-2, -3, 10),
                      false, false, WedgeRelation.WEDGE_IS_DISJOINT);
        [Fact]
        public void Test_DisjointWithEqualityOnTheOtherSide() =>
            TestWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(0, 5, 10),
                      new S2Point(4, -5, 10), new S2Point(-2, 3, 10),
                      false, false, WedgeRelation.WEDGE_IS_DISJOINT);
        [Fact]
        public void Test_DisjointWithEqualityOnBothSides() =>
            TestWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(4, -5, 10),
                      new S2Point(4, -5, 10), new S2Point(-2, 3, 10),
                      false, false, WedgeRelation.WEDGE_IS_DISJOINT);
        [Fact]
        public void Test_BContainsAWithEqualityOnOneSide() =>
            TestWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(1, -5, 10),
                      new S2Point(2, 1, 10), new S2Point(-1, -1, 10),
                      false, true, WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
        [Fact]
        public void Test_BContainsAWithEequalityOnTheOtherSide() =>
            TestWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(1, -5, 10),
                      new S2Point(-2, 1, 10), new S2Point(1, -5, 10),
                      false, true, WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);

        private static void TestWedge(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2, bool contains, bool intersects, WedgeRelation wedgeRelation)
        {
            a0 = a0.Normalized;
            ab1 = ab1.Normalized;
            a2 = a2.Normalized;
            b0 = b0.Normalized;
            b2 = b2.Normalized;
            Assert.True(contains == S2WedgeRelations.WedgeContains(a0, ab1, a2, b0, b2));
            Assert.True(intersects == S2WedgeRelations.WedgeIntersects(a0, ab1, a2, b0, b2));
            Assert.True(wedgeRelation == S2WedgeRelations.GetWedgeRelation(a0, ab1, a2, b0, b2));
        }
    }
}
