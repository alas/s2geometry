using System.Collections.Generic;
using Xunit;

namespace S2Geometry
{
    public class S2ShapeIndexRegionTests
    {
        // Pad by at least twice the maximum error for reliable results.
        private static readonly double kPadding = 2 * (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kIntersectsRectErrorUVDist);

        [Fact]
        public void Test_S2ShapeIndexRegion_GetCapBound() {
            var id = S2CellId.FromDebugString("3/0123012301230123012301230123");

            // Add a polygon that is slightly smaller than the cell being tested.
            MutableS2ShapeIndex index = new();
            index.Add(NewPaddedCell(id, -kPadding));
            S2Cap cell_bound = new S2Cell(id).GetCapBound();
            S2Cap index_bound = index.MakeS2ShapeIndexRegion().GetCapBound();
            Assert.True(index_bound.Contains(cell_bound));

            // Note that S2CellUnion.GetCapBound returns a slightly larger bound than
            // S2Cell.GetBound even when the cell union consists of a single S2CellId.
            Assert.True(index_bound.RadiusAngle <= 1.00001 * cell_bound.RadiusAngle);
        }

        [Fact]
        public void Test_S2ShapeIndexRegion_GetRectBound() {
            var id = S2CellId.FromDebugString("3/0123012301230123012301230123");

            // Add a polygon that is slightly smaller than the cell being tested.
            MutableS2ShapeIndex index = new();
            index.Add(NewPaddedCell(id, -kPadding));
            S2LatLngRect cell_bound = new S2Cell(id).GetRectBound();
            S2LatLngRect index_bound = index.MakeS2ShapeIndexRegion().GetRectBound();
            Assert.Equal(index_bound, cell_bound);
        }

        [Fact]
        public void Test_S2ShapeIndexRegion_GetCellUnionBoundMultipleFaces() {
            var ids = new List<S2CellId>{ MakeCellId("3/00123"), MakeCellId("2/11200013") };
            MutableS2ShapeIndex index = new MutableS2ShapeIndex();
            foreach (var id in ids) index.Add(NewPaddedCell(id, -kPadding));
            var covering = new List<S2CellId>();
            index.MakeS2ShapeIndexRegion().GetCellUnionBound(covering);
            ids.Sort();
            Assert.Equal(ids, covering);
        }

        [Fact]
        public void Test_S2ShapeIndexRegion_GetCellUnionBoundOneFace() {
            // This tests consists of 3 pairs of S2CellIds.  Each pair is located within
            // one of the children of face 5, namely the cells 5/0, 5/1, and 5/3.
            // We expect GetCellUnionBound to compute the smallest cell that bounds the
            // pair on each face.
            S2CellId[] input = {
    MakeCellId("5/010"), MakeCellId("5/0211030"),
    MakeCellId("5/110230123"), MakeCellId("5/11023021133"),
    MakeCellId("5/311020003003030303"), MakeCellId("5/311020023"),
  };
            S2CellId[] expected = {
    MakeCellId("5/0"), MakeCellId("5/110230"), MakeCellId("5/3110200")
  };
            MutableS2ShapeIndex index = new MutableS2ShapeIndex();
            foreach (var id in input) {
                // Add each shape 3 times to ensure that the S2ShapeIndex subdivides.
                for (int copy = 0; copy < 3; ++copy) {
                    index.Add(NewPaddedCell(id, -kPadding));
                }
            }
            var actual = new List<S2CellId>();
            index.MakeS2ShapeIndexRegion().GetCellUnionBound(actual);
            Assert.Equal(expected, actual);
        }

        [Fact]
        public void Test_S2ShapeIndexRegion_ContainsCellMultipleShapes()
        {
            Assert.True(false); //TODO

            var id = S2CellId.FromDebugString("3/0123012301230123012301230123");

            // Add a polygon that is slightly smaller than the cell being tested.
            MutableS2ShapeIndex index = new();
            index.Add(NewPaddedCell(id, -kPadding));
            Assert.False(index.MakeS2ShapeIndexRegion().Contains(new S2Cell(id)));

            // Add a second polygon that is slightly larger than the cell being tested.
            // Note that Contains() should return true if *any* shape contains the cell.
            index.Add(NewPaddedCell(id, kPadding));
            Assert.True(index.MakeS2ShapeIndexRegion().Contains(new S2Cell(id)));

            // Verify that all children of the cell are also contained.
            for (var child = id.ChildBegin(); child != id.ChildEnd(); child = child.Next) {
                Assert.True(index.MakeS2ShapeIndexRegion().Contains(new S2Cell(child)));
            }
        }

        [Fact]
        public void Test_S2ShapeIndexRegion_IntersectsShrunkenCell() {
            var target = S2CellId.FromDebugString("3/0123012301230123012301230123");

            // Add a polygon that is slightly smaller than the cell being tested.
            MutableS2ShapeIndex index = new();
            index.Add(NewPaddedCell(target, -kPadding));
            var region = index.MakeS2ShapeIndexRegion();

            // Check that the index intersects the cell itself, but not any of the
            // neighboring cells.
            Assert.True(region.MayIntersect(new S2Cell(target)));
            var nbrs = new List<S2CellId>();
            target.AppendAllNeighbors(target.Level, nbrs);
            foreach (var id in nbrs) {
                Assert.False(region.MayIntersect(new S2Cell(id)));
            }
        }

        [Fact]
        public void Test_S2ShapeIndexRegion_IntersectsExactCell() {
            var target = S2CellId.FromDebugString("3/0123012301230123012301230123");

            // Adds a polygon that exactly follows a cell boundary.
            MutableS2ShapeIndex index = new();
            index.Add(NewPaddedCell(target, 0.0));
            var region = index.MakeS2ShapeIndexRegion();

            // Check that the index intersects the cell and all of its neighbors.
            var ids = new List<S2CellId>{ target };
            target.AppendAllNeighbors(target.Level, ids);
            foreach (S2CellId id in ids) {
                Assert.True(region.MayIntersect(new S2Cell(id)));
            }
        }

        private static S2CellId MakeCellId(string str)
        {
            return S2CellId.FromDebugString(str);
        }

        private static S2Shape NewPaddedCell(S2CellId id, double padding_uv)
        {
            var ij = new int[2];
            var face = id.ToFaceIJOrientation(out ij[0], out ij[1], out _, false);
            var uv = S2CellId.IJLevelToBoundUV(ij, id.Level).Expanded(padding_uv);
            var vertices = new S2Point[4];
            for (int i = 0; i < 4; ++i)
            {
                vertices[i] = S2Coords.FaceUVtoXYZ(face, uv.GetVertex(i)).Normalized;
            }
            return new S2LaxLoopShape(vertices);
        }
    }
}
