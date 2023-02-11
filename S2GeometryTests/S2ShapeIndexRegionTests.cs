namespace S2Geometry;

public class S2ShapeIndexRegionTests
{
    private static S2CellId MakeCellId(string str) => S2CellIdUtils.FromDebugString(str);

    // Pad by at least twice the maximum error for reliable results.
    private static readonly double kPadding = 2 * (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kIntersectsRectErrorUVDist);

    private static S2Shape NewPaddedCell(S2CellId id, double padding_uv)
    {
        var ij = new int[2];
        var face = id.ToFaceIJOrientation(out ij[0], out ij[1], out _, false);
        var uv = S2CellId.IJLevelToBoundUV(ij, id.Level()).Expanded(padding_uv);
        var vertices = new S2Point[4];
        for (int i = 0; i < 4; ++i)
        {
            vertices[i] = S2.FaceUVtoXYZ(face, uv.GetVertex(i)).Normalize();
        }
        return new S2LaxLoopShape(vertices);
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_GetCapBound()
    {
        var id = S2CellIdUtils.FromDebugString("3/0123012301230123012301230123");

        // Add a polygon that is slightly smaller than the cell being tested.
        MutableS2ShapeIndex index = new();
        index.Add(NewPaddedCell(id, -kPadding));
        S2Cap cell_bound = new S2Cell(id).GetCapBound();
        S2Cap index_bound = index.MakeS2ShapeIndexRegion().GetCapBound();
        Assert.True(index_bound.Contains(cell_bound));

        // Note that S2CellUnion.GetCapBound returns a slightly larger bound than
        // S2Cell.GetBound even when the cell union consists of a single S2CellId.
        Assert.True(index_bound.RadiusAngle() <= 1.00001 * cell_bound.RadiusAngle());
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_GetRectBound()
    {
        var id = S2CellIdUtils.FromDebugString("3/0123012301230123012301230123");

        // Add a polygon that is slightly smaller than the cell being tested.
        MutableS2ShapeIndex index = new();
        index.Add(NewPaddedCell(id, -kPadding));
        S2LatLngRect cell_bound = new S2Cell(id).GetRectBound();
        S2LatLngRect index_bound = index.MakeS2ShapeIndexRegion().GetRectBound();
        Assert.Equal(index_bound, cell_bound);
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_GetCellUnionBoundMultipleFaces()
    {
        var ids = new List<S2CellId> { MakeCellId("3/00123"), MakeCellId("2/11200013") };
        MutableS2ShapeIndex index = new();
        foreach (var id in ids) index.Add(NewPaddedCell(id, -kPadding));
        var covering = new List<S2CellId>();
        index.MakeS2ShapeIndexRegion().GetCellUnionBound(covering);
        ids.Sort();
        Assert.Equal(ids, covering);
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_GetCellUnionBoundOneFace()
    {
        // This tests consists of 3 pairs of S2CellIds.  Each pair is located within
        // one of the children of face 5, namely the cells 5/0, 5/1, and 5/3.
        // We expect GetCellUnionBound to compute the smallest cell that bounds the
        // pair on each face.
        S2CellId[] input =
        {
            MakeCellId("5/010"), MakeCellId("5/0211030"),
            MakeCellId("5/110230123"), MakeCellId("5/11023021133"),
            MakeCellId("5/311020003003030303"), MakeCellId("5/311020023"),
        };
        S2CellId[] expected =
        {
            MakeCellId("5/0"), MakeCellId("5/110230"), MakeCellId("5/3110200")
        };
        MutableS2ShapeIndex index = new();
        foreach (var id in input)
        {
            // Add each shape 3 times to ensure that the S2ShapeIndex subdivides.
            for (int copy = 0; copy < 3; ++copy)
            {
                index.Add(NewPaddedCell(id, -kPadding));
            }
        }
        List<S2CellId> actual = new();
        index.MakeS2ShapeIndexRegion().GetCellUnionBound(actual);
        Assert.Equal(expected, actual);
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_ContainsCellMultipleShapes()
    {
        Assert.True(false); //TODO

        var id = S2CellIdUtils.FromDebugString("3/0123012301230123012301230123");

        // Add a polygon that is slightly smaller than the cell being tested.
        MutableS2ShapeIndex index = new();
        index.Add(NewPaddedCell(id, -kPadding));
        Assert.False(index.MakeS2ShapeIndexRegion().Contains(new S2Cell(id)));

        // Add a second polygon that is slightly larger than the cell being tested.
        // Note that Contains() should return true if *any* shape contains the cell.
        index.Add(NewPaddedCell(id, kPadding));
        Assert.True(index.MakeS2ShapeIndexRegion().Contains(new S2Cell(id)));

        // Verify that all children of the cell are also contained.
        for (var child = id.ChildBegin(); child != id.ChildEnd(); child = child.Next())
        {
            Assert.True(index.MakeS2ShapeIndexRegion().Contains(new S2Cell(child)));
        }
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_IntersectsShrunkenCell()
    {
        var target = S2CellIdUtils.FromDebugString("3/0123012301230123012301230123");

        // Add a polygon that is slightly smaller than the cell being tested.
        MutableS2ShapeIndex index = new();
        index.Add(NewPaddedCell(target, -kPadding));
        var region = index.MakeS2ShapeIndexRegion();

        // Check that the index intersects the cell itself, but not any of the
        // neighboring cells.
        Assert.True(region.MayIntersect(new S2Cell(target)));
        var nbrs = new List<S2CellId>();
        target.AppendAllNeighbors(target.Level(), nbrs);
        foreach (var id in nbrs)
        {
            Assert.False(region.MayIntersect(new S2Cell(id)));
        }
    }

    [Fact]
    internal void Test_S2ShapeIndexRegion_IntersectsExactCell()
    {
        var target = S2CellIdUtils.FromDebugString("3/0123012301230123012301230123");

        // Adds a polygon that exactly follows a cell boundary.
        MutableS2ShapeIndex index = new();
        index.Add(NewPaddedCell(target, 0.0));
        var region = index.MakeS2ShapeIndexRegion();

        // Check that the index intersects the cell and all of its neighbors.
        List<S2CellId> ids = new(){ target };
        target.AppendAllNeighbors(target.Level(), ids);
        foreach (S2CellId id in ids)
        {
            Assert.True(region.MayIntersect(new S2Cell(id)));
        }
    }

    // Tests that VisitIntersectingShapes() produces results that are consistent
    // with MayIntersect() and Contains() for the given S2ShapeIndex.  It tests
    // all cells in the given index, all ancestors of those cells, and a randomly
    // chosen subset of descendants of those cells.
    internal class VisitIntersectingShapesTest
    {
        internal VisitIntersectingShapesTest(S2ShapeIndex index)
        {
            index_ = index;
            iter_ = new(index);
            region_ = new(index);
            // Create an S2ShapeIndex for each shape in the original index, so that we
            // can use MayIntersect() and Contains() to determine the status of
            // individual shapes.
            for (int s = 0; s < index_.NumShapeIds(); ++s)
            {
                var shape_index = new MutableS2ShapeIndex();
                shape_index.Add(new S2WrappedShape(index_.Shape(s)));
                shape_indexes_.Add(shape_index);
            }
        }

        internal void Run()
        {
            for (S2CellId id = S2CellId.Begin(0);
                 id != S2CellId.End(0); id = id.Next())
            {
                TestCell(new S2Cell(id));
            }
        }

        private void TestCell(S2Cell target) {
            // Indicates whether each shape that intersects "target" also contains it.
            Dictionary<int, bool> shape_contains = new();
            Assert.True(region_.VisitIntersectingShapes(
                target, (S2Shape shape, bool contains_target) => {
                    // Verify that each shape is visited at most once.
                    Assert.False(shape_contains.ContainsKey(shape.Id));
                    shape_contains[shape.Id] = contains_target;
                    return true;
                }));
            for (int s = 0; s < index_.NumShapeIds(); ++s) {
                var shape_region = shape_indexes_[s].MakeS2ShapeIndexRegion();
                if (!shape_region.MayIntersect(target)) {
                    Assert.False(shape_contains.ContainsKey(s));
                } else {
                    Assert.Equal(shape_contains[s], shape_region.Contains(target));
                }
            }
            var (cellRelation, pos) = index_.LocateCell(target.Id);
            iter_.SetPosition(pos);
            switch (cellRelation)
            {
                case S2CellRelation.DISJOINT:
                    return;

                case S2CellRelation.SUBDIVIDED:
                    {
                        S2Cell[] children = new S2Cell[4];
                        Assert.True(target.Subdivide(children));
                        foreach (var child in children) {
                            TestCell(child);
                        }
                        return;
                    }

                case S2CellRelation.INDEXED:
                    {
                        // We check a few random descendant cells by continuing randomly down
                        // one branch of the tree for a few levels.
                        if (target.IsLeaf() || S2Testing.Random.OneIn(3)) return;
                        TestCell(new S2Cell(target.Id.Child(S2Testing.Random.Uniform(4))));
                        return;
                    }
            }
        }

        private readonly S2ShapeIndex index_;
        private readonly S2ShapeIndex.Enumerator iter_;
        private readonly S2ShapeIndexRegion<S2ShapeIndex> region_;
        private readonly List<MutableS2ShapeIndex> shape_indexes_ = new();
    }

    [Fact]
    internal void Test_VisitIntersectingShapes_Points()
    {
        List<S2Point> vertices = new();
        for (int i = 0; i < 100; ++i)
        {
            vertices.Add(S2Testing.RandomPoint());
        }
        MutableS2ShapeIndex index = new();
        index.Add(new S2PointVectorShape(vertices.ToArray()));
        new VisitIntersectingShapesTest(index).Run();
    }

    [Fact]
    internal void Test_VisitIntersectingShapes_Polylines() {
        MutableS2ShapeIndex index = new();
        S2Cap center_cap=new(new S2Point(1, 0, 0), S1Angle.FromRadians(0.5));
        for (int i = 0; i < 50; ++i)
        {
            S2Point center = S2Testing.SamplePoint(center_cap);
            S2Point[] vertices;
            if (S2Testing.Random.OneIn(10))
            {
                vertices = new[]{ center, center};  // Try a few degenerate polylines.
            }
            else
            {
                vertices = S2Testing.MakeRegularPoints(
                    center, S1Angle.FromRadians(S2Testing.Random.RandDouble()),
                    S2Testing.Random.Uniform(20) + 3);
            }
            index.Add(new S2LaxPolylineShape(vertices));
        }
        new VisitIntersectingShapesTest(index).Run();
    }

    [Fact]
    internal void Test_VisitIntersectingShapes_Polygons() {
        MutableS2ShapeIndex index = new();
        S2Cap center_cap = new(new(1, 0, 0), S1Angle.FromRadians(0.5));
        S2Testing.Fractal fractal = new();
        for (int i = 0; i < 10; ++i)
        {
            fractal.SetLevelForApproxMaxEdges(3 * 64);
            S2Point center = S2Testing.SamplePoint(center_cap);
            index.Add(new S2Loop.Shape(
                fractal.MakeLoop(S2Testing.GetRandomFrameAt(center),
                S1Angle.FromRadians(S2Testing.Random.RandDouble()))));
        }
        // Also add a big polygon containing most of the polygons above to ensure
        // that we test containment of cells that are ancestors of index cells.
        index.Add(NewPaddedCell(S2CellId.FromFace(0), 0));
        new VisitIntersectingShapesTest(index).Run();
    }
}
