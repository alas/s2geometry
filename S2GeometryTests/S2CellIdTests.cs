namespace S2Geometry;

public class S2CellIdTests
{
    private const int kMaxExpandLevel = 3;
    private const int kMaxWalkLevel = 8;
    private const double kCellSize = 1.0 / (1 << kMaxWalkLevel);

    [Fact]
    internal void Test_S2CellId_DefaultConstructor()
    {
        S2CellId id = new();
        Assert.Equal(0UL, id.Id);
        Assert.False(id.IsValid());
    }

    [Fact]
    internal void Test_S2CellId_S2CellIdHash()
    {
        Assert.Equal(
            GetCellId(0, 90).GetHashCode(),
            GetCellId(0, 90).GetHashCode());
    }

    [Fact]
    internal void Test_S2CellId_FaceDefinitions()
    {
        Assert.Equal(0UL, GetCellId(0, 0).Face());
        Assert.Equal(1UL, GetCellId(0, 90).Face());
        Assert.Equal(2UL, GetCellId(90, 0).Face());
        Assert.Equal(3UL, GetCellId(0, 180).Face());
        Assert.Equal(4UL, GetCellId(0, -90).Face());
        Assert.Equal(5UL, GetCellId(-90, 0).Face());
    }

    [Fact]
    internal void Test_S2CellId_FromFace()
    {
        for (int face = 0; face < 6; ++face)
        {
            Assert.Equal(S2CellId.FromFacePosLevel(face, 0, 0), S2CellId.FromFace(face));
        }
    }

    [Fact]
    internal void Test_S2CellId_ParentChildRelationships()
    {
        S2CellId id = S2CellId.FromFacePosLevel(3, 0x12345678, S2.kMaxCellLevel - 4);
        Assert.True(id.IsValid());
        Assert.Equal(3UL, id.Face());
        Assert.Equal(0x12345700UL, id.Pos());
        Assert.Equal(S2.kMaxCellLevel - 4, id.Level());
        Assert.False(id.IsLeaf());

        Assert.Equal(0x12345610UL, id.ChildBegin(id.Level() + 2).Pos());
        Assert.Equal(0x12345640UL, id.ChildBegin().Pos());
        Assert.Equal(0x12345400UL, id.Parent().Pos());
        Assert.Equal(0x12345000UL, id.Parent(id.Level() - 2).Pos());

        // Check ordering of children relative to parents.
        Assert.True(id.ChildBegin() < id);
        Assert.True(id.ChildEnd() > id);
        Assert.Equal(id.ChildEnd(), id.ChildBegin().Next().Next().Next().Next());
        Assert.Equal(id.RangeMin(), id.ChildBegin(S2.kMaxCellLevel));
        Assert.Equal(id.RangeMax().Next(), id.ChildEnd(S2.kMaxCellLevel));

        // Check that cells are represented by the position of their center
        // along the Hilbert curve.
        Assert.Equal(2 * id.Id, id.RangeMin().Id + id.RangeMax().Id);
    }

    [Fact]
    internal void Test_S2CellId_SentinelRangeMinMax()
    {
        Assert.Equal(S2CellId.Sentinel, S2CellId.Sentinel.RangeMin());
        Assert.Equal(S2CellId.Sentinel, S2CellId.Sentinel.RangeMax());
    }

    [Fact]
    internal void Test_S2CellId_CenterSiTi()
    {
        S2CellId id = S2CellId.FromFacePosLevel(3, 0x12345678,
                                                 S2.kMaxCellLevel);
        // Check that the (si, ti) coordinates of the center end in a
        // 1 followed by (30 - level) 0s.

        // Leaf level, 30.
        id.CenterSiTi(out var si, out var ti);
        Assert.Equal(1 << 0, si & 1);
        Assert.Equal(1 << 0, ti & 1);

        // Level 29.
        id.Parent(S2.kMaxCellLevel - 1).CenterSiTi(out si, out ti);
        Assert.Equal(1 << 1, si & 3);
        Assert.Equal(1 << 1, ti & 3);

        // Level 28.
        id.Parent(S2.kMaxCellLevel - 2).CenterSiTi(out si, out ti);
        Assert.Equal(1 << 2, si & 7);
        Assert.Equal(1 << 2, ti & 7);

        // Level 20.
        id.Parent(S2.kMaxCellLevel - 10).CenterSiTi(out si, out ti);
        Assert.Equal(1 << 10, si & ((1 << 11) - 1));
        Assert.Equal(1 << 10, ti & ((1 << 11) - 1));

        // Level 10.
        id.Parent(S2.kMaxCellLevel - 20).CenterSiTi(out si, out ti);
        Assert.Equal(1 << 20, si & ((1 << 21) - 1));
        Assert.Equal(1 << 20, ti & ((1 << 21) - 1));

        // Level 0.
        id.Parent(0).CenterSiTi(out si, out ti);
        Assert.Equal(1 << 30, si & ((1U << 31) - 1));
        Assert.Equal(1 << 30, ti & ((1U << 31) - 1));
    }

    [Fact]
    internal void Test_S2CellId_Wrapping()
    {
        // Check wrapping from beginning of Hilbert curve to end and vice versa.
        Assert.Equal(S2CellId.End(0).Prev(), S2CellId.Begin(0).PrevWrap());

        Assert.Equal(S2CellId.FromFacePosLevel(
                      5, ~0UL >> S2CellId.kFaceBits, S2.kMaxCellLevel),
                  S2CellId.Begin(S2.kMaxCellLevel).PrevWrap());
        Assert.Equal(S2CellId.FromFacePosLevel(
                      5, ~0UL >> S2CellId.kFaceBits, S2.kMaxCellLevel),
                  S2CellId.Begin(S2.kMaxCellLevel).AdvanceWrap(-1));

        Assert.Equal(S2CellId.Begin(4), S2CellId.End(4).Prev().NextWrap());
        Assert.Equal(S2CellId.Begin(4), S2CellId.End(4).Advance(-1).AdvanceWrap(1));

        Assert.Equal(S2CellId.FromFacePosLevel(0, 0, S2.kMaxCellLevel),
                  S2CellId.End(S2.kMaxCellLevel).Prev().NextWrap());
        Assert.Equal(S2CellId.FromFacePosLevel(0, 0, S2.kMaxCellLevel),
                  S2CellId.End(S2.kMaxCellLevel).Advance(-1).AdvanceWrap(1));
    }

    [Fact]
    internal void Test_S2CellId_Advance()
    {
        S2CellId id = S2CellId.FromFacePosLevel(3, 0x12345678,
                                                 S2.kMaxCellLevel - 4);
        // Check basic properties of advance().
        Assert.Equal(S2CellId.End(0), S2CellId.Begin(0).Advance(7));
        Assert.Equal(S2CellId.End(0), S2CellId.Begin(0).Advance(12));
        Assert.Equal(S2CellId.Begin(0), S2CellId.End(0).Advance(-7));
        Assert.Equal(S2CellId.Begin(0), S2CellId.End(0).Advance(-12000000));
        int num_level_5_cells = 6 << (2 * 5);
        Assert.Equal(S2CellId.End(5).Advance(500 - num_level_5_cells),
                  S2CellId.Begin(5).Advance(500));
        Assert.Equal(id.Next().ChildBegin(S2.kMaxCellLevel),
                  id.ChildBegin(S2.kMaxCellLevel).Advance(256));
        Assert.Equal(S2CellId.FromFacePosLevel(5, 0, S2.kMaxCellLevel),
                  S2CellId.FromFacePosLevel(1, 0, S2.kMaxCellLevel)
                  .Advance((long)(4UL << (2 * S2.kMaxCellLevel))));

        // Check basic properties of advance_wrap().
        Assert.Equal(S2CellId.FromFace(1), S2CellId.Begin(0).AdvanceWrap(7));
        Assert.Equal(S2CellId.Begin(0), S2CellId.Begin(0).AdvanceWrap(12));
        Assert.Equal(S2CellId.FromFace(4), S2CellId.FromFace(5).AdvanceWrap(-7));
        Assert.Equal(S2CellId.Begin(0), S2CellId.Begin(0).AdvanceWrap(-12000000));
        Assert.Equal(S2CellId.Begin(5).AdvanceWrap(6644),
                  S2CellId.Begin(5).AdvanceWrap(-11788));
        Assert.Equal(id.Next().ChildBegin(S2.kMaxCellLevel),
                  id.ChildBegin(S2.kMaxCellLevel).AdvanceWrap(256));
        Assert.Equal(S2CellId.FromFacePosLevel(1, 0, S2.kMaxCellLevel),
                  S2CellId.FromFacePosLevel(5, 0, S2.kMaxCellLevel)
                  .AdvanceWrap((long)(2UL << (2 * S2.kMaxCellLevel))));
    }

    [Fact]
    internal void Test_S2CellId_DistanceFromBegin()
    {
        Assert.Equal(6, S2CellId.End(0).DistanceFromBegin());
        Assert.Equal(6 * (1L << (2 * S2.kMaxCellLevel)),
                  S2CellId.End(S2.kMaxCellLevel).DistanceFromBegin());

        Assert.Equal(0, S2CellId.Begin(0).DistanceFromBegin());
        Assert.Equal(0, S2CellId.Begin(S2.kMaxCellLevel).DistanceFromBegin());

        S2CellId id = S2CellId.FromFacePosLevel(3, 0x12345678, S2.kMaxCellLevel - 4);
        Assert.Equal(id, S2CellId.Begin(id.Level()).Advance(id.DistanceFromBegin()));
    }

    [Fact]
    internal void Test_S2CellId_MaximumTile()
    {
        // This method is tested more thoroughly in s2cell_union_test.cc.
        for (int iter = 0; iter < 1000; ++iter)
        {
            S2CellId id = S2Testing.GetRandomCellId(10);

            // Check that "limit" is returned for tiles at or beyond "limit".
            Assert.Equal(id, id.MaximumTile(id));
            Assert.Equal(id, id.Child(0).MaximumTile(id));
            Assert.Equal(id, id.Child(1).MaximumTile(id));
            Assert.Equal(id, id.Next().MaximumTile(id));
            Assert.Equal(id.Child(0), id.MaximumTile(id.Child(0)));

            // Check that the tile size is increased when possible.
            Assert.Equal(id, id.Child(0).MaximumTile(id.Next()));
            Assert.Equal(id, id.Child(0).MaximumTile(id.Next().Child(0)));
            Assert.Equal(id, id.Child(0).MaximumTile(id.Next().Child(1).Child(0)));
            Assert.Equal(id, id.Child(0).Child(0).MaximumTile(id.Next()));
            Assert.Equal(id, id.Child(0).Child(0).Child(0).MaximumTile(id.Next()));

            // Check that the tile size is decreased when necessary.
            Assert.Equal(id.Child(0), id.MaximumTile(id.Child(0).Next()));
            Assert.Equal(id.Child(0), id.MaximumTile(id.Child(0).Next().Child(0)));
            Assert.Equal(id.Child(0), id.MaximumTile(id.Child(0).Next().Child(1)));
            Assert.Equal(id.Child(0).Child(0),
                      id.MaximumTile(id.Child(0).Child(0).Next()));
            Assert.Equal(id.Child(0).Child(0).Child(0),
                      id.MaximumTile(id.Child(0).Child(0).Child(0).Next()));

            // Check that the tile size is otherwise unchanged.
            Assert.Equal(id, id.MaximumTile(id.Next()));
            Assert.Equal(id, id.MaximumTile(id.Next().Child(0)));
            Assert.Equal(id, id.MaximumTile(id.Next().Child(1).Child(0)));
        }
    }

    [Fact]
    internal void Test_S2CellId_GetCommonAncestorLevel()
    {
        // Two identical cell ids.
        Assert.Equal(0, S2CellId.FromFace(0)
            .CommonAncestorLevel(S2CellId.FromFace(0)));
        Assert.Equal(30, S2CellId.FromFace(0).ChildBegin(30).CommonAncestorLevel(S2CellId.FromFace(0).ChildBegin(30)));

        // One cell id is a descendant of the other.
        Assert.Equal(0, S2CellId.FromFace(0).ChildBegin(30).CommonAncestorLevel(S2CellId.FromFace(0)));
        Assert.Equal(0, S2CellId.FromFace(5).CommonAncestorLevel(S2CellId.FromFace(5).ChildEnd(30).Prev()));

        // Two cells that have no common ancestor.
        Assert.Equal(-1, S2CellId.FromFace(0).CommonAncestorLevel(S2CellId.FromFace(5)));
        Assert.Equal(-1, S2CellId.FromFace(2).ChildBegin(30).CommonAncestorLevel(S2CellId.FromFace(3).ChildEnd(20)));

        // Two cells that have a common ancestor distinct from both of them.
        Assert.Equal(8, S2CellId.FromFace(5).ChildBegin(9).Next().ChildBegin(15).CommonAncestorLevel(
                      S2CellId.FromFace(5).ChildBegin(9).ChildBegin(20)));
        Assert.Equal(1, S2CellId.FromFace(0).ChildBegin(2).ChildBegin(30).
                  CommonAncestorLevel(
                      S2CellId.FromFace(0).ChildBegin(2).Next().ChildBegin(5)));
    }

    [Fact]
    internal void Test_S2CellId_Inverses()
    {
        // Check the conversion of random leaf cells to S2LatLngs and back.
        for (int i = 0; i < 200000; ++i)
        {
            S2CellId id = S2Testing.GetRandomCellId(S2.kMaxCellLevel);
            Assert.True(id.IsLeaf());
            Assert.Equal(S2.kMaxCellLevel, id.Level());
            S2LatLng center = id.ToLatLng();
            Assert.Equal(id.Id, new S2CellId(center).Id);
        }
    }

    [Fact]
    internal void Test_S2CellId_Tokens()
    {
        // Test random cell ids at all levels.
        for (int i = 0; i < 10000; ++i)
        {
            var id = S2Testing.GetRandomCellId();
            var token = id.ToToken();
            Assert.True(token.Length <= 16);
            Assert.Equal(id, S2CellId.FromToken(token));
            Assert.Equal(id, S2CellId.FromToken(token));
        }
        // Check that invalid cell ids can be encoded, and round-trip is
        // the identity operation.
        string token2 = S2CellId.None.ToToken();
        Assert.Equal(S2CellId.None, S2CellId.FromToken(token2));
        Assert.Equal(S2CellId.None, S2CellId.FromToken(token2));

        // Sentinel is invalid.
        token2 = S2CellId.Sentinel.ToToken();
        Assert.Equal(S2CellId.FromToken(token2), S2CellId.Sentinel);
        Assert.Equal(S2CellId.FromToken(token2),
                  S2CellId.Sentinel);

        // Check an invalid face.
        token2 = S2CellId.FromFace(7).ToToken();
        Assert.Equal(S2CellId.FromToken(token2), S2CellId.FromFace(7));
        Assert.Equal(S2CellId.FromToken(token2),
                  S2CellId.FromFace(7));

        // Check that supplying tokens with non-alphanumeric characters
        // returns S2CellId.None.
        Assert.Equal(S2CellId.None, S2CellId.FromToken("876b e99"));
        Assert.Equal(S2CellId.None, S2CellId.FromToken("876bee99\n"));
        Assert.Equal(S2CellId.None, S2CellId.FromToken("876[ee99"));
        Assert.Equal(S2CellId.None, S2CellId.FromToken(" 876bee99"));
    }

    [Fact]
    internal void Test_S2CellId_EncodeDecode()
    {
        S2CellId id = new(0x7837423);
        Encoder encoder = new();
        id.Encode(encoder);
        var decoder = encoder.Decoder();
        var (success, decoded_id) = S2CellId.Decode(decoder);
        Assert.True(success);
        Assert.Equal(id, decoded_id);
    }

    [Fact]
    internal void Test_S2CellId_EncodeDecodeNoneCell()
    {
        S2CellId none_id = S2CellId.None;
        Encoder encoder = new();
        none_id.Encode(encoder);
        var decoder = encoder.Decoder();
        var (success, decoded_id) = S2CellId.Decode(decoder);
        Assert.True(success);
        Assert.Equal(none_id, decoded_id);
    }

    [Fact]
    internal void Test_S2CellId_DecodeFailsWithTruncatedBuffer()
    {
        S2CellId id = new(0x7837423);
        Encoder encoder = new();
        id.Encode(encoder);

        // Truncate encoded buffer.
        var decoder = encoder.Decoder(2);
        var (success, _) = S2CellId.Decode(decoder);
        Assert.False(success);
    }

    [Fact]
    internal void Test_S2CellId_Containment()
    {
        // Test contains() and intersects().
        var parent_map = new Dictionary<S2CellId, S2CellId>();
        var cells = new List<S2CellId>();
        for (int face = 0; face < 6; ++face)
        {
            ExpandCell(S2CellId.FromFace(face), cells, parent_map);
        }
        foreach (var end_id in cells)
        {
            foreach (var begin_id in cells)
            {
                bool contained = true;
                for (var id = begin_id; id != end_id; id = parent_map[id])
                {
                    if (!parent_map.ContainsKey(id))
                    {
                        contained = false;
                        break;
                    }
                }
                Assert.Equal(contained, end_id.Contains(begin_id));
                Assert.Equal(contained,
                          begin_id >= end_id.RangeMin() &&
                          begin_id <= end_id.RangeMax());
                Assert.Equal(end_id.Intersects(begin_id),
                          end_id.Contains(begin_id) || begin_id.Contains(end_id));
            }
        }
    }

    [Fact]
    internal void Test_S2CellId_Continuity()
    {
        // Make sure that sequentially increasing cell ids form a continuous
        // path over the surface of the sphere, i.e. there are no
        // discontinuous jumps from one region to another.

        double max_dist = S2.kMaxEdge.GetValue(kMaxWalkLevel);
        S2CellId end = S2CellId.End(kMaxWalkLevel);
        S2CellId id = S2CellId.Begin(kMaxWalkLevel);
        for (; id != end; id = id.Next())
        {
            Assert.True(id.ToPointRaw().Angle(id.NextWrap().ToPointRaw()) <= max_dist);
            Assert.Equal(id.NextWrap(), id.AdvanceWrap(1));
            Assert.Equal(id, id.NextWrap().AdvanceWrap(-1));

            // Check that the ToPointRaw() returns the center of each cell
            // in (s,t) coordinates.
            S2.XYZtoFaceUV(id.ToPointRaw(), out var u, out var v);
            Assert2.Near(Math.IEEERemainder(S2.UVtoST(u), 0.5 * kCellSize), 0.0, S2.DoubleError);
            Assert2.Near(Math.IEEERemainder(S2.UVtoST(v), 0.5 * kCellSize), 0.0, S2.DoubleError);
        }
    }

    [Fact]
    internal void Test_S2CellId_Coverage()
    {
        // Make sure that random points on the sphere can be represented to the
        // expected level of accuracy, which in the worst case is Math.Sqrt(2/3) times
        // the maximum arc length between the points on the sphere associated with
        // adjacent values of "i" or "j".  (It is Math.Sqrt(2/3) rather than 1/2 because
        // the cells at the corners of each face are stretched -- they have 60 and
        // 120 degree angles.)

        double max_dist = 0.5 * S2.kMaxDiag.GetValue(S2.kMaxCellLevel);
        for (int i = 0; i < 1000000; ++i)
        {
            S2Point p = S2Testing.RandomPoint();
            S2Point q = new S2CellId(p).ToPointRaw();
            Assert.True(p.Angle(q) <= max_dist);
        }
    }

    [Fact]
    internal void Test_S2CellId_Neighbors()
    {
        // Check the edge neighbors of face 1.
        var out_faces = new[] { 5, 3, 2, 0 };
        var face_nbrs = new S2CellId[4];
        S2CellId.FromFace(1).EdgeNeighbors(face_nbrs);
        for (int i = 0; i < 4; ++i)
        {
            Assert.True(face_nbrs[i].IsFace());
            Assert.Equal(out_faces[i], (int)face_nbrs[i].Face());
        }

        // Check the edge neighbors of the corner cells at all levels.  This case is
        // trickier because it requires projecting onto adjacent faces.
        const int kMaxIJ = S2CellId.kMaxSize - 1;
        for (int level = 1; level <= S2.kMaxCellLevel; ++level)
        {
            S2CellId id2 = S2CellId.FromFaceIJ(1, 0, 0).Parent(level);
            var nbrs2 = new S2CellId[4];
            id2.EdgeNeighbors(nbrs2);
            // These neighbors were determined manually using the face and axis
            // relationships defined in s2coords.cc.
            int size_ij = S2CellId.SizeIJ(level);
            Assert.Equal(S2CellId.FromFaceIJ(5, kMaxIJ, kMaxIJ).Parent(level), nbrs2[0]);
            Assert.Equal(S2CellId.FromFaceIJ(1, size_ij, 0).Parent(level), nbrs2[1]);
            Assert.Equal(S2CellId.FromFaceIJ(1, 0, size_ij).Parent(level), nbrs2[2]);
            Assert.Equal(S2CellId.FromFaceIJ(0, kMaxIJ, 0).Parent(level), nbrs2[3]);
        }

        // Check the vertex neighbors of the center of face 2 at level 5.
        var nbrs = new List<S2CellId>();
        new S2CellId(new S2Point(0, 0, 1)).AppendVertexNeighbors(5, nbrs);
        nbrs.Sort();
        for (int i = 0; i < 4; ++i)
        {
            Assert.Equal(S2CellId.FromFaceIJ(
                2, (1 << 29) - ((i < 2) ? 1 : 0), (1 << 29) - ((i == 0 || i == 3) ? 1 : 0))
                .Parent(5), nbrs[i]);
        }
        nbrs.Clear();

        // Check the vertex neighbors of the corner of faces 0, 4, and 5.
        S2CellId id1 = S2CellId.FromFacePosLevel(0, 0, S2.kMaxCellLevel);
        id1.AppendVertexNeighbors(0, nbrs);
        nbrs.Sort();
        Assert.Equal(3, nbrs.Count);
        Assert.Equal(S2CellId.FromFace(0), nbrs[0]);
        Assert.Equal(S2CellId.FromFace(4), nbrs[1]);
        Assert.Equal(S2CellId.FromFace(5), nbrs[2]);

        // Check that AppendAllNeighbors produces results that are consistent
        // with AppendVertexNeighbors for a bunch of random cells.
        for (var i = 0; i < 1000; ++i)
        {
            S2CellId id2 = S2Testing.GetRandomCellId();
            if (id2.IsLeaf()) id2 = id2.Parent();

            // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell ids,
            // so it's not reasonable to use large values of "diff".
            int max_diff = Math.Min(5, S2.kMaxCellLevel - id2.Level() - 1);
            int level = id2.Level() + S2Testing.Random.Uniform(max_diff + 1);
            TestAllNeighbors(id2, level);
        }
    }

    [Fact]
    internal void Test_S2CellId_ExpandedByDistanceUV()
    {
        const double max_dist_degrees = 10;
        for (var iter = 0; iter < 100; ++iter)
        {
            S2CellId id = S2Testing.GetRandomCellId();
            double dist_degrees = S2Testing.Random.UniformDouble(-max_dist_degrees, max_dist_degrees);
            TestExpandedByDistanceUV(id, S1Angle.FromDegrees(dist_degrees));
        }
    }

    [Fact]
    internal void Test_S2CellId_ToString()
    {
        Assert.Equal("3/", S2CellId.FromFace(3).ToString());
        Assert.Equal("4/000000000000000000000000000000",
                  S2CellId.FromFace(4).RangeMin().ToString());
        Assert.Equal("Invalid: 0000000000000000", S2CellId.None.ToString());
    }

    [Fact]
    internal void Test_S2CellId_FromDebugString()
    {
        Assert.Equal(S2CellId.FromFace(3), S2CellId.FromDebugString("3/"));
        Assert.Equal(S2CellId.FromFace(0).Child(2).Child(1),
                  S2CellId.FromDebugString("0/21"));
        Assert.Equal(S2CellId.FromFace(4).RangeMin(),
                  S2CellId.FromDebugString("4/000000000000000000000000000000"));
        Assert.Equal(S2CellId.None,
                  S2CellId.FromDebugString("4/0000000000000000000000000000000"));
        Assert.Equal(S2CellId.None, S2CellId.FromDebugString(""));
        Assert.Equal(S2CellId.None, S2CellId.FromDebugString("7/"));
        Assert.Equal(S2CellId.None, S2CellId.FromDebugString(" /"));
        Assert.Equal(S2CellId.None, S2CellId.FromDebugString("3:0"));
        Assert.Equal(S2CellId.None, S2CellId.FromDebugString("3/ 12"));
        Assert.Equal(S2CellId.None, S2CellId.FromDebugString("3/1241"));
    }

    [Fact]
    internal void Test_S2CellId_OutputOperator()
    {
        S2CellId cell = new(0xbb04000000000000UL);
        Assert.Equal("5/31200", cell.ToString());
    }

    private static S2CellId GetCellId(double lat_degrees, double lng_degrees)
        => new(S2LatLng.FromDegrees(lat_degrees, lng_degrees));

    // Returns a random point on the boundary of the given rectangle.
    private static R2Point SampleBoundary(R2Rect rect)
    {
        var uv = new double[2];
        int d = S2Testing.Random.Uniform(2);
        uv[d] = S2Testing.Random.UniformDouble(rect[d][0], rect[d][1]);
        uv[1 - d] = S2Testing.Random.OneIn(2) ? rect[1 - d][0] : rect[1 - d][1];
        return R2Point.FromCoords(uv);
    }

    // Returns the closest point to "uv" on the boundary of "rect".
    private static R2Point ProjectToBoundary(R2Point uv, R2Rect rect)
    {
        double du0 = Math.Abs(uv[0] - rect[0][0]);
        double du1 = Math.Abs(uv[0] - rect[0][1]);
        double dv0 = Math.Abs(uv[1] - rect[1][0]);
        double dv1 = Math.Abs(uv[1] - rect[1][1]);
        double dmin = new[] { du0, du1, dv0, dv1 }.Min();
        if (du0 == dmin) return new R2Point(rect[0][0], rect[1].Project(uv[1]));
        if (du1 == dmin) return new R2Point(rect[0][1], rect[1].Project(uv[1]));
        if (dv0 == dmin) return new R2Point(rect[0].Project(uv[0]), rect[1][0]);
        Assert.Equal(dmin, dv1); // Bug in ProjectToBoundary
        return new R2Point(rect[0].Project(uv[0]), rect[1][1]);
    }

    private static void TestExpandedByDistanceUV(S2CellId id, S1Angle distance)
    {
        R2Rect bound = id.BoundUV();
        R2Rect expanded = S2CellId.ExpandedByDistanceUV(bound, distance);
        for (int iter = 0; iter < 100; ++iter)
        {
            // Choose a point on the boundary of the rectangle.
            int face = S2Testing.Random.Uniform(6);
            R2Point center_uv = SampleBoundary(bound);
            S2Point center = S2.FaceUVtoXYZ(face, center_uv).Normalize();

            // Now sample a point from a disc of radius (2 * distance).
            S2Point p = S2Testing.SamplePoint(new S2Cap(center, 2 * S1Angle.FromRadians(distance.Abs())));

            // Find the closest point on the boundary to the sampled point.
            if (!S2.FaceXYZtoUV(face, p, out var uv)) continue;
            R2Point closest_uv = ProjectToBoundary(uv, bound);
            S2Point closest = S2.FaceUVtoXYZ(face, closest_uv).Normalize();
            S1Angle actual_dist = new(p, closest);

            if (distance >= S1Angle.Zero)
            {
                // "expanded" should contain all points in the original bound, and also
                // all points within "distance" of the boundary.
                if (bound.Contains(uv) || actual_dist < distance)
                {
                    Assert.True(expanded.Contains(uv));
                }
            }
            else
            {
                // "expanded" should not contain any points within "distance" of the
                // original boundary.
                if (actual_dist < -distance)
                {
                    Assert.False(expanded.Contains(uv));
                }
            }
        }
    }

    private static void ExpandCell(S2CellId parent, List<S2CellId> cells, Dictionary<S2CellId, S2CellId> parent_map)
    {
        cells.Add(parent);
        if (parent.Level() == kMaxExpandLevel) return;
        int face = parent.ToFaceIJOrientation(out _, out _, out var orientation, true);
        Assert.Equal((int)parent.Face(), face);

        var child = parent.ChildBegin();
        for (int pos = 0; child != parent.ChildEnd(); child = child.Next(), ++pos)
        {
            parent_map[child] = parent;
            // Do some basic checks on the children.
            Assert.Equal(child, parent.Child(pos));
            Assert.Equal(pos, child.ChildPosition());
            // Test child_position(level) on all the child's ancestors.
            for (var ancestor = child; ancestor.Level() >= 1; ancestor = parent_map[ancestor])
            {
                Assert.Equal(child.ChildPosition(ancestor.Level()), ancestor.ChildPosition());
            }
            Assert.Equal(pos, child.ChildPosition(child.Level()));
            Assert.Equal(parent.Level() + 1, child.Level());
            Assert.False(child.IsLeaf());
            Assert.Equal(face, child.ToFaceIJOrientation(out _, out _, out var child_orientation, true));
            Assert.Equal(orientation ^ S2.kPosToOrientation[pos], child_orientation);
            ExpandCell(child, cells, parent_map);
        }
    }

    private static void TestAllNeighbors(S2CellId id, int level)
    {
        Assert.True(level >= id.Level());
        Assert.True(level < S2.kMaxCellLevel);

        // We compute AppendAllNeighbors, and then add in all the children of "id"
        // at the given level.  We then compare this against the result of finding
        // all the vertex neighbors of all the vertices of children of "id" at the
        // given level.  These should give the same result.
        var all = new List<S2CellId>();
        var expected = new List<S2CellId>();

        id.AppendAllNeighbors(level, all);
        S2CellId end = id.ChildEnd(level + 1);
        for (S2CellId c = id.ChildBegin(level + 1); c != end; c = c.Next())
        {
            all.Add(c.Parent());
            c.AppendVertexNeighbors(level, expected);
        }
        // Sort the results and eliminate duplicates.
        all = new SortedSet<S2CellId>(all).ToList();
        expected = new SortedSet<S2CellId>(expected).ToList();
        Assert.Equal(expected, all);
    }
}
