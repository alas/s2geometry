namespace S2Geometry;

public class S2CoordsTests
{
    [Fact]
    internal void Test_TraversalOrder()
    {
        for (int r = 0; r < 4; ++r)
        {
            for (int i = 0; i < 4; ++i)
            {
                // Check consistency with respect to swapping axes.
                Assert.Equal(S2.kIJtoPos[r][i],
                          S2.kIJtoPos[r ^ S2.kSwapMask][SwapAxes(i)]);
                Assert.Equal(S2.kPosToIJ[r][i],
                          SwapAxes(S2.kPosToIJ[r ^ S2.kSwapMask][i]));

                // Check consistency with respect to reversing axis directions.
                Assert.Equal(S2.kIJtoPos[r][i],
                          S2.kIJtoPos[r ^ S2.kInvertMask][InvertBits(i)]);
                Assert.Equal(S2.kPosToIJ[r][i],
                          InvertBits(S2.kPosToIJ[r ^ S2.kInvertMask][i]));

                // Check that the two tables are inverses of each other.
                Assert.Equal(S2.kIJtoPos[r][S2.kPosToIJ[r][i]], i);
                Assert.Equal(S2.kPosToIJ[r][S2.kIJtoPos[r][i]], i);
            }
        }
    }

    [Fact]
    internal void Test_ST_UV_Conversions()
    {
        // Check boundary conditions.
        for (double s = 0; s <= 1; s += 0.5)
        {
            /*volatile*/
            double u = S2.STtoUV(s);
            Assert.Equal(u, 2 * s - 1);
        }
        for (double u = -1; u <= 1; ++u)
        {
            /*volatile*/
            double s = S2.UVtoST(u);
            Assert.Equal(s, 0.5 * (u + 1));
        }
        // Check that UVtoST and STtoUV are inverses.
        for (double x = 0; x <= 1; x += 0.0001)
        {
            Assert2.Near(S2.UVtoST(S2.STtoUV(x)), x, S2.DoubleError);
            Assert2.Near(S2.STtoUV(S2.UVtoST(2 * x - 1)), 2 * x - 1, S2.DoubleError);
        }
    }

    [Fact]
    internal void Test_FaceUVtoXYZ()
    {
        // Check that each face appears exactly once.
        S2Point sum = S2Point.Empty;
        for (int face = 0; face < 6; ++face)
        {
            S2Point center = S2.FaceUVtoXYZ(face, 0, 0);
            Assert.Equal(S2.GetNorm(face), center);
            Assert.Equal(1, Math.Abs(center[center.LargestAbsComponent()]));
            sum += center.Fabs();
        }
        Assert.Equal(sum, new S2Point(2, 2, 2));

        // Check that each face has a right-handed coordinate system.
        for (int face = 0; face < 6; ++face)
        {
            Assert.Equal(1, S2.GetUAxis(face).CrossProd(S2.GetVAxis(face))
                      .DotProd(S2.FaceUVtoXYZ(face, 0, 0)));
        }

        // Check that the Hilbert curves on each face combine to form a
        // continuous curve over the entire cube.
        for (int face = 0; face < 6; ++face)
        {
            // The Hilbert curve on each face starts at (-1,-1) and terminates
            // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
            int sign = ((face & S2.kSwapMask) != 0) ? -1 : 1;
            Assert.Equal(S2.FaceUVtoXYZ(face, sign, -sign),
                      S2.FaceUVtoXYZ((face + 1) % 6, -1, -1));
        }
    }

    [Fact]
    internal void Test_FaceXYZtoUVW()
    {
        for (int face = 0; face < 6; ++face)
        {
            Assert.Equal(new S2Point(0, 0, 0), S2.FaceXYZtoUVW(face, S2Point.Empty));
            Assert.Equal(new S2Point(1, 0, 0), S2.FaceXYZtoUVW(face, S2.GetUAxis(face)));
            Assert.Equal(new S2Point(-1, 0, 0), S2.FaceXYZtoUVW(face, -S2.GetUAxis(face)));
            Assert.Equal(new S2Point(0, 1, 0), S2.FaceXYZtoUVW(face, S2.GetVAxis(face)));
            Assert.Equal(new S2Point(0, -1, 0), S2.FaceXYZtoUVW(face, -S2.GetVAxis(face)));
            Assert.Equal(new S2Point(0, 0, 1), S2.FaceXYZtoUVW(face, S2.GetNorm(face)));
            Assert.Equal(new S2Point(0, 0, -1), S2.FaceXYZtoUVW(face, -S2.GetNorm(face)));
        }
    }

    [Fact]
    internal void Test_XYZToFaceSiTi()
    {
        // Check the conversion of random cells to center points and back.
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            for (int i = 0; i < 1000; ++i)
            {
                S2CellId id = S2Testing.GetRandomCellId(level);
                int actual_level = S2.XYZtoFaceSiTi(id.ToPoint(), out var face, out var si, out var ti);
                Assert.Equal(level, actual_level);
                S2CellId actual_id =
                    S2CellId.FromFaceIJ(face, (int)(si / 2), (int)(ti / 2)).Parent(level);
                Assert.Equal(id, actual_id);

                // Now test a point near the cell center but not equal to it.
                S2Point p_moved = id.ToPoint() + new S2Point(1e-13, 1e-13, 1e-13);
                actual_level = S2.XYZtoFaceSiTi(p_moved, out var face_moved, out var si_moved, out var ti_moved);
                Assert.Equal(-1, actual_level);
                Assert.Equal(face, face_moved);
                Assert.Equal(si, si_moved);
                Assert.Equal(ti, ti_moved);

                // Finally, test some random (si,ti) values that may be at different
                // levels, or not at a valid level at all (for example, si == 0).
                int face_random = S2Testing.Random.Uniform(S2CellId.kNumFaces);
                uint si_random, ti_random;
                uint mask = ~0U << (S2.kMaxCellLevel - level);
                do
                {
                    si_random = S2Testing.Random.Rand32() & mask;
                    ti_random = S2Testing.Random.Rand32() & mask;
                } while (si_random > S2.kMaxSiTi || ti_random > S2.kMaxSiTi);
                S2Point p_random = S2.FaceSiTitoXYZ(face_random, si_random, ti_random);
                actual_level = S2.XYZtoFaceSiTi(p_random, out face, out si, out ti);
                if (face != face_random)
                {
                    // The chosen point is on the edge of a top-level face cell.
                    Assert.Equal(-1, actual_level);
                    Assert.True(si == 0 || si == S2.kMaxSiTi ||
                                ti == 0 || ti == S2.kMaxSiTi);
                }
                else
                {
                    Assert.Equal(si_random, si);
                    Assert.Equal(ti_random, ti);
                    if (actual_level >= 0)
                    {
                        Assert.Equal(p_random, S2CellId.FromFaceIJ(face, (int)(si / 2), (int)(ti / 2)).Parent(actual_level).ToPoint());
                    }
                }
            }
        }
    }

    [Fact]
    internal void Test_UVNorms()
    {
        // Check that GetUNorm and GetVNorm compute right-handed normals for
        // an edge in the increasing U or V direction.
        for (int face = 0; face < 6; ++face)
        {
            for (double x = -1; x <= 1; x += 1 / 1024.0)
            {
                Assert2.DoubleEqual(S2.FaceUVtoXYZ(face, x, -1)
                    .CrossProd(S2.FaceUVtoXYZ(face, x, 1))
                    .Angle(S2.GetUNorm(face, x)), 0);
                Assert2.DoubleEqual(S2.FaceUVtoXYZ(face, -1, x)
                    .CrossProd(S2.FaceUVtoXYZ(face, 1, x))
                    .Angle(S2.GetVNorm(face, x)), 0);
            }
        }
    }

    [Fact]
    internal void Test_UVWAxis()
    {
        for (int face = 0; face < 6; ++face)
        {
            // Check that axes are consistent with FaceUVtoXYZ.
            Assert.Equal(S2.FaceUVtoXYZ(face, 1, 0) - S2.FaceUVtoXYZ(face, 0, 0),
                      S2.GetUAxis(face));
            Assert.Equal(S2.FaceUVtoXYZ(face, 0, 1) - S2.FaceUVtoXYZ(face, 0, 0),
                      S2.GetVAxis(face));
            Assert.Equal(S2.FaceUVtoXYZ(face, 0, 0), S2.GetNorm(face));

            // Check that every face coordinate frame is right-handed.
            Assert.Equal(1, S2.GetUAxis(face).CrossProd(S2.GetVAxis(face))
                      .DotProd(S2.GetNorm(face)));

            // Check that GetUVWAxis is consistent with GetUAxis, GetVAxis, GetNorm.
            Assert.Equal(S2.GetUAxis(face), S2.GetUVWAxis(face, 0));
            Assert.Equal(S2.GetVAxis(face), S2.GetUVWAxis(face, 1));
            Assert.Equal(S2.GetNorm(face), S2.GetUVWAxis(face, 2));
        }
    }

    [Fact]
    internal void Test_UVWFace()
    {
        // Check that GetUVWFace is consistent with GetUVWAxis.
        for (int face = 0; face < 6; ++face)
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                Assert.Equal(S2.GetFace(-S2.GetUVWAxis(face, axis)),
                          S2.GetUVWFace(face, axis, 0));
                Assert.Equal(S2.GetFace(S2.GetUVWAxis(face, axis)),
                          S2.GetUVWFace(face, axis, 1));
            }
        }
    }

    private static int SwapAxes(int ij)
    {
        return ((ij >> 1) & 1) + ((ij & 1) << 1);
    }

    private static int InvertBits(int ij)
    {
        return ij ^ 3;
    }
}
