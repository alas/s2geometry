// S2PaddedCell represents an S2Cell whose (u,v)-range has been expanded on
// all sides by a given amount of "padding".  Unlike S2Cell, its methods and
// representation are optimized for clipping edges against S2Cell boundaries
// to determine which cells are intersected by a given set of edges.
//
// This class is intended to be copied by value as desired.

namespace S2Geometry;

public class S2PaddedCell
{
    public S2CellId Id { get; }
    public double Padding { get; }
    public int Level { get; }

    // Return the bound for this cell (including padding).
    public R2Rect Bound { get; }

    // The rectangle in (u,v)-space that belongs to all four padded children.
    //
    // Return the "middle" of the padded cell, defined as the rectangle that
    // belongs to all four children.
    //
    // Note that this method is *not* thread-safe, because the return value is
    // computed on demand and cached.  (It is expected that this class will be
    // mainly useful in the context of single-threaded recursive algorithms.)
    public R2Rect Middle
    {
        get
        {
            // We compute this field lazily because it is not needed the majority of the
            // time (i.e., for cells where the recursion terminates).
            if (!middle_.HasValue)
            {
                int ij_size = S2CellId.SizeIJ(Level);
                double u = S2.STtoUV(S2.SiTitoST((uint)(2 * ij_lo_[0] + ij_size)));
                double v = S2.STtoUV(S2.SiTitoST((uint)(2 * ij_lo_[1] + ij_size)));
                middle_ = new R2Rect(new R1Interval(u - Padding, u + Padding),
                                     new R1Interval(v - Padding, v + Padding));
            }
            return middle_.Value;
        }
        set { middle_ = value; }
    }
    private R2Rect? middle_;

    private readonly int[] ij_lo_;      // Minimum (i,j)-coordinates of this cell, before padding
    private readonly int orientation_;  // Hilbert curve orientation of this cell (see s2coords.h)

    #region Constructors

    // Construct an S2PaddedCell for the given cell id and padding.
    public S2PaddedCell(S2CellId id, double padding)
    {
        ij_lo_ = new int[2];
        Id = id;
        Padding = padding;
        if (Id.IsFace())
        {
            // Fast path for constructing a top-level face (the most common case).
            double limit = 1 + padding;
            Bound = new R2Rect(new R1Interval(-limit, limit),
                                new R1Interval(-limit, limit));
            middle_ = new R2Rect(new R1Interval(-padding, padding),
                                 new R1Interval(-padding, padding));
            ij_lo_[0] = ij_lo_[1] = 0;
            orientation_ = (int)(Id.Face() & 1);
            Level = 0;
        }
        else
        {
            var ij = new int[2];
            id.ToFaceIJOrientation(out ij[0], out ij[1], out orientation_, true);
            Level = id.Level();
            Bound = S2CellId.IJLevelToBoundUV(ij, Level).Expanded(padding);
            int ij_size = S2CellId.SizeIJ(Level);
            ij_lo_[0] = ij[0] & -ij_size;
            ij_lo_[1] = ij[1] & -ij_size;
        }
    }

    // Construct the child of "parent" with the given (i,j) index.  The four
    // child cells have indices of (0,0), (0,1), (1,0), (1,1), where the i and j
    // indices correspond to increasing u- and v-values respectively.
    public S2PaddedCell(S2PaddedCell parent, int i, int j)
    {
        ij_lo_ = new int[2];
        Padding = parent.Padding;
        Level = parent.Level + 1;
        // Compute the position and orientation of the child incrementally from the
        // orientation of the parent.
        int pos = S2.kIJtoPos[parent.orientation_][2 * i + j];
        Id = parent.Id.Child(pos);
        int ij_size = S2CellId.SizeIJ(Level);
        ij_lo_[0] = parent.ij_lo_[0] + i * ij_size;
        ij_lo_[1] = parent.ij_lo_[1] + j * ij_size;
        orientation_ = parent.orientation_ ^ S2.kPosToOrientation[pos];
        // For each child, one corner of the bound is taken directly from the parent
        // while the diagonally opposite corner is taken from middle().
        var middle = parent.Middle;
        var x = new double[] { parent.Bound[0][0], parent.Bound[0][1] };
        var y = new double[] { parent.Bound[1][0], parent.Bound[1][1] };
        x[1 - i] = middle[0][1 - i];
        y[1 - j] = middle[1][1 - j];
        Bound = new R2Rect(new R1Interval(x), new R1Interval(y));
    }

    #endregion

    // Return the (i,j) coordinates for the child cell at the given traversal
    // position.  The traversal position corresponds to the order in which child
    // cells are visited by the Hilbert curve.
    public void GetChildIJ(int pos, out int i, out int j)
    {
        int ij = S2.kPosToIJ[orientation_][pos];
        i = ij >> 1;
        j = ij & 1;
    }

    // Return the smallest cell that contains all descendants of this cell whose
    // bounds intersect "rect".  For algorithms that use recursive subdivision
    // to find the cells that intersect a particular object, this method can be
    // used to skip all the initial subdivision steps where only one child needs
    // to be expanded.
    //
    // Note that this method is not the same as returning the smallest cell that
    // contains the intersection of this cell with "rect".  Because of the
    // padding, even if one child completely contains "rect" it is still
    // possible that a neighboring child also intersects "rect".
    //
    // REQUIRES: bound().Intersects(rect)
    public S2CellId ShrinkToFit(R2Rect rect)
    {
        MyDebug.Assert(Bound.Intersects(rect));

        // Quick rejection test: if "rect" contains the center of this cell along
        // either axis, then no further shrinking is possible.
        int ij_size = S2CellId.SizeIJ(Level);
        if (Level == 0)
        {
            // Fast path (most calls to this function start with a face cell).
            if (rect[0].Contains(0) || rect[1].Contains(0)) return Id;
        }
        else
        {
            if (rect[0].Contains(S2.STtoUV(S2.SiTitoST((uint)(2 * ij_lo_[0] + ij_size)))) ||
                rect[1].Contains(S2.STtoUV(S2.SiTitoST((uint)(2 * ij_lo_[1] + ij_size)))))
            {
                return Id;
            }
        }
        // Otherwise we expand "rect" by the given padding() on all sides and find
        // the range of coordinates that it spans along the i- and j-axes.  We then
        // compute the highest bit position at which the min and max coordinates
        // differ.  This corresponds to the first cell level at which at least two
        // children intersect "rect".

        // Increase the padding to compensate for the error in S2Coords.UVtoST().
        // (The constant below is a provable upper bound on the additional error.)
        R2Rect padded = rect.Expanded(Padding + 1.5 * S2.DoubleEpsilon);
        var ij_min = new int[2];  // Min i- or j- coordinate spanned by "padded"
        var ij_xor = new int[2];  // XOR of the min and max i- or j-coordinates
        for (int d = 0; d < 2; ++d)
        {
            ij_min[d] = Math.Max(ij_lo_[d], S2.STtoIJ(S2.UVtoST(padded[d][0])));
            int ij_max = Math.Min(ij_lo_[d] + ij_size - 1,
                             S2.STtoIJ(S2.UVtoST(padded[d][1])));
            ij_xor[d] = ij_min[d] ^ ij_max;
        }
        // Compute the highest bit position where the two i- or j-endpoints differ,
        // and then choose the cell level that includes both of these endpoints.  So
        // if both pairs of endpoints are equal we choose kMaxLevel; if they differ
        // only at bit 0, we choose (kMaxLevel - 1), and so on.
        var level_msb = ((ij_xor[0] | ij_xor[1]) << 1) + 1;
        var level = S2.kMaxCellLevel - BitsUtils.FindMSBSetNonZero((uint)level_msb);
        if (level <= Level) return Id;
        return S2CellId.FromFaceIJ((int)Id.Face(), ij_min[0], ij_min[1]).Parent(level);
    }

    // Return the center of this cell.
    public S2Point GetCenter()
    {
        var ij_size = S2CellId.SizeIJ(Level);
        var si = (uint)(2 * ij_lo_[0] + ij_size);
        var ti = (uint)(2 * ij_lo_[1] + ij_size);
        return S2.FaceSiTitoXYZ((int)Id.Face(), si, ti).Normalize();
    }

    // Return the vertex where the S2 space-filling curve enters this cell.
    public S2Point GetEntryVertex()
    {
        // The curve enters at the (0,0) vertex unless the axis directions are
        // reversed, in which case it enters at the (1,1) vertex.
        uint i = (uint)ij_lo_[0];
        uint j = (uint)ij_lo_[1];
        if ((orientation_ & S2.kInvertMask) != 0)
        {
            int ij_size = S2CellId.SizeIJ(Level);
            i = (uint)(i + ij_size);
            j = (uint)(j + ij_size);
        }
        return S2.FaceSiTitoXYZ((int)Id.Face(), 2 * i, 2 * j).Normalize();
    }

    // Return the vertex where the S2 space-filling curve exits this cell.
    public S2Point GetExitVertex()
    {
        // The curve exits at the (1,0) vertex unless the axes are swapped or
        // inverted but not both, in which case it exits at the (0,1) vertex.
        uint i = (uint)ij_lo_[0];
        uint j = (uint)ij_lo_[1];
        int ij_size = S2CellId.SizeIJ(Level);
        if (orientation_ == 0 || orientation_ == S2.kSwapMask + S2.kInvertMask)
        {
            i = (uint)(i + ij_size);
        }
        else
        {
            j = (uint)(j + ij_size);
        }
        return S2.FaceSiTitoXYZ((int)Id.Face(), 2 * i, 2 * j).Normalize();
    }
}
