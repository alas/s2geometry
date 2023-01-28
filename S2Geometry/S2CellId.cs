// An S2CellId is a 64-bit unsigned integer that uniquely identifies a
// cell in the S2 cell decomposition.  It has the following format:
//
//   id = [face][face_pos]
//
//   face:     a 3-bit number (range 0..5) encoding the cube face.
//
//   face_pos: a 61-bit number encoding the position of the center of this
//             cell along the Hilbert curve over this face (see the Wiki
//             pages for details).
//
// Sequentially increasing cell ids follow a continuous space-filling curve
// over the entire sphere.  They have the following properties:
//
//  - The id of a cell at level k consists of a 3-bit face number followed
//    by k bit pairs that recursively select one of the four children of
//    each cell.  The next bit is always 1, and all other bits are 0.
//    Therefore, the level of a cell is determined by the position of its
//    lowest-numbered bit that is turned on (for a cell at level k, this
//    position is 2 * (kMaxLevel - k).)
//
//  - The id of a parent cell is at the midpoint of the range of ids spanned
//    by its children (or by its descendants at any level).
//
// Leaf cells are often used to represent points on the unit sphere, and
// this class provides methods for converting directly between these two
// representations.  For cells that represent 2D regions rather than
// discrete point, it is better to use the S2Cell class.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.

namespace S2Geometry;

using System.Globalization;

public readonly record struct S2CellId(UInt64 Id) : IComparable<S2CellId>, IDecoder<S2CellId>
{
    #region Fields, Constants

    // Returns an invalid cell id.
    public static readonly S2CellId None = new(0);
    // Returns an invalid cell id guaranteed to be larger than any
    // valid cell id.  Useful for creating indexes.
    public static readonly S2CellId Sentinel = new(~0UL);
    // Although only 60 bits are needed to represent the index of a leaf cell, the
    // extra position bit lets us encode each cell as its Hilbert curve position
    // at the cell center, which is halfway along the portion of the Hilbert curve
    // that fills that cell.
    public const int kFaceBits = 3;
    public const int kNumFaces = 6;
    public const int kMaxLevel =
        S2.kMaxCellLevel;  // Valid levels: 0..kMaxLevel
    public const int kPosBits = 2 * S2.kMaxCellLevel + 1;
    public const int kMaxSize = 1 << S2.kMaxCellLevel;
    private const int kLookupBits = 4;
    private const double kScale = 1.0 / kMaxSize;
    private const double kLimit = 1.0 + S2.DoubleEpsilon;
    // This is the offset required to wrap around from the beginning of the
    // Hilbert curve to the end or vice versa; see NextWrap and PrevWrap.
    // SWIG doesn't understand UInt64{}, so use static_cast.
    private const UInt64 kWrapOffset = (UInt64)kNumFaces << kPosBits;

    #endregion

    #region Constructors

    // Construct a leaf cell containing the given point "p".  Usually there is
    // exactly one such cell, but for points along the edge of a cell, any
    // adjacent cell may be (deterministically) chosen.  This is because
    // S2CellIds are considered to be closed sets.  The returned cell will
    // always contain the given point, i.e.
    //
    //   S2Cell(S2CellId(p)).Contains(p)
    //
    // is always true.  The point "p" does not need to be normalized.
    //
    // If instead you want every point to be contained by exactly one S2Cell,
    // you will need to convert the S2CellIds to S2Loops (which implement point
    // containment this way).
    public S2CellId(S2Point p) : this(FromS2Point(p).Id) { }

    // Construct a leaf cell containing the given normalized S2LatLng.
    public S2CellId(S2LatLng ll) : this(ll.ToPoint()) { }

    #endregion

    #region Factories

    // Return the cell corresponding to a given S2 cube face.
    public static S2CellId FromFace(int face)
    {
        return new S2CellId(((UInt64)face << kPosBits) + LowestOnBitForLevel(0));
    }

    // Return a cell given its face (range 0..5), Hilbert curve position within
    // that face (an unsigned integer with S2CellId.kPosBits bits), and level
    // (range 0..kMaxLevel).  The given position will be modified to correspond
    // to the Hilbert curve position at the center of the returned cell.  This
    // is a static function rather than a constructor in order to indicate what
    // the arguments represent.
    public static S2CellId FromFacePosLevel(int face, UInt64 pos, int level)
    {
        var cell = new S2CellId(((UInt64)face << kPosBits) + (pos | 1));
        return cell.Parent(level);
    }

    public static S2CellId FromS2Point(S2Point p)
    {
        int face = S2.XYZtoFaceUV(p, out double u, out double v);
        int i = S2.STtoIJ(S2.UVtoST(u));
        int j = S2.STtoIJ(S2.UVtoST(v));
        return FromFaceIJ(face, i, j);
    }

    #endregion

    #region S2CellId

    // Return the direction vector corresponding to the center of the given
    // cell.  The vector returned by ToPointRaw is not necessarily unit length.
    // This method returns the same result as S2Cell.GetCenter().
    //
    // The maximum directional error in ToPoint() (compared to the exact
    // mathematical result) is 1.5 * S2Constants.DoubleEpsilon radians, and the maximum length
    // error is 2 * S2Constants.DoubleEpsilon (the same as Normalize).
    public S2Point ToPoint() => ToPointRaw().Normalize();
    public S2Point ToPointRaw()
    {
        int face = CenterSiTi(out int si, out int ti);
        return S2.FaceSiTitoXYZ(face, (uint)si, (uint)ti);
    }

    // Return the center of the cell in (s,t) coordinates (see s2coords.h).
    public R2Point CenterST()
    {
        CenterSiTi(out int si, out int ti);
        return new R2Point(S2.SiTitoST((uint)si), S2.SiTitoST((uint)ti));
    }

    // Return the edge length of this cell in (s,t)-space.
    public double SizeST() => SizeST(Level());

    // Return the edge length in (s,t)-space of cells at the given level.
    public static double SizeST(int level) => S2.IJtoSTMin(SizeIJ(level));

    // Return the bounds of this cell in (s,t)-space.
    public R2Rect BoundST()
    {
        double size = SizeST();
        return R2Rect.FromCenterSize(CenterST(), new R2Point(size, size));
    }

    // Return the center of the cell in (u,v) coordinates (see s2coords.h).
    // Note that the center of the cell is defined as the point at which it is
    // recursively subdivided into four children; in general, it is not at the
    // midpoint of the (u,v) rectangle covered by the cell.
    public R2Point CenterUV()
    {
        CenterSiTi(out int si, out int ti);
        return new R2Point(
            S2.STtoUV(S2.SiTitoST((uint)si)),
            S2.STtoUV(S2.SiTitoST((uint)ti)));
    }

    // Return the bounds of this cell in (u,v)-space.
    public R2Rect BoundUV()
    {
        var ij = new int[2];
        ToFaceIJOrientation(out ij[0], out ij[1], out _, false);
        return IJLevelToBoundUV(ij, Level());
    }

    // Expand a rectangle in (u,v)-space so that it contains all points within
    // the given distance of the boundary, and return the smallest such
    // rectangle.  If the distance is negative, then instead shrink this
    // rectangle so that it excludes all points within the given absolute
    // distance of the boundary.
    //
    // Distances are measured *on the sphere*, not in (u,v)-space.  For example,
    // you can use this method to expand the (u,v)-bound of an S2CellId so that
    // it contains all points within 5km of the original cell.  You can then
    // test whether a point lies within the expanded bounds like this:
    //
    //   R2Point uv;
    //   if (S2Coords.FaceXYZtoUV(face, point, &uv) && bound.Contains(uv)) { ... }
    //
    // Limitations:
    //
    //  - Because the rectangle is drawn on one of the six cube-face planes
    //    (i.e., {x,y,z} = +/-1), it can cover at most one hemisphere.  This
    //    limits the maximum amount that a rectangle can be expanded.  For
    //    example, S2CellId bounds can be expanded safely by at most 45 degrees
    //    (about 5000 km on the Earth's surface).
    //
    //  - The implementation is not exact for negative distances.  The resulting
    //    rectangle will exclude all points within the given distance of the
    //    boundary but may be slightly smaller than necessary.
    public static R2Rect ExpandedByDistanceUV(R2Rect uv, S1Angle distance)
    {
        // Expand each of the four sides of the rectangle just enough to include all
        // points within the given distance of that side.  (The rectangle may be
        // expanded by a different amount in (u,v)-space on each side.)
        double u0 = uv[0][0], u1 = uv[0][1], v0 = uv[1][0], v1 = uv[1][1];
        double max_u = Math.Max(Math.Abs(u0), Math.Abs(u1));
        double max_v = Math.Max(Math.Abs(v0), Math.Abs(v1));
        double sin_dist = Math.Sin(distance.Radians);
        var ep1 = ExpandEndpoint(u0, max_v, -sin_dist);
        var ep2 = ExpandEndpoint(u1, max_v, sin_dist);
        var ep3 = ExpandEndpoint(v0, max_u, -sin_dist);
        var ep4 = ExpandEndpoint(v1, max_u, sin_dist);
        return new R2Rect(new R1Interval(ep1, ep2), new R1Interval(ep3, ep4));
    }

    // Return the (face, si, ti) coordinates of the center of the cell.  Note
    // that although (si,ti) coordinates span the range [0,2**31] in general,
    // the cell center coordinates are always in the range [1,2**31-1] and
    // therefore can be represented using a signed 32-bit integer.
    public int CenterSiTi(out int psi, out int pti)
    {
        // First we compute the discrete (i,j) coordinates of a leaf cell contained
        // within the given cell.  Given that cells are represented by the Hilbert
        // curve position corresponding at their center, it turns out that the cell
        // returned by ToFaceIJOrientation is always one of two leaf cells closest
        // to the center of the cell (unless the given cell is a leaf cell itself,
        // in which case there is only one possibility).
        //
        // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
        // jmin) be the coordinates of its lower left-hand corner, the leaf cell
        // returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
        // (imin + s/2 - 1, jmin + s/2 - 1).  The first case is the one we want.
        // We can distinguish these two cases by looking at the low bit of "i" or
        // "j".  In the second case the low bit is one, unless s == 2 (i.e. the
        // level just above leaf cells) in which case the low bit is zero.
        //
        // In the code below, the expression ((i ^ (int(Id) >> 2)) & 1) is true
        // if we are in the second case described above.
        int face = ToFaceIJOrientation(out int i, out int j, out _, false);
        int delta = IsLeaf() ? 1 : (((i ^ ((int)Id >> 2)) & 1) != 0) ? 2 : 0;

        // Note that (2 * {i,j} + delta) will never overflow a 32-bit integer.
        psi = 2 * i + delta;
        pti = 2 * j + delta;
        return face;
    }

    // Return the S2LatLng corresponding to the center of the given cell.
    public S2LatLng ToLatLng()
    {
        return new S2LatLng(ToPointRaw());
    }

    // Return true if id() represents a valid cell.
    //
    // All methods require IsValid to be true unless otherwise specified
    // (although not all methods enforce this).
    public bool IsValid() => Face() < kNumFaces && (LowestOnBit() & 0x1555555555555555UL) != 0;

    // Which cube face this cell belongs to, in the range 0..5.
    public ulong Face() => Id >> kPosBits;

    // The position of the cell center along the Hilbert curve over this face,
    // in the range 0..(2**kPosBits-1).
    public ulong Pos() => Id & (~0UL >> kFaceBits);

    // Return the subdivision level of the cell (range 0..kMaxLevel).
    public int Level()
    {
        // We can't just Assert.True(IsValid) because we want Level to be
        // defined for end-iterators, i.e. S2CellId.End(kLevel).  However there is
        // no good way to define None.Level, so we do prohibit that.
        MyDebug.Assert(Id != 0);

        // A special case for leaf cells is not worthwhile.
        // Fast path for leaf cells.
        if (IsLeaf())
        {
            return S2.kMaxCellLevel;
        }
        var x = (uint)Id;
        var level = -1;
        if (x != 0)
        {
            level += 16;
        }
        else
        {
            x = (uint)(Id >> 32);
        }
        // We only need to look at even-numbered bits to determine the
        // level of a valid cell id.
        x &= (uint)-x; // Get lowest bit.
        if ((x & 0x00005555) != 0)
        {
            level += 8;
        }
        if ((x & 0x00550055) != 0)
        {
            level += 4;
        }
        if ((x & 0x05050505) != 0)
        {
            level += 2;
        }
        if ((x & 0x11111111) != 0)
        {
            level += 1;
        }
        MyDebug.Assert(level >= 0 && level <= S2.kMaxCellLevel);
        return level;
    }

    // Return the edge length of this cell in (i,j)-space.
    public int SizeIJ() => SizeIJ(Level());

    // Like the above, but return the size of cells at the given level.
    public static int SizeIJ(int level)
    {
        return 1 << (S2.kMaxCellLevel - level);
    }

    // Return true if this is a leaf cell (more efficient than checking
    // whether Level == kMaxLevel).
    public bool IsLeaf() => (Id & 1) != 0;

    // Return true if this is a top-level face cell (more efficient than
    // checking whether Level == 0).
    public bool IsFace() => (Id & (LowestOnBitForLevel(0) - 1)) == 0;

    // Return the child position (0..3) of this cell within its parent.
    // REQUIRES: Level >= 1.
    // No need for a special implementation; the compiler optimizes this well.
    public int ChildPosition() => ChildPosition(Level());

    // Return the child position (0..3) of this cell's ancestor at the given
    // level within its parent.  For example, child_position(1) returns the
    // position of this cell's level-1 ancestor within its top-level face cell.
    // REQUIRES: 1 <= level <= this.Level.
    public int ChildPosition(int level)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(level >= 1);
        MyDebug.Assert(level <= Level());
        return (int)(Id >> (2 * (S2.kMaxCellLevel - level) + 1)) & 3;
    }

    // These methods return the range of cell ids that are contained within this
    // cell (including itself).  The range is *inclusive* (i.e. test using >=
    // and <=) and the return values of both methods are valid leaf cell ids.
    // In other words, a.contains(b) if and only if
    //
    //     (b >= a.RangeMin && b <= a.RangeMax)
    //
    // If you want to iterate through all the descendants of this cell at a
    // particular level, use ChildBegin(level) and ChildEnd(level) instead.
    // Also see maximum_tile(), which can be used to iterate through a range of
    // cells using S2CellIds at different levels that are as large as possible.
    //
    // If you need to convert the range to a semi-open interval [min, limit)
    // (e.g., in order to use a key-value store that only supports semi-open
    // range queries), do not attempt to define "limit" as range_max.Next.
    // The problem is that leaf S2CellIds are 2 units apart, so the semi-open
    // interval [min, limit) includes an additional value (range_max.Id + 1)
    // which happens to be a valid S2CellId about one-third of the time and
    // is *never* contained by this cell.  (It always corresponds to a cell that
    // is larger than this one.)  You can define "limit" as (range_max.Id + 1)
    // if necessary (which is not always a valid S2CellId but can still be used
    // with FromToken/ToToken), or you can convert RangeMax to the key space
    // of your key-value store and define "limit" as Successor(key).
    //
    // Note that Sentinel().RangeMin == Sentinel.RangeMax == Sentinel().
    public S2CellId RangeMin() => new(Id - (LowestOnBit() - 1));

    public S2CellId RangeMax() => new(Id + (LowestOnBit() - 1));

    // Return true if the given cell is contained within this one.
    public bool Contains(S2CellId other)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(other.IsValid());
        return other >= RangeMin() && other <= RangeMax();
    }

    // Return true if the given cell intersects this one.
    public bool Intersects(S2CellId other)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(other.IsValid());
        return other.RangeMin() <= RangeMax() && other.RangeMax() >= RangeMin();
    }

    // Return the cell at the previous level or at the given level (which must
    // be less than or equal to the current level).
    public S2CellId Parent()
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(!IsFace());
        UInt64 new_lsb = LowestOnBit() << 2;
        return new S2CellId((Id & (~new_lsb + 1)) | new_lsb);
    }

    public S2CellId Parent(int level)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(level >= 0);
        MyDebug.Assert(level <= Level());
        UInt64 new_lsb = LowestOnBitForLevel(level);
        return new S2CellId((Id & (~new_lsb + 1)) | new_lsb);
    }

    // Return the immediate child of this cell at the given traversal order
    // position (in the range 0 to 3).  This cell must not be a leaf cell.
    public S2CellId Child(int position)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(!IsLeaf());
        // To change the level, we need to move the least-significant bit two
        // positions downward.  We do this by subtracting (4 * new_lsb) and adding
        // new_lsb.  Then to advance to the given child cell, we add
        // (2 * position * new_lsb).
        UInt64 new_lsb = LowestOnBit() >> 2;
        return new S2CellId(Id + (2 * (UInt64)position + 1 - 4) * new_lsb);
    }

    // Iterator-style methods for traversing the immediate children of a cell or
    // all of the children at a given level (greater than or equal to the current
    // level).  Note that the end value is exclusive, just like standard STL
    // iterators, and may not even be a valid cell id.  You should iterate using
    // code like this:
    //
    //   for(S2CellId c = id.ChildBegin(); c != id.ChildEnd(); c = c.Next)
    //     ...
    //
    // The convention for advancing the iterator is "c = c.Next" rather
    // than "++c" to avoid possible confusion with incrementing the
    // underlying 64-bit cell id.
    public S2CellId ChildBegin()
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(!IsLeaf());
        UInt64 old_lsb = LowestOnBit();
        return new S2CellId(Id - old_lsb + (old_lsb >> 2));
    }

    public S2CellId ChildBegin(int level)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(level >= Level());
        MyDebug.Assert(level <= S2.kMaxCellLevel);
        return new S2CellId(Id - LowestOnBit() + LowestOnBitForLevel(level));
    }

    public S2CellId ChildEnd()
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(!IsLeaf());
        UInt64 old_lsb = LowestOnBit();
        return new S2CellId(Id + old_lsb + (old_lsb >> 2));
    }

    public S2CellId ChildEnd(int level)
    {
        MyDebug.Assert(IsValid());
        MyDebug.Assert(level >= Level());
        MyDebug.Assert(level <= S2.kMaxCellLevel);
        return new(Id + LowestOnBit() + LowestOnBitForLevel(level));
    }

    // Return the next/previous cell at the same level along the Hilbert curve.
    // Works correctly when advancing from one face to the next, but
    // does *not* wrap around from the last face to the first or vice versa.
    public S2CellId Next() => new(Id + (LowestOnBit() << 1));

    public S2CellId Prev() => new(Id - (LowestOnBit() << 1));

    // This method advances or retreats the indicated number of steps along the
    // Hilbert curve at the current level, and returns the new position.  The
    // position is never advanced past End() or before Begin().
    public S2CellId Advance(Int64 steps)
    {
        if (steps == 0) return this;

        // We clamp the number of steps if necessary to ensure that we do not
        // advance past the End() or before the Begin() of this level.  Note that
        // min_steps and max_steps always fit in a signed 64-bit integer.

        int step_shift = 2 * (S2.kMaxCellLevel - Level()) + 1;
        if (steps < 0)
        {
            Int64 min_steps = -(Int64)(Id >> step_shift);
            if (steps < min_steps) steps = min_steps;
        }
        else
        {
            Int64 max_steps = (Int64)((kWrapOffset + LowestOnBit() - Id) >> step_shift);
            if (steps > max_steps) steps = max_steps;
        }
        // If steps is negative, then shifting it left has undefined behavior.
        // Cast to UInt64 for a 2's complement answer.
        return new S2CellId(Id + ((UInt64)(steps) << step_shift));
    }

    // Returns the number of steps that this cell is from Begin(Level). The
    // return value is always non-negative.
    public Int64 DistanceFromBegin()
    {
        int step_shift = 2 * (S2.kMaxCellLevel - Level()) + 1;
        return (Int64)(Id >> step_shift);
    }

    // Like Next and Prev, but these methods wrap around from the last face
    // to the first and vice versa.  They should *not* be used for iteration in
    // conjunction with ChildBegin(), ChildEnd(), Begin(), or End().  The
    // input must be a valid cell id.
    public S2CellId NextWrap()
    {
        MyDebug.Assert(IsValid());
        S2CellId n = Next();
        if (n.Id < kWrapOffset) return n;
        return new S2CellId(n.Id - kWrapOffset);
    }

    public S2CellId PrevWrap()
    {
        MyDebug.Assert(IsValid());
        S2CellId p = Prev();
        if (p.Id < kWrapOffset) return p;
        return new S2CellId(p.Id + kWrapOffset);
    }

    // This method advances or retreats the indicated number of steps along the
    // Hilbert curve at the current level, and returns the new position.  The
    // position wraps between the first and last faces as necessary.  The input
    // must be a valid cell id.
    public S2CellId AdvanceWrap(Int64 steps)
    {
        MyDebug.Assert(IsValid());
        if (steps == 0) return this;

        int step_shift = 2 * (S2.kMaxCellLevel - Level()) + 1;
        if (steps < 0)
        {
            Int64 min_steps = -(Int64)(Id >> step_shift);
            if (steps < min_steps)
            {
                Int64 step_wrap = (Int64)(kWrapOffset >> step_shift);
                steps %= step_wrap;
                if (steps < min_steps) steps += step_wrap;
            }
        }
        else
        {
            // Unlike advance(), we don't want to return End(level).
            Int64 max_steps = (Int64)((kWrapOffset - Id) >> step_shift);
            if (steps > max_steps)
            {
                Int64 step_wrap = (Int64)(kWrapOffset >> step_shift);
                steps %= step_wrap;
                if (steps > max_steps) steps -= step_wrap;
            }
        }
        return new S2CellId(Id + ((UInt64)(steps) << step_shift));
    }

    // Return the largest cell with the same RangeMin and such that
    // RangeMax < limit.RangeMin.  Returns "limit" if no such cell exists.
    // This method can be used to generate a small set of S2CellIds that covers
    // a given range (a "tiling").  This example shows how to generate a tiling
    // for a semi-open range of leaf cells [start, limit):
    //
    //   for (S2CellId id = start.maximum_tile(limit);
    //        id != limit; id = id.Next.maximum_tile(limit)) { ... }
    //
    // Note that in general the cells in the tiling will be of different sizes;
    // they gradually get larger (near the middle of the range) and then
    // gradually get smaller (as "limit" is approached).
    public S2CellId MaximumTile(S2CellId limit)
    {
        S2CellId id = this;
        S2CellId start = id.RangeMin();
        if (start >= limit.RangeMin()) return limit;

        if (id.RangeMax() >= limit)
        {
            // The cell is too large.  Shrink it.  Note that when generating coverings
            // of S2CellId ranges, this loop usually executes only once.  Also because
            // id.RangeMin < limit.RangeMin, we will always exit the loop by the
            // time we reach a leaf cell.
            do { id = id.Child(0); } while (id.RangeMax() >= limit);
            return id;
        }
        // The cell may be too small.  Grow it if necessary.  Note that generally
        // this loop only iterates once.
        while (!id.IsFace())
        {
            S2CellId parent = id.Parent();
            if (parent.RangeMin() != start || parent.RangeMax() >= limit) break;
            id = parent;
        }
        return id;
    }

    // Returns the level of the lowest common ancestor of this cell and "other",
    // i.e. the maximum level where this->parent(level) == other.parent(level).
    // Note that this definition also covers the situation where this cell is a
    // descendant of "other" or vice versa, or the two cells are the same,
    // since this->parent(this->level()) == *this.
    //
    // Returns -1 if the two cells do not have any common ancestor (i.e., they
    // are from different faces).
    public int CommonAncestorLevel(S2CellId other)
    {
        // Basically we find the first bit position at which the two S2CellIds
        // differ and convert that to a level.  The Math.Max() below is necessary for the
        // case where one S2CellId is a descendant of the other.
        UInt64 bits = Math.Max(Id ^ other.Id, Math.Max(LowestOnBit(), other.LowestOnBit()));
        MyDebug.Assert(bits != 0, "Because LowestOnBit is non-zero.");

        // Compute the position of the most significant bit, and then map the bit
        // position as follows:
        // {0} . 30, {1,2} . 29, {3,4} . 28, ... , {59,60} . 0, {61,62,63} . -1.
        return Math.Max(60 - BitsUtils.FindMSBSetNonZero64(bits), -1) >> 1;
    }

    // C++ Iterator-style methods for traversing all the cells along the Hilbert
    // curve at a given level (across all 6 faces of the cube).  Note that the
    // end value is exclusive (just like standard STL iterators), and is not a
    // valid cell id.
    public static S2CellId Begin(int level)
        => FromFace(0).ChildBegin(level);

    public static S2CellId End(int level)
        => FromFace(5).ChildEnd(level);

    // Methods to encode and decode cell ids to compact text strings suitable
    // for display or indexing.  Cells at lower levels (i.e. larger cells) are
    // encoded into fewer characters.  The maximum token length is 16.
    //
    // Tokens preserve ordering, i.e. ToToken(x) < ToToken(y) iff x < y.
    //
    // ToToken() returns a string by value for convenience; the compiler
    // does this without intermediate copying in most cases.
    //
    // These methods guarantee that FromToken(ToToken(x)) == x even when
    // "x" is an invalid cell id.  All tokens are alphanumeric strings.
    // FromToken() returns None for malformed inputs.
    public string ToToken()
    {
        // Simple implementation: print the id in hex without trailing zeros.
        // Using hex has the advantage that the tokens are case-insensitive, all
        // characters are alphanumeric, no characters require any special escaping
        // in queries for most indexing systems, and it's easy to compare cell
        // tokens against the feature ids of the corresponding features.
        //
        // Using base 64 would produce slightly shorter tokens, but for typical cell
        // sizes used during indexing (up to level 15 or so) the average savings
        // would be less than 2 bytes per cell which doesn't seem worth it.

        // "0" with trailing 0s stripped is the empty string, which is not a
        // reasonable token.  Encode as "X".
        if (Id == 0) return "X";

        return $"{Id:X16}".TrimEnd('0');
    }

    public static S2CellId FromToken(string token)
    {
        if (token.Length > 16) return None;

        var allHexa = token.All(t =>
        {
            return '0' <= t && t <= '9' || 'a' <= t && t <= 'f' || 'A' <= t && t <= 'F';
        });
        if (!allHexa) return None;

        token = (token + "0000000000000000")[..16];
        var success = UInt64.TryParse(token, NumberStyles.HexNumber, CultureInfo.InvariantCulture, out var id);
        if (!success) return None;

        return new S2CellId(id);
    }

    // Return the four cells that are adjacent across the cell's four edges.
    // Neighbors are returned in the order defined by S2Cell.GetEdge.  All
    // neighbors are guaranteed to be distinct.
    public void EdgeNeighbors(S2CellId[] neighbors)
    {
        int level = Level();
        int size = SizeIJ(level);
        int face = ToFaceIJOrientation(out int i, out int j, out _, false);

        // Edges 0, 1, 2, 3 are in the down, right, up, left directions.
        neighbors[0] = FromFaceIJSame(face, i, j - size, j - size >= 0)
            .Parent(level);
        neighbors[1] = FromFaceIJSame(face, i + size, j, i + size < kMaxSize)
            .Parent(level);
        neighbors[2] = FromFaceIJSame(face, i, j + size, j + size < kMaxSize)
            .Parent(level);
        neighbors[3] = FromFaceIJSame(face, i - size, j, i - size >= 0)
            .Parent(level);
    }

    // Return the neighbors of closest vertex to this cell at the given level,
    // by appending them to "output".  Normally there are four neighbors, but
    // the closest vertex may only have three neighbors if it is one of the 8
    // cube vertices.
    //
    // Requires: level < this.Level, so that we can determine which vertex is
    // closest (in particular, level == kMaxLevel is not allowed).
    public void AppendVertexNeighbors(int level, List<S2CellId> cellsId)
    {
        // "level" must be strictly less than this cell's level so that we can
        // determine which vertex this cell is closest to.
        MyDebug.Assert(level < Level());
        int face = ToFaceIJOrientation(out int i, out int j, out _, false);

        // Determine the i- and j-offsets to the closest neighboring cell in each
        // direction.  This involves looking at the next bit of "i" and "j" to
        // determine which quadrant of this.parent(level) this cell lies in.
        int halfsize = SizeIJ(level + 1);
        int size = halfsize << 1;
        bool isame, jsame;
        int ioffset, joffset;
        if ((i & halfsize) != 0)
        {
            ioffset = size;
            isame = (i + size) < kMaxSize;
        }
        else
        {
            ioffset = -size;
            isame = (i - size) >= 0;
        }
        if ((j & halfsize) != 0)
        {
            joffset = size;
            jsame = (j + size) < kMaxSize;
        }
        else
        {
            joffset = -size;
            jsame = (j - size) >= 0;
        }

        cellsId.Add(Parent(level));
        cellsId.Add(FromFaceIJSame(face, i + ioffset, j, isame).Parent(level));
        cellsId.Add(FromFaceIJSame(face, i, j + joffset, jsame).Parent(level));
        // If i- and j- edge neighbors are *both* on a different face, then this
        // vertex only has three neighbors (it is one of the 8 cube vertices).
        if (isame || jsame)
        {
            cellsId.Add(FromFaceIJSame(face, i + ioffset, j + joffset, isame && jsame)
                .Parent(level));
        }
    }

    // Append all neighbors of this cell at the given level to "output".  Two
    // cells X and Y are neighbors if their boundaries intersect but their
    // interiors do not.  In particular, two cells that intersect at a single
    // point are neighbors.  Note that for cells adjacent to a face vertex, the
    // same neighbor may be appended more than once.
    //
    // REQUIRES: nbr_level >= this.Level.
    public void AppendAllNeighbors(int nbr_level, List<S2CellId> cellsId)
    {
        MyDebug.Assert(nbr_level >= Level());
        int face = ToFaceIJOrientation(out int i, out int j, out _, false);

        // Find the coordinates of the lower left-hand leaf cell.  We need to
        // normalize (i,j) to a known position within the cell because nbr_level
        // may be larger than this cell's level.
        int size = SizeIJ();
        i &= -size;
        j &= -size;

        int nbr_size = SizeIJ(nbr_level);
        MyDebug.Assert(nbr_size <= size);

        // We compute the top-bottom, left-right, and diagonal neighbors in one
        // pass.  The loop test is at the end of the loop to avoid 32-bit overflow.
        for (int k = -nbr_size; ; k += nbr_size)
        {
            bool same_face;
            if (k < 0)
            {
                same_face = (j + k >= 0);
            }
            else if (k >= size)
            {
                same_face = (j + k < kMaxSize);
            }
            else
            {
                same_face = true;
                // Top and bottom neighbors.
                cellsId.Add(FromFaceIJSame(face, i + k, j - nbr_size,
                                                 j - size >= 0).Parent(nbr_level));
                cellsId.Add(FromFaceIJSame(face, i + k, j + size,
                                                 j + size < kMaxSize).Parent(nbr_level));
            }
            // Left, right, and diagonal neighbors.
            cellsId.Add(FromFaceIJSame(face, i - nbr_size, j + k,
                                             same_face && i - size >= 0)
                              .Parent(nbr_level));
            cellsId.Add(FromFaceIJSame(face, i + size, j + k,
                                             same_face && i + size < kMaxSize)
                              .Parent(nbr_level));
            if (k >= size) break;
        }
    }

    /////////////////////////////////////////////////////////////////////
    // Low-level methods.

    // Return a leaf cell given its cube face (range 0..5) and
    // i- and j-coordinates (see s2coords.h).
    public static S2CellId FromFaceIJ(int face, int i, int j)
    {
        // Initialization if not done yet
        MaybeInit();

        // Optimization notes:
        //  - Non-overlapping bit fields can be combined with either "+" or "|".
        //    Generally "+" seems to produce better code, but not always.

        // Note that this value gets shifted one bit to the left at the end
        // of the function.
        UInt64 n = (UInt64)face << (kPosBits - 1);

        // Alternating faces have opposite Hilbert curve orientations; this
        // is necessary in order for all faces to have a right-handed
        // coordinate system.
        UInt64 bits = (UInt64)(face & S2.kSwapMask);

        // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
        // curve position.  The lookup table transforms a 10-bit key of the form
        // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
        // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
        // Hilbert curve orientation respectively.
        var i2 = (UInt64)i;
        var j2 = (UInt64)j;
        for (var k = 7; k >= 0; k--)
        {
            UInt64 mask = (1UL << kLookupBits) - 1UL;
            bits += ((i2 >> (k * kLookupBits)) & mask) << (kLookupBits + 2);
            bits += ((j2 >> (k * kLookupBits)) & mask) << 2;
            bits = (ulong)LookupPos[bits];
            n |= (bits >> 2) << (k * 2 * kLookupBits);
            bits &= (S2.kSwapMask | S2.kInvertMask);
        }

        return new S2CellId(n * 2 + 1);
    }

    // Return the (face, i, j) coordinates for the leaf cell corresponding to
    // this cell id.  Since cells are represented by the Hilbert curve position
    // at the center of the cell, the returned (i,j) for non-leaf cells will be
    // a leaf cell adjacent to the cell center.  If "orientation" is non-null,
    // also return the Hilbert curve orientation for the current cell.
    public int ToFaceIJOrientation(out int pi, out int pj, out int orientation, bool setOrientation)
    {
        // Initialization if not done yet
        MaybeInit();

        int i = 0, j = 0;
        var face = (int)Face();
        int bits = face & S2.kSwapMask;

        // Each iteration maps 8 bits of the Hilbert curve position into
        // 4 bits of "i" and "j".  The lookup table transforms a key of the
        // form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
        // letters [ijpo] represents bits of "i", "j", the Hilbert curve
        // position, and the Hilbert curve orientation respectively.
        //
        // On the first iteration we need to be careful to clear out the bits
        // representing the cube face.

        for (var k = 7; k >= 0; k--)
        {
            int nbits = (k == 7)
                ? (S2.kMaxCellLevel - 7 * kLookupBits)
                : kLookupBits;
            bits += ((int)(Id >> (k * 2 * kLookupBits + 1))
                        & ((1 << (2 * nbits)) - 1)) << 2;
            bits = LookupIj[bits];
            i += (bits >> (kLookupBits + 2)) << (k * kLookupBits);
            j += ((bits >> 2) & ((1 << kLookupBits) - 1)) << (k * kLookupBits);
            bits &= (S2.kSwapMask | S2.kInvertMask);
        }

        pi = i;
        pj = j;

        if (setOrientation)
        {
            // The position of a non-leaf cell at level "n" consists of a prefix of
            // 2*n bits that identifies the cell, followed by a suffix of
            // 2*(kMaxLevel-n)+1 bits of the form 10*.  If n==kMaxLevel, the suffix is
            // just "1" and has no effect.  Otherwise, it consists of "10", followed
            // by (kMaxLevel-n-1) repetitions of "00", followed by "0".  The "10" has
            // no effect, while each occurrence of "00" has the effect of reversing
            // the kSwapMask bit.
            MyDebug.Assert(0 == S2.kPosToOrientation[2]);
            MyDebug.Assert(S2.kSwapMask == S2.kPosToOrientation[0]);
            if ((LowestOnBit() & 0x1111111111111110UL) != 0)
            {
                bits ^= S2.kSwapMask;
            }
            orientation = bits;
        }
        else
        {
            orientation = 0;
        }
        return face;
    }

    // Return the lowest-numbered bit that is on for this cell id, which is
    // equal to (UInt64{1} << (2 * (kMaxLevel - level))).  So for example,
    // a.LowestOnBit <= b.LowestOnBit if and only if a.Level >= b.Level, but the
    // first test is more efficient.
    public ulong LowestOnBit() => Id & (~Id + 1);

    // Return the lowest-numbered bit that is on for cells at the given level.
    public static UInt64 LowestOnBitForLevel(int level) => 1UL << (2 * (S2.kMaxCellLevel - level));

    // Return the bound in (u,v)-space for the cell at the given level containing
    // the leaf cell with the given (i,j)-coordinates.
    public static R2Rect IJLevelToBoundUV(int[] ij, int level)
    {
        var r = new R1Interval[2];
        int cell_size = SizeIJ(level);
        for (int d = 0; d < 2; ++d)
        {
            int ij_lo = ij[d] & -cell_size;
            int ij_hi = ij_lo + cell_size;
            r[d] = new R1Interval(
                S2.STtoUV(S2.IJtoSTMin(ij_lo)),
                S2.STtoUV(S2.IJtoSTMin(ij_hi)));
        }
        return new R2Rect(r[0], r[1]);
    }

    // Given a face and a point (i,j) where either i or j is outside the valid
    // range [0..kMaxSize-1], this function first determines which neighboring
    // face "contains" (i,j), and then returns the leaf cell on that face which
    // is adjacent to the given face and whose distance from (i,j) is minimal.
    private static S2CellId FromFaceIJWrap(int face, int i, int j)
    {
        // Convert i and j to the coordinates of a leaf cell just beyond the
        // boundary of this face.  This prevents 32-bit overflow in the case
        // of finding the neighbors of a face cell.
        i = Math.Max(-1, Math.Min(kMaxSize, i));
        j = Math.Max(-1, Math.Min(kMaxSize, j));

        // We want to wrap these coordinates onto the appropriate adjacent face.
        // The easiest way to do this is to convert the (i,j) coordinates to (x,y,z)
        // (which yields a point outside the normal face boundary), and then call
        // S2Coords.XYZtoFaceUV() to project back onto the correct face.
        //
        // The code below converts (i,j) to (si,ti), and then (si,ti) to (u,v) using
        // the linear projection (u=2*s-1 and v=2*t-1).  (The code further below
        // converts back using the inverse projection, s=0.5*(u+1) and t=0.5*(v+1).
        // Any projection would work here, so we use the simplest.)  We also clamp
        // the (u,v) coordinates so that the point is barely outside the
        // [-1,1]x[-1,1] face rectangle, since otherwise the reprojection step
        // (which divides by the new z coordinate) might change the other
        // coordinates enough so that we end up in the wrong leaf cell.

        // The arithmetic below is designed to avoid 32-bit integer overflows.
        MyDebug.Assert(0 == kMaxSize % 2);
        double u = Math.Max(-kLimit, Math.Min(kLimit, kScale * (2 * (i - kMaxSize / 2) + 1)));
        double v = Math.Max(-kLimit, Math.Min(kLimit, kScale * (2 * (j - kMaxSize / 2) + 1)));

        // Find the leaf cell coordinates on the adjacent face, and convert
        // them to a cell id at the appropriate level.
        face = S2.XYZtoFaceUV(S2.FaceUVtoXYZ(face, u, v), out u, out v);
        return FromFaceIJ(face, S2.STtoIJ(0.5 * (u + 1)), S2.STtoIJ(0.5 * (v + 1)));
    }

    // Inline helper function that calls FromFaceIJ if "same_face" is true,
    // or FromFaceIJWrap if "same_face" is false.
    private static S2CellId FromFaceIJSame(int face, int i, int j, bool same_face) => same_face
            ? FromFaceIJ(face, i, j)
            : FromFaceIJWrap(face, i, j);

    // The following lookup tables are used to convert efficiently between an
    // (i,j) cell index and the corresponding position along the Hilbert curve.
    // "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
    // orientation of the current cell into 8 bits representing the order in which
    // that subcell is visited by the Hilbert curve, plus 2 bits indicating the
    // new orientation of the Hilbert curve within that subcell.  (Cell
    // orientations are represented as combination of kSwapMask and kInvertMask.)
    //
    // "lookup_ij" is an inverted table used for mapping in the opposite
    // direction.
    //
    // We also experimented with looking up 16 bits at a time (14 bits of position
    // plus 2 of orientation) but found that smaller lookup tables gave better
    // performance.  (2KB fits easily in the primary cache.)

    //private UInt16 lookup_pos = [1 << (2 * kLookupBits + 2)];
    //private UInt16 lookup_ij = [1 << (2 * kLookupBits + 2)];
    private static readonly int[] LookupPos = new int[1 << (2 * kLookupBits + 2)];
    private static readonly int[] LookupIj = new int[1 << (2 * kLookupBits + 2)];

    private static void InitLookupCell(int level, int i, int j, int orig_orientation, int pos, int orientation)
    {
        if (level == kLookupBits)
        {
            int ij = (i << kLookupBits) + j;
            LookupPos[(ij << 2) + orig_orientation] = (pos << 2) + orientation;
            LookupIj[(pos << 2) + orig_orientation] = (ij << 2) + orientation;
        }
        else
        {
            level++;
            i <<= 1;
            j <<= 1;
            pos <<= 2;
            var r = S2.kPosToIJ[orientation];
            InitLookupCell(level, i + (r[0] >> 1), j + (r[0] & 1), orig_orientation,
                           pos, orientation ^ S2.kPosToOrientation[0]);
            InitLookupCell(level, i + (r[1] >> 1), j + (r[1] & 1), orig_orientation,
                           pos + 1, orientation ^ S2.kPosToOrientation[1]);
            InitLookupCell(level, i + (r[2] >> 1), j + (r[2] & 1), orig_orientation,
                           pos + 2, orientation ^ S2.kPosToOrientation[2]);
            InitLookupCell(level, i + (r[3] >> 1), j + (r[3] & 1), orig_orientation,
                           pos + 3, orientation ^ S2.kPosToOrientation[3]);
        }
    }

    // This is a helper function for ExpandedByDistanceUV().
    //
    // Given an edge of the form (u,v0)-(u,v1), let max_v = Math.Max(abs(v0), abs(v1)).
    // This method returns a new u-coordinate u' such that the distance from the
    // line u=u' to the given edge (u,v0)-(u,v1) is exactly the given distance
    // (which is specified as the sine of the angle corresponding to the distance).
    private static double ExpandEndpoint(double u, double max_v, double sin_dist)
    {
        // This is based on solving a spherical right triangle, similar to the
        // calculation in S2Cap.GetRectBound.
        double sin_u_shift = sin_dist * Math.Sqrt((1 + u * u + max_v * max_v) / (1 + u * u));
        double cos_u_shift = Math.Sqrt(1 - sin_u_shift * sin_u_shift);
        // The following is an expansion of tan(atan(u) + asin(sin_u_shift)).
        return (cos_u_shift * u + sin_u_shift) / (cos_u_shift - sin_u_shift * u);
    }

    private static void MaybeInit()
    {
        if (!callOnceFlag)
        {
            InitLookupCell(0, 0, 0, 0, 0, 0);
            InitLookupCell(0, 0, 0, S2.kSwapMask, 0, S2.kSwapMask);
            InitLookupCell(0, 0, 0, S2.kInvertMask, 0, S2.kInvertMask);
            InitLookupCell(0, 0, 0, S2.kSwapMask | S2.kInvertMask, 0, S2.kSwapMask | S2.kInvertMask);
            callOnceFlag = true;
        }
    }
    private static bool callOnceFlag = false;

    #endregion

    #region IEncoder

    // Use encoder to generate a serialized representation of this cell id.
    // Can also encode an invalid cell.
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(sizeof(UInt64));  // A single UInt64.
        encoder.Put64(Id);
    }

    // Decodes an S2CellId encoded by Encode(). Returns true on success.
    public static (bool, S2CellId) Decode(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(UInt64)) return (false, default);
        var cell = new S2CellId(decoder.Get64());
        return (true, cell);
    }

    #endregion

    #region IComparable

    // When S2CellId is used as a key in one of the absl::btree container types,
    // indicate that linear rather than binary search should be used.  This is
    // much faster when the comparison function is cheap.
    //typedef std::true_type absl_btree_prefer_linear_node_search;

    public int CompareTo(S2CellId other)
    {
        if (Id < other.Id) return -1;
        if (Id > other.Id) return 1;
        return 0;
    }

    public static bool operator <(S2CellId x, S2CellId y) => x.Id < y.Id;
    public static bool operator >(S2CellId x, S2CellId y) => x.Id > y.Id;
    public static bool operator <=(S2CellId x, S2CellId y) => x.Id <= y.Id;
    public static bool operator >=(S2CellId x, S2CellId y) => x.Id >= y.Id;

    #endregion

    #region Object

    // Creates a human readable debug string.  Used for << and available for
    // direct usage as well.  The format is "f/dd..d" where "f" is a digit in
    // the range [0-5] representing the S2CellId face, and "dd..d" is a string
    // of digits in the range [0-3] representing each child's position with
    // respect to its parent.  (Note that the latter string may be empty.)
    //
    // For example "4/" represents S2CellId.FromFace(4), and "3/02" represents
    // S2CellId.FromFace(3).child(0).child(2).
    public override string ToString()
    {
        if (!IsValid())
        {
            return $"Invalid: {Id:X16}";
        }
        var output = new System.Text.StringBuilder();
        output.AppendFormat("{0}/", Face());
        for (int current_level = 1; current_level <= Level(); current_level++)
        {
            output.Append("0123"[ChildPosition(current_level)]);
        }
        return output.ToString();
    }

    #endregion
}
