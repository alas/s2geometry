// This file contains documentation of the various coordinate systems used
// throughout the library.  Most importantly, S2 defines a framework for
// decomposing the unit sphere into a hierarchy of "cells".  Each cell is a
// quadrilateral bounded by four geodesics.  The top level of the hierarchy is
// obtained by projecting the six faces of a cube onto the unit sphere, and
// lower levels are obtained by subdividing each cell into four children
// recursively.  Cells are numbered such that sequentially increasing cells
// follow a continuous space-filling curve over the entire sphere.  The
// transformation is designed to make the cells at each level fairly uniform
// in size.
//
//
////////////////////////// S2Cell Decomposition /////////////////////////
//
// The following methods define the cube-to-sphere projection used by
// the S2Cell decomposition.
//
// In the process of converting a latitude-longitude pair to a 64-bit cell
// id, the following coordinate systems are used:
//
//  (id)
//    An S2CellId is a 64-bit encoding of a face and a Hilbert curve position
//    on that face.  The Hilbert curve position implicitly encodes both the
//    position of a cell and its subdivision level (see s2cell_id.h).
//
//  (face, i, j)
//    Leaf-cell coordinates.  "i" and "j" are integers in the range
//    [0,(2**30)-1] that identify a particular leaf cell on the given face.
//    The (i, j) coordinate system is right-handed on each face, and the
//    faces are oriented such that Hilbert curves connect continuously from
//    one face to the next.
//
//  (face, s, t)
//    Cell-space coordinates.  "s" and "t" are real numbers in the range
//    [0,1] that identify a point on the given face.  For example, the point
//    (s, t) = (0.5, 0.5) corresponds to the center of the top-level face
//    cell.  This point is also a vertex of exactly four cells at each
//    subdivision level greater than zero.
//
//  (face, si, ti)
//    Discrete cell-space coordinates.  These are obtained by multiplying
//    "s" and "t" by 2**31 and rounding to the nearest unsigned integer.
//    Discrete coordinates lie in the range [0,2**31].  This coordinate
//    system can represent the edge and center positions of all cells with
//    no loss of precision (including non-leaf cells).  In binary, each
//    coordinate of a level-k cell center ends with a 1 followed by
//    (30 - k) 0s.  The coordinates of its edges end with (at least)
//    (31 - k) 0s.
//
//  (face, u, v)
//    Cube-space coordinates in the range [-1,1].  To make the cells at each
//    level more uniform in size after they are projected onto the sphere,
//    we apply a nonlinear transformation of the form u=f(s), v=f(t).
//    The (u, v) coordinates after this transformation give the actual
//    coordinates on the cube face (modulo some 90 degree rotations) before
//    it is projected onto the unit sphere.
//
//  (face, u, v, w)
//    Per-face coordinate frame.  This is an extension of the (face, u, v)
//    cube-space coordinates that adds a third axis "w" in the direction of
//    the face normal.  It is always a right-handed 3D coordinate system.
//    Cube-space coordinates can be converted to this frame by setting w=1,
//    while (u,v,w) coordinates can be projected onto the cube face by
//    dividing by w, i.e. (face, u/w, v/w).
//
//  (x, y, z)
//    Direction vector (S2Point).  Direction vectors are not necessarily unit
//    length, and are often chosen to be points on the biunit cube
//    [-1,+1]x[-1,+1]x[-1,+1].  They can be be normalized to obtain the
//    corresponding point on the unit sphere.
//
//  (lat, lng)
//    Latitude and longitude (S2LatLng).  Latitudes must be between -90 and
//    90 degrees inclusive, and longitudes must be between -180 and 180
//    degrees inclusive.
//
// Note that the (i, j), (s, t), (si, ti), and (u, v) coordinate systems are
// right-handed on all six faces.

// The name "S2" is derived from the mathematical
// symbol for the two-dimensional unit sphere (note that the "2" refers to the
// dimension of the surface, not the space it is embedded in).

// We have implemented three different projections from cell-space (s,t) to
// cube-space (u,v): linear, quadratic, and tangent.  They have the following
// tradeoffs:
//
//   Linear - This is the fastest transformation, but also produces the least
//   uniform cell sizes.  Cell areas vary by a factor of about 5.2, with the
//   largest cells at the center of each face and the smallest cells in
//   the corners.
//
//   Tangent - Transforming the coordinates via atan() makes the cell sizes
//   more uniform.  The areas vary by a maximum ratio of 1.4 as opposed to a
//   maximum ratio of 5.2.  However, each call to atan() is about as expensive
//   as all of the other calculations combined when converting from points to
//   cell ids, i.e. it reduces performance by a factor of 3.
//
//   Quadratic - This is an approximation of the tangent projection that
//   is much faster and produces cells that are almost as uniform in size.
//   It is about 3 times faster than the tangent projection for converting
//   cell ids to points or vice versa.  Cell areas vary by a maximum ratio of
//   about 2.1.
//
// Here is a table comparing the cell uniformity using each projection.  "Area
// ratio" is the maximum ratio over all subdivision levels of the largest cell
// area to the smallest cell area at that level, "edge ratio" is the maximum
// ratio of the longest edge of any cell to the shortest edge of any cell at
// the same level, and "diag ratio" is the ratio of the longest diagonal of
// any cell to the shortest diagonal of any cell at the same level.  "ToPoint"
// and "FromPoint" are the times in microseconds required to convert cell ids
// to and from points (unit vectors) respectively.  "ToPointRaw" is the time
// to convert to a non-unit-length vector, which is all that is needed for
// some purposes.
//
//               Area    Edge    Diag   ToPointRaw  ToPoint  FromPoint
//              Ratio   Ratio   Ratio             (microseconds)
// -------------------------------------------------------------------
// Linear:      5.200   2.117   2.959      0.020     0.087     0.085
// Tangent:     1.414   1.414   1.704      0.237     0.299     0.258
// Quadratic:   2.082   1.802   1.932      0.033     0.096     0.108
//
// The worst-case cell aspect ratios are about the same with all three
// projections.  The maximum ratio of the longest edge to the shortest edge
// within the same cell is about 1.4 and the maximum ratio of the diagonals
// within the same cell is about 1.7.
//
// This data was produced using s2cell_test and s2cell_id_test.

#define S2_QUADRATIC_PROJECTION
//#define S2_LINEAR_PROJECTION
//#define S2_TAN_PROJECTION

namespace S2Geometry;

using System.Collections.ObjectModel;

public static partial class S2
{
    // Convert an s- or t-value to the corresponding u- or v-value.  This is
    // a non-linear transformation from [-1,1] to [-1,1] that attempts to
    // make the cell sizes more uniform.
    public static double STtoUV(double s)
    {
#if S2_LINEAR_PROJECTION
            return 2 * s - 1;
#elif S2_TAN_PROJECTION
            // Unfortunately, tan(S2Constants.M_PI_4) is slightly less than 1.0.  This isn't due to
            // a flaw in the implementation of tan(), it's because the derivative of
            // tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
            // point numbers on either side of the infinite-precision value of pi/4 have
            // tangents that are slightly below and slightly above 1.0 when rounded to
            // the nearest double-precision result.

            s = tan(S2Constants.M_PI_2 * s - S2Constants.M_PI_4);
            return s + (1.0 / (int64{1} << 53)) * s;
#elif S2_QUADRATIC_PROJECTION
        if (s >= 0.5) return (1 / 3.0) * (4 * s * s - 1);
        else return (1 / 3.0) * (1 - 4 * (1 - s) * (1 - s));
#else
            throw new NotImplementedException("Unknown value for S2_PROJECTION");
#endif
    }

    // The inverse of the STtoUV transformation.  Note that it is not always
    // true that UVtoST(STtoUV(x)) == x due to numerical errors.
    public static double UVtoST(double u)
    {
#if S2_LINEAR_PROJECTION
            return 0.5 * (u + 1);
#elif S2_TAN_PROJECTION
            volatile double a = atan(u);
            return (2 * S2Constants.M_1_PI) * (a + S2Constants.M_PI_4);
#elif S2_QUADRATIC_PROJECTION
        if (u >= 0) return 0.5 * Math.Sqrt(1 + 3 * u);
        else return 1 - 0.5 * Math.Sqrt(1 - 3 * u);
#else
            throw new NotImplementedException("Unknown value for S2_PROJECTION");
#endif
    }

    // Convert the i- or j-index of a leaf cell to the minimum corresponding s-
    // or t-value contained by that cell.  The argument must be in the range
    // [0..2**30], i.e. up to one position beyond the normal range of valid leaf
    // cell indices.
    public static double IJtoSTMin(int i)
    {
        Debug.Assert(i >= 0 && i <= S2.kLimitIJ);
        return (1.0 / S2.kLimitIJ) * i;
    }

    // Return the i- or j-index of the leaf cell containing the given
    // s- or t-value.  If the argument is outside the range spanned by valid
    // leaf cell indices, return the index of the closest valid leaf cell (i.e.,
    // return values are clamped to the range of valid leaf cell indices).
    public static int STtoIJ(double s)
    {
        return (int)Math.Max(0, Math.Min(S2.kLimitIJ - 1,
            (long)Math.Round(S2.kLimitIJ * s - 0.5)));
    }

    // Convert an si- or ti-value to the corresponding s- or t-value.
    public static double SiTitoST(uint si)
    {
        Debug.Assert(si <= S2.kMaxSiTi);
        return (1.0 / S2.kMaxSiTi) * si;
    }

    // Return the si- or ti-coordinate that is nearest to the given s- or
    // t-value.  The result may be outside the range of valid (si,ti)-values.
    public static uint STtoSiTi(double s)
    {
        // kMaxSiTi == 2^31, so the result doesn't fit in an Int32 when s == 1.
        return (uint)Math.Round(s * S2.kMaxSiTi);
    }

    // Convert (face, u, v) coordinates to a direction vector (not
    // necessarily unit length).
    public static S2Point FaceUVtoXYZ(int face, double u, double v)
    {
        return face switch
        {
            0 => new S2Point(1, u, v),
            1 => new S2Point(-u, 1, v),
            2 => new S2Point(-u, -v, 1),
            3 => new S2Point(-1, -v, -u),
            4 => new S2Point(v, -1, -u),
            _ => new S2Point(v, u, -1),
        };
    }
    public static S2Point FaceUVtoXYZ(int face, R2Point uv)
    {
        return FaceUVtoXYZ(face, uv[0], uv[1]);
    }

    // If the dot product of p with the given face normal is positive,
    // set the corresponding u and v values (which may lie outside the range
    // [-1,1]) and return true.  Otherwise return false.
    public static bool FaceXYZtoUV(int face, S2Point p, out double pu, out double pv)
    {
        pu = 0;
        pv = 0;
        if (face < 3)
        {
            if (p[face] <= 0) return false;
        }
        else
        {
            if (p[face - 3] >= 0) return false;
        }
        ValidFaceXYZtoUV(face, p, out pu, out pv);
        return true;
    }

    public static bool FaceXYZtoUV(int face, S2Point p, out R2Point puv)
    {
        var result = FaceXYZtoUV(face, p, out var x, out var y);
        puv = new R2Point(x, y);
        return result;
    }

    // Given a *valid* face for the given point p (meaning that dot product
    // of p with the face normal is positive), return the corresponding
    // u and v values (which may lie outside the range [-1,1]).
    public static void ValidFaceXYZtoUV(int face, S2Point p, out double pu, out double pv)
    {
        Debug.Assert(p.DotProd(GetNorm(face)) > 0);
        switch (face)
        {
            case 0: pu = p[1] / p[0]; pv = p[2] / p[0]; break;
            case 1: pu = -p[0] / p[1]; pv = p[2] / p[1]; break;
            case 2: pu = -p[0] / p[2]; pv = -p[1] / p[2]; break;
            case 3: pu = p[2] / p[0]; pv = p[1] / p[0]; break;
            case 4: pu = p[2] / p[1]; pv = -p[0] / p[1]; break;
            default: pu = -p[1] / p[2]; pv = -p[0] / p[2]; break;
        }
    }
    public static void ValidFaceXYZtoUV(int face, S2Point p, out R2Point puv)
    {
        ValidFaceXYZtoUV(face, p, out var x, out var y);
        puv = new R2Point(x, y);
    }

    // Transform the given point P to the (u,v,w) coordinate frame of the given
    // face (where the w-axis represents the face normal).
    public static S2Point FaceXYZtoUVW(int face, S2Point p)
    {
        // The result coordinates are simply the dot products of P with the (u,v,w)
        // axes for the given face (see kFaceUVWAxes).
        return face switch
        {
            0 => new S2Point(p.Y, p.Z, p.X),
            1 => new S2Point(-p.X, p.Z, p.Y),
            2 => new S2Point(-p.X, -p.Y, p.Z),
            3 => new S2Point(-p.Z, -p.Y, -p.X),
            4 => new S2Point(-p.Z, p.X, -p.Y),
            _ => new S2Point(p.Y, p.X, -p.Z),
        };
    }

    // Return the face containing the given direction vector.  (For points on
    // the boundary between faces, the result is arbitrary but repeatable.)
    public static int GetFace(S2Point p)
    {
        int face = p.LargestAbsComponent();
        if (p[face] < 0) face += 3;
        return face;
    }

    // Convert a direction vector (not necessarily unit length) to
    // (face, u, v) coordinates.
    public static int XYZtoFaceUV(S2Point p, out double pu, out double pv)
    {
        int face = GetFace(p);
        ValidFaceXYZtoUV(face, p, out pu, out pv);
        return face;
    }
    public static int XYZtoFaceUV(S2Point p, out R2Point puv)
    {
        var res = XYZtoFaceUV(p, out var x, out var y);
        puv = new R2Point(x, y);
        return res;
    }

    // Convert a direction vector (not necessarily unit length) to
    // (face, si, ti) coordinates and, if p is exactly equal to the center of a
    // cell, return the level of this cell (-1 otherwise).
    public static int XYZtoFaceSiTi(S2Point p, out int face, out uint si, out uint ti)
    {
        face = XYZtoFaceUV(p, out var u, out var v);
        si = STtoSiTi(UVtoST(u));
        ti = STtoSiTi(UVtoST(v));
        // If the levels corresponding to si,ti are not equal, then p is not a cell
        // center.  The si,ti values 0 and kMaxSiTi need to be handled specially
        // because they do not correspond to cell centers at any valid level; they
        // are mapped to level -1 by the code below.
        int level = S2.kMaxCellLevel - BitsUtils.FindLSBSetNonZero(si | S2.kMaxSiTi);
        if (level < 0 || level != S2.kMaxCellLevel - BitsUtils.FindLSBSetNonZero(ti | S2.kMaxSiTi))
        {
            return -1;
        }
        Debug.Assert(level <= S2.kMaxCellLevel);
        // In infinite precision, this test could be changed to ST == SiTi. However,
        // due to rounding errors, UVtoST(XYZtoFaceUV(FaceUVtoXYZ(STtoUV(...)))) is
        // not idempotent. On the other hand, center_raw is computed exactly the same
        // way p was originally computed (if it is indeed the center of an S2Cell):
        // the comparison can be exact.
        var center = FaceSiTitoXYZ(face, si, ti).Normalize();
        return p == center ? level : -1;
    }

    // Convert (face, si, ti) coordinates to a direction vector (not necessarily
    // unit length).
    public static S2Point FaceSiTitoXYZ(int face, uint si, uint ti)
    {
        double u = STtoUV(SiTitoST(si));
        double v = STtoUV(SiTitoST(ti));
        return FaceUVtoXYZ(face, u, v);
    }

    // Return the right-handed normal (not necessarily unit length) for an
    // edge in the direction of the positive v-axis at the given u-value on
    // the given face.  (This vector is perpendicular to the plane through
    // the sphere origin that contains the given edge.)
    public static S2Point GetUNorm(int face, double u)
    {
        return face switch
        {
            0 => new S2Point(u, -1, 0),
            1 => new S2Point(1, u, 0),
            2 => new S2Point(1, 0, u),
            3 => new S2Point(-u, 0, 1),
            4 => new S2Point(0, -u, 1),
            _ => new S2Point(0, -1, -u),
        };
    }

    // Return the right-handed normal (not necessarily unit length) for an
    // edge in the direction of the positive u-axis at the given v-value on
    // the given face.
    public static S2Point GetVNorm(int face, double v)
    {
        return face switch
        {
            0 => new S2Point(-v, 0, 1),
            1 => new S2Point(0, -v, 1),
            2 => new S2Point(0, -1, -v),
            3 => new S2Point(v, -1, 0),
            4 => new S2Point(1, v, 0),
            _ => new S2Point(1, 0, v),
        };
    }

    // Return the unit-length normal, u-axis, or v-axis for the given face.
    public static S2Point GetNorm(int face)
    {
        return GetUVWAxis(face, 2);
    }
    public static S2Point GetUAxis(int face)
    {
        return GetUVWAxis(face, 0);
    }
    public static S2Point GetVAxis(int face)
    {
        return GetUVWAxis(face, 1);
    }

    // Return the given axis of the given face (u=0, v=1, w=2).
    public static S2Point GetUVWAxis(int face, int axis)
    {
        return new S2Point(kFaceUVWAxes[face][axis]);
    }

    // With respect to the (u,v,w) coordinate system of a given face, return the
    // face that lies in the given direction (negative=0, positive=1) of the
    // given axis (u=0, v=1, w=2).  For example, GetUVWFace(4, 0, 1) returns the
    // face that is adjacent to face 4 in the positive u-axis direction.
    public static int GetUVWFace(int face, int axis, bool direction) => GetUVWFace(face, axis, direction ? 1 : 0);
    public static int GetUVWFace(int face, int axis, int direction)
    {
        Debug.Assert(face >= 0 && face <= 5);
        Debug.Assert(axis >= 0 && axis <= 2);
        Debug.Assert(new[] { 0, 1 }.Contains(direction));
        return kFaceUVWFaces[face][axis][direction];
    }

    #region Tables

    // The canonical Hilbert traversal order looks like an inverted 'U':
    // the subcells are visited in the order (0,0), (0,1), (1,1), (1,0).
    // The following tables encode the traversal order for various
    // orientations of the Hilbert curve (axes swapped and/or directions
    // of the axes reversed).

    // kIJtoPos[orientation][ij] . pos
    //
    // Given a cell orientation and the (i,j)-index of a subcell (0=(0,0),
    // 1=(0,1), 2=(1,0), 3=(1,1)), return the order in which this subcell is
    // visited by the Hilbert curve (a position in the range [0..3]).
    public static readonly ReadOnlyCollection<ReadOnlyCollection<int>> kIJtoPos = new List<ReadOnlyCollection<int>>
    {
        // (0,0) (0,1) (1,0) (1,1)
        new List<int>{ 0, 1, 3, 2 }.AsReadOnly(),  // canonical order
        new List<int>{ 0, 3, 1, 2 }.AsReadOnly(),  // axes swapped
        new List<int>{ 2, 3, 1, 0 }.AsReadOnly(),  // bits inverted
        new List<int>{ 2, 1, 3, 0 }.AsReadOnly(),  // swapped & inverted
    }.AsReadOnly();

    // kPosToIJ[orientation][pos] . ij
    //
    // Return the (i,j) index of the subcell at the given position 'pos' in the
    // Hilbert curve traversal order with the given orientation.  This is the
    // inverse of the previous table:
    //
    //   kPosToIJ[r][kIJtoPos[r][ij]] == ij
    public static readonly ReadOnlyCollection<ReadOnlyCollection<int>> kPosToIJ = new List<ReadOnlyCollection<int>>
    {
        // 0  1  2  3
        new List<int>{ 0, 1, 3, 2 }.AsReadOnly(),    // canonical order:    (0,0), (0,1), (1,1), (1,0)
        new List<int>{ 0, 2, 3, 1 }.AsReadOnly(),    // axes swapped:       (0,0), (1,0), (1,1), (0,1)
        new List<int>{ 3, 2, 0, 1 }.AsReadOnly(),    // bits inverted:      (1,1), (1,0), (0,0), (0,1)
        new List<int>{ 3, 1, 0, 2 }.AsReadOnly(),    // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
    }.AsReadOnly();

    // kPosToOrientation[pos] . orientation_modifier
    //
    // Return a modifier indicating how the orientation of the child subcell
    // with the given traversal position [0..3] is related to the orientation
    // of the parent cell.  The modifier should be XOR-ed with the parent
    // orientation to obtain the curve orientation in the child.
    public static readonly ReadOnlyCollection<int> kPosToOrientation = new List<int>
    {
        S2.kSwapMask,
        0,
        0,
        S2.kInvertMask + S2.kSwapMask,
    }.AsReadOnly();

    // The precomputed neighbors of each face (see GetUVWFace).
    private static readonly ReadOnlyCollection<ReadOnlyCollection<ReadOnlyCollection<int>>> kFaceUVWFaces = new List<ReadOnlyCollection<ReadOnlyCollection<int>>> {
            new List<ReadOnlyCollection<int>> {
                new List<int>{ 4, 1 }.AsReadOnly(),
                new List<int>{ 5, 2 }.AsReadOnly(),
                new List<int>{ 3, 0 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<int>> {
                new List<int>{ 0, 3 }.AsReadOnly(),
                new List<int>{ 5, 2 }.AsReadOnly(),
                new List<int>{ 4, 1 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<int>>{
                new List<int>{ 0, 3 }.AsReadOnly(),
                new List<int>{ 1, 4 }.AsReadOnly(),
                new List<int>{ 5, 2 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<int>>{
                new List<int>{ 2, 5 }.AsReadOnly(),
                new List<int>{ 1, 4 }.AsReadOnly(),
                new List<int>{ 0, 3 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<int>>{
                new List<int>{ 2, 5 }.AsReadOnly(),
                new List<int>{ 3, 0 }.AsReadOnly(),
                new List<int>{ 1, 4 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<int>>{
                new List<int>{ 4, 1 }.AsReadOnly(),
                new List<int>{ 3, 0 }.AsReadOnly(),
                new List<int>{ 2, 5 }.AsReadOnly()
            }.AsReadOnly()
        }.AsReadOnly();

    // The U,V,W axes for each face.
    private static readonly ReadOnlyCollection<ReadOnlyCollection<ReadOnlyCollection<double>>> kFaceUVWAxes = new List<ReadOnlyCollection<ReadOnlyCollection<double>>>
    {
            new List<ReadOnlyCollection<double>>{
                new List<double>{ 0,  1,  0 }.AsReadOnly(),
                new List<double>{ 0,  0,  1 }.AsReadOnly(),
                new List<double>{ 1,  0,  0 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<double>>{
                new List<double>{-1,  0,  0 }.AsReadOnly(),
                new List<double>{ 0,  0,  1 }.AsReadOnly(),
                new List<double>{ 0,  1,  0 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<double>>{
                new List<double>{-1,  0,  0 }.AsReadOnly(),
                new List<double>{ 0, -1,  0 }.AsReadOnly(),
                new List<double>{ 0,  0,  1 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<double>>{
                new List<double>{ 0,  0, -1 }.AsReadOnly(),
                new List<double>{ 0, -1,  0 }.AsReadOnly(),
                new List<double>{-1,  0,  0 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<double>>{
                new List<double>{ 0,  0, -1 }.AsReadOnly(),
                new List<double>{ 1,  0,  0 }.AsReadOnly(),
                new List<double>{ 0, -1,  0 }.AsReadOnly()
            }.AsReadOnly(),
            new List<ReadOnlyCollection<double>>{
                new List<double>{ 0,  1,  0 }.AsReadOnly(),
                new List<double>{ 1,  0,  0 }.AsReadOnly(),
                new List<double>{ 0,  0, -1 }.AsReadOnly()
            }.AsReadOnly()
    }.AsReadOnly();

    #endregion
}
