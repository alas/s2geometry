// This class is a stopgap measure that allows some of the S2 spherical
// geometry machinery to be applied to planar geometry.  An S2R2Rect
// represents a closed axis-aligned rectangle in the (x,y) plane (an R2Rect),
// but it also happens to be a subtype of S2Region, which means that you can
// use an S2RegionCoverer to approximate it as a collection of S2CellIds.
//
// With respect to the S2Cell decomposition, an S2R2Rect is interpreted as a
// region of (s,t)-space on face 0.  In particular, the rectangle [0,1]x[0,1]
// corresponds to the S2CellId that covers all of face 0.  This means that
// only rectangles that are subsets of [0,1]x[0,1] can be approximated using
// the S2RegionCoverer interface.
//
// The S2R2Rect class is also a convenient way to find the (s,t)-region
// covered by a given S2CellId (see the FromCell and FromCellId methods).
//
// TODO(ericv): If the geometry library is extended to have better support
// for planar geometry, then this class should no longer be necessary.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator, however it is
// not a "plain old datatype" (POD) because it has virtual functions.

namespace S2Geometry;

public readonly record struct S2R2Rect(R2Rect Rect) : IS2Region<S2R2Rect>
{
    #region Fields, Constants

    // The canonical empty rectangle.  Use IsEmpty to test for empty
    // rectangles, since they have more than one representation.
    public static readonly S2R2Rect Empty = new(R2Rect.Empty);

    #endregion

    #region Constructors

    // Construct a rectangle from the given lower-left and upper-right points.
    public S2R2Rect(R2Point lo, R2Point hi)
        : this(new R2Rect(lo, hi)) { }

    // Construct a rectangle from the given intervals in x and y.  The two
    // intervals must either be both empty or both non-empty.
    public S2R2Rect(R1Interval x, R1Interval y)
        : this(new R2Rect(x, y)) { }

    #endregion

    #region Factories

    // Construct a rectangle that corresponds to the boundary of the given cell
    // is (s,t)-space.  Such rectangles are always a subset of [0,1]x[0,1].
    public static S2R2Rect FromCell(S2Cell cell)
    {
        // S2Cells have a more efficient GetSizeST() method than S2CellIds.
        double size = cell.SizeST();
        return FromCenterSize(cell.Id.CenterST(), new R2Point(size, size));
    }

    public static S2R2Rect FromCellId(S2CellId id)
    {
        double size = id.SizeST();
        return FromCenterSize(id.CenterST(), new R2Point(size, size));
    }

    // Construct a rectangle from a center point and size in each dimension.
    // Both components of size should be non-negative, i.e. this method cannot
    // be used to create an empty rectangle.
    public static S2R2Rect FromCenterSize(R2Point center, R2Point size) =>
        new(R2Rect.FromCenterSize(center, size));

    // Convenience method to construct a rectangle containing a single point.
    public static S2R2Rect FromPoint(R2Point p) => new(R2Rect.FromPoint(p));

    // Convenience method to construct the minimal bounding rectangle containing
    // the two given points.  This is equivalent to starting with an empty
    // rectangle and calling AddPoint() twice.  Note that it is different than
    // the S2R2Rect(lo, hi) constructor, where the first point is always
    // used as the lower-left corner of the resulting rectangle.
    public static S2R2Rect FromPointPair(R2Point p1, R2Point p2) => new(R2Rect.FromPointPair(p1, p2));

    #endregion

    #region S2R2Rect

    // Methods that allow the S2R2Rect to be accessed as a vector.
    public R1Interval this[int i] => Rect[i];

    // Accessor methods.
    public R1Interval X() => Rect.X;

    public R1Interval Y() => Rect.Y;

    public R2Point Lo() => Rect.Lo();

    public R2Point Hi() => Rect.Hi();

    // Return true if the rectangle is valid, which essentially just means
    // that if the bound for either axis is empty then both must be.
    public bool IsValid() => Rect.IsValid();

    // Return true if the rectangle is empty, i.e. it contains no points at all.
    public bool IsEmpty() => Rect.IsEmpty();

    // Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.
    // Vertex 0 is in the lower-left corner.  For convenience, the argument is
    // reduced modulo 4 to the range [0..3].
    public R2Point Vertex(int k) => Rect.GetVertex(k);

    // Return the vertex in direction "i" along the x-axis (0=left, 1=right) and
    // direction "j" along the y-axis (0=down, 1=up).  Equivalently, return the
    // vertex constructed by selecting endpoint "i" of the x-interval (0=lo,
    // 1=hi) and vertex "j" of the y-interval.
    public R2Point Vertex(int i, int j) => Rect.GetVertex(i, j);

    // Return the center of the rectangle in (x,y)-space
    // (in general this is not the center of the region on the sphere).
    public R2Point Center() => Rect.GetCenter();

    // Return the width and height of this rectangle in (x,y)-space.  Empty
    // rectangles have a negative width and height.
    public R2Point GetSize() => Rect.GetSize();

    // Return true if the rectangle contains the given point.  Note that
    // rectangles are closed regions, i.e. they contain their boundary.
    public bool Contains(R2Point p) => Rect.Contains(p);

    // Return true if and only if the given point is contained in the interior
    // of the region (i.e. the region excluding its boundary).
    public bool InteriorContains(R2Point p) => Rect.InteriorContains(p);

    // Return true if and only if the rectangle contains the given other
    // rectangle.
    public bool Contains(S2R2Rect other) => Rect.Contains(other.Rect);

    // Return true if and only if the interior of this rectangle contains all
    // points of the given other rectangle (including its boundary).
    public bool InteriorContains(S2R2Rect other) => Rect.InteriorContains(other.Rect);

    // Return true if this rectangle and the given other rectangle have any
    // points in common.
    public bool Intersects(S2R2Rect other) => Rect.Intersects(other.Rect);

    // Return true if and only if the interior of this rectangle intersects
    // any point (including the boundary) of the given other rectangle.
    public bool InteriorIntersects(S2R2Rect other) => Rect.InteriorIntersects(other.Rect);

    // Increase the size of the bounding rectangle to include the given point.
    // The rectangle is expanded by the minimum amount possible.
    public S2R2Rect AddPoint(R2Point p) => new(Rect.AddPoint(p));

    // Return the closest point in the rectangle to the given point "p".
    // The rectangle must be non-empty.
    public R2Point Project(R2Point p) => Rect.Project(p);

    // Return a rectangle that has been expanded on each side in the x-direction
    // by margin.x(), and on each side in the y-direction by margin.y().  If
    // either margin is negative, then shrink the interval on the corresponding
    // sides instead.  The resulting rectangle may be empty.  Any expansion of
    // an empty rectangle remains empty.
    public S2R2Rect Expanded(R2Point margin) => new(Rect.Expanded(margin));
    public S2R2Rect Expanded(double margin) => new(Rect.Expanded(margin));

    // Return the smallest rectangle containing the union of this rectangle and
    // the given rectangle.
    public S2R2Rect Union(S2R2Rect other) => new(Rect.Union(other.Rect));

    // Return the smallest rectangle containing the intersection of this
    // rectangle and the given rectangle.
    public S2R2Rect Intersection(S2R2Rect other) => new(Rect.Intersection(other.Rect));

    // Return true if the x- and y-intervals of the two rectangles are the same
    // up to the given tolerance (see r1interval.h for details).
    public bool ApproxEquals(S2R2Rect other) => ApproxEquals(other, S1Angle.FromRadians(S2.DoubleError));
    public bool ApproxEquals(S2R2Rect other, S1Angle max_error) => Rect.ApproxEquals(other.Rect, max_error.Radians);

    // Return the unit-length S2Point corresponding to the given point "p" in
    // the (s,t)-plane.  "p" need not be restricted to the range [0,1]x[0,1].
    public static S2Point ToS2Point(R2Point p) =>
        S2.FaceUVtoXYZ(0, S2.STtoUV(p.X), S2.STtoUV(p.Y)).Normalize();

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        if (IsEmpty()) return S2Cap.Empty;

        // The rectangle is a convex polygon on the sphere, since it is a subset of
        // one cube face.  Its bounding cap is also a convex region on the sphere,
        // and therefore we can bound the rectangle by just bounding its vertices.
        // We use the rectangle's center in (s,t)-space as the cap axis.  This
        // doesn't yield the minimal cap but it's pretty close.
        S2Cap cap = S2Cap.FromPoint(ToS2Point(Center()));
        for (int k = 0; k < 4; ++k)
        {
            cap = cap.AddPoint(ToS2Point(Vertex(k)));
        }
        return cap;
    }

    public S2LatLngRect GetRectBound()
    {
        // This is not very tight but hopefully good enough.
        return GetCapBound().GetRectBound();
    }

    public bool Contains(S2Point p)
    {
        if (S2.GetFace(p) != 0) return false;
        S2.ValidFaceXYZtoUV(0, p, out var u, out var v);
        return Contains(new R2Point(S2.UVtoST(u), S2.UVtoST(v)));
    }
    public bool Contains(S2Cell cell)
    {
        if (cell.Face != 0) return false;
        return Contains(FromCell(cell));
    }
    public bool MayIntersect(S2Cell cell)
    {
        if (cell.Face != 0) return false;
        return Intersects(FromCell(cell));
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2R2Rect(Rect);

    #endregion

    #region Object

    public override string ToString() => $"[Lo{Rect.Lo()}, Hi{Rect.Hi()}]";

    #endregion
}
