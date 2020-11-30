using System;

namespace S2Geometry
{
    /// <summary>
    /// An R2Rect represents a closed axis-aligned rectangle in the (x,y) plane.
    /// </summary>
    public readonly struct R2Rect : IEquatable<R2Rect>
    {
        #region Fields, Constants

        public readonly R1Interval X;
        public readonly R1Interval Y;

        /// <summary>
        /// The canonical empty rectangle.  Use IsEmpty to test for empty
        /// rectangles, since they have more than one representation.
        /// </summary>
        public static readonly R2Rect Empty = new R2Rect(R1Interval.Empty, R1Interval.Empty);

        #endregion

        #region Constructors

        /// <summary>
        /// Construct a rectangle from the given lower-left and upper-right points.
        /// </summary>
        public R2Rect(R2Point lo, R2Point hi)
        {
            X = new R1Interval(lo.X, hi.X);
            Y = new R1Interval(lo.Y, hi.Y);
            Assert.True(IsValid);
        }

        /// <summary>
        /// Construct a rectangle from the given intervals in x and y.  The two
        /// intervals must either be both empty or both non-empty.
        /// </summary>
        public R2Rect(R1Interval x, R1Interval y)
        {
            X = x;
            Y = y;
            Assert.True(IsValid);
        }

        #endregion

        #region Factories

        /// <summary>
        /// Construct a rectangle from a center point and size in each dimension.
        /// Both components of size should be non-negative, i.e. this method cannot
        /// be used to create an empty rectangle.
        /// </summary>
        public static R2Rect FromCenterSize(R2Point center, R2Point size)
        {
            var i1 = new R1Interval(center.X - 0.5 * size.X, center.X + 0.5 * size.X);
            var i2 = new R1Interval(center.Y - 0.5 * size.Y, center.Y + 0.5 * size.Y);
            return new R2Rect(i1, i2);
        }

        /// <summary>
        /// Convenience method to construct the minimal bounding rectangle containing
        /// the two given points.  This is equivalent to starting with an empty
        /// rectangle and calling AddPoint() twice.  Note that it is different than
        /// the R2Rect(lo, hi) constructor, where the first point is always
        /// used as the lower-left corner of the resulting rectangle.
        /// </summary>
        public static R2Rect FromPointPair(R2Point p1, R2Point p2)
        {
            var i1 = R1Interval.FromPointPair(p1.X, p2.X);
            var i2 = R1Interval.FromPointPair(p1.Y, p2.Y);
            return new R2Rect(i1, i2);
        }

        /// <summary>
        /// Convenience method to construct a rectangle containing a single point.
        /// </summary>
        public static R2Rect FromPoint(R2Point p) => new R2Rect(p, p);

        #endregion

        #region R2Rect

        public R2Point Lo => new R2Point(X.Lo, Y.Lo);

        public R2Point Hi => new R2Point(X.Hi, Y.Hi);

        public R1Interval this[int index] => index switch
        {
            0 => X,
            1 => Y,
            _ => throw new ArgumentOutOfRangeException(nameof(index))
        };

        /// <summary>
        /// Return true if the rectangle is valid, which essentially just means
        /// that if the bound for either axis is empty then both must be.
        /// </summary>
        /// <remarks>The x/y ranges must either be both empty or both non-empty.</remarks>
        public bool IsValid => X.IsEmpty == Y.IsEmpty;

        /// <summary>
        /// Return true if the rectangle is empty, i.e. it contains no points at all.
        /// </summary>
        public bool IsEmpty => X.IsEmpty;

        /// <summary>
        /// Return the center of the rectangle in (x,y)-space.
        /// </summary>
        public R2Point Center => new R2Point(X.Center, Y.Center);

        /// <summary>
        /// Return the width and height of this rectangle in (x,y)-space.  Empty
        /// rectangles have a negative width and height.
        /// </summary>
        public R2Point Size => new R2Point(X.Length, Y.Length);

        /// <summary>
        /// Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.
        /// Vertex 0 is in the lower-left corner.  For convenience, the argument is
        /// reduced modulo 4 to the range [0..3].
        /// </summary>
        public R2Point GetVertex(int k)
        {
            // Twiddle bits to return the points in CCW order (lower left, lower right,
            // upper right, upper left).
            var j = (k >> 1) & 1;
            var l = j ^ (k & 1);
            var vx = l == 0 ? X.Lo : X.Hi;
            var vy = j == 0 ? Y.Lo : Y.Hi;
            return new R2Point(vx, vy);
        }

        /// <summary>
        /// Return the vertex in direction "i" along the x-axis (0=left, 1=right) and
        /// direction "j" along the y-axis (0=down, 1=up).  Equivalently, return the
        /// vertex constructed by selecting endpoint "i" of the x-interval (0=lo,
        /// 1=hi) and vertex "j" of the y-interval.
        /// </summary>
        public R2Point GetVertex(int i, int j)
        {
            var vx = i == 0 ? X.Lo : X.Hi;
            var vy = j == 0 ? Y.Lo : Y.Hi;
            return new R2Point(vx, vy);
        }

        /// <summary>
        /// Return true if and only if the rectangle contains the given other
        /// rectangle.
        /// </summary>
        public bool Contains(R2Rect other) => X.Contains(other.X) && Y.Contains(other.Y);

        /// <summary>
        /// Return true if the rectangle contains the given point.  Note that
        /// rectangles are closed regions, i.e. they contain their boundary.
        /// </summary>
        public bool Contains(R2Point p) => X.Contains(p.X) && Y.Contains(p.Y);

        /// <summary>
        /// Return true if and only if the interior of this rectangle con
        /// points of the given other rectangle (including its boundary).tains all
        /// </summary>
        public bool InteriorContains(R2Rect other) => X.InteriorContains(other.X) && Y.InteriorContains(other.Y);

        /// <summary>
        /// Return true if and only if the given point is contained in the interior
        /// of the region (i.e. the region excluding its boundary).
        /// </summary>
        public bool InteriorContains(R2Point p) => X.InteriorContains(p.X) && Y.InteriorContains(p.Y);

        /// <summary>
        /// Return true if this rectangle and the given other rectangle have any
        /// points in common.
        /// </summary>
        public bool Intersects(R2Rect other) => X.Intersects(other.X) && Y.Intersects(other.Y);

        /// <summary>
        /// Return true if and only if the interior of this rectangle intersects
        /// any point (including the boundary) of the given other rectangle.
        /// </summary>
        public bool InteriorIntersects(R2Rect other) => X.InteriorIntersects(other.X) && Y.InteriorIntersects(other.Y);

        /// <summary>
        /// Expand the rectangle to include the given point.  The rectangle is
        /// expanded by the minimum amount possible.
        /// </summary>
        public R2Rect AddPoint(R2Point p) => new R2Rect(R1Interval.AddPoint(X, p.X), R1Interval.AddPoint(Y, p.Y));

        /// <summary>
        /// Expand the rectangle to include the given other rectangle.  This is the
        /// same as replacing the rectangle by the union of the two rectangles, but
        /// is somewhat more efficient.
        /// </summary>
        public R2Rect AddRect(R2Rect other) => new R2Rect(R1Interval.AddInterval(X, other.X), R1Interval.AddInterval(Y, other.Y));

        /// <summary>
        /// Return the closest point in the rectangle to the given point "p".
        /// The rectangle must be non-empty.
        /// </summary>
        public R2Point Project(R2Point p) => new R2Point(X.Project(p.X), Y.Project(p.Y));

        /// <summary>
        /// Return a rectangle that has been expanded on each side in the x-direction
        /// by margin.x(), and on each side in the y-direction by margin.y().  If
        /// either margin is empty, then shrink the interval on the corresponding
        /// sides instead.  The resulting rectangle may be empty.  Any expansion of
        /// an empty rectangle remains empty.
        /// </summary>
        public R2Rect Expanded(R2Point margin)
        {
            var xx = X.Expanded(margin.X);
            var yy = Y.Expanded(margin.Y);
            if (xx.IsEmpty || yy.IsEmpty) return Empty;
            return new R2Rect(xx, yy);
        }

        /// <summary>
        /// Return a rectangle that has been expanded on each side in the x-direction
        /// by margin.x(), and on each side in the y-direction by margin.y().  If
        /// either margin is empty, then shrink the interval on the corresponding
        /// sides instead.  The resulting rectangle may be empty.  Any expansion of
        /// an empty rectangle remains empty.
        /// </summary>
        public R2Rect Expanded(double margin) => Expanded(new R2Point(margin, margin));

        /// <summary>
        /// Return the smallest rectangle containing the union of this rectangle and
        /// the given rectangle.
        /// </summary>
        public R2Rect Union(R2Rect other) => new R2Rect(X.Union(other.X), Y.Union(other.Y));

        /// <summary>
        /// Return the smallest rectangle containing the intersection of this
        /// rectangle and the given rectangle.
        /// </summary>
        public R2Rect Intersection(R2Rect other)
        {
            var xx = X.Intersection(other.X);
            var yy = Y.Intersection(other.Y);
            if (xx.IsEmpty || yy.IsEmpty) return Empty;
            return new R2Rect(xx, yy);
        }

        /// <summary>
        /// Return true if the x- and y-intervals of the two rectangles are the same
        /// up to the given tolerance (see r1interval.h for details).
        /// </summary>
        public bool ApproxEquals(R2Rect other, double max_error = S2Constants.DoubleError)
            => X.ApproxEquals(other.X, max_error) && Y.ApproxEquals(other.Y, max_error);

        #endregion

        #region IEquatable

        /// <summary>
        /// Return true if two rectangles contains the same set of points.
        /// </summary>
        public bool Equals(R2Rect other) => X == other.X && Y == other.Y;

        public override bool Equals(object obj) => obj is R2Rect rect && Equals(rect);

        public override int GetHashCode() => HashCode.Combine(X, Y);

        /// <summary>
        /// Return true if two rectangles contain the same set of points.
        /// </summary>
        public static bool operator ==(R2Rect left, R2Rect right) => Equals(left, right);

        /// <summary>
        /// Return true if two rectangles do not contain the same set of points.
        /// </summary>
        public static bool operator !=(R2Rect left, R2Rect right) => !Equals(left, right);

        #endregion
        
        #region Object

        public override string ToString() => $"[Lo{Lo}, Hi{Hi}]"; 

        #endregion
    }
}
