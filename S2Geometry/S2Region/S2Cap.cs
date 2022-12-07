// S2Cap represents a disc-shaped region defined by a center and radius.
// Technically this shape is called a "spherical cap" (rather than disc)
// because it is not planar; the cap represents a portion of the sphere that
// has been cut off by a plane.  The boundary of the cap is the circle defined
// by the intersection of the sphere and the plane.  For containment purposes,
// the cap is a closed set, i.e. it contains its boundary.
//
// For the most part, you can use a spherical cap wherever you would use a
// disc in planar geometry.  The radius of the cap is measured along the
// surface of the sphere (rather than the straight-line distance through the
// interior).  Thus a cap of radius Pi/2 is a hemisphere, and a cap of radius
// Pi covers the entire sphere.
//
// A cap can also be defined by its center point and height.  The height is
// simply the distance from the center point to the cutoff plane.  There is
// also support for empty and full caps, which contain no points and all
// points respectively.
//
// This class is intended to be copied by value as desired.  It uses the
// default copy constructor and assignment operator, however it is not a
// "plain old datatype" (POD) because it has virtual functions.

namespace S2Geometry;

public readonly record struct S2Cap : IS2Region<S2Cap>, IDecoder<S2Cap>
{
    #region Fields, Constants

    public S2Point Center { get; init; }
    public S1ChordAngle Radius { get; init; }

    // an empty cap, i.e. a cap that contains no points.
    public static readonly S2Cap Empty = new(new S2Point(1, 0, 0), S1ChordAngle.Negative);

    // Return a full cap, i.e. a cap that contains all points.
    public static readonly S2Cap Full = new(new S2Point(1, 0, 0), S1ChordAngle.Straight);

    #endregion

    #region Constructors

    // Constructs a cap with the given center and radius.  A negative radius
    // yields an empty cap; a radius of 180 degrees or more yields a full cap
    // (containing the entire sphere).  "center" should be unit length.
    public S2Cap(S2Point center, S1Angle radius)
    {
        Center = center;
        Radius = new S1ChordAngle(S1Angle.Min(radius, S1Angle.FromRadians(Math.PI)));
        // The "min" calculation above is necessary to handle S1Angle.Infinity().
        System.Diagnostics.Debug.Assert(IsValid());
    }

    // Constructs a cap where the angle is expressed as an S1ChordAngle.  This
    // constructor is more efficient than the one above.
    public S2Cap(S2Point center, S1ChordAngle radius)
    {
        Center = center;
        Radius = radius;
        System.Diagnostics.Debug.Assert(IsValid());
    }

    #endregion

    #region Factories

    // Convenience function that creates a cap containing a single point.  This
    // method is more efficient that the S2Cap(center, radius) constructor.
    public static S2Cap FromPoint(S2Point center)
    {
        return new S2Cap(center, S1ChordAngle.Zero);
    }

    // Returns a cap with the given center and height (see comments above).  A
    // negative height yields an empty cap; a height of 2 or more yields a full
    // cap.  "center" should be unit length.
    public static S2Cap FromCenterHeight(S2Point center, double height)
    {
        return new S2Cap(center, S1ChordAngle.FromLength2(2 * height));
    }

    // Return a cap with the given center and surface area.  Note that the area
    // can also be interpreted as the solid angle subtended by the cap (because
    // the sphere has unit radius).  A negative area yields an empty cap; an
    // area of 4*Pi or more yields a full cap.  "center" should be unit length.
    public static S2Cap FromCenterArea(S2Point center, double area)
    {
        return new S2Cap(center, S1ChordAngle.FromLength2(area / Math.PI));
    }

    #endregion

    #region S2Cap

    // Returns the height of the cap, i.e. the distance from the center point to
    // the cutoff plane.
    public double Height() => 0.5 * Radius.Length2;

    // Return the cap radius as an S1Angle.  (Note that the cap angle is stored
    // internally as an S1ChordAngle, so this method requires a trigonometric
    // operation and may yield a slightly different result than the value passed
    // to the (S2Point, S1Angle) constructor.)
    public S1Angle RadiusAngle() => Radius.ToAngle();

    // Return the area of the cap.
    public double Area() => S2.M_2_PI * Math.Max(0.0, Height());

    // Return the true centroid of the cap multiplied by its surface area (see
    // s2centroids.h for details on centroids). The result lies on the ray from
    // the origin through the cap's center, but it is not unit length. Note that
    // if you just want the "surface centroid", i.e. the normalized result, then
    // it is much simpler just to call center().
    //
    // The reason for multiplying the result by the cap area is to make it
    // easier to compute the centroid of more complicated shapes.  The centroid
    // of a union of disjoint regions can be computed simply by adding their
    // GetCentroid() results. Caveat: for caps that contain a single point
    // (i.e., zero radius), this method always returns the origin (0, 0, 0).
    // This is because shapes with no area don't affect the centroid of a
    // union whose total area is positive.
    public S2Point Centroid()
    {
        // From symmetry, the centroid of the cap must be somewhere on the line
        // from the origin to the center of the cap on the surface of the sphere.
        // When a sphere is divided into slices of constant thickness by a set of
        // parallel planes, all slices have the same surface area. This implies
        // that the radial component of the centroid is simply the midpoint of the
        // range of radial distances spanned by the cap. That is easily computed
        // from the cap height.
        if (IsEmpty()) return S2Point.Empty;
        double r = 1.0 - 0.5 * Height();
        return r * Area() * Center;
    }

    // We allow negative heights (to represent empty caps) but heights are
    // normalized so that they do not exceed 2.
    public bool IsValid() => Center.IsUnitLength() && Radius.Length2 <= 4;

    // Return true if the cap is empty, i.e. it contains no points.
    public bool IsEmpty() => Radius.IsNegative();

    // Return true if the cap is full, i.e. it contains all points.
    public bool IsFull() => Radius.Length2 == 4;

    // Return the complement of the interior of the cap.  A cap and its
    // complement have the same boundary but do not share any interior points.
    // The complement operator is not a bijection because the complement of a
    // singleton cap (containing a single point) is the same as the complement
    // of an empty cap.
    public S2Cap Complement()
    {
        // The complement of a full cap is an empty cap, not a singleton.
        // Also make sure that the complement of an empty cap is full.
        if (IsFull()) return Empty;
        if (IsEmpty()) return Full;
        return new S2Cap(-Center, S1ChordAngle.FromLength2(4 - Radius.Length2));
    }

    // Return true if and only if this cap contains the given other cap
    // (in a set containment sense, e.g. every cap contains the empty cap).
    public bool Contains(S2Cap other)
    {
        if (IsFull() || other.IsEmpty()) return true;
        return Radius >= new S1ChordAngle(Center, other.Center) + other.Radius;
    }

    // Return true if and only if this cap intersects the given other cap,
    // i.e. whether they have any points in common.
    public bool Intersects(S2Cap other)
    {
        if (IsEmpty() || other.IsEmpty()) return false;
        return Radius + other.Radius >= new S1ChordAngle(Center, other.Center);
    }

    // Return true if and only if the interior of this cap intersects the
    // given other cap.  (This relationship is not symmetric, since only
    // the interior of this cap is used.)
    public bool InteriorIntersects(S2Cap other)
    {
        // Make sure this cap has an interior and the other cap is non-empty.
        if (Radius.Length2 <= 0 || other.IsEmpty()) return false;
        return Radius + other.Radius > new S1ChordAngle(Center, other.Center);
    }

    // Return true if and only if the given point is contained in the interior
    // of the cap (i.e. the cap excluding its boundary).  "p" should be be a
    // unit-length vector.
    public bool InteriorContains(S2Point p)
    {
        System.Diagnostics.Debug.Assert(p.IsUnitLength());
        return IsFull() || new S1ChordAngle(Center, p) < Radius;
    }

    // Increase the cap height if necessary to include the given point.  If the
    // cap is empty then the center is set to the given point, but otherwise the
    // center is not changed.  "p" should be a unit-length vector.
    public S2Cap AddPoint(S2Point p)
    {
        // Compute the squared chord length, then convert it into a height.
        System.Diagnostics.Debug.Assert(p.IsUnitLength());
        var newCenter = Center;
        var newRadius = S1ChordAngle.Zero;
        if (IsEmpty())
        {
            newCenter = p;
        }
        else
        {
            // After calling cap.AddPoint(p), cap.Contains(p) must be true.  However
            // we don't need to do anything special to achieve this because Contains()
            // does exactly the same distance calculation that we do here.
            newRadius = S1ChordAngle.Max(Radius, new S1ChordAngle(Center, p));
        }
        return new S2Cap(newCenter, newRadius);
    }

    // Increase the cap height if necessary to include "other".  If the current
    // cap is empty it is set to the given other cap.
    public S2Cap AddCap(S2Cap other)
    {
        if (IsEmpty())
        {
            return other;
        }
        else if (!other.IsEmpty())
        {
            // We round up the distance to ensure that the cap is actually contained.
            // TODO(ericv): Do some error analysis in order to guarantee this.
            var dist = new S1ChordAngle(Center, other.Center) + other.Radius;
            var newRadius = S1ChordAngle.Max(Radius, dist.PlusError(S2.DoubleEpsilon * dist.Length2));
            return new S2Cap(Center, newRadius);
        }
        else
        {
            return this;
        }
    }

    // Return a cap that contains all points within a given distance of this
    // cap.  Note that any expansion of the empty cap is still empty.
    public S2Cap Expanded(S1Angle distance)
    {
        System.Diagnostics.Debug.Assert(distance.Radians >= 0);
        if (IsEmpty()) return Empty;
        return new S2Cap(Center, Radius + new S1ChordAngle(distance));
    }

    // Return the smallest cap which encloses this cap and "other".
    public S2Cap Union(S2Cap other)
    {
        if (Radius < other.Radius)
        {
            return other.Union(this);
        }
        if (IsFull() || other.IsEmpty())
        {
            return this;
        }
        // This calculation would be more efficient using S1ChordAngles.
        var this_radius = RadiusAngle();
        var other_radius = other.RadiusAngle();
        S1Angle distance = new(Center, other.Center);
        if (this_radius >= distance + other_radius)
        {
            return this;
        }

        var result_radius = 0.5 * (distance + this_radius + other_radius);
        var result_center = S2.GetPointOnLine(
            Center,
            other.Center,
            0.5 * (distance - this_radius + other_radius));
        return new S2Cap(result_center, result_radius);
    }

    // Here are some useful relationships between the cap height (h), the cap
    // radius (r), the maximum chord length from the cap's center (d), and the
    // radius of cap's base (a).
    //
    //     h = 1 - cos(r)
    //       = 2 * sin^2(r/2)
    //   d^2 = 2 * h
    //       = a^2 + h^2

    // Return true if the cap intersects "cell", given that the cap does contain
    // any of the cell vertices (supplied in "vertices", an array of length 4).
    private bool Intersects(S2Cell cell, S2Point[] vertices)
    {
        // Return true if this cap intersects any point of 'cell' excluding its
        // vertices (which are assumed to already have been checked).

        // If the cap is a hemisphere or larger, the cell and the complement of the
        // cap are both convex.  Therefore since no vertex of the cell is contained,
        // no other interior point of the cell is contained either.
        if (Radius >= S1ChordAngle.Right) return false;

        // We need to check for empty caps due to the center check just below.
        if (IsEmpty()) return false;

        // Optimization: return true if the cell contains the cap center.  (This
        // allows half of the edge checks below to be skipped.)
        if (cell.Contains(Center)) return true;

        // At this point we know that the cell does not contain the cap center,
        // and the cap does not contain any cell vertex.  The only way that they
        // can intersect is if the cap intersects the interior of some edge.

        double sin2_angle = Radius.Sin2();
        for (int k = 0; k < 4; ++k)
        {
            S2Point edge = cell.EdgeRaw(k);
            double dot = Center.DotProd(edge);
            if (dot > 0)
            {
                // The center is in the interior half-space defined by the edge.  We don't
                // need to consider these edges, since if the cap intersects this edge
                // then it also intersects the edge on the opposite side of the cell
                // (because we know the center is not contained with the cell).
                continue;
            }
            // The Norm2() factor is necessary because "edge" is not normalized.
            if (dot * dot > sin2_angle * edge.Norm2())
            {
                return false;  // Entire cap is on the exterior side of this edge.
            }
            // Otherwise, the great circle containing this edge intersects
            // the interior of the cap.  We just need to check whether the point
            // of closest approach occurs between the two edge endpoints.
            var dir = edge.CrossProd(Center);
            if (dir.DotProd(vertices[k]) < 0 && dir.DotProd(vertices[(k + 1) & 3]) > 0)
                return true;
        }
        return false;
    }

    // Return true if the cap center and height differ by at most "max_error"
    // from the given cap "other".
    public bool ApproxEquals(S2Cap other)
    {
        return ApproxEquals(other, S1Angle.FromRadians(1e-14));
    }

    public bool ApproxEquals(S2Cap other, S1Angle max_error_angle)
    {
        double max_error = max_error_angle.Radians;
        double r2 = Radius.Length2;
        double other_r2 = other.Radius.Length2;
        return (S2.ApproxEquals(Center, other.Center, max_error_angle) &&
                Math.Abs(r2 - other_r2) <= max_error) ||
               (IsEmpty() && other_r2 <= max_error) ||
               (other.IsEmpty() && r2 <= max_error) ||
               (IsFull() && other_r2 >= 2 - max_error) ||
               (other.IsFull() && r2 >= 2 - max_error);
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        return this;
    }

    public S2LatLngRect GetRectBound()
    {
        if (IsEmpty()) return S2LatLngRect.Empty;

        // Convert the center to a (lat,lng) pair, and compute the cap angle.
        var center_ll = new S2LatLng(Center);
        var cap_angle = RadiusAngle().Radians;

        var all_longitudes = false;
        var lat = new double[2];
        var lng = new double[2];
        lng[0] = -Math.PI;
        lng[1] = Math.PI;

        // Check whether cap includes the south pole.
        lat[0] = center_ll.LatRadians - cap_angle;
        if (lat[0] <= -S2.M_PI_2)
        {
            lat[0] = -S2.M_PI_2;
            all_longitudes = true;
        }
        // Check whether cap includes the north pole.
        lat[1] = center_ll.LatRadians + cap_angle;
        if (lat[1] >= S2.M_PI_2)
        {
            lat[1] = S2.M_PI_2;
            all_longitudes = true;
        }
        if (!all_longitudes)
        {
            // Compute the range of longitudes covered by the cap.  We use the law
            // of sines for spherical triangles.  Consider the triangle ABC where
            // A is the north pole, B is the center of the cap, and C is the point
            // of tangency between the cap boundary and a line of longitude.  Then
            // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
            // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
            // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
            // minus the latitude).  This formula also works for negative latitudes.
            //
            // The formula for sin(a) follows from the relationship h = 1 - cos(a).

            var sin_a = Math.Sin(Radius.Radians());
            var sin_c = Math.Cos(center_ll.LatRadians);
            if (sin_a <= sin_c)
            {
                var angle_A = Math.Asin(sin_a / sin_c);
                lng[0] = Math.IEEERemainder(center_ll.LngRadians - angle_A, S2.M_2_PI);
                lng[1] = Math.IEEERemainder(center_ll.LngRadians + angle_A, S2.M_2_PI);
            }
        }
        return new S2LatLngRect(new R1Interval(lat[0], lat[1]), new S1Interval(lng[0], lng[1]));
    }

    // Computes a covering of the S2Cap.  In general the covering consists of at
    // most 4 cells except for very large caps, which may need up to 6 cells.
    // The output is not sorted.
    public void GetCellUnionBound(List<S2CellId> cellsId)
    {
        // TODO(ericv): The covering could be made quite a bit tighter by mapping
        // the cap to a rectangle in (i,j)-space and finding a covering for that.
        cellsId.Clear();

        // Find the maximum level such that the cap contains at most one cell vertex
        // and such that S2CellId.AppendVertexNeighbors() can be called.
        int level = S2.kMinWidth.GetLevelForMinValue(RadiusAngle().Radians) - 1;

        // If level < 0, then more than three face cells are required.
        if (level < 0)
        {
            for (int face = 0; face < 6; ++face)
            {
                cellsId.Add(S2CellId.FromFace(face));
            }
        }
        else
        {
            // The covering consists of the 4 cells at the given level that share the
            // cell vertex that is closest to the cap center.
            new S2CellId(Center).AppendVertexNeighbors(level, cellsId);
        }
    }

    public bool Contains(S2Cell cell)
    {
        // If the cap does not contain all cell vertices, return false.
        // We check the vertices before taking the Complement() because we can't
        // accurately represent the complement of a very small cap (a height
        // of 2-epsilon is rounded off to 2).
        var vertices = new S2Point[4];
        for (int k = 0; k < 4; ++k)
        {
            vertices[k] = cell.Vertex(k);
            if (!Contains(vertices[k])) return false;
        }
        // Otherwise, return true if the complement of the cap does not intersect
        // the cell.  (This test is slightly conservative, because technically we
        // want Complement().InteriorIntersects() here.)
        return !Complement().Intersects(cell, vertices);
    }

    public bool MayIntersect(S2Cell cell)
    {
        // If the cap contains any cell vertex, return true.
        var vertices = new S2Point[4];
        for (int k = 0; k < 4; ++k)
        {
            vertices[k] = cell.Vertex(k);
            if (Contains(vertices[k])) return true;
        }
        return Intersects(cell, vertices);
    }

    /// <param name="p">should be a unit-length vector.</param>
    public bool Contains(S2Point p)
    {
        System.Diagnostics.Debug.Assert(p.IsUnitLength());
        return new S1ChordAngle(Center, p) <= Radius;
    }

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2Cap to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(4 * sizeof(double));

        encoder.PutDouble(Center.X);
        encoder.PutDouble(Center.Y);
        encoder.PutDouble(Center.Z);
        encoder.PutDouble(Radius.Length2);

        System.Diagnostics.Debug.Assert(encoder.Avail() >= 0);
    }

    // Decodes an S2Cap encoded with Encode().  Returns true on success.
    public static (bool, S2Cap) Decode(Decoder d)
    {
        if (d.Avail() < 4 * sizeof(double)) return (false, default);

        var x = d.GetDouble();
        var y = d.GetDouble();
        var z = d.GetDouble();
        var cen = new S2Point(x, y, z);
        var rad = S1ChordAngle.FromLength2(d.GetDouble());

        var result = new S2Cap(cen, rad);

#if s2debug
        System.Diagnostics.Debug.Assert(result.IsValid()); // Invalid S2Cap
#endif

        return (true, result);
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2Cap(Center, Radius);

    #endregion

    #region IEquatable

    /// <summary>
    /// Return true if two caps are identical.
    /// </summary>
    public bool Equals(S2Cap? other)
    {
        if (other is null) return false;

        return (Center == other.Value.Center && Radius == other.Value.Radius) ||
            (IsEmpty() && other.Value.IsEmpty()) ||
            (IsFull() && other.Value.IsFull());
    }

    public override int GetHashCode() => HashCode.Combine(Center, Radius);


    #endregion

    #region Object

    public override string ToString() => $"[Center={Center}, Radius={RadiusAngle()}]";

    #endregion
}
