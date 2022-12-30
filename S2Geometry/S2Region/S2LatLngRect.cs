// An S2LatLngRect represents a closed latitude-longitude rectangle.  It is
// capable of representing the empty and full rectangles as well as single
// points.  Note that the latitude-longitude space is considered to have a
// *cylindrical* topology rather than a spherical one, i.e. the poles have
// multiple lat/lng representations.  An S2LatLngRect may be defined so that
// includes some representations of a pole but not others.  Use the
// PolarClosure() method if you want to expand a rectangle so that it contains
// all possible representations of any contained poles.
//
// Because S2LatLngRect uses S1Interval to store the longitude range,
// longitudes of -180 degrees are treated specially.  Except for empty
// and full longitude spans, -180 degree longitudes will turn into +180
// degrees.  This sign flip causes lng_lo() to be greater than lng_hi(),
// indicating that the rectangle will wrap around through -180 instead of
// through +179. Thus the math is consistent within the library, but the sign
// flip can be surprising, especially when working with map projections where
// -180 and +180 are at opposite ends of the flattened map.  See the comments
// on S1Interval for more details.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator, however it is
// not a "plain old datatype" (POD) because it has virtual functions.

namespace S2Geometry;

public readonly record struct S2LatLngRect : IS2Region<S2LatLngRect>, IDecoder<S2LatLngRect>
{
    #region Fields, Constants

    public R1Interval Lat { get; init; }
    public S1Interval Lng { get; init; }

    // The canonical empty and full rectangles, as derived from the Empty
    // and Full R1 and S1 Intervals.

    // Empty: lat_lo=1, lat_hi=0, lng_lo=Pi, lng_hi=-Pi (radians)
    public static readonly S2LatLngRect Empty = new(R1Interval.Empty, S1Interval.Empty);

    // The full allowable range of latitudes and longitudes.
    public static readonly R1Interval FullLat = new(-S2.M_PI_2, S2.M_PI_2);
    public static readonly S1Interval FullLng = S1Interval.Full;

    // Full: lat_lo=-Pi/2, lat_hi=Pi/2, lng_lo=-Pi, lng_hi=Pi (radians)
    public static readonly S2LatLngRect Full = new(FullLat, FullLng);

    #endregion

    #region Constructors

    // Construct a rectangle from minimum and maximum latitudes and longitudes.
    // If lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude
    // line. Both points must be normalized, with lo.lat() <= hi.lat().
    // The rectangle contains all the points p such that 'lo' <= p <= 'hi',
    // where '<=' is defined in the obvious way.
    public S2LatLngRect(S2LatLng lo, S2LatLng hi)
    {
        Lat = new R1Interval(lo.LatRadians, hi.LatRadians);
        Lng = new S1Interval(lo.LngRadians, hi.LngRadians);
        if (!IsValid()) Debug.WriteLine($"Invalid rect: lo={lo}, hi={hi}");
    }

    // Construct a rectangle from latitude and longitude intervals.  The two
    // intervals must either be both empty or both non-empty, and the latitude
    // interval must not extend outside [-90, +90] degrees.
    // Note that both intervals (and hence the rectangle) are closed.
    public S2LatLngRect(R1Interval lat, S1Interval lng)
    {
        Lat = lat;
        Lng = lng;
        if (!IsValid()) Debug.WriteLine($"Invalid rect: lat={lat}, lng={lng}");
    }

    #endregion

    #region Factories

    // Construct a rectangle of the given size centered around the given point.
    // "center" needs to be normalized, but "size" does not.  The latitude
    // interval of the result is clamped to [-90,90] degrees, and the longitude
    // interval of the result is Full() if and only if the longitude size is
    // 360 degrees or more.  Examples of clamping (in degrees):
    //
    //   center=(80,170),  size=(40,60)   . lat=[60,90],   lng=[140,-160]
    //   center=(10,40),   size=(210,400) . lat=[-90,90],  lng=[-180,180]
    //   center=(-90,180), size=(20,50)   . lat=[-90,-80], lng=[155,-155]
    public static S2LatLngRect FromCenterSize(S2LatLng center, S2LatLng size)
    {
        return FromPoint(center).Expanded(0.5 * size);
    }

    // Construct a rectangle containing a single (normalized) point.
    public static S2LatLngRect FromPoint(S2LatLng p)
    {
        if (!p.IsValid()) Debug.WriteLine($"Invalid S2LatLng: {p}");

        return new S2LatLngRect(p, p);
    }

    // Construct the minimal bounding rectangle containing the two given
    // normalized points.  This is equivalent to starting with an empty
    // rectangle and calling AddPoint() twice.  Note that it is different than
    // the S2LatLngRect(lo, hi) constructor, where the first point is always
    // used as the lower-left corner of the resulting rectangle.
    public static S2LatLngRect FromPointPair(S2LatLng p1, S2LatLng p2)
    {
        if (!p1.IsValid()) Debug.WriteLine($"Invalid S2LatLng 1: {p1}");
        if (!p2.IsValid()) Debug.WriteLine($"Invalid S2LatLng 2: {p2}");

        return new S2LatLngRect(
            R1Interval.FromPointPair(p1.LatRadians, p2.LatRadians),
            S1Interval.FromPointPair(p1.LngRadians, p2.LngRadians));
    }

    #endregion

    #region S2LatLngRect

    // Accessor methods.
    public S1Angle LatLo() => S1Angle.FromRadians(Lat.Lo);

    public S1Angle LatHi() => S1Angle.FromRadians(Lat.Hi);

    public S1Angle LngLo() => S1Angle.FromRadians(Lng.Lo);

    public S1Angle LngHi() => S1Angle.FromRadians(Lng.Hi);

    public S2LatLng Lo() => new(LatLo(), LngLo());

    public S2LatLng Hi() => new(LatHi(), LngHi());

    // Returns true if the rectangle is valid, which essentially just means
    // that the latitude bounds do not exceed Pi/2 in absolute value and
    // the longitude bounds do not exceed Pi in absolute value.  Also, if
    // either the latitude or longitude bound is empty then both must be.
    public bool IsValid()
    {
        // The lat/lng ranges must either be both empty or both non-empty.
        return Math.Abs(Lat.Lo) <= S2.M_PI_2 &&
               Math.Abs(Lat.Hi) <= S2.M_PI_2 &&
               Lng.IsValid() &&
               Lat.IsEmpty() == Lng.IsEmpty();
    }

    // Returns true if the rectangle is empty, i.e. it contains no points at all.
    public bool IsEmpty() => Lat.IsEmpty();

    // Returns true if the rectangle is full, i.e. it contains all points.
    public bool IsFull() => Lat == FullLat && Lng.IsFull();

    // Returns true if the rectangle is a point, i.e. lo() == hi()
    public bool IsPoint() => Lat.Lo == Lat.Hi && Lng.Lo == Lng.Hi;

    // Returns true if lng_.Lo > lng_.Hi, i.e. the rectangle crosses
    // the 180 degree longitude line.
    public bool IsInverted() => Lng.IsInverted();

    // Returns the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order
    // (lower left, lower right, upper right, upper left).  For convenience, the
    // argument is reduced modulo 4 to the range [0..3].
    public S2LatLng Vertex(int k)
    {
        // Twiddle bits to return the points in CCW order (lower left, lower right,
        // upper right, upper left).
        int i = (k >> 1) & 1;
        return S2LatLng.FromRadians(Lat[i], Lng[i ^ (k & 1)]);
    }

    // Returns the center of the rectangle in latitude-longitude space
    // (in general this is not the center of the region on the sphere).
    public S2LatLng Center()
    {
        return S2LatLng.FromRadians(Lat.GetCenter(), Lng.GetCenter());
    }

    // Returns the width and height of this rectangle in latitude-longitude
    // space.  Empty rectangles have a negative width and height.
    public S2LatLng Size()
    {
        return S2LatLng.FromRadians(Lat.GetLength(), Lng.GetLength());
    }

    // Returns the surface area of this rectangle on the unit sphere.
    public double Area()
    {
        if (IsEmpty()) return 0.0;
        // This is the size difference of the two spherical caps, multiplied by
        // the longitude ratio.
        return Lng.GetLength() * (Math.Sin(LatHi().Radians) - Math.Sin(LatLo().Radians));
    }

    // Returns the true centroid of the rectangle multiplied by its surface area
    // (see s2centroids.h for details on centroids).  The result is not unit
    // length, so you may want to normalize it.  Note that in general the
    // centroid is *not* at the center of the rectangle, and in fact it may not
    // even be contained by the rectangle.  (It is the "center of mass" of the
    // rectangle viewed as subset of the unit sphere, i.e. it is the point in
    // space about which this curved shape would rotate.)
    //
    // The reason for multiplying the result by the rectangle area is to make it
    // easier to compute the centroid of more complicated shapes.  The centroid
    // of a union of disjoint regions can be computed simply by adding their
    // GetCentroid() results.
    public S2Point Centroid()
    {
        // When a sphere is divided into slices of constant thickness by a set of
        // parallel planes, all slices have the same surface area.  This implies
        // that the z-component of the centroid is simply the midpoint of the
        // z-interval spanned by the S2LatLngRect.
        //
        // Similarly, it is easy to see that the (x,y) of the centroid lies in the
        // plane through the midpoint of the rectangle's longitude interval.  We
        // only need to determine the distance "d" of this point from the z-axis.
        //
        // Let's restrict our attention to a particular z-value.  In this z-plane,
        // the S2LatLngRect is a circular arc.  The centroid of this arc lies on a
        // radial line through the midpoint of the arc, and at a distance from the
        // z-axis of
        //
        //     r * (sin(alpha) / alpha)
        //
        // where r = Math.Sqrt(1-z^2) is the radius of the arc, and "alpha" is half of
        // the arc length (i.e., the arc covers longitudes [-alpha, alpha]).
        //
        // To find the centroid distance from the z-axis for the entire rectangle,
        // we just need to integrate over the z-interval.  This gives
        //
        //    d = Integrate[Math.Sqrt(1-z^2)*sin(alpha)/alpha, z1..z2] / (z2 - z1)
        //
        // where [z1, z2] is the range of z-values covered by the rectangle.  This
        // simplifies to
        //
        //    d = sin(alpha)/(2*alpha*(z2-z1))*(z2*r2 - z1*r1 + theta2 - theta1)
        //
        // where [theta1, theta2] is the latitude interval, z1=sin(theta1),
        // z2=sin(theta2), r1=cos(theta1), and r2=cos(theta2).
        //
        // Finally, we want to return not the centroid itself, but the centroid
        // scaled by the area of the rectangle.  The area of the rectangle is
        //
        //    A = 2 * alpha * (z2 - z1)
        //
        // which fortunately appears in the denominator of "d".

        if (IsEmpty()) return S2Point.Empty;
        var z1 = Math.Sin(LatLo().Radians);
        var z2 = Math.Sin(LatHi().Radians);
        var r1 = Math.Cos(LatLo().Radians);
        var r2 = Math.Cos(LatHi().Radians);
        var alpha = 0.5 * Lng.GetLength();
        var r = Math.Sin(alpha) * (r2 * z2 - r1 * z1 + Lat.GetLength());
        var lng = Lng.GetCenter();
        var z = alpha * (z2 + z1) * (z2 - z1);  // scaled by the area
        return new S2Point(r * Math.Cos(lng), r * Math.Sin(lng), z);
    }

    // More efficient version of Contains() that accepts a S2LatLng rather than
    // an S2Point.  The argument must be normalized.
    public bool Contains(S2LatLng ll)
    {
        if (!ll.IsValid()) Debug.WriteLine($"Invalid S2LatLng in S2LatLngRect.Contains: {ll}");

        return (Lat.Contains(ll.LatRadians) &&
                Lng.Contains(ll.LngRadians));
    }

    // Returns true if and only if the given point is contained in the interior
    // of the region (i.e. the region excluding its boundary).  The point 'p'
    // does not need to be normalized.
    public bool InteriorContains(S2Point p)
    {
        return InteriorContains(new S2LatLng(p));
    }

    // More efficient version of InteriorContains() that accepts a S2LatLng
    // rather than an S2Point.  The argument must be normalized.
    public bool InteriorContains(S2LatLng ll)
    {
        if (!ll.IsValid()) Debug.WriteLine($"Invalid S2LatLng in S2LatLngRect.InteriorContains: {ll}");

        return (Lat.InteriorContains(ll.LatRadians) &&
                Lng.InteriorContains(ll.LngRadians));
    }

    // Returns true if and only if the rectangle contains the given other
    // rectangle.
    public bool Contains(S2LatLngRect other)
    {
        return Lat.Contains(other.Lat) && Lng.Contains(other.Lng);
    }

    // Returns true if and only if the interior of this rectangle contains all
    // points of the given other rectangle (including its boundary).
    public bool InteriorContains(S2LatLngRect other)
    {
        return (Lat.InteriorContains(other.Lat) &&
                Lng.InteriorContains(other.Lng));
    }

    // Returns true if this rectangle and the given other rectangle have any
    // points in common.
    public bool Intersects(S2LatLngRect other)
    {
        return Lat.Intersects(other.Lat) && Lng.Intersects(other.Lng);
    }

    // Returns true if this rectangle intersects the given cell.  (This is an
    // exact test and may be fairly expensive, see also MayIntersect below.)
    public bool Intersects(S2Cell cell)
    {
        // First we eliminate the cases where one region completely contains the
        // other.  Once these are disposed of, then the regions will intersect
        // if and only if their boundaries intersect.

        if (IsEmpty()) return false;
        if (Contains(cell.CenterRaw())) return true;
        if (cell.Contains(Center().ToPoint())) return true;

        // Quick rejection test (not required for correctness).
        if (!Intersects(cell.GetRectBound())) return false;

        // Precompute the cell vertices as points and latitude-longitudes.  We also
        // check whether the S2Cell contains any corner of the rectangle, or
        // vice-versa, since the edge-crossing tests only check the edge interiors.

        var cell_v = new S2Point[4];
        var cell_ll = new S2LatLng[4];
        for (int i = 0; i < 4; ++i)
        {
            cell_v[i] = cell.Vertex(i);  // Must be normalized.
            cell_ll[i] = new S2LatLng(cell_v[i]);
            if (Contains(cell_ll[i])) return true;
            if (cell.Contains(Vertex(i).ToPoint())) return true;
        }

        // Now check whether the boundaries intersect.  Unfortunately, a
        // latitude-longitude rectangle does not have straight edges -- two edges
        // are curved, and at least one of them is concave.

        for (int i = 0; i < 4; ++i)
        {
            S1Interval edge_lng = S1Interval.FromPointPair(
                cell_ll[i].LngRadians, cell_ll[(i + 1) & 3].LngRadians);
            if (!Lng.Intersects(edge_lng)) continue;

            S2Point a = cell_v[i];
            S2Point b = cell_v[(i + 1) & 3];
            if (edge_lng.Contains(Lng.Lo))
            {
                if (IntersectsLngEdge(a, b, Lat, Lng.Lo)) return true;
            }
            if (edge_lng.Contains(Lng.Hi))
            {
                if (IntersectsLngEdge(a, b, Lat, Lng.Hi)) return true;
            }
            if (IntersectsLatEdge(a, b, Lat.Lo, Lng)) return true;
            if (IntersectsLatEdge(a, b, Lat.Hi, Lng)) return true;
        }
        return false;
    }

    // Returns true if and only if the interior of this rectangle intersects
    // any point (including the boundary) of the given other rectangle.
    public bool InteriorIntersects(S2LatLngRect other)
    {
        return (Lat.InteriorIntersects(other.Lat) &&
                Lng.InteriorIntersects(other.Lng));
    }

    // Returns true if the boundary of this rectangle intersects the given
    // geodesic edge (v0, v1).
    public bool BoundaryIntersects(S2Point v0, S2Point v1)
    {
        if (IsEmpty()) return false;
        if (!Lng.IsFull())
        {
            if (IntersectsLngEdge(v0, v1, Lat, Lng.Lo)) return true;
            if (IntersectsLngEdge(v0, v1, Lat, Lng.Hi)) return true;
        }
        if (Lat.Lo != -S2.M_PI_2 && IntersectsLatEdge(v0, v1, Lat.Lo, Lng))
        {
            return true;
        }
        if (Lat.Hi != S2.M_PI_2 && IntersectsLatEdge(v0, v1, Lat.Hi, Lng))
        {
            return true;
        }
        return false;
    }

    // Increase the size of the bounding rectangle to include the given point.
    // The rectangle is expanded by the minimum amount possible.  The S2LatLng
    // argument must be normalized.
    public S2LatLngRect AddPoint(S2Point p)
    {
        return AddPoint(new S2LatLng(p));
    }
    public S2LatLngRect AddPoint(S2LatLng ll)
    {
        if (!ll.IsValid()) Debug.WriteLine($"Invalid S2LatLng in S2LatLngRect.AddPoint: {ll}");

        return new S2LatLngRect(
            R1Interval.AddPoint(Lat, ll.LatRadians),
            S1Interval.AddPoint(Lng, ll.LngRadians));
    }

    // Returns a rectangle that has been expanded by margin.lat() on each side in
    // the latitude direction, and by margin.lng() on each side in the longitude
    // direction.  If either margin is negative, then shrinks the rectangle on
    // the corresponding sides instead.  The resulting rectangle may be empty.
    //
    // As noted above, the latitude-longitude space has the topology of a
    // cylinder.  Longitudes "wrap around" at +/-180 degrees, while latitudes
    // are clamped to range [-90, 90].  This means that any expansion (positive
    // or negative) of the full longitude range remains full (since the
    // "rectangle" is actually a continuous band around the cylinder), while
    // expansion of the full latitude range remains full only if the margin is
    // positive.
    //
    // If either the latitude or longitude interval becomes empty after
    // expansion by a negative margin, the result is empty.
    //
    // Note that if an expanded rectangle contains a pole, it may not contain
    // all possible lat/lng representations of that pole (see header above).
    // Use the PolarClosure() method if you do not want this behavior.
    //
    // If you are trying to grow a rectangle by a certain *distance* on the
    // sphere (e.g. 5km), use the ExpandedByDistance() method instead.
    public S2LatLngRect Expanded(S2LatLng margin)
    {
        var lat = Lat.Expanded(margin.LatRadians);
        var lng = Lng.Expanded(margin.LngRadians);
        if (lat.IsEmpty() || lng.IsEmpty()) return Empty;
        return new S2LatLngRect(lat.Intersection(FullLat), lng);
    }

    // If the rectangle does not include either pole, returns it unmodified.
    // Otherwise expands the longitude range to Full() so that the rectangle
    // contains all possible representations of the contained pole(s).
    public S2LatLngRect PolarClosure()
    {
        if (Lat.Lo == -S2.M_PI_2 || Lat.Hi == S2.M_PI_2)
        {
            return new S2LatLngRect(Lat, S1Interval.Full);
        }
        return this;
    }

    // Returns the smallest rectangle containing the union of this rectangle and
    // the given rectangle.
    public S2LatLngRect Union(S2LatLngRect other)
    {
        return new S2LatLngRect(Lat.Union(other.Lat), Lng.Union(other.Lng));
    }

    // Returns the smallest rectangle containing the intersection of this
    // rectangle and the given rectangle.  Note that the region of intersection
    // may consist of two disjoint rectangles, in which case a single rectangle
    // spanning both of them is returned.
    public S2LatLngRect Intersection(S2LatLngRect other)
    {
        var lat = Lat.Intersection(other.Lat);
        var lng = Lng.Intersection(other.Lng);
        if (lat.IsEmpty() || lng.IsEmpty())
        {
            // The lat/lng ranges must either be both empty or both non-empty.
            return Empty;
        }
        return new S2LatLngRect(lat, lng);
    }

    // Expands this rectangle so that it contains all points within the given
    // distance of the boundary, and return the smallest such rectangle.  If the
    // distance is negative, then instead shrinks this rectangle so that it
    // excludes all points within the given absolute distance of the boundary,
    // and returns the largest such rectangle.
    //
    // Unlike Expanded(), this method treats the rectangle as a set of points on
    // the sphere, and measures distances on the sphere.  For example, you can
    // use this method to find a rectangle that contains all points within 5km
    // of a given rectangle.  Because this method uses the topology of the
    // sphere, note the following:
    //
    //  - The full and empty rectangles have no boundary on the sphere.  Any
    //    expansion (positive or negative) of these rectangles leaves them
    //    unchanged.
    //
    //  - Any rectangle that covers the full longitude range does not have an
    //    east or west boundary, therefore no expansion (positive or negative)
    //    will occur in that direction.
    //
    //  - Any rectangle that covers the full longitude range and also includes
    //    a pole will not be expanded or contracted at that pole, because it
    //    does not have a boundary there.
    //
    //  - If a rectangle is within the given distance of a pole, the result will
    //    include the full longitude range (because all longitudes are present
    //    at the poles).
    //
    // Expansion and contraction are defined such that they are inverses whenever
    // possible, i.e.
    //
    //   rect.ExpandedByDistance(x).ExpandedByDistance(-x) == rect
    //
    // (approximately), so long as the first operation does not cause a
    // rectangle boundary to disappear (i.e., the longitude range newly becomes
    // full or empty, or the latitude range expands to include a pole).
    public S2LatLngRect ExpandedByDistance(S1Angle distance)
    {
        if (distance >= S1Angle.Zero)
        {
            // The most straightforward approach is to build a cap centered on each
            // vertex and take the union of all the bounding rectangles (including the
            // original rectangle; this is necessary for very large rectangles).

            // TODO(ericv): Update this code to use an algorithm like the one below.
            var radius = new S1ChordAngle(distance);
            S2LatLngRect r = this;
            for (int k = 0; k < 4; ++k)
            {
                r = r.Union(new S2Cap(Vertex(k).ToPoint(), radius).GetRectBound());
            }
            return r;
        }
        else
        {
            // Shrink the latitude interval unless the latitude interval contains a pole
            // and the longitude interval is full, in which case the rectangle has no
            // boundary at that pole.
            var lo = Lat.Lo <= FullLat.Lo && Lng.IsFull() ? FullLat.Lo : Lat.Lo - distance.Radians;
            var hi = Lat.Hi >= FullLat.Hi && Lng.IsFull() ? FullLat.Hi : Lat.Hi + distance.Radians;
            var lat_result = new R1Interval(lo, hi);

            if (lat_result.IsEmpty())
            {
                return Empty;
            }

            // Maximum absolute value of a latitude in lat_result. At this latitude,
            // the cap occupies the largest longitude interval.
            double max_abs_lat = Math.Max(-lat_result.Lo, lat_result.Hi);

            // Compute the largest longitude interval that the cap occupies. We use the
            // law of sines for spherical triangles. For the details, see the comment in
            // S2Cap.GetRectBound().
            //
            // When sin_a >= sin_c, the cap covers all the latitude.
            double sin_a = Math.Sin(-distance.Radians);
            double sin_c = Math.Cos(max_abs_lat);
            double max_lng_margin = sin_a < sin_c ? Math.Asin(sin_a / sin_c) : S2.M_PI_2;

            var lng_result = Lng.Expanded(-max_lng_margin);
            if (lng_result.IsEmpty())
            {
                return Empty;
            }
            return new S2LatLngRect(lat_result, lng_result);
        }
    }

    // Returns the minimum distance (measured along the surface of the sphere) to
    // the given S2LatLngRect. Both S2LatLngRects must be non-empty.
    public S1Angle GetDistance(S2LatLngRect other)
    {
        var a = this;
        var b = other;
        Debug.Assert(!a.IsEmpty());
        Debug.Assert(!b.IsEmpty());

        // First, handle the trivial cases where the longitude intervals overlap.
        if (a.Lng.Intersects(b.Lng))
        {
            if (a.Lat.Intersects(b.Lat))
                return S1Angle.FromRadians(0);  // Intersection between a and b.

            // We found an overlap in the longitude interval, but not in the latitude
            // interval. This means the shortest path travels along some line of
            // longitude connecting the high-latitude of the lower rect with the
            // low-latitude of the higher rect.
            S1Angle lo, hi;
            if (a.Lat.Lo > b.Lat.Hi)
            {
                lo = b.LatHi();
                hi = a.LatLo();
            }
            else
            {
                lo = a.LatHi();
                hi = b.LatLo();
            }
            return hi - lo;
        }

        // The longitude intervals don't overlap. In this case, the closest points
        // occur somewhere on the pair of longitudinal edges which are nearest in
        // longitude-space.
        S1Angle a_lng, b_lng;
        S1Interval lo_hi = S1Interval.FromPointPair(a.Lng.Lo, b.Lng.Hi);
        S1Interval hi_lo = S1Interval.FromPointPair(a.Lng.Hi, b.Lng.Lo);
        if (lo_hi.GetLength() < hi_lo.GetLength())
        {
            a_lng = a.LngLo();
            b_lng = b.LngHi();
        }
        else
        {
            a_lng = a.LngHi();
            b_lng = b.LngLo();
        }

        // The shortest distance between the two longitudinal segments will include at
        // least one segment endpoint. We could probably narrow this down further to a
        // single point-edge distance by comparing the relative latitudes of the
        // endpoints, but for the sake of clarity, we'll do all four point-edge
        // distance tests.
        var a_lo = new S2LatLng(a.LatLo(), a_lng).ToPoint();
        var a_hi = new S2LatLng(a.LatHi(), a_lng).ToPoint();
        var b_lo = new S2LatLng(b.LatLo(), b_lng).ToPoint();
        var b_hi = new S2LatLng(b.LatHi(), b_lng).ToPoint();

        return new[]
        {
                S2.GetDistance(a_lo, b_lo, b_hi),
                S2.GetDistance(a_hi, b_lo, b_hi),
                S2.GetDistance(b_lo, a_lo, a_hi),
                S2.GetDistance(b_hi, a_lo, a_hi),
            }.Min();
    }

    // Returns the minimum distance (measured along the surface of the sphere)
    // from a given point to the rectangle (both its boundary and its interior).
    // The latlng must be valid.
    public S1Angle GetDistance(S2LatLng p)
    {
        // The algorithm here is the same as in GetDistance(S2LatLngRect), only
        // with simplified calculations.
        S2LatLngRect a = this;
        if (a.IsEmpty()) Debug.WriteLine($"Empty S2LatLngRect in GetDistance: {a}");
        if (!p.IsValid()) Debug.WriteLine($"Invalid S2LatLng in S2LatLngRect.GetDistance: {p}");

        if (a.Lng.Contains(p.LngRadians))
        {
            return S1Angle.FromRadians(Math.Max(0.0, Math.Max(p.LatRadians - a.Lat.Hi,
                                                 a.Lat.Lo - p.LatRadians)));
        }

        var interval = new S1Interval(a.Lng.Hi, a.Lng.ComplementCenter());
        double a_lng;
        if (interval.Contains(p.LngRadians))
        {
            a_lng = a.Lng.Hi;
        }
        else
        {
            a_lng = a.Lng.Lo;
        }
        S2Point lo = S2LatLng.FromRadians(a.Lat.Lo, a_lng).ToPoint();
        S2Point hi = S2LatLng.FromRadians(a.Lat.Hi, a_lng).ToPoint();
        return S2.GetDistance(p.ToPoint(), lo, hi);
    }

    // Returns the (directed or undirected) Hausdorff distance (measured along the
    // surface of the sphere) to the given S2LatLngRect. The directed Hausdorff
    // distance from rectangle A to rectangle B is given by
    //     h(A, B) = max_{p in A} min_{q in B} d(p, q).
    // The Hausdorff distance between rectangle A and rectangle B is given by
    //     H(A, B) = max{h(A, B), h(B, A)}.
    public S1Angle GetDirectedHausdorffDistance(S2LatLngRect other)
    {
        if (IsEmpty())
        {
            return S1Angle.FromRadians(0);
        }
        if (other.IsEmpty())
        {
            return S1Angle.FromRadians(Math.PI);  // maximum possible distance on S2
        }

        double lng_distance = Lng.GetDirectedHausdorffDistance(other.Lng);
        Debug.Assert(lng_distance >= 0);
        return GetDirectedHausdorffDistance(lng_distance, Lat, other.Lat);
    }
    public S1Angle GetHausdorffDistance(S2LatLngRect other)
    {
        return S1Angle.Max(
            GetDirectedHausdorffDistance(other),
            other.GetDirectedHausdorffDistance(this));
    }

    // Returns true if the latitude and longitude intervals of the two rectangles
    // are the same up to the given tolerance (see r1interval.h and s1interval.h
    // for details).
    public bool ApproxEquals(S2LatLngRect other)
    {
        return ApproxEquals(other, S1Angle.FromRadians(S2.DoubleError));
    }
    public bool ApproxEquals(S2LatLngRect other, S1Angle max_error)
    {
        return (Lat.ApproxEquals(other.Lat, max_error.Radians) &&
                Lng.ApproxEquals(other.Lng, max_error.Radians));
    }

    // ApproxEquals() with separate tolerances for latitude and longitude.
    public bool ApproxEquals(S2LatLngRect other, S2LatLng max_error)
    {
        return (Lat.ApproxEquals(other.Lat, max_error.LatRadians) &&
                Lng.ApproxEquals(other.Lng, max_error.LngRadians));
    }

    // Returns true if the edge AB intersects the given edge of constant
    // longitude.
    public static bool IntersectsLngEdge(S2Point a, S2Point b, R1Interval lat, double lng)
    {
        // Return true if the segment AB intersects the given edge of constant
        // longitude.  The nice thing about edges of constant longitude is that
        // they are straight lines on the sphere (geodesics).

        var p1 = S2LatLng.FromRadians(lat.Lo, lng).ToPoint();
        var p2 = S2LatLng.FromRadians(lat.Hi, lng).ToPoint();
        return S2.CrossingSign(a, b, p1, p2) > 0;
    }

    // Returns true if the edge AB intersects the given edge of constant
    // latitude.  Requires the vectors to have unit length.
    public static bool IntersectsLatEdge(S2Point a, S2Point b, double lat, S1Interval lng)
    {
        // Return true if the segment AB intersects the given edge of constant
        // latitude.  Unfortunately, lines of constant latitude are curves on
        // the sphere.  They can intersect a straight edge in 0, 1, or 2 points.
        Debug.Assert(a.IsUnitLength());
        Debug.Assert(b.IsUnitLength());

        // First, compute the normal to the plane AB that points vaguely north.
        var z = S2.RobustCrossProd(a, b).Normalize();
        if (z[2] < 0) z = -z;

        // Extend this to an orthonormal frame (x,y,z) where x is the direction
        // where the great circle through AB achieves its maximium latitude.
        var y = S2.RobustCrossProd(z, new S2Point(0, 0, 1)).Normalize();
        var x = y.CrossProd(z);
        Debug.Assert(x.IsUnitLength());
        Debug.Assert(x[2] >= 0);

        // Compute the angle "theta" from the x-axis (in the x-y plane defined
        // above) where the great circle intersects the given line of latitude.
        double sin_lat = Math.Sin(lat);
        if (Math.Abs(sin_lat) >= x[2])
        {
            return false;  // The great circle does not reach the given latitude.
        }
        Debug.Assert(x[2] > 0);
        double cos_theta = sin_lat / x[2];
        double sin_theta = Math.Sqrt(1 - cos_theta * cos_theta);
        double theta = Math.Atan2(sin_theta, cos_theta);

        // The candidate intersection points are located +/- theta in the x-y
        // plane.  For an intersection to be valid, we need to check that the
        // intersection point is contained in the interior of the edge AB and
        // also that it is contained within the given longitude interval "lng".

        // Compute the range of theta values spanned by the edge AB.
        S1Interval ab_theta = S1Interval.FromPointPair(
            Math.Atan2(a.DotProd(y), a.DotProd(x)),
            Math.Atan2(b.DotProd(y), b.DotProd(x)));

        if (ab_theta.Contains(theta))
        {
            // Check if the intersection point is also in the given "lng" interval.
            S2Point isect = x * cos_theta + y * sin_theta;
            if (lng.Contains(Math.Atan2(isect[1], isect[0]))) return true;
        }
        if (ab_theta.Contains(-theta))
        {
            // Check if the intersection point is also in the given "lng" interval.
            S2Point isect = x * cos_theta - y * sin_theta;
            if (lng.Contains(Math.Atan2(isect[1], isect[0]))) return true;
        }
        return false;
    }

    // Return the directed Hausdorff distance from one longitudinal edge spanning
    // latitude range 'a_lat' to the other longitudinal edge spanning latitude
    // range 'b_lat', with their longitudinal difference given by 'lng_diff'.
    private static S1Angle GetDirectedHausdorffDistance(double lng_diff, R1Interval a, R1Interval b)
    {
        // By symmetry, we can assume a's longtitude is 0 and b's longtitude is
        // lng_diff. Call b's two endpoints b_lo and b_hi. Let H be the hemisphere
        // containing a and delimited by the longitude line of b. The Voronoi diagram
        // of b on H has three edges (portions of great circles) all orthogonal to b
        // and meeting at b_lo cross b_hi.
        // E1: (b_lo, b_lo cross b_hi)
        // E2: (b_hi, b_lo cross b_hi)
        // E3: (-b_mid, b_lo cross b_hi), where b_mid is the midpoint of b
        //
        // They subdivide H into three Voronoi regions. Depending on how longitude 0
        // (which contains edge a) intersects these regions, we distinguish two cases:
        // Case 1: it intersects three regions. This occurs when lng_diff <= S2Constants.M_PI_2.
        // Case 2: it intersects only two regions. This occurs when lng_diff > S2Constants.M_PI_2.
        //
        // In the first case, the directed Hausdorff distance to edge b can only be
        // realized by the following points on a:
        // A1: two endpoints of a.
        // A2: intersection of a with the equator, if b also intersects the equator.
        //
        // In the second case, the directed Hausdorff distance to edge b can only be
        // realized by the following points on a:
        // B1: two endpoints of a.
        // B2: intersection of a with E3
        // B3: farthest point from b_lo to the interior of D, and farthest point from
        //     b_hi to the interior of U, if any, where D (resp. U) is the portion
        //     of edge a below (resp. above) the intersection point from B2.

        Debug.Assert(lng_diff >= 0);
        Debug.Assert(lng_diff <= Math.PI);

        if (lng_diff == 0)
        {
            return S1Angle.FromRadians(a.GetDirectedHausdorffDistance(b));
        }

        // Assumed longtitude of b.
        double b_lng = lng_diff;
        // Two endpoints of b.
        S2Point b_lo = S2LatLng.FromRadians(b.Lo, b_lng).ToPoint();
        S2Point b_hi = S2LatLng.FromRadians(b.Hi, b_lng).ToPoint();

        // Handling of each case outlined at the top of the function starts here.


        // Cases A1 and B1.
        S2Point a_lo = S2LatLng.FromRadians(a.Lo, 0).ToPoint();
        S2Point a_hi = S2LatLng.FromRadians(a.Hi, 0).ToPoint();
        var max_distance = S2.GetDistance(a_lo, b_lo, b_hi);
        max_distance = S1Angle.Max(max_distance, S2.GetDistance(a_hi, b_lo, b_hi));

        if (lng_diff <= S2.M_PI_2)
        {
            // Case A2.
            if (a.Contains(0) && b.Contains(0))
            {
                max_distance = S1Angle.Max(max_distance, S1Angle.FromRadians(lng_diff));
            }
        }
        else
        {
            // Case B2.
            S2Point p = GetBisectorIntersection(b, b_lng);
            double p_lat = S2LatLng.Latitude(p).Radians;
            if (a.Contains(p_lat))
            {
                max_distance = S1Angle.Max(max_distance, new S1Angle(p, b_lo));
            }

            // Case B3.
            if (p_lat > a.Lo)
            {
                max_distance = S1Angle.Max(max_distance, GetInteriorMaxDistance(
                    new R1Interval(a.Lo, Math.Min(p_lat, a.Hi)), b_lo));
            }
            if (p_lat < a.Hi)
            {
                max_distance = S1Angle.Max(max_distance, GetInteriorMaxDistance(
                    new R1Interval(Math.Max(p_lat, a.Lo), a.Hi), b_hi));
            }
        }

        return max_distance;
    }

    // Return max distance from a point b to the segment spanning latitude range
    // a_lat on longitude 0, if the max occurs in the interior of a_lat. Otherwise
    // return -1.
    private static S1Angle GetInteriorMaxDistance(R1Interval a_lat, S2Point b)
    {
        // Longitude 0 is in the y=0 plane. b.x() >= 0 implies that the maximum
        // does not occur in the interior of a_lat.
        if (a_lat.IsEmpty() || b.X >= 0) return S1Angle.FromRadians(-1);

        // Project b to the y=0 plane. The antipodal of the normalized projection is
        // the point at which the maxium distance from b occurs, if it is contained
        // in a_lat.
        var intersection_point = new S2Point(-b.X, 0, -b.Z).Normalize();
        if (a_lat.InteriorContains(
            S2LatLng.Latitude(intersection_point).Radians))
        {
            return new S1Angle(b, intersection_point);
        }
        else
        {
            return S1Angle.FromRadians(-1);
        }
    }

    // Return the intersection of longitude 0 with the bisector of an edge
    // on longitude 'lng' and spanning latitude range 'lat'.
    private static S2Point GetBisectorIntersection(R1Interval lat, double lng)
    {
        lng = Math.Abs(lng);
        double lat_center = lat.GetCenter();
        // A vector orthogonal to the bisector of the given longitudinal edge.
        S2LatLng ortho_bisector;
        if (lat_center >= 0)
        {
            ortho_bisector = S2LatLng.FromRadians(lat_center - S2.M_PI_2, lng);
        }
        else
        {
            ortho_bisector = S2LatLng.FromRadians(-lat_center - S2.M_PI_2, lng - Math.PI);
        }
        // A vector orthogonal to longitude 0.
        var ortho_lng = new S2Point(0, -1, 0);
        return S2.RobustCrossProd(ortho_lng, ortho_bisector.ToPoint());
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        // We consider two possible bounding caps, one whose axis passes
        // through the center of the lat-long rectangle and one whose axis
        // is the north or south pole.  We return the smaller of the two caps.

        if (IsEmpty()) return S2Cap.Empty;

        double pole_z, pole_angle;
        if (Lat.Lo + Lat.Hi < 0)
        {
            // South pole axis yields smaller cap.
            pole_z = -1;
            pole_angle = S2.M_PI_2 + Lat.Hi;
        }
        else
        {
            pole_z = 1;
            pole_angle = S2.M_PI_2 - Lat.Lo;
        }
        // Ensure that the bounding cap is conservative taking into account errors
        // in the arithmetic above and the S1Angle/S1ChordAngle conversion.
        S2Cap pole_cap = new(new S2Point(0, 0, pole_z),
            S1Angle.FromRadians((1 + 2 * S2.DoubleEpsilon) * pole_angle));

        // For bounding rectangles that span 180 degrees or less in longitude, the
        // maximum cap size is achieved at one of the rectangle vertices.  For
        // rectangles that are larger than 180 degrees, we punt and always return a
        // bounding cap centered at one of the two poles.
        if (Lng.GetLength() < S2.M_2_PI)
        {
            var mid_cap = new S2Cap(Center().ToPoint(), S1Angle.Zero);
            for (int k = 0; k < 4; ++k)
            {
                mid_cap = mid_cap.AddPoint(Vertex(k).ToPoint());
            }
            if (mid_cap.Height() < pole_cap.Height())
                return mid_cap;
        }
        return pole_cap;
    }

    public S2LatLngRect GetRectBound() => this;

    public bool Contains(S2Cell cell)
    {
        // A latitude-longitude rectangle contains a cell if and only if it contains
        // the cell's bounding rectangle.  This test is exact from a mathematical
        // point of view, assuming that the bounds returned by S2Cell.GetRectBound()
        // are tight.  However, note that there can be a loss of precision when
        // converting between representations -- for example, if an S2Cell is
        // converted to a polygon, the polygon's bounding rectangle may not contain
        // the cell's bounding rectangle.  This has some slightly unexpected side
        // effects; for instance, if one creates an S2Polygon from an S2Cell, the
        // polygon will contain the cell, but the polygon's bounding box will not.
        return Contains(cell.GetRectBound());
    }

    // This test is cheap but is NOT exact.  Use Intersects() if you want a more
    // accurate and more expensive test.  Note that when this method is used by
    // an S2RegionCoverer, the accuracy isn't all that important since if a cell
    // may intersect the region then it is subdivided, and the accuracy of this
    // method goes up as the cells get smaller.
    public bool MayIntersect(S2Cell cell)
    {
        // This test is cheap but is NOT exact (see s2latlng_rect.h).
        return Intersects(cell.GetRectBound());
    }

    // The point 'p' does not need to be normalized.
    public bool Contains(S2Point p) => Contains(new S2LatLng(p));

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2LatLngRect(Lat, Lng);

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2LatLngRect to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        encoder.Ensure(40);  // sufficient

        encoder.Put8(S2.kCurrentLosslessEncodingVersionNumber);
        encoder.PutDouble(Lat.Lo);
        encoder.PutDouble(Lat.Hi);
        encoder.PutDouble(Lng.Lo);
        encoder.PutDouble(Lng.Hi);

        Debug.Assert(encoder.Avail() >= 0);
    }

    // Decodes an S2LatLngRect encoded with Encode().  Returns true on success.
    public static (bool, S2LatLngRect) Decode(Decoder decoder)
    {
        if (decoder.Avail() < sizeof(byte) + 4 * sizeof(double))
            return (false, default);

        byte version = decoder.Get8();
        if (version > S2.kCurrentLosslessEncodingVersionNumber)
            return (false, default);

        double lat_lo = decoder.GetDouble();
        double lat_hi = decoder.GetDouble();
        var lat = new R1Interval(lat_lo, lat_hi);
        double lng_lo = decoder.GetDouble();
        double lng_hi = decoder.GetDouble();
        var lng = new S1Interval(lng_lo, lng_hi);
        var result = new S2LatLngRect(lat, lng);
        if (!result.IsValid())
        {
#if s2debug
            Debug.WriteLine($"Invalid result in S2LatLngRect.Decode: {result}");
#endif
            return (false, default);
        }

        return (true, result);
    }

    #endregion

    #region Object

    public override string ToString()
    {
        return $"[Lo{Lo()}, Hi{Hi()}]";
    }

    #endregion
} 
