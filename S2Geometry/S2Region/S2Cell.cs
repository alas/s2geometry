// An S2Cell is an S2Region object that represents a cell.  Unlike S2CellIds,
// it supports efficient containment and intersection tests.  However, it is
// also a more expensive representation (currently 48 bytes rather than 8).

// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator, however it is
// not a "plain old datatype" (POD) because it has virtual functions.

namespace S2Geometry;

public readonly record struct S2Cell : IS2Region<S2Cell>, IComparable<S2Cell>, IDecoder<S2Cell>
{
    #region Fields, Constants

    public readonly S2CellId Id { get; init; }
    public readonly int Face { get; init; }
    public readonly int Level { get; init; }
    public readonly int Orientation { get; init; }
    // The bounds of this cell in (u,v)-space.
    public readonly R2Rect BoundUV { get; init; }

    // Since S2Cells are copied by value, the following assertion is a reminder
    // not to add fields unnecessarily.  An S2Cell currently consists of 43 data
    // bytes, one vtable pointer, plus alignment overhead.  This works out to 48
    // bytes on 32 bit architectures and 56 bytes on 64 bit architectures.
    //
    // The expression below rounds up (43 + sizeof(byte[])) to the nearest
    // multiple of sizeof(byte[]).
    // static_assert(sizeof(S2Cell) <= ((43+2*sizeof(byte[])-1) & -sizeof(byte[])), "S2Cell is getting bloated");

    // The 4 cells around the equator extend to +/-45 degrees latitude at the
    // midpoints of their top and bottom edges.  The two cells covering the
    // poles extend down to +/-35.26 degrees at their vertices.  The maximum
    // error in this calculation is 0.5 * S2Constants.DoubleEpsilon.
    private static readonly double kPoleMinLat = Math.Asin(Math.Sqrt(1.0 / 3)) - 0.5 * S2.DoubleEpsilon;

    #endregion

    #region Constructors

    // An S2Cell always corresponds to a particular S2CellId.  The other
    // constructors are just convenience methods.
    public S2Cell(S2CellId id)
    {
        Id = id;
        var ij = new int[2];
        Face = (sbyte)id.ToFaceIJOrientation(out ij[0], out ij[1], out int orientation, true);
        Orientation = orientation;  // Compress int to a byte.
        Level = id.Level();
        BoundUV = S2CellId.IJLevelToBoundUV(ij, Level);
    }

    // Convenience constructors.  The S2LatLng must be normalized.
    public S2Cell(S2Point p) : this(new S2CellId(p)) { }
    public S2Cell(S2LatLng ll) : this(new S2CellId(ll)) { }

    private S2Cell(S2CellId id, int face, int level, int orientation, R2Rect uv)
    {
        Id = id;
        Face = face;
        Level = level;
        Orientation = orientation;
        BoundUV = uv;
    }

    #endregion

    #region Factories

    // Returns the cell corresponding to the given S2 cube face.
    public static S2Cell FromFace(int face) => new(S2CellId.FromFace(face));

    // Returns a cell given its face (range 0..5), Hilbert curve position within
    // that face (an unsigned integer with S2CellId.kPosBits bits), and level
    // (range 0..kMaxLevel).  The given position will be modified to correspond
    // to the Hilbert curve position at the center of the returned cell.  This
    // is a static function rather than a constructor in order to indicate what
    // the arguments represent.
    public static S2Cell FromFacePosLevel(int face, UInt64 pos, int level)
        => new(S2CellId.FromFacePosLevel(face, pos, level));

    #endregion

    #region S2Cell

    public bool IsLeaf() => Level == S2.kMaxCellLevel;

    // These are equivalent to the S2CellId methods, but have a more efficient
    // implementation since the level has been precomputed.
    public int SizeIJ() => S2CellId.SizeIJ(Level);

    public double SizeST() => S2CellId.SizeST(Level);

    // Returns the k-th vertex of the cell (k = 0,1,2,3).  Vertices are returned
    // in CCW order (lower left, lower right, upper right, upper left in the UV
    // plane).  The points returned by GetVertexRaw are not normalized.
    // For convenience, the argument is reduced modulo 4 to the range [0..3].
    public S2Point Vertex(int k) => VertexRaw(k).Normalize();
    public S2Point VertexRaw(int k) => S2.FaceUVtoXYZ(Face, BoundUV.GetVertex(k));

    // Returns the inward-facing normal of the great circle passing through the
    // edge from vertex k to vertex k+1 (mod 4).  The normals returned by
    // GetEdgeRaw are not necessarily unit length.  For convenience, the
    // argument is reduced modulo 4 to the range [0..3].
    public S2Point Edge(int k) => EdgeRaw(k).Normalize();
    public S2Point EdgeRaw(int k) => (k & 3) switch
    {
        0 => S2.GetVNorm(Face, BoundUV[1][0]), // Bottom
        1 => S2.GetUNorm(Face, BoundUV[0][1]), // Right
        2 => -S2.GetVNorm(Face, BoundUV[1][1]),// Top
        _ => -S2.GetUNorm(Face, BoundUV[0][0]),// Left
    };

    // If this is not a leaf cell, sets children[0..3] to the four children of
    // this cell (in traversal order) and return true.  Otherwise returns false.
    // This method is equivalent to the following:
    //
    // for (pos=0, id=ChildBegin(); id != ChildEnd(); id = id.Next, ++pos)
    //   children[pos] = S2Cell(id);
    //
    // except that it is more than two times faster.
    public bool Subdivide(S2Cell[] children)
    {
        //Debug.Assert(children.Length == 4);

        // This function is equivalent to just iterating over the child cell ids
        // and calling the S2Cell constructor, but it is about 2.5 times faster.

        if (Id.IsLeaf()) return false;

        // Compute the cell midpoint in uv-space.
        var uv_mid = Id.CenterUV();

        // Create four children with the appropriate bounds.
        var id = Id.ChildBegin();
        for (int pos = 0; pos < 4; ++pos, id = id.Next())
        {
            // We want to split the cell in half in "u" and "v".  To decide which
            // side to set equal to the midpoint value, we look at cell's (i,j)
            // position within its parent.  The index for "i" is in bit 1 of ij.
            int ij = S2.kPosToIJ[Orientation][pos];
            int i = ij >> 1;
            int j = ij & 1;
            var tmp = new double[2];
            tmp[i] = BoundUV[0][i];
            tmp[1 - i] = uv_mid[0];
            var x = new R1Interval(tmp[0], tmp[1]);
            tmp[j] = BoundUV[1][j];
            tmp[1 - j] = uv_mid[1];
            var y = new R1Interval(tmp[0], tmp[1]);
            var uv = new R2Rect(x, y);
            children[pos] = new S2Cell(id, Face, Level + 1, Orientation ^ S2.kPosToOrientation[pos], uv);
        }
        return true;
    }

    // Returns the direction vector corresponding to the center in (s,t)-space of
    // the given cell.  This is the point at which the cell is divided into four
    // subcells; it is not necessarily the centroid of the cell in (u,v)-space
    // or (x,y,z)-space.  The point returned by GetCenterRaw is not necessarily
    // unit length.
    public S2Point Center() => CenterRaw().Normalize();
    public S2Point CenterRaw() => Id.ToPointRaw();

    // Returns the average area for cells at the given level.
    public static double AverageArea(int level) => S2.kAvgArea.GetValue(level);

    // Returns the average area of cells at this level in steradians.  This is
    // accurate to within a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is
    // extremely cheap to compute.
    public double AverageArea() => AverageArea(Level);

    // Returns the approximate area of this cell in steradians.  This method is
    // accurate to within 3% percent for all cell sizes and accurate to within
    // 0.1% for cells at level 5 or higher (i.e. squares 350km to a side or
    // smaller on the Earth's surface).  It is moderately cheap to compute.
    public double ApproxArea()
    {
        // All cells at the first two levels have the same area.
        if (Level < 2) return AverageArea(Level);

        // First, compute the approximate area of the cell when projected
        // perpendicular to its normal.  The cross product of its diagonals gives
        // the normal, and the length of the normal is twice the projected area.
        double flat_area = 0.5 * (Vertex(2) - Vertex(0)).
                           CrossProd(Vertex(3) - Vertex(1)).Norm();

        // Now, compensate for the curvature of the cell surface by pretending
        // that the cell is shaped like a spherical cap.  The ratio of the
        // area of a spherical cap to the area of its projected disc turns out
        // to be 2 / (1 + Math.Sqrt(1 - r*r)) where "r" is the radius of the disc.
        // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
        // Here we set Pi*r*r == flat_area to find the equivalent disc.
        return flat_area * 2 / (1 + Math.Sqrt(1 - Math.Min(S2.M_1_PI * flat_area, 1.0)));
    }

    // Returns the area of this cell as accurately as possible.  This method is
    // more expensive but it is accurate to 6 digits of precision even for leaf
    // cells (whose area is approximately 1e-18).
    public double ExactArea()
    {
        // There is a straightforward mathematical formula for the exact surface
        // area (based on 4 calls to asin), but as the cell size gets small this
        // formula has too much cancellation error.  So instead we compute the area
        // as the sum of two triangles (which is very accurate at all cell levels).
        S2Point v0 = Vertex(0);
        S2Point v1 = Vertex(1);
        S2Point v2 = Vertex(2);
        S2Point v3 = Vertex(3);
        return S2.Area(v0, v1, v2) + S2.Area(v0, v2, v3);
    }

    // Returns the distance from the cell to the given point.  Returns zero if
    // the point is inside the cell.
    public S1ChordAngle Distance(S2Point target)
    {
        return DistanceInternal(target, true /*to_interior*/);
    }

    // Return the distance from the cell boundary to the given point.
    public S1ChordAngle BoundaryDistance(S2Point target)
    {
        return DistanceInternal(target, false /*to_interior*/);
    }

    // Returns the maximum distance from the cell (including its interior) to the
    // given point.
    public S1ChordAngle MaxDistance(S2Point target)
    {
        // First check the 4 cell vertices.  If all are within the hemisphere
        // centered around target, the max distance will be to one of these vertices.
        S2Point target_uvw = S2.FaceXYZtoUVW(Face, target);
        S1ChordAngle max_dist = new[]
        {
                VertexChordDist(target_uvw, 0, 0),
                VertexChordDist(target_uvw, 1, 0),
                VertexChordDist(target_uvw, 0, 1),
                VertexChordDist(target_uvw, 1, 1),
            }.Max();

        if (max_dist <= S1ChordAngle.Right)
        {
            return max_dist;
        }

        // Otherwise, find the minimum distance d_min to the antipodal point and the
        // maximum distance will be Pi - d_min.
        return S1ChordAngle.Straight - Distance(-target);
    }

    // Returns the minimum distance from the cell to the given edge AB.  Returns
    // zero if the edge intersects the cell interior.
    public S1ChordAngle Distance(S2Point a, S2Point b)
    {
        // Possible optimizations:
        //  - Currently the (cell vertex, edge endpoint) distances are computed
        //    twice each, and the length of AB is computed 4 times.
        //  - To fix this, refactor GetDistance(target) so that it skips calculating
        //    the distance to each cell vertex.  Instead, compute the cell vertices
        //    and distances in this function, and add a low-level UpdateMinDistance
        //    that allows the XA, XB, and AB distances to be passed in.
        //  - It might also be more efficient to do all calculations in UVW-space,
        //    since this would involve transforming 2 points rather than 4.

        // First, check the minimum distance to the edge endpoints A and B.
        // (This also detects whether either endpoint is inside the cell.)
        var min_dist = S1ChordAngle.Min(Distance(a), Distance(b));
        if (min_dist == S1ChordAngle.Zero) return min_dist;

        // Otherwise, check whether the edge crosses the cell boundary.
        // Note that S2EdgeCrosser needs pointers to vertices.
        var v = new S2Point[4];
        for (int i = 0; i < 4; ++i)
        {
            v[i] = Vertex(i);
        }
        var crosser = new S2EdgeCrosser(a, b, v[3]);
        for (int i = 0; i < 4; ++i)
        {
            if (crosser.CrossingSign(v[i]) >= 0)
            {
                return S1ChordAngle.Zero;
            }
        }
        // Finally, check whether the minimum distance occurs between a cell vertex
        // and the interior of the edge AB.  (Some of this work is redundant, since
        // it also checks the distance to the endpoints A and B again.)
        //
        // Note that we don't need to check the distance from the interior of AB to
        // the interior of a cell edge, because the only way that this distance can
        // be minimal is if the two edges cross (already checked above).
        for (int i = 0; i < 4; ++i)
        {
            S2.UpdateMinDistance(v[i], a, b, ref min_dist);
        }
        return min_dist;
    }

    // Returns the maximum distance from the cell (including its interior) to the
    // given edge AB.
    public S1ChordAngle MaxDistance(S2Point a, S2Point b)
    {
        // If the maximum distance from both endpoints to the cell is less than Pi/2
        // then the maximum distance from the edge to the cell is the maximum of the
        // two endpoint distances.
        var max_dist = S1ChordAngle.Max(MaxDistance(a), MaxDistance(b));
        if (max_dist <= S1ChordAngle.Right)
        {
            return max_dist;
        }

        return S1ChordAngle.Straight - Distance(-a, -b);
    }

    // Returns the distance from the cell to the given cell.  Returns zero if
    // one cell contains the other.
    public S1ChordAngle Distance(S2Cell target)
    {
        // If the cells intersect, the distance is zero.  We use the (u,v) ranges
        // rather S2CellId.intersects() so that cells that share a partial edge or
        // corner are considered to intersect.
        if (Face == target.Face && BoundUV.Intersects(target.BoundUV))
        {
            return S1ChordAngle.Zero;
        }

        // Otherwise, the minimum distance always occurs between a vertex of one
        // cell and an edge of the other cell (including the edge endpoints).  This
        // represents a total of 32 possible (vertex, edge) pairs.
        //
        // TODO(ericv): This could be optimized to be at least 5x faster by pruning
        // the set of possible closest vertex/edge pairs using the faces and (u,v)
        // ranges of both cells.
        var va = new S2Point[4];
        var vb = new S2Point[4];
        for (int i = 0; i < 4; ++i)
        {
            va[i] = Vertex(i);
            vb[i] = target.Vertex(i);
        }
        var min_dist = S1ChordAngle.Infinity;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                S2.UpdateMinDistance(va[i], vb[j], vb[(j + 1) & 3], ref min_dist);
                S2.UpdateMinDistance(vb[i], va[j], va[(j + 1) & 3], ref min_dist);
            }
        }
        return min_dist;
    }

    // Returns the maximum distance from the cell (including its interior) to the
    // given target cell.
    public S1ChordAngle MaxDistance(S2Cell target)
    {
        // Need to check the antipodal target for intersection with the cell. If it
        // intersects, the distance is S1ChordAngle.Straight.
        if (Face == OppositeFace(target.Face) &&
            BoundUV.Intersects(OppositeUV(target.BoundUV)))
        {
            return S1ChordAngle.Straight;
        }

        // Otherwise, the maximum distance always occurs between a vertex of one
        // cell and an edge of the other cell (including the edge endpoints).  This
        // represents a total of 32 possible (vertex, edge) pairs.
        //
        // TODO(user): When the maximum distance is at most Pi/2, the maximum is
        // always attained between a pair of vertices, and this could be made much
        // faster by testing each vertex pair once rather than the current 4 times.
        var va = new S2Point[4];
        var vb = new S2Point[4];
        for (int i = 0; i < 4; ++i)
        {
            va[i] = Vertex(i);
            vb[i] = target.Vertex(i);
        }
        S1ChordAngle max_dist = S1ChordAngle.Negative;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                S2.UpdateMaxDistance(va[i], vb[j], vb[(j + 1) & 3], ref max_dist);
                S2.UpdateMaxDistance(vb[i], va[j], va[(j + 1) & 3], ref max_dist);
            }
        }
        return max_dist;
    }
    private static int OppositeFace(int face) => face >= 3 ? face - 3 : face + 3;

    // The antipodal UV is the transpose of the original UV, interpreted within
    // the opposite face.
    private static R2Rect OppositeUV(R2Rect uv) => new(uv[1], uv[0]);

    // Returns the latitude or longitude of the cell vertex given by (i,j),
    // where "i" and "j" are either 0 or 1.
    private double Latitude(int i, int j)
    {
        var p = S2.FaceUVtoXYZ(Face, BoundUV[0][i], BoundUV[1][j]);
        return S2LatLng.Latitude(p).Radians;
    }
    private double Longitude(int i, int j)
    {
        S2Point p = S2.FaceUVtoXYZ(Face, BoundUV[0][i], BoundUV[1][j]);
        return S2LatLng.Longitude(p).Radians;
    }

    // Return the squared chord distance from point P to corner vertex (i,j).
    private S1ChordAngle VertexChordDist(S2Point p, int i, int j)
    {
        var vertex = new S2Point(BoundUV[0][i], BoundUV[1][j], 1).Normalize();
        return new S1ChordAngle(p, vertex);
    }

    // Given a point P and either the lower or upper edge of the S2Cell (specified
    // by setting "v_end" to 0 or 1 respectively), return true if P is closer to
    // the interior of that edge than it is to either endpoint.
    private bool UEdgeIsClosest(S2Point target, int v_end)
    {
        var u0 = BoundUV[0][0];
        var u1 = BoundUV[0][1];
        var v = BoundUV[1][v_end];
        // These are the normals to the planes that are perpendicular to the edge
        // and pass through one of its two endpoints.
        var dir0 = new S2Point(v * v + 1, -u0 * v, -u0);
        var dir1 = new S2Point(v * v + 1, -u1 * v, -u1);
        return target.DotProd(dir0) > 0 && target.DotProd(dir1) < 0;
    }

    // Given a point P and either the left or right edge of the S2Cell (specified
    // by setting "u_end" to 0 or 1 respectively), return true if P is closer to
    // the interior of that edge than it is to either endpoint.
    private bool VEdgeIsClosest(S2Point target, int u_end)
    {
        var v0 = BoundUV[1][0];
        var v1 = BoundUV[1][1];
        var u = BoundUV[0][u_end];
        // See comments above.
        var dir0 = new S2Point(-u * v0, u * u + 1, -v0);
        var dir1 = new S2Point(-u * v1, u * u + 1, -v1);
        return target.DotProd(dir0) > 0 && target.DotProd(dir1) < 0;
    }

    // Returns the distance from the given point to the interior of the cell if
    // "to_interior" is true, and to the boundary of the cell otherwise.
    private S1ChordAngle DistanceInternal(S2Point target_xyz, bool to_interior)
    {
        // All calculations are done in the (u,v,w) coordinates of this cell's face.
        S2Point target = S2.FaceXYZtoUVW(Face, target_xyz);

        // Compute dot products with all four upward or rightward-facing edge
        // normals.  "dirIJ" is the dot product for the edge corresponding to axis
        // I, endpoint J.  For example, dir01 is the right edge of the S2Cell
        // (corresponding to the upper endpoint of the u-axis).
        double dir00 = target[0] - target[2] * BoundUV[0][0];
        double dir01 = target[0] - target[2] * BoundUV[0][1];
        double dir10 = target[1] - target[2] * BoundUV[1][0];
        double dir11 = target[1] - target[2] * BoundUV[1][1];
        bool inside = true;
        if (dir00 < 0)
        {
            inside = false;  // Target is to the left of the cell
            if (VEdgeIsClosest(target, 0)) return EdgeDistance(-dir00, BoundUV[0][0]);
        }
        if (dir01 > 0)
        {
            inside = false;  // Target is to the right of the cell
            if (VEdgeIsClosest(target, 1)) return EdgeDistance(dir01, BoundUV[0][1]);
        }
        if (dir10 < 0)
        {
            inside = false;  // Target is below the cell
            if (UEdgeIsClosest(target, 0)) return EdgeDistance(-dir10, BoundUV[1][0]);
        }
        if (dir11 > 0)
        {
            inside = false;  // Target is above the cell
            if (UEdgeIsClosest(target, 1)) return EdgeDistance(dir11, BoundUV[1][1]);
        }
        if (inside)
        {
            if (to_interior) return S1ChordAngle.Zero;
            // Although you might think of S2Cells as rectangles, they are actually
            // arbitrary quadrilaterals after they are projected onto the sphere.
            // Therefore the simplest approach is just to find the minimum distance to
            // any of the four edges.
            return new[]
                {
                        EdgeDistance(-dir00, BoundUV[0][0]),
                        EdgeDistance(dir01, BoundUV[0][1]),
                        EdgeDistance(-dir10, BoundUV[1][0]),
                        EdgeDistance(dir11, BoundUV[1][1]),
                    }.Min();
        }
        // Otherwise, the closest point is one of the four cell vertices.  Note that
        // it is *not* trivial to narrow down the candidates based on the edge sign
        // tests above, because (1) the edges don't meet at right angles and (2)
        // there are points on the far side of the sphere that are both above *and*
        // below the cell, etc.
        return new[]
            {
                    VertexChordDist(target, 0, 0),
                    VertexChordDist(target, 1, 0),
                    VertexChordDist(target, 0, 1),
                    VertexChordDist(target, 1, 1),
                }.Min();
    }

    // Given the dot product of a point P with the normal of a u- or v-edge at the
    // given coordinate value, return the distance from P to that edge.
    private static S1ChordAngle EdgeDistance(double dirIJ, double uv)
    {
        // Let P by the target point and let R be the closest point on the given
        // edge AB.  The desired distance PR can be expressed as PR^2 = PQ^2 + QR^2
        // where Q is the point P projected onto the plane through the great circle
        // through AB.  We can compute the distance PQ^2 perpendicular to the plane
        // from "dirIJ" (the dot product of the target point P with the edge
        // normal) and the squared length the edge normal (1 + uv**2).
        var pq2 = (dirIJ * dirIJ) / (1 + uv * uv);

        // We can compute the distance QR as (1 - OQ) where O is the sphere origin,
        // and we can compute OQ^2 = 1 - PQ^2 using the Pythagorean theorem.
        // (This calculation loses accuracy as angle POQ approaches Pi/2.)
        var qr = 1 - Math.Sqrt(1 - pq2);
        return S1ChordAngle.FromLength2(pq2 + qr * qr);
    }

    #endregion

    #region S2Region

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public S2Cap GetCapBound()
    {
        // Use the cell center in (u,v)-space as the cap axis.  This vector is
        // very close to GetCenter() and faster to compute.  Neither one of these
        // vectors yields the bounding cap with minimal surface area, but they
        // are both pretty close.
        //
        // It's possible to show that the two vertices that are furthest from
        // the (u,v)-origin never determine the maximum cap size (this is a
        // possible future optimization).

        var center = S2.FaceUVtoXYZ(Face, BoundUV.GetCenter()).Normalize();
        var cap = S2Cap.FromPoint(center);
        for (int k = 0; k < 4; ++k)
        {
            cap = cap.AddPoint(Vertex(k));
        }
        return cap;
    }

    public S2LatLngRect GetRectBound()
    {
        if (Level > 0)
        {
            // Except for cells at level 0, the latitude and longitude extremes are
            // attained at the vertices.  Furthermore, the latitude range is
            // determined by one pair of diagonally opposite vertices and the
            // longitude range is determined by the other pair.
            //
            // We first determine which corner (i,j) of the cell has the largest
            // absolute latitude.  To maximize latitude, we want to find the point in
            // the cell that has the largest absolute z-coordinate and the smallest
            // absolute x- and y-coordinates.  To do this we look at each coordinate
            // (u and v), and determine whether we want to minimize or maximize that
            // coordinate based on the axis direction and the cell's (u,v) quadrant.
            double u = BoundUV[0][0] + BoundUV[0][1];
            double v = BoundUV[1][0] + BoundUV[1][1];
            int i = (S2.GetUAxis(Face)[2] == 0 ? (u < 0) : (u > 0)) ? 1 : 0;
            int j = (S2.GetVAxis(Face)[2] == 0 ? (v < 0) : (v > 0)) ? 1 : 0;
            var lat = R1Interval.FromPointPair(Latitude(i, j),
                                                       Latitude(1 - i, 1 - j));
            var lng = S1Interval.FromPointPair(Longitude(i, 1 - j),
                                                       Longitude(1 - i, j));

            // We grow the bounds slightly to make sure that the bounding rectangle
            // contains S2LatLng(P) for any point P inside the loop L defined by the
            // four *normalized* vertices.  Note that normalization of a vector can
            // change its direction by up to 0.5 * S2Constants.DoubleEpsilon radians, and it is not
            // enough just to add Normalize() calls to the code above because the
            // latitude/longitude ranges are not necessarily determined by diagonally
            // opposite vertex pairs after normalization.
            //
            // We would like to bound the amount by which the latitude/longitude of a
            // contained point P can exceed the bounds computed above.  In the case of
            // longitude, the normalization error can change the direction of rounding
            // leading to a maximum difference in longitude of 2 * S2Constants.DoubleEpsilon.  In
            // the case of latitude, the normalization error can shift the latitude by
            // up to 0.5 * S2Constants.DoubleEpsilon and the other sources of error can cause the
            // two latitudes to differ by up to another 1.5 * S2Constants.DoubleEpsilon, which also
            // leads to a maximum difference of 2 * S2Constants.DoubleEpsilon.
            return new S2LatLngRect(lat, lng)
                .Expanded(S2LatLng.FromRadians(2 * S2.DoubleEpsilon, 2 * S2.DoubleEpsilon))
                .PolarClosure();
        }

        // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
        Debug.Assert(((Face < 3) ? 1 : -1) == S2.GetNorm(Face)[Face % 3]);
        var bound = Face switch
        {
            0 => new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4),
                                  new S1Interval(-S2.M_PI_4, S2.M_PI_4)),
            1 => new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4),
                                  new S1Interval(S2.M_PI_4, 3 * S2.M_PI_4)),
            2 => new S2LatLngRect(new R1Interval(kPoleMinLat, S2.M_PI_2), S1Interval.Full),
            3 => new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4),
                                  new S1Interval(3 * S2.M_PI_4, -3 * S2.M_PI_4)),
            4 => new S2LatLngRect(new R1Interval(-S2.M_PI_4, S2.M_PI_4),
                                  new S1Interval(-3 * S2.M_PI_4, -S2.M_PI_4)),
            _ => new S2LatLngRect(new R1Interval(-S2.M_PI_2, -kPoleMinLat), S1Interval.Full),
        };
        // Finally, we expand the bound to account for the error when a point P is
        // converted to an S2LatLng to test for containment.  (The bound should be
        // large enough so that it contains the computed S2LatLng of any contained
        // point, not just the infinite-precision version.)  We don't need to expand
        // longitude because longitude is calculated via a single call to atan2(),
        // which is guaranteed to be semi-monotonic.  (In fact the Gnu implementation
        // is also correctly rounded, but we don't even need that here.)
        return bound.Expanded(S2LatLng.FromRadians(S2.DoubleEpsilon, 0));
    }

    public bool Contains(S2Cell cell)
    {
        return Id.Contains(cell.Id);
    }

    public bool MayIntersect(S2Cell cell)
    {
        return Id.Intersects(cell.Id);
    }

    // Returns true if the cell contains the given point "p".  Note that unlike
    // S2Loop/S2Polygon, S2Cells are considered to be closed sets.  This means
    // that points along an S2Cell edge (or at a vertex) belong to the adjacent
    // cell(s) as well.
    //
    // If instead you want every point to be contained by exactly one S2Cell,
    // you will need to convert the S2Cells to S2Loops (which implement point
    // containment this way).
    //
    // The point "p" does not need to be normalized.
    public bool Contains(S2Point p)
    {
        // We can't just call XYZtoFaceUV, because for points that lie on the
        // boundary between two faces (i.e. u or v is +1/-1) we need to return
        // true for both adjacent cells.
        if (!S2.FaceXYZtoUV(Face, p, out R2Point uv)) return false;

        // Expand the (u,v) bound to ensure that
        //
        //   S2Cell(S2CellId(p)).Contains(p)
        //
        // is always true.  To do this, we need to account for the error when
        // converting from (u,v) coordinates to (s,t) coordinates.  At least in the
        // case of S2_QUADRATIC_PROJECTION, the total error is at most S2Constants.DoubleEpsilon.
        return BoundUV.Expanded(S2.DoubleEpsilon).Contains(uv);
    }

    #endregion

    #region IEncoder

    // Appends a serialized representation of the S2Cell to "encoder".
    //
    // REQUIRES: "encoder" uses the default constructor, so that its buffer
    //           can be enlarged as necessary by calling Ensure(int).
    public void Encode(Encoder encoder, CodingHint hint = CodingHint.COMPACT)
    {
        Id.Encode(encoder, hint);
    }

    // Decodes an S2Cell encoded with Encode().  Returns true on success.
    public static (bool, S2Cell) Decode(Decoder decoder)
    {
        var (success, id) = S2CellId.Decode(decoder);
        if (!success)
            return (false, default);

        return (true, new S2Cell(id));
    }

    #endregion

    #region ICustomCloneable

    public object CustomClone() => new S2Cell(Id);

    #endregion

    #region IEquatable

    public bool Equals(S2Cell cell) => Id.Equals(cell.Id);

    public override int GetHashCode() => Id.GetHashCode();

    #endregion

    #region IComparable

    public int CompareTo(S2Cell other)
    {
        if (Id < other.Id) return -1;
        if (Id > other.Id) return 1;
        return 0;
    }

    public static bool operator <(S2Cell x, S2Cell y) => x.Id < y.Id;
    public static bool operator >(S2Cell x, S2Cell y) => x.Id > y.Id;
    public static bool operator <=(S2Cell x, S2Cell y) => x.Id <= y.Id;
    public static bool operator >=(S2Cell x, S2Cell y) => x.Id >= y.Id;

    #endregion
}
