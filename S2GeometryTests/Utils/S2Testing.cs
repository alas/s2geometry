namespace S2Geometry;

using System.Text;

public static class Assert2
{
    public static void Equal(double a, double b, double maxRatio = S2.DoubleError)
    {
        if (a == b) return;

        var full = a != 0 ? a : b;
        var dif = a - b;
        var difRatio = Math.Abs(dif / full);
        Assert.True(difRatio <= maxRatio);
    }

    /// <summary>
    /// Verifies that the two double values are approximately equal, to within 4 ULPs from each other.
    /// 
    /// Replacement for EXPECT_DOUBLE_EQ and ASSERT_DOUBLE_EQ
    /// </summary>
    public static void DoubleEqual(double value1, double value2, long units = 4L)
    {
        var lValue1 = BitConverter.DoubleToInt64Bits(value1);
        var lValue2 = BitConverter.DoubleToInt64Bits(value2);

        // If the signs are different, return false except for +0 and -0.
        if (((ulong)lValue1 >> 63) != ((ulong)lValue2 >> 63))
        {
            Assert.True(value1 == value2);
            return;
        }

        var diff = Math.Abs(lValue1 - lValue2);

        Assert.True(diff <= units);
    }

    /// <summary>
    /// Verifies that the two float values are approximately equal, to within 4 ULPs from each other.
    /// 
    /// Replacement for EXPECT_FLOAT_EQ and ASSERT_FLOAT_EQ
    /// </summary>
    public static void FloatEqual(float value1, float value2)
    {
        var lValue1 = BitConverter.ToInt32(BitConverter.GetBytes(value1), 0);
        var lValue2 = BitConverter.ToInt32(BitConverter.GetBytes(value2), 0);

        // If the signs are different, return false except for +0 and -0.
        if ((lValue1 >> 31) != (lValue2 >> 31))
            Assert.True(value1 == value2);

        var diff = Math.Abs(lValue1 - lValue2);

        Assert.True(diff <= 4L);
    }
    public static void FloatEqual(float value1, double value2) =>
        FloatEqual(value1, (float)value2);

    /// <summary>
    /// Verifies that the difference between value1 and value2 does not exceed the absolute error bound abs_error.
    /// 
    /// Replacement for EXPECT_NEAR
    /// </summary>
    public static void Near(double value1, double value2, double absoluteError = 0.0001) =>
        Assert.True(Math.Abs(value1 - value2) <= absoluteError);
}

// This class defines various static functions that are useful for writing
// unit tests.
public static class S2Testing
{
    public static class StrPoints
    {
        // Given a string where each character "ch" represents a vertex (such as
        // "abac"), returns a vector of S2Points of the form (ch, 0, 0).  Note that
        // these points are not unit length and therefore are not suitable for general
        // use, however they are useful for testing certain functions.
        public static S2PointLoopSpan StrToPoints(string str)
        {
            return Encoding.ASCII.GetBytes(str)
                .Select(t => new S2Point(t, 0, 0))
                .ToList();
        }

        public static string PointsToStr(IEnumerable<S2Point> points)
        {
            return Encoding.ASCII.GetString(
                points.Select(t => (byte)t[0]).ToArray());
        }
    }

    // A deterministically-seeded random number generator.
    //
    // Functions in this class return random numbers that are as good as random()
    // is.  The results are reproducible since the seed is deterministic.  This
    // class is *NOT* thread-safe; it is only intended for testing purposes.
    public static class Random
    {
        // Currently this class is based on random(), therefore it makes no sense to
        // make a copy.
        public static int RandomSeed { get; set; } = 1;
        private static System.Random rnd = new(RandomSeed);

        // Reset the generator state using the given seed.
        public static void Reset(int seed)
        {
            rnd = new System.Random(seed);
        }

        public static int Next(int max)
        {
            lock (rnd) return rnd.Next(max);
        }

        /// <summary>
        /// Return a 64-bit unsigned integer whose lowest "num_bits" are random, and
        /// whose other bits are zero.
        /// </summary>
        public static UInt64 GetBits(int num_bits)
        {
            Assert.True(num_bits >= 0);
            Assert.True(num_bits <= 64);

            // This code uses random(), which returns an integer in the range
            // from 0 to (2^31)-1 inclusive (i.e. all of the lower 31 bits are
            // in play, and the 32nd bit and higher are 0) regardless of whether
            // its return type (long) is larger than 32 bits.  See
            //
            // www.gnu.org/software/libc/manual/html_node/BSD-Random.html#BSD-Random
            //
            // Note that at some point the manual page in linux claimed that the range
            // is 0 to RAND_MAX as defined in stdlib.h.  RAND_MAX however is part only
            // of the ISO rand() interface.  At least as of glibc-2.21, rand() is
            // simply an alias for random().  On other systems, rand() may differ,
            // but random() should always adhere to the behavior specified in BSD.
            const int RAND_BITS = 31;

            UInt64 result = 0;
            for (int bits = 0; bits < num_bits; bits += RAND_BITS)
            {
                result = (result << RAND_BITS) + (ulong)(uint)Next(int.MaxValue);
            }
            if (num_bits < 64)
            {  // Not legal to shift by full bitwidth of type
                result &= (1UL << num_bits) - 1;
            }
            return result;
        }

        /// <summary>
        /// Return a uniformly distributed 64-bit unsigned integer.
        /// </summary>
        public static UInt64 Rand64()
        {
            return GetBits(64);
        }

        /// <summary>
        /// Return a uniformly distributed 32-bit unsigned integer.
        /// </summary>
        public static UInt32 Rand32()
        {
            return (UInt32)GetBits(32);
        }

        /// <summary>
        /// Return a uniformly distributed "double" in the range [0,1).  Note that
        /// the values returned are all multiples of 2**-53, which means that not all
        /// possible values in this range are returned.
        /// </summary>
        public static double RandDouble()
        {
            int NUM_BITS = 53;
            return MathUtils.Ldexp(GetBits(NUM_BITS), -NUM_BITS);
        }

        /// <summary>
        /// Return a uniformly distributed integer in the range [0,n).
        /// </summary>
        /// <param name="n"></param>
        public static Int32 Uniform(Int32 n)
        {
            Assert.True(n > 0);
            return (Int32)(RandDouble() * n);
        }

        // Return a uniformly distributed "double" in the range [min, limit).
        public static double UniformDouble(double min, double limit)
        {
            Assert.True(min < limit);
            return min + RandDouble() * (limit - min);
        }

        // Return true with probability 1 in n.
        public static bool OneIn(Int32 n)
        {
            return Uniform(n) == 0;
        }

        // Skewed: pick "base" uniformly from range [0,max_log] and then
        // return "base" random bits.  The effect is to pick a number in the
        // range [0,2^max_log-1] with bias towards smaller numbers.
        public static Int32 Skewed(int max_log)
        {
            Assert.True(max_log >= 0);
            Assert.True(max_log <= 31);
            Int32 base_ = Uniform(max_log + 1);
            return (int)(GetBits(31) & ((1UL << base_) - 1));
        }
    }

    // Append the vertices of "loop" to "vertices".
    public static void AppendLoopVertices(S2Loop loop, List<S2Point> vertices)
    {
        var base_ = loop.Vertices;
        vertices.AddRange(base_);
    }

    // Returns a vector of points shaped as a regular polygon with
    // num_vertices vertices, all on a circle of the specified angular
    // radius around the center.  The radius is the actual distance from
    // the center to the circle along the sphere.
    //
    // If you want to construct a regular polygon, try this:
    //   S2Polygon polygon(S2Loop.MakeRegularLoop(center, radius, num_vertices));
    public static S2Point[] MakeRegularPoints(S2Point center, S1Angle radius, int num_vertices)
    {
        var loop = S2Loop.MakeRegularLoop(center, radius, num_vertices);
        var points = new S2Point[loop.NumVertices];
        for (int i = 0; i < loop.NumVertices; i++)
        {
            points[i] = loop.Vertex(i);
        }
        return points;
    }

    // Convert a distance on the Earth's surface to an angle.
    // Do not use these methods in non-testing code; use s2earth.h instead.
    public static S1Angle MetersToAngle(double meters)
    {
        return KmToAngle(0.001 * meters);
    }

    public static S1Angle KmToAngle(double km)
    {
        return S1Angle.FromRadians(km / S2Earth.RadiusKm);
    }

    // Convert an area in steradians (as returned by the S2 area methods) to
    // square meters or square kilometers.
    public static double AreaToMeters2(double steradians)
    {
        return 1e6 * AreaToKm2(steradians);
    }

    public static double AreaToKm2(double steradians)
    {
        return steradians * S2Earth.RadiusKm * S2Earth.RadiusKm;
    }

    // The Dump*() functions are for use within a debugger.  They are similar to
    // the corresponding s2textformat::ToString() functions except that they
    // prefix their output with a label and they don't require default arguments
    // or constructing absl::Span objects (which gdb doesn't know how to do).
    public static string Dump(S2Point p)
    {
        return "S2Point: " + p.ToDebugString() + Environment.NewLine;
    }

    public static string Dump(S2Point[] points)
    {
        return "S2Polygon: " + points.ToDebugString() + Environment.NewLine;
    }
    
    public static string Dump(S2Loop loop)
    {
        return "S2Polygon: " + loop.ToDebugString() + Environment.NewLine;
    }

    public static string Dump(S2Polyline polyline)
    {
        return "S2Polyline: " + polyline.ToDebugString() + Environment.NewLine;
    }

    public static string Dump(S2Polygon polygon)
    {
        return "S2Polygon: " + polygon.ToDebugString() + Environment.NewLine;
    }

    public static string Dump(S2LaxPolylineShape polyline)
    {
        return "S2Polyline: " + polyline.ToDebugString() + Environment.NewLine;
    }
    
    public static string Dump(S2LaxPolygonShape polygon)
    {
        return "S2Polygon: " + polygon.ToDebugString() + Environment.NewLine;
    }

    // Outputs the contents of an S2ShapeIndex in human-readable form.
    public static string Dump(S2ShapeIndex index)
    {
        var sb = new StringBuilder();
        sb.AppendLine("S2ShapeIndex: " + index);
        sb.AppendLine("  " + index.ToDebugString());
        foreach (var it in index.GetNewEnumerable())
        {
            sb.AppendLine("  id: " + it.Item1);
            var cell = it.Item2;
            for (var s = 0; s < cell.NumClipped(); ++s)
            {
                var clipped = cell.Clipped(s);
                sb.Append("    shape_id " + clipped.ShapeId + ": ");
                for (var e = 0; e < clipped.NumEdges; ++e)
                {
                    if (e > 0) sb.Append(", ");
                    sb.Append(clipped.Edge(e));
                }
                sb.AppendLine();
            }
        }
        return sb.ToString();
    }

    // Return a random unit-length vector.
    public static S2Point RandomPoint()
    {
        // The order of evaluation of function arguments is unspecified,
        // so we may not just call S2Point with three RandDouble-based args.
        // Use temporaries to induce sequence points between calls.
        double x = Random.UniformDouble(-1, 1);
        double y = Random.UniformDouble(-1, 1);
        double z = Random.UniformDouble(-1, 1);
        return new S2Point(x, y, z).Normalize();
    }

    // Return a right-handed coordinate frame (three orthonormal vectors).
    public static void GetRandomFrame(out S2Point x, out S2Point y, out S2Point z)
    {
        z = RandomPoint();
        GetRandomFrameAt(z, out x, out y);
    }

    public static S2PointS2Point GetRandomFrame()
    {
        return GetRandomFrameAt(RandomPoint());
    }

    // Given a unit-length z-axis, compute x- and y-axes such that (x,y,z) is a
    // right-handed coordinate frame (three orthonormal vectors).
    public static void GetRandomFrameAt(S2Point z, out S2Point x, out S2Point y)
    {
        x = z.CrossProd(RandomPoint()).Normalize();
        y = z.CrossProd(x).Normalize();
    }

    public static S2PointS2Point GetRandomFrameAt(S2Point z)
    {
        GetRandomFrameAt(z, out var x, out var y);
        return new S2PointS2Point(x, y, z);
    }

    // Return a random cell id at the given level or at a randomly chosen
    // level.  The distribution is uniform over the space of cell ids,
    // but only approximately uniform over the surface of the sphere.
    public static S2CellId GetRandomCellId(int level)
    {
        var face = Random.Uniform(S2CellId.kNumFaces);
        var pos = Random.Rand64() & ((1UL + S2CellId.kPosBits) - 1);
        return S2CellId.FromFacePosLevel(face, pos, level);
    }

    public static S2CellId GetRandomCellId()
    {
        return GetRandomCellId(Random.Uniform(S2.kMaxCellLevel + 1));
    }

    // Return a cap with a random axis such that the log of its area is
    // uniformly distributed between the logs of the two given values.
    // (The log of the cap angle is also approximately uniformly distributed.)
    public static S2Cap GetRandomCap(double min_area, double max_area)
    {
        double cap_area = max_area * Math.Pow(min_area / max_area, Random.RandDouble());
        Assert.True(cap_area >= min_area);
        Assert.True(cap_area <= max_area);

        // The surface area of a cap is 2*Pi times its height.
        return S2Cap.FromCenterArea(RandomPoint(), cap_area);
    }

    // Return a polygon with the specified center, number of concentric loops
    // and vertices per loop.
    public static void ConcentricLoopsPolygon(S2Point center, int num_loops, int num_vertices_per_loop, out S2Polygon polygon)
    {
        var m = S2.GetFrame(center);
        var loops = new List<S2Loop>();
        for (int li = 0; li < num_loops; ++li)
        {
            var vertices = new List<S2Point>();
            double radius = 0.005 * (li + 1) / num_loops;
            double radian_step = S2.M_2_PI / num_vertices_per_loop;
            for (int vi = 0; vi < num_vertices_per_loop; ++vi)
            {
                double angle = vi * radian_step;
                var p = new S2Point(radius * Math.Cos(angle), radius * Math.Sin(angle), 1);
                vertices.Add(S2.FromFrame(m, p.Normalize()));
            }
            loops.Add(new S2Loop(vertices));
        }
        polygon = new S2Polygon(loops);
    }

    // Return a point chosen uniformly at random (with respect to area)
    // from the given cap.
    public static S2Point SamplePoint(S2Cap cap)
    {
        // We consider the cap axis to be the "z" axis.  We choose two other axes to
        // complete the coordinate frame.

        var m = S2.GetFrame(cap.Center);

        // The surface area of a spherical cap is directly proportional to its
        // height.  First we choose a random height, and then we choose a random
        // point along the circle at that height.

        double h = Random.RandDouble() * cap.Height();
        double theta = S2.M_2_PI * Random.RandDouble();
        double r = Math.Sqrt(h * (2 - h));  // Radius of circle.

        // The result should already be very close to unit-length, but we might as
        // well make it accurate as possible.
        return S2.FromFrame(m, new S2Point(Math.Cos(theta) * r, Math.Sin(theta) * r, 1 - h))
               .Normalize();
    }

    // Return a point chosen uniformly at random (with respect to area on the
    // sphere) from the given latitude-longitude rectangle.
    public static S2Point SamplePoint(S2LatLngRect rect)
    {
        // First choose a latitude uniformly with respect to area on the sphere.
        double sin_lo = Math.Sin(rect.Lat.Lo);
        double sin_hi = Math.Sin(rect.Lat.Hi);
        double lat = Math.Asin(Random.UniformDouble(sin_lo, sin_hi));

        // Now choose longitude uniformly within the given range.
        double lng = rect.Lng.Lo + Random.RandDouble() * rect.Lng.GetLength();
        return S2LatLng.FromRadians(lat, lng).Normalized().ToPoint();
    }

    // Checks that "covering" completely covers the given region.  If
    // "check_tight" is true, also checks that it does not contain any cells
    // that do not intersect the given region.  ("id" is only used internally.)
    public static void CheckCovering(IS2Region region, S2CellUnion covering, bool check_tight, S2CellId? id = null)
    {
        id ??= new S2CellId();
        if (!id.Value.IsValid())
        {
            for (int face = 0; face < 6; ++face)
            {
                CheckCovering(region, covering, check_tight, S2CellId.FromFace(face));
            }
            return;
        }

        if (!region.MayIntersect(new S2Cell(id.Value)))
        {
            // If region does not intersect id, then neither should the covering.
            if (check_tight) Assert.True(!covering.Intersects(id.Value));

        }
        else if (!covering.Contains(id.Value))
        {
            // The region may intersect id, but we can't assert that the covering
            // intersects id because we may discover that the region does not actually
            // intersect upon further subdivision.  (MayIntersect is not exact.)
            Assert.True(!region.Contains(new S2Cell(id.Value)));
            Assert.True(!id.Value.IsLeaf());
            S2CellId end = id.Value.ChildEnd();
            S2CellId child;
            for (child = id.Value.ChildBegin(); child != end; child = child.Next())
            {
                CheckCovering(region, covering, check_tight, child);
            }
        }
    }

    // A simple class that generates "Koch snowflake" fractals (see Wikipedia
    // for an introduction).  There is an option to control the fractal
    // dimension (between 1.0 and 2.0); values between 1.02 and 1.50 are
    // reasonable simulations of various coastlines.  The default dimension
    // (about 1.26) corresponds to the standard Koch snowflake.  (The west coast
    // of Britain has a fractal dimension of approximately 1.25.)
    //
    // The fractal is obtained by starting with an equilateral triangle and
    // recursively subdividing each edge into four segments of equal length.
    // Therefore the shape at level "n" consists of 3*(4**n) edges.  Multi-level
    // fractals are also supported: if you set min_level() to a non-negative
    // value, then the recursive subdivision has an equal probability of
    // stopping at any of the levels between the given min and max (inclusive).
    // This yields a fractal where the perimeter of the original triangle is
    // approximately equally divided between fractals at the various possible
    // levels.  If there are k distinct levels {min,..,max}, the expected number
    // of edges at each level "i" is approximately 3*(4**i)/k.
    public class Fractal
    {
        // You must call set_max_level() or SetLevelForApproxMaxEdges() before
        // calling MakeLoop().
        public Fractal()
        {
            _maxLevel = -1;
            min_level_arg_ = -1;
            min_level_ = -1;
            FractalDimension = Math.Log(4) / Math.Log(3); /* standard Koch curve */
            edge_fraction_ = 0;
            offset_fraction_ = 0;
            ComputeOffsets();
        }

        // Set the maximum subdivision level for the fractal (see above).
        // REQUIRES: max_level >= 0
        public int MaxLevel
        {
            get => _maxLevel;
            set
            {
                Assert.True(value >= 0);
                _maxLevel = value;
                ComputeMinLevel();
            }
        }
        private int _maxLevel;

        // Set the minimum subdivision level for the fractal (see above).  The
        // default value of -1 causes the min and max levels to be the same.  A
        // min_level of 0 should be avoided since this creates a significant
        // chance that none of the three original edges will be subdivided at all.
        //
        // DEFAULT: max_level()
        public int MinLevel
        {
            get => min_level_arg_;
            set
            {
                Assert.True(value >= -1);
                min_level_arg_ = value;
                ComputeMinLevel();
            }
        }
        private int min_level_;      // Actual min level (depends on max_level_)
        private int min_level_arg_;  // Value set by user

        private void ComputeMinLevel()
        {
            if (min_level_arg_ >= 0 && min_level_arg_ <= MaxLevel)
            {
                min_level_ = min_level_arg_;
            }
            else
            {
                min_level_ = MaxLevel;
            }
        }

        // Set the fractal dimension.  The default value of approximately 1.26
        // corresponds to the stardard Koch curve.  The value must lie in the
        // range [1.0, 2.0).
        //
        // DEFAULT: log(4) / log(3) ~= 1.26
        public double FractalDimension
        {
            get => _fractalDimension;
            set
            {
                Assert.True(value >= 1.0);
                Assert.True(value < 2);
                _fractalDimension = value;
                ComputeOffsets();
            }
        }
        private double _fractalDimension;

        private void ComputeOffsets()
        {
            edge_fraction_ = Math.Pow(4.0, -1.0 / FractalDimension);
            offset_fraction_ = Math.Sqrt(edge_fraction_ - 0.25);
        }

        // Set the min and/or max level to produce approximately the given number
        // of edges.  (The values are rounded to a nearby value of 3*(4**n).)
        public void SetLevelForApproxMinEdges(int min_edges)
        {
            // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
            MinLevel = ((int)Math.Round(0.5 * Math.Log2(min_edges / 3)));
        }

        public void SetLevelForApproxMaxEdges(int max_edges)
        {
            // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
            MaxLevel = ((int)Math.Round(0.5 * Math.Log2(max_edges / 3)));
        }

        // Return a lower bound on ratio (Rmin / R), where "R" is the radius
        // passed to MakeLoop() and "Rmin" is the minimum distance from the
        // fractal boundary to its center, where all distances are measured in the
        // tangent plane at the fractal's center.  This can be used to inscribe
        // another geometric figure within the fractal without intersection.
        public double MinRadiusRactor()
        {
            // The minimum radius is attained at one of the vertices created by the
            // first subdivision step as long as the dimension is not too small (at
            // least kMinDimensionForMinRadiusAtLevel1, see below).  Otherwise we fall
            // back on the incircle radius of the original triangle, which is always a
            // lower bound (and is attained when dimension = 1).
            //
            // The value below was obtained by letting AE be an original triangle edge,
            // letting ABCDE be the corresponding polyline after one subdivision step,
            // and then letting BC be tangent to the inscribed circle at the center of
            // the fractal O.  This gives rise to a pair of similar triangles whose edge
            // length ratios can be used to solve for the corresponding "edge fraction".
            // This method is slightly conservative because it is computed using planar
            // rather than spherical geometry.  The value below is equal to
            // -log(4)/log((2 + cbrt(2) - cbrt(4))/6).
            double kMinDimensionForMinRadiusAtLevel1 = 1.0852230903040407;
            if (FractalDimension >= kMinDimensionForMinRadiusAtLevel1)
            {
                return Math.Sqrt(1 + 3 * edge_fraction_ * (edge_fraction_ - 1));
            }
            return 0.5;
        }

        // Return the ratio (Rmax / R), where "R" is the radius passed to
        // MakeLoop() and "Rmax" is the maximum distance from the fractal boundary
        // to its center, where all distances are measured in the tangent plane at
        // the fractal's center.  This can be used to inscribe the fractal within
        // some other geometric figure without intersection.
        public double MaxRadiusFactor()
        {
            // The maximum radius is always attained at either an original triangle
            // vertex or at a middle vertex from the first subdivision step.
            return Math.Max(1.0, offset_fraction_ * Math.Sqrt(3) + 0.5);
        }

        private void GetR2Vertices(List<R2Point> vertices)
        {
            // The Koch "snowflake" consists of three Koch curves whose initial edges
            // form an equilateral triangle.
            var v0 = new R2Point(1.0, 0.0);
            var v1 = new R2Point(-0.5, Math.Sqrt(3) / 2);
            var v2 = new R2Point(-0.5, -Math.Sqrt(3) / 2);
            GetR2VerticesHelper(v0, v1, 0, vertices);
            GetR2VerticesHelper(v1, v2, 0, vertices);
            GetR2VerticesHelper(v2, v0, 0, vertices);
        }

        // Given the two endpoints (v0,v4) of an edge, recursively subdivide the edge
        // to the desired level, and insert all vertices of the resulting curve up to
        // but not including the endpoint "v4".
        private void GetR2VerticesHelper(R2Point v0, R2Point v4, int level, List<R2Point> vertices)
        {
            if (level >= min_level_ && Random.OneIn(MaxLevel - level + 1))
            {
                // Stop subdivision at this level.
                vertices.Add(v0);
                return;
            }
            // Otherwise compute the intermediate vertices v1, v2, and v3.
            var dir = v4 - v0;
            var v1 = v0 + edge_fraction_ * dir;
            var v2 = 0.5 * (v0 + v4) - offset_fraction_ * dir.GetOrtho();
            var v3 = v4 - edge_fraction_ * dir;

            // And recurse on the four sub-edges.
            GetR2VerticesHelper(v0, v1, level + 1, vertices);
            GetR2VerticesHelper(v1, v2, level + 1, vertices);
            GetR2VerticesHelper(v2, v3, level + 1, vertices);
            GetR2VerticesHelper(v3, v4, level + 1, vertices);
        }

        // Return a fractal loop centered around the z-axis of the given
        // coordinate frame, with the first vertex in the direction of the
        // positive x-axis.  In order to avoid self-intersections, the fractal is
        // generated by first drawing it in a 2D tangent plane to the unit sphere
        // (touching at the fractal's center point) and then projecting the edges
        // onto the sphere.  This has the side effect of shrinking the fractal
        // slightly compared to its nominal radius.
        public S2Loop MakeLoop(S2PointS2Point frame, S1Angle nominal_radius)
        {
            List<R2Point> r2vertices = new();
            GetR2Vertices(r2vertices);
            List<S2Point> vertices = new();
            var r = nominal_radius.Radians;
            foreach (var v in r2vertices)
            {
                S2Point p = new(v[0] * r, v[1] * r, 1);
                vertices.Add(S2.FromFrame(frame, p).Normalize());
            }
            return new S2Loop(vertices);
        }

        // The ratio of the sub-edge length to the original edge length at each
        // subdivision step.
        private double edge_fraction_;

        // The distance from the original edge to the middle vertex at each
        // subdivision step, as a fraction of the original edge length.
        private double offset_fraction_;
    }
}
