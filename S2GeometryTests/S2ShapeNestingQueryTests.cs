namespace S2Geometry;

public class S2ShapeNestingQueryTests
{
    #region Data structures and utility methods

    internal record struct RingSpec(
        S2LatLng Center,
        double RadiusDeg,
        // Should we reverse the vertex order?
        bool Reverse = false);

    // Builds an S2LaxPolygonShape out of nested rings according to given ring
    // specs.  Each spec specifies a ring center and a radius in degrees.  The
    // radius will be made positive if it's not already, and crossing the poles will
    // cause a S2_CHECK failure.
    //
    // Rings are generated in counter-clockwise orientation around their center
    // by default, if reverse is specified, then the orientation becomes clockwise.
    internal static S2LaxPolygonShape RingShape(int verticesPerLoop, List<RingSpec> ringSpecs)
    {
        var radian_step = 2 * S2.M_PI / verticesPerLoop;

        List<List<S2Point>> loops = [];
        foreach (var spec in ringSpecs)
        {
            // Check that we're not in reach of the poles.
            var radius = Math.Abs(spec.RadiusDeg);
            Assert.True(spec.Center.Lat().GetDegrees() + radius < +90);
            Assert.True(spec.Center.Lat().GetDegrees() - radius > -90);

            List<S2Point> vertices = new(verticesPerLoop);
            for (int i = 0; i < verticesPerLoop; ++i)
            {
                double angle = i * radian_step;
                S2LatLng pnt = S2LatLng.FromDegrees(radius * Math.Sin(angle), radius * Math.Cos(angle));
                vertices.Add((spec.Center + pnt).Normalized().ToPoint());
            }

            if (spec.Reverse)
            {
                vertices.Reverse();
            }

            loops.Add(vertices);
        }

        return new S2LaxPolygonShape(loops);
    }

    // Specify a circular arc about a center point.  The arc has thickness and
    // extends from the starting angle to the ending angle.
    //
    // If offset is specified, the result array of points is rotated by that amount
    // to change which vertex is first.
    //
    // By default, the resulting shape is oriented CCW around its center point, but
    // if reverse is set, then the points are reversed and its orientation changes
    // to CW.
    internal record struct ArcSpec(
        S2LatLng Center,
        double RadiusDeg,
        double Thickness,
        double StartDeg,
        double EndDeg,
        // If non-zero, rotate ring vertices by this many points.  The total "real"
        // shift will be offset % vertices_per_loop.
        int Offset,
        // Should we reverse order of vertices?
        bool Reverse);

    // Builds an S2LaxPolygonShape from one or more `ArcSpec`s.  Each spec yields an
    // arc on a circle made to have the specified thickness. The inner and outer
    // edges have their ends connected with a butt cap.  The arcs are generated in
    // counter-clockwise order around their center by default.  If reverse is given,
    // then the orientation becomes clockwise.
    internal static S2LaxPolygonShape ArcShape(int vertices_per_loop, List<ArcSpec> specs)
    {
        var deg2rad = (double degrees) => { return S2.M_PI / 180.0 * degrees; };

        List<List<S2Point>> loops = [];
        foreach (var spec in specs)
        {
            double start_rad = deg2rad(spec.StartDeg);
            double end_rad = deg2rad(spec.EndDeg);

            Assert.True(start_rad < end_rad);
            Assert.True(spec.RadiusDeg > 0);
            Assert.True(spec.Thickness > 0);
            Assert.True(vertices_per_loop % 2 == 0);

            double radius_inner = spec.RadiusDeg - spec.Thickness;
            double radius_outer = spec.RadiusDeg + spec.Thickness;
            double radian_step =
                (end_rad - start_rad) / (vertices_per_loop / 2 - 1);

            // Don't allow arcs to go over the poles to avoid the singularity there
            // which makes the resulting logic much more complex.
            Assert.True(spec.Center.Lat().GetDegrees() + spec.RadiusDeg + spec.Thickness <
                        +90);
            Assert.True(spec.Center.Lat().GetDegrees() - spec.RadiusDeg - spec.Thickness >
                        -90);

            // Generate inner and outer edges at the same time with implied butt joint.
            List<S2Point> vertices = new(vertices_per_loop);
            for (int i = 0; i < vertices_per_loop / 2; ++i)
            {
                double angle = start_rad + i * radian_step;
                double sina = Math.Sin(angle);
                double cosa = Math.Cos(angle);

                S2LatLng pnt0 =
                    S2LatLng.FromDegrees(radius_outer * sina, radius_outer * cosa);
                S2LatLng pnt1 =
                    S2LatLng.FromDegrees(radius_inner * sina, radius_inner * cosa);

                vertices[i] = (spec.Center + pnt0).Normalized().ToPoint();
                vertices[vertices_per_loop - i - 1] =
                    (spec.Center + pnt1).Normalized().ToPoint();
            }

            if (spec.Offset != 0)
            {
                int offset = spec.Offset % vertices_per_loop;
                vertices.RotateInPlace(offset);
            }

            if (spec.Reverse)
            {
                vertices.Reverse();
            }

            loops.Add(vertices);
        }
        return new S2LaxPolygonShape(loops);
    }

    #endregion

    [Fact]
    internal void Test_S2ShapeNestingQuery_OneChainAlwaysShell()
    {
        const int kNumEdges = 100;

        MutableS2ShapeIndex index = [];
        int id = index.Add(RingShape(kNumEdges, [new(S2LatLng.FromDegrees(0.0, 0.0), 1.0)]));

        S2ShapeNestingQuery query = new(index);
        var relations = query.ComputeShapeNesting(id);
        Assert.Equal(relations.Count, 1);
        Assert.True(relations[0].IsShell());
        Assert.False(relations[0].IsHole());
        Assert.True(relations[0].Parent < 0);
        Assert.Equal(relations[0].GetHoles().Count, 0);
    }

    [Fact]
    internal void Test_S2ShapeNestingQuery_TwoChainsFormPair()
    {
        const int kNumEdges = 100;
        var kCenter = S2LatLng.FromDegrees(0.0, 0.0);

        {
            // Nested rings, like a donut.
            MutableS2ShapeIndex index = [];
            int id = index.Add(
                RingShape(kNumEdges, [new(kCenter, 1.0, false), new(kCenter, 0.5, true)]));

            S2ShapeNestingQuery query = new(index);
            var relations = query.ComputeShapeNesting(id);

            // First chain should be a shell and the second a hole.
            Assert.Equal(relations.Count, 2);
            Assert.True(relations[0].IsShell());
            Assert.True(relations[1].IsHole());
            Assert.False(relations[0].IsHole());
            Assert.False(relations[1].IsShell());

            // First chain should have no parent and one hole, which is second chain.
            Assert.True(relations[0].Parent < 0);
            Assert.Equal(relations[0].GetHoles().Count, 1);
            Assert.Equal(relations[0].GetHoles()[0], 1);

            // Second chain should have one parent, which is chain zero, and no holes.
            Assert.Equal(relations[1].Parent, 0);
            Assert.Equal(relations[1].GetHoles().Count, 0);
        }

        {
            // Swapping ring ordering shouldn't change anything.
            MutableS2ShapeIndex index = [];
            int id = index.Add(
                RingShape(kNumEdges, [new(kCenter, 0.5, true), new(kCenter, 1.0, false)]));

            S2ShapeNestingQuery query = new(index);
            var relations = query.ComputeShapeNesting(id);

            // First chain should be a shell and the second a hole.
            Assert.Equal(relations.Count, 2);
            Assert.True(relations[0].IsShell());
            Assert.True(relations[1].IsHole());
            Assert.False(relations[0].IsHole());
            Assert.False(relations[1].IsShell());

            // First chain should have no parent and one hole, which is second chain.
            Assert.True(relations[0].Parent < 0);
            Assert.Equal(relations[0].GetHoles().Count, 1);
            Assert.Equal(relations[0].GetHoles()[0], 1);

            // Second chain should have one parent, which is chain zero, and no holes.
            Assert.Equal(relations[1].Parent, 0);
            Assert.Equal(relations[1].GetHoles().Count, 0);
        }

        {
            // If we reverse the vertex order of the rings.  We should end up with two
            // shells since the hole and shell don't face each other.
            MutableS2ShapeIndex index = [];
            int id = index.Add(
                RingShape(kNumEdges, [new(kCenter, 1.0, true), new(kCenter, 0.5, false)]));

            S2ShapeNestingQuery query = new(index);
            var relations = query.ComputeShapeNesting(id);

            // Both chains should be a shell with no holes
            Assert.Equal(relations.Count, 2);
            for (int i = 0; i < 2; ++i)
            {
                Assert.True(relations[i].IsShell());
                Assert.False(relations[i].IsHole());
                Assert.True(relations[i].Parent < 0);
                Assert.Equal(relations[i].GetHoles().Count, 0);
            }
        }
    }

    [Fact]
    internal void Test_S2ShapeNestingQuery_CanSetDatumShellOption()
    {
        const int kNumEdges = 100;
        var kCenter = S2LatLng.FromDegrees(0.0, 0.0);

        // Nested rings, like a donut.
        MutableS2ShapeIndex index = [];
        int id = index.Add(
            RingShape(kNumEdges, [new(kCenter, 1.0, false), new(kCenter, 0.5, true)]));

        // We should be able to override the default datum shell strategy.
        S2ShapeNestingQuery.Options options = new()
        {
            DatumStrategy = (S2Shape s) => 1
        };
        S2ShapeNestingQuery query = new(index, options);

        var relations = query.ComputeShapeNesting(id);

        // Second chain should be a shell and the first a hole.
        Assert.Equal(relations.Count, 2);
        Assert.True(relations[1].IsShell());
        Assert.True(relations[0].IsHole());
        Assert.False(relations[1].IsHole());
        Assert.False(relations[0].IsShell());
    }

    [Fact]
    internal void Test_S2ShapeNestingQuery_ShellCanHaveMultipleHoles()
    {
        const int kNumEdges = 16;

        // A ring with four holes in it like a shirt button.
        MutableS2ShapeIndex index = [];
        int id = index.Add(
            RingShape(kNumEdges, [
            new(S2LatLng.FromDegrees(0.5, 0.5), 2.0),
            new(S2LatLng.FromDegrees(1.0, 0.5), 0.25, true),
            new(S2LatLng.FromDegrees(0.0, 0.5), 0.25, true),
            new(S2LatLng.FromDegrees(0.5, 1.0), 0.25, true),
            new(S2LatLng.FromDegrees(0.5, 0.0), 0.25, true)
            ]));

        S2ShapeNestingQuery query = new(index);
        var relations = query.ComputeShapeNesting(id);

        // First chain should be a shell and have four holes.
        Assert.Equal(relations.Count, 5);
        Assert.True(relations[0].IsShell());
        Assert.False(relations[0].IsHole());
        Assert.True(relations[0].Parent < 0);
        Assert.Equal(relations[0].GetHoles().Count, 4);

        // Chain zero should have the others as holes, and the others should have
        // chain zero as their parent.
        for (int i = 1; i < 5; ++i)
        {
            Assert.Equal(relations[0].GetHoles()[i - 1], i);

            Assert.True(relations[i].IsHole());
            Assert.False(relations[i].IsShell());
            Assert.Equal(relations[i].Parent, 0);
            Assert.Equal(relations[i].GetHoles().Count, 0);
        }
    }

    [Fact]
    internal void Test_S2ShapeNestingQuery_ExactPathIsIrrelevant()
    {
        const int kNumEdges = 32;
        var kCenter = S2LatLng.FromDegrees(0.0, 0.0);

        // The path we take from the datum shell to the inner shell shouldn't matter
        // for the final classification.  So build a set of nested rings that are open
        // on one end (like a 'C'), so that they're highly concave.  Shift the datum
        // ring and other rings a point at a time to get coverage on the various
        // permutations.
        for (int offset0 = 0; offset0 < kNumEdges; ++offset0)
        {
            for (int offset1 = 0; offset1 < kNumEdges; ++offset1)
            {
                //S2_VLOG(1) << "Offset (" << offset0 << "," << offset1 << ")";

                MutableS2ShapeIndex index = [];
                int id = index.Add(ArcShape(
                    //            center radius thickness start end  offset  reverse
                    kNumEdges, [
                new(kCenter, 0.3, 0.15, -240.0, 60.0, offset0, false),
                new(kCenter, 0.3, 0.05, -230.0, 50.0, offset1, true),
                new(kCenter, 1.0, 0.15, -85.0, 265.0, offset1, false),
                new(kCenter, 1.0, 0.05, -80.0, 260.0, offset1, true)
                ]));

                S2ShapeNestingQuery query = new(index);
                var relations = query.ComputeShapeNesting(id);

                Assert.Equal(relations.Count, 4);
                Assert.True(relations[0].IsShell());
                Assert.True(relations[1].IsHole());
                Assert.Equal(relations[1].Parent, 0);
                Assert.True(relations[2].IsShell());
                Assert.True(relations[3].IsHole());
                Assert.Equal(relations[3].Parent, 2);
            }
        }
    }

    [Theory]
    // Test even/odd number of rings, outer ring is first.
    [InlineData(31, 0, false)]
    [InlineData(32, 0, false)]
    [InlineData(31, 0, true)]
    [InlineData(32, 0, true)]

            // Test even/odd number of rings, last ring is first.
    [InlineData(31, 30, true)]
    [InlineData(32, 31, true)]

            // Tof rings, intermediate ring is first.
    [InlineData(31, 31 / 13, true)]
    [InlineData(32, 32 / 13, true)]
    [InlineData(31, 31 / 3, true)]
    [InlineData(32, 32 / 3, true)]
    internal void Test_NestingTest_NestedChainsPartitionCorrectly(
        int depth,        // How many nested loops to generate
        int firstChain,   // Which nested loop is the first loop in the list
        bool shuffle)    // If true, shuffle the loop order (other than first loop)
    {
        const int kNumEdges = 16;
        var kCenter = S2LatLng.FromDegrees(0.0, 0.0);

        List<RingSpec> rings =
        [
            new RingSpec(kCenter, 2.0 / (firstChain + 1), firstChain % 2 == 1)
        ];

        for (int i = 0; i < depth; ++i)
        {
            if (i == firstChain)
            {
                continue;
            }
            rings.Add(new RingSpec(kCenter, 2.0 / (i + 1), i % 2 == 1));
        }

        if (shuffle)
        {
            S2Testing.Random.Reset(S2Testing.Random.RandomSeed);
            //gtl.legacy_random_shuffle(rings.begin() + 1, rings.end(), S2Testing.rnd);
        }

        MutableS2ShapeIndex index = [];
        int id = index.Add(RingShape(kNumEdges, rings));

        S2ShapeNestingQuery query = new(index);
        var relations = query.ComputeShapeNesting(id);
        Assert.Equal(relations.Count, depth);

        // In the case of the outer ring being the first chain, and no shuffling,
        // then the outer chain should be a shell and then we alternate hole shell
        // as we move inwards.
        if (firstChain == 0 && !shuffle)
        {
            Assert.True(relations[0].IsShell());
            Assert.Equal(relations[0].GetHoles().Count, 1);
            Assert.Equal(relations[0].GetHoles()[0], 1);

            for (int chain = 1; chain < depth; chain++)
            {
                if (chain % 2 == 1)
                {
                    // We expect this chain to be a hole
                    Assert.True(relations[chain].IsHole());
                    Assert.False(relations[chain].IsShell());
                    Assert.Equal(relations[chain].Parent, chain - 1);
                }
                else
                {
                    // We expect this chain to be a shell
                    Assert.False(relations[chain].IsHole());
                    Assert.True(relations[chain].IsShell());
                    Assert.Equal(relations[chain].Parent, -1);
                }
            }
        }

        // We should always be able to divide the set of chains into pairs of shells
        // and holes, possibly with one shell left over.
        int num_shells = 0;
        int num_holes = 0;
        for (int chain = 0; chain < depth; chain++)
        {
            if (relations[chain].IsShell())
            {
                num_shells++;
                foreach (var child in relations[chain].GetHoles())
                {
                    Assert.Equal(relations[child].Parent, chain);
                }
            }

            if (relations[chain].IsHole())
            {
                num_holes++;
                int parent = relations[chain].Parent;
                var holes = relations[parent].GetHoles();
                //Assert.NE(std.find(holes.begin(), holes.end(), chain), holes.end());
            }
        }

        // Everything is a hole or a shell.
        Assert.Equal(num_holes + num_shells, depth);
    }
}
