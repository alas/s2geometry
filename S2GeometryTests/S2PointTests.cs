namespace S2Geometry;

public class S2PointTests
{
    [Fact]
    internal void Test_S2Point_HashSpreads()
    {
        int kTestPoints = 1 << 16;
        List<int> set = new();
        SortedSet<S2Point> points = new();
        S2Point base_ = new(1, 1, 1);
        for (var i = 0; i < kTestPoints; ++i)
        {
            // All points in a tiny cap to test avalanche property of hash
            // function (the cap would be of radius 1mm on Earth (4*10^9/2^35).
            S2Point perturbed = base_ + S2Testing.RandomPoint() / (1UL << 35);
            perturbed = perturbed.Normalize();
            set.Add(perturbed.GetHashCode());
            points.Add(perturbed);
        }
        // A real collision is extremely unlikely.
        Assert.Equal(0, kTestPoints - points.Count);
        // Allow a few for the hash.
        Assert.True(10 >= (kTestPoints - set.Count));
    }

    [Fact]
    internal void Test_S2Point_IsAVector()
    {
        // Check that we haven't accidentally increased the size or modified
        // the alignment of S2Point.
        Assert.Equal(typeof(S2Point).SizeOf(), typeof(Vector3<double>).SizeOf());
        Assert.Equal(typeof(S2Point).AlignOf(), typeof(Vector3<double>).AlignOf());
    }

    /*[Fact]
    internal void Test_S2Point_CoderWorks()
    {
        S2Testing.Random.Reset(S2Testing.Random.RandomSeed);

        S2Point point = S2Testing.RandomPoint();
        S2Error error;
        S2Point decoded = S2Coder_Testing.RoundTrip(new S2Point_Coder(), point, out error);
        Assert.True(error.IsOk());
        Assert.Equal(decoded, point);
    }*/

    [Fact]
    internal void Test_S2Point_SubtractionWorks()
    {
        S2Point a = new(1, 2, 3);
        S2Point b = new(1, 1, 1);
        a -= b;
        Assert.Equal(a, new S2Point(0, 1, 2));
    }

    [Fact]
    internal void Test_S2Point_ElementWiseDivisionWorks()
    {
        S2Point a = new(4, 8, 16);
        S2Point b = new(2, 2, 2);
        Assert.Equal(a.DivComponents(b), new S2Point(2, 4, 8));
    }

    [Fact]
    internal void Test_S2Point_SqrtWorks()
    {
        S2Point a = new(4, 9, 16);
        Assert.Equal(a.Sqrt(), new S2Point(2, 3, 4));
    }

    [Fact]
    internal void Test_S2Point_FloorWorks()
    {
        S2Point a = new(1.4, 1.5, 1.6);
        Assert.Equal(a.Floor(), new S2Point(1, 1, 1));
    }

    [Fact]
    internal void Test_S2Point_CeilWorks()
    {
        S2Point a = new(1.4, 1.5, 1.6);
        Assert.Equal(a.Ceiling(), new S2Point(2, 2, 2));
    }

    [Fact]
    internal void Test_S2Point_FRoundWorks()
    {
        S2Point a = new(1.4, 1.5, 1.6);
        Assert.Equal(a.Round(), new S2Point(1, 2, 2));
    }
}
