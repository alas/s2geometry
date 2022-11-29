namespace S2Geometry;

public class S2PointTests
{
    [Fact]
    public void Test_S2_GoodFastHashSpreads()
    {
        int kTestPoints = 1 << 16;
        var set = new List<int>();
        var points = new SortedSet<S2Point>();
        var base_ = new S2Point(1, 1, 1);
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
        Assert.True(10 >= kTestPoints - set.Count);
    }
}
