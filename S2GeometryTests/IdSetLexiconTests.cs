namespace S2Geometry;

using IdSet = List<Int32>;
using Seq = List<Int32>;

public class IdSetLexiconTests
{
    [Fact]
    internal void Test_IdSetLexicon_EmptySet()
    {
        IdSetLexicon lexicon = new();
        ExpectIdSet([], lexicon.IdSet_(lexicon.Add([])));
    }

    [Fact]
    internal void Test_IdSetLexicon_SingletonSets()
    {
        IdSetLexicon lexicon = new();
        Assert.Equal(5, lexicon.Add([5]));
        Assert.Equal(0, lexicon.Add([0, 0]));
        Assert.Equal(1, IdSetLexicon.AddSingleton(1));
        var m = int.MaxValue;
        Assert.Equal(m, lexicon.Add([m]));

        ExpectIdSet([0], lexicon.IdSet_(0));
        ExpectIdSet([1], lexicon.IdSet_(1));
        ExpectIdSet([5], lexicon.IdSet_(5));
        ExpectIdSet([m], lexicon.IdSet_(m));
    }

    [Fact]
    internal void Test_IdSetLexicon_SetsAreSorted()
    {
        IdSetLexicon lexicon = new();
        Assert.Equal(~0, lexicon.Add([2, 5]));
        Assert.Equal(~1, lexicon.Add([3, 2, 5]));
        Assert.Equal(~0, lexicon.Add([5, 2]));
        Assert.Equal(~1, lexicon.Add([5, 3, 2, 5]));

        ExpectIdSet([2, 5], lexicon.IdSet_(~0));
        ExpectIdSet([2, 3, 5], lexicon.IdSet_(~1));
    }

    [Fact]
    internal void Test_IdSetLexicon_Clear()
    {
        IdSetLexicon lexicon = new();
        Assert.Equal(~0, lexicon.Add([1, 2]));
        Assert.Equal(~1, lexicon.Add([3, 4]));
        lexicon.Clear();
        Assert.Equal(~0, lexicon.Add([3, 4]));
        Assert.Equal(~1, lexicon.Add([1, 2]));
    }

    private static void ExpectIdSet(IdSet expected, IdSet actual)
    {
        Assert.Equal(expected.Count, actual.Count);
        Assert.True(expected.SequenceEqual(actual));
    }
}
