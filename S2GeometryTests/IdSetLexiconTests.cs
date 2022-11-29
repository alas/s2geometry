namespace S2Geometry;

using IdSet = List<Int32>;
using Seq = List<Int32>;

public class IdSetLexiconTests
{
    [Fact]
    public void Test_IdSetLexicon_EmptySet()
    {
        IdSetLexicon lexicon = new();
        ExpectIdSet(new IdSet { }, lexicon.IdSet_(lexicon.Add(new Seq { })));
    }

    [Fact]
    public void Test_IdSetLexicon_SingletonSets()
    {
        IdSetLexicon lexicon = new();
        Assert.Equal(5, lexicon.Add(new Seq { 5 }));
        Assert.Equal(0, lexicon.Add(new Seq { 0, 0 }));
        Assert.Equal(1, IdSetLexicon.AddSingleton(1));
        var m = int.MaxValue;
        Assert.Equal(m, lexicon.Add(new Seq { m }));

        ExpectIdSet(new IdSet { 0 }, lexicon.IdSet_(0));
        ExpectIdSet(new IdSet { 1 }, lexicon.IdSet_(1));
        ExpectIdSet(new IdSet { 5 }, lexicon.IdSet_(5));
        ExpectIdSet(new IdSet { m }, lexicon.IdSet_(m));
    }

    [Fact]
    public void Test_IdSetLexicon_SetsAreSorted()
    {
        IdSetLexicon lexicon = new();
        Assert.Equal(~0, lexicon.Add(new Seq { 2, 5 }));
        Assert.Equal(~1, lexicon.Add(new Seq { 3, 2, 5 }));
        Assert.Equal(~0, lexicon.Add(new Seq { 5, 2 }));
        Assert.Equal(~1, lexicon.Add(new Seq { 5, 3, 2, 5 }));

        ExpectIdSet(new IdSet { 2, 5 }, lexicon.IdSet_(~0));
        ExpectIdSet(new IdSet { 2, 3, 5 }, lexicon.IdSet_(~1));
    }

    [Fact]
    public void Test_IdSetLexicon_Clear()
    {
        IdSetLexicon lexicon = new();
        Assert.Equal(~0, lexicon.Add(new Seq { 1, 2 }));
        Assert.Equal(~1, lexicon.Add(new Seq { 3, 4 }));
        lexicon.Clear();
        Assert.Equal(~0, lexicon.Add(new Seq { 3, 4 }));
        Assert.Equal(~1, lexicon.Add(new Seq { 1, 2 }));
    }

    private static void ExpectIdSet(IdSet expected, IdSet actual)
    {
        Assert.Equal(expected.Count, actual.Count);
        Assert.True(expected.SequenceEqual(actual));
    }
}
