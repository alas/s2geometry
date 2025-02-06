namespace S2Geometry;

using Seq = List<long>;

public class SequenceLexiconTests
{
    [Fact]
    internal void Test_SequenceLexicon_int64()
    {
        SequenceLexicon<long> lex = new();
        Assert.Equal(0, lex.Add(new Seq { }));
        Assert.Equal(1, lex.Add(new Seq { 5 }));
        Assert.Equal(0, lex.Add(new Seq { }));
        Assert.Equal(2, lex.Add(new Seq { 5, 5 }));
        Assert.Equal(3, lex.Add(new Seq { 5, 0, -3 }));
        Assert.Equal(1, lex.Add(new Seq { 5 }));
        Assert.Equal(4, lex.Add(new Seq { 0x7fffffffffffffff }));
        Assert.Equal(3, lex.Add(new Seq { 5, 0, -3 }));
        Assert.Equal(0, lex.Add(new Seq { }));
        Assert.Equal(5, lex.Count());
        Assert.Equal([], lex.Sequence(0));
        Assert.Equal([5], lex.Sequence(1));
        Assert.Equal([5, 5], lex.Sequence(2));
        Assert.Equal([5, 0, -3], lex.Sequence(3));
        Assert.Equal([0x7fffffffffffffff], lex.Sequence(4));
    }

    [Fact]
    internal void Test_SequenceLexicon_Clear()
    {
        SequenceLexicon<long> lex = new();
        Assert.Equal(0, lex.Add(new Seq { 1 }));
        Assert.Equal(1, lex.Add(new Seq { 2 }));
        lex.Clear();
        Assert.Equal(0, lex.Add(new Seq { 2 }));
        Assert.Equal(1, lex.Add(new Seq { 1 }));
    }

    [Fact]
    internal void Test_SequenceLexicon_CopyConstructor()
    {
        var original = new SequenceLexicon<long>();
        Assert.Equal(0, original.Add(new Seq { 1, 2 }));
        var lex = original;
        Assert.Equal(1, lex.Add(new Seq { 3, 4 }));
        ExpectSequence([1, 2], lex.Sequence(0));
        ExpectSequence([3, 4], lex.Sequence(1));
    }

    [Fact]
    internal void Test_SequenceLexicon_MoveConstructor()
    {
        var original = new SequenceLexicon<long>();
        Assert.Equal(0, original.Add(new Seq { 1, 2}));
        var lex = original;
        Assert.Equal(1, lex.Add(new Seq { 3, 4}));
        ExpectSequence([1, 2], lex.Sequence(0));
        ExpectSequence([3, 4], lex.Sequence(1));
    }

    [Fact]
    internal void Test_SequenceLexicon_CopyAssignmentOperator()
    {
        SequenceLexicon<long> original = new();
        Assert.Equal(0, original.Add(new Seq { 1, 2}));
        SequenceLexicon<long> lex = new();
        Assert.Equal(0, lex.Add(new Seq { 3, 4}));
        Assert.Equal(1, lex.Add(new Seq { 5, 6}));
        lex = original;
        Assert.Equal(1, lex.Add(new Seq { 7, 8}));
        ExpectSequence([1, 2], lex.Sequence(0));
        ExpectSequence([7, 8], lex.Sequence(1));
    }

    [Fact]
    internal void Test_SequenceLexicon_MoveAssignmentOperator()
    {
        SequenceLexicon<long> original = new();
        Assert.Equal(0, original.Add(new Seq { 1, 2}));
        SequenceLexicon<long> lex = new();
        Assert.Equal(0, lex.Add(new Seq { 3, 4}));
        Assert.Equal(1, lex.Add(new Seq { 5, 6}));
        lex = original;
        Assert.Equal(1, lex.Add(new Seq { 7, 8}));
        ExpectSequence([1, 2], lex.Sequence(0));
        ExpectSequence([7, 8], lex.Sequence(1));
    }

    private static void ExpectSequence<T>(List<T> expected, List<T> actual) {
        Assert.Equal(expected.Count, actual.Count);
        Assert.True(expected.SequenceEqual(actual));
    }
}
