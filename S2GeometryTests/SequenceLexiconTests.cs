using System;
using System.Collections.Generic;
using System.Linq;
using Xunit;
using Seq = System.Collections.Generic.List<System.Int64>;

namespace S2Geometry
{
    public class SequenceLexiconTests
    {
        [Fact]
        public void Test_SequenceLexicon_int64()
        {
            SequenceLexicon<Int64> lex = new();
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
            Assert.Equal(new Seq { }, lex.Sequence(0));
            Assert.Equal(new Seq { 5 }, lex.Sequence(1));
            Assert.Equal(new Seq { 5, 5 }, lex.Sequence(2));
            Assert.Equal(new Seq { 5, 0, -3 }, lex.Sequence(3));
            Assert.Equal(new Seq { 0x7fffffffffffffff }, lex.Sequence(4));
        }

        [Fact]
        public void Test_SequenceLexicon_Clear()
        {
            SequenceLexicon<Int64> lex = new();
            Assert.Equal(0, lex.Add(new Seq { 1 }));
            Assert.Equal(1, lex.Add(new Seq { 2 }));
            lex.Clear();
            Assert.Equal(0, lex.Add(new Seq { 2 }));
            Assert.Equal(1, lex.Add(new Seq { 1 }));
        }

        [Fact]
        public void Test_SequenceLexicon_CopyConstructor()
        {
            var original = new SequenceLexicon<Int64>();
            Assert.Equal(0, original.Add(new Seq { 1, 2 }));
            var lex = original;
            Assert.Equal(1, lex.Add(new Seq { 3, 4 }));
            ExpectSequence(new Seq { 1, 2 }, lex.Sequence(0));
            ExpectSequence(new Seq { 3, 4 }, lex.Sequence(1));
        }

        [Fact]
        public void Test_SequenceLexicon_MoveConstructor()
        {
            var original = new SequenceLexicon<Int64>();
            Assert.Equal(0, original.Add(new Seq { 1, 2}));
            var lex = original;
            Assert.Equal(1, lex.Add(new Seq { 3, 4}));
            ExpectSequence(new Seq { 1, 2}, lex.Sequence(0));
            ExpectSequence(new Seq { 3, 4}, lex.Sequence(1));
        }

        [Fact]
        public void Test_SequenceLexicon_CopyAssignmentOperator()
        {
            SequenceLexicon<Int64> original = new();
            Assert.Equal(0, original.Add(new Seq { 1, 2}));
            SequenceLexicon<Int64> lex = new();
            Assert.Equal(0, lex.Add(new Seq { 3, 4}));
            Assert.Equal(1, lex.Add(new Seq { 5, 6}));
            lex = original;
            Assert.Equal(1, lex.Add(new Seq { 7, 8}));
            ExpectSequence(new Seq { 1, 2}, lex.Sequence(0));
            ExpectSequence(new Seq { 7, 8}, lex.Sequence(1));
        }

        [Fact]
        public void Test_SequenceLexicon_MoveAssignmentOperator()
        {
            SequenceLexicon<Int64> original = new();
            Assert.Equal(0, original.Add(new Seq { 1, 2}));
            SequenceLexicon<Int64> lex = new();
            Assert.Equal(0, lex.Add(new Seq { 3, 4}));
            Assert.Equal(1, lex.Add(new Seq { 5, 6}));
            lex = original;
            Assert.Equal(1, lex.Add(new Seq { 7, 8}));
            ExpectSequence(new Seq { 1, 2}, lex.Sequence(0));
            ExpectSequence(new Seq { 7, 8}, lex.Sequence(1));
        }

        private static void ExpectSequence<T>(List<T> expected, List<T> actual) {
            Assert.Equal(expected.Count, actual.Count);
            Assert.True(expected.SequenceEqual(actual));
        }
    }
}
