using System;
using System.Collections.Generic;

namespace S2Geometry
{
    // SequenceLexicon is a class for compactly representing sequences of values
    // (e.g., tuples).  It automatically eliminates duplicates, and maps the
    // remaining sequences to sequentially increasing integer ids.  See also
    // ValueLexicon and IdSetLexicon.
    //
    // Each distinct sequence is mapped to a 32-bit integer.  The space used for
    // each sequence is approximately 11 bytes plus the memory needed to represent
    // the sequence elements.  For example, a sequence of three "double"s would
    // need about 11 + 3*8 = 35 bytes.  Note also that sequences are referred to
    // using 32-bit ids rather than 64-bit pointers.
    //
    // This class has the same thread-safety properties as "string": const methods
    // are thread safe, and non-const methods are not thread safe.
    //
    // Example usage:
    //
    //   SequenceLexicon<string> lexicon;
    //   vector<string> pets {"cat", "dog", "parrot"};
    //   uint pets_id = lexicon.Add(pets);
    //   S2_CHECK_EQ(pets_id, lexicon.Add(pets));
    //   string values;
    //   for (const auto& pet : lexicon.sequence(pets_id)) {
    //     values += pet;
    //   }
    //   S2_CHECK_EQ("catdogparrot", values);
    //
    public class SequenceLexicon<T> where T : IComparable<T>
    {
        private readonly List<T> values_ = new();
        private readonly List<int> begins_ = new() { 0 };
        private readonly Dictionary<ulong, int> id_set_ = new();

        public SequenceLexicon() { }

        // Clears all data from the lexicon.
        public void Clear()
        {
            values_.Clear();
            begins_.Clear();
            id_set_.Clear();
            begins_.Add(0);
        }

        // Add the given sequence of values to the lexicon if it is not already
        // present, and return its integer id.  Ids are assigned sequentially
        // starting from zero.
        public int Add(IEnumerable<T> values)
        {
            var hash = HashMix.GetHash(values);
            if (!id_set_.ContainsKey(hash))
            {
                values_.AddRange(values);
                begins_.Add(values_.Count);
                var id = begins_.Count - 2;
                id_set_.Add(hash, id);
                return id;
            }
            else
            {
                return id_set_[hash];
            }
        }

        public int Count()
        {
            return begins_.Count - 1;
        }

        // Return the value sequence with the given id.  This method can be used
        // with range-based for loops as follows:
        //   for (const auto& value : lexicon.sequence(id)) { ... }
        public List<T> Sequence(int id)
        {
            List<T> result = new();
            for (var i = begins_[id]; i < begins_[id + 1]; i++)
            {
                result.Add(values_[i]);
            }
            return result;
        }
    }
}
