// IdSetLexicon is a class for compactly representing sets of non-negative
// integers such as array indices ("id sets").  It is especially suitable when
// either (1) there are many duplicate sets, or (2) there are many singleton
// or empty sets.  See also ValueLexicon and SequenceLexicon.
//
// Each distinct id set is mapped to a 32-bit integer.  Empty and singleton
// sets take up no additional space whatsoever; the set itself is represented
// by the unique id assigned to the set. Sets of size 2 or more occupy about
// 11 bytes per set plus 4 bytes per element (as compared to 24 bytes per set
// plus 4 bytes per element for std.vector).  Duplicate sets are
// automatically eliminated.  Note also that id sets are referred to using
// 32-bit integers rather than 64-bit pointers.
//
// This class is especially useful in conjunction with ValueLexicon<T>.  For
// example, suppose that you want to label objects with a set of strings.  You
// could use a ValueLexicon<string> to map the strings to "label ids" (32-bit
// integers), and then use IdSetLexicon to map each set of labels to a "label
// set id".  Each reference to that label set then takes up only 4 bytes.
//
// Example usage:
//
//   ValueLexicon<string> labels_;
//   IdSetLexicon label_sets_;
//
//   int32 GetLabelSet(const vector<string>& label_strings) {
//     vector<int32> label_ids;
//     for (const var& str : label_strings) {
//       label_ids.push_back(labels_.Add(str));
//     }
//     return label_sets_.Add(label_ids);
//   }
//
//   int label_set_id = GetLabelSet(...);
//   for (var id : label_sets_.id_set(label_set_id)) {
//     S2_LOG(INFO) << id;
//   }
//
// This class is similar to SequenceLexicon, except:
//
// 1. Empty and singleton sets are represented implicitly; they use no space.
// 2. Sets are represented rather than sequences; the ordering of values is
//    not important and duplicates are removed.
// 3. The values must be 32-bit non-negative integers (only).

namespace S2Geometry;

public class IdSetLexicon
{
    public const int kEmptySetId = int.MinValue;
    private readonly SequenceLexicon<int> id_sets_ = new();

    // Add the given set of integers to the lexicon if it is not already
    // present, and return the unique id for this set.  The values are automatically sorted and
    // duplicates are removed.  Returns a signed integer representing this set.
    public int Add(List<int> ids)
    {
        if (!ids.Any())
        {
            // Empty sets have a special id chosen not to conflict with other ids.
            return kEmptySetId;
        }
        else if (ids.Count == 1)
        {
            // Singleton sets are represented by their element.
            return ids.First();
        }
        else
        {
            var ss = new SortedSet<int>(ids).ToList();

            // After eliminating duplicates, we may now have a singleton.
            if (ss.Count == 1) return ss.First();

            // Non-singleton sets are represented by the bitwise complement of the id
            // returned by SequenceLexicon.
            return ~id_sets_.Add(ss);
        }
    }

    public List<int> IdSet_(int set_id)
    {
        if (set_id >= 0)
        {
            return new List<int> { set_id };
        }
        else if (set_id == kEmptySetId)
        {
            return new List<int> { };
        }
        else
        {
            return id_sets_.Sequence(~set_id);
        }
    }

    public void Clear() => id_sets_.Clear();

    public static int AddSingleton(int id)
    {
        System.Diagnostics.Debug.Assert(id >= 0);
        System.Diagnostics.Debug.Assert(id <= int.MaxValue);
        // Singleton sets are represented by their element.
        return id;
    }
}
