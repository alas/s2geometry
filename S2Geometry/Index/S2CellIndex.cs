// S2CellIndex stores a collection of (cell_id, label) pairs.  The S2CellIds
// may be overlapping or contain duplicate values.  For example, an
// S2CellIndex could store a collection of S2CellUnions, where each
// S2CellUnion has its own label.
//
// Labels are 32-bit non-negative integers, and are typically used to map the
// results of queries back to client data structures.  Labels other than
// integers can be supported by using a ValueLexicon, which maintains a set of
// distinct labels and maps them to sequentially numbered integers.  For
// example, the following code uses strings as labels:
//
//   ValueLexicon<string> my_label_lexicon;
//   string label_str = ...;
//   cell_index.Add(cell_id, my_label_lexicon.Add(label_str));
//   ...
//   Int32 label = ...;
//   string label_str = my_label_lexicon.value(label);
//
// To build an S2CellIndex, call Add() for each (cell_id, label) pair, and
// then call the Build() method.  For example:
//
//   S2CellId[] contents = ...;
//   for (int i = 0; i < contents.Count; ++i) {
//     index.Add(contents[i], i /*label*/);
//   }
//   index.Build();
//
// There is also a convenience method that adds an S2CellUnion:
//
//     index.Add(cell_union, label);
//
// Note that the index is not dynamic; the contents of the index cannot be
// changed once it has been built.
//
// There are several options for retrieving data from the index.  The simplest
// is to use a built-in method such as GetIntersectingLabels (which returns
// the labels of all cells that intersect a given target S2CellUnion):
//
//   Label[] labels = index.GetIntersectingLabels(target_union);
//
// Alternatively, you can use an external class such as S2ClosestCellQuery,
// which computes the cell(s) that are closest to a given target geometry.
// For example, here is how to find all cells that are closer than
// "distance_limit" to a given target point:
//
//   S2ClosestCellQuery query(out index);
//   query.Options().set_max_distance(distance_limit);
//   S2ClosestCellQuery.PointTarget target(target_point);
//   foreach (var& result in query.FindClosestCells(out target)) {
//     // result.distance() is the distance to the target.
//     // result.cell_id() is the indexed S2CellId.
//     // result.label() is the integer label associated with the S2CellId.
//     DoSomething(target_point, result);
//   }
//
// Finally, you can access the index contents directly.  Internally, the index
// consists of a set of non-overlapping leaf cell ranges that subdivide the
// sphere and such that each range intersects a particular set of (cell_id,
// label) pairs.  Data is accessed using the following iterator types:
//
//   RangeIterator:
//    - used to seek and iterate over the non-overlapping leaf cell ranges.
//   NonEmptyRangeIterator:
//    - like RangeIterator, but skips ranges whose contents are empty.
//   ContentsIterator:
//    - iterates over the (cell_id, label) pairs that intersect a given range.
//
// Note that these are low-level, efficient types intended mainly for
// implementing new query classes.  Most clients should use either the
// built-in methods such as VisitIntersectingCells and GetIntersectingLabels,
// or a helper such as S2ClosestCellQuery or S2Closest*Query.CellUnionTarget.

using System.Collections;
using System.Diagnostics.CodeAnalysis;

namespace S2Geometry;

public class S2CellIndex
{
    // A special label indicating that ContentsIterator.done() is true.
    private const int kDoneContents = -1;

    // Default constructor.
    public S2CellIndex()
    {
        range_nodes_ = new List<RangeNode>();
        cell_tree_ = new List<CellNode>();
    }

    // Returns the number of (cell_id, label) pairs in the index.
    public int NumCells() => cell_tree_.Count;

    // Adds the given (cell_id, label) pair to the index.  Note that the index
    // is not valid until Build() is called.
    //
    // The S2CellIds in the index may overlap (including duplicate values).
    // Duplicate (cell_id, label) pairs are also allowed, although be aware that
    // S2ClosestCellQuery will eliminate such duplicates anyway.
    //
    // REQUIRES: cell_id.IsValid
    public void Add(S2CellId cell_id, Int32 label)
    {
        System.Diagnostics.Debug.Assert(cell_id.IsValid());
        System.Diagnostics.Debug.Assert(label >= 0);
        cell_tree_.Add(new CellNode(cell_id, label, -1));
    }

    // Convenience function that adds a collection of cells with the same label.
    public void Add(S2CellUnion cell_ids, Int32 label)
    {
        foreach (S2CellId cell_id in cell_ids)
        {
            Add(cell_id, label);
        }
    }

    // Constructs the index.  This method may only be called once.  No iterators
    // may be used until the index is built.
    public void Build()
    {
        // To build the cell tree and leaf cell ranges, we maintain a stack of
        // (cell_id, label) pairs that contain the current leaf cell.  This class
        // represents an instruction to push or pop a (cell_id, label) pair.
        //
        // If label >= 0, the (cell_id, label) pair is pushed on the stack.
        // If cell_id == S2CellId.Sentinel(), a pair is popped from the stack.
        // Otherwise the stack is unchanged but a RangeNode is still emitted.


        var deltas = new List<DeltaBuild>(2 * cell_tree_.Count + 2);
        // Create two deltas for each (cell_id, label) pair: one to add the pair to
        // the stack (at the start of its leaf cell range), and one to remove it from
        // the stack (at the end of its leaf cell range).
        foreach (var node in cell_tree_)
        {
            deltas.Add(new DeltaBuild(node.CellId.RangeMin(), node.CellId, node.Label));
            deltas.Add(new DeltaBuild(node.CellId.RangeMax().Next(), S2CellId.Sentinel, -1));
        }
        // We also create two special deltas to ensure that a RangeNode is emitted at
        // the beginning and end of the S2CellId range.
        deltas.Add(new DeltaBuild(S2CellId.Begin(S2.kMaxCellLevel), S2CellId.None, -1));
        deltas.Add(new DeltaBuild(S2CellId.End(S2.kMaxCellLevel), S2CellId.None, -1));
        deltas.Sort();

        // Now walk through the deltas to build the leaf cell ranges and cell tree
        // (which is essentially a permanent form of the "stack" described above).
        cell_tree_.Clear();
        range_nodes_.Capacity = deltas.Count;
        int contents = -1;
        for (int i = 0; i < deltas.Count;)
        {
            var start_id = deltas[i].StartId;
            // Process all the deltas associated with the current start_id.
            for (; i < deltas.Count && deltas[i].StartId == start_id; ++i)
            {
                if (deltas[i].Label >= 0)
                {
                    cell_tree_.Add(new CellNode(deltas[i].CellId, deltas[i].Label, contents));
                    contents = cell_tree_.Count - 1;
                }
                else if (deltas[i].CellId == S2CellId.Sentinel)
                {
                    contents = cell_tree_[contents].Parent;
                }
            }
            range_nodes_.Add(new RangeNode(start_id, contents));
        }
    }

    // Clears the index so that it can be re-used.
    public void Clear()
    {
        cell_tree_.Clear();
        range_nodes_.Clear();
    }

    // A function that is called with each (cell_id, label) pair to be visited.
    // The function may return false in order to indicate that no further
    // (cell_id, label) pairs are needed.
    public delegate bool CellVisitor(S2CellId cell_id, Int32 label);

    // Visits all (cell_id, label) pairs in the given index that intersect the
    // given S2CellUnion "target", terminating early if the given CellVisitor
    // function returns false (in which case VisitIntersectingCells returns false
    // as well).  Each (cell_id, label) pair in the index is visited at most
    // once.  (If the index contains duplicates, then each copy is visited.)
    public bool VisitIntersectingCells(S2CellUnion target, CellVisitor visitor)
    {
        if (!target.CellIds.Any()) return true;

        var contents = new ContentsEnumerator(this);

        var targetIndex = 0;
        var targetLimit = target.CellIds.Count;
        var rangeIndex = 0;
        do
        {
            var limitId = range_nodes_[rangeIndex + 1].StartId;
            if (limitId <= target.CellIds[targetIndex].RangeMin())
            {
                // Only seek when necessary.
                var rangeTarget = target.CellIds[targetIndex].RangeMin();
                rangeIndex = range_nodes_.GetUpperBound(new RangeNode(rangeTarget, -1)) - 1;
            }
            for (; range_nodes_[rangeIndex].StartId <= target.CellIds[targetIndex].RangeMax(); rangeIndex++)
            {
                var rNode = range_nodes_[rangeIndex];
                contents.StartUnion(rNode);
                while (contents.MoveNext())
                {
                    if (!visitor(contents.Current.CellId, contents.Current.Label))
                    {
                        return false;
                    }
                }
            }
            // Check whether the next target cell is also contained by the leaf cell
            // range that we just processed.  If so, we can skip over all such cells
            // using binary search.  This speeds up benchmarks by between 2x and 10x
            // when the average number of intersecting cells is small (< 1).
            if (++targetIndex < targetLimit && target.CellIds[targetIndex].RangeMax() < range_nodes_[rangeIndex].StartId)
            {
                // Skip to the first target cell that extends past the previous range.
                targetIndex = target.CellIds.GetLowerBound(range_nodes_[rangeIndex].StartId, targetIndex + 1, targetLimit);
                if (target.CellIds[targetIndex - 1].RangeMax() >= range_nodes_[rangeIndex].StartId) --targetIndex;
            }
        } while (targetIndex < targetLimit);
        return true;
    }

    // Convenience function that returns the labels of all indexed cells that
    // intersect the given S2CellUnion "target".
    public SortedSet<Int32> GetIntersectingLabels(S2CellUnion target)
    {
        var labels = new SortedSet<Int32>();
        VisitIntersectingCells(target, (S2CellId cell_id, Int32 label) =>
        {
            labels.Add(label);
            return true;
        });
        return labels;
    }

    // A tree of (cell_id, label) pairs such that if X is an ancestor of Y, then
    // X.cell_id contains Y.cell_id.  The contents of a given range of leaf
    // cells can be represented by pointing to a node of this tree.
    private readonly List<CellNode> cell_tree_;

    // The last element of range_nodes_ is a sentinel value, which is necessary
    // in order to represent the range covered by the previous element.
    private readonly List<RangeNode> range_nodes_;

    // Represents a node in the (cell_id, label) tree.  Cells are organized in a
    // tree such that the ancestors of a given node contain that node.
    public readonly struct CellNode
    {
        public readonly S2CellId CellId { get; init; }
        public readonly Int32 Parent { get; init; }
        public readonly Int32 Label { get; init; }

        public static readonly CellNode Zero = new(S2CellId.None, kDoneContents, -1);

        public CellNode(S2CellId cell_id, Int32 label, Int32 parent)
        { CellId = cell_id; Label = label; Parent = parent; }
    }

    // An iterator that visits the entire set of indexed (cell_id, label) pairs
    // in an unspecified order.
    public IEnumerable<CellNode> GetCellEnumerable() => cell_tree_.AsEnumerable();

    #region RangeNode IEnumerator

    // An IEnumerator that seeks and iterates over a set of non-overlapping leaf
    // cell ranges that cover the entire sphere.  The indexed (s2cell_id, label)
    // pairs that intersect the current leaf cell range can be visited using
    // ContentsIterator (see below).  The last element is a sentinel value.
    public RangeNodeEnumerator GetRangeNodeEnumerator() => new(range_nodes_, -1);

    public class RangeNodeEnumerator : IReversableEnumerator<RangeNode>, ICustomCloneable
    {
        private readonly List<RangeNode> range_nodes_;
        protected int Position;

        // Initializes a RangeIterator for the given S2CellIndex.  The iterator is
        // initially *unpositioned*; you must call a positioning method such as
        // Begin() or Seek() before accessing its contents.
        public RangeNodeEnumerator(S2CellIndex index)
        {
            range_nodes_ = index.range_nodes_;
            Position = -1;
            System.Diagnostics.Debug.Assert(range_nodes_.Any()); // Call Build() first.
        }

        public RangeNodeEnumerator(List<RangeNode> rangeNodes, int position)
        {
            range_nodes_ = rangeNodes;
            Position = position;
            System.Diagnostics.Debug.Assert(range_nodes_.Any()); // Call Build() first.
        }

        public virtual bool MoveNext()
        {
            Position++;
            return !IsAtSentinel();
        }

        // If advancing the iterator "n" times would leave it positioned on a
        // valid range, does so and returns true.  Otherwise leaves the iterator
        // unmodified and returns false.
        public bool Advance(int n)
        {
            // Note that the last element of range_nodes_ is a sentinel value.
            if (n >= range_nodes_.Count - 1 - Position) return false;
            Position += n;
            return true;
        }

        // Positions the iterator so that done() is true.
        // Note that the last element of range_nodes_ is a sentinel value.
        public void Finish() => Position = range_nodes_.Count - 1;

        /// <remarks>The last item in range_nodes_ is a Sentinel</remarks>
        protected bool IsAtSentinel() => Position >= (range_nodes_.Count - 1);

        public void Reset() => Position = -1;
        public bool Done() => IsAtSentinel();

        public void SetPosition(int position) => Position = position;

        public S2CellId GetLimitId() => range_nodes_[Position + 1].StartId;

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        public virtual bool MovePrevious()
        {
            if (Position == 0) return false;

            Position--;
            return Position >= 0;
        }

        // Positions the iterator at the range containing "target". (Such a range
        // always exists as long as the target is a valid leaf cell.)
        //
        // REQUIRES: target.IsLeaf
        public void Seek(S2CellId target)
        {
            System.Diagnostics.Debug.Assert(target.IsLeaf());
            Position = range_nodes_.GetUpperBound(new RangeNode(target, -1)) - 1;
        }

        object IEnumerator.Current => Current;

        public RangeNode Current
        {
            get
            {
                try
                {
                    return range_nodes_[Position];
                }
                catch (IndexOutOfRangeException)
                {
                    throw new InvalidOperationException();
                }
            }
        }

        public void Dispose() { GC.SuppressFinalize(this); }

        public object CustomClone() => new RangeNodeEnumerator(range_nodes_, Position);
    }

    #endregion

    #region NonEmpty RangeNode IEnumerator

    // Like RangeNodeEnumerable, but only visits leaf cell ranges that overlap at
    // least one (cell_id, label) pair.
    public NonEmptyRangeEnumerator GetNERNEnum() => new(this);

    // Like RangeNodeEnumerable, but only visits leaf cell ranges that overlap at
    // least one (cell_id, label) pair.
    public class NonEmptyRangeEnumerator : RangeNodeEnumerator
    {
        public NonEmptyRangeEnumerator(S2CellIndex index) : base(index) { }

        // Advances the iterator to the next non-empty range of leaf cells.
        public override bool MoveNext()
        {
            do
            {
                Position++;
            } while (Current.IsEmpty && !IsAtSentinel());

            return !IsAtSentinel();
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous non-empty entry and
        // returns true.
        public override bool MovePrevious()
        {
            while (Position > 0)
            {
                Position--;
                if (!Current.IsEmpty) return true;
            }
            // Return to original position.
            if (Current.IsEmpty && !IsAtSentinel())
            {
                MoveNext();
            }
            return false;
        }
    }

    #endregion

    #region Contents IEnumerator

    // An iterator that visits the (cell_id, label) pairs that cover a set of
    // leaf cell ranges (see RangeIterator).  Note that when multiple leaf cell
    // ranges are visited, this class only guarantees that each result will be
    // reported at least once, i.e. duplicate values may be suppressed.  If you
    // want duplicate values to be reported again, be sure to call Clear() first.
    //
    // [In particular, the implementation guarantees that when multiple leaf
    // cell ranges are visited in monotonically increasing order, then each
    // (cell_id, label) pair is reported exactly once.]
    public class ContentsEnumerator : IEnumerator<CellNode>
    {
        #region Fields, Constants

        // A pointer to the cell tree itself (owned by the S2CellIndex).
        private readonly List<CellNode> cell_tree_;

        // The value of it.start_id() from the previous call to StartUnion().
        // This is used to check whether these values are monotonically
        // increasing.
        private S2CellId prev_start_id_;

        // The maximum index within the cell_tree_ vector visited during the
        // previous call to StartUnion().  This is used to eliminate duplicate
        // values when StartUnion() is called multiple times.
        private Int32 node_cutoff_;

        // The maximum index within the cell_tree_ vector visited during the
        // current call to StartUnion().  This is used to update node_cutoff_.
        private Int32 next_node_cutoff_;

        public CellNode Current { get; private set; }

        object IEnumerator.Current => Current;

        #endregion

        #region Constructors

        // Convenience constructor that calls Init().
        public ContentsEnumerator(S2CellIndex index)
        {
            cell_tree_ = index.cell_tree_;
            Reset();
        }

        #endregion

        #region IEnumerator

        // Clears all state with respect to which range(s) have been visited.
        public void Reset()
        {
            Current = CellNode.Zero;
            prev_start_id_ = S2CellId.None;
            node_cutoff_ = -1;
            next_node_cutoff_ = -1;
        }

        // Advances the iterator to the next (cell_id, label) pair covered by the
        // current leaf cell range.
        // REQUIRES: !done()
        public bool MoveNext()
        {
            if (Current.Parent <= node_cutoff_)
            {
                // We have already processed this node and its ancestors.
                node_cutoff_ = next_node_cutoff_;
                return false;
            }
            else
            {
                Current = cell_tree_[Current.Parent];
                return true;
            }
        }

        public void Dispose() { GC.SuppressFinalize(this); }

        #endregion

        #region ContentsEnumerator

        // Positions the ContentsIterator at the first (cell_id, label) pair that
        // covers the given leaf cell range.  Note that when multiple leaf cell
        // ranges are visited using the same ContentsIterator, duplicate values
        // may be suppressed.  If you don't want this behavior, call Clear() first.
        public bool StartUnion(RangeNode rn)
        {
            if (rn.StartId < prev_start_id_)
            {
                node_cutoff_ = -1;  // Can't automatically eliminate duplicates.
            }
            prev_start_id_ = rn.StartId;

            // TODO(ericv): Since RangeNode only uses 12 of its 16 bytes, we could add a
            // "label" field without using any extra space.  Then we could store a leaf
            // node of cell_tree_ directly in each RangeNode, where the cell_id is
            // implicitly defined as the one that covers the current leaf cell range.
            // This would save quite a bit of space; e.g. if the given cells are
            // non-overlapping, then cell_tree_ would be empty (since every node is a
            // leaf node and could therefore be stored directly in a RangeNode).  It
            // would also be faster because cell_tree_ would rarely be accessed.
            if (rn.Contents <= node_cutoff_)
            {
                return false;
            }
            else
            {
                Current = cell_tree_[rn.Contents];
            }

            // When visiting ancestors, we can stop as soon as the node index is smaller
            // than any previously visited node index.  Because indexes are assigned
            // using a preorder traversal, such nodes are guaranteed to have already
            // been reported.
            next_node_cutoff_ = rn.Contents;
            return true;
        }

        #endregion
    }

    #endregion

    // Convenience class that represents a (cell_id, label) pair.
    public readonly struct LabelledCell : IEquatable<LabelledCell>, IComparable<LabelledCell>
    {
        #region Fields, Constants

        public readonly S2CellId CellId { get; init; }
        public readonly Int32 Label { get; init; }

        public static readonly LabelledCell Zero = new(S2CellId.None, -1);

        #endregion

        #region Constructors

        public LabelledCell(S2CellId cell_id, Int32 label)
        { CellId = cell_id; Label = label; }

        #endregion

        #region IEquatable

        public bool Equals(LabelledCell y) => CellId == y.CellId && Label == y.Label;
        public override bool Equals(object? obj) => obj is LabelledCell cell && Equals(cell);
        public override int GetHashCode() => HashCode.Combine(CellId.Id, Label);
        public static bool operator ==(LabelledCell x, LabelledCell y) => Equals(x, y);
        public static bool operator !=(LabelledCell x, LabelledCell y) => !Equals(x, y);

        #endregion

        #region IComparable

        public static bool operator <(LabelledCell x, LabelledCell y)
        {
            if (x.CellId < y.CellId) return true;
            if (y.CellId < x.CellId) return false;
            return x.Label < y.Label;
        }

        public static bool operator >(LabelledCell x, LabelledCell y)
        {
            if (x.CellId > y.CellId) return true;
            if (y.CellId > x.CellId) return false;
            return x.Label > y.Label;
        }
        public int CompareTo([AllowNull] LabelledCell other)
        {
            if (CellId.CompareTo(other.CellId) != 0)
                return CellId.CompareTo(other.CellId);

            return Label.CompareTo(other.Label);
        }

        #endregion

        #region Object

        public override string ToString() => $"({CellId}, {Label})";

        #endregion
    }


    private readonly struct DeltaBuild : IComparable<DeltaBuild>
    {
        public readonly S2CellId StartId { get; init; }
        public readonly S2CellId CellId { get; init; }
        public readonly int Label { get; init; }

        public DeltaBuild(S2CellId start_id, S2CellId cell_id, Int32 label)
        { StartId = start_id; CellId = cell_id; Label = label; }

        // Deltas are sorted first by start_id, then in reverse order by cell_id,
        // and then by label.  This is necessary to ensure that (1) larger cells
        // are pushed on the stack before smaller cells, and (2) cells are popped
        // off the stack before any new cells are added.
        public int CompareTo(DeltaBuild other)
        {
            if (StartId.CompareTo(other.StartId) != 0)
            {
                return StartId.CompareTo(other.StartId);
            }
            else if (other.CellId.CompareTo(CellId) != 0)
            {
                return other.CellId.CompareTo(CellId);
            }
            return Label.CompareTo(other.Label);
        }
    }

    // A RangeNode represents a range of leaf S2CellIds.  The range starts at
    // "start_id" (a leaf cell) and ends at the "start_id" field of the next
    // RangeNode.  "contents" points to the node of cell_tree_ representing the
    // cells that overlap this range.
    public readonly struct RangeNode : IComparable<RangeNode>
    {
        public readonly S2CellId StartId { get; init; }  // First leaf cell contained by this range.
        public readonly int Contents { get; init; }     // Contents of this node (an index within cell_tree_).

        public RangeNode(S2CellId start_id, int contents)
        {
            StartId = start_id; Contents = contents;
        }

        // Comparison operator needed for upper_bound().
        public int CompareTo(RangeNode other) => StartId.CompareTo(other.StartId);

        // Returns true if no (s2cell_id, label) pairs intersect this range.
        // Also returns true if it is the last (sentinel).
        public bool IsEmpty => Contents == kDoneContents;

        public override string ToString() => $"{Contents} - {StartId}";
    }
}
