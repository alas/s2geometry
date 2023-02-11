// S2CellRangeIterator is a wrapper around an S2CellIterator that caches the
// range_min() and range_max() that each cell covers as it iterates.  This lets
// us define range based methods such as SeekTo() and Locate() efficiently.
//
// Computing the range_max() and range_min() of a cell isn't expensive but it's
// not free either, so we extend the S2CellIterator interface instead of
// integrating this functionality there, allowing the user to pay only for what
// they use.
//
// An S2CellRangeIterator wraps an S2CellIterator, but is also itself an
// S2CellIterator and thus can be used anywhere one is required.
//
// Note(Alas): S2CellRangeEnumerator is a decorator that wraps S2CellEnumerator

namespace S2Geometry;

internal static class S2CellRangeEnumeratorFactory
{
    // Builds a new S2CellRangeIterator from an index, supporting type inference.
    //
    // We may wish to provide overloads for other types in the future, so we
    // disqualify this function from overload resolution using std::enable_if when
    // the type isn't an S2ShapeIndex.
    //
    // The index must live for the duration of the iterator, so we take it by const
    // pointer instead of reference to avoid binding to temporaries.
    //
    //template<typename IndexType, typename std::enable_if<S2ShapeIndex::ImplementedBy<IndexType>{}, bool>::type = true>
    public static S2CellRangeEnumerator<S2ShapeIndexCell> Make(S2ShapeIndex index) =>
        new(index.GetNewEnumerator(S2ShapeIndex.InitialPosition.BEGIN));

    // Builds a new S2CellRangeIterator from an S2CellIterator, explicitly supports
    // type inference for the Iterator parameter.
    //template<typename Iterator, typename std::enable_if<S2CellIterator::ImplementedBy<Iterator>{}, bool>::type = true>
    public static S2CellRangeEnumerator<T> Make<T>(S2CellEnumerator<T> enumerator) where T : S2ShapeIndex =>
        new(enumerator);
}

internal sealed class S2CellRangeEnumerator<T> : S2CellEnumerator<T>
{
    #region Fields and Properties

    public S2CellEnumerator<T> Enumerator { get; }

    // The min and max leaf cell ids covered by the current cell.  If done() is
    // true, these methods return a value larger than any valid cell id.
    public S2CellId RangeMin { get; private set; }
    public S2CellId RangeMax { get; private set; }

    #endregion

    #region Constructor

    // Construct a new S2CellRangeIterator positioned at the beginning.
    public S2CellRangeEnumerator(S2CellEnumerator<T> enumerator)
    {
        Enumerator = enumerator;
        Refresh();
    }

    #endregion

    // The current S2CellId and cell contents.
    public override S2CellId Id => Enumerator.Id;

    // Queries the relationship between two range iterators.  Returns -1 if this
    // iterator's current position entirely precedes the other iterator's current
    // position, +1 if it entirely follows, and 0 if they overlap.
    public int Relation(S2CellRangeEnumerator<T> b)
    {
        if (RangeMax < b.RangeMin) return -1;
        if (RangeMin > b.RangeMax) return +1;
        return 0;
    }

    #region S2CellIterator API

    public override bool MoveNext()
    {
        var status = Enumerator.MoveNext();
        Refresh();
        return status;
    }
    public override bool MovePrevious()
    {
        bool status = Enumerator.MovePrevious();
        Refresh();
        return status;
    }
    public override void Seek(S2CellId target)
    {
        Enumerator.Seek(target);
        Refresh();
    }
    public override void Finish()
    {
        Enumerator.Finish();
        Refresh();
    }
    public override bool Done() => Enumerator.Done();
    public override bool Locate(S2Point target)
    {
        bool status = Enumerator.Locate(target);
        Refresh();
        return status;
    }
    public override S2CellRelation Locate(S2CellId target)
    {
        // Let T be the target cell id, let I = Seek(T.range_min()) and let Prev(I) be
        // the predecessor of I.  If T contains any index cells, then T contains I.
        // Similarly, if T is contained by an index cell, then the containing cell is
        // either I or Prev(I).  We test for containment by comparing the ranges of
        // leaf cells spanned by T, I, and Prev(I).
        Seek(target.RangeMin());
        if (!Done())
        {
            // The target is contained by the cell we landed on, so it's indexed.
            if (Id >= target && RangeMin <= target)
            {
                return S2CellRelation.INDEXED;
            }

            // The cell we landed on is contained by the target, so it's subdivided.
            if (Id <= target.RangeMax())
            {
                return S2CellRelation.SUBDIVIDED;
            }
        }

        // Otherwise check the previous cell (if it exists).  If it contains the
        // target then it's indexed, otherwise the target cell is disjoint.
        if (MovePrevious() && RangeMax >= target)
        {
            return S2CellRelation.INDEXED;
        }
        return S2CellRelation.DISJOINT;
    }

    public override void SetPosition(int position)
    {
        Enumerator.SetPosition(position);
        Refresh();
    }

    public override void Reset() => throw new NotImplementedException();

    public override void Dispose() => throw new NotImplementedException();

    #endregion

    // The same as above, but uses another S2CellRangeIterator as the target.
    //
    // Convenience re-implementation of the above function, see it for details.
    public S2CellRelation Locate(S2CellRangeEnumerator<T> target)
    {
        Seek(target.RangeMin);
        if (!Done())
        {
            // The target is contained by the cell we landed on, so it's indexed.
            if (Id >= target.Id && RangeMin <= target.Id)
            {
                return S2CellRelation.INDEXED;
            }

            // The cell we landed on is contained by the target, so it's subdivided.
            if (Id <= target.RangeMax)
            {
                return S2CellRelation.SUBDIVIDED;
            }
        }

        // Otherwise check the previous cell (if it exists).  If it contains the
        // target then it's indexed, otherwise the target cell is disjoint.
        if (MovePrevious() && RangeMax >= target.Id)
        {
            return S2CellRelation.INDEXED;
        }
        return S2CellRelation.DISJOINT;
    }

    // Position the iterator at the first cell that overlaps or follows
    // "target", i.e. such that range_max() >= target.range_min().
    public void SeekTo(S2CellRangeEnumerator<T> target)
    {
        Seek(target.RangeMin);

        // If the current cell does not overlap "target", it is possible that the
        // previous cell is the one we are looking for.  This can only happen when
        // the previous cell contains "target" but has a smaller S2CellId.
        if (Done() || RangeMin > target.RangeMax)
        {
            if (MovePrevious() && RangeMax < target.Id)
            {
                MoveNext();
            }
        }
        Refresh();
    }

    // Position the iterator at the first cell that follows "target", i.e. the
    // first cell such that range_min() > target.range_max().
    public void SeekBeyond(S2CellRangeEnumerator<T> target)
    {
        Seek(target.RangeMax.Next());
        if (!Done() && RangeMin <= target.RangeMax)
        {
            MoveNext();
        }
        Refresh();
    }

    // Updates internal state after the iterator has been repositioned.
    //
    // This method is inline, but is only called by non-inline methods defined in
    // this file.  Putting the definition here enforces this requirement.
    private void Refresh()
    {
        if (Done())
        {
            RangeMin = S2CellId.Sentinel.RangeMin();
            RangeMax = S2CellId.Sentinel.RangeMax();
        }
        else
        {
            RangeMin = Id.RangeMin();
            RangeMax = Id.RangeMax();
        }
    }
}
