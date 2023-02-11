// An abstract base class for iterators over any sorted collection keyed by
// S2CellId.
//
// This is intentionally opaque to the type of the value we might be mapping to,
// which can be anything (index cells, points, integers, other S2CellIds, etc),
// and instead only defines the common subset of functionality for positioning
// the iterator and querying the current cell id.
//
// This class is generally not used directly via pointer but instead to type
// check that a given template parameter implements the interface.
//
// Several default implementations for methods like SeekTo and SeekPast are
// given as static methods.  Inheritors should call them directly when defining
// their overloads of the API and mark their methods final to ensure
// de-virtualization when using sub-classes directly.
//
// A canonical implementation for SeekTo thus might look like:
//
//   class MyIterator : public S2CellIterator {
//    public:
//     void Locate(const S2Point& point) final {
//       return LocateImpl(*this, point);
//     }
//   };
//
// This will ensure code that uses MyIterator directly and calls Locate will
// directly use the methods of MyIterator instead of calling through the vtable.

namespace S2Geometry;

using System.Collections;

public abstract class S2CellEnumerator<T> : IReversableEnumerator<T>
{
    // Returns the current S2CellId that the iterator is positioned at.  This
    // function should be cheap to call (ideally directly returning the current
    // value).  When the iterator is done, this should return S2CellId::Sentinel.
    public abstract S2CellId Id { get; }

    public virtual T Current => throw new NotImplementedException();

    object IEnumerator.Current => throw new NotImplementedException();

    // Return true if the iterator has reached the end of the input. This function
    // should be cheap to call to check if iteration has ended.
    public abstract bool Done();

    // Positions the iterator at the next value.  Must not be called when done()
    // is true.
    public abstract bool MoveNext();

    // Positions the iterator at the previous value.  Returns false if the
    // iterator is already at the start.
    public abstract bool MovePrevious();

    // Seeks the iterator to the first cell with id() >= target or the end
    // of the iterator if no such cell exists.
    public abstract void Seek(S2CellId target);

    // Positions the iterator at the cell containing target and returns true. If
    // no such cell exists, return false and leave the iterator in an undefined
    // (but valid) state.
    public abstract bool Locate(S2Point target);

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and returns
    // INDEXED.  Otherwise if T contains one or more (smaller) index cells, it
    // positions the iterator at the first such cell I and returns SUBDIVIDED.
    // Otherwise it returns DISJOINT and leaves the iterator in an undefined
    // (but valid) state.
    public abstract S2CellRelation Locate(S2CellId target);

    // Positions the iterator past the last value.  After calling this function,
    // the done() method should return true.
    public abstract void Finish();

    public abstract void SetPosition(int position);

    /// <summary>
    /// Begin
    /// </summary>
    public abstract void Reset();

    public abstract void Dispose();

    protected bool LocateImpl<TEnumerator>(TEnumerator iter, S2Point point)
        where TEnumerator : S2CellEnumerator<T>
    {
        //static_assert(S2CellIterator::ImplementedBy<Iterator>{ }, "Iterator must implement the S2CellIterator API.");

        // Let I = Seek(T), where T is the leaf cell containing the target point, and
        // let Prev(I) be the predecessor of I.  If T is contained by an index cell,
        // then the containing cell is either I or Prev(I).  We test for containment
        // by comparing the ranges of leaf cells spanned by T, I, and Prev(I).
        S2CellId target = new(point);

        iter.Seek(target);
        if (!iter.Done() && iter.Id.RangeMin() <= target)
        {
            return true;
        }

        if (iter.MovePrevious() && iter.Id.RangeMax() >= target)
        {
            return true;
        }
        return false;
    }

    protected S2CellRelation LocateImpl<TEnumerator>(TEnumerator iter, S2CellId target)
        where TEnumerator : S2CellEnumerator<T>
    {
        //static_assert(S2CellIterator::ImplementedBy<Iterator>{ }, "Iterator must implement the S2CellIterator API.");

        // Let T be the target cell id, let I = Seek(T.range_min()) and let Prev(I) be
        // the predecessor of I.  If T contains any index cells, then T contains I.
        // Similarly, if T is contained by an index cell, then the containing cell is
        // either I or Prev(I).  We test for containment by comparing the ranges of
        // leaf cells spanned by T, I, and Prev(I).
        iter.Seek(target.RangeMin());
        if (!iter.Done())
        {
            // The target is contained by the cell we landed on, so it's indexed.
            if (iter.Id >= target && iter.Id.RangeMin() <= target)
            {
                return S2CellRelation.INDEXED;
            }

            // The cell we landed on is contained by the target, so it's subdivided.
            if (iter.Id <= target.RangeMax())
            {
                return S2CellRelation.SUBDIVIDED;
            }
        }

        // Otherwise check the previous cell (if it exists).  If it contains the
        // target then it's indexed, otherwise the target cell is disjoint.
        if (iter.MovePrevious() && iter.Id.RangeMax() >= target)
        {
            return S2CellRelation.INDEXED;
        }
        return S2CellRelation.DISJOINT;
    }
}

// Possible relationships between two S2CellIds in an index.
public enum S2CellRelation
{
    INDEXED,       // Target is contained by an index cell
    SUBDIVIDED,    // Target is subdivided into one or more index cells
    DISJOINT       // Target does not intersect any index cells
};
