using System.Collections;

namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // RangeEnumerator is a wrapper over S2ShapeIndex.Iterator with extra methods
    // that are useful for merging the contents of two or more S2ShapeIndexes.
    public class RangeEnumerator : IReversableEnumerator<S2ShapeIndexIdCell>
    {
        private readonly S2ShapeIndex _index;
        private readonly S2ShapeIndex.Enumerator _it;
        // The min and max leaf cell ids covered by the current cell.  If done() is
        // true, these methods return a value larger than any valid cell id.
        public S2CellId RangeMin { get; private set; }
        public S2CellId RangeMax { get; private set; }

        // Construct a new RangeIterator positioned at the first cell of the index.
        public RangeEnumerator(S2ShapeIndex index)
        {
            _index = index;
            _it = new S2ShapeIndex.Enumerator(index);
            //Refresh();
        }

        // The current S2CellId and cell contents.
        public S2CellId Id => _it.Id;
        public S2ShapeIndexCell Cell => _it.Cell;

        // Various other convenience methods for the current cell.
        public S2ClippedShape Clipped() { return Cell.Clipped(0); }
        public int NumEdges() { return Clipped().NumEdges; }
        public bool ContainsCenter() { return Clipped().ContainsCenter; }

        public S2ShapeIndexIdCell Current => _it.Current;
        object IEnumerator.Current => _it.Current;

        public bool MoveNext() { var has = _it.MoveNext(); Refresh(); return has; }
        public bool MovePrevious() { var has = _it.MovePrevious(); Refresh(); return has; }
        public void Reset() => _it.Reset();
        public bool Done() => _it.Done();
        public void SetPosition(int position) => _it.SetPosition(position);
        public void Dispose() { GC.SuppressFinalize(this); }

        // Position the iterator at the first cell that overlaps or follows
        // "target", i.e. such that RangeMax >= target.RangeMin.
        public void SeekTo(RangeEnumerator target)
        {
            _it.SetPosition(_index.SeekCell(target.RangeMin).pos);
            // If the current cell does not overlap "target", it is possible that the
            // previous cell is the one we are looking for.  This can only happen when
            // the previous cell contains "target" but has a smaller S2CellId.
            if (_it.Done() || _it.Id.RangeMin() > target.RangeMax)
            {
                if (_it.MovePrevious() && _it.Id.RangeMax() < target.Id) _it.MoveNext();
            }
            Refresh();
        }

        // Position the iterator at the first cell that follows "target", i.e. the
        // first cell such that RangeMin > target.RangeMax.
        public void SeekBeyond(RangeEnumerator target)
        {
            _it.SetPosition(_index.SeekCell(target.RangeMax.Next()).pos);
            if (!_it.Done() && _it.Id.RangeMin() <= target.RangeMax)
            {
                _it.MoveNext();
            }
            Refresh();
        }

        // Updates internal state after the iterator has been repositioned.
        private void Refresh()
        {
            RangeMin = Id.RangeMin();
            RangeMax = Id.RangeMax();
        }
    }
}
