namespace S2Geometry;

using System.Collections;

// S2PointSpan represents a view of an S2Point array.  It is used to pass
// vertex arrays to functions that don't care about the actual array type
// (e.g. std::vector<S2Point> or S2Point[]).
//
// NOTE: S2PointSpan has an implicit constructor from any container type with
// data() and size() methods (such as std::vector and std::array).  Therefore
// you can use such containers as arguments for any S2PointSpan parameter.
public class S2PointSpan : IList<S2Point>
{
    private readonly IList<S2Point> Path;
    public S2PointSpan() => Path = [];
    public S2PointSpan(int capacity) => Path = new List<S2Point>(capacity);
    public S2PointSpan(IList<S2Point> path) => Path = path;

    public virtual S2Point this[int index] { get => Path[index]; set => Path[index] = value; }

    public virtual int Count => Path.Count;

    public bool IsReadOnly => Path.IsReadOnly;

    public void Add(S2Point item) => Path.Add(item);

    public void Clear() => Path.Clear();

    public bool Contains(S2Point item) => Path.Contains(item);

    public void CopyTo(S2Point[] array, int arrayIndex) => Path.CopyTo(array, arrayIndex);

    public IEnumerator<S2Point> GetEnumerator() => Path.GetEnumerator();

    public int IndexOf(S2Point item) => Path.IndexOf(item);

    public void Insert(int index, S2Point item) => Path.Insert(index, item);

    public bool Remove(S2Point item) => Path.Remove(item);

    public void RemoveAt(int index) => Path.RemoveAt(index);

    IEnumerator IEnumerable.GetEnumerator() => Path.GetEnumerator();

    public static implicit operator S2PointSpan(S2Point[] path) => new(path);
    public static implicit operator S2PointSpan(List<S2Point> path) => new(path);
}

// Like S2PointSpan, except that operator[] maps index values in the range
// [n, 2*n-1] to the range [0, n-1] by subtracting n (where n == size()).
// In other words, two full copies of the vertex array are available.  (This
// is a compromise between convenience and efficiency, since computing the
// index modulo "n" is surprisingly expensive.)
//
// This property is useful for implementing algorithms where the elements of
// the span represent the vertices of a loop.
public class S2PointLoopSpan : S2PointSpan
{
    public S2PointLoopSpan() : base() { }
    public S2PointLoopSpan(int capacity) : base(capacity) { }
    public S2PointLoopSpan(IList<S2Point> path) : base(path) { }

    // Like operator[], but allows index values in the range [0, 2*size()-1]
    // where each index i >= size() is mapped to i - size().
    public override S2Point this[int i]
    {
        get
        {
            MyDebug.Assert(i >= 0);
            MyDebug.Assert(i < 2 * this.Count);

            int j = i - this.Count;
            return base[j < 0 ? i : j];
        }
    }

    public static implicit operator S2PointLoopSpan(S2Point[] path) => new(path);
    public static implicit operator S2PointLoopSpan(List<S2Point> path) => new(path);
}
