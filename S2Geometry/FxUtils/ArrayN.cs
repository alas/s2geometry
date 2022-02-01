using System.Collections;

namespace S2Geometry;

// Ensure that we don't usually need to allocate memory.
public class Array16<T> : ArrayN<T> { public Array16() : base() { } public Array16(Func<T> getDefault) : base(getDefault) { } public override int Size { get { return 16; } } }
public class Array6<T> : ArrayN<T> { public Array6() : base() { } public Array6(Func<T> getDefault) : base(getDefault) { } public override int Size { get { return 6; } } }
public class Array2<T> : ArrayN<T> { public Array2() : base() { } public Array2(Func<T> getDefault) : base(getDefault) { } public override int Size { get { return 2; } } }
public abstract class ArrayN<T> : IEnumerable<T>
{
    private readonly T[] arr;
    public int Count = 0;
    public abstract int Size { get; }

    public ArrayN() => arr = new T[Size];

    public ArrayN(Func<T> getDefault) => arr = new T[Size].Fill(getDefault);

    public T this[int i]
    {
        get
        {
            Assert.True(i < Size);
            return arr[i];
        }
        set
        {
            Assert.True(i < Size);
            arr[i] = value;
        }
    }

    public void Clear() => Count = 0;

    public void Add(T item)
    {
        Assert.True(Count < Size);
        arr[Count++] = item;
    }

    public void Reserve(int index)
    {
        Assert.True(index + 1 <= Size);
        Assert.True(index <= Count);
        if (Count < (index + 1)) Count++;
    }

    public IEnumerator<T> GetEnumerator() => ((IEnumerable<T>)arr).GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}
