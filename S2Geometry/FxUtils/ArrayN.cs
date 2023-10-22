namespace S2Geometry;

using System.Collections;

// Ensure that we don't usually need to allocate memory.
public class Array16<T> : ArrayN<T> { public Array16() : base() { } public Array16(Func<T> getDefault) : base(getDefault) { } public override int Size => 16; }
public class Array8<T> : ArrayN<T> { public Array8() : base() { } public Array8(Func<T> getDefault) : base(getDefault) { } public override int Size => 8; }
public class Array6<T> : ArrayN<T> { public Array6() : base() { } public Array6(Func<T> getDefault) : base(getDefault) { } public override int Size => 6; }
public class Array4<T> : ArrayN<T> { public Array4() : base() { } public Array4(Func<T> getDefault) : base(getDefault) { } public override int Size => 4; }
public class Array2<T> : ArrayN<T> { public Array2() : base() { } public Array2(Func<T> getDefault) : base(getDefault) { } public override int Size => 2; }
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
            MyDebug.Assert(i < Size);
            return arr[i];
        }
        set
        {
            MyDebug.Assert(i < Size);
            arr[i] = value;
        }
    }

    public void Clear() => Count = 0;

    public void Add(T item)
    {
        MyDebug.Assert(Count < Size);
        arr[Count++] = item;
    }

    public void Reserve(int index)
    {
        MyDebug.Assert(index + 1 <= Size);
        MyDebug.Assert(index <= Count);
        if (Count < (index + 1)) Count++;
    }

    public IEnumerator<T> GetEnumerator() => new ArrayNEnumerator(this);

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    private class ArrayNEnumerator(ArrayN<T> arrayn) : IEnumerator<T>
    {
        private int position = -1;

        public T Current => arrayn[position];

        object IEnumerator.Current => Current;

        public void Dispose() { }
        public bool MoveNext() => ++position >= 0 && position < arrayn.Size;
        public void Reset() => position = -1;
    }
}
