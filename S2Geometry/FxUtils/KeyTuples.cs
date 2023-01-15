namespace S2Geometry;

/// <summary>
/// A tuple where Item1 is the identifier and has 2 Items in total
/// </summary>
public readonly record struct KeyData<T, U>(T Item1, U Item2)
    : IComparable<KeyData<T, U>>
     where T : IEquatable<T>, IComparable<T>
{
    public override string ToString() => $"{Item1}:{Item2}";

    public bool Equals(KeyData<T, U> other) => Equals(Item1, other.Item1);
    public override int GetHashCode() => Item1.GetHashCode();

    public int CompareTo(KeyData<T, U> other) => Item1.CompareTo(other.Item1);
    public static bool operator <(KeyData<T, U> x, KeyData<T, U> y) => x.CompareTo(y) < 0;
    public static bool operator >(KeyData<T, U> x, KeyData<T, U> y) => x.CompareTo(y) > 0;
    public static bool operator <=(KeyData<T, U> x, KeyData<T, U> y) => x.CompareTo(y) <= 0;
    public static bool operator >=(KeyData<T, U> x, KeyData<T, U> y) => x.CompareTo(y) >= 0;
}
/// <summary>
/// A tuple where Item1 is the identifier and has 3 Items in total
/// </summary>
public readonly record struct KeyData2<T, U, V>(T Item1, U Item2, V Item3)
    : IComparable<KeyData2<T, U, V>>
     where T : IEquatable<T>, IComparable<T>
{
    public override string ToString() => $"{Item1}:{Item2}:{Item3}";

    public bool Equals(KeyData2<T, U, V> other) => Equals(Item1, other.Item1);
    public override int GetHashCode() => Item1.GetHashCode();

    public int CompareTo(KeyData2<T, U, V> other) => Item1.CompareTo(other.Item1);
    public static bool operator <(KeyData2<T, U, V> x, KeyData2<T, U, V> y) => x.CompareTo(y) < 0;
    public static bool operator >(KeyData2<T, U, V> x, KeyData2<T, U, V> y) => x.CompareTo(y) > 0;
    public static bool operator <=(KeyData2<T, U, V> x, KeyData2<T, U, V> y) => x.CompareTo(y) <= 0;
    public static bool operator >=(KeyData2<T, U, V> x, KeyData2<T, U, V> y) => x.CompareTo(y) >= 0;
}
/// <summary>
/// KeyData2 with types tailored for S2PointIndex
/// </summary>
public readonly record struct TreeNode<T>(S2CellId Id, S2Point Point, T Data)
    : IComparable<TreeNode<T>>
{
    public override string ToString() => $"{Id}:{Point}:{Data}";

    public bool Equals(TreeNode<T> other) => Equals(Id, other.Id);
    public override int GetHashCode() => Id.GetHashCode();

    public int CompareTo(TreeNode<T> other) => Id.CompareTo(other.Id);
    public static bool operator <(TreeNode<T> x, TreeNode<T> y) => x.Id.CompareTo(y.Id) < 0;
    public static bool operator >(TreeNode<T> x, TreeNode<T> y) => x.Id.CompareTo(y.Id) > 0;
    public static bool operator <=(TreeNode<T> x, TreeNode<T> y) => x.Id.CompareTo(y.Id) <= 0;
    public static bool operator >=(TreeNode<T> x, TreeNode<T> y) => x.Id.CompareTo(y.Id) >= 0;
}
/// <summary>
/// A tuple where Item1 is the identifier and has 2 Items in total, Item1 is reverse compared
/// </summary>
public readonly record struct ReverseKeyData<T, U>(T Item1, U Item2)
    : IComparable<ReverseKeyData<T, U>>
     where T : IEquatable<T>, IComparable<T>
{
    public override string ToString() => $"{Item1}:{Item2}";

    public bool Equals(ReverseKeyData<T, U> other) => Equals(Item1, other.Item1);
    public override int GetHashCode() => Item1.GetHashCode();

    public int CompareTo(ReverseKeyData<T, U> other) => other.Item1.CompareTo(Item1);
    public static bool operator <(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y) => y.CompareTo(x) < 0;
    public static bool operator >(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y) => y.CompareTo(x) > 0;
    public static bool operator <=(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y) => y.CompareTo(x) <= 0;
    public static bool operator >=(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y) => y.CompareTo(x) >= 0;
}
/// <summary>
/// A tuple where Item1 is the identifier and has 3 Items in total, Item1 is reverse compared
/// </summary>
public readonly record struct ReverseKeyData2<T, U, V>(T Item1, U Item2, V Item3)
    : IComparable<ReverseKeyData2<T, U, V>>
     where T : IEquatable<T>, IComparable<T>
{
    public override string ToString() => $"{Item1}:{Item2}:{Item3}";

    public bool Equals(ReverseKeyData2<T, U, V> other) => Equals(Item1, other.Item1);
    public override int GetHashCode() => Item1.GetHashCode();

    public int CompareTo(ReverseKeyData2<T, U, V> other) => other.Item1.CompareTo(Item1);
    public static bool operator <(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y) => y.CompareTo(x) < 0;
    public static bool operator >(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y) => y.CompareTo(x) > 0;
    public static bool operator <=(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y) => y.CompareTo(x) <= 0;
    public static bool operator >=(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y) => y.CompareTo(x) >= 0;
}
/// <summary>
/// A tuple where Item1 and Item2 are the identifier and has 2 Items in total
/// </summary>
public readonly record struct KeyKey<T, U>(T Item1, U Item2)
    : IComparable<KeyKey<T, U>>
     where T : IEquatable<T>, IComparable<T>
     where U : IEquatable<U>, IComparable<U>
{
    public override string ToString() => $"{Item1}:{Item2}";

    public int CompareTo(KeyKey<T, U> other)
    {
        var c = Item1.CompareTo(other.Item1);
        if (c != 0) return c;

        return Item2.CompareTo(other.Item2);
    }
    public static bool operator <(KeyKey<T, U> x, KeyKey<T, U> y) => x.CompareTo(y) < 0;
    public static bool operator >(KeyKey<T, U> x, KeyKey<T, U> y) => x.CompareTo(y) > 0;
    public static bool operator <=(KeyKey<T, U> x, KeyKey<T, U> y) => x.CompareTo(y) <= 0;
    public static bool operator >=(KeyKey<T, U> x, KeyKey<T, U> y) => x.CompareTo(y) >= 0;
}
/// <summary>
/// A tuple where Item1(Int32) and Item2(Int32) are the identifier and has 2 Items in total
/// </summary>
public readonly record struct Int32Int32(Int32 Item1, Int32 Item2)
    : IComparable<Int32Int32>
{
    public override string ToString() => $"{Item1}:{Item2}";

    public int CompareTo(Int32Int32 other)
    {
        if (Item1 < other.Item1) return -1;
        if (Item1 > other.Item1) return 1;
        if (Item2 < other.Item2) return -1;
        if (Item2 > other.Item2) return 1;
        return 0;
    }

    public static bool operator <(Int32Int32 x, Int32Int32 y) => x.CompareTo(y) < 0;
    public static bool operator >(Int32Int32 x, Int32Int32 y) => x.CompareTo(y) > 0;
    public static bool operator <=(Int32Int32 x, Int32Int32 y) => x.CompareTo(y) <= 0;
    public static bool operator >=(Int32Int32 x, Int32Int32 y) => x.CompareTo(y) >= 0;
}
/// <summary>
/// A tuple where Item1(S2CellId) and Item2(Int32) are the identifier and has 2 Items in total
/// </summary>
public readonly record struct S2CellIdInt32(S2CellId Item1, Int32 Item2)
    : IComparable<S2CellIdInt32>
{
    public override string ToString() => $"{Item1}:{Item2}";

    public int CompareTo(S2CellIdInt32 other)
    {
        if (Item1.Id < other.Item1.Id) return -1;
        if (Item1.Id > other.Item1.Id) return 1;
        if (Item2 < other.Item2) return -1;
        if (Item2 > other.Item2) return 1;
        return 0;
    }

    public static bool operator <(S2CellIdInt32 x, S2CellIdInt32 y) => x.CompareTo(y) < 0;
    public static bool operator >(S2CellIdInt32 x, S2CellIdInt32 y) => x.CompareTo(y) > 0;
    public static bool operator <=(S2CellIdInt32 x, S2CellIdInt32 y) => x.CompareTo(y) <= 0;
    public static bool operator >=(S2CellIdInt32 x, S2CellIdInt32 y) => x.CompareTo(y) >= 0;
}
/// <summary>
/// A tuple where Item1, Item2 and Item3 are the identifier and has 3 Items in total
/// </summary>
public readonly record struct Key3<T, U, V>(T Item1, U Item2, V Item3)
    : IComparable<Key3<T, U, V>>
     where T : IEquatable<T>, IComparable<T>
     where U : IEquatable<U>, IComparable<U>
     where V : IEquatable<V>, IComparable<V>
{
    public override string ToString() => $"{Item1}:{Item2}:{Item3}";

    public int CompareTo(Key3<T, U, V> other)
    {
        var c = Item1.CompareTo(other.Item1);
        if (c != 0) return c;

        c = Item2.CompareTo(other.Item2);
        if (c != 0) return c;

        return Item3.CompareTo(other.Item3);
    }
    public static bool operator <(Key3<T, U, V> x, Key3<T, U, V> y) => x.CompareTo(y) < 0;
    public static bool operator >(Key3<T, U, V> x, Key3<T, U, V> y) => x.CompareTo(y) > 0;
    public static bool operator <=(Key3<T, U, V> x, Key3<T, U, V> y) => x.CompareTo(y) <= 0;
    public static bool operator >=(Key3<T, U, V> x, Key3<T, U, V> y) => x.CompareTo(y) >= 0;
}
