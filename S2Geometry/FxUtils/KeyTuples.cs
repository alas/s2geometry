using System;
using System.Diagnostics.CodeAnalysis;

namespace S2Geometry
{
    /// <summary>
    /// A tuple where Item1 is the identifier and has 2 Items in total
    /// </summary>
    public readonly struct KeyData<T, U>
        : IEquatable<KeyData<T, U>>, IComparable<KeyData<T, U>>
         where T : IEquatable<T>, IComparable<T>
    {
        public readonly T Item1;
        public readonly U Item2;
        public KeyData(T item1, U item2) { Item1 = item1; Item2 = item2; }
        public override string ToString()
        {
            return $"{Item1}:{Item2}";
        }
        public override int GetHashCode()
        {
            return Item1.GetHashCode();
        }
        public int CompareTo([AllowNull] KeyData<T, U> other)
        {
            return Item1.CompareTo(other.Item1);
        }
        public bool Equals(KeyData<T, U> other)
        {
            return Equals(Item1, other.Item1);
        }
        public override bool Equals(object obj)
        {
            return obj is KeyData<T, U> id && Equals(Item1, id.Item1);
        }
        public static bool operator ==(KeyData<T, U> x, KeyData<T, U> y)
        {
            return Equals(x.Item1, y.Item1);
        }
        public static bool operator !=(KeyData<T, U> x, KeyData<T, U> y)
        {
            return !Equals(x.Item1, y.Item1);
        }
        public static bool operator <(KeyData<T, U> x, KeyData<T, U> y)
        {
            return x.Item1.CompareTo(y.Item1) < 0;
        }
        public static bool operator >(KeyData<T, U> x, KeyData<T, U> y)
        {
            return x.Item1.CompareTo(y.Item1) > 0;
        }
        public static bool operator <=(KeyData<T, U> x, KeyData<T, U> y)
        {
            return x.Item1.CompareTo(y.Item1) <= 0;
        }
        public static bool operator >=(KeyData<T, U> x, KeyData<T, U> y)
        {
            return x.Item1.CompareTo(y.Item1) >= 0;
        }
    }
    /// <summary>
    /// A tuple where Item1 is the identifier and has 3 Items in total
    /// </summary>
    public readonly struct KeyData2<T, U, V> : IEquatable<KeyData2<T, U, V>>, IComparable<KeyData2<T, U, V>>
         where T : IEquatable<T>, IComparable<T>
    {
        public readonly T Item1;
        public readonly U Item2;
        public readonly V Item3;
        public KeyData2(T item1, U item2, V item3) { Item1 = item1; Item2 = item2; Item3 = item3; }

        public override string ToString()
        {
            return $"{Item1}:{Item2}:{Item3}";
        }
        public override int GetHashCode()
        {
            return Item1.GetHashCode();
        }
        public int CompareTo([AllowNull] KeyData2<T, U, V> other)
        {
            return Item1.CompareTo(other.Item1);
        }
        public bool Equals(KeyData2<T, U, V> other)
        {
            return Equals(Item1, other.Item1);
        }
        public override bool Equals(object obj)
        {
            return obj is KeyData2<T, U, V> id && Equals(Item1, id.Item1);
        }
        public static bool operator ==(KeyData2<T, U, V> x, KeyData2<T, U, V> y)
        {
            return Equals(x.Item1, y.Item1);
        }
        public static bool operator !=(KeyData2<T, U, V> x, KeyData2<T, U, V> y)
        {
            return !Equals(x.Item1, y.Item1);
        }
        public static bool operator <(KeyData2<T, U, V> x, KeyData2<T, U, V> y)
        {
            return x.Item1.CompareTo(y.Item1) < 0;
        }
        public static bool operator >(KeyData2<T, U, V> x, KeyData2<T, U, V> y)
        {
            return x.Item1.CompareTo(y.Item1) > 0;
        }
        public static bool operator <=(KeyData2<T, U, V> x, KeyData2<T, U, V> y)
        {
            return x.Item1.CompareTo(y.Item1) <= 0;
        }
        public static bool operator >=(KeyData2<T, U, V> x, KeyData2<T, U, V> y)
        {
            return x.Item1.CompareTo(y.Item1) >= 0;
        }
    }
    /// <summary>
    /// A tuple where Item1 is the identifier and has 2 Items in total, Item1 is reverse compared
    /// </summary>
    public readonly struct ReverseKeyData<T, U>
        : IEquatable<ReverseKeyData<T, U>>, IComparable<ReverseKeyData<T, U>>
         where T : IEquatable<T>, IComparable<T>
    {
        public readonly T Item1;
        public readonly U Item2;
        public ReverseKeyData(T item1, U item2) { Item1 = item1; Item2 = item2; }
        public override string ToString()
        {
            return $"{Item1}:{Item2}";
        }
        public override int GetHashCode()
        {
            return Item1.GetHashCode();
        }
        public int CompareTo([AllowNull] ReverseKeyData<T, U> other)
        {
            return other.Item1.CompareTo(Item1);
        }
        public bool Equals(ReverseKeyData<T, U> other)
        {
            return Equals(Item1, other.Item1);
        }
        public override bool Equals(object obj)
        {
            return obj is ReverseKeyData<T, U> id && Equals(Item1, id.Item1);
        }
        public static bool operator ==(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y)
        {
            return Equals(x.Item1, y.Item1);
        }
        public static bool operator !=(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y)
        {
            return !Equals(x.Item1, y.Item1);
        }
        public static bool operator <(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y)
        {
            return y.Item1.CompareTo(x.Item1) < 0;
        }
        public static bool operator >(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y)
        {
            return y.Item1.CompareTo(x.Item1) > 0;
        }
        public static bool operator <=(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y)
        {
            return y.Item1.CompareTo(x.Item1) <= 0;
        }
        public static bool operator >=(ReverseKeyData<T, U> x, ReverseKeyData<T, U> y)
        {
            return y.Item1.CompareTo(x.Item1) >= 0;
        }
    }
    /// <summary>
    /// A tuple where Item1 is the identifier and has 3 Items in total, Item1 is reverse compared
    /// </summary>
    public readonly struct ReverseKeyData2<T, U, V> : IEquatable<ReverseKeyData2<T, U, V>>, IComparable<ReverseKeyData2<T, U, V>>
         where T : IEquatable<T>, IComparable<T>
    {
        public readonly T Item1;
        public readonly U Item2;
        public readonly V Item3;
        public ReverseKeyData2(T item1, U item2, V item3) { Item1 = item1; Item2 = item2; Item3 = item3; }

        public override string ToString()
        {
            return $"{Item1}:{Item2}:{Item3}";
        }
        public override int GetHashCode()
        {
            return Item1.GetHashCode();
        }
        public int CompareTo([AllowNull] ReverseKeyData2<T, U, V> other)
        {
            return other.Item1.CompareTo(Item1);
        }
        public bool Equals(ReverseKeyData2<T, U, V> other)
        {
            return Equals(Item1, other.Item1);
        }
        public override bool Equals(object obj)
        {
            return obj is ReverseKeyData2<T, U, V> id && Equals(Item1, id.Item1);
        }
        public static bool operator ==(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y)
        {
            return Equals(x.Item1, y.Item1);
        }
        public static bool operator !=(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y)
        {
            return !Equals(x.Item1, y.Item1);
        }
        public static bool operator <(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y)
        {
            return y.Item1.CompareTo(x.Item1) < 0;
        }
        public static bool operator >(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y)
        {
            return y.Item1.CompareTo(x.Item1) > 0;
        }
        public static bool operator <=(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y)
        {
            return y.Item1.CompareTo(x.Item1) <= 0;
        }
        public static bool operator >=(ReverseKeyData2<T, U, V> x, ReverseKeyData2<T, U, V> y)
        {
            return y.Item1.CompareTo(x.Item1) >= 0;
        }
    }
    /// <summary>
    /// A tuple where Item1 and Item2 are the identifier and has 2 Items in total
    /// </summary>
    public readonly struct KeyKey<T, U>
        : IEquatable<KeyKey<T, U>>, IComparable<KeyKey<T, U>>
         where T : IEquatable<T>, IComparable<T>
         where U : IEquatable<U>, IComparable<U>
    {
        public readonly T Item1;
        public readonly U Item2;
        public KeyKey(T item1, U item2) { Item1 = item1; Item2 = item2; }
        public override string ToString()
        {
            return $"{Item1}:{Item2}";
        }
        public override int GetHashCode()
        {
            return HashCode.Combine(Item1, Item2);
        }
        public int CompareTo([AllowNull] KeyKey<T, U> other)
        {
            if (Item1.CompareTo(other.Item1) != 0)
                return Item1.CompareTo(other.Item1);

            return Item2.CompareTo(other.Item2);
        }
        public bool Equals(KeyKey<T, U> other)
        {
            return Equals(Item1, other.Item1) && Equals(Item2, other.Item2);
        }
        public override bool Equals(object obj)
        {
            return obj is KeyKey<T, U> id && Equals(id);
        }
        public static bool operator ==(KeyKey<T, U> x, KeyKey<T, U> y)
        {
            return Equals(x, y);
        }
        public static bool operator !=(KeyKey<T, U> x, KeyKey<T, U> y)
        {
            return !Equals(x, y);
        }
        public static bool operator <(KeyKey<T, U> x, KeyKey<T, U> y)
        {
            return x.CompareTo(y) < 0;
        }
        public static bool operator >(KeyKey<T, U> x, KeyKey<T, U> y)
        {
            return x.CompareTo(y) > 0;
        }
        public static bool operator <=(KeyKey<T, U> x, KeyKey<T, U> y)
        {
            return x.CompareTo(y) <= 0;
        }
        public static bool operator >=(KeyKey<T, U> x, KeyKey<T, U> y)
        {
            return x.CompareTo(y) >= 0;
        }
    }
}
