using System;
using System.Collections.Generic;
using System.Linq;

namespace S2Geometry
{
    public static class LinqUtils
    {
        /// <summary>
        /// Discard the last n items
        /// </summary>
        /// <remarks>from: https://stackoverflow.com/questions/1779129/how-to-take-all-but-the-last-element-in-a-sequence-using-linq</remarks>
        public static IEnumerable<T> SkipLastN<T>(this IEnumerable<T> source, int n)
        {
            var it = source.GetEnumerator();
            bool hasRemainingItems;
            var cache = new Queue<T>(n + 1);

            do
            {
                if (hasRemainingItems = it.MoveNext())
                {
                    cache.Enqueue(it.Current);
                    if (cache.Count > n)
                        yield return cache.Dequeue();
                }
            } while (hasRemainingItems);
        }

        public static T[] Fill<T>(this T[] arr, T value)
        {
            Array.Fill(arr, value);
            return arr;
        }

        public static T[] Fill<T>(this T[] arr, Func<T> getDefault)
        {
            for (int i = 0; i < arr.Length; i++)
            {
                arr[i] = getDefault();
            }
            return arr;
        }

        public static IList<T> Fill<T>(this IList<T> arr, T value, int count)
        {
            for (int i = 0; i < count; i++)
            {
                arr.Add(value);
            }
            return arr;
        }

        public static IList<T> Fill<T>(this IList<T> arr, Func<T> getNewValue, int count)
        {
            for (int i = 0; i < count; i++)
            {
                arr.Add(getNewValue());
            }
            return arr;
        }

        public static void Iota(this List<int> arr, int value, int count)
        {
            if (count == 0) return;

            arr.Capacity = count;

            arr.Add(value);
            for (int i = 1; i < count; i++)
            {
                arr.Add(++value);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <remarks>from: https://stackoverflow.com/questions/16323386/fast-efficient-way-to-get-index-of-minimum-value-in-listt/16323480</remarks>
        public static int IndexOfMin<T>(this IList<T> self) where T : IComparable<T>
        {
            if (self == null)
            {
                throw new ArgumentNullException(nameof(self));
            }

            if (self.Count == 0)
            {
                throw new ArgumentException("List is empty.", nameof(self));
            }

            T min = self[0];
            int minIndex = 0;

            for (int i = 1; i < self.Count; i++)
            {
                if (self[i].CompareTo(min) < 0)
                {
                    min = self[i];
                    minIndex = i;
                }
            }

            return minIndex;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <remarks>from: https://stackoverflow.com/questions/8094867/good-gethashcode-override-for-list-of-foo-objects-respecting-the-order</remarks>
        public static int GetSequenceHashCode<T>(this IEnumerable<T> list)
        {
            if (list == null) return 0;
            const int seedValue = 0x2D2816FE;
            const int primeNumber = 397;
            return list.Aggregate(seedValue, (current, item) => (current * primeNumber)
                + (Equals(item, default(T)) ? 0 : item.GetHashCode()));
        }

        public static T[] DeepCopy<T>(this T[] arr) where T : ICloneable
        {
            var newArr = new T[arr.Length];
            for (var i = 0; i < arr.Length; i++)
            {
                newArr[i] = (T)arr[i].Clone();
            }
            return newArr;
        }

        public static List<T> DeepCopy<T>(this List<T> arr) where T : ICloneable
        {
            return arr.ConvertAll(t => (T)t.Clone());
        }
        
        #region lower_bound

        /// <summary>
        /// replacement for C++ lower_bound using indexes instead of iterators
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="arr"></param>
        /// <param name="first">index of the first element in the range to examine</param>
        /// <param name="last">index of the last element in the range to examine</param>
        /// <param name="value"></param>
        /// <returns>index of the first element that is not less than value, or last if no such element is found.</returns>
        public static int GetLowerBound<T>(this List<T> arr, T value, int first, int last, IComparer<T> comp)
        {
            var i = arr.BinarySearch(first, last - first, value, comp);
            if (i < 0) return ~i;

            while (i > 0 && comp.Compare(arr[i], value) == 0) i--;

            return i + 1;
        }

        public static int GetLowerBound<T>(this List<T> arr, T value, IComparer<T> comp)
        {
            var i = comp == null ? arr.BinarySearch(value) : arr.BinarySearch(0, arr.Count, value, comp);
            if (i < 0) return ~i;

            while (i > 0 && comp.Compare(arr[i], value) == 0) i--;

            return i + 1;
        }

        public static int GetLowerBound<T>(this List<T> arr, T value, int first, int last) where T : IComparable<T>
        {
            var i = arr.BinarySearch(first, last - first, value, null);
            if (i < 0) return ~i;

            while (i > 0 && value.CompareTo(arr[i]) == 0) i--;

            return i + 1;
        }

        public static int GetLowerBound<T>(this List<T> arr, T value) where T : IComparable<T>
        {
            var i = arr.BinarySearch(value);
            if (i < 0) return ~i;

            while (i > 0 && value.CompareTo(arr[i]) == 0) i--;

            return i + 1;
        }

        public static int GetLowerBound<T>(this T[] arr, T value, int first, int last) where T : IComparable<T>
        {
            var i = Array.BinarySearch(arr, first, last - first, value, null);
            if (i < 0) return ~i;

            while (i > 0 && arr[i].CompareTo(value) == 0) i--;

            return i + 1;
        }

        public static int GetLowerBound<T>(this T[] arr, T value) where T : IComparable<T>
        {
            var i = Array.BinarySearch(arr, value);
            if (i < 0) return ~i;

            while (i > 0 && arr[i].CompareTo(value) == 0) i--;

            return i + 1;
        }

        #endregion

        #region upper_bound

        public static int GetUpperBound<T>(this List<T> arr, T value, int first, int last) where T : IComparable<T>
        {
            var i = arr.BinarySearch(first, last - first, value, null);
            if (i < 0) return ~i;

            while (i < arr.Count && arr[i].CompareTo(value) == 0) i++;

            return i;
        }

        public static int GetUpperBound<T>(this List<T> arr, T value) where T : IComparable<T>
        {
            var i = arr.BinarySearch(value);
            if (i < 0) return ~i;

            while (i < arr.Count && arr[i].CompareTo(value) == 0) i++;

            return i;
        }

        public static int GetUpperBound<T>(this List<T> arr, T value, IComparer<T> comp)
        {
            var i = arr.BinarySearch(value, comp);
            if (i < 0) return ~i;

            while (i < arr.Count && comp.Compare(arr[i], value) == 0) i++;

            return i;
        }

        public static int GetUpperBound<T>(this T[] arr, T value) where T : IComparable<T>
        {
            var i = Array.BinarySearch(arr, value);
            if (i < 0) return ~i;

            while (i < arr.Length && arr[i].CompareTo(value) == 0) i++;

            return i;
        }

        #endregion

        /// <summary>
        /// Determines if string array is sorted from A -> Z
        /// </summary>
        /// <remarks>from: https://www.dotnetperls.com/issorted</remarks>
        public static bool IsSorted<T>(this IList<T> arr) where T : IComparable<T>
        {
            for (int i = 1; i < arr.Count; i++)
            {
                if (arr[i - 1].CompareTo(arr[i]) > 0) // If previous is bigger, return false
                {
                    return false;
                }
            }
            return true;
        }

        public static bool IsSorted<T>(this IList<T> arr, Func<T, T, int> comp)
        {
            for (int i = 1; i < arr.Count; i++)
            {
                if (comp(arr[i - 1], arr[i]) > 0) // If previous is bigger, return false
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Determines if string array is sorted from Z -> A
        /// </summary>
        /// <remarks>from: https://www.dotnetperls.com/issorted</remarks>
        public static bool IsSortedDescending<T>(this IList<T> arr) where T : IComparable<T>
        {
            for (int i = arr.Count - 2; i >= 0; i--)
            {
                if (arr[i].CompareTo(arr[i + 1]) < 0) // If previous is smaller, return false
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Inserts item into list in order, List is alreay sorted
        /// </summary>
        /// <remarks>from: https://stackoverflow.com/questions/12172162/how-to-insert-item-into-list-in-order</remarks>
        public static void AddSorted<T>(this List<T> @this, T item) where T : IComparable<T>
        {
            if (@this.Count == 0)
            {
                @this.Add(item);
                return;
            }
            if (@this[^1].CompareTo(item) <= 0)
            {
                @this.Add(item);
                return;
            }
            if (@this[0].CompareTo(item) >= 0)
            {
                @this.Insert(0, item);
                return;
            }
            int index = @this.BinarySearch(item);
            if (index < 0) index = ~index;
            @this.Insert(index, item);
        }
        public static bool AddSortedUnique<T>(this List<T> @this, T item) where T : IComparable<T>
        {
            if (@this.Count == 0)
            {
                @this.Add(item);
                return true;
            }
            if (@this[^1].CompareTo(item) < 0)
            {
                @this.Add(item);
                return true;
            }
            if (@this[0].CompareTo(item) > 0)
            {
                @this.Insert(0, item);
                return true;
            }
            int index = @this.BinarySearch(item);
            if (index < 0)
            {
                @this.Insert(~index, item);
                return true;
            }

            return false;
        }

        /// <summary>
        /// Inserts item into list in order, List is alreay sorted
        /// </summary>
        /// <remarks>from: https://stackoverflow.com/questions/12172162/how-to-insert-item-into-list-in-order</remarks>
        public static void AddSorted<T>(this List<T> @this, T item, IComparer<T> comp)
        {
            if (@this.Count == 0)
            {
                @this.Add(item);
                return;
            }
            if (comp.Compare(@this[^1], item) <= 0)
            {
                @this.Add(item);
                return;
            }
            if (comp.Compare(@this[0], item) >= 0)
            {
                @this.Insert(0, item);
                return;
            }
            int index = @this.BinarySearch(0, @this.Count, item, comp);
            if (index < 0) index = ~index;
            @this.Insert(index, item);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <remarks>from: https://stackoverflow.com/questions/12231569/is-there-in-c-sharp-a-method-for-listt-like-resize-in-c-for-vectort</remarks>
        public static void Resize<T>(this List<T> list, int sz, Func<T> getValue)
        {
            int cur = list.Count;
            if (sz < cur)
            {
                list.RemoveRange(sz, cur - sz);
            }
            else if (sz > cur)
            {
                if (sz > list.Capacity)//this bit is purely an optimisation, to avoid multiple automatic capacity changes.
                    list.Capacity = sz;
                for (var i = 0; i < (sz - cur); i++)
                    list.Add(getValue());
            }
        }
        public static void Resize<T>(this List<T> list, int sz, T c = default)
        {
            int cur = list.Count;
            if (sz < cur)
            {
                list.RemoveRange(sz, cur - sz);
            }
            else if (sz > cur)
            {
                if (sz > list.Capacity)//this bit is purely an optimisation, to avoid multiple automatic capacity changes.
                    list.Capacity = sz;
                list.AddRange(Enumerable.Repeat(c, sz - cur));
            }
        }

        /// <summary>
        /// Determines whether a System.Collections.Generic.List\<TObject\> is a subset of the specified collection.
        /// http://stackoverflow.com/questions/332973/linq-check-whether-an-array-is-a-subset-of-another
        /// </summary>
        /// <param name="list"></param>
        /// <param name="listToFind"></param>
        /// <returns></returns>
        public static bool IsSubsetOf<T>(this IEnumerable<T> coll1, IEnumerable<T> coll2)
        {
            bool isSubset = !coll1.Except(coll2).Any();
            return isSubset;
        }

        public static void RotateInPlace<T>(this List<T> lst, int count)
        {
            if (count == 0) return;

            var tmp = lst.Take(count).ToList();
            lst.RemoveRange(0, count);
            lst.AddRange(tmp);
        }

        public static void RotateInPlace<T>(this T[] lst, int count)
        {
            if (count == 0) return;

            var tmp = lst.Take(count).ToArray();
            Array.Copy(lst, count, lst, 0, lst.Length - count);
            Array.Copy(tmp, 0, lst, lst.Length - count, count);
        }

        /// <summary>
        /// from: https://stackoverflow.com/questions/1287567/is-using-random-and-orderby-a-good-shuffle-algorithm
        /// </summary>
        public static IEnumerable<T> Shuffle<T>(this List<T> source, Func<int, int> random)
        {
            for (int i = source.Count - 1; i >= 0; i--)
            {
                // Swap element "i" with a random earlier element it (or itself)
                // ... except we don't really need to swap it fully, as we can
                // return it immediately, and afterwards it's irrelevant.
                int swapIndex = random(i + 1);
                yield return source[swapIndex];
                source[swapIndex] = source[i];
            }
        }

        public static void Swap<T>(List<T> a, List<T> b)
        {
            var tmp = a.ToList();
            a.Clear();
            a.AddRange(b);
            b.Clear();
            b.AddRange(tmp);
        }

        /// <summary>
        /// https://stackoverflow.com/questions/16192906/net-dictionary-get-or-create-new
        /// </summary>
        public static TValue GetOrCreate<TKey, TValue>(this IDictionary<TKey, TValue> dict, TKey key, Func<TValue> createNew)
        {
            if (!dict.TryGetValue(key, out var val))
            {
                val = createNew();
                dict.Add(key, val);
            }

            return val;
        }

        public static T GetRemIndex<T>(this T[] arr, int index)
        {
            return arr[index % arr.Length];
        }
    }
}
