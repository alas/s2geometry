// S2PointIndex maintains an index of points sorted by leaf S2CellId.  Each
// point can optionally store auxiliary data such as an integer or pointer.
// This can be used to map results back to client data structures.
//
// The class supports adding or removing points dynamically, and provides a
// seekable iterator interface for navigating the index.
//
// You can use this class in conjunction with S2ClosestPointQuery to find the
// closest index points to a given query point.  For example:
//
// void Test(S2Point[] index_points,
//           S2Point[] target_points) {
//   // The template argument allows auxiliary data to be attached to each
//   // point (in this case, the array index).
//   S2PointIndex<int> index;
//   for (int i = 0; i < index_points.size(); ++i) {
//     index.Add(index_points[i], i);
//   }
//   S2ClosestPointQuery<int> query(out index);
//   query.Options().set_max_results(5);
//   for (S2Point target_point : target_points) {
//     S2ClosestPointQueryPointTarget target(target_point);
//     foreach (var& result in query.FindClosestPoints(out target)) {
//       // The Result class contains the following methods:
//       //   distance() is the distance to the target.
//       //   point() is the indexed point.
//       //   data() is the auxiliary data.
//       DoSomething(target_point, result);
//     }
//   }
// }
//
// The Data argument defaults to an empty class, which uses no additional
// space beyond the S2Point itself.  In this case the Data argument is
// required.  For example:
//
//   S2PointIndex<> index;
//   index.Add(point);
//
// Points can be added or removed from the index at any time by calling Add()
// or Remove().  However when the index is modified, you must call Init() on
// each iterator before using it again (or simply create a new iterator).
//
//   index.Add(new_point, 123456);
//   it.Init(out index);
//   it.Seek(target.RangeMin);
//
// You can also access the index directly using the iterator interface.  For
// example, here is how to iterate through all the points in a given S2CellId
// "target_id":
//
//   S2PointIndex<int>.Iterator it(out index);
//   it.Seek(target_id.RangeMin);
//   for (; !it.done() && it.id() <= target_id.RangeMax; it.Next()) {
//     DoSomething(it.id(), it.point(), it.data());
//   }
//
// TODO(ericv): Consider adding an S2PointIndexRegion class, which could be
// used to efficiently compute coverings of a collection of S2Points.
//
// REQUIRES: "Data" has default and copy constructors.
// REQUIRES: "Data" has operator== and operator<.

namespace S2Geometry;

using System.Collections;
using System.Runtime.InteropServices;

public class S2PointIndex<Data> : IEnumerable<TreeNode<Data>>
{
    // Default constructor.
    public S2PointIndex() { }

    // Returns the number of points in the index.
    public int NumPoints()
    {
        return map_.Count;
    }

    // Convenience function for the case when Data is an empty class.
    public void Add(S2Point point)
    {
        //Assert.True(typeof(Data) == void); // Data must be empty
        Add(point, default);
    }
    // Adds the given point to the index.  Invalidates all iterators.
    public void Add(S2Point point, Data data)
    {
        var id = new S2CellId(point);
        var tup = new TreeNode<Data>(id, point, data);
        map_.AddSorted(tup);
    }

    // Convenience function for the case when Data is an empty class.
    public void Remove(S2Point point)
    {
        //Assert.True(typeof(Data) == void); // Data must be empty
        Remove(point, default);
    }
    // Removes the given point from the index.  Both the "point" and "data"
    // fields must match the point to be removed.  Returns false if the given
    // point was not present.  Invalidates all iterators.
    public bool Remove(S2Point point, Data data)
    {
        var id = new S2CellId(point);
        var limit = map_.Count;
        var point_data = (point, data);
        var init = map_.GetLowerBound(new(id, point, data));
        for (var i = init; i < limit && map_[i].Id == id; i++)
        {
            var cur = (map_[i].Point, map_[i].Data);
            if (Equals(cur, point_data))
            {
                map_.RemoveAt(i);
                return true;
            }
        }
        return false;
    }

    // Resets the index to its original empty state.  Invalidates all iterators.
    public void Clear()
    {
        map_.Clear();
    }

    // Returns the number of bytes currently occupied by the index.
    public int SpaceUsed()
    {
        var sizeoftreenode = SizeHelper.SizeOf(typeof(List<TreeNode<Data>>));
        return SizeHelper.SizeOf(this) - sizeoftreenode + SizeHelper.SizeOf(map_);
    }

    public int Seek(S2CellId target)
    {
        return map_.GetLowerBound(new(target, default, default));
    }

    public TreeNode<Data>? Get(int pos)
    {
        return map_.Count > pos ? map_[pos] : null;
    }

    public IEnumerator<TreeNode<Data>> GetEnumerator()
    {
        return map_.GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    private readonly List<TreeNode<Data>> map_ = []; // gtl.btree_multimap<S2CellId, PointData>;
}
