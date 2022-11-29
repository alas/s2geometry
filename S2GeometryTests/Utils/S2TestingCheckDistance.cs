namespace S2Geometry;

public static class S2TestingCheckDistance<Id, Distance> where Distance : IEquatable<Distance>, IComparable<Distance>, IDistance<Distance> where Id : IEquatable<Id>
{
    public delegate void WriteLineDelegate(string s);

    // Compare two sets of "closest" items, where "expected" is computed via brute
    // force (i.e., considering every possible candidate) and "actual" is computed
    // using a spatial data structure.  Here "max_size" is a bound on the maximum
    // number of items, "max_distance" is a limit on the distance to any item, and
    // "max_error" is the maximum error allowed when selecting which items are
    // closest (see S2ClosestEdgeQuery.Options.max_error).
    public static bool CheckDistanceResults(List<(Distance, Id)> expected, List<(Distance, Id)> actual,
                    int max_size, Distance max_distance, S1ChordAngle max_error, WriteLineDelegate writeline)
    {
        return CheckResultSet(actual, expected, max_size, max_distance, max_error, kMaxPruningError, "Missing", writeline) & /*not &&*/
               CheckResultSet(expected, actual, max_size, max_distance, max_error, S1ChordAngle.Zero, "Extra", writeline);
    }
    // This is a conservative bound on the error in computing the distance from
    // the target geometry to an S2Cell.  Such errors can cause candidates to be
    // pruned from the result set even though they may be slightly closer.
    private static readonly S1ChordAngle kMaxPruningError = S1ChordAngle.FromRadians(S2.DoubleError);

    // Check that result set "x" contains all the expected results from "y", and
    // does not include any duplicate results.
    private static bool CheckResultSet(List<(Distance, Id)> x, List<(Distance, Id)> y,
                    int max_size, Distance max_distance, S1ChordAngle max_error, S1ChordAngle max_pruning_error,
                    string label, WriteLineDelegate writeline)
    {
        // Results should be sorted by distance, but not necessarily then by Id.
        Assert.True(x.IsSorted(((Distance, Id) x, (Distance, Id) y) =>
        {
            return x.Item1.CompareTo(y.Item1);
        }));

        // Result set X should contain all the items from Y whose distance is less
        // than "limit" computed below.
        Distance limit = Distance.Zero;
        if (x.Count < max_size)
        {
            // Result set X was not limited by "max_size", so it should contain all
            // the items up to "max_distance", except that a few items right near the
            // distance limit may be missed because the distance measurements used for
            // pruning S2Cells are not conservative.
            if (Equals(max_distance, Distance.Infinity))
            {
                limit = max_distance;
            }
            else
            {
                limit = max_distance - max_pruning_error;
            }
        }
        else if (x.Any())
        {
            // Result set X contains only the closest "max_size" items, to within a
            // tolerance of "max_error + max_pruning_error".
            limit = (x.Last().Item1 - max_error) - max_pruning_error;
        }

        bool result = true;
        foreach (var yp in y)
        {
            // Note that this test also catches duplicate values.
            int count = x.Count(((Distance, Id) xp) =>
            {
                return Equals(xp.Item2, yp.Item2);
            });
            if (yp.Item1 < limit && count != 1)
            {
                result = false;
                writeline($"{(count > 1 ? "Duplicate" : label)} distance = {yp.Item1}, id = {yp.Item2}");
            }
        }

        return result;
    }
}


