namespace S2Geometry;

public class R1IntervalTests
{
    [Fact]
    internal void Test_R1Interval_TestBasic() {
        // Constructors and accessors.
        R1Interval unit = new(0, 1);
        R1Interval negunit = new(-1, 0);
        Assert.Equal(0, unit.Lo);
        Assert.Equal(1, unit.Hi);
        Assert.Equal(-1, negunit[0]);
        Assert.Equal(0, negunit[1]);
        R1Interval ten = new(0, 10);
        Assert.Equal(new R1Interval(0, 10), ten);
        ten = new R1Interval(-10, ten.Hi);
        Assert.Equal(new R1Interval(-10, 10), ten);
        ten = new R1Interval(ten.Lo, 0);
        Assert.True(new R1Interval(-10, 0) == ten);
        ten = new R1Interval(0, 10);
        Assert.Equal(new R1Interval(0, 10), ten);

        // IsEmpty
        R1Interval half = new(0.5, 0.5);
        Assert.False(unit.IsEmpty());
        Assert.False(half.IsEmpty());
        R1Interval empty = R1Interval.Empty;
        Assert.True(empty.IsEmpty());

        // == and !=
        Assert.True(empty == R1Interval.Empty);
        Assert.True(unit == new R1Interval(0, 1));
        Assert.True(unit != empty);
        Assert.True(new R1Interval(1, 2) != new R1Interval(1, 3));

        // Check that the default R1Interval is identical to Empty().
        R1Interval default_empty = R1Interval.Empty;
        Assert.True(default_empty.IsEmpty());
        Assert.Equal(empty.Lo, default_empty.Lo);
        Assert.Equal(empty.Hi, default_empty.Hi);

        // GetCenter(), GetLength()
        Assert.Equal(0.5, unit.GetCenter());
        Assert.Equal(0.5, half.GetCenter());
        Assert.Equal(1.0, negunit.GetLength());
        Assert.Equal(0, half.GetLength());
        Assert.True(empty.GetLength() < 0);

        // Contains(double), InteriorContains(double)
        Assert.True(unit.Contains(0.5));
        Assert.True(unit.InteriorContains(0.5));
        Assert.True(unit.Contains(0));
        Assert.False(unit.InteriorContains(0));
        Assert.True(unit.Contains(1));
        Assert.False(unit.InteriorContains(1));

        // Contains(R1Interval), InteriorContains(R1Interval)
        // Intersects(R1Interval), InteriorIntersects(R1Interval)
        TestIntervalOps(empty, empty, "TTFF");
        TestIntervalOps(empty, unit, "FFFF");
        TestIntervalOps(unit, half, "TTTT");
        TestIntervalOps(unit, unit, "TFTT");
        TestIntervalOps(unit, empty, "TTFF");
        TestIntervalOps(unit, negunit, "FFTF");
        TestIntervalOps(unit, new R1Interval(0, 0.5), "TFTT");
        TestIntervalOps(half, new R1Interval(0, 0.5), "FFTF");

        // AddPoint()
        R1Interval r = empty;
        r = R1Interval.AddPoint(r, 5);
        Assert.Equal(5, r.Lo);
        Assert.Equal(5, r.Hi);
        r = R1Interval.AddPoint(r, -1);
        Assert.Equal(-1, r.Lo);
        Assert.Equal(5, r.Hi);
        r = R1Interval.AddPoint(r, 0);
        Assert.Equal(-1, r.Lo);
        Assert.Equal(5, r.Hi);

        // Project()
        Assert.Equal(0.3, new R1Interval(0.1, 0.4).Project(0.3));
        Assert.Equal(0.1, new R1Interval(0.1, 0.4).Project(-7.0));
        Assert.Equal(0.4, new R1Interval(0.1, 0.4).Project(0.6));

        // FromPointPair()
        Assert.Equal(new R1Interval(4, 4), R1Interval.FromPointPair(4, 4));
        Assert.Equal(new R1Interval(-2, -1), R1Interval.FromPointPair(-1, -2));
        Assert.Equal(new R1Interval(-5, 3), R1Interval.FromPointPair(-5, 3));

        // Expanded()
        Assert.Equal(empty, empty.Expanded(0.45));
        Assert.Equal(new R1Interval(-0.5, 1.5), unit.Expanded(0.5));
        Assert.Equal(new R1Interval(0.5, 0.5), unit.Expanded(-0.5));
        Assert.True(unit.Expanded(-0.51).IsEmpty());
        Assert.True(unit.Expanded(-0.51).Expanded(0.51).IsEmpty());

        // Union(), Intersection()
        Assert.Equal(new R1Interval(99, 100), new R1Interval(99, 100).Union(empty));
        Assert.Equal(new R1Interval(99, 100), empty.Union(new R1Interval(99, 100)));
        Assert.True(new R1Interval(5, 3).Union(new R1Interval(0, -2)).IsEmpty());
        Assert.True(new R1Interval(0, -2).Union(new R1Interval(5, 3)).IsEmpty());
        Assert.Equal(unit, unit.Union(unit));
        Assert.Equal(new R1Interval(-1, 1), unit.Union(negunit));
        Assert.Equal(new R1Interval(-1, 1), negunit.Union(unit));
        Assert.Equal(unit, half.Union(unit));
        Assert.Equal(half, unit.Intersection(half));
        Assert.Equal(new R1Interval(0, 0), unit.Intersection(negunit));
        Assert.True(negunit.Intersection(half).IsEmpty());
        Assert.True(unit.Intersection(empty).IsEmpty());
        Assert.True(empty.Intersection(unit).IsEmpty());
    }

    [Fact]
    internal void Test_R1Interval_ApproxEquals()
    {
        // Choose two values kLo and kHi such that it's okay to shift an endpoint by
        // kLo (i.e., the resulting interval is equivalent) but not by kHi.
        const double kLo = 4 * S2.DoubleEpsilon;  // < max_error default
        const double kHi = 6 * S2.DoubleEpsilon;  // > max_error default

        // Empty intervals.
        R1Interval empty = R1Interval.Empty;
        Assert.True(empty.ApproxEquals(empty));
        Assert.True(new R1Interval(0, 0).ApproxEquals(empty));
        Assert.True(empty.ApproxEquals(new R1Interval(0, 0)));
        Assert.True(new R1Interval(1, 1).ApproxEquals(empty));
        Assert.True(empty.ApproxEquals(new R1Interval(1, 1)));
        Assert.False(empty.ApproxEquals(new R1Interval(0, 1)));
        Assert.True(empty.ApproxEquals(new R1Interval(1, 1 + 2 * kLo)));
        Assert.False(empty.ApproxEquals(new R1Interval(1, 1 + 2 * kHi)));

        // Singleton intervals.
        Assert.True(new R1Interval(1, 1).ApproxEquals(new R1Interval(1, 1)));
        Assert.True(new R1Interval(1, 1).ApproxEquals(new R1Interval(1 - kLo, 1 - kLo)));
        Assert.True(new R1Interval(1, 1).ApproxEquals(new R1Interval(1 + kLo, 1 + kLo)));
        Assert.False(new R1Interval(1, 1).ApproxEquals(new R1Interval(1 - kHi, 1)));
        Assert.False(new R1Interval(1, 1).ApproxEquals(new R1Interval(1, 1 + kHi)));
        Assert.True(new R1Interval(1, 1).ApproxEquals(new R1Interval(1 - kLo, 1 + kLo)));
        Assert.False(new R1Interval(0, 0).ApproxEquals(new R1Interval(1, 1)));

        // Other intervals.
        Assert.True(new R1Interval(1 - kLo, 2 + kLo).ApproxEquals(new R1Interval(1, 2)));
        Assert.True(new R1Interval(1 + kLo, 2 - kLo).ApproxEquals(new R1Interval(1, 2)));
        Assert.False(new R1Interval(1 - kHi, 2 + kLo).ApproxEquals(new R1Interval(1, 2)));
        Assert.False(new R1Interval(1 + kHi, 2 - kLo).ApproxEquals(new R1Interval(1, 2)));
        Assert.False(new R1Interval(1 - kLo, 2 + kHi).ApproxEquals(new R1Interval(1, 2)));
        Assert.False(new R1Interval(1 + kLo, 2 - kHi).ApproxEquals(new R1Interval(1, 2)));
    }

    private static void TestIntervalOps(R1Interval x, R1Interval y, string expected)
    {
        // Test all of the interval operations on the given pair of intervals.
        // "expected" is a sequence of "T" and "F" characters corresponding to
        // the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.

        Assert.Equal(expected[0] == 'T', x.Contains(y));
        Assert.Equal(expected[1] == 'T', x.InteriorContains(y));
        Assert.Equal(expected[2] == 'T', x.Intersects(y));
        Assert.Equal(expected[3] == 'T', x.InteriorIntersects(y));

        Assert.Equal(x.Contains(y), x.Union(y) == x);
        Assert.Equal(x.Intersects(y), !x.Intersection(y).IsEmpty());

        R1Interval z = x;
        z = R1Interval.AddInterval(z, y);
        Assert.Equal(x.Union(y), z);
    }
} 
