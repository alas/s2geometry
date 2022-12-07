// Most of the R2Rect methods have trivial implementations in terms of the
// R1Interval class, so most of the testing is done in that unit test.

namespace S2Geometry;

public class R2RectTests
{
    [Fact]
    internal void Test_R2Rect_EmptyRectangles()
    {
        // Test basic properties of empty rectangles.
        R2Rect empty = R2Rect.Empty;
        Assert.True(empty.IsValid());
        Assert.True(empty.IsEmpty());
        Assert.Equal(empty, empty);
    }

    [Fact]
    internal void Test_R2Rect_ConstructorsAndAccessors()
    {
        // Check various constructors and accessor methods.
        R2Rect r = new(new R2Point(0.1, 0), new R2Point(0.25, 1));
        Assert.Equal(0.1, r.X.Lo);
        Assert.Equal(0.25, r.X.Hi);
        Assert.Equal(0.0, r.Y.Lo);
        Assert.Equal(1.0, r.Y.Hi);

        Assert.Equal(0.1, r[0][0]);
        Assert.Equal(0.25, r[0][1]);
        Assert.Equal(0.0, r[1][0]);
        Assert.Equal(1.0, r[1][1]);

        Assert.Equal(new R1Interval(0.1, 0.25), r.X);
        Assert.Equal(new R1Interval(0, 1), r.Y);

        Assert.Equal(new R1Interval(0.1, 0.25), r[0]);
        Assert.Equal(new R1Interval(0, 1), r[1]);

        r = new R2Rect(new R1Interval(3, 4), new R1Interval(5, 6));
        Assert.Equal(new R1Interval(3, 4), r[0]);
        Assert.Equal(new R1Interval(5, 6), r[1]);

        Assert.Equal(r, r);
        Assert.NotEqual(r, R2Rect.Empty);

        R2Rect r2 = new();
        Assert.True(r2.IsEmpty());
        Assert.Equal(r2, R2Rect.Empty);
    }

    [Fact]
    internal void Test_R2Rect_FromCenterSize()
    {
        // FromCenterSize()
        Assert.True(R2Rect.FromCenterSize(new R2Point(0.3, 0.5), new R2Point(0.2, 0.4)).
                    ApproxEquals(new R2Rect(new R2Point(0.2, 0.3), new R2Point(0.4, 0.7))));
        Assert.True(R2Rect.FromCenterSize(new R2Point(1, 0.1), new R2Point(0, 2)).
                    ApproxEquals(new R2Rect(new R2Point(1, -0.9), new R2Point(1, 1.1))));
    }

    [Fact]
    internal void Test_R2Rect_FromPoint()
    {
        // FromPoint(), FromPointPair()
        R2Rect d1 = new(new R2Point(0.1, 0), new R2Point(0.25, 1));
        Assert.Equal(new R2Rect(d1.Lo(), d1.Lo()), R2Rect.FromPoint(d1.Lo()));
        Assert.Equal(new R2Rect(new R2Point(0.15, 0.3), new R2Point(0.35, 0.9)),
                  R2Rect.FromPointPair(new R2Point(0.15, 0.9), new R2Point(0.35, 0.3)));
        Assert.Equal(new R2Rect(new R2Point(0.12, 0), new R2Point(0.83, 0.5)),
                  R2Rect.FromPointPair(new R2Point(0.83, 0), new R2Point(0.12, 0.5)));
    }

    [Fact]
    internal void Test_R2Rect_SimplePredicates()
    {
        // GetCenter(), GetVertex(), Contains(R2Point), InteriorContains(R2Point).
        R2Point sw1 = new(0, 0.25);
        R2Point ne1 = new(0.5, 0.75);
        R2Rect r1 = new(sw1, ne1);

        Assert.Equal(new R2Point(0.25, 0.5), r1.GetCenter());
        Assert.Equal(new R2Point(0, 0.25), r1.GetVertex(0));
        Assert.Equal(new R2Point(0.5, 0.25), r1.GetVertex(1));
        Assert.Equal(new R2Point(0.5, 0.75), r1.GetVertex(2));
        Assert.Equal(new R2Point(0, 0.75), r1.GetVertex(3));
        Assert.True(r1.Contains(new R2Point(0.2, 0.4)));
        Assert.False(r1.Contains(new R2Point(0.2, 0.8)));
        Assert.False(r1.Contains(new R2Point(-0.1, 0.4)));
        Assert.False(r1.Contains(new R2Point(0.6, 0.1)));
        Assert.True(r1.Contains(sw1));
        Assert.True(r1.Contains(ne1));
        Assert.False(r1.InteriorContains(sw1));
        Assert.False(r1.InteriorContains(ne1));

        // Make sure that GetVertex() returns vertices in CCW order.
        for (int k = 0; k < 4; ++k)
        {
            R2Point a = r1.GetVertex(k - 1);
            R2Point b = r1.GetVertex(k);
            R2Point c = r1.GetVertex(k + 1);
            Assert.True((b - a).GetOrtho().DotProd(c - a) > 0);
        }
    }

    [Fact]
    internal void Test_R2Rect_IntervalOperations()
    {
        // Contains(R2Rect), InteriorContains(R2Rect),
        // Intersects(), InteriorIntersects(), Union(), Intersection().
        //
        // Much more testing of these methods is done in s1interval_test
        // and r1interval_test.

        R2Rect empty = R2Rect.Empty;
        R2Point sw1 = new(0, 0.25);
        R2Point ne1 = new(0.5, 0.75);
        R2Rect r1 = new(sw1, ne1);
        R2Rect r1_mid = new(new R2Point(0.25, 0.5), new R2Point(0.25, 0.5));
        R2Rect r_sw1 = new(sw1, sw1);
        R2Rect r_ne1 = new(ne1, ne1);

        TestIntervalOps(r1, r1_mid, "TTTT", r1, r1_mid);
        TestIntervalOps(r1, r_sw1, "TFTF", r1, r_sw1);
        TestIntervalOps(r1, r_ne1, "TFTF", r1, r_ne1);

        Assert.Equal(new R2Rect(new R2Point(0, 0.25), new R2Point(0.5, 0.75)), r1);
        TestIntervalOps(r1, new R2Rect(new R2Point(0.45, 0.1), new R2Point(0.75, 0.3)), "FFTT",
                        new R2Rect(new R2Point(0, 0.1), new R2Point(0.75, 0.75)),
                        new R2Rect(new R2Point(0.45, 0.25), new R2Point(0.5, 0.3)));
        TestIntervalOps(r1, new R2Rect(new R2Point(0.5, 0.1), new R2Point(0.7, 0.3)), "FFTF",
                        new R2Rect(new R2Point(0, 0.1), new R2Point(0.7, 0.75)),
                        new R2Rect(new R2Point(0.5, 0.25), new R2Point(0.5, 0.3)));
        TestIntervalOps(r1, new R2Rect(new R2Point(0.45, 0.1), new R2Point(0.7, 0.25)), "FFTF",
                        new R2Rect(new R2Point(0, 0.1), new R2Point(0.7, 0.75)),
                        new R2Rect(new R2Point(0.45, 0.25), new R2Point(0.5, 0.25)));

        TestIntervalOps(new R2Rect(new R2Point(0.1, 0.2), new R2Point(0.1, 0.3)),
                        new R2Rect(new R2Point(0.15, 0.7), new R2Point(0.2, 0.8)), "FFFF",
                        new R2Rect(new R2Point(0.1, 0.2), new R2Point(0.2, 0.8)),
                        empty);

        // Check that the intersection of two rectangles that overlap in x but not y
        // is valid, and vice versa.
        TestIntervalOps(new R2Rect(new R2Point(0.1, 0.2), new R2Point(0.4, 0.5)),
                        new R2Rect(new R2Point(0, 0), new R2Point(0.2, 0.1)), "FFFF",
                        new R2Rect(new R2Point(0, 0), new R2Point(0.4, 0.5)), empty);
        TestIntervalOps(new R2Rect(new R2Point(0, 0), new R2Point(0.1, 0.3)),
                        new R2Rect(new R2Point(0.2, 0.1), new R2Point(0.3, 0.4)), "FFFF",
                        new R2Rect(new R2Point(0, 0), new R2Point(0.3, 0.4)), empty);
    }

    [Fact]
    internal void Test_R2Rect_AddPoint()
    {
        // AddPoint()
        R2Point sw1 = new(0, 0.25);
        R2Point ne1 = new(0.5, 0.75);
        R2Rect r1 = new(sw1, ne1);

        R2Rect r2 = R2Rect.Empty;
        r2 = r2.AddPoint(new R2Point(0, 0.25));
        r2 = r2.AddPoint(new R2Point(0.5, 0.25));
        r2 = r2.AddPoint(new R2Point(0, 0.75));
        r2 = r2.AddPoint(new R2Point(0.1, 0.4));
        Assert.Equal(r1, r2);
    }

    [Fact]
    internal void Test_R2Rect_Project()
    {
        R2Rect r1 = new(new R1Interval(0, 0.5), new R1Interval(0.25, 0.75));

        Assert.Equal(new R2Point(0, 0.25), r1.Project(new R2Point(-0.01, 0.24)));
        Assert.Equal(new R2Point(0, 0.48), r1.Project(new R2Point(-5.0, 0.48)));
        Assert.Equal(new R2Point(0, 0.75), r1.Project(new R2Point(-5.0, 2.48)));
        Assert.Equal(new R2Point(0.19, 0.75), r1.Project(new R2Point(0.19, 2.48)));
        Assert.Equal(new R2Point(0.5, 0.75), r1.Project(new R2Point(6.19, 2.48)));
        Assert.Equal(new R2Point(0.5, 0.53), r1.Project(new R2Point(6.19, 0.53)));
        Assert.Equal(new R2Point(0.5, 0.25), r1.Project(new R2Point(6.19, -2.53)));
        Assert.Equal(new R2Point(0.33, 0.25), r1.Project(new R2Point(0.33, -2.53)));
        Assert.Equal(new R2Point(0.33, 0.37), r1.Project(new R2Point(0.33, 0.37)));
    }

    [Fact]
    internal void Test_R2Rect_Expanded()
    {
        // Expanded()
        Assert.True(R2Rect.Empty.Expanded(new R2Point(0.1, 0.3)).IsEmpty());
        Assert.True(R2Rect.Empty.Expanded(new R2Point(-0.1, -0.3)).IsEmpty());
        Assert.True(new R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                    Expanded(new R2Point(0.1, 0.3)).
                    ApproxEquals(new R2Rect(new R2Point(0.1, 0.1), new R2Point(0.4, 1.0))));
        Assert.True(new R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                    Expanded(new R2Point(-0.1, 0.3)).IsEmpty());
        Assert.True(new R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                    Expanded(new R2Point(0.1, -0.2)).IsEmpty());
        Assert.True(new R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).
                    Expanded(new R2Point(0.1, -0.1)).
                    ApproxEquals(new R2Rect(new R2Point(0.1, 0.5), new R2Point(0.4, 0.6))));
        Assert.True(new R2Rect(new R2Point(0.2, 0.4), new R2Point(0.3, 0.7)).Expanded(0.1).
                    ApproxEquals(new R2Rect(new R2Point(0.1, 0.3), new R2Point(0.4, 0.8))));
    }

    private static void TestIntervalOps(R2Rect x, R2Rect y, string expected_rexion,
        R2Rect expected_union, R2Rect expected_intersection)
    {
        // Test all of the interval operations on the given pair of intervals.
        // "expected_rexion" is a sequence of "T" and "F" characters corresponding
        // to the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.

        Assert.Equal(expected_rexion[0] == 'T', x.Contains(y));
        Assert.Equal(expected_rexion[1] == 'T', x.InteriorContains(y));
        Assert.Equal(expected_rexion[2] == 'T', x.Intersects(y));
        Assert.Equal(expected_rexion[3] == 'T', x.InteriorIntersects(y));

        Assert.Equal(x.Union(y) == x, x.Contains(y));
        Assert.Equal(!x.Intersection(y).IsEmpty(), x.Intersects(y));

        Assert.Equal(expected_union, x.Union(y));
        Assert.Equal(expected_intersection, x.Intersection(y));

        R2Rect r = x;
        r = r.AddRect(y);
        Assert.Equal(expected_union, r);
        if (y.GetSize() == new R2Point(0, 0))
        {
            r = x;
            r = r.AddPoint(y.Lo());
            Assert.Equal(expected_union, r);
        }
    }
}
