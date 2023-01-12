namespace S2Geometry;

public class S1IntervalTestsBase
{
    #region Fields, Constants

    private readonly S1Interval empty, full;
    private readonly S1Interval zero, pi2, pi, mipi, mipi2;
    private readonly S1Interval quad1, quad2, quad3, quad4;
    private readonly S1Interval quad12, quad23, quad34, quad41;
    private readonly S1Interval quad123, quad234, quad341, quad412;
    private readonly S1Interval mid12, mid23, mid34, mid41;

    #endregion

    #region Constructor

    // Create some standard intervals to use in the tests.  These include the
    // empty and full intervals, intervals containing a single point, and
    // intervals spanning one or more "quadrants" which are numbered as follows:
    //    quad1 == [0, Pi/2]
    //    quad2 == [Pi/2, Pi]
    //    quad3 == [-Pi, -Pi/2]
    //    quad4 == [-Pi/2, 0]
    public S1IntervalTestsBase()
    {
        empty = S1Interval.Empty;
        full = S1Interval.Full;
        // Single-point intervals:
        zero = new S1Interval(0, 0);
        pi2 = new S1Interval(S2.M_PI_2, S2.M_PI_2);
        pi = new S1Interval(Math.PI, Math.PI);
        mipi = new S1Interval(-Math.PI, -Math.PI);  // Same as "pi" after normalization.
        mipi2 = new S1Interval(-S2.M_PI_2, -S2.M_PI_2);
        // Single quadrants:
        quad1 = new S1Interval(0, S2.M_PI_2);
        quad2 = new S1Interval(S2.M_PI_2, -Math.PI);
        quad3 = new S1Interval(Math.PI, -S2.M_PI_2);
        quad4 = new S1Interval(-S2.M_PI_2, 0);
        // Quadrant pairs:
        quad12 = new S1Interval(0, -Math.PI);
        quad23 = new S1Interval(S2.M_PI_2, -S2.M_PI_2);
        quad34 = new S1Interval(-Math.PI, 0);
        quad41 = new S1Interval(-S2.M_PI_2, S2.M_PI_2);
        // Quadrant triples:
        quad123 = new S1Interval(0, -S2.M_PI_2);
        quad234 = new S1Interval(S2.M_PI_2, 0);
        quad341 = new S1Interval(Math.PI, S2.M_PI_2);
        quad412 = new S1Interval(-S2.M_PI_2, -Math.PI);
        // Small intervals around the midpoints between quadrants, such that
        // the center of each interval is offset slightly CCW from the midpoint.
        mid12 = new S1Interval(S2.M_PI_2 - 0.01, S2.M_PI_2 + 0.02);
        mid23 = new S1Interval(Math.PI - 0.01, -Math.PI + 0.02);
        mid34 = new S1Interval(-S2.M_PI_2 - 0.01, -S2.M_PI_2 + 0.02);
        mid41 = new S1Interval(-0.01, 0.02);
    }

    #endregion

    [Fact]
    internal void Test_S1IntervalTestBase_ConstructorsAndAccessors()
    {
        // Spot-check the constructors and accessors.
        Assert.Equal(0, quad12.Lo);
        Assert.Equal(quad12.Hi, Math.PI);
        Assert.Equal(quad34[0], Math.PI);
        Assert.Equal(0, quad34[1]);
        Assert.Equal(quad34.Bounds(), new R2Point(Math.PI, 0).Bounds());
        Assert.Equal(pi.Lo, Math.PI);
        Assert.Equal(pi.Hi, Math.PI);

        // Check that [-Pi, -Pi] is normalized to [Pi, Pi].
        Assert.Equal(mipi.Lo, Math.PI);
        Assert.Equal(mipi.Hi, Math.PI);
        Assert.Equal(quad23.Lo, S2.M_PI_2);
        Assert.Equal(quad23.Hi, -S2.M_PI_2);

        // Check that the default S1Interval is identical to Empty().
        S1Interval default_empty = S1Interval.Empty;
        Assert.True(default_empty.IsValid());
        Assert.True(default_empty.IsEmpty());
        Assert.Equal(empty.Lo, default_empty.Lo);
        Assert.Equal(empty.Hi, default_empty.Hi);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_SimplePredicates()
    {
        // is_valid(), IsEmpty, IsFull, IsInverted
        Assert.True(zero.IsValid() && !zero.IsEmpty() && !zero.IsFull());
        Assert.True(empty.IsValid() && empty.IsEmpty() && !empty.IsFull());
        Assert.True(empty.IsInverted());
        Assert.True(full.IsValid() && !full.IsEmpty() && full.IsFull());
        Assert.True(!quad12.IsEmpty() && !quad12.IsFull() && !quad12.IsInverted());
        Assert.True(!quad23.IsEmpty() && !quad23.IsFull() && quad23.IsInverted());
        Assert.True(pi.IsValid() && !pi.IsEmpty() && !pi.IsInverted());
        Assert.True(mipi.IsValid() && !mipi.IsEmpty() && !mipi.IsInverted());
    }

    [Fact]
    internal void Test_S1IntervalTestBase_AlmostEmptyOrFull()
    {
        // Test that rounding errors don't cause intervals that are almost empty or
        // full to be considered empty or full.  The following value is the greatest
        // representable value less than Pi.
        double kAlmostPi = Math.PI - 2 * S2.DoubleEpsilon;
        Assert.False(new S1Interval(-kAlmostPi, Math.PI).IsFull());
        Assert.False(new S1Interval(-Math.PI, kAlmostPi).IsFull());
        Assert.False(new S1Interval(Math.PI, -kAlmostPi).IsEmpty());
        Assert.False(new S1Interval(kAlmostPi, -Math.PI).IsEmpty());
    }

    [Fact]
    internal void Test_S1IntervalTestBase_GetCenter()
    {
        Assert.Equal(quad12.GetCenter(), S2.M_PI_2);
        Assert2.DoubleEqual(new S1Interval(3.1, 2.9).GetCenter(), 3.0 - Math.PI);
        Assert2.DoubleEqual(new S1Interval(-2.9, -3.1).GetCenter(), Math.PI - 3.0);
        Assert2.DoubleEqual(new S1Interval(2.1, -2.1).GetCenter(), Math.PI);
        Assert.Equal(pi.GetCenter(), Math.PI);
        Assert.Equal(mipi.GetCenter(), Math.PI);
        Assert.Equal(Math.Abs(quad23.GetCenter()), Math.PI);
        Assert2.DoubleEqual(quad123.GetCenter(), 0.75 * Math.PI);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_GetLength()
    {
        Assert.Equal(quad12.GetLength(), Math.PI);
        Assert.Equal(0, pi.GetLength());
        Assert.Equal(0, mipi.GetLength());
        Assert2.DoubleEqual(quad123.GetLength(), 1.5 * Math.PI);
        Assert.Equal(Math.Abs(quad23.GetLength()), Math.PI);
        Assert.Equal(full.GetLength(), S2.M_2_PI);
        Assert.True(empty.GetLength() < 0);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_Complement()
    {
        Assert.True(empty.Complement().IsFull());
        Assert.True(full.Complement().IsEmpty());
        Assert.True(pi.Complement().IsFull());
        Assert.True(mipi.Complement().IsFull());
        Assert.True(zero.Complement().IsFull());
        Assert.True(quad12.Complement().ApproxEquals(quad34));
        Assert.True(quad34.Complement().ApproxEquals(quad12));
        Assert.True(quad123.Complement().ApproxEquals(quad4));
    }

    [Fact]
    internal void Test_S1IntervalTestBase_Contains()
    {
        // Contains(double), InteriorContains(double)
        Assert.True(!empty.Contains(0) && !empty.Contains(Math.PI) &&
                    !empty.Contains(-Math.PI));
        Assert.True(!empty.InteriorContains(Math.PI) && !empty.InteriorContains(-Math.PI));
        Assert.True(full.Contains(0) && full.Contains(Math.PI) && full.Contains(-Math.PI));
        Assert.True(full.InteriorContains(Math.PI) && full.InteriorContains(-Math.PI));
        Assert.True(quad12.Contains(0) && quad12.Contains(Math.PI) &&
                    quad12.Contains(-Math.PI));
        Assert.True(quad12.InteriorContains(S2.M_PI_2) && !quad12.InteriorContains(0));
        Assert.True(!quad12.InteriorContains(Math.PI) &&
                    !quad12.InteriorContains(-Math.PI));
        Assert.True(quad23.Contains(S2.M_PI_2) && quad23.Contains(-S2.M_PI_2));
        Assert.True(quad23.Contains(Math.PI) && quad23.Contains(-Math.PI));
        Assert.True(!quad23.Contains(0));
        Assert.True(!quad23.InteriorContains(S2.M_PI_2) &&
                    !quad23.InteriorContains(-S2.M_PI_2));
        Assert.True(quad23.InteriorContains(Math.PI) && quad23.InteriorContains(-Math.PI));
        Assert.True(!quad23.InteriorContains(0));
        Assert.True(pi.Contains(Math.PI) && pi.Contains(-Math.PI) && !pi.Contains(0));
        Assert.True(!pi.InteriorContains(Math.PI) && !pi.InteriorContains(-Math.PI));
        Assert.True(mipi.Contains(Math.PI) && mipi.Contains(-Math.PI) && !mipi.Contains(0));
        Assert.True(!mipi.InteriorContains(Math.PI) && !mipi.InteriorContains(-Math.PI));
        Assert.True(zero.Contains(0) && !zero.InteriorContains(0));
    }

    [Fact]
    internal void Test_S1IntervalTestBase_IntervalOps()
    {
        // Contains(S1Interval), InteriorContains(S1Interval),
        // Intersects(), InteriorIntersects(), Union(), Intersection()
        TestIntervalOps(empty, empty, "TTFF", empty, empty);
        TestIntervalOps(empty, full, "FFFF", full, empty);
        TestIntervalOps(empty, zero, "FFFF", zero, empty);
        TestIntervalOps(empty, pi, "FFFF", pi, empty);
        TestIntervalOps(empty, mipi, "FFFF", mipi, empty);

        TestIntervalOps(full, empty, "TTFF", full, empty);
        TestIntervalOps(full, full, "TTTT", full, full);
        TestIntervalOps(full, zero, "TTTT", full, zero);
        TestIntervalOps(full, pi, "TTTT", full, pi);
        TestIntervalOps(full, mipi, "TTTT", full, mipi);
        TestIntervalOps(full, quad12, "TTTT", full, quad12);
        TestIntervalOps(full, quad23, "TTTT", full, quad23);

        TestIntervalOps(zero, empty, "TTFF", zero, empty);
        TestIntervalOps(zero, full, "FFTF", full, zero);
        TestIntervalOps(zero, zero, "TFTF", zero, zero);
        TestIntervalOps(zero, pi, "FFFF", new S1Interval(0, Math.PI), empty);
        TestIntervalOps(zero, pi2, "FFFF", quad1, empty);
        TestIntervalOps(zero, mipi, "FFFF", quad12, empty);
        TestIntervalOps(zero, mipi2, "FFFF", quad4, empty);
        TestIntervalOps(zero, quad12, "FFTF", quad12, zero);
        TestIntervalOps(zero, quad23, "FFFF", quad123, empty);

        TestIntervalOps(pi2, empty, "TTFF", pi2, empty);
        TestIntervalOps(pi2, full, "FFTF", full, pi2);
        TestIntervalOps(pi2, zero, "FFFF", quad1, empty);
        TestIntervalOps(pi2, pi, "FFFF", new S1Interval(S2.M_PI_2, Math.PI), empty);
        TestIntervalOps(pi2, pi2, "TFTF", pi2, pi2);
        TestIntervalOps(pi2, mipi, "FFFF", quad2, empty);
        TestIntervalOps(pi2, mipi2, "FFFF", quad23, empty);
        TestIntervalOps(pi2, quad12, "FFTF", quad12, pi2);
        TestIntervalOps(pi2, quad23, "FFTF", quad23, pi2);

        TestIntervalOps(pi, empty, "TTFF", pi, empty);
        TestIntervalOps(pi, full, "FFTF", full, pi);
        TestIntervalOps(pi, zero, "FFFF", new S1Interval(Math.PI, 0), empty);
        TestIntervalOps(pi, pi, "TFTF", pi, pi);
        TestIntervalOps(pi, pi2, "FFFF", new S1Interval(S2.M_PI_2, Math.PI), empty);
        TestIntervalOps(pi, mipi, "TFTF", pi, pi);
        TestIntervalOps(pi, mipi2, "FFFF", quad3, empty);
        TestIntervalOps(pi, quad12, "FFTF", new S1Interval(0, Math.PI), pi);
        TestIntervalOps(pi, quad23, "FFTF", quad23, pi);

        TestIntervalOps(mipi, empty, "TTFF", mipi, empty);
        TestIntervalOps(mipi, full, "FFTF", full, mipi);
        TestIntervalOps(mipi, zero, "FFFF", quad34, empty);
        TestIntervalOps(mipi, pi, "TFTF", mipi, mipi);
        TestIntervalOps(mipi, pi2, "FFFF", quad2, empty);
        TestIntervalOps(mipi, mipi, "TFTF", mipi, mipi);
        TestIntervalOps(mipi, mipi2, "FFFF", new S1Interval(-Math.PI, -S2.M_PI_2), empty);
        TestIntervalOps(mipi, quad12, "FFTF", quad12, mipi);
        TestIntervalOps(mipi, quad23, "FFTF", quad23, mipi);

        TestIntervalOps(quad12, empty, "TTFF", quad12, empty);
        TestIntervalOps(quad12, full, "FFTT", full, quad12);
        TestIntervalOps(quad12, zero, "TFTF", quad12, zero);
        TestIntervalOps(quad12, pi, "TFTF", quad12, pi);
        TestIntervalOps(quad12, mipi, "TFTF", quad12, mipi);
        TestIntervalOps(quad12, quad12, "TFTT", quad12, quad12);
        TestIntervalOps(quad12, quad23, "FFTT", quad123, quad2);
        TestIntervalOps(quad12, quad34, "FFTF", full, quad12);

        TestIntervalOps(quad23, empty, "TTFF", quad23, empty);
        TestIntervalOps(quad23, full, "FFTT", full, quad23);
        TestIntervalOps(quad23, zero, "FFFF", quad234, empty);
        TestIntervalOps(quad23, pi, "TTTT", quad23, pi);
        TestIntervalOps(quad23, mipi, "TTTT", quad23, mipi);
        TestIntervalOps(quad23, quad12, "FFTT", quad123, quad2);
        TestIntervalOps(quad23, quad23, "TFTT", quad23, quad23);
        TestIntervalOps(quad23, quad34, "FFTT", quad234, new S1Interval(-Math.PI, -S2.M_PI_2));

        TestIntervalOps(quad1, quad23, "FFTF", quad123, new S1Interval(S2.M_PI_2, S2.M_PI_2));
        TestIntervalOps(quad2, quad3, "FFTF", quad23, mipi);
        TestIntervalOps(quad3, quad2, "FFTF", quad23, pi);
        TestIntervalOps(quad2, pi, "TFTF", quad2, pi);
        TestIntervalOps(quad2, mipi, "TFTF", quad2, mipi);
        TestIntervalOps(quad3, pi, "TFTF", quad3, pi);
        TestIntervalOps(quad3, mipi, "TFTF", quad3, mipi);

        TestIntervalOps(quad12, mid12, "TTTT", quad12, mid12);
        TestIntervalOps(mid12, quad12, "FFTT", quad12, mid12);

        S1Interval quad12eps = new(quad12.Lo, mid23.Hi);
        S1Interval quad2hi = new(mid23.Lo, quad12.Hi);
        TestIntervalOps(quad12, mid23, "FFTT", quad12eps, quad2hi);
        TestIntervalOps(mid23, quad12, "FFTT", quad12eps, quad2hi);

        // This test checks that the union of two disjoint intervals is the smallest
        // interval that contains both of them.  Note that the center of "mid34" is
        // slightly CCW of -Pi/2 so that there is no ambiguity about the result.
        S1Interval quad412eps = new(mid34.Lo, quad12.Hi);
        TestIntervalOps(quad12, mid34, "FFFF", quad412eps, empty);
        TestIntervalOps(mid34, quad12, "FFFF", quad412eps, empty);

        S1Interval quadeps12 = new(mid41.Lo, quad12.Hi);
        S1Interval quad1lo = new(quad12.Lo, mid41.Hi);
        TestIntervalOps(quad12, mid41, "FFTT", quadeps12, quad1lo);
        TestIntervalOps(mid41, quad12, "FFTT", quadeps12, quad1lo);

        S1Interval quad2lo = new(quad23.Lo, mid12.Hi);
        S1Interval quad3hi = new(mid34.Lo, quad23.Hi);
        S1Interval quadeps23 = new(mid12.Lo, quad23.Hi);
        S1Interval quad23eps = new(quad23.Lo, mid34.Hi);
        S1Interval quadeps123 = new(mid41.Lo, quad23.Hi);
        TestIntervalOps(quad23, mid12, "FFTT", quadeps23, quad2lo);
        TestIntervalOps(mid12, quad23, "FFTT", quadeps23, quad2lo);
        TestIntervalOps(quad23, mid23, "TTTT", quad23, mid23);
        TestIntervalOps(mid23, quad23, "FFTT", quad23, mid23);
        TestIntervalOps(quad23, mid34, "FFTT", quad23eps, quad3hi);
        TestIntervalOps(mid34, quad23, "FFTT", quad23eps, quad3hi);
        TestIntervalOps(quad23, mid41, "FFFF", quadeps123, empty);
        TestIntervalOps(mid41, quad23, "FFFF", quadeps123, empty);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_AddPoint()
    {
        Assert.Equal(zero, S1Interval.AddPoint(empty, 0));
        Assert.Equal(pi, S1Interval.AddPoint(empty, Math.PI));
        Assert.Equal(mipi, S1Interval.AddPoint(empty, -Math.PI));
        Assert.Equal(pi, S1Interval.AddPoints(empty, Math.PI, -Math.PI));
        Assert.Equal(mipi, S1Interval.AddPoints(empty, -Math.PI, Math.PI));
        Assert.Equal(mid12, S1Interval.AddPoints(empty, mid12.Lo, mid12.Hi));
        Assert.Equal(mid23, S1Interval.AddPoints(empty, mid23.Lo, mid23.Hi));
        Assert.Equal(quad123, S1Interval.AddPoints(quad1, -0.9 * Math.PI, -S2.M_PI_2));
        Assert.True(S1Interval.AddPoint(full, 0).IsFull());
        Assert.True(S1Interval.AddPoint(full, Math.PI).IsFull());
        Assert.True(S1Interval.AddPoint(full, -Math.PI).IsFull());
    }

    [Fact]
    internal void Test_S1IntervalTestBase_Project()
    {
        S1Interval r = new(-Math.PI, -Math.PI);
        Assert.Equal(Math.PI, r.Project(-Math.PI));
        Assert.Equal(Math.PI, r.Project(0));
        r = new S1Interval(0, Math.PI);
        Assert.Equal(0.1, r.Project(0.1));
        Assert.Equal(0, r.Project(-S2.M_PI_2 + S2.DoubleError));
        Assert.Equal(Math.PI, r.Project(-S2.M_PI_2 - S2.DoubleError));
        r = new S1Interval(Math.PI - 0.1, -Math.PI + 0.1);
        Assert.Equal(Math.PI, r.Project(Math.PI));
        Assert.Equal(Math.PI - 0.1, r.Project(S2.DoubleError));
        Assert.Equal(-Math.PI + 0.1, r.Project(-S2.DoubleError));
        Assert.Equal(0, S1Interval.Full.Project(0));
        Assert.Equal(Math.PI, S1Interval.Full.Project(Math.PI));
        Assert.Equal(Math.PI, S1Interval.Full.Project(-Math.PI));
    }

    [Fact]
    internal void Test_S1IntervalTestBase_FromPointPair()
    {
        Assert.Equal(S1Interval.FromPointPair(-Math.PI, Math.PI), pi);
        Assert.Equal(S1Interval.FromPointPair(Math.PI, -Math.PI), pi);
        Assert.Equal(S1Interval.FromPointPair(mid34.Hi, mid34.Lo), mid34);
        Assert.Equal(S1Interval.FromPointPair(mid23.Lo, mid23.Hi), mid23);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_Expanded()
    {
        Assert.Equal(empty.Expanded(1), empty);
        Assert.Equal(full.Expanded(1), full);
        Assert.Equal(zero.Expanded(1), new S1Interval(-1, 1));
        Assert.Equal(mipi.Expanded(0.01), new S1Interval(Math.PI - 0.01, -Math.PI + 0.01));
        Assert.Equal(pi.Expanded(27), full);
        Assert.Equal(pi.Expanded(S2.M_PI_2), quad23);
        Assert.Equal(pi2.Expanded(S2.M_PI_2), quad12);
        Assert.Equal(mipi2.Expanded(S2.M_PI_2), quad34);

        Assert.Equal(empty.Expanded(-1), empty);
        Assert.Equal(full.Expanded(-1), full);
        Assert.Equal(quad123.Expanded(-27), empty);
        Assert.Equal(quad234.Expanded(-27), empty);
        Assert.Equal(quad123.Expanded(-S2.M_PI_2), quad2);
        Assert.Equal(quad341.Expanded(-S2.M_PI_2), quad4);
        Assert.Equal(quad412.Expanded(-S2.M_PI_2), quad1);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_ApproxEquals()
    {
        // Choose two values kLo and kHi such that it's okay to shift an endpoint by
        // kLo (i.e., the resulting interval is equivalent) but not by kHi.
        const double kLo = 4 * S2.DoubleEpsilon;  // < max_error default
        const double kHi = 6 * S2.DoubleEpsilon;  // > max_error default

        // Empty intervals.
        Assert.True(empty.ApproxEquals(empty));
        Assert.True(zero.ApproxEquals(empty) && empty.ApproxEquals(zero));
        Assert.True(pi.ApproxEquals(empty) && empty.ApproxEquals(pi));
        Assert.True(mipi.ApproxEquals(empty) && empty.ApproxEquals(mipi));
        Assert.False(empty.ApproxEquals(full));
        Assert.True(empty.ApproxEquals(new S1Interval(1, 1 + 2 * kLo)));
        Assert.False(empty.ApproxEquals(new S1Interval(1, 1 + 2 * kHi)));
        Assert.True(new S1Interval(Math.PI - kLo, -Math.PI + kLo).ApproxEquals(empty));

        // Full intervals.
        Assert.True(full.ApproxEquals(full));
        Assert.False(full.ApproxEquals(empty));
        Assert.False(full.ApproxEquals(zero));
        Assert.False(full.ApproxEquals(pi));
        Assert.True(full.ApproxEquals(new S1Interval(kLo, -kLo)));
        Assert.False(full.ApproxEquals(new S1Interval(2 * kHi, 0)));
        Assert.True(new S1Interval(-Math.PI + kLo, Math.PI - kLo).ApproxEquals(full));
        Assert.False(new S1Interval(-Math.PI, Math.PI - 2 * kHi).ApproxEquals(full));

        // Singleton intervals.
        Assert.True(pi.ApproxEquals(pi) && mipi.ApproxEquals(pi));
        Assert.True(pi.ApproxEquals(new S1Interval(Math.PI - kLo, Math.PI - kLo)));
        Assert.False(pi.ApproxEquals(new S1Interval(Math.PI - kHi, Math.PI - kHi)));
        Assert.True(pi.ApproxEquals(new S1Interval(Math.PI - kLo, -Math.PI + kLo)));
        Assert.False(pi.ApproxEquals(new S1Interval(Math.PI - kHi, -Math.PI)));
        Assert.False(zero.ApproxEquals(pi));
        Assert.True(pi.Union(mid12).Union(zero).ApproxEquals(quad12));
        Assert.True(quad2.Intersection(quad3).ApproxEquals(pi));
        Assert.True(quad3.Intersection(quad2).ApproxEquals(pi));

        // Intervals whose corresponding endpoints are nearly the same but where the
        // endpoints are in opposite order (i.e., inverted intervals).
        Assert.False(new S1Interval(0, kLo).ApproxEquals(new S1Interval(kLo, 0)));
        Assert.False(new S1Interval(Math.PI - 0.5 * kLo, -Math.PI + 0.5 * kLo).
                     ApproxEquals(new S1Interval(-Math.PI + 0.5 * kLo, Math.PI - 0.5 * kLo)));

        // Other intervals.
        Assert.True(new S1Interval(1 - kLo, 2 + kLo).ApproxEquals(new S1Interval(1, 2)));
        Assert.True(new S1Interval(1 + kLo, 2 - kLo).ApproxEquals(new S1Interval(1, 2)));
        Assert.True(new S1Interval(2 - kLo, 1 + kLo).ApproxEquals(new S1Interval(2, 1)));
        Assert.True(new S1Interval(2 + kLo, 1 - kLo).ApproxEquals(new S1Interval(2, 1)));
        Assert.False(new S1Interval(1 - kHi, 2 + kLo).ApproxEquals(new S1Interval(1, 2)));
        Assert.False(new S1Interval(1 + kHi, 2 - kLo).ApproxEquals(new S1Interval(1, 2)));
        Assert.False(new S1Interval(2 - kHi, 1 + kLo).ApproxEquals(new S1Interval(2, 1)));
        Assert.False(new S1Interval(2 + kHi, 1 - kLo).ApproxEquals(new S1Interval(2, 1)));
        Assert.False(new S1Interval(1 - kLo, 2 + kHi).ApproxEquals(new S1Interval(1, 2)));
        Assert.False(new S1Interval(1 + kLo, 2 - kHi).ApproxEquals(new S1Interval(1, 2)));
        Assert.False(new S1Interval(2 - kLo, 1 + kHi).ApproxEquals(new S1Interval(2, 1)));
        Assert.False(new S1Interval(2 + kLo, 1 - kHi).ApproxEquals(new S1Interval(2, 1)));
    }

    [Fact]
    internal void Test_S1IntervalTestBase_OperatorEquals()
    {
        Assert.Equal(empty, empty);
        Assert.Equal(full, full);
        Assert.NotEqual(full, empty);
    }

    [Fact]
    internal void Test_S1IntervalTestBase_GetDirectedHausdorffDistance()
    {
        Assert2.Near(0.0, empty.GetDirectedHausdorffDistance(empty));
        Assert2.Near(0.0, empty.GetDirectedHausdorffDistance(mid12));
        Assert2.Near(Math.PI, mid12.GetDirectedHausdorffDistance(empty));

        Assert.Equal(0.0, quad12.GetDirectedHausdorffDistance(quad123));
        S1Interval in_ = new(3.0, -3.0);  // an interval whose complement center is 0.
        Assert2.Near(3.0, new S1Interval(-0.1, 0.2).GetDirectedHausdorffDistance(in_));
        Assert2.Near(3.0 - 0.1, new S1Interval(0.1, 0.2).GetDirectedHausdorffDistance(in_));
        Assert2.Near(3.0 - 0.1, new S1Interval(-0.2, -0.1).GetDirectedHausdorffDistance(in_));
    }

    [Fact]
    internal void Test_S1IntervalTestBase_Quad41()
    {
        Assert.Equal(quad41, new S1Interval(-S2.M_PI_2, S2.M_PI_2));
    }

    private static void TestIntervalOps(S1Interval x, S1Interval y, string expected_relation, S1Interval expected_union, S1Interval expected_intersection)
    {
        // Test all of the interval operations on the given pair of intervals.
        // "expected_relation" is a sequence of "T" and "F" characters corresponding
        // to the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.

        Assert.Equal(x.Contains(y), expected_relation[0] == 'T');
        Assert.Equal(x.InteriorContains(y), expected_relation[1] == 'T');
        Assert.Equal(x.Intersects(y), expected_relation[2] == 'T');
        Assert.Equal(x.InteriorIntersects(y), expected_relation[3] == 'T');

        // bounds() returns a reference to a member variable, so we need to
        // make a copy when invoking it on a temporary object.
        Assert.Equal(R2Point.FromCoords(x.Union(y).Bounds()).Bounds(), expected_union.Bounds());
        Assert.Equal(R2Point.FromCoords(x.Intersection(y).Bounds()).Bounds(),
                  expected_intersection.Bounds());

        Assert.Equal(x.Contains(y), x.Union(y) == x);
        Assert.Equal(x.Intersects(y), !x.Intersection(y).IsEmpty());

        if (y.Lo == y.Hi)
        {
            var r = S1Interval.AddPoint(x, y.Lo);
            Assert.Equal(r.Bounds(), expected_union.Bounds());
        }
    }
}
