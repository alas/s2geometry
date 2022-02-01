namespace S2Geometry;

internal static class Assert
{
    public static void Equal(double a, double b, double maxRatio = S2.DoubleError)
    {
        if (a == b) return;

        var full = a != 0 ? a : b;
        var dif = a - b;
        var difRatio = Math.Abs(dif / full);
        if (!(difRatio <= maxRatio)) throw new AssertionException("Assertion Failed");
    }
    public static void Equal<T>(T a, T b) where T : IEquatable<T>
    {
        if (!Equals(a, b)) throw new AssertionException("Assertion Failed");
    }
    public static void NotEqual<T>(T a, T b) where T : IEquatable<T>
    {
        if (Equals(a, b)) throw new AssertionException("Assertion Failed");
    }
    public static void Near(double a, double b, double absolute_range = 0.0001)
    {
        if (!(Math.Abs(a) - Math.Abs(b) <= absolute_range)) throw new AssertionException("Assertion Failed");
    }
    public static void True(bool expression)
    {
        if (!expression) throw new AssertionException("Assertion Failed");
    }
}

public class AssertionException : Exception
{
    public AssertionException(string message) : base(message) { }
}
