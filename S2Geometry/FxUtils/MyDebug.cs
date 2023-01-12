namespace S2Geometry;

using System.Diagnostics;

internal static class MyDebug
{
    internal static void Assert(bool statement)
    {
        if (!statement) throw new InvalidDataException();
    }

    [Conditional("DEBUG")]
    internal static void WriteLine(string? message) => Debug.WriteLine(message);

    [Conditional("DEBUG")]
    internal static void Write(string? message) => Debug.Write(message);
}
