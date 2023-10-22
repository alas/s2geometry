namespace S2Geometry;

using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;

internal static class MyDebug
{
    [Conditional("DEBUG")]
    internal static void Assert([DoesNotReturnIf(false)] bool statement, string? message = null)
    {
        if (!statement) throw new InvalidDataException(message ?? string.Empty);
    }

    [Conditional("DEBUG")]
    internal static void WriteLine(string? message) => Debug.WriteLine(message);

    [Conditional("DEBUG")]
    internal static void Write(string? message) => Debug.Write(message);
}
