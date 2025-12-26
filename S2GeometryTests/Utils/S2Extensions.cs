namespace S2Geometry;

internal static class S2Extensions
{
    public static double[] Bounds(this R2Point me) => [me.X, me.Y];

    public static double[] Bounds(this S1Interval me) => [me.Lo, me.Hi];

    /// <summary>
    /// Exports the latitude and longitude in degrees, separated by a comma.
    /// e.g. "94.518000,150.300000"
    /// </summary>
    public static string ToStringInDegrees(this S2LatLng me)
    {
        var pt = me.Normalized();
        return FormattableString.Invariant($"{pt.LatDegrees():F6},{pt.LngDegrees():F6}");
    }
}
