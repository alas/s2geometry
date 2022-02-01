namespace S2Geometry;

public static partial class S2
{
    public const double M_1_PI = 1.0 / Math.PI;
    /// <summary>
    /// Constant used to convert radians to degrees
    /// </summary>
    public const double M_180_PI = 180.0 / Math.PI;
    public const double M_PI = Math.PI;
    public const double M_PI_180 = Math.PI / 180.0;
    /// <summary>
    /// 1.570
    /// </summary>
    public const double M_PI_2 = Math.PI / 2.0;
    public const double M_PI_4 = Math.PI / 4.0;
    public const double M_2_PI = Math.PI * 2.0;
    public const double M_4_PI = Math.PI * 4.0;
    public const double M_N2_PI = Math.PI * -2.0;
    /// <summary>
    /// 1.414
    /// </summary>
    public static readonly double M_SQRT2 = Math.Sqrt(2.0);
    public static readonly double M_SQRT1_2 = 1.0 / M_SQRT2;

    /// <summary>
    /// This is the number of levels needed to specify a leaf cell.  This
    /// constant is defined here so that the S2Metrics.Metric class and the conversion
    /// functions in S2Coords can be implemented without using S2CellId.
    /// </summary>
    public const int kMaxCellLevel = 30;

    /// <summary>
    /// The maximum index of a valid leaf cell plus one.  The range of valid leaf
    /// cell indices is [0..kLimitIJ-1].
    /// </summary>
    public const int kLimitIJ = 1 << kMaxCellLevel;  // == S2CellId.kMaxSize

    /// <summary>
    /// The maximum value of an si- or ti-coordinate.  The range of valid (si,ti)
    /// values is [0..kMaxSiTi].
    /// </summary>
    public const uint kMaxSiTi = 1U << (kMaxCellLevel + 1);

    // Together these flags define a cell orientation.  If 'kSwapMask'
    // is true, then canonical traversal order is flipped around the
    // diagonal (i.e. i and j are swapped with each other).  If
    // 'kInvertMask' is true, then the traversal order is rotated by 180
    // degrees (i.e. the bits of i and j are inverted, or equivalently,
    // the axis directions are reversed).
    public const int kSwapMask = 0x01;
    public const int kInvertMask = 0x02;

    public const byte kCurrentLosslessEncodingVersionNumber = 1;

    #region C++ Double Constants

    // C++ constants do not map well to C# constants
    // do not use double.Epsilon to avoid confusion, it is orders of
    // magnitude smaller than C++ DBL_EPSILON

    public const double DoubleError = 1e-15;

    // Difference between 1.0 and the minimum double greater than 1.0
    // (std.numeric_limits<double>.epsilon / DBL_EPSILON)
    public const double DoubleEpsilon = 2.2204460492503131e-16;

    // Minimum positive normalized value
    // (std.numeric_limits<double>.min() / DBL_MIN)
    public const double DoubleMinNorm = 2.22507e-308;

    #endregion
}
