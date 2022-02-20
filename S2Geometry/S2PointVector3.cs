namespace S2Geometry;

// A simple class to handle 3x3 matrices
// The aim of this class is to be able to manipulate 3x3 matrices
// and 3D vectors as naturally as possible and make calculations
// readable.
// For that reason, the operators +, -, * are overloaded.
// (Reading a = a + b*2 - c is much easier to read than
// a = Sub(Add(a, Mul(b,2)),c)   )
//
// Please be careful about overflows when using those matrices wth integer types
// The calculations are carried with double. eg : if you are using uint8 as the
// base type, all values will be modulo 256.
// This feature is necessary to use the class in a more general framework with
// double != plain old data type.
//
// Note(Alas): sadly I can't cleanly make it generic due to limitations in c# so I only added support for double and S2Point
public readonly record struct S2PointVector3
{
    #region Fields, Constants

    private readonly double M00 { get; }
    private readonly double M01 { get; }
    private readonly double M02 { get; }
    private readonly double M10 { get; }
    private readonly double M11 { get; }
    private readonly double M12 { get; }
    private readonly double M20 { get; }
    private readonly double M21 { get; }
    private readonly double M22 { get; }

    #endregion

    #region Constructors

    public S2PointVector3(
        double m00, double m01, double m02,
        double m10, double m11, double m12,
        double m20, double m21, double m22)
    {
        M00 = m00; M01 = m01; M02 = m02;
        M10 = m10; M11 = m11; M12 = m12;
        M20 = m20; M21 = m21; M22 = m22;
    }

    public S2PointVector3(S2Point a, S2Point b, S2Point c) : this(
        a[0], a[1], a[2],
        b[0], b[1], b[2],
        c[0], c[1], c[2])
    {
    }

    #endregion

    #region S2PointVector3

    // Return matrix element (i,j) with 0<=i<=2 0<=j<=2
    public S2Point Row(int index) => index switch
    {
        0 => new S2Point(M00, M01, M02),
        1 => new S2Point(M10, M11, M12),
        2 => new S2Point(M20, M21, M22),
        _ => throw new ArgumentException(null, nameof(index)),
    };

    public S2Point Col(int index) => index switch
    {
        0 => new S2Point(M00, M10, M20),
        1 => new S2Point(M01, M11, M21),
        2 => new S2Point(M02, M12, M22),
        _ => throw new ArgumentException(null, nameof(index)),
    };

    // Return the determinant of the matrix
    public double Det() =>
        M00 * M11 * M22
        + M01 * M12 * M20
        + M02 * M10 * M21
        - M20 * M11 * M02
        - M21 * M12 * M00
        - M22 * M10 * M01;

    // Return the trace of the matrix
    public double Trace() => M00 + M11 + M22;

    public double GetIJ(int i, int j)
    {
        var v = Row(i);
        return j switch
        {
            0 => v[0],
            1 => v[1],
            2 => v[2],
            _ => throw new ArgumentException(null, nameof(j)),
        };
    }

    // Return the transposed matrix
    public S2PointVector3 Transpose() => new(
        M00, M10, M20,
        M01, M11, M21,
        M02, M12, M22);

    // Return the transposed of the matrix of the cofactors
    // (Useful for inversion for example)
    public S2PointVector3 ComatrixTransposed() => new(
        M11 * M22 - M21 * M12,
        M21 * M02 - M01 * M22,
        M01 * M12 - M11 * M02,

        M12 * M20 - M22 * M10,
        M22 * M00 - M02 * M20,
        M02 * M10 - M12 * M00,

        M10 * M21 - M20 * M11,
        M20 * M01 - M00 * M21,
        M00 * M11 - M10 * M01);

    // Matrix inversion
    public S2PointVector3 Inverse()
    {
        var det = Det();
        System.Diagnostics.Debug.Assert(det != 0); // Can't inverse. Determinant = 0.
        return (1 / det) * ComatrixTransposed();
    }

    // Return true is one of the elements of the matrix is NaN
    public bool IsNaN()
    {
        for (var i = 0; i < 3; ++i)
        {
            for (var j = 0; j < 3; ++j)
            {
                if (double.IsNaN(Row(i)[j]))
                {
                    return true;
                }
            }
        }
        return false;
    }

    #endregion

    #region Operators

    // Matrix addition
    public static S2PointVector3 operator +(S2PointVector3 left, S2PointVector3 right) => new(
        left.M00 + right.M00,
        left.M01 + right.M01,
        left.M02 + right.M02,

        left.M10 + right.M10,
        left.M11 + right.M11,
        left.M12 + right.M12,

        left.M20 + right.M20,
        left.M21 + right.M21,
        left.M22 + right.M22);

    // Matrix subtration
    public static S2PointVector3 operator -(S2PointVector3 left, S2PointVector3 right) => new(
        left.M00 - right.M00,
        left.M01 - right.M01,
        left.M02 - right.M02,

        left.M10 - right.M10,
        left.M11 - right.M11,
        left.M12 - right.M12,

        left.M20 - right.M20,
        left.M21 - right.M21,
        left.M22 - right.M22);

    // Matrix multiplication by a scalar
    public static S2PointVector3 operator *(S2PointVector3 left, double k) => new(
        left.M00 * k,
        left.M01 * k,
        left.M02 * k,

        left.M10 * k,
        left.M11 * k,
        left.M12 * k,

        left.M20 * k,
        left.M21 * k,
        left.M22 * k);

    // Matrix multiplication by a scalar
    public static S2PointVector3 operator *(double k, S2PointVector3 m) => m * k;

    // Matrix multiplication
    public static S2PointVector3 operator *(S2PointVector3 left, S2PointVector3 right)
    {
        return new S2PointVector3(
            left.M00 * right.M00 + left.M01 * right.M10 + left.M02 * right.M20,
            left.M00 * right.M01 + left.M01 * right.M11 + left.M02 * right.M21,
            left.M00 * right.M02 + left.M01 * right.M12 + left.M02 * right.M22,

            left.M10 * right.M00 + left.M11 * right.M10 + left.M12 * right.M20,
            left.M10 * right.M01 + left.M11 * right.M11 + left.M12 * right.M21,
            left.M10 * right.M02 + left.M11 * right.M12 + left.M12 * right.M22,

            left.M20 * right.M00 + left.M21 * right.M10 + left.M22 * right.M20,
            left.M20 * right.M01 + left.M21 * right.M11 + left.M22 * right.M21,
            left.M20 * right.M02 + left.M21 * right.M12 + left.M22 * right.M22);
    }

    // Multiplication of a matrix by a vector
    public static S2Point operator *(S2PointVector3 left, S2Point v) => new(
            left.M00 * v[0] + left.M01 * v[1] + left.M02 * v[2],
            left.M10 * v[0] + left.M11 * v[1] + left.M12 * v[2],
            left.M20 * v[0] + left.M21 * v[1] + left.M22 * v[2]);

    public static S2Point operator *(S2Point v, S2PointVector3 m) => m * v;

    #endregion

    #region Object

    public override string ToString() =>
$@"[{M00} {M01} {M02} 
 {M10} {M11} {M12} 
 {M20} {M21} {M22} ]";

    #endregion
}
