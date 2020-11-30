using System;
namespace S2Geometry
{
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
    // sadly i can't cleanly make it generic in c# so I only make it for double and S2Point
    public readonly struct S2PointVector3 : IEquatable<S2PointVector3>
    {
        #region Fields, Constants

        private readonly double m_00;
        private readonly double m_01;
        private readonly double m_02;
        private readonly double m_10;
        private readonly double m_11;
        private readonly double m_12;
        private readonly double m_20;
        private readonly double m_21;
        private readonly double m_22;

        #endregion

        #region Constructors

        // Constructor explicitly setting the values of all the coefficient of
        // the matrix
        public S2PointVector3(double m00, double m01, double m02,
                         double m10, double m11, double m12,
                         double m20, double m21, double m22)
        {
            m_00 = m00;
            m_01 = m01;
            m_02 = m02;

            m_10 = m10;
            m_11 = m11;
            m_12 = m12;

            m_20 = m20;
            m_21 = m21;
            m_22 = m22;
        }

        public S2PointVector3(S2Point a, S2Point b, S2Point c)
        {
            m_00 = a[0];
            m_01 = a[1];
            m_02 = a[2];

            m_10 = b[0];
            m_11 = b[1];
            m_12 = b[2];

            m_20 = c[0];
            m_21 = c[1];
            m_22 = c[2];
        }

        #endregion

        #region S2PointVector3

        // Return matrix element (i,j) with 0<=i<=2 0<=j<=2
        public S2Point this[int index]
        {
            get
            {
                return index switch
                {
                    0 => new S2Point(m_00, m_01, m_02),
                    1 => new S2Point(m_10, m_11, m_12),
                    2 => new S2Point(m_20, m_21, m_22),
                    _ => throw new ArgumentException(null, nameof(index)),
                };
            }
        }

        public S2Point Col(int index)
        {
            return index switch
            {
                0 => new S2Point(m_00, m_10, m_20),
                1 => new S2Point(m_01, m_11, m_21),
                2 => new S2Point(m_02, m_12, m_22),
                _ => throw new ArgumentException(null, nameof(index)),
            };
        }

        // Return the determinant of the matrix
        public double Det() =>
            m_00 * m_11 * m_22
            + m_01 * m_12 * m_20
            + m_02 * m_10 * m_21
            - m_20 * m_11 * m_02
            - m_21 * m_12 * m_00
            - m_22 * m_10 * m_01;

        // Return the trace of the matrix
        public double Trace() => m_00 + m_11 + m_22;

        public double GetIJ(int i, int j)
        {
            var v = this[i];
            return j switch
            {
                0 => v[0],
                1 => v[1],
                2 => v[2],
                _ => throw new ArgumentException(null, nameof(j)),
            };
        }

        // Return the transposed matrix
        public S2PointVector3 Transpose() => new S2PointVector3(
            m_00, m_10, m_20,
            m_01, m_11, m_21,
            m_02, m_12, m_22);

        // Return the transposed of the matrix of the cofactors
        // (Useful for inversion for example)
        public S2PointVector3 ComatrixTransposed() => new S2PointVector3(
            m_11 * m_22 - m_21 * m_12,
            m_21 * m_02 - m_01 * m_22,
            m_01 * m_12 - m_11 * m_02,

            m_12 * m_20 - m_22 * m_10,
            m_22 * m_00 - m_02 * m_20,
            m_02 * m_10 - m_12 * m_00,

            m_10 * m_21 - m_20 * m_11,
            m_20 * m_01 - m_00 * m_21,
            m_00 * m_11 - m_10 * m_01);

        // Matrix inversion
        public S2PointVector3 Inverse()
        {
            var det = Det();
            Assert.True(det != 0); // Can't inverse. Determinant = 0.
            return (1 / det) * ComatrixTransposed();
        }

        // Return true is one of the elements of the matrix is NaN
        public bool IsNaN()
        {
            for (var i = 0; i < 3; ++i)
            {
                for (var j = 0; j < 3; ++j)
                {
                    if (double.IsNaN(this[i][j]))
                    {
                        return true;
                    }
                }
            }
            return false;
        } 
        
        #endregion

        #region IEquatable

        public override int GetHashCode()
        {
            var hash = HashCode.Combine(
                Math.Abs(m_00),
                Math.Abs(m_01),
                Math.Abs(m_02),
                Math.Abs(m_10),
                Math.Abs(m_11),
                Math.Abs(m_12));
            return HashCode.Combine(hash,
                Math.Abs(m_20),
                Math.Abs(m_21),
                Math.Abs(m_22));
        }

        public bool Equals(S2PointVector3 other) =>
            m_00 == other.m_00 &&
            m_01 == other.m_01 &&
            m_02 == other.m_02 &&
            m_10 == other.m_10 &&
            m_11 == other.m_11 &&
            m_12 == other.m_12 &&
            m_20 == other.m_20 &&
            m_21 == other.m_21 &&
            m_22 == other.m_22;

        public override bool Equals(object obj) => obj is S2PointVector3 mat && Equals(mat);

        public static bool operator ==(S2PointVector3 left, S2PointVector3 right) => Equals(left, right);

        public static bool operator !=(S2PointVector3 left, S2PointVector3 right) => !Equals(left, right); 

        #endregion

        #region Operators

        // Matrix addition
        public static S2PointVector3 operator +(S2PointVector3 left, S2PointVector3 right) => new S2PointVector3(
            left.m_00 + right.m_00,
            left.m_01 + right.m_01,
            left.m_02 + right.m_02,

            left.m_10 + right.m_10,
            left.m_11 + right.m_11,
            left.m_12 + right.m_12,

            left.m_20 + right.m_20,
            left.m_21 + right.m_21,
            left.m_22 + right.m_22);

        // Matrix subtration
        public static S2PointVector3 operator -(S2PointVector3 left, S2PointVector3 right) => new S2PointVector3(
            left.m_00 - right.m_00,
            left.m_01 - right.m_01,
            left.m_02 - right.m_02,

            left.m_10 - right.m_10,
            left.m_11 - right.m_11,
            left.m_12 - right.m_12,

            left.m_20 - right.m_20,
            left.m_21 - right.m_21,
            left.m_22 - right.m_22);

        // Matrix multiplication by a scalar
        public static S2PointVector3 operator *(S2PointVector3 left, double k) => new S2PointVector3(
            left.m_00 * k,
            left.m_01 * k,
            left.m_02 * k,

            left.m_10 * k,
            left.m_11 * k,
            left.m_12 * k,

            left.m_20 * k,
            left.m_21 * k,
            left.m_22 * k);

        // Matrix multiplication by a scalar
        public static S2PointVector3 operator *(double k, S2PointVector3 m) => m * k;

        // Matrix multiplication
        public static S2PointVector3 operator *(S2PointVector3 left, S2PointVector3 right)
        {
            return new S2PointVector3(
                left.m_00 * right.m_00 + left.m_01 * right.m_10 + left.m_02 * right.m_20,
                left.m_00 * right.m_01 + left.m_01 * right.m_11 + left.m_02 * right.m_21,
                left.m_00 * right.m_02 + left.m_01 * right.m_12 + left.m_02 * right.m_22,

                left.m_10 * right.m_00 + left.m_11 * right.m_10 + left.m_12 * right.m_20,
                left.m_10 * right.m_01 + left.m_11 * right.m_11 + left.m_12 * right.m_21,
                left.m_10 * right.m_02 + left.m_11 * right.m_12 + left.m_12 * right.m_22,

                left.m_20 * right.m_00 + left.m_21 * right.m_10 + left.m_22 * right.m_20,
                left.m_20 * right.m_01 + left.m_21 * right.m_11 + left.m_22 * right.m_21,
                left.m_20 * right.m_02 + left.m_21 * right.m_12 + left.m_22 * right.m_22);
        }

        // Multiplication of a matrix by a vector
        public static S2Point operator *(S2PointVector3 left, S2Point v)
        {
            return new S2Point(
                left.m_00 * v[0] + left.m_01 * v[1] + left.m_02 * v[2],
                left.m_10 * v[0] + left.m_11 * v[1] + left.m_12 * v[2],
                left.m_20 * v[0] + left.m_21 * v[1] + left.m_22 * v[2]);
        }

        public static S2Point operator *(S2Point v, S2PointVector3 m) => m * v;

        #endregion

        #region Object

        public override string ToString() =>
$@"[{m_00} {m_01} {m_02} 
 {m_10} {m_11} {m_12} 
 {m_20} {m_21} {m_22} ]"; 
        
        #endregion
    }
}
