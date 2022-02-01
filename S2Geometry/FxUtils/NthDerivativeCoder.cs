//   far-far-superior implementation courtesy of amc@google.com (Adam Costello)
//
// Nth Derivative Coding
//   (In signal processing disciplines, this is known as N-th Delta Coding.)
//
// Good for varint coding integer sequences with polynomial trends.
//
// Instead of coding a sequence of values directly, code its nth-order discrete
// derivative.  Overflow in integer addition and subtraction makes this a
// lossless transform.
//
//                                       constant     linear      quadratic
//                                        trend       trend         trend
//                                      /        \  /        \  /           \_
// input                               |0  0  0  0  1  2  3  4  9  16  25  36
// 0th derivative(identity)            |0  0  0  0  1  2  3  4  9  16  25  36
// 1st derivative(delta coding)        |   0  0  0  1  1  1  1  5   7   9  11
// 2nd derivative(linear prediction)   |      0  0  1  0  0  0  4   2   2   2
//                                      -------------------------------------
//                                      0  1  2  3  4  5  6  7  8   9  10  11
//                                                  n in sequence
//
// Higher-order codings can break even or be detrimental on other sequences.
//
//                                           random            oscillating
//                                      /               \  /                  \_
// input                               |5  9  6  1   8  8  2 -2   4  -4   6  -6
// 0th derivative(identity)            |5  9  6  1   8  8  2 -2   4  -4   6  -6
// 1st derivative(delta coding)        |   4 -3 -5   7  0 -6 -4   6  -8  10 -12
// 2nd derivative(linear prediction)   |     -7 -2  12 -7 -6  2  10 -14  18 -22
//                                      ---------------------------------------
//                                      0  1  2  3  4   5  6  7   8   9  10  11
//                                                  n in sequence
//
// Note that the nth derivative isn't available until sequence item n.  Earlier
// values are coded at lower order.  For the above table, read 5 4 -7 -2 12 ...
//
// A caveat on class usage.  Encode() and Decode() share state.  Using both
// without a Reset() in-between probably doesn't make sense.

namespace S2Geometry;

public class NthDerivativeCoder
{
    // #if (~0 != -1)
    // #error Sorry, this code needs twos complement integers.
    // #endif

    // range of supported Ns: [ N_MIN, N_MAX ]
    private const int N_MIN = 0;
    private const int N_MAX = 10;

    // Initialize a new NthDerivativeCoder of the given N.
    public NthDerivativeCoder(int n)
    {
        N = n;
        if (n < N_MIN || n > N_MAX)
        {
            // Unsupported N, Using 0 instead.
            N = 0;
        }
        Reset();
    }

    // Encode the next value in the sequence.  Don't mix with Decode() calls.
    public Int32 Encode(Int32 k)
    {
        for (int i = 0; i < m_; ++i)
        {
            UInt32 delta = (UInt32)((UInt32)k - memory_[i]);
            memory_[i] = k;
            k = (Int32)delta;
        }
        if (m_ < N)
            memory_[m_++] = k;
        return k;
    }

    // Decode the next value in the sequence.  Don't mix with Encode() calls.
    public Int32 Decode(Int32 k)
    {
        if (m_ < N)
            m_++;
        for (int i = m_ - 1; i >= 0; --i)
            k = memory_[i] = (Int32)(memory_[i] + (UInt32)k);
        return k;
    }

    // Reset state.
    public void Reset()
    {
        for (int i = 0; i < N; ++i)
            memory_[i] = 0;
        m_ = 0;
    }

    // accessors
    public int N { get; }

    private int m_;  // the derivative order in which to code the next value(ramps to n_)
    private readonly Int32[] memory_ = new Int32[N_MAX];  // value memory. [0] is oldest
} 
