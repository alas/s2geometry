using System.Collections.Generic;

namespace S2Geometry
{
    // Author: jyrki@google.com (Jyrki Alakuijala)
    //
    // Fast mixing of hash values -- not strong enough for fingerprinting.
    // May change from time to time.
    //
    // Values given are expected to be hashes from good hash functions.
    // What constitutes a good hash may depend on your application. As a rule of
    // thumb, if std::hash<int> is strong enough for your hashing need if
    // your data were just ints, it will most likely be the correct choice
    // for a mixed hash of data members. HashMix does one round of multiply and
    // rotate mixing, so you get some additional collision avoidance guarantees
    // compared to just using std::hash<int> directly.
    //
    // Possible use:
    //
    // struct Xyzzy {
    //   int x;
    //   int y;
    //   string id;
    // };
    //
    // #ifndef SWIG
    // template<> struct XyzzyHash<Xyzzy> {
    //   size_t operator()(const Xyzzy& c) const {
    //     HashMix mix(hash<int>()(c.x));
    //     mix.Mix(hash<int>()(c.y));
    //     mix.Mix(GoodFastHash<string>()(c.id));
    //     return mix.get();
    //   }
    // }
    // #endif
    //
    // HashMix is a lower level interface than std::hash<std::tuple<>>.
    // Use std::hash<std::tuple<>> instead of HashMix where appropriate.
    public static class HashMix
    {
        private const ulong kMul = 0xdc3eb94af8ab4c93UL;
        public static ulong GetHash<T>(IEnumerable<T> values)
        {
            ulong mix = 1;
            foreach (var value in values)
            {
                // Multiplicative hashing will mix bits better in the msb end ...
                mix *= kMul;
                // ... and rotating will move the better mixed msb-bits to lsb-bits.
                mix = ((mix << 19) | (mix >> (64 - 19))) + (ulong)(uint)value.GetHashCode();
            }
            return mix;
        }
    }
}
