
// The S2 library defines extra validity checks throughout the code that can
// optionally be enabled or disabled.  By default, these validity checks are
// enabled in debug-mode builds (including fastbuild) and disabled in
// optimized builds.
//
// There are two ways to change the default behavior:
//
//  - The command line --s2debug flag, which changes the global default.
//
//  - The S2Debug enum, which allows validity checks to be enabled or disabled
//    for specific objects (e.g., an S2Polygon).
//
// If you want to intentionally create invalid geometry (e.g., in a test), the
// S2Debug enum is preferable.  For example, to create an invalid S2Polygon,
// you can do this:
//
//   S2Polygon invalid;
//   invalid.set_s2debug_override(S2Debug::DISABLE);
//
// There is also a convenience constructor:
//
//   vector<unique_ptr<S2Loop>> loops = ...;
//   S2Polygon invalid(loops, S2Debug::DISABLE);
//
// There are a few checks that cannot be disabled this way (e.g., internal
// functions that require S2Points to be unit length).  If you absolutely need
// to disable these checks, you can set FLAGS_s2debug for the duration of a
// specific test like this:
//
// TEST(MyClass, InvalidGeometry) {
//   FLAGS_s2debug = false;  // Automatically restored between tests
//   ...
// }


// Command line flag that enables extra validity checking throughout the S2
// code.  It is turned on by default in debug-mode builds.
//S2_DECLARE_bool(s2debug);

namespace S2Geometry;

// Class that allows the --s2debug validity checks to be enabled or disabled
// for specific objects (e.g., see S2Polygon).
public enum S2Debug : byte
{
    ALLOW,    // Validity checks are controlled by --s2debug
    DISABLE   // No validity checks even when --s2debug is true
}
