using Xunit;

namespace S2Geometry
{
    public class S2PolylineSimplifierTests
    {
        [Fact]
        public void Test_S2PolylineSimplifier_Reuse() {
            // Check that Init() can be called more than once.
            S1ChordAngle radius = new S1ChordAngle(S1Angle.FromDegrees(10));
            var s = new S2PolylineSimplifier(new S2Point(1, 0, 0));
            Assert.True(s.TargetDisc(new S2Point(1, 1, 0).Normalized, radius));
            Assert.True(s.TargetDisc(new S2Point(1, 1, 0.1).Normalized, radius));
            Assert.False(s.Extend(new S2Point(1, 1, 0.4).Normalized));

            // s.Init(S2Point(0, 1, 0));
            Assert.True(s.TargetDisc(new S2Point(1, 1, 0.3).Normalized, radius));
            Assert.True(s.TargetDisc(new S2Point(1, 1, 0.2).Normalized, radius));
            Assert.False(s.Extend(new S2Point(1, 1, 0).Normalized));
        }

        [Fact]
        public void Test_S2PolylineSimplifier_NoConstraints() {
            // Noraints, dst == src.
            CheckSimplify("0:1", "0:1", "", "", System.Array.Empty<bool>(), 0, true);

            // Noraints, dst != src.
            CheckSimplify("0:1", "1:0", "", "", System.Array.Empty<bool>(), 0, true);

            // Noraints, (src, dst) longer than 90 degrees (not supported).
            CheckSimplify("0:0", "0:91", "", "", System.Array.Empty<bool>(), 0, false);
        }

        [Fact]
        public void Test_S2PolylineSimplifier_TargetOnePoint() {
            // Three points on a straight line.  In theory zero tolerance should work,
            // but in practice there are floating point errors.
            CheckSimplify("0:0", "0:2", "0:1", "", System.Array.Empty<bool>(), 1e-10, true);

            // Three points where the middle point is too far away.
            CheckSimplify("0:0", "0:2", "1:1", "", System.Array.Empty<bool>(), 0.9, false);

            // A target disc that contains the source vertex.
            CheckSimplify("0:0", "0:2", "0:0.1", "", System.Array.Empty<bool>(), 1.0, true);

            // A target disc that contains the destination vertex.
            CheckSimplify("0:0", "0:2", "0:2.1", "", System.Array.Empty<bool>(), 1.0, true);
        }

        [Fact]
        public void Test_S2PolylineSimplifier_AvoidOnePoint() {
            // Three points on a straight line, attempting to avoid the middle point.
            CheckSimplify("0:0", "0:2", "", "0:1", new[] { true}, 1e-10, false);

            // Three points where the middle point can be successfully avoided.
            CheckSimplify("0:0", "0:2", "", "1:1", new[] { true}, 0.9, true);

            // Three points where the middle point is on the left, but where the client
            // requires the point to be on the right of the edge.
            CheckSimplify("0:0", "0:2", "", "1:1", new[] { false}, 1e-10, false);
        }

        [Fact]
        public void Test_S2PolylineSimplifier_TargetAndAvoid() {
            // Target several points that are separated from the proposed edge by about
            // 0.7 degrees, and avoid several points that are separated from the
            // proposed edge by about 1.4 degrees.
            CheckSimplify("0:0", "10:10", "2:3, 4:3, 7:8",
                          "4:2, 7:5, 7:9", new[]{ true, true, false}, 1.0, true);

            // The same example, but one point to be targeted is 1.4 degrees away.
            CheckSimplify("0:0", "10:10", "2:3, 4:6, 7:8",
                          "4:2, 7:5, 7:9", new[]{ true, true, false}, 1.0, false);

            // The same example, but one point to be avoided is 0.7 degrees away.
            CheckSimplify("0:0", "10:10", "2:3, 4:3, 7:8",
                          "4:2, 6:5, 7:9", new[]{ true, true, false}, 1.0, false);
        }

        [Fact]
        public void Test_S2PolylineSimplifier_Precision() {
            // This is a rough upper bound on both the error in constructing the disc
            // locations (i.e., S2.InterpolateAtDistance, etc), and also on the
            // padding that S2PolylineSimplifier uses to ensure that its results are
            // conservative (i.e., the error calculated by GetSemiwidth).
            S1Angle kMaxError = S1Angle.FromRadians(25 * S2Constants.DoubleEpsilon);

            // We repeatedly generate a random edge.  We then target several discs that
            // barely overlap the edge, and avoid several discs that barely miss the
            // edge.  About half the time, we choose one disc and make it slightly too
            // large or too small so that targeting fails.
            int kIters = 1000;  // Passes with 1 million iterations.
            for (int iter = 0; iter < kIters; ++iter) {
                S2Testing.Random.Reset(iter + 1);  // Easier to reproduce a specific case.
                S2Point src = S2Testing.RandomPoint();
                var simplifier = new S2PolylineSimplifier(src);
                S2Point dst = S2EdgeDistances.InterpolateAtDistance(
                    S1Angle.FromRadians(S2Testing.Random.RandDouble()),
                    src, S2Testing.RandomPoint());
                S2Point n = S2PointUtil.RobustCrossProd(src, dst).Normalized;

                // If bad_disc >= 0, then we make targeting fail for that disc.
                int kNumDiscs = 5;
                int bad_disc = S2Testing.Random.Uniform(2 * kNumDiscs) - kNumDiscs;
                for (int i = 0; i < kNumDiscs; ++i) {
                    double f = S2Testing.Random.RandDouble();
                    S2Point a = ((1 - f) * src + f * dst).Normalized;
                    S1Angle r = S1Angle.FromRadians(S2Testing.Random.RandDouble());
                    bool on_left = S2Testing.Random.OneIn(2);
                    S2Point x = S2EdgeDistances.InterpolateAtDistance(r, a, on_left ? n : -n);
                    // We grow the radius slightly if we want to target the disc and shrink
                    // it otherwise, *unless* we want targeting to fail for this disc, in
                    // which case these actions are reversed.
                    bool avoid = S2Testing.Random.OneIn(2);
                    bool grow_radius = (avoid == (i == bad_disc));
                    var radius = new S1ChordAngle(grow_radius ? r + kMaxError : r - kMaxError);
                    if (avoid) {
                        simplifier.AvoidDisc(x, radius, on_left);
                    } else {
                        simplifier.TargetDisc(x, radius);
                    }
                }
                // The result is true iff all the discraints were satisfiable.
                Assert.Equal(bad_disc < 0, simplifier.Extend(dst));
            }
        }

        private void CheckSimplify(string src, string dst,
                           string target, string avoid,
                           bool[] disc_on_left,
                           double radius_degrees, bool expected_result)
        {
            S1ChordAngle radius = new S1ChordAngle(S1Angle.FromDegrees(radius_degrees));
            var s = new S2PolylineSimplifier(S2TextFormat.MakePointOrDie(src));
            foreach (S2Point p in S2TextFormat.ParsePointsOrDie(target))
            {
                s.TargetDisc(p, radius);
            }
            int i = 0;
            foreach (S2Point p in S2TextFormat.ParsePointsOrDie(avoid))
            {
                s.AvoidDisc(p, radius, disc_on_left[i++]);
            }
            Assert.Equal(expected_result, s.Extend(S2TextFormat.MakePointOrDie(dst)));
        }
    }
}
