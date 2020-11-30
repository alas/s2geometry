using System;
namespace S2Geometry
{
    // This class allows edges to be efficiently tested for intersection with a
    // given fixed edge AB.  It is especially efficient when testing for
    // intersection with an edge chain connecting vertices v0, v1, v2, ...
    //
    // Example usage:
    //
    //   void CountIntersections(S2Point a, S2Point b,
    //                           (S2Point, S2Point)[] edges) {
    //     int count = 0;
    //     S2EdgeCrosser crosser(out a, &b);
    //     foreach (var& edge in edges) {
    //       if (crosser.CrossingSign(out edge.Item1, &edge.Item2) >= 0) {
    //         ++count;
    //       }
    //     }
    //     return count;
    //   }
    //
    // This class expects that the client already has all the necessary vertices
    // stored in memory, so that this class can refer to them with pointers and
    // does not need to make its own copies.  If this is not the case (e.g., you
    // want to pass temporary objects as vertices), see S2CopyingEdgeCrosser.
    public class S2EdgeCrosser
    {
        // Default constructor; must be followed by a call to Init().
        public S2EdgeCrosser() {}

        // Convenience constructor that calls Init() with the given fixed edge AB.
        // The arguments "a" and "b" must point to values that persist for the
        // lifetime of the S2EdgeCrosser object (or until the next Init() call).
        public S2EdgeCrosser(S2Point a, S2Point b)
        {
            A = a;
            B = b;
            a_cross_b_ = A.CrossProd(B);
            have_tangents_ = false;
            c_ = S2Point.Empty;
            Assert.True(a.IsUnitLength);
            Assert.True(b.IsUnitLength);
        }

        // Accessors for the constructor arguments.
        public S2Point A { get; private set; }
        public S2Point B { get; private set; }

        // Initialize the S2EdgeCrosser with the given fixed edge AB.  The arguments
        // "a" and "b" must point to values that persist for the lifetime of the
        // S2EdgeCrosser object (or until the next Init() call).
        public void Init(S2Point a, S2Point b)
        {
            A = a;
            B = b;
            a_cross_b_ = a.CrossProd(B);
            have_tangents_ = false;
            c_ = S2Point.Empty;
        }

        // This function determines whether the edge AB intersects the edge CD.
        // Returns +1 if AB crosses CD at a point that is interior to both edges.
        // Returns  0 if any two vertices from different edges are the same.
        // Returns -1 otherwise.
        //
        // Note that if an edge is degenerate (A == B or C == D), the return value
        // is 0 if two vertices from different edges are the same and -1 otherwise.
        //
        // Properties of CrossingSign:
        //
        //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
        //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
        //  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
        //  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
        //
        // This function implements an exact, consistent perturbation model such
        // that no three points are ever considered to be collinear.  This means
        // that even if you have 4 points A, B, C, D that lie exactly in a line
        // (say, around the equator), C and D will be treated as being slightly to
        // one side or the other of AB.  This is done in a way such that the
        // results are always consistent (see S2Pred.Sign).
        //
        // Note that if you want to check an edge against a chain of other edges,
        // it is slightly more efficient to use the single-argument version of
        // CrossingSign below.
        //
        // The arguments must point to values that persist until the next call.
        public int CrossingSign(S2Point c, S2Point d)
        {
            if (c != c_) RestartAt(c);
            return CrossingSign(d);
        }

        // This method extends the concept of a "crossing" to the case where AB
        // and CD have a vertex in common.  The two edges may or may not cross,
        // according to the rules defined in VertexCrossing() below.  The rules
        // are designed so that point containment tests can be implemented simply
        // by counting edge crossings.  Similarly, determining whether one edge
        // chain crosses another edge chain can be implemented by counting.
        //
        // Returns true if CrossingSign(c, d) > 0, or AB and CD share a vertex
        // and VertexCrossing(a, b, c, d) returns true.
        //
        // The arguments must point to values that persist until the next call.
        public bool EdgeOrVertexCrossing(S2Point c, S2Point d)
        {
            if (c != c_) RestartAt(c);
            return EdgeOrVertexCrossing(d);
        }

        ///////////////////////// Edge Chain Methods ///////////////////////////
        //
        // You don't need to use these unless you're trying to squeeze out every
        // last drop of performance.  Essentially all you are saving is a test
        // whether the first vertex of the current edge is the same as the second
        // vertex of the previous edge.  Example usage:
        //
        //   S2Point[] chain;
        //   crosser.RestartAt(out chain[0]);
        //   for (int i = 1; i < chain.size(); ++i) {
        //     if (crosser.EdgeOrVertexCrossing(out chain[i])) { ++count; }
        //   }

        // Convenience constructor that uses AB as the fixed edge, and C as the
        // first vertex of the vertex chain (equivalent to calling RestartAt(c)).
        //
        // The arguments must point to values that persist until the next call.
        public S2EdgeCrosser(S2Point a, S2Point b, S2Point c)
        {
            A = a;
            B = b;
            a_cross_b_ = A.CrossProd(B);
            have_tangents_ = false;
            Assert.True(a.IsUnitLength);
            Assert.True(b.IsUnitLength);
            RestartAt(c);
        }

        // Call this method when your chain 'jumps' to a new place.
        // The argument must point to a value that persists until the next call.
        public void RestartAt(S2Point c)
        {
            Assert.True(c.IsUnitLength);
            c_ = c;
            acb_ = -S2Pred.TriageSign(A, B, c_, a_cross_b_);
        }

        // Like CrossingSign above, but uses the last vertex passed to one of
        // the crossing methods (or RestartAt) as the first vertex of the current
        // edge.
        //
        // The argument must point to a value that persists until the next call.
        public int CrossingSign(S2Point d)
        {
            Assert.True(d.IsUnitLength);
            // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
            // all be oriented the same way (CW or CCW).  We keep the orientation of ACB
            // as part of our state.  When each new point D arrives, we compute the
            // orientation of BDA and check whether it matches ACB.  This checks whether
            // the points C and D are on opposite sides of the great circle through AB.

            // Recall that TriageSign is invariant with respect to rotating its
            // arguments, i.e. ABD has the same orientation as BDA.
            int bda = S2Pred.TriageSign(A, B, d, a_cross_b_);
            if (acb_ == -bda && bda != 0)
            {
                // The most common case -- triangles have opposite orientations.  Save the
                // current vertex D as the next vertex C, and also save the orientation of
                // the new triangle ACB (which is opposite to the current triangle BDA).
                c_ = d;
                acb_ = -bda;
                return -1;
            }
            bda_ = bda;
            return CrossingSignInternal(d);
        }

        // Like EdgeOrVertexCrossing above, but uses the last vertex passed to one
        // of the crossing methods (or RestartAt) as the first vertex of the
        // current edge.
        //
        // The argument must point to a value that persists until the next call.
        public bool EdgeOrVertexCrossing(S2Point d)
        {
            // We need to copy c_ since it is clobbered by CrossingSign().
            var c = c_;
            int crossing = CrossingSign(d);
            if (crossing < 0) return false;
            if (crossing > 0) return true;
            return S2EdgeCrossings.VertexCrossing(A, B, c, d);
        }

        // Returns the last vertex of the current edge chain being tested, i.e. the
        // C vertex that will be used to construct the edge CD when one of the
        // methods above is called.
        public S2Point C => c_;
        // These functions handle the "slow path" of CrossingSign().
        private int CrossingSignInternal(S2Point d)
        {
            // Compute the actual result, and then save the current vertex D as the next
            // vertex C, and save the orientation of the next triangle ACB (which is
            // opposite to the current triangle BDA).
            int result = CrossingSignInternal2(d);
            c_ = d;
            acb_ = -bda_;
            return result;
        }
        private int CrossingSignInternal2(S2Point d)
        {
            // At this point, a very common situation is that A,B,C,D are four points on
            // a line such that AB does not overlap CD.  (For example, this happens when
            // a line or curve is sampled finely, or when geometry is constructed by
            // computing the union of S2CellIds.)  Most of the time, we can determine
            // that AB and CD do not intersect by computing the two outward-facing
            // tangents at A and B (parallel to AB) and testing whether AB and CD are on
            // opposite sides of the plane perpendicular to one of these tangents.  This
            // is moderately expensive but still much cheaper than S2Pred.ExpensiveSign.
            if (!have_tangents_) {
                S2Point norm = S2PointUtil.RobustCrossProd(A, B).Normalized;
                a_tangent_ = A.CrossProd(norm);
                b_tangent_ = norm.CrossProd(B);
                have_tangents_ = true;
            }

            if ((c_.DotProd(a_tangent_) > kError && d.DotProd(a_tangent_) > kError) ||
                (c_.DotProd(b_tangent_) > kError && d.DotProd(b_tangent_) > kError)) {
                return -1;
            }

            // Otherwise, eliminate the cases where two vertices from different edges
            // are equal.  (These cases could be handled in the code below, but we would
            // rather avoid calling ExpensiveSign whenever possible.)
            if (A == c_ || A == d || B == c_ || B == d) return 0;

            // Eliminate cases where an input edge is degenerate.  (Note that in most
            // cases, if CD is degenerate then this method is not even called because
            // acb_ and bda have different signs.)
            if (A == B || c_ == d) return -1;

            // Otherwise it's time to break out the big guns.
            if (acb_ == 0) acb_ = -S2Pred.ExpensiveSign(A, B, c_);
            Assert.True(acb_ != 0);
            if (bda_ == 0) bda_ = S2Pred.ExpensiveSign(A, B, d);
            Assert.True(bda_ != 0);
            if (bda_ != acb_) return -1;

            var c_cross_d = c_.CrossProd(d);
            int cbd = -S2Pred.Sign(c_, d, B, c_cross_d);
            Assert.True(cbd != 0);
            if (cbd != acb_) return -1;
            int dac = S2Pred.Sign(c_, d, A, c_cross_d);
            Assert.True(dac != 0);
            return (dac != acb_) ? -1 : 1;
        }

        // The error in RobustCrossProd() is insignificant.  The maximum error in
        // the call to CrossProd() (i.e., the maximum norm of the error vector) is
        // (0.5 + 1/Math.Sqrt(3)) * S2Constants.DoubleEpsilon.  The maximum error in each call to
        // DotProd() below is S2Constants.DoubleEpsilon.  (There is also a small relative error
        // term that is insignificant because we are comparing the result against a
        // constant that is very close to zero.)
        private static readonly double kError = (1.5 + 1 / Math.Sqrt(3)) * S2Constants.DoubleEpsilon;
       

        // Used internally by S2CopyingEdgeCrosser.  Updates "c_" only.
        internal void SetC(S2Point c) { c_ = c; }

        private S2Point a_cross_b_;

        // To reduce the number of calls to S2Pred.ExpensiveSign(), we compute an
        // outward-facing tangent at A and B if necessary.  If the plane
        // perpendicular to one of these tangents separates AB from CD (i.e., one
        // edge on each side) then there is no intersection.
        private bool have_tangents_;  // True if the tangents have been computed.
        private S2Point a_tangent_;   // Outward-facing tangent at A.
        private S2Point b_tangent_;   // Outward-facing tangent at B.

        // The fields below are updated for each vertex in the chain.
        internal S2Point c_;       // Previous vertex in the vertex chain.
        private int acb_;                // The orientation of triangle ACB.

        // The field below is a temporary used by CrossingSignInternal().
        private int bda_;                // The orientation of triangle BDA.
    }

    // S2CopyingEdgeCrosser is exactly like S2EdgeCrosser, except that it makes its
    // own copy of all arguments so that they do not need to persist between
    // calls.  This is less efficient, but makes it possible to use points that
    // are generated on demand and cannot conveniently be stored by the client.
    public class S2CopyingEdgeCrosser
    {
        // These methods are all exactly like S2EdgeCrosser, except that the
        // arguments can be temporaries.
        public S2CopyingEdgeCrosser() {}
        public S2CopyingEdgeCrosser(S2Point a, S2Point b)
        {
            A = a;
            B = b;
            C = S2Point.Empty;
            crosser_ = new S2EdgeCrosser(A, B);
        }
        public S2Point A { get; private set; }
        public S2Point B { get; private set; }
        public S2Point C { get; private set; }
        public void Init(S2Point a, S2Point b)
        {
            A = a;
            B = b;
            C = S2Point.Empty;
            crosser_.Init(A, B);
        }
        public int CrossingSign(S2Point c, S2Point d)
        {
            if (c != C || crosser_.c_ == S2Point.Empty) RestartAt(c);
            return CrossingSign(d);
        }
        public bool EdgeOrVertexCrossing(S2Point c, S2Point d)
        {
            if (c != C || crosser_.c_ == S2Point.Empty) RestartAt(c);
            return EdgeOrVertexCrossing(d);
        }
        public S2CopyingEdgeCrosser(S2Point a, S2Point b, S2Point c)
        {
            A = a;
            B = b;
            C = c;
            crosser_ = new S2EdgeCrosser(A, B, c);
        }
        public void RestartAt(S2Point c)
        {
            C = c;
            crosser_.RestartAt(C);
        }
        public int CrossingSign(S2Point d)
        {
            int result = crosser_.CrossingSign(d);
            C = d;
            crosser_.SetC(C);
            return result;
        }
        public bool EdgeOrVertexCrossing(S2Point d)
        {
            bool result = crosser_.EdgeOrVertexCrossing(d);
            C = d;
            crosser_.SetC(C);
            return result;
        }

        // TODO(ericv): It would be more efficient to implement S2CopyingEdgeCrosser
        // directly rather than as a wrapper around S2EdgeCrosser.
        private readonly S2EdgeCrosser crosser_;
    }
}
