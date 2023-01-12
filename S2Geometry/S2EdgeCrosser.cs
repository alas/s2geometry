// This file defines two classes S2EdgeCrosser and S2CopyingEdgeCrosser that
// allow edges to be efficiently tested for intersection with a given fixed
// edge AB.  They are especially efficient when testing for intersection with
// an edge chain connecting vertices v0, v1, v2, ...
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
// The example above uses S2EdgeCrosser, which requires that the client
// already has all the necessary vertices stored in memory so that this class
// can refer to them with pointers and does not need to make its own copies.
// If this is not the case (e.g., you want to pass temporary objects as
// vertices) then you should use S2CopyingEdgeCrosser, which has exactly the
// same API except that vertices are passed by const reference and do not need
// to persist.
//
// The class below is instantiated twice:
//
//  - For S2EdgeCrosser, ArgType is (const S2Point*) and all points must be
//    stored persistently by the client.
//
//  - For S2CopyingEdgeCrosser, ArgType is (const S2Point&) and points may
//    be temporary objects (since this class makes its own copies).
//
// Note that S2EdgeCrosser is 5-10% faster in real applications when its
// requirements can be met.  (Also note that simple benchmarks are not
// sufficient to measure this performance difference; it seems to have
// something to do with cache performance.)

// Explicitly instantiate the classes we need so that the methods above can be
// omitted from the .h file (and to reduce compilation time).
//using S2EdgeCrosserBase_PointerRep = S2EdgeCrosserBase<S2Point_PointerRep>;
//using S2EdgeCrosserBase_ValueRep = S2EdgeCrosserBase<S2Point_ValueRep>;

// S2EdgeCrosser implements the API above by using (const S2Point *) as the
// argument type and requiring that all points are stored persistently by the
// client.  If this is not the case, use S2CopyingEdgeCrosser (below).
//using S2EdgeCrosser = S2EdgeCrosserBase<S2::internal::S2Point_PointerRep>;

// S2CopyingEdgeCrosser is exactly like S2EdgeCrosser except that it makes its
// own copy of all arguments so that they do not need to persist between
// calls.  This is less efficient, but makes it possible to use points that
// are generated on demand and cannot conveniently be stored by the client.
//using S2CopyingEdgeCrosser = S2EdgeCrosserBase<S2::internal::S2Point_ValueRep>;

// Let the compiler know that these classes are explicitly instantiated in the
// .cc file; this helps to reduce compilation time.
//extern template class S2EdgeCrosserBase<S2::internal::S2Point_PointerRep>;
//extern template class S2EdgeCrosserBase<S2::internal::S2Point_ValueRep>;

global using S2EdgeCrosser = S2Geometry.S2EdgeCrosserBase;
global using S2CopyingEdgeCrosser = S2Geometry.S2EdgeCrosserBase;

namespace S2Geometry;

public class S2EdgeCrosserBase
{
    // Accessors for the constructor arguments.
    //
    // const S2Point* S2EdgeCrosser::a();
    // const S2Point* S2EdgeCrosser::b();
    // const S2Point& S2CopyingEdgeCrosser::a();
    // const S2Point& S2CopyingEdgeCrosser::b();
    public S2Point A { get; private set; }
    public S2Point B { get; private set; }

    // Default constructor; must be followed by a call to Init().
    public S2EdgeCrosserBase() { }

    // Convenience constructor that calls Init() with the given fixed edge AB.
    //
    // For S2EdgeCrosser (only), the arguments "a" and "b" must point to values
    // that persist for the lifetime of the S2EdgeCrosser object.
    //
    // S2EdgeCrosser(const S2Point* a, const S2Point* b);
    // S2CopyingEdgeCrosser(const S2Point& a, const S2Point& b);
    public S2EdgeCrosserBase(S2Point a, S2Point b)
    {
        A = a;
        B = b;
        a_cross_b_ = A.CrossProd(B);
        have_tangents_ = false;
        C = S2Point.Empty;
        MyDebug.Assert(a.IsUnitLength());
        MyDebug.Assert(b.IsUnitLength());
    }

    // Initialize the object with the given fixed edge AB.
    //
    // For S2EdgeCrosser (only), the arguments "a" and "b" must point to values
    // that persist for the lifetime of the S2EdgeCrosser object.
    //
    // void S2EdgeCrosser::Init(const S2Point* a, const S2Point* b);
    // void S2CopyingEdgeCrosser::Init(const S2Point& a, const S2Point& b);
    public void Init(S2Point a, S2Point b)
    {
        A = a;
        B = b;
        a_cross_b_ = a.CrossProd(B);
        have_tangents_ = false;
        C = S2Point.Empty;
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
    // For S2EdgeCrosser (only), the arguments must point to values that persist
    // until the next call.
    //
    // int S2EdgeCrosser::CrossingSign(const S2Point* c, const S2Point* d);
    // int S2CopyingEdgeCrosser::CrossingSign(const S2Point& c, const S2Point& d);
    public int CrossingSign(S2Point c, S2Point d)
    {
        if (c != C) RestartAt(c);
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
    // For S2EdgeCrosser (only), the arguments must point to values that persist
    // until the next call.
    //
    // bool S2EdgeCrosser::EdgeOrVertexCrossing(const S2Point* c,
    //                                          const S2Point* d);
    // bool S2CopyingEdgeCrosser::EdgeOrVertexCrossing(const S2Point& c,
    //                                                 const S2Point& d);
    public bool EdgeOrVertexCrossing(S2Point c, S2Point d)
    {
        if (c != C) RestartAt(c);
        return EdgeOrVertexCrossing(d);
    }

    // Like EdgeOrVertexCrossing() but returns -1 if AB crosses CD from left to
    // right, +1 if AB crosses CD from right to left, and 0 otherwise.  This
    // implies that if CD bounds some region according to the "interior is on
    // the left" rule, this function returns -1 when AB exits the region and +1
    // when AB enters.
    //
    // This method allows computing the change in winding number from point A to
    // point B by summing the signed edge crossings of AB with the edges of the
    // loop(s) used to define the winding number.
    //
    // The arguments must point to values that persist until the next call.
    //
    // int S2EdgeCrosser::SignedEdgeOrVertexCrossing(const S2Point* c,
    //                                               const S2Point* d);
    // int S2CopyingEdgeCrosser::SignedEdgeOrVertexCrossing(const S2Point& c,
    //                                                      const S2Point& d);
    public int SignedEdgeOrVertexCrossing(S2Point c, S2Point d)
    {
        if (c != C) RestartAt(c);
        return SignedEdgeOrVertexCrossing(d);
    }

    // If the preceding call to CrossingSign() returned +1 (indicating that the
    // edge crosses the edge CD), this method returns -1 if AB crossed CD from
    // left to right and +1 if AB crossed CD from right to left.  Otherwise its
    // return value is undefined.
    //
    // When AB crosses CD, the crossing sign is Sign(ABC).  S2EdgeCrosser doesn't
    // store this, but it does store the sign of the *next* triangle ACB.  These
    // two values happen to be the same.
    public int last_interior_crossing_sign() => acb_;
       

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
    // For S2EdgeCrosser (only), the arguments must point to values that persist
    // until the next call.
    //
    // S2EdgeCrosser(S2Point const* a, S2Point const* b, S2Point const* c);
    // S2CopyingEdgeCrosser(S2Point const& a, S2Point const& b, S2Point const& c);
    public S2EdgeCrosserBase(S2Point a, S2Point b, S2Point c)
    {
        A = a;
        B = b;
        a_cross_b_ = A.CrossProd(B);
        have_tangents_ = false;
        MyDebug.Assert(a.IsUnitLength());
        MyDebug.Assert(b.IsUnitLength());
        RestartAt(c);
    }

    // Call this method when your chain 'jumps' to a new place.
    //
    // For S2EdgeCrosser (only), the argument must point to a value that
    // persists until the next call.
    //
    // void S2EdgeCrosser::RestartAt(S2Point const* c);
    // void S2CopyingEdgeCrosser::RestartAt(S2Point const& c);
    public void RestartAt(S2Point c)
    {
        MyDebug.Assert(c.IsUnitLength());
        C = c;
        acb_ = -S2Pred.TriageSign(A, B, C, a_cross_b_);
    }

    // Like CrossingSign above, but uses the last vertex passed to one of the
    // crossing methods (or RestartAt) as the first vertex of the current edge.
    //
    // For S2EdgeCrosser (only), the argument must point to a value that
    // persists until the next call.
    //
    // int S2EdgeCrosser::CrossingSign(S2Point const* d);
    // int S2CopyingEdgeCrosser::CrossingSign(S2Point const& d);
    public int CrossingSign(S2Point d)
    {
        MyDebug.Assert(d.IsUnitLength());
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
            C = d;
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
    // For S2EdgeCrosser (only), the argument must point to a value that
    // persists until the next call.
    //
    // bool S2EdgeCrosser::EdgeOrVertexCrossing(S2Point const* d);
    // bool S2CopyingEdgeCrosser::EdgeOrVertexCrossing(S2Point const& d);
    public bool EdgeOrVertexCrossing(S2Point d)
    {
        // We need to copy c_ since it is clobbered by CrossingSign().
        var c = C;
        int crossing = CrossingSign(d);
        if (crossing < 0) return false;
        if (crossing > 0) return true;
        return S2.VertexCrossing(A, B, c, d);
    }

    // Like EdgeOrVertexCrossing above, but uses the last vertex passed to one
    // of the crossing methods (or RestartAt) as the first vertex of the
    // current edge.
    //
    // For S2EdgeCrosser (only), the argument must point to a value that
    // persists until the next call.
    //
    // int S2EdgeCrosser::SignedEdgeOrVertexCrossing(S2Point const* d);
    // int S2CopyingEdgeCrosser::SignedEdgeOrVertexCrossing(S2Point const& d);
    public int SignedEdgeOrVertexCrossing(S2Point d)
    {
        // We need to copy c_ since it is clobbered by CrossingSign().
        var c = C;
        var crossing = CrossingSign(d);
        if (crossing < 0) return 0;
        if (crossing > 0) return last_interior_crossing_sign();
        return S2.SignedVertexCrossing(A, B, c, d);
    }

    // Returns the last vertex of the current edge chain being tested, i.e.
    // the C vertex that will be used to construct the edge CD when one of the
    // methods above is called.
    //
    // const S2Point* S2EdgeCrosser::c();
    // const S2Point& S2CopyingEdgeCrosser::c();
    public S2Point C { get; private set; } // Previous vertex in the vertex chain.
    // These functions handle the "slow path" of CrossingSign().
    private int CrossingSignInternal(S2Point d)
    {
        // Compute the actual result, and then save the current vertex D as the next
        // vertex C, and save the orientation of the next triangle ACB (which is
        // opposite to the current triangle BDA).
        int result = CrossingSignInternal2(d);
        C = d;
        acb_ = -bda_;
        return result;
    }
    private int CrossingSignInternal2(S2Point d)
    {
        // At this point it is still very likely that CD does not cross AB.  Two
        // common situations are (1) CD crosses the great circle through AB but does
        // not cross AB itself, or (2) A,B,C,D are four points on a line such that
        // AB does not overlap CD.  For example, the latter happens when a line or
        // curve is sampled finely, or when geometry is constructed by computing the
        // union of S2CellIds.
        //
        // Most of the time, we can determine that AB and CD do not intersect by
        // computing the two outward-facing tangents at A and B (parallel to AB) and
        // testing whether AB and CD are on opposite sides of the plane perpendicular
        // to one of these tangents.  This is somewhat expensive but still much
        // cheaper than s2pred::ExpensiveSign.
        if (!have_tangents_)
        {
            S2Point norm = S2.RobustCrossProd(A, B);
            a_tangent_ = A.CrossProd(norm);
            b_tangent_ = norm.CrossProd(B);
            have_tangents_ = true;
        }

        if ((C.DotProd(a_tangent_) > kError && d.DotProd(a_tangent_) > kError) ||
            (C.DotProd(b_tangent_) > kError && d.DotProd(b_tangent_) > kError))
        {
            return -1;
        }

        // Otherwise, eliminate the cases where two vertices from different edges
        // are equal.  (These cases could be handled in the code below, but we would
        // rather avoid calling ExpensiveSign whenever possible.)
        if (A == C || A == d || B == C || B == d) return 0;

        // Eliminate cases where an input edge is degenerate.  (Note that in most
        // cases, if CD is degenerate then this method is not even called because
        // acb_ and bda have different signs.)
        if (A == B || C == d) return -1;

        // Otherwise it's time to break out the big guns.
        if (acb_ == 0) acb_ = -S2Pred.ExpensiveSign(A, B, C);
        MyDebug.Assert(acb_ != 0);
        if (bda_ == 0) bda_ = S2Pred.ExpensiveSign(A, B, d);
        MyDebug.Assert(bda_ != 0);
        if (bda_ != acb_) return -1;

        var c_cross_d = C.CrossProd(d);
        int cbd = -S2Pred.Sign(C, d, B, c_cross_d);
        MyDebug.Assert(cbd != 0);
        if (cbd != acb_) return -1;
        int dac = S2Pred.Sign(C, d, A, c_cross_d);
        MyDebug.Assert(dac != 0);
        return (dac != acb_) ? -1 : 1;
    }

    // The error in RobustCrossProd() is insignificant.  The maximum error in
    // the call to CrossProd() (i.e., the maximum norm of the error vector) is
    // (0.5 + 1/Math.Sqrt(3)) * S2Constants.DoubleEpsilon.  The maximum error in each call to
    // DotProd() below is S2Constants.DoubleEpsilon.  (There is also a small relative error
    // term that is insignificant because we are comparing the result against a
    // constant that is very close to zero.)
    private static readonly double kError = (1.5 + 1 / Math.Sqrt(3)) * S2.DoubleEpsilon;

    private S2Point a_cross_b_;

    // To reduce the number of calls to S2Pred.ExpensiveSign(), we compute an
    // outward-facing tangent at A and B if necessary.  If the plane
    // perpendicular to one of these tangents separates AB from CD (i.e., one
    // edge on each side) then there is no intersection.
    private bool have_tangents_;    // True if the tangents have been computed.
    private S2Point a_tangent_;     // Outward-facing tangent at A.
    private S2Point b_tangent_;     // Outward-facing tangent at B.

    // The fields below are updated for each vertex in the chain.  acb_ is
    // initialized to avoid undefined behavior in the case where the edge chain
    // starts with the invalid point (0, 0, 0).
    private int acb_ = 0;                // The orientation of triangle ACB.

    // The field below is a temporary used by CrossingSignInternal().
    private int bda_;               // The orientation of triangle BDA.
}
