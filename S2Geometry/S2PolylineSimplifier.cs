// This is a helper class for simplifying polylines.  It allows you to compute
// a maximal edge that intersects a sequence of discs, and that optionally
// avoids a different sequence of discs.  The results are conservative in that
// the edge is guaranteed to intersect or avoid the specified discs using
// exact arithmetic (see s2predicates.h).
//
// Note that S2Builder can also simplify polylines and supports more features
// (e.g., snapping to S2CellId centers), so it is only recommended to use this
// class if S2Builder does not meet your needs.
//
// Here is a simple example showing how to simplify a polyline into a sequence
// of edges that stay within "max_error" of the original edges:
//
//   S2Point[] v = { ... };
//   S2PolylineSimplifier simplifier;
//   simplifier.Init(v[0]);
//   for (int i = 1; i < v.size(); ++i) {
//     if (!simplifier.Extend(v[i])) {
//       OutputEdge(simplifier.src(), v[i-1]);
//       simplifier.Init(v[i-1]);
//     }
//     simplifier.TargetDisc(v[i], max_error);
//   }
//   OutputEdge(simplifer.src(), v.Last());
//
// Note that the points targeted by TargetDisc do not need to be the same as
// the candidate endpoints passed to Extend.  So for example, you could target
// the original vertices of a polyline, but only consider endpoints that are
// snapped to E7 coordinates or S2CellId centers.
//
// Please be aware that this class works by maintaining a range of acceptable
// angles (bearings) from the start vertex to the hypothetical destination
// vertex.  It does not keep track of distances to any of the discs to be
// targeted or avoided.  Therefore to use this class correctly,raints
// should be added in increasing order of distance.  (The actual requirement
// is slightly weaker than this, which is why it is not enforced, but
// basically you should only call TargetDisc() and AvoidDisc() with arguments
// that you want torain the immediately following call to Extend().)

namespace S2Geometry;

public class S2PolylineSimplifier
{
    /// <summary>
    /// Source vertex of the output edge.
    /// </summary>
    public S2Point Source { get; }              // Output edge source vertex.
    private readonly S2Point x_dir_, y_dir_;    // Orthonormal frame for mapping vectors to angles.
    private S1Interval window_;                 // Allowable range of angles for the output edge.

    // We store the discs to avoid individually until TargetDisc() is first
    // called with a disc that does not contain the source vertex.  At that time
    // all such discs are processed by using them to constrain "window_", and
    // this vector is cleared.
    private readonly List<RangeToAvoid> ranges_to_avoid_ = [];

    // Starts a new simplified edge at "src".
    public S2PolylineSimplifier(S2Point src)
    {
        Source = src;
        window_ = S1Interval.Full;
        ranges_to_avoid_.Clear();

        // Precompute basis vectors for the tangent space at "src".  This is similar
        // to GetFrame() except that we don't normalize the vectors.  As it turns
        // out, the two basis vectors below have the same magnitude (up to the
        // length error in S2Point.Normalize).

        // Find the index of the component whose magnitude is smallest.
        var tmp = Source.Fabs();
        int i = tmp[0] < tmp[1]
            ? (tmp[0] < tmp[2] ? 0 : 2)
            : (tmp[1] < tmp[2] ? 1 : 2);

        // We define the "y" basis vector as the cross product of "src" and the
        // basis vector for axis "i".  Let "j" and "k" be the indices of the other
        // two components in cyclic order.
        int j = i == 2 ? 0 : i + 1, k = i == 0 ? 2 : i - 1;
        var np = new double[3];
        np[i] = 0;
        np[j] = Source[k];
        np[k] = -Source[j];
        y_dir_ = new S2Point(np);

        // Compute the cross product of "y_dir" and "src".  We write out the cross
        // product here mainly for documentation purposes; it also happens to save a
        // few multiplies because unfortunately the optimizer does *not* get rid of
        // multiplies by zero (since these multiplies propagate NaN, for example).
        np[i] = Source[j] * Source[j] + Source[k] * Source[k];
        np[j] = -Source[j] * Source[i];
        np[k] = -Source[k] * Source[i];
        x_dir_ = new S2Point(np);
    }

    // Returns true if the edge (src, dst) satisfies all of the targeting
    // requirements so far.  Returns false if the edge would be longer than
    // 90 degrees (such edges are not supported).
    public bool Extend(S2Point dst)
    {
        // We limit the maximum edge length to 90 degrees in order to simplify the
        // error bounds.  (The error gets arbitrarily large as the edge length
        // approaches 180 degrees.)
        if (new S1ChordAngle(Source, dst) > S1ChordAngle.Right) return false;

        // Otherwise check whether this vertex is in the acceptable angle range.
        double dir = GetDirection(dst);
        if (!window_.Contains(dir)) return false;

        // Also check any angles ranges to avoid that have not been processed yet.
        foreach (var range in ranges_to_avoid_)
        {
            if (range.Interval.Contains(dir)) return false;
        }
        return true;
    }

    // Requires that the output edge must pass through the given disc.
    public bool TargetDisc(S2Point point, S1ChordAngle radius)
    {
        // Shrink the target interval by the maximum error from all sources.  This
        // guarantees that the output edge will intersect the given disc.
        double semiwidth = GetSemiwidth(point, radius, -1 /*round down*/);
        if (semiwidth >= Math.PI)
        {
            // The target disc contains "src", so there is nothing to do.
            return true;
        }
        if (semiwidth < 0)
        {
            window_ = S1Interval.Empty;
            return false;
        }
        // Otherwise compute the angle interval corresponding to the target disc and
        // intersect it with the current window.
        double center = GetDirection(point);
        S1Interval target = S1Interval.FromPoint(center).Expanded(semiwidth);
        window_ = window_.Intersection(target);

        // If there are any angle ranges to avoid, they can be processed now.
        foreach (var range in ranges_to_avoid_)
        {
            AvoidRange(range.Interval, range.OnLeft);
        }
        ranges_to_avoid_.Clear();

        return !window_.IsEmpty();
    }

    // Requires that the output edge must avoid the given disc.  "disc_on_left"
    // specifies whether the disc must be to the left or right of the output
    // edge AB.  (This feature allows the simplified edge to preserve the
    // topology of the original polyline with respect to other nearby points.)
    //
    // More precisely, let AB be the output edge, P be the center of the disc,
    // and r be its radius.  Then this method ensures that
    //
    //   (1) Distance(AB, P) > r, and
    //   (2) if DotProd(AB, AP) > 0, then Sign(ABP) > 0 iff disc_on_left is true.
    //
    // The second condition says that "disc_on_left" has an effect if and only
    // if P is not behind the source vertex A with respect to the direction AB.
    //
    // If your input is a polyline, you can compute "disc_on_left" as follows.
    // Let the polyline be ABCDE and assume that it already avoids a set of
    // points X_i.  Suppose that you have aleady added ABC to the simplifier, and
    // now want to extend the edge chain to D.  First find the X_i that are near
    // the edge CD, then discard the ones such that AX_i <= AC or AX_i >= AD
    // (since these points have either already been considered or aren't
    // relevant yet).  Now X_i is to the left of the polyline if and only if
    // s2pred::OrderedCCW(A, D, X_i, C) (in other words, if X_i is to the left of
    // the angle wedge ACD).  Note that simply testing s2pred::Sign(C, D, X_i)
    // or s2pred::Sign(A, D, X_i) does not handle all cases correctly.
    public bool AvoidDisc(S2Point point, S1ChordAngle radius, bool disc_on_left)
    {
        // Expand the interval by the maximum error from all sources.  This
        // guarantees that the final output edge will avoid the given disc.
        double semiwidth = GetSemiwidth(point, radius, 1 /*round up*/);
        if (semiwidth >= Math.PI)
        {
            // The disc to avoid contains "src", so it can't be avoided.
            window_ = S1Interval.Empty;
            return false;
        }
        // Compute the disallowed range of angles: the angle subtended by the disc
        // on one side, and 90 degrees on the other (to satisfy "disc_on_left").
        double center = GetDirection(point);
        double dleft = disc_on_left ? S2.M_PI_2 : semiwidth;
        double dright = disc_on_left ? semiwidth : S2.M_PI_2;
        S1Interval avoid_interval = new(Math.IEEERemainder(center - dright, 2 * S2.M_PI),
                            Math.IEEERemainder(center + dleft, 2 * S2.M_PI));

        if (window_.IsFull())
        {
            // Discs to avoid can't be processed until window_ is reduced to at most
            // 180 degrees by a call to TargetDisc().  Save it for later.
            ranges_to_avoid_.Add(new RangeToAvoid(avoid_interval, disc_on_left));
            return true;
        }
        AvoidRange(avoid_interval, disc_on_left);
        return !window_.IsEmpty();
    }

    private void AvoidRange(S1Interval avoid_interval, bool disc_on_left)
    {
        // If "avoid_interval" is a proper subset of "window_", then in theory the
        // result should be two intervals.  One interval points towards the given
        // disc and pass on the correct side of it, while the other interval points
        // away from the disc.  However the latter interval never contains an
        // acceptable output edge direction (as long as this class is being used
        // correctly) and can be safely ignored.  This is true because (1) "window_"
        // is not full, which means that it contains at least one vertex of the input
        // polyline and is at most 180 degrees in length, and (2) "disc_on_left" is
        // computed with respect to the next edge of the input polyline, which means
        // that the next input vertex is either inside "avoid_interval" or somewhere
        // in the 180 degrees to its right/left according to "disc_on_left", which
        // means that it cannot be contained by the subinterval that we ignore.
        MyDebug.Assert(!window_.IsFull());
        if (window_.Contains(avoid_interval))
        {
            if (disc_on_left)
            {
                window_ = new S1Interval(window_.Lo, avoid_interval.Lo);
            }
            else 
            {
                window_ = new S1Interval(avoid_interval.Hi, window_.Hi);
            }
        }
        else
        {
            window_ = window_.Intersection(avoid_interval.Complement());
        }
    }

    private double GetDirection(S2Point p)
    {
        return Math.Atan2(p.DotProd(y_dir_), p.DotProd(x_dir_));
    }

    // Computes half the angle in radians subtended from the source vertex by a
    // disc of radius "r" centered at "p", rounding the result conservatively up
    // or down according to whether round_direction is +1 or -1.  (So for example,
    // if round_direction == +1 then the return value is an upper bound on the
    // true result.)
    private double GetSemiwidth(S2Point p, S1ChordAngle r, int round_direction)
    {
        const double DBL_ERR = 0.5 * S2.DoubleEpsilon;

        // Using spherical trigonometry,
        //
        //   sin(semiwidth) = sin(r) / sin(a)
        //
        // where "a" is the angle between "src" and "p".  Rather than measuring
        // these angles, instead we measure the squared chord lengths through the
        // interior of the sphere (i.e., Cartersian distance).  Letting "r2" be the
        // squared chord distance corresponding to "r", and "a2" be the squared
        // chord distance corresponding to "a", we use the relationships
        //
        //    sin^2(r) = r2 (1 - r2 / 4)
        //    sin^2(a) = d2 (1 - d2 / 4)
        //
        // which follow from the fact that r2 = (2 * sin(r / 2)) ^ 2, etc.

        // "a2" has a relative error up to 5 * DBL_ERR, plus an absolute error of up
        // to 64 * DBL_ERR^2 (because "src" and "p" may differ from unit length by
        // up to 4 * DBL_ERR).  We can correct for the relative error later, but for
        // the absolute error we use "round_direction" to account for it now.
        var r2 = r.Length2;
        var a2 = new S1ChordAngle(Source, p).Length2;
        a2 -= 64 * DBL_ERR * DBL_ERR * round_direction;
        if (a2 <= r2) return Math.PI;  // The given disc contains "src".

        double sin2_r = r2 * (1 - 0.25 * r2);
        double sin2_a = a2 * (1 - 0.25 * a2);
        double semiwidth = Math.Asin(Math.Sqrt(sin2_r / sin2_a));

        // We compute bounds on the errors from all sources:
        //
        //   - The call to GetSemiwidth (this call).
        //   - The call to GetDirection that computes the center of the interval.
        //   - The call to GetDirection in Extend that tests whether a given point
        //     is an acceptable destination vertex.
        //
        // Summary of the errors in GetDirection:
        //
        // y_dir_ has no error.
        //
        // x_dir_ has a relative error of DBL_ERR in two components, a relative
        // error of 2 * DBL_ERR in the other component, plus an overall relative
        // length error of 4 * DBL_ERR (compared to y_dir_) because "src" is assumed
        // to be normalized only to within the tolerances of S2Point.Normalize().
        //
        // p.DotProd(y_dir_) has a relative error of 1.5 * DBL_ERR and an
        // absolute error of 1.5 * DBL_ERR * y_dir_.Norm.
        //
        // p.DotProd(x_dir_) has a relative error of 5.5 * DBL_ERR and an absolute
        // error of 3.5 * DBL_ERR * y_dir_.Norm (noting that x_dir_ and y_dir_
        // have the same length to within a relative error of 4 * DBL_ERR).
        //
        // It's possible to show by taking derivatives that these errors can affect
        // the angle atan2(y, x) by up 7.093 * DBL_ERR radians.  Rounding up and
        // including the call to atan2 gives a final error bound of 10 * DBL_ERR.
        //
        // Summary of the errors in GetSemiwidth:
        //
        // The distance a2 has a relative error of 5 * DBL_ERR plus an absolute
        // error of 64 * DBL_ERR^2 because the points "src" and "p" may differ from
        // unit length (by up to 4 * DBL_ERR).  We have already accounted for the
        // absolute error above, leaving only the relative error.
        //
        // sin2_r has a relative error of 2 * DBL_ERR.
        //
        // sin2_a has a relative error of 12 * DBL_ERR assuming that a2 <= 2,
        // i.e. distance(src, p) <= 90 degrees.  (The relative error gets
        // arbitrarily larger as this distance approaches 180 degrees.)
        //
        // semiwidth has a relative error of 17 * DBL_ERR.
        //
        // Finally, (center +/- semiwidth) has a rounding error of up to 4 * DBL_ERR
        // because in theory, the result magnitude may be as large as 1.5 * Math.PI
        // which is larger than 4.0.  This gives a total error of:
        double error = (2 * 10 + 4) * DBL_ERR + 17 * DBL_ERR * semiwidth;
        return semiwidth + round_direction * error;
    }

    // Unfortunately, the discs to avoid cannot be processed until the direction
    // of the output edge is constrained to lie within an S1Interval of at most
    // 180 degrees.  This happens only when the first target disc is added that
    // does not contain the source vertex.  Until that time we simply store all
    // the discs as ranges of directions to avoid.
    private readonly record struct RangeToAvoid(
                    S1Interval Interval,    // Range of directions to avoid.
                    bool OnLeft);          // Is this disc to the left of the output edge?
}
