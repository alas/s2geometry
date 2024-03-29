// S2DistanceTarget represents a geometric object to which distances are
// measured.  For example, there are subtypes for measuring distances to a
// point, an edge, or to an S2ShapeIndex (an arbitrary collection of
// geometry).  S2DistanceTarget objects are provided for the benefit of
// classes that measure distances and/or find nearby geometry, such as
// S2ClosestEdgeQuery, S2ClosestPointQuery, and S2ClosestCellQuery.
//
// Implementations do *not* need to be thread-safe.  They may cache data or
// allocate temporary data structures in order to improve performance.  For
// this reason, S2DistanceTarget objects are typically passed as pointers
// rather than as references.
//
// The Distance template argument is used to represent distances.  Usually
// this type is a thin wrapper around S1ChordAngle, but another distance type
// may be substituted as long as it implements the API below.  This can be
// used to change the comparison function (e.g., to find the furthest edges
// from the target), to get more accuracy, or to measure non-spheroidal
// distances (e.g., using the WGS84 ellipsoid).
//
// The Distance concept is as follows:
//
// class Distance {
//  public:
//   // Default and copy constructors, assignment operator:
//   Distance();
//   Distance(Distance&);
//   Distance& operator=(Distance&);
//
//   // Factory methods:
//   static Distance Zero();      // Returns a zero distance.
//   static Distance Infinity();  // Larger than any valid distance.
//   static Distance Negative();  // Smaller than any valid distance.
//
//   // Comparison operators:
//   friend bool operator==(Distance x, Distance y);
//   friend bool operator<(Distance x, Distance y);
//
//   // Delta represents the positive difference between two distances.
//   // It is used together with operator-() to implement Options.max_error().
//   // Typically Distance.Delta is simply S1ChordAngle.
//   class Delta {
//    public:
//     Delta();
//     Delta(Delta&);
//     Delta& operator=(Delta&);
//     friend bool operator==(Delta x, Delta y);
//     static Delta Zero();
//   };
//
//   // Subtraction operator.  Note that the second argument represents a
//   // delta between two distances.  This distinction is important for
//   // classes that compute maximum distances (e.g., S2FurthestEdgeQuery).
//   friend Distance operator-(Distance x, Delta delta);
//
//   // Method that returns an upper bound on the S1ChordAngle corresponding
//   // to this Distance (needed to implement Options.max_distance
//   // efficiently).  For example, if Distance measures WGS84 ellipsoid
//   // distance then the corresponding angle needs to be 0.56% larger.
//   S1ChordAngle GetChordAngleBound();
// };

namespace S2Geometry;

public abstract class S2DistanceTarget<Distance> where Distance : IEquatable<Distance>, IComparable<Distance>, IDistance<Distance>
{
    // Returns an S2Cap that bounds the set of points whose distance to the
    // target is Distance.Zero().
    public abstract S2Cap GetCapBound();

    // If the distance to the point "p" is less than "min_dist", then updates "min_dist" and
    // returns true.  Otherwise returns false.
    public abstract bool UpdateMinDistance(S2Point p, ref Distance min_dist);

    // If the distance to the edge (v0, v1) is less than "min_dist", then
    // updates "min_dist" and returns true.  Otherwise returns false.
    public abstract bool UpdateMinDistance(S2Point v0, S2Point v1, ref Distance min_dist);

    // If the distance to the given S2Cell (including its interior) is less
    // than "min_dist", then updates "min_dist" and returns true.  Otherwise
    // returns false.
    public abstract bool UpdateMinDistance(S2Cell cell, ref Distance min_dist);

    // Finds all polygons in the given "query_index" that completely contain a
    // connected component of the target geometry.  (For example, if the
    // target consists of 10 points, this method finds polygons that contain
    // any of those 10 points.)  For each such polygon, "visitor" is called
    // with the S2Shape of the polygon along with a point of the target
    // geometry that is contained by that polygon.
    //
    // Optionally, any polygon that intersects the target geometry may also be
    // returned.  In other words, this method returns all polygons that
    // contain any connected component of the target, along with an arbitrary
    // subset of the polygons that intersect the target.
    //
    // For example, suppose that "query_index" contains two abutting polygons
    // A and B.  If the target consists of two points "a" contained by A and
    // "b" contained by B, then both A and B are returned.  But if the target
    // consists of the edge "ab", then any subset of {A, B} could be returned
    // (because both polygons intersect the target but neither one contains
    // the edge "ab").
    //
    // If "visitor" returns false, this method terminates early and returns
    // false as well.  Otherwise returns true.
    //
    // NOTE(ericv): This method exists only for the purpose of implementing
    // S2ClosestEdgeQuery::Options::include_interiors() efficiently.  Its API is
    // unlikely to be useful for other purposes.
    //
    // CAVEAT: Containing shapes may be visited more than once.
    public abstract bool VisitContainingShapes(S2ShapeIndex query_index, ShapeVisitor visitor);

    public delegate bool ShapeVisitor(S2Shape containing_shape, S2Point target_point);

    // Specifies that whenever one of the UpdateMinDistance() methods above
    // returns "true", the returned distance is allowed to be up to "max_error"
    // larger than the true minimum distance.  In other words, it gives this
    // target object permission to terminate its distance calculation as soon as
    // it has determined that (1) the minimum distance is less than "min_dist"
    // and (2) the best possible further improvement is less than "max_error".
    //
    // If the target takes advantage of "max_error" to optimize its distance
    // calculation, this method must return "true".  (Most target types can use
    // the default implementation which simply returns false.)
    internal virtual S1ChordAngle MaxError { get; set; }

    // The following method is provided as a convenience for classes that
    // compute distances to a collection of indexed geometry, such as
    // S2ClosestPointQuery, S2ClosestEdgeQuery, and S2ClosestCellQuery.  It
    // returns the maximum number of indexed objects for which it is faster to
    // compute the distance by brute force (e.g., by testing every edge) rather
    // than by using an index.  (The appropriate value is different for each
    // index type and can be estimated for a given (distance target, index type)
    // pair by running benchmarks.)
    //
    // By default this method returns -1, indicating that it is not implemented.
    public virtual int MaxBruteForceIndexSize => -1;
}
