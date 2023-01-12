// The bulk of this file consists of "tests" that attempt to construct worst
// cases for the various constants used in S2CellIdSnapFunction and
// IntLatLngSnapFunction implementations.  For all of theseants I have
// done hand analysis of the planar configurations, but sometimes the
// spherical case is slightly better or worse because of the spherical
// distortion.

namespace S2Geometry;

public class S2BuilderUtil_SnapFunctionsTests
{
    private static readonly S2CellId kSearchRootId = S2CellId.FromFace(0);
    private static readonly S2CellId kSearchFocusId = S2CellId.FromFace(0).Child(3);
    private readonly ITestOutputHelper _logger;

    public S2BuilderUtil_SnapFunctionsTests(ITestOutputHelper logger) { _logger = logger; }

    [Fact]
    internal void Test_S2CellIdSnapFunction_LevelToFromSnapRadius()
    {
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            S1Angle radius = S2CellIdSnapFunction.MinSnapRadiusForLevel(level);
            Assert.Equal(level, S2CellIdSnapFunction.LevelForMaxSnapRadius(radius));
            Assert.Equal(Math.Min(level + 1, S2.kMaxCellLevel),
                      S2CellIdSnapFunction.LevelForMaxSnapRadius(0.999 * radius));
        }
        Assert.Equal(0,
                  S2CellIdSnapFunction.LevelForMaxSnapRadius(
                      S1Angle.FromRadians(5)));
        Assert.Equal(S2.kMaxCellLevel,
                  S2CellIdSnapFunction.LevelForMaxSnapRadius(
                      S1Angle.FromRadians(1e-30)));
    }

    [Fact]
    internal void Test_S2CellIdSnapFunction_SnapPoint()
    {
        for (int iter = 0; iter < 1000; ++iter)
        {
            for (int level = 0; level <= S2.kMaxCellLevel; ++level)
            {
                // This checks that points are snapped to the correct level, since
                // S2CellId centers at different levels are always different.
                S2CellIdSnapFunction f = new(level);
                S2Point p = S2Testing.GetRandomCellId(level).ToPoint();
                Assert.Equal(p, f.SnapPoint(p));
            }
        }
    }

    [Fact]
    internal void Test_IntLatLngSnapFunction_ExponentToFromSnapRadius()
    {
        for (int exponent = IntLatLngSnapFunction.kMinExponent;
             exponent <= IntLatLngSnapFunction.kMaxExponent; ++exponent)
        {
            S1Angle radius = IntLatLngSnapFunction.MinSnapRadiusForExponent(exponent);
            Assert.Equal(exponent,
                      IntLatLngSnapFunction.ExponentForMaxSnapRadius(radius));
            Assert.Equal(Math.Min(exponent + 1, IntLatLngSnapFunction.kMaxExponent),
                      IntLatLngSnapFunction.ExponentForMaxSnapRadius(0.999 * radius));
        }
        Assert.Equal(IntLatLngSnapFunction.kMinExponent,
                  IntLatLngSnapFunction.ExponentForMaxSnapRadius(
                      S1Angle.FromRadians(5)));
        Assert.Equal(IntLatLngSnapFunction.kMaxExponent,
                  IntLatLngSnapFunction.ExponentForMaxSnapRadius(
                      S1Angle.FromRadians(1e-30)));
    }

    [Fact]
    internal void Test_IntLatLngSnapFunction_SnapPoint()
    {
        for (int iter = 0; iter < 1000; ++iter)
        {
            // Test that IntLatLngSnapFunction does not modify points that were
            // generated using the S2LatLng.From{E5,E6,E7} methods.  This ensures
            // that both functions are using bitwise-compatible conversion methods.
            S2Point p = S2Testing.RandomPoint();
            S2LatLng ll = new(p);
            S2Point p5 = S2LatLng.FromE5(ll.Lat().E5(), ll.Lng().E5()).ToPoint();
            Assert.Equal(p5, new IntLatLngSnapFunction(5).SnapPoint(p5));
            S2Point p6 = S2LatLng.FromE6(ll.Lat().E6(), ll.Lng().E6()).ToPoint();
            Assert.Equal(p6, new IntLatLngSnapFunction(6).SnapPoint(p6));
            S2Point p7 = S2LatLng.FromE7(ll.Lat().E7(), ll.Lng().E7()).ToPoint();
            Assert.Equal(p7, new IntLatLngSnapFunction(7).SnapPoint(p7));

            // Make sure that we're not snapping using some lower exponent.
            S2Point p7not6 = S2LatLng.FromE7(10 * ll.Lat().E6() + 1,
                                              10 * ll.Lng().E6() + 1).ToPoint();
            Assert.NotEqual(p7not6, new IntLatLngSnapFunction(6).SnapPoint(p7not6));
        }
    }

    [Fact]
    internal void Test_S2CellIdSnapFunction_MinVertexSeparationSnapRadiusRatio()
    {
        // The purpose of this "test" is to compute a lower bound to the fraction
        // (min_vertex_separation() / snap_radius()).  Essentially this involves
        // searching for two adjacent cells A and B such when one of the corner
        // vertices of B is snapped to the center of B, the distance to the center
        // of A decreases as much as possible.  In other words, we want the ratio
        //
        //   distance(center(A), center(B)) / distance(center(A), vertex(B))
        //
        // to be as small as possible.  We do this by considering one cell level at
        // a time, and remembering the cells that had the lowest ratios.  When we
        // proceed from one level to the next, we consider all the children of those
        // cells and keep the best ones.
        //
        // The reason we can restrict the search to children of cells at the
        // previous level is that the ratio above is essentially a function of the
        // local distortions created by projecting the S2 cube space onto the
        // sphere.  These distortions change smoothly over the sphere, so by keeping
        // a fairly large number of candidates ("num_to_keep"), we are essentially
        // keeping all the neighbors of the optimal cell as well.
        double best_score = 1e10;
        List<S2CellId> best_cells = new();
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            double score = GetS2CellIdMinVertexSeparation(level, best_cells);
            best_score = Math.Min(best_score, score);
        }
        _logger.WriteLine($"min_vertex_sep / snap_radius ratio: {best_score:f15}");
    }

    [Fact]
    internal void Test_S2CellIdSnapFunction_MinEdgeVertexSeparationForLevel()
    {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for any
        // snap radius at each level.
        double score = GetS2CellIdMinEdgeSeparation("min_sep_for_level",
            (int level, S1Angle edge_sep, S1Angle min_snap_radius, S1Angle max_snap_radius)
            => edge_sep.Radians / S2.kMinDiag.GetValue(level));
        _logger.WriteLine($"min_edge_vertex_sep / kMinDiag ratio: {score:f15}");
    }
    
    [Fact]
    internal void Test_S2CellIdSnapFunction_MinEdgeVertexSeparationAtMinSnapRadius()
    {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for the
        // special case where the minimum snap radius is being used.
        double score = GetS2CellIdMinEdgeSeparation("min_sep_at_min_radius",
            (int level, S1Angle edge_sep, S1Angle min_snap_radius, S1Angle max_snap_radius) =>
            {
                double min_radius_at_level = S2.kMaxDiag.GetValue(level) / 2;
                return (min_snap_radius.Radians <= (1 + 1e-10) * min_radius_at_level) ?
                  (edge_sep.Radians / S2.kMinDiag.GetValue(level)) : 100.0;
            });
        _logger.WriteLine($"min_edge_vertex_sep / kMinDiag at MinSnapRadiusForLevel: {score:f15}");
    }

    [Fact]
    internal void Test_S2CellIdSnapFunction_MinEdgeVertexSeparationSnapRadiusRatio()
    {
        // Computes the minimum edge separation expressed as a fraction of the
        // maximum snap radius that could yield that edge separation.
        double score = GetS2CellIdMinEdgeSeparation("min_sep_snap_radius_ratio",
            (int level, S1Angle edge_sep, S1Angle min_snap_radius, S1Angle max_snap_radius) =>
            {
                return edge_sep.Radians / max_snap_radius.Radians;
            });
        _logger.WriteLine($"min_edge_vertex_sep / snap_radius ratio: {score:f15}");
    }

    [Fact]
    internal void Test_IntLatLngSnapFunction_MinVertexSeparationSnapRadiusRatio()
    {
        double best_score = 1e10;
        List<IntLatLng> best_configs = new();
        Int64 scale = 18L;
        for (int lat0 = 0; lat0 <= 9; ++lat0)
        {
            best_configs.Add(new(lat0, 0L));
        }
        for (int exp = 0; exp <= 10; ++exp, scale *= 10)
        {
            double score = GetLatLngMinVertexSeparation(scale, 10 * scale, best_configs);
            best_score = Math.Min(best_score, score);
        }
        _logger.WriteLine($"min_vertex_sep / snap_radius ratio: {best_score:f15}");
    }

    [Fact]
    internal void Test_IntLatLngSnapFunction_MinEdgeVertexSeparationForLevel()
    {
        // Computes the minimum edge separation (as a fraction of kMinDiag) for any
        // snap radius at each level.
        double score = GetLatLngMinEdgeSeparation("min_sep_for_level",
            (Int64 scale, S1Angle edge_sep, S1Angle max_snap_radius) =>
            {
                double e_unit = Math.PI / scale;
                return edge_sep.Radians / e_unit;
            });
        _logger.WriteLine($"min_edge_vertex_sep / e_unit ratio: {score:f15}");
    }

    [Fact]
    internal void Test_IntLatLngSnapFunction_MinEdgeVertexSeparationSnapRadiusRatio()
    {
        // Computes the minimum edge separation expressed as a fraction of the
        // maximum snap radius that could yield that edge separation.
        double score = GetLatLngMinEdgeSeparation("min_sep_snap_radius_ratio",
            (Int64 scale, S1Angle edge_sep, S1Angle max_snap_radius) =>
            {
                return edge_sep.Radians / max_snap_radius.Radians;
            });
        _logger.WriteLine($"min_edge_vertex_sep / snap_radius ratio: {score:f15}");
    }

    private static S1Angle GetMaxVertexDistance(S2Point p, S2CellId id)
    {
        S2Cell cell = new(id);
        return new[]
        {
    new S1Angle(p, cell.Vertex(0)),
    new S1Angle(p, cell.Vertex(1)),
    new S1Angle(p, cell.Vertex(2)),
    new S1Angle(p, cell.Vertex(3)),
}.Max();
    }

    // Helper function that computes the vertex separation between "id0" and its
    // neighbors.
    private static void UpdateS2CellIdMinVertexSeparation(S2CellId id0, List<(double, S2CellId)> scores)
    {
        S2Point site0 = id0.ToPoint();
        List<S2CellId> nbrs = new();
        id0.AppendAllNeighbors(id0.Level(), nbrs);
        foreach (S2CellId id1 in nbrs)
        {
            S2Point site1 = id1.ToPoint();
            S1Angle vertex_sep = new(site0, site1);
            S1Angle max_snap_radius = GetMaxVertexDistance(site0, id1);
            Assert.True(max_snap_radius >=
                      S2CellIdSnapFunction.MinSnapRadiusForLevel(id0.Level()));
            double r = vertex_sep / max_snap_radius;
            scores.Add((r, id0));
        }
    }

    private double GetS2CellIdMinVertexSeparation(int level, List<S2CellId> best_cells)
    {
        // The worst-case separation ratios always occur when the snap_radius is not
        // much larger than the minimum, since this allows the site spacing to be
        // reduced by as large a fraction as possible.
        //
        // For the minimum vertex separation ratio, we choose a site and one of its
        // 8-way neighbors, then look at the ratio of the distance to the center of
        // that neighbor to the distance to the furthest corner of that neighbor
        // (which is the largest possible snap radius for this configuration).
        List<(double ratio, S2CellId)> scores = new();
        if (level == 0)
        {
            UpdateS2CellIdMinVertexSeparation(kSearchRootId, scores);
        }
        else
        {
            foreach (var parent in best_cells)
            {
                for (var id0 = parent.ChildBegin();
                     id0 != parent.ChildEnd(); id0 = id0.Next())
                {
                    UpdateS2CellIdMinVertexSeparation(id0, scores);
                }
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and keep
        // the best "num_to_keep" of them to seed the next round.
        scores = new SortedSet<(double, S2CellId)>(scores).ToList();
        best_cells.Clear();
        int num_to_keep = 300;
        int num_to_print = 1;
        foreach (var (ratio, id) in scores)
        {
            if (--num_to_print >= 0)
            {
                R2Point uv = id.CenterUV();
                _logger.WriteLine($"Level {level:d2}: min_vertex_sep_ratio = {ratio:f15} u={uv[0]:f6} v={uv[1]:f6} {id.ToToken()}");
            }
            if (kSearchFocusId.Contains(id) || id.Contains(kSearchFocusId))
            {
                var inserted = best_cells.AddSortedUnique(id);
                if (inserted && --num_to_keep <= 0) break;
            }
        }
        return scores[0].ratio;
    }

    private static S1Angle GetCircumRadius(S2Point a, S2Point b, S2Point c)
    {
        // We return this value is the circumradius is very large.
        S1Angle kTooBig = S1Angle.FromRadians(Math.PI);
        double turn_angle = S2.TurnAngle(a, b, c);
        if (Math.Abs(Math.IEEERemainder(turn_angle, Math.PI)) < 1e-2) return kTooBig;

        double a2 = (b - c).Norm2();
        double b2 = (c - a).Norm2();
        double c2 = (a - b).Norm2();
        if (a2 > 2 || b2 > 2 || c2 > 2) return kTooBig;
        double ma = a2 * (b2 + c2 - a2);
        double mb = b2 * (c2 + a2 - b2);
        double mc = c2 * (a2 + b2 - c2);
        S2Point p = (ma * a + mb * b + mc * c) / (ma + mb + mc);
        return new S1Angle(p, a);
    }

    private static List<S2CellId> GetNeighbors(S2CellId id)
    {
        int kNumLayers = 2;
        List<S2CellId> nbrs = new()
        {
            id
        };
        for (int layer = 0; layer < kNumLayers; ++layer)
        {
            List<S2CellId> new_nbrs = new();
            foreach (var nbr in nbrs)
            {
                nbr.AppendAllNeighbors(id.Level(), new_nbrs);
            }
            nbrs.AddRange(new_nbrs);
            nbrs.Remove(id);
            nbrs = new SortedSet<S2CellId>(nbrs).ToList();
        }
        return nbrs;
    }

    // S2CellIdMinEdgeSeparationFunction defines an objective function that will
    // be optimized by GetS2CellIdMinEdgeSeparation() by finding worst-case
    // configurations of S2CellIds.  We use this to find the worst cases under
    // various conditions (e.g., when the minimum snap radius at a given level is
    // being used).  The objective function is called for a specific configuration
    // of vertices that are snapped at the given S2CellId level.  "edge_sep" is
    // the edge-vertex distance that is achieved by this configuration, and
    // "min_snap_radius" and "max_snap_radius" are the minimum and maximum snap
    // radii for which this configuration is valid (i.e., where the desired
    // snapping will take place).
    private delegate double S2CellIdMinEdgeSeparationFunction(int level, S1Angle edge_sep,
        S1Angle min_snap_radius, S1Angle max_snap_radius);

    // Returns the minimum value of the given objective function over sets of
    // nearby vertices that are designed to minimize the edge-vertex separation
    // when an edge is snapped.
    private double GetS2CellIdMinEdgeSeparation(
        string label, S2CellIdMinEdgeSeparationFunction objective,
        int level, List<S2CellId> best_cells)
    {
        // To find minimum edge separations, we choose a cell ("id0") and two nearby
        // cells ("id1" and "id2"), where "nearby" is defined by GetNeighbors().
        // Let "site0", "site1", and "site2" be the centers of these cells.  The
        // idea is to consider an input edge E that intersects the Voronoi regions
        // of "site1" and "site2" (and therefore snaps to an edge E' between these
        // sites) but does not not intersect the Voronoi region of "site0" (and
        // therefore can't be snapped to site0).  The goal is to search for snapped
        // edges E' that approach site0 as closely as possible.
        //
        // To do this, we first compute the circumradius of the three cell centers
        // ("site0", "site1", and "site2"); this is the minimum snap radius in order
        // for it to be possible to construct an edge E that snaps to "site1" and
        // "site2" but not to "site0".  We also compute the distance from "site0" to
        // the snapped edge.  Next we find the corner vertex of "id1" and "id2" that
        // is furthest from "site0"; the smaller of these two distances is the
        // maximum snap radius such that "site1" and "site2" can be chosen as
        // sites after choosing "site0".  If the maximum is less than the minimum,
        // then this configuration is rejected; otherwise we evaluate the given
        // objective function and keep the configurations that result in the
        // smallest values.
        //
        // The optimization process works by keeping track of the set of S2CellIds
        // that yielded the best results at the previous level, and exploring all
        // the nearby neighbor combinations of the children of those cells at the
        // next level.  In order to get better coverage, we keep track of the best
        // score and configuration (i.e. the two neighboring cells "id1" and "id2")
        // for each initial cell "id0".
        Dictionary<S2CellId, double> best_scores = new();
        Dictionary<S2CellId, (S2CellId, S2CellId)> best_configs = new();
        foreach (var parent in best_cells)
        {
            for (var id0 = parent.ChildBegin(level);
                 id0 != parent.ChildEnd(level); id0 = id0.Next())
            {
                S2Point site0 = id0.ToPoint();
                var nbrs = GetNeighbors(id0);
                foreach (S2CellId id1 in nbrs)
                {
                    S2Point site1 = id1.ToPoint();
                    S1Angle max_v1 = GetMaxVertexDistance(site0, id1);
                    foreach (S2CellId id2 in nbrs)
                    {
                        if (id2 <= id1) continue;
                        S2Point site2 = id2.ToPoint();
                        S1Angle min_snap_radius = GetCircumRadius(site0, site1, site2);
                        if (min_snap_radius > SnapFunction.kMaxSnapRadius)
                        {
                            continue;
                        }
                        // Note that it is only the original points *before* snapping that
                        // need to be at least "snap_radius" away from "site0".  The points
                        // after snapping ("site1" and "site2") may be closer.
                        S1Angle max_v2 = GetMaxVertexDistance(site0, id2);
                        S1Angle max_snap_radius = new[] { max_v1, max_v2 }.Min();
                        if (min_snap_radius > max_snap_radius) continue;
                        Assert.True(max_snap_radius >= S2CellIdSnapFunction.MinSnapRadiusForLevel(level));

                        // This is a valid configuration, so evaluate it.
                        S1Angle edge_sep = S2.GetDistance(site0, site1, site2);
                        double score = objective(level, edge_sep, min_snap_radius, max_snap_radius);
                        var best_score = best_scores[id0];
                        if (best_score == 0 || best_score > score)
                        {
                            best_score = score;
                            best_configs[id0] = (id1, id2);
                        }
                    }
                }
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and
        // generate a set of candidates for the next round by generating all the
        // 8-way neighbors of the best candidates, and keeping up to"num_to_keep" of
        // them.  The results vary slightly according to how many candidates we
        // keep, but the variations are much smaller than the conservative
        // assumptions made by the S2CellIdSnapFunction implementation.
        int num_to_keep = 100;
#if DEBUG
        num_to_keep = 20;
#endif
        int num_to_print = 3;
        List<(double ratio, S2CellId)> sorted = new();
        foreach (var entry in best_scores)
        {
            sorted.Add((entry.Value, entry.Key));
        }
        sorted.Sort();
        best_cells.Clear();
        _logger.WriteLine($"Level {level}:");
        foreach (var (ratio, id) in sorted)
        {
            if (--num_to_print >= 0)
            {
                R2Point uv = id.CenterUV();
                var (c1, c2) = best_configs[id];
                _logger.WriteLine($"  {label} = {ratio:f15} u={uv[0]:f7.4} v={uv[1]:f7.4} {id.ToToken()} {c1.ToToken()} {c2.ToToken()}");
            }
            List<S2CellId> nbrs = new(1)
            {
                id
            };
            id.AppendAllNeighbors(id.Level(), nbrs);
            foreach (var nbr in nbrs)
            {
                // The S2Cell hierarchy has many regions that are symmetrical.  We can
                // eliminate most of the "duplicates" by restricting the search to cells
                // in kS2CellIdFocus.
                if (kSearchFocusId.Contains(nbr) || nbr.Contains(kSearchFocusId))
                {
                    var inserted = best_cells.AddSortedUnique(nbr);
                    if (inserted && --num_to_keep <= 0)
                    {
                        return sorted[0].ratio;
                    }
                }
            }
        }
        return sorted[0].ratio;
    }

    private double GetS2CellIdMinEdgeSeparation(
        string label, S2CellIdMinEdgeSeparationFunction objective)
    {
        double best_score = 1e10;
        List<S2CellId> best_cells = new()
        {
            kSearchRootId
        };
        for (int level = 0; level <= S2.kMaxCellLevel; ++level)
        {
            double score = GetS2CellIdMinEdgeSeparation(label, objective, level, best_cells);
            best_score = Math.Min(best_score, score);
        }
        return best_score;
    }

    private static bool IsValid(IntLatLng ll, Int64 scale)
    {
        // A coordinate value of "scale" corresponds to 180 degrees.
        return (Math.Abs(ll.Lat) <= scale / 2 && Math.Abs(ll.Lng) <= scale);
    }

    private static bool HasValidVertices(IntLatLng ll, Int64 scale)
    {
        // Like IsValid, but excludes latitudes of 90 and longitudes of 180.
        // A coordinate value of "scale" corresponds to 180 degrees.
        return (Math.Abs(ll.Lat) < scale / 2 && Math.Abs(ll.Lng) < scale);
    }

    private static IntLatLng Rescale(IntLatLng ll, double scale_factor)
    {
        return new(
            Convert.ToInt64(scale_factor * ll.Lat),
            Convert.ToInt64(scale_factor * ll.Lng));
    }

    private static S2Point ToPoint(IntLatLng ll, Int64 scale)
    {
        return S2LatLng.FromRadians(
            ll.Lat * (Math.PI / scale),
            ll.Lng * (Math.PI / scale)).ToPoint();
    }

    private static S2Point GetVertex(IntLatLng ll, Int64 scale, int i)
    {
        // Return the points in CCW order starting from the lower left.
        int dlat = (i == 0 || i == 3) ? -1 : 1;
        int dlng = (i == 0 || i == 1) ? -1 : 1;
        return ToPoint(new(2 * ll.Lat + dlat, 2 * ll.Lng + dlng), 2 * scale);
    }

    private static S1Angle GetMaxVertexDistance(S2Point p, IntLatLng ll, Int64 scale)
    {
        return new[]
        {
    new S1Angle(p, GetVertex(ll, scale, 0)),
    new S1Angle(p, GetVertex(ll, scale, 1)),
    new S1Angle(p, GetVertex(ll, scale, 2)),
    new S1Angle(p, GetVertex(ll, scale, 3)),
}.Max();
    }

    private double GetLatLngMinVertexSeparation(Int64 old_scale, Int64 scale, List<IntLatLng> best_configs)
    {
        // The worst-case separation ratios always occur when the snap_radius is not
        // much larger than the minimum, since this allows the site spacing to be
        // reduced by as large a fraction as possible.
        //
        // For the minimum vertex separation ratio, we choose a site and one of its
        // 8-way neighbors, then look at the ratio of the distance to the center of
        // that neighbor to the distance to the furthest corner of that neighbor
        // (which is the largest possible snap radius for this configuration).
        S1Angle min_snap_radius_at_scale = S1Angle.FromRadians(S2.M_SQRT1_2 * Math.PI / scale);
        List<(double ratio, IntLatLng ll)> scores = new();
        double scale_factor = (double)scale / old_scale;
        foreach (var parent in best_configs)
        {
            var (lat, lng) = Rescale(parent, scale_factor);
            for (long dlat0 = -7; dlat0 <= 7; ++dlat0)
            {
                IntLatLng ll0 = new(lat + dlat0, lng);
                if (!IsValid(ll0, scale) || ll0.Lat < 0) continue;
                S2Point site0 = ToPoint(ll0, scale);
                for (int dlat1 = 0; dlat1 <= 2; ++dlat1)
                {
                    for (int dlng1 = 0; dlng1 <= 5; ++dlng1)
                    {
                        IntLatLng ll1 = new(ll0.Lat + dlat1, ll0.Lng + dlng1);
                        if (ll1 == ll0 || !HasValidVertices(ll1, scale)) continue;
                        S1Angle max_snap_radius = GetMaxVertexDistance(site0, ll1, scale);
                        if (max_snap_radius < min_snap_radius_at_scale) continue;
                        S2Point site1 = ToPoint(ll1, scale);
                        S1Angle vertex_sep = new(site0, site1);
                        double r = vertex_sep / max_snap_radius;
                        scores.AddSortedUnique((r, ll0));
                    }
                }
            }
        }

        best_configs.Clear();
        int num_to_keep = 100;
        int num_to_print = 1;
        foreach (var (ratio, ll) in scores)
        {
            if (--num_to_print >= 0)
            {
                var ds = ToPoint(ll, scale).ToDebugString();
                _logger.WriteLine($"Scale {scale:14}: min_vertex_sep_ratio = {ratio:f15},{ds}");
            }
            var inserted = best_configs.AddSortedUnique(ll);
            if (inserted && --num_to_keep <= 0) break;
        }
        return scores[0].ratio;
    }

    private delegate double LatLngMinEdgeSeparationFunction(Int64 scale, S1Angle edge_sep, S1Angle max_snap_radius);

    private double GetLatLngMinEdgeSeparation(
        string label, LatLngMinEdgeSeparationFunction objective,
        Int64 scale, List<LatLngConfig> best_configs)
    {
        var min_snap_radius_at_scale = S1Angle.FromRadians(S2.M_SQRT1_2 * Math.PI / scale);
        List<LatLngConfigScore> scores = new();
        for (var i = 0; i < best_configs.Count; i++)
        {
            var (scale_, ll0_, ll1_, ll2_) = best_configs[i];
            // To reduce duplicates, we require that site0 always has longitude 0.
            Assert.Equal(0, ll0_.Lng);
            double scale_factor = (double)scale / scale_;
            ll0_ = Rescale(ll0_, scale_factor);
            ll1_ = Rescale(ll1_, scale_factor);
            ll2_ = Rescale(ll2_, scale_factor);
            for (int dlat0 = -1; dlat0 <= 1; ++dlat0)
            {
                IntLatLng ll0 = new(ll0_.Lat + dlat0, ll0_.Lng);
                // To reduce duplicates, we require that site0.latitude >= 0.
                if (!IsValid(ll0, scale) || ll0.Lat < 0) continue;

                var site0 = ToPoint(ll0, scale);
                for (int dlat1 = -1; dlat1 <= 1; ++dlat1)
                {
                    for (int dlng1 = -2; dlng1 <= 2; ++dlng1)
                    {
                        IntLatLng ll1 = new(ll1_.Lat + dlat0 + dlat1, dlng1);
                        if (ll1 == ll0 || !HasValidVertices(ll1, scale)) continue;
                        // Only consider neighbors within 2 latitude units of site0.
                        if (Math.Abs(ll1.Lat - ll0.Lat) > 2) continue;

                        var site1 = ToPoint(ll1, scale);
                        var max_v1 = GetMaxVertexDistance(site0, ll1, scale);
                        for (int dlat2 = -1; dlat2 <= 1; ++dlat2)
                        {
                            for (int dlng2 = -2; dlng2 <= 2; ++dlng2)
                            {
                                IntLatLng ll2 = new(ll2_.Lat + dlat0 + dlat2, ll2_.Lng + dlng2);
                                if (!HasValidVertices(ll2, scale)) continue;
                                // Only consider neighbors within 2 latitude units of site0.
                                if (Math.Abs(ll2.Lat - ll0.Lat) > 2) continue;
                                // To reduce duplicates, we require ll1 < ll2 lexicographically
                                // and site2.longitude >= 0.  (It's *not* okay to
                                // require site1.longitude >= 0, because then some configurations
                                // with site1.latitude == site2.latitude would be missed.)
                                if (ll2.CompareTo(ll1) <= 0 || ll2.Lng < 0) continue;

                                S2Point site2 = ToPoint(ll2, scale);
                                S1Angle min_snap_radius = GetCircumRadius(site0, site1, site2);
                                if (min_snap_radius > SnapFunction.kMaxSnapRadius)
                                {
                                    continue;
                                }
                                // Only the original points *before* snapping that need to be at
                                // least "snap_radius" away from "site0".  The points after
                                // snapping ("site1" and "site2") may be closer.
                                S1Angle max_v2 = GetMaxVertexDistance(site0, ll2, scale);
                                S1Angle max_snap_radius = new[] { max_v1, max_v2 }.Min();
                                if (min_snap_radius > max_snap_radius) continue;
                                if (max_snap_radius < min_snap_radius_at_scale) continue;

                                // This is a valid configuration, so evaluate it.
                                S1Angle edge_sep = S2.GetDistance(site0, site1, site2);
                                double score = objective(scale, edge_sep, max_snap_radius);
                                LatLngConfig config = new(scale, ll0, ll1, ll2);
                                scores.Add(new(score, config));
                            }
                        }
                    }
                }
            }
        }
        // Now sort the entries, print out the "num_to_print" best ones, and keep
        // the best "num_to_keep" of them to seed the next round.
        scores = new SortedSet<LatLngConfigScore>(scores).ToList();
        best_configs.Clear();
        int num_to_keep = 200;
#if DEBUG
        num_to_keep = 50;
#endif
        int num_to_print = 3;
        _logger.WriteLine($"Scale {scale}:");
        foreach (var (ratio, config) in scores)
        {
            var scale2 = config.Scale;
            if (--num_to_print >= 0)
            {
                var s1 = ToPoint(config.LL0, scale2).ToDebugString();
                var s2 = ToPoint(config.LL1, scale2).ToDebugString();
                var s3 = ToPoint(config.LL2, scale2).ToDebugString();
                _logger.WriteLine($"  {label} = {ratio:f15} {s1} {s2} {s3}");
            }
            // Optional: filter the candidates to concentrate on a specific region
            // (e.g., the north pole).
            best_configs.Add(config);
            if (--num_to_keep <= 0) break;
        }
        return scores[0].Score;
    }

    private double GetLatLngMinEdgeSeparation(
        string label, LatLngMinEdgeSeparationFunction objective)
    {
        double best_score = 1e10;
        List<LatLngConfig> best_configs = new();
        Int64 scale = 6L;  // Initially points are 30 degrees apart.
        int max_lng = (int)scale;
        int max_lat = (int)(scale / 2);
        for (int lat0 = 0; lat0 <= max_lat; ++lat0)
        {
            for (int lat1 = lat0 - 2; lat1 <= Math.Min(max_lat, lat0 + 2); ++lat1)
            {
                for (int lng1 = 0; lng1 <= max_lng; ++lng1)
                {
                    for (int lat2 = lat1; lat2 <= Math.Min(max_lat, lat0 + 2); ++lat2)
                    {
                        for (int lng2 = 0; lng2 <= max_lng; ++lng2)
                        {
                            IntLatLng ll0 = new(lat0, 0);
                            IntLatLng ll1 = new(lat1, lng1);
                            IntLatLng ll2 = new(lat2, lng2);
                            if (ll2.CompareTo(ll1) <= 0) continue;

                            best_configs.Add(new(scale, ll0, ll1, ll2));
                        }
                    }
                }
            }
        }
        Int64 target_scale = 180L;
        for (int exp = 0; exp <= 10; ++exp, target_scale *= 10)
        {
            while (scale < target_scale)
            {
                scale = Math.Min((Int64)(1.8 * scale), target_scale);
                double score = GetLatLngMinEdgeSeparation(label, objective, scale,
                                                          best_configs);
                if (scale == target_scale)
                {
                    best_score = Math.Min(best_score, score);
                }
            }
        }
        return best_score;
    }

    // A scaled S2LatLng with integer coordinates, similar to E7 coordinates,
    // except that the scale is variable (see LatLngConfig below).
    private record IntLatLng(long Lat, long Lng) : IComparable<IntLatLng>
    {
        public int CompareTo(IntLatLng? other) 
        {
            if (other is null) return 1;

            var c = Lat.CompareTo(other.Lat);
            if (c != 0) return c;

            return Lng.CompareTo(other.Lng);
        }
    }
    // A triple of scaled S2LatLng coordinates.  The coordinates are multiplied by
    // (Math.PI / scale) to convert them to radians.
    private record LatLngConfig(long Scale, IntLatLng LL0, IntLatLng LL1, IntLatLng LL2) : IComparable<LatLngConfig>
    {
        public int CompareTo(LatLngConfig? other)
        {
            if (other is null) return 1;

            var c = LL0.CompareTo(other.LL0);
            if (c != 0) return c;

            c = LL1.CompareTo(other.LL1);
            if (c != 0) return c;

            return LL2.CompareTo(other.LL2);
        }
    }
    private record LatLngConfigScore(double Score, LatLngConfig LLC);
}
