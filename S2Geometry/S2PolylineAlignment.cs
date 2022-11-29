// This library provides code to compute vertex alignments between S2Polylines.
//
// A vertex "alignment" or "warp" between two polylines is a matching between
// pairs of their vertices. Users can imagine pairing each vertex from
// S2Polyline `a` with at least one other vertex in S2Polyline `b`. The "cost"
// of an arbitrary alignment is defined as the summed value of the squared
// chordal distance between each pair of points in the warp path. An "optimal
// alignment" for a pair of polylines is defined as the alignment with least
// cost. Note: optimal alignments are not necessarily unique. The standard way
// of computing an optimal alignment between two sequences is the use of the
// `Dynamic Timewarp` algorithm.
//
// We provide three methods for computing (via Dynamic Timewarp) the optimal
// alignment between two S2Polylines. These methods are performance-sensitive,
// and have been reasonably optimized for space- and time- usage. On modern
// hardware, it is possible to compute exact alignments between 4096x4096
// element polylines in ~70ms, and approximate alignments much more quickly.
//
// The results of an alignment operation are captured in a VertexAlignment
// object. In particular, a VertexAlignment keeps track of the total cost of
// alignment, as well as the warp path (a sequence of pairs of indices into each
// polyline whose vertices are linked together in the optimal alignment)
//
// For a worked example, consider the polylines
//
// a = [(1, 0), (5, 0), (6, 0), (9, 0)] and
// b = [(2, 0), (7, 0), (8, 0)].
//
// The "cost matrix" between these two polylines (using squared chordal
// distance, .Norm2, as our distance function) looks like this:
//
//        (2, 0)  (7, 0)  (8, 0)
// (1, 0)     1      36      49
// (5, 0)     9       4       9
// (6, 0)    16       1       4
// (9, 0)    49       4       1
//
// The Dynamic Timewarp DP table for this cost matrix has cells defined by
//
// table[i][j] = cost(i,j) + Math.Min(table[i-1][j-1], table[i][j-1], table[i-1, j])
//
//        (2, 0)  (7, 0)  (8, 0)
// (1, 0)     1      37      86
// (5, 0)    10       5      14
// (6, 0)    26       6       9
// (9, 0)    75      10       7
//
// Starting at the bottom right corner of the DP table, we can work our way
// backwards to the upper left corner  to recover the reverse of the warp path:
// (3, 2) . (2, 1) . (1, 1) . (0, 0). The VertexAlignment produced containing
// this has alignment_cost = 7 and warp_path = {(0, 0), (1, 1), (2, 1), (3, 2)}.
//
// We also provide methods for performing alignment of multiple sequences. These
// methods return a single, representative polyline from a non-empty collection
// of polylines, for various definitions of "representative."
//
// GetMedoidPolyline() returns a new polyline (point-for-point-equal to some
// existing polyline from the collection) that minimizes the summed vertex
// alignment cost to all other polylines in the collection.
//
// GetConsensusPolyline() returns a new polyline (unlikely to be present in the
// input collection) that represents a "weighted consensus" polyline. This
// polyline is constructed iteratively using the Dynamic Timewarp Barycenter
// Averaging algorithm of F. Petitjean, A. Ketterlin, and P. Gancarski, which
// can be found here:
// https://pdfs.semanticscholar.org/a596/8ca9488199291ffe5473643142862293d69d.pdf

namespace S2Geometry;

using System.Text;

using WarpPath = List<(int start, int end)>;

public static class S2PolylineAlignment
{
    #region VertexAlignment

    // GetExactVertexAlignment takes two non-empty polylines as input, and returns
    // the VertexAlignment corresponding to the optimal alignment between them. This
    // method is quadratic O(A*B) in both space and time complexity.
    public static VertexAlignment GetExactVertexAlignment(S2Polyline a, S2Polyline b)
    {
        int a_n = a.NumVertices();
        int b_n = b.NumVertices();
        System.Diagnostics.Debug.Assert(a_n > 0); // A is empty polyline.
        System.Diagnostics.Debug.Assert(b_n > 0); // B is empty polyline.
        var w = new Window(new (int start, int end)[a_n].Fill((0, b_n)));
        return DynamicTimewarp(a, b, w);
    }

    // GetExactVertexAlignmentCost takes two non-empty polylines as input, and
    // returns the *cost* of their optimal alignment. A standard, traditional
    // dynamic timewarp algorithm can output both a warp path and a cost, but
    // requires quadratic space to reconstruct the path by walking back through the
    // Dynamic Programming cost table. If all you really need is the warp cost (i.e.
    // you're inducing a similarity metric between S2Polylines, or something
    // equivalent), you can overwrite the DP table and useant space -
    // O(max(A,B)). This method provides that space-efficiency optimization.
    public static double GetExactVertexAlignmentCost(S2Polyline a, S2Polyline b)
    {
        // This is the constant-space implementation of Dynamic Timewarp that can
        // compute the alignment cost, but not the warp path.

        int a_n = a.NumVertices();
        int b_n = b.NumVertices();
        System.Diagnostics.Debug.Assert(a_n > 0); // A is empty polyline.
        System.Diagnostics.Debug.Assert(b_n > 0); // B is empty polyline.
        var cost = new double[b_n].Fill(DOUBLE_MAX);
        double left_diag_min_cost = 0;
        for (int row = 0; row < a_n; ++row)
        {
            for (int col = 0; col < b_n; ++col)
            {
                double up_cost = cost[col];
                cost[col] = Math.Min(left_diag_min_cost, up_cost) +
                            (a.Vertex(row) - b.Vertex(col)).Norm2();
                left_diag_min_cost = Math.Min(cost[col], up_cost);
            }
            left_diag_min_cost = DOUBLE_MAX;
        }
        return cost.Last();
    }

    // GetApproxVertexAlignment takes two non-empty polylines `a` and `b` as input,
    // and a `radius` paramater GetApproxVertexAlignment (quickly) computes an
    // approximately optimal vertex alignment of points between polylines `a` and
    // `b` by implementing the algorithm described in `FastDTW: Toward Accurate
    // Dynamic Time Warping in Linear Time and Space` by Stan Salvador and Philip
    // Chan. Details can be found below:
    //
    // https://pdfs.semanticscholar.org/05a2/0cde15e172fc82f32774dd0cf4fe5827cad2.pdf
    //
    // The `radius` parameter controls the distance we search outside of the
    // projected warp path during the refining step. Smaller values of `radius`
    // correspond to a smaller search window, and therefore distance computation on
    // fewer cells, which leads to a faster (but worse) approximation.
    // This method is O(max(A, B)) in both space and time complexity.
    public static VertexAlignment GetApproxVertexAlignment(S2Polyline a, S2Polyline b, int radius)
    {
        // Determined experimentally, through benchmarking, as about the points at
        // which ExactAlignment is faster than ApproxAlignment, so we use these as
        // our switchover points to exact computation mode.
        const int kSizeSwitchover = 32;
        const double kDensitySwitchover = 0.85;
        int a_n = a.NumVertices();
        int b_n = b.NumVertices();
        System.Diagnostics.Debug.Assert(a_n > 0); // A is empty polyline.
        System.Diagnostics.Debug.Assert(b_n > 0); // B is empty polyline.
        System.Diagnostics.Debug.Assert(radius >= 0); // Radius is negative.

        // If we've hit the point where doing a full, direct solve is guaranteed to
        // be faster, then terminate the recursion and do that.
        if (a_n - radius < kSizeSwitchover || b_n - radius < kSizeSwitchover)
        {
            return GetExactVertexAlignment(a, b);
        }

        // If we've hit the point where the window will be probably be so full that we
        // might as well compute an exact solution, then terminate recursion to do so.
        if (Math.Max(a_n, b_n) * (2 * radius + 1) > a_n * b_n * kDensitySwitchover)
        {
            return GetExactVertexAlignment(a, b);
        }

        // Otherwise, shrink the input polylines, recursively compute the vertex
        // alignment using this method, and then compute the final alignment using
        // the projected alignment `proj` on an upsampled, dilated window.
        var a_half = HalfResolution(a);
        var b_half = HalfResolution(b);
        var proj = GetApproxVertexAlignment(a_half, b_half, radius);
        var w = new Window(proj.WarpPath_).Upsample(a_n, b_n).Dilate(radius);
        return DynamicTimewarp(a, b, w);
    }

    // A convience overload for GetApproxVertexAlignment which computes and uses
    // suggested default parameter of radius = Math.Max(a.size(), b.size())^0.25
    public static VertexAlignment GetApproxVertexAlignment(S2Polyline a, S2Polyline b)
    {
        int max_length = Math.Max(a.NumVertices(), b.NumVertices());
        int radius = (int)Math.Pow(max_length, 0.25);
        return GetApproxVertexAlignment(a, b, radius);
    }

    public readonly struct VertexAlignment
    {
        // `alignment_cost` represents the sum of the squared chordal distances
        // between each pair of vertices in the warp path. Specifically,
        // cost = sum_{(i, j) \in path} (a.vertex(i) - b.vertex(j)).Norm2;
        // This means that the units of alignment_cost are "squared distance". This is
        // an optimization to avoid the (expensive) atan computation of the true
        // spherical angular distance between the points, as well as an unnecessary
        // square root. All we need to compute vertex alignment is a metric that
        // satisifies the triangle inequality, and squared chordal distance works as
        // well as spherical S1Angle distance for this purpose.
        public readonly double AlignmentCost;

        // Each entry (i, j) of `warp_path` represents a pairing between vertex
        // a.vertex(i) and vertex b.vertex(j) in the optimal alignment.
        // The warp_path is defined in forward order, such that the result of
        // aligning polylines `a` and `b` is always a warp_path with warp_path.First()
        // = {0,0} and warp_path.Last() = {a.NumVertices - 1, b.NumVertices - 1}
        // Note that this DOES NOT define an alignment from a point sequence to an
        // edge sequence. That functionality may come at a later date.
        public readonly WarpPath WarpPath_;

        public VertexAlignment(double cost, WarpPath path)
        {
            AlignmentCost = cost;
            WarpPath_ = path;
        }
    }

    // Perform dynamic timewarping by filling in the DP table on cells that are
    // inside our search window. For an exact (all-squares) evaluation, this
    // incurs bounds checking overhead - we don't need to ensure that we're inside
    // the appropriate cells in the window, because it's guaranteed. Structuring
    // the program to reuse code for both the EXACT and WINDOWED cases by
    // abstracting EXACT as a window with full-covering strides is done for
    // maintainability reasons. One potential optimization here might be to overload
    // this function to skip bounds checking when the window is full.
    //
    // As a note of general interest, the Dynamic Timewarp algorithm as stated here
    // prefers shorter warp paths, when two warp paths might be equally costly. This
    // is because it favors progressing in the sequences simultaneously due to the
    // equal weighting of a diagonal step in the cost table with a horizontal or
    // vertical step. This may be counterintuitive, but represents the standard
    // implementation of this algorithm. TODO(user) - future implementations could
    // allow weights on the lookup costs to mitigate this.
    //
    // This is the hottest routine in the whole package, please be careful to
    // profile any future changes made here.
    //
    // This method takes time proportional to the number of cells in the window,
    // which can range from O(max(a, b)) cells (best) to O(a*b) cells (worst)
    private static VertexAlignment DynamicTimewarp(S2Polyline a, S2Polyline b, Window w)
    {
        var rows = a.NumVertices();
        var cols = b.NumVertices();
        var costs = new double[rows][].Fill(() => new double[cols]);

        (int start, int end) curr;
        var prev = ColumnStrideAll();

        for (int row1 = 0; row1 < rows; ++row1)
        {
            curr = w.GetColumnStride(row1);
            for (var col1 = curr.start; col1 < curr.end; ++col1)
            {
                var d_cost = BoundsCheckedTableCost(row1 - 1, col1 - 1, prev, costs);
                var u_cost = BoundsCheckedTableCost(row1 - 1, col1 - 0, prev, costs);
                var l_cost = BoundsCheckedTableCost(row1 - 0, col1 - 1, curr, costs);
                costs[row1][col1] = new[] { d_cost, u_cost, l_cost }.Min() + (a.Vertex(row1) - b.Vertex(col1)).Norm2();
            }
            prev = curr;
        }

        // Now we walk back through the cost table and build up the warp path.
        // Somewhat surprisingly, it is faster to recover the path this way than it
        // is to save the comparisons from the computation we *already did* to get the
        // direction we came from. The author speculates that this behavior is
        // assignment-cost-related: to persist direction, we have to do extra
        // stores/loads of "directional" information, and the extra assignment cost
        // this incurs is larger than the cost to simply redo the comparisons.
        // It's probably worth revisiting this assumption in the future.
        // As it turns out, the following code ends up effectively free.
        var warp_path = new WarpPath(Math.Max(a.NumVertices(), b.NumVertices()));

        int row = a.NumVertices() - 1;
        int col = b.NumVertices() - 1;
        curr = w.GetCheckedColumnStride(row);
        prev = w.GetCheckedColumnStride(row - 1);
        while (row >= 0 && col >= 0)
        {
            warp_path.Add((row, col));
            double d_cost = BoundsCheckedTableCost(row - 1, col - 1, prev, costs);
            double u_cost = BoundsCheckedTableCost(row - 1, col - 0, prev, costs);
            double l_cost = BoundsCheckedTableCost(row - 0, col - 1, curr, costs);
            if (d_cost <= u_cost && d_cost <= l_cost)
            {
                row -= 1;
                col -= 1;
                curr = w.GetCheckedColumnStride(row);
                prev = w.GetCheckedColumnStride(row - 1);
            }
            else if (u_cost <= l_cost)
            {
                row -= 1;
                curr = w.GetCheckedColumnStride(row);
                prev = w.GetCheckedColumnStride(row - 1);
            }
            else
            {
                col -= 1;
            }
        }
        warp_path.Reverse();
        return new VertexAlignment(costs.Last().Last(), warp_path);
    }

    private static double BoundsCheckedTableCost(int row, int col, (int start, int end) stride, double[][] table)
    {
        if (row < 0 && col < 0)
        {
            return 0.0;
        }
        else if (row < 0 || col < 0 || !InRange(col, stride))
        {
            return DOUBLE_MAX;
        }
        else
        {
            return table[row][col];
        }
    }

    private const double DOUBLE_MAX = double.MaxValue;

    // Reduce the number of vertices of polyline `in` by selecting every other
    // vertex for inclusion in a new polyline. Specifically, we take even-index
    // vertices [0, 2, 4,...]. For an even-length polyline, the last vertex is not
    // selected. For an odd-length polyline, the last vertex is selected.
    // Constructs and returns a new S2Polyline in linear time.
    public static S2Polyline HalfResolution(S2Polyline input)
    {
        int n = input.NumVertices();
        var vertices = new List<S2Point>();
        for (int i = 0; i < n; i += 2)
        {
            vertices.Add(input.Vertex(i));
        }
        return new S2Polyline(vertices.ToArray());
    }

    // A ColumnStride is a [start, end) range of columns in a search window.
    // It enables us to lazily fill up our CostTable structures by providing bounds
    // checked access for reads. We also use them to keep track of structured,
    // sparse Window matrices by tracking start and end columns for each row.
    public static bool InRange(int index, (int start, int end) columnStride)
    {
        return columnStride.start <= index && index < columnStride.end;
    }

    // Returns a ColumnStride where InRange evaluates to `true` for all
    // non-negative inputs less than int.MaxValue;
    public static (int start, int end) ColumnStrideAll()
    {
        return (-1, int.MaxValue);
    }

    #endregion

    #region Medoid

    // We use some of the symmetry of our metric to avoid computing all N^2
    // alignments. Specifically, because cost_fn(a, b) = cost_fn(b, a), and
    // cost_fn(a, a) = 0, we can compute only the lower triangle of cost matrix
    // and then mirror it across the diagonal to save on cost_fn invocations.
    public static int GetMedoidPolyline(S2Polyline[] polylines, MedoidOptions options)
    {
        int num_polylines = polylines.Length;
        bool approx = options.Approx;
        System.Diagnostics.Debug.Assert(num_polylines > 0);

        // costs[i] stores total cost of aligning [i] with all other polylines.
        var costs = new double[num_polylines].Fill(0.0);
        for (int i = 0; i < num_polylines; ++i)
        {
            for (int j = i + 1; j < num_polylines; ++j)
            {
                double cost = CostFn(polylines[i], polylines[j], approx);
                costs[i] += cost;
                costs[j] += cost;
            }
        }
        return costs.IndexOfMin();
    }

    // Helper methods for GetMedoidPolyline and GetConsensusPolyline to var-select
    // appropriate cost function / alignment functions.
    private static double CostFn(S2Polyline a, S2Polyline b, bool approx)
    {
        return approx ? GetApproxVertexAlignment(a, b).AlignmentCost
                      : GetExactVertexAlignmentCost(a, b);
    }

    // GetMedoidPolyline returns the index `p` of a "medoid" polyline from a
    // non-empty collection of `polylines` such that
    //
    // sum_{all j in `polylines`} VertexAlignmentCost(p, j) is minimized.
    //
    // In the case of a tie for minimal summed alignment cost, we return the lowest
    // index - this tie is guaranteed to happen in the two-polyline-input case.
    //
    // ASYMPTOTIC BEHAVIOR:
    // Computation may require up to (N^2 - N) / 2 alignment cost function
    // evaluations, for N input polylines. For polylines of length U, V, the
    // alignment cost function evaluation is O(U+V) if options.approx = true and
    // O(U*V) if options.approx = false.
    public class MedoidOptions
    {
        // If options.approx = false, we compute vertex alignment costs exactly.
        // If options.approx = true, we use approximate vertex alignment
        // computation, called with the default radius parameter.
        public bool Approx { get; set; } = true;
    }

    #endregion

    #region Consensus

    // Allocates and returns a new "consensus" polyline from a
    // non-empty collection of polylines. We iteratively apply Dynamic Timewarp
    // Barycenter Averaging to an initial `seed` polyline, which improves the
    // consensus alignment quality each iteration. For implementation details, see
    //
    // https://pdfs.semanticscholar.org/a596/8ca9488199291ffe5473643142862293d69d.pdf
    //
    // The returned polyline from this method is unlikely to be point-for-point
    // equal to an input polyline, whereas a polyline returned from
    // GetMedoidPolyline() is guaranteed to match an input polyline point-for-point.
    // NOTE: the number of points in our returned consensus polyline is always equal
    // to the number of points in the initial seed, which is implementation-defined.
    // If the collection of polylines has a large resolution distribution, it might
    // be a good idea to reinterpolate them to have about the same number of points.
    // In practice, this doesn't seem to matter, but is probably worth noting.
    //
    // ASYMPTOTIC BEHAVIOR:
    // Seeding this algorithm requires O(1) vertex alignments if seed_medoid =
    // false, and O(N^2) vertex alignments if seed_medoid = true. Once the seed
    // polyline is chosen, computing the consensus polyline requires at most
    // (iteration_cap)*N vertex alignments. For polylines of length U, V, the
    // alignment cost function evaluation is O(U+V) if options.approx = true, and
    // O(U*V) if options.approx = false.
    public static S2Polyline GetConsensusPolyline(S2Polyline[] polylines, ConsensusOptions options)
    {
        // Implements Iterative Dynamic Timewarp Barycenter Averaging algorithm from
        //
        // https://pdfs.semanticscholar.org/a596/8ca9488199291ffe5473643142862293d69d.pdf
        //
        // Algorithm:
        // Initialize consensus sequence with either the medoid or an arbitrary
        // element (chosen here to be the first element in the input collection).
        // While the consensus polyline `consensus` hasn't converged and we haven't
        // exceeded our iteration cap:
        //   For each polyline `p` in the input,
        //     Compute vertex alignment from the current consensus to `p`.
        //     For each (c_index, p_index) pair in the warp path,
        //       Add the S2Point pts.vertex(p_index) to S2Point consensus[c_index]
        //   Normalize (compute centroid) of each consensus point.
        //   Determine if consensus is converging; if no vertex has moved or we've hit
        //   the iteration cap, halt.
        //
        //  This algorithm takes O(iteration_cap * num_polylines) pairwise alignments.

        int num_polylines = polylines.Length;
        System.Diagnostics.Debug.Assert(num_polylines > 0);
        bool approx = options.Approx;

        // Seed a consensus polyline, either arbitrarily with first element, or with
        // the medoid. If seeding with medoid, inherit approx parameter from options.
        int seed_index = 0;
        if (options.SeedMedoid)
        {
            var medoid_options = new MedoidOptions { Approx = approx };
            seed_index = GetMedoidPolyline(polylines, medoid_options);
        }
        var consensus = (S2Polyline)polylines[seed_index].CustomClone();
        int num_consensus_vertices = consensus.NumVertices();
        System.Diagnostics.Debug.Assert(num_consensus_vertices > 1);

        var converged = false;
        int iterations = 0;
        while (!converged && iterations < options.IterationCap)
        {
            var points = new Dictionary<int, S2Point>();
            foreach (var polyline in polylines)
            {
                var alignment = AlignmentFn(consensus, polyline, approx);
                foreach (var (start, end) in alignment.WarpPath_)
                {
                    if (points.ContainsKey(start))
                    {
                        points[start] = points[start] + polyline.Vertex(end);
                    }
                    else
                    {
                        points.Add(start, polyline.Vertex(end));
                    }
                }
            }

            ++iterations;
            var new_consensus = new S2Polyline(points
                .OrderBy(t => t.Key)
                .Select(t => t.Value.Normalize())
                .ToArray());
            converged = new_consensus.ApproxEquals(consensus);
            consensus = new_consensus;
        }
        return consensus;
    }

    private static VertexAlignment AlignmentFn(S2Polyline a, S2Polyline b, bool approx)
    {
        return approx ? GetApproxVertexAlignment(a, b)
                      : GetExactVertexAlignment(a, b);
    }

    public class ConsensusOptions
    {
        // If options.approx = false, vertex alignments are computed with
        // GetExactVertexAlignment. If options.approx = true, vertex alignments are
        // computed with GetApproxVertexAlignment, called with default radius
        // parameter.
        public bool Approx { get; set; } = true;

        // If options.seed_medoid = true, we seed the consensus polyline with the
        // medoid of the collection. This is a more expensive approach, but may result
        // in higher quality consensus sequences by avoiding bad arbitrary initial
        // seeds. Seeding with the medoid will incur up to (N^2 - N) / 2 evaluations
        // of the vertex alignment function. If options.seed_medoid = false, we seed
        // the consensus polyline by taking an arbitrary element from the collection.
        public bool SeedMedoid { get; set; } = false;

        // options.iteration_cap controls the maximum number of DBA refining steps we
        // apply to the initial seed.
        public int IterationCap { get; set; } = 5;
    }

    #endregion

    // A Window is a sparse binary matrix with specific structuralraints
    // on allowed element configurations. It is used in this library to represent
    // "search windows" for windowed dynamic timewarping.
    //
    // Valid Windows require the following structural conditions to hold:
    // 1) All rows must consist of a single contiguous stride of `true` values.
    // 2) All strides are greater than zero length (i.e. no empty rows).
    // 3) The index of the first `true` column in a row must be at least as
    //    large as the index of the first `true` column in the previous row.
    // 4) The index of the last `true` column in a row must be at least as large
    //    as the index of the last `true` column in the previous row.
    // 5) strides[0].start = 0 (the first cell is always filled).
    // 6) strides[n_rows-1].end = n_cols (the last cell is filled).
    //
    // Example valid strided_masks (* = filled, . = unfilled)
    //   0 1 2 3 4 5
    // 0 * * * . . .
    // 1 . * * * . .
    // 2 . * * * . .
    // 3 . . * * * *
    // 4 . . * * * *
    //   0 1 2 3 4 5
    // 0 * * * * . .
    // 1 . * * * * .
    // 2 . . * * * .
    // 3 . . . . * *
    // 4 . . . . . *
    //   0 1 2 3 4 5
    // 0 * * . . . .
    // 1 . * . . . .
    // 2 . . * * * .
    // 3 . . . . . *
    // 4 . . . . . *
    //
    // Example invalid strided_masks:
    //   0 1 2 3 4 5
    // 0 * * * . * * <-- more than one continuous run
    // 1 . * * * . .
    // 2 . * * * . .
    // 3 . . * * * *
    // 4 . . * * * *
    //   0 1 2 3 4 5
    // 0 * * * . . .
    // 1 . * * * . .
    // 2 . * * * . .
    // 3 * * * * * * <-- start index not monotonically increasing
    // 4 . . * * * *
    //   0 1 2 3 4 5
    // 0 * * * . . .
    // 1 . * * * * .
    // 2 . * * * . . <-- end index not monotonically increasing
    // 3 . . * * * *
    // 4 . . * * * *
    //   0 1 2 3 4 5
    // 0 . * . . . . <-- does not fill upper left corner
    // 1 . * . . . .
    // 2 . * . . . .
    // 3 . * * * . .
    // 4 . . * * * *
    public class Window
    {
        // Construct a Window from a non-empty list of column strides.
        public Window((int start, int end)[] strides)
        {
            System.Diagnostics.Debug.Assert(strides.Any()); // Cannotruct empty window.
            System.Diagnostics.Debug.Assert(strides[0].start == 0); // First element of start_cols is non-zero.
            strides_ = strides;
            rows_ = strides.Length;
            cols_ = strides.Last().end;
            System.Diagnostics.Debug.Assert(IsValid); // Constructor validity check fail.
        }

        // Construct a Window from a non-empty sequence of warp path index pairs.
        public Window(WarpPath warp_path)
        {
            System.Diagnostics.Debug.Assert(warp_path.Any()); // Cannot construct window from empty warp path.
            System.Diagnostics.Debug.Assert(warp_path.First() == (0, 0)); // Must start at (0, 0).
            rows_ = warp_path.Last().start + 1;
            System.Diagnostics.Debug.Assert(rows_ > 0); // Must have at least one row.
            cols_ = warp_path.Last().end + 1;
            System.Diagnostics.Debug.Assert(cols_ > 0); // Must have at least one column.
            strides_ = new (int start, int end)[rows_];

            int prev_row = 0;
            int curr_row = 0;
            int stride_start = 0;
            int stride_stop = 0;
            foreach (var (start, end) in warp_path)
            {
                curr_row = start;
                if (curr_row > prev_row)
                {
                    strides_[prev_row] = (stride_start, stride_stop);
                    stride_start = end;
                    prev_row = curr_row;
                }
                stride_stop = end + 1;
            }
            System.Diagnostics.Debug.Assert(curr_row == rows_ - 1);
            strides_[rows_ - 1] = (stride_start, stride_stop);
            System.Diagnostics.Debug.Assert(IsValid); // Constructor validity check fail.
        }

        // Return the (not-bounds-checked) stride for this row.
        public (int start, int end) GetColumnStride(int row)
        {
            return strides_[row];
        }

        // Return the (bounds-checked) stride for this row.
        // If row < 0, returns ColumnStride.All()
        public (int start, int end) GetCheckedColumnStride(int row)
        {
            return (row > 0) ? strides_[row] : ColumnStrideAll();
        }

        // Return a new, larger Window that is an upscaled version of this window
        // Used by ApproximateAlignment window expansion step.
        public Window Upsample(int new_rows, int new_cols)
        {
            System.Diagnostics.Debug.Assert(new_rows >= rows_); // Upsampling: New_rows < current_rows
            System.Diagnostics.Debug.Assert(new_cols >= cols_); // Upsampling: New_cols < current_cols
            var row_scale = ((double)new_rows) / rows_;
            var col_scale = ((double)new_cols) / cols_;
            var new_strides = new (int start, int end)[new_rows];
            for (int row = 0; row < new_rows; ++row)
            {
                var (start, end) = strides_[(int)((row + 0.5) / row_scale)];
                new_strides[row] = (
                    start: (int)(col_scale * start + 0.5),
                    end: (int)(col_scale * end + 0.5));
            }
            return new Window(new_strides);
        }

        // Return a new, equal-size Window by dilating this window with a square
        // structuring element with half-length `radius`. Radius = 1 corresponds to
        // a 3x3 square morphological dilation.
        // Used by ApproximateAlignment window expansion step.
        public Window Dilate(int radius)
        {
            // This code takes advantage of the fact that the dilation window is square to
            // ensure that we can compute the stride for each output row in constant time.
            // TODO (mrdmnd): a potential optimization might be to combine this method and
            // the Upsample method into a single "Expand" method. For the sake of
            // testing, I haven't done that here, but I think it would be fairly
            // straightforward to do so. This method generally isn't very expensive so it
            // feels unnecessary to combine them.

            System.Diagnostics.Debug.Assert(radius >= 0); // Negative dilation radius.
            var new_strides = new (int start, int end)[rows_];
            int prev_row, next_row;
            for (int row = 0; row < rows_; ++row)
            {
                prev_row = Math.Max(0, row - radius);
                next_row = Math.Min(row + radius, rows_ - 1);
                new_strides[row] = (
                    start: Math.Max(0, strides_[prev_row].start - radius),
                    end: Math.Min(strides_[next_row].end + radius, cols_));
            }
            return new Window(new_strides);
        }

        // Return a string representation of this window.
        // implemented primarily for testing purposes.
        public override string ToString()
        {
            var sb = new StringBuilder();
            for (int row = 0; row < rows_; ++row)
            {
                for (int col = 0; col < cols_; ++col)
                {
                    sb.Append(InRange(col, strides_[row]) ? " *" : " .");
                }
                sb.AppendLine();
            }
            return sb.ToString();
        }

        // Valid Windows require the following structural conditions to hold:
        // 1) All rows must consist of a single contiguous stride of `true` values.
        // 2) All strides are greater than zero length (i.e. no empty rows).
        // 3) The index of the first `true` column in a row must be at least as
        //    large as the index of the first `true` column in the previous row.
        // 4) The index of the last `true` column in a row must be at least as large
        //    as the index of the last `true` column in the previous row.
        // 5) strides[0].start = 0 (the first cell is always filled).
        // 6) strides[n_rows-1].end = n_cols (the last cell is filled).
        //
        // Returns true if this window's data represents a valid window.
        private bool IsValid
        {
            get
            {
                if (rows_ <= 0 || cols_ <= 0 || strides_.First().start != 0 ||
                    strides_.Last().end != cols_)
                {
                    return false;
                }

                var prev = (start: -1, end: -1);
                foreach (var curr in strides_)
                {
                    if (curr.end <= curr.start || curr.start < prev.start ||
                        curr.end < prev.end)
                    {
                        return false;
                    }
                    prev = curr;
                }
                return true;
            }
        }

        private readonly int rows_;
        private readonly int cols_;
        private readonly (int start, int end)[] strides_;
    }
}

