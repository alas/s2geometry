// On a sphere, polygon hierarchy is ambiguous.  If you imagine two chains
// around a sphere at +/- 10 degrees latitude, then we're completely justified
// to consider either one a shell, with the other being a hole in that shell.
//
// Given this ambiguity, we need to specify a strategy to choose a chain to be a
// shell by definition (which we call the datum shell). Other chains are then
// classified relative to that chain.  Often, the first chain in a shape will be
// a shell by construction, so use that as our default strategy.
//
// To specify the datum chain, we just need a function that takes an S2Shape and
// produces a chain id to use.  We use a function pointer to avoid the copying
// issues associated with std::function.  Non-capturing lambdas may still be
// used with a function pointer, but this may change to an `AnyInvocable` in the
// future.

namespace S2Geometry;

using S2DatumStrategy = Func<S2Shape, System.Int32>;

// `S2ShapeNestingQuery` defines a query to determine the relationships between
// chains in a shape.  Chains are either shells or holes.  Shells have no parent
// and may have zero or more holes associated with them.  Holes belong to a
// single parent shell and have zero holes of their own.  Polygon interiors are
// always on the left of their boundaries, so shells and holes face each other
// in the sense of containing the interior between them.
//
// It's always possible to reverse how we interpret shells and holes due to the
// geometry being on a sphere; so to construct a consistent relationship we need
// to have a strategy to identify a chain to use as a 'datum' shell to which all
// the other chains are classified relatively.  The strategy used for this is
// configurable through the `S2ShapeNestingQueryOptions` and defaults to
// FIRST_CHAIN, selecting the first chain in the shape as the datum.
//
// The `ShapeNesting`() function determines the relationship between shells and
// returns a vector of `ChainRelation` instances.  The `ChainRelation`s are in
// 1:1 correspondence with the chain_id in the shape, e.g. chain 1's relations
// are always given by `ShapeNesting()[1]`.
//
// In the common case of a shell having one or no holes, we store the hole
// relationships efficiently so that no memory allocation is required when
// creating a `ChainRelation`.  Complex geometry may exceed those limits and
// allocate however.
//
// Restrictions:
//   `S2ShapeNestingQuery` is only meaningful for S2ShapeIndex instances
//   containing 2D geometry.
//
//   The query currently doesn't handle any sort of degeneracy in the underlying
//   geometry.
//
public class S2ShapeNestingQuery
{
    public S2ShapeIndex Index { get; set; }
    public Options Options_ { get; set; }

    public S2ShapeNestingQuery(S2ShapeIndex index, Options? options = null)
    {
        Index = index;
        Options_ = options ?? new();
    }

    // Evaluates the relationships between chains in a shape.  This determines
    // which chains are shells, and which are holes associated with a parent
    // shell.
    //
    // The returned `ChainRelation` instances are in 1:1 correspondence with the
    // chains in the shape, i.e. chain id 3 responds to `result[3]`.
    //
    // Finds three consecutive vertices that aren't degenerate and stores them in
    // the given array.  Returns false if there aren't three consecutive
    // non-degenerate points.
    public List<ChainRelation> ComputeShapeNesting(int shape_id)
    {
        var shape = Index.Shape(shape_id);
        if (shape is null || shape.NumChains() == 0)
        {
            return new();
        }
        MyDebug.Assert(shape.Dimension() == 2);

        int num_chains = shape.NumChains();

        // A single chain is always a shell, with no holes.
        if (num_chains == 1)
        {
            return new() { ChainRelation.MakeShell() };
        }

        // Sets to track possible parents and possible children for each chain.
        var parents = new Bitmap64[num_chains].Fill(() => new Bitmap64(num_chains, false));
        var children = new Bitmap64[num_chains].Fill(() => new Bitmap64(num_chains, false));

        // We'll compute edge crossings along a line segment from the datum shell to a
        // random point on the other chains.  This choice is arbitrary, so we'll use
        // the first vertex of edge 1 so we can easily get the next and previous
        // points to check for orientation.
        var datum_shell = Options_.DatumStrategy(shape);
        var vertices = new S2Point[]
        {
            shape.ChainEdge(datum_shell, 0).V0,
            shape.ChainEdge(datum_shell, 1).V0,
            shape.ChainEdge(datum_shell, 2).V0,
        };
        var start_point = vertices[1];

        S2CrossingEdgeQuery crossing_query = new(Index);
        List<S2ShapeUtil.ShapeEdge> edges=new();
        for (int chain = 0; chain < num_chains; ++chain)
        {
            if (chain == datum_shell)
            {
                continue;
            }
            //S2_VLOG(1) << "Processing chain " << chain;

            // Find a close point on the target chain out of 4 equally spaced ones.
            int end_idx = ClosestOfNPoints(start_point, shape, chain, 4);
            S2Point end_point = shape.ChainEdge(chain, end_idx).V0;

            // We need to know whether we're inside the datum shell at the end, so we
            // need to properly seed its starting state.  If we start by entering the
            // datum shell's interior _and_ end by arriving from the target chain's
            // interior, we set it to true.
            //
            // As we cross edges from the datum to the target chain the total number of
            // datum shell _or_ target chain edges we'll cross is either even or odd.
            // Each of these edges toggles our "insideness" relative to the datum shell,
            if (S2Pred.OrderedCCW(vertices[2], end_point, vertices[0], start_point))
            {
                //S2_VLOG(1) << "  Edge starts into interior of datum chain";
                parents[chain].Set(datum_shell, true);
                children[datum_shell].Set(chain, true);
            }

            // Arriving from the interior of the target chain?
            S2Point next = NextChainEdge(shape, chain, end_idx).V0;
            S2Point prev = PrevChainEdge(shape, chain, end_idx).V0;
            if (S2Pred.OrderedCCW(next, start_point, prev, end_point))
            {
                //S2_VLOG(1) << "  Edge ends from interior of target chain";
                parents[chain].Set(chain, true);
            }
            //S2_VLOG(2) << "    Initial set: " << parents[chain].ToString(8);

            // Query all the edges crossed by the line from the datum shell to a point
            // on this chain.  Only look at edges that belong to the requested shape.
            // Using INTERIOR here will avoid returning the two edges on the datum and
            // target shells that are touched by the endpoints of our line segment.
            crossing_query.GetCrossingEdges(start_point, end_point, shape,
                                            CrossingType.INTERIOR,
                                            edges);

            // Walk through the intersected chains and toggle corresponding bits.
            foreach (var edge in edges)
            {
                var other_chain = shape.GetChainPosition(edge.Id.EdgeId).ChainId;

                parents[chain].Toggle(other_chain);
                if (other_chain != chain)
                {
                    children[other_chain].Toggle(chain);
                }
                //S2_VLOG(1) << "  Crosses chain " << other_chain;
                //S2_VLOG(2) << "    Parent set: " << parents[chain].ToString(8);
            }

            // Now set the final state.  Remove the target chain from its own parent set
            // to make following logic simpler.  The datum shell is a potential parent
            // if both parent and target chain bits are set.
            parents[chain].Set(datum_shell, parents[chain].Get(datum_shell) &&
                                                parents[chain].Get(chain));
            parents[chain].Set(chain, false);
        }

        /*if (S2_VLOG_IS_ON(2))
        {
            S2_LOG(INFO) << "Current parent set";
            for (int chain = 0; chain<num_chains; ++chain)
            {
                S2_LOG(INFO) << "  " << absl::StrFormat("%2d", chain) << ": "
                          << parents[chain].ToString(8);
        }
        }*/

        // Look at each chain with a single parent and remove the parent from any of
        // its child chains.  This enforces the constraint that if A is a parent of B
        // and B is a parent of C, then A shouldn't directly be a parent of C.
        for (int current_chain = 0; current_chain < num_chains; ++current_chain)
        {
            if (parents[current_chain].GetOnesCount() != 1)
            {
                continue;
            }

            int parent_chain;
            parents[current_chain].FindFirstSetBit(out parent_chain);

            int next_chain = current_chain;
            int child = 0;
            for (; children[current_chain].FindNextSetBit(ref child); child++)
            {
                if (parents[child].Get(parent_chain))
                {
                    parents[child].Set(parent_chain, false);

                    // If this chain has a single parent now, we have to process it as well,
                    // so if we've already passed it in the outer loop, we have to back up.
                    if (parents[child].GetOnesCount() == 1 && child < next_chain)
                    {
                        next_chain = child;
                    }
                }
            }

            /*S2_VLOG(1) << "Chain " << current_chain << " has one parent";
            if (S2_VLOG_IS_ON(2))
            {
                S2_LOG(INFO) << "  Parent set now:";
                for (int chain = 0; chain < num_chains; ++chain)
                {
                    S2_LOG(INFO) << "  " << absl::StrFormat("%2d", chain) << ": "
                              << parents[chain].ToString(8);
                }
            }*/

            // Backup current chain so next loop increment sets it properly
            if (next_chain != current_chain)
            {
                current_chain = next_chain - 1;
            }
        }

        // Each chain now points to its immediate parent.  Scan through and set child
        // to point to parent and vice-versa.
        var relations = new List<ChainRelation>().Fill(() => new(), num_chains);
        for (int chain = 0; chain < num_chains; ++chain)
        {
            MyDebug.Assert(parents[chain].GetOnesCount() <= 1);

            if (parents[chain].FindFirstSetBit(out int parent))
            {
                relations[chain].Parent = parent;
                relations[parent].AddHole(chain);
            }
        }

        // Detach any chains that are even depth from their parent and make them
        // shells.  This is effectively implementing the even/odd rule.
        for (int chain = 0; chain < num_chains; ++chain)
        {
            int depth = -1, current = chain;
            do
            {
                ++depth;
                current = relations[current].Parent;
            } while (current >= 0 && depth < num_chains);
            MyDebug.Assert(depth < num_chains);

            if (depth != 0 && (depth % 2 == 0))
            {
                relations[chain].ClearParent();
            }
        }

        return relations;
    }

    // Takes N equally spaced points from the given chain of the shape and finds
    // the one closest to the target point, returning its index.
    private static int ClosestOfNPoints(S2Point target, S2Shape shape,
                                   int chain, int num_points)
    {
        int chain_len = shape.GetChain(chain).Length;

        // If chain_len is < num_points, we still want to use whatever points there
        // are in the chain, so max 1 the minimum step size and we'll modulo to get
        // back into bounds for the chain.
        int step = Math.Max(1, chain_len / num_points);

        double min_dist2 = double.PositiveInfinity;
        int closest_idx = 0;
        for (int i = 0; i < num_points; ++i)
        {
            int idx = (i * step) % chain_len;
            S2Point point = shape.ChainEdge(chain, idx).V0;
            double dist2 = (target - point).Norm2();
            if (dist2 < min_dist2)
            {
                min_dist2 = dist2;
                closest_idx = idx;
            }
        }
        return closest_idx;
    }

    // Returns the next edge of a particular chain, handling index wrap around.
    private static S2Shape.Edge NextChainEdge(S2Shape shape, int chain, int edge)
    {
        return shape.ChainEdge(chain, (edge + 1) % shape.GetChain(chain).Length);
    }

    // Returns the previous edge of a particular chain, handling index wrap around.
    private static S2Shape.Edge PrevChainEdge(S2Shape shape, int chain, int edge)
    {
        int index = edge - 1;
        if (index < 0)
        {
            index = shape.GetChain(chain).Length - 1;
        }
        return shape.ChainEdge(chain, index);
    }

    // Class to model options for the query, passed to `Init()`.
    public class Options
    {
        public S2DatumStrategy DatumStrategy { get; set; }

        public Options() { DatumStrategy = GetFirstChain; }

        // Returns the first chain id in a shape (always zero), used as the
        // default datum strategy.
        public static int GetFirstChain(S2Shape shape) { return 0; }
    }

    // `ChainRelation` models the parent/child relationship between chains in a
    // shape and chain classification as a shell or a hole.
    public class ChainRelation
    {
        // So that the query can add Holes to the chain before returning it.
        //friend class S2ShapeNestingQuery;

        // `absl::InlinedVector` doesn't round up its inline capacity even if it
        // wouldn't make the object any bigger.  It stores a pointer and a capacity
        // with the inlined data in a union:
        //
        //   [ pointer | capacity ]
        //   [        T[N]        ]
        //
        // We pay for the size of the pointer and capacity regardless, so it makes
        // sense for us to scale the number of reserved elements to take advantage
        // of this.  On 64-bit systems we can store 4 int32s and on 32-bit we can
        // store 2.

        private readonly ArrayN<Int32> holes_;

        // Returns the id of the parent chain of the associated chain.  Chains that
        // are shells don't have a parent and have a parent id of -1.
        public Int32 Parent { get; set; }

        public ChainRelation(Int32 parent = -1)
        {
            Parent=parent;
            if (IntPtr.Size == 4)
            {
                //32 bits
                holes_ = new Array2<Int32>();
            }
            else if (IntPtr.Size == 8)
            {
                //64 bits
                holes_ = new Array4<Int32>();
            }
            else
            {
                throw new NotImplementedException("only 32/64 bits supported");
            }
        }

        // Builds a `ChainRelation` that's a shell with given holes.
        public static ChainRelation MakeShell(List<Int32>? holes = null)
        {
            ChainRelation relation = new();
            holes ??= new();
            foreach (var chain in holes)
            {
                relation.AddHole(chain);
            }
            return relation;
        }

        // Builds a `ChainRelation` that's a hole with given parent.
        public static ChainRelation MakeHole(Int32 parent)
        {
            return new ChainRelation(parent);
        }

        // Returns true if the associated chain is a shell.  Otherwise the chain
        // is a hole in its parent chain.
        public bool IsShell() { return Parent < 0; }

        // Returns true if the associated chain is a hole.  This is true iff the
        // chain is not a shell.
        public bool IsHole() { return !IsShell(); }

        // Return number of holes in the associated chain.
        public int GetNumHoles() { return holes_.Count; }

        // Returns a read only view over the hole ids for the associated chain.
        public List<Int32> GetHoles() { return new List<Int32>(holes_); }

        // Adds the given id as a hole of this chain.  The `ParentId` of the hole
        // should therefore be the id of this chain as an invariant.
        internal void AddHole(Int32 id) { holes_.Add(id); }

        // Clears the parent of the associated chain (thus marking it as a shell).
        public void ClearParent() { Parent = -1; }
    }
}
