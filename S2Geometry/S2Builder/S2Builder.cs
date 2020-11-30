// Internal flag intended to be set from within a debugger.
// #define s2builder_verbose = false;

using System;
using System.Collections.Generic;
using System.Linq;
using S2Geometry.S2BuilderUtil;
using S2Geometry.S2ShapeUtil;

#if s2builder_verbose
using System.Text;
#endif

namespace S2Geometry
{
    // These types define the output vertices and edges.
    using OutputEdge = KeyKey<int, int>;            // Defines an output edge.   
    using LayerEdgeId = KeyKey<int, int>;           // Identifies an output edge in a particular layer.   
    using InputEdge = KeyKey<int, int>;             // Defines an input edge.
    using InputVertexKey = KeyKey<S2CellId, int>;   // Sort key for prioritizing input vertices.

    // S2Builder is a tool for assembling polygonal geometry from edges.  Here are
    // some of the things it is designed for:
    //
    // 1. Building polygons, polylines, and polygon meshes from unsorted
    //    collections of edges.
    //
    // 2. Snapping geometry to discrete representations (such as S2CellId centers
    //    or E7 lat/lng coordinates) while preserving the input topology and with
    //    guaranteed error bounds.
    //
    // 3. Simplifying geometry (e.g. for indexing, display, or storage).
    //
    // 4. Importing geometry from other formats, including repairing geometry
    //    that has errors.
    //
    // 5. As a tool for implementing more complex operations such as polygon
    //    intersections and unions.
    //
    // The implementation is based on the framework of "snap rounding".  Unlike
    // most snap rounding implementations, S2Builder defines edges as geodesics on
    // the sphere (straight lines) and uses the topology of the sphere (i.e.,
    // there are no "seams" at the poles or 180th meridian).  The algorithm is
    // designed to be 100% robust for arbitrary input geometry.  It offers the
    // following properties:
    //
    //   - Guaranteed bounds on how far input vertices and edges can move during
    //     the snapping process (i.e., at most the given "snap_radius").
    //
    //   - Guaranteed minimum separation between edges and vertices other than
    //     their endpoints (similar to the goals of Iterated Snap Rounding).  In
    //     other words, edges that do not intersect in the output are guaranteed
    //     to have a minimum separation between them.
    //
    //   - Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the
    //     input already meets the output criteria then it will not be modified.
    //
    //   - Preservation of the input topology (up to the creation of
    //     degeneracies).  This means that there exists a continuous deformation
    //     from the input to the output such that no vertex crosses an edge.  In
    //     other words, self-intersections won't be created, loops won't change
    //     orientation, etc.
    //
    //   - The ability to snap to arbitrary discrete point sets (such as S2CellId
    //     centers, E7 lat/lng points on the sphere, or simply a subset of the
    //     input vertices), rather than being limited to an integer grid.
    //
    // Here are some of its other features:
    //
    //  - It can handle both directed and undirected edges.  Undirected edges can
    //    be useful for importing data from other formats, e.g. where loops have
    //    unspecified orientations.
    //
    //  - It can eliminate self-intersections by finding all edge pairs that cross
    //    and adding a new vertex at each intersection point.
    //
    //  - It can simplify polygons to within a specified tolerance.  For example,
    //    if two vertices are close enough they will be merged, and if an edge
    //    passes nearby a vertex then it will be rerouted through that vertex.
    //    Optionally, it can also detect nearly straight chains of short edges and
    //    replace them with a single long edge, while maintaining the same
    //    accuracy, separation, and topology guarantees ("simplify_edge_chains").
    //
    //  - It supports many different output types through the concept of "layers"
    //    (polylines, polygons, polygon meshes, etc).  You can build multiple
    //    layers at once in order to ensure that snapping does not create
    //    intersections between different objects (for example, you can simplify a
    //    set of contour lines without the risk of having them cross each other).
    //
    //  - It supports edge labels, which allow you to attach arbitrary information
    //    to edges and have it preserved during the snapping process.  (This can
    //    also be achieved using layers, at a coarser level of granularity.)
    //
    // Caveats:
    //
    //  - Because S2Builder only works with edges, it cannot distinguish between
    //    the empty and full polygons.  If your application can generate both the
    //    empty and full polygons, you must implement logic outside of this class.
    //
    // Example showing how to snap a polygon to E7 coordinates:
    //
    //  using s2builderutil.IntLatLngSnapFunction;
    //  S2Builder builder(S2Builder.Options(IntLatLngSnapFunction(7)));
    //  S2Polygon output;
    //  builder.StartLayer(new s2builderutil.S2PolygonLayer(output));
    //  builder.AddPolygon(input);
    //  S2Error error;
    //  if (!builder.Build(out error)) {
    //    S2_LOG(ERROR) << error;
    //    ...
    //  }
    //
    // The algorithm is based on the idea of choosing a set of sites and computing
    // their "limited radius Voronoi diagram", which is obtained by intersecting
    // each Voronoi region with a disc of fixed radius (the "snap radius")
    // centered around the corresponding site.
    //
    // For each input edge, we then determine the sequence of Voronoi regions
    // crossed by that edge, and snap the edge to the corresponding sequence of
    // sites.  (In other words, each input edge is replaced by an edge chain.)
    //
    // The sites are chosen by starting with the set of input vertices, optionally
    // snapping them to discrete point set (such as S2CellId centers or lat/lng E7
    // coordinates), and then choosing a subset such that no two sites are closer
    // than the given "snap_radius".  Note that the sites do not need to be spaced
    // regularly -- their positions are completely arbitrary.
    //
    // Rather than computing the full limited radius Voronoi diagram, instead we
    // compute on demand the sequence of Voronoi regions intersected by each edge.
    // We do this by first finding all the sites that are within "snap_radius" of
    // the edge, sorting them by distance from the edge origin, and then using an
    // incremental algorithm.
    //
    // We implement the minimum edge-vertex separation property by snapping all
    // the input edges and checking whether any site (the "site to avoid") would
    // then be too close to the edge.  If so we add another site (the "separation
    // site") along the input edge, positioned so that the new snapped edge is
    // guaranteed to be far enough away from the site to avoid.  We then find all
    // the input edges that are within "snap_radius" of the new site, and resnap
    // those edges.  (It is very rare that we need to add a separation site, even
    // when sites are densely spaced.)
    //
    // Idempotency is implemented by explicitly checking whether the input
    // geometry already meets the output criteria.  This is not as sophisticated
    // as Stable Snap Rounding (Hershberger); I have worked out the details and it
    // is possible to adapt those ideas here, but it would make the implementation
    // significantly more complex.
    //
    // The only way that different output layers interact is in the choice of
    // Voronoi sites:
    //
    //  - Vertices from all layers contribute to the initial selection of sites.
    //
    //  - Edges in any layer that pass too close to a site can cause new sites to
    //    be added (which affects snapping in all layers).
    //
    //  - Simplification can be thought of as removing sites.  A site can be
    //    removed only if the snapped edges stay within the error bounds of the
    //    corresponding input edges in all layers.
    //
    // Otherwise all layers are processed independently.  For example, sibling
    // edge pairs can only cancel each other within a single layer (if desired).
    public partial class S2Builder
    {
        public class Options
        {
            public Options() => SnapFunction = new IdentitySnapFunction(S1Angle.Zero);

            // Convenience constructor that calls set_snap_function().
            public Options(SnapFunction snapFunction) => SnapFunction = (SnapFunction)snapFunction.Clone();

            public Options(Options options)
            {
                SnapFunction = (SnapFunction)options.SnapFunction.Clone();
                SplitCrossingEdges = options.SplitCrossingEdges;
                SimplifyEdgeChains = options.SimplifyEdgeChains;
                Idempotent = options.Idempotent;
            }

            // Sets the desired snap function.  The snap function is copied
            // internally, so you can safely pass a temporary object.
            //
            // Note that if your input data includes vertices that were created using
            // S2EdgeCrossings.GetIntersection(), then you should use a "snap_radius" of
            // at least S2.kIntersectionSnapRadius, e.g. by calling
            //
            //  options.set_snap_function(s2builderutil.IdentitySnapFunction(
            //      S2.kIntersectionSnapRadius));
            //
            // DEFAULT: s2builderutil.IdentitySnapFunction(S1Angle.Zero())
            // [This does no snapping and preserves all input vertices exactly.]
            public SnapFunction SnapFunction { get; set; }

            // If true, then detect all pairs of crossing edges and eliminate them by
            // adding a new vertex at their intersection point.
            //
            // When this option is true, the effective snap_radius() for edges is
            // increased by S2EdgeCrossings.kIntersectionError to take into account the
            // additional error when computing intersection points.  In other words,
            // edges may move by up to snap_radius() + S2EdgeCrossings.kIntersectionError.
            //
            // Undirected edges should always be used when the output is a polygon,
            // since splitting a directed loop at a self-intersection converts it into
            // two loops that don't define a consistent interior according to the
            // "interior is on the left" rule.  (On the other hand, it is fine to use
            // directed edges when defining a polygon *mesh* because in that case the
            // input consists of sibling edge pairs.)
            //
            // Self-intersections can also arise when importing data from a 2D
            // projection.  You can minimize this problem by subdividing the input
            // edges so that the S2 edges (which are geodesics) stay close to the
            // original projected edges (which are curves on the sphere).  This can
            // be done using s2builderutil.EdgeSplitter(), for example.
            //
            // DEFAULT: false
            public bool SplitCrossingEdges { get; set; } = false;

            // If true, then simplify the output geometry by replacing nearly straight
            // chains of short edges with a single long edge.
            //
            // The combined effect of snapping and simplifying will not change the
            // input by more than the guaranteed tolerances (see the list documented
            // with the SnapFunction class).  For example, simplified edges are
            // guaranteed to pass within snap_radius() of the *original* positions of
            // all vertices that were removed from that edge.  This is a much tighter
            // guarantee than can be achieved by snapping and simplifying separately.
            //
            // However, note that this option does not guarantee idempotency.  In
            // other words, simplifying geometry that has already been simplified once
            // may simplify it further.  (This is unavoidable, since tolerances are
            // measured with respect to the original geometry, which is no longer
            // available when the geometry is simplified a second time.)
            //
            // When the output consists of multiple layers, simplification is
            // guaranteed to be consistent: for example, edge chains are simplified in
            // the same way across layers, and simplification preserves topological
            // relationships between layers (e.g., no crossing edges will be created).
            // Note that edge chains in different layers do not need to be identical
            // (or even have the same number of vertices, etc) in order to be
            // simplified together.  All that is required is that they are close
            // enough together so that the same simplified edge can meet all of their
            // individual snapping guarantees.
            //
            // Note that edge chains are approximated as parametric curves rather than
            // point sets.  This means that if an edge chain backtracks on itself (for
            // example, ABCDEFEDCDEFGH) then such backtracking will be preserved to
            // within snap_radius() (for example, if the preceding point were all in a
            // straight line then the edge chain would be simplified to ACFCFH, noting
            // that C and F have degree > 2 and therefore can't be simplified away).
            //
            // Simplified edges are assigned all labels associated with the edges of
            // the simplified chain.
            //
            // For this option to have any effect, a SnapFunction with a non-zero
            // snap_radius() must be specified.  Also note that vertices specified
            // using ForceVertex are never simplified away.
            //
            // DEFAULT: false
            public bool SimplifyEdgeChains
            {
                get => _simplifyEdgeChains;
                set
                {
                    _simplifyEdgeChains = value;

                    // Simplification requires a non-zero snap radius, and while it might be
                    // possible to do some simplifying without snapping, it is much simpler to
                    // always snap (even if the input geometry already meets the other output
                    // requirements).  We need to compute edge_sites_ in order to avoid
                    // approaching non-incident vertices too closely, for example.
                    Idempotent = false;
                }
            }
            private bool _simplifyEdgeChains = false;

            // If true, then snapping occurs only when the input geometry does not
            // already meet the S2Builder output guarantees (see the SnapFunction
            // class description for details).  This means that if all input vertices
            // are at snapped locations, all vertex pairs are separated by at least
            // min_vertex_separation(), and all edge-vertex pairs are separated by at
            // least min_edge_vertex_separation(), then no snapping is done.
            //
            // If false, then all vertex pairs and edge-vertex pairs closer than
            // "snap_radius" will be considered for snapping.  This can be useful, for
            // example, if you know that your geometry contains errors and you want to
            // make sure that features closer together than "snap_radius" are merged.
            //
            // This option is automatically turned off by simplify_edge_chains(),
            // since simplifying edge chains is never guaranteed to be idempotent.
            //
            // DEFAULT: true
            public bool Idempotent { get; set; } = true;
        }

        // An S2Shape used to represent the entire collection of S2Builder input edges.
        // Vertices are specified as indices into a vertex vector to save space.
        public sealed class VertexIdEdgeVectorShape : S2Shape
        {
            // Requires that "edges" isant for the lifetime of this object.
            public VertexIdEdgeVectorShape(List<InputEdge> edges, List<S2Point> vertices)
            {
                edges_ = edges; vertices_ = vertices;
            }

            public S2Point Vertex0(int e) => Vertex(edges_[e].Item1);
            public S2Point Vertex1(int e) { return Vertex(edges_[e].Item2); }

            // S2Shape interface:
            public override int NumEdges => edges_.Count;
            public override Edge GetEdge(int e)
            {
                return new(vertices_[edges_[e].Item1], vertices_[edges_[e].Item2]);
            }
            public override int Dimension() { return 1; }
            public override ReferencePoint GetReferencePoint()
            {
                return ReferencePoint.FromContained(false);
            }
            public override int NumChains() { return edges_.Count; }
            public override Chain GetChain(int i) { return new Chain(i, 1); }
            public override Edge ChainEdge(int i, int j) { return GetEdge(i); }
            public override ChainPosition GetChainPosition(int e)
            {
                return new ChainPosition(e, 0);
            }

            private S2Point Vertex(int i) { return vertices_[i]; }

            private readonly List<InputEdge> edges_;
            private readonly List<S2Point> vertices_;
        }

        // A class that encapsulates the state needed for simplifying edge chains.
        public class EdgeChainSimplifier
        {
            // The graph "g" contains all edges from all layers.  "edge_layers"
            // indicates the original layer for each edge.  "site_vertices" is a map
            // from SiteId to the set of InputVertexIds that were snapped to that site.
            // "layer_edges" and "layer_input_edge_ids" are output arguments where the
            // simplified edge chains will be placed.  The input and output edges are
            // not sorted.
            public EdgeChainSimplifier(
              S2Builder builder, Graph g,
              List<int> edge_layers,
              List<List<int>> site_vertices,
              List<List<OutputEdge>> layer_edges,
              List<List<int>> layer_input_edge_ids,
              IdSetLexicon input_edge_id_set_lexicon)
            {
                builder_ = builder; g_ = g; in_ = new Graph.VertexInMap(g); out_ = new Graph.VertexOutMap(g); edge_layers_ = edge_layers;
                site_vertices_ = site_vertices; layer_edges_ = layer_edges;
                layer_input_edge_ids_ = layer_input_edge_ids;
                input_edge_id_set_lexicon_ = input_edge_id_set_lexicon;
                layer_begins_ = builder_.layer_begins_;
                is_interior_ = new List<bool>(g.NumVertices); used_ = new List<bool>(g.NumEdges);
                new_edges_.Capacity = g.NumEdges;
                new_input_edge_ids_.Capacity = g.NumEdges;
                new_edge_layers_.Capacity = g.NumEdges;
            }

            public void Run()
            {
                // Determine which vertices can be interior vertices of an edge chain.
                for (var v = 0; v < g_.NumVertices; ++v)
                {
                    is_interior_[v] = IsInterior(v);
                }
                // Attempt to simplify all edge chains that start from a non-interior
                // vertex.  (This takes care of all chains except loops.)
                for (var e = 0; e < g_.NumEdges; ++e)
                {
                    if (used_[e]) continue;
                    var edge = g_.GetEdge(e);
                    if (is_interior_[edge.Item1]) continue;
                    if (!is_interior_[edge.Item2])
                    {
                        OutputEdge(e);  // An edge between two non-interior vertices.
                    }
                    else
                    {
                        SimplifyChain(edge.Item1, edge.Item2);
                    }
                }
                // If there are any edges left, they form one or more disjoint loops where
                // all vertices are interior vertices.
                //
                // TODO(ericv): It would be better to start from the edge with the smallest
                // min_input_edge_id(), since that would make the output more predictable
                // for testing purposes.  It also means that we won't create an edge that
                // spans the start and end of a polyline if the polyline is snapped into a
                // loop.  (Unfortunately there are pathological examples that prevent us
                // from guaranteeing this in general, e.g. there could be two polylines in
                // different layers that snap to the same loop but start at different
                // positions.  In general we only consider input edge ids to be a hint
                // towards the preferred output ordering.)
                for (var e = 0; e < g_.NumEdges; ++e)
                {
                    if (used_[e]) continue;
                    var edge = g_.GetEdge(e);
                    if (edge.Item1 == edge.Item2)
                    {
                        // Note that it is safe to output degenerate edges as we go along,
                        // because this vertex has at least one non-degenerate outgoing edge and
                        // therefore we will (or just did) start an edge chain here.
                        OutputEdge(e);
                    }
                    else
                    {
                        SimplifyChain(edge.Item1, edge.Item2);
                    }
                }

                // Finally, copy the output edges into the appropriate layers.  They don't
                // need to be sorted because the input edges were also unsorted.
                for (var e = 0; e < new_edges_.Count; ++e)
                {
                    var layer = new_edge_layers_[e];
                    layer_edges_[layer].Add(new_edges_[e]);
                    layer_input_edge_ids_[layer].Add(new_input_edge_ids_[e]);
                }
            }

            // A helper class for determining whether a vertex can be an interior vertex
            // of a simplified edge chain.  Such a vertex must be adjacent to exactly two
            // vertices (across all layers combined), and in each layer the number of
            // incoming edges from one vertex must equal the number of outgoing edges to
            // the other vertex (in both directions).  Furthermore the vertex cannot have
            // any degenerate edges in a given layer unless it has at least one
            // non-degenerate edge in that layer as well.  (Note that usually there will
            // not be any degenerate edges at all, since most layer types discard them.)
            //
            // The last condition is necessary to prevent the following: suppose that one
            // layer has a chain ABC and another layer has a degenerate edge BB (with no
            // other edges incident to B).  Then we can't simplify ABC to AC because there
            // would be no suitable replacement for the edge BB (since the input edge that
            // mapped to BB can't be replaced by any of the edges AA, AC, or CC without
            // moving further than snap_radius).
            public class InteriorVertexMatcher
            {
                // Checks whether "v0" can be an interior vertex of an edge chain.
                public InteriorVertexMatcher(int v0)
                {
                    v0_ = v0; v1_ = -1; v2_ = -1; n0_ = 0; n1_ = 0; n2_ = 0; excess_out_ = 0;
                    too_many_endpoints_ = false;
                }

                // Starts analyzing the edges of a new layer.
                public void StartLayer()
                {
                    excess_out_ = n0_ = n1_ = n2_ = 0;
                }

                // This method should be called for each edge incident to "v0" in a given
                // layer.  (For degenerate edges, it should be called twice.)
                public void Tally(int v, bool outgoing)
                {
                    excess_out_ += outgoing ? 1 : -1;  // outdegree - indegree
                    if (v == v0_)
                    {
                        ++n0_;  // Counts both endpoints of each degenerate edge.
                    }
                    else
                    {
                        // We keep track of the total number of edges (incoming or outgoing)
                        // connecting v0 to up to two adjacent vertices.
                        if (v1_ < 0) v1_ = v;
                        if (v1_ == v)
                        {
                            ++n1_;
                        }
                        else
                        {
                            if (v2_ < 0) v2_ = v;
                            if (v2_ == v)
                            {
                                ++n2_;
                            }
                            else
                            {
                                too_many_endpoints_ = true;
                            }
                        }
                    }
                }

                // This method should be called after processing the edges for each layer.
                // It returns true if "v0" is an interior vertex based on the edges so far.
                public bool Matches()
                {
                    // We check that there are the same number of incoming and outgoing edges
                    // in each direction by verifying that (1) indegree(v0) == outdegree(v0)
                    // and (2) the total number of edges (incoming and outgoing) to "v1" and
                    // "v2" are equal.  We also check the condition on degenerate edges that
                    // is documented above.
                    return (!too_many_endpoints_ && excess_out_ == 0 && n1_ == n2_ &&
                            (n0_ == 0 || n1_ > 0));
                }

                private readonly int v0_;
                private int v1_;
                private int v2_;
                private int n0_, n1_, n2_;
                private int excess_out_;           // outdegree(v0) - indegree(v0)
                private bool too_many_endpoints_;  // Have we seen more than two adjacent vertices?
            }

            // Copies the given edge to the output and marks it as used.
            private void OutputEdge(int e)
            {
                new_edges_.Add(g_.GetEdge(e));
                new_input_edge_ids_.Add(g_.InputEdgeIdSetId(e));
                new_edge_layers_.Add(edge_layers_[e]);
                used_[e] = true;
            }
            // Returns the layer than a given input edge belongs to.
            private int InputEdgeLayer(int id)
            {
                // NOTE(ericv): If this method shows up in profiling, the result could be
                // stored with each edge (i.e., edge_layers_ and new_edge_layers_).
                Assert.True(id >= 0);
                return layer_begins_.GetUpperBound(id) - 1;
            }
            // Returns true if VertexId "v" can be an interior vertex of a simplified edge
            // chain.  (See the InteriorVertexMatcher class for what this implies.)
            private bool IsInterior(int v)
            {
                // Check a few simple prerequisites.
                if (out_.Degree(v) == 0) return false;
                if (out_.Degree(v) != in_.Degree(v)) return false;
                if (v < builder_.num_forced_sites_) return false;  // Keep forced vertices.

                // Sort the edges so that they are grouped by layer.
                var edges = new List<int>();
                foreach (var e in out_.EdgeIds(v)) edges.Add(e);
                foreach (var e in in_.EdgeIds(v)) edges.Add(e);
                edges.Sort(new GraphEdgeComp(edge_layers_));
                // Now feed the edges in each layer to the InteriorVertexMatcher.
                var matcher = new InteriorVertexMatcher(v);
                for (var e = 0; e < edges.Count; ) {
                    var layer = edge_layers_[e];
                    matcher.StartLayer();
                    for (; e < edges.Count && edge_layers_[e] == layer; ++e)
                    {
                        var edge = g_.GetEdge(edges[e]);
                        if (edge.Item1 == v) matcher.Tally(edge.Item2, true /*outgoing*/);
                        if (edge.Item2 == v) matcher.Tally(edge.Item1, false /*outgoing*/);
                    }
                    if (!matcher.Matches()) return false;
                }
                return true;
            }
            private class GraphEdgeComp : IComparer<int>
            {
                private readonly List<int> edge_layers_;
                public GraphEdgeComp(List<int> edge_layers) { edge_layers_ = edge_layers; }
                public int Compare(int x, int y)
                {
                    return edge_layers_[x].CompareTo(edge_layers_[y]);
                }
            }
            // Follows the edge chain starting with (v0, v1) until either we find a
            // non-interior vertex or we return to the original vertex v0.  At each vertex
            // we simplify a subchain of edges that is as long as possible.
            private void SimplifyChain(int v0, int v1)
            {
                // Avoid allocating "chain" each time by reusing it.
                var chain = tmp_vertices_;
                var vstart = v0;
                bool done;
                do
                {
                    // Simplify a subchain of edges starting (v0, v1).
                    var simplifier = new S2PolylineSimplifier(g_.Vertex(v0));
                    AvoidSites(v0, v0, v1, simplifier);
                    chain.Add(v0);
                    do
                    {
                        chain.Add(v1);
                        done = !is_interior_[v1] || v1 == vstart;
                        if (done) break;

                        // Attempt to extend the chain to the next vertex.
                        var vprev = v0;
                        v0 = v1;
                        v1 = FollowChain(vprev, v0);
                    } while (TargetInputVertices(v0, simplifier) &&
                             AvoidSites(chain[0], v0, v1, simplifier) &&
                             simplifier.Extend(g_.Vertex(v1)));

                    if (chain.Count == 2)
                    {
                        OutputAllEdges(chain[0], chain[1]);  // Could not simplify.
                    }
                    else
                    {
                        MergeChain(chain);
                    }
                    // Note that any degenerate edges that were not merged into a chain are
                    // output by EdgeChainSimplifier.Run().
                    chain.Clear();
                } while (!done);
            }
            // Given an edge (v0, v1) where v1 is an interior vertex, returns the (unique)
            // next vertex in the edge chain.
            private int FollowChain(int v0, int v1)
            {
                Assert.True(is_interior_[v1]);
                foreach (var e in out_.EdgeIds(v1))
                {
                    var v = g_.GetEdge(e).Item2;
                    if (v != v0 && v != v1) return v;
                }
                throw new ApplicationException("Could not find next edge in edge chain");
            }
            // Copies all input edges between v0 and v1 (in both directions) to the output.
            private void OutputAllEdges(int v0, int v1)
            {
                foreach (var e in out_.EdgeIds(v0, v1)) OutputEdge(e);
                foreach (var e in out_.EdgeIds(v1, v0)) OutputEdge(e);
            }
            // Ensures that the simplified edge passes within "edge_snap_radius" of all
            // the *input* vertices that snapped to the given vertex "v".
            private bool TargetInputVertices(int v, S2PolylineSimplifier simplifier)
            {
                foreach (var i in site_vertices_[v])
                {
                    if (!simplifier.TargetDisc(builder_.input_vertices_[i],
                                                builder_.edge_snap_radius_ca_))
                    {
                        return false;
                    }
                }
                return true;
            }
            // Given the starting vertex v0 and last edge (v1, v2) of an edge chain,
            // restricts the allowable range of angles in order to ensure that all sites
            // near the edge (v1, v2) are avoided by at least min_edge_vertex_separation.
            private bool AvoidSites(int v0, int v1, int v2, S2PolylineSimplifier simplifier)
            {
                var p0 = g_.Vertex(v0);
                var p1 = g_.Vertex(v1);
                var p2 = g_.Vertex(v2);
                var r1 = new S1ChordAngle(p0, p1);
                var r2 = new S1ChordAngle(p0, p2);

                // The distance from the start of the edge chain must increase monotonically
                // for each vertex, since we don't want to simplify chains that backtrack on
                // themselves (we want a parametric approximation, not a geometric one).
                if (r2 < r1) return false;

                // We also limit the maximum edge length in order to guarantee that the
                // simplified edge stays with max_edge_deviation() of all the input edges
                // that snap to it.
                if (r2 >= builder_.min_edge_length_to_split_ca_) return false;

                // Otherwise it is sufficient to consider the nearby sites (edge_sites_) for
                // a single input edge that snapped to (v1, v2) or (v2, v1).  This is
                // because each edge has a list of all sites within (max_edge_deviation +
                // min_edge_vertex_separation), and since the output edge is within
                // max_edge_deviation of all input edges, this list includes all sites
                // within min_edge_vertex_separation of the output edge.
                //
                // Usually there is only one edge to choose from, but it's not much more
                // effort to choose the edge with the shortest list of edge_sites_.
                var best = -1;
                var edge_sites = builder_.edge_sites_;
                foreach (var e in out_.EdgeIds(v1, v2))
                {
                    foreach (var id in g_.InputEdgeIds(e))
                    {
                        if (best < 0 || edge_sites[id].Count < edge_sites[best].Count)
                            best = id;
                    }
                }
                foreach (var e in out_.EdgeIds(v2, v1))
                {
                    foreach (var id in g_.InputEdgeIds(e))
                    {
                        if (best < 0 || edge_sites[id].Count < edge_sites[best].Count)
                            best = id;
                    }
                }
                Assert.True(best >= 0);  // Because there is at least one outgoing edge.

                foreach (var v in edge_sites[best])
                {
                    // This test is optional since these sites are excluded below anyway.
                    if (v == v0 || v == v1 || v == v2) continue;

                    // We are only interested in sites whose distance from "p0" is in the
                    // range (r1, r2).  Sites closer than "r1" have already been processed,
                    // and sites further than "r2" aren't relevant yet.
                    var p = g_.Vertex(v);
                    var r = new S1ChordAngle(p0, p);
                    if (r <= r1 || r >= r2) continue;

                    // We need to figure out whether this site is to the left or right of the
                    // edge chain.  For the first edge this is easy.  Otherwise, since we are
                    // only considering sites in the radius range (r1, r2), we can do this by
                    // checking whether the site is to the left of the wedge (p0, p1, p2).
                    var disc_on_left = (v1 == v0) ? (S2Pred.Sign(p1, p2, p) > 0)
                                                   : S2Pred.OrderedCCW(p0, p2, p, p1);
                    if (!simplifier.AvoidDisc(p, builder_.min_edge_site_separation_ca_,
                                               disc_on_left))
                    {
                        return false;
                    }
                }
                return true;
            }
            // Given the vertices in a simplified edge chain, adds the corresponding
            // simplified edge(s) to the output.  Note that (1) the edge chain may exist
            // in multiple layers, (2) the edge chain may exist in both directions, (3)
            // there may be more than one copy of an edge chain (in either direction)
            // within a single layer.
            private void MergeChain(List<int> vertices)
            {
                // Suppose that all interior vertices have M outgoing edges and N incoming
                // edges.  Our goal is to group the edges into M outgoing chains and N
                // incoming chains, and then replace each chain by a single edge.
                var merged_input_ids = new List<List<int>>();
                var degenerate_ids = new List<int>();
                int num_out;  // Edge count in the outgoing direction.
                for (var i = 1; i < vertices.Count; ++i)
                {
                    var v0 = vertices[i - 1];
                    var v1 = vertices[i];
                    var out_edges = out_.EdgeIds(v0, v1);
                    var in_edges = out_.EdgeIds(v1, v0);
                    if (i == 1)
                    {
                        // Allocate space to store the input edge ids associated with each edge.
                        num_out = out_edges.Count;
                        merged_input_ids.Capacity = num_out + in_edges.Count;
                        foreach (var ids in merged_input_ids)
                        {
                            ids.Capacity = vertices.Count - 1;
                        }
                    }
                    else
                    {
                        // For each interior vertex, we build a list of input edge ids
                        // associated with degenerate edges.  Each input edge ids will be
                        // assigned to one of the output edges later.  (Normally there are no
                        // degenerate edges at all since most layer types don't want them.)
                        Assert.True(is_interior_[v0]);
                        foreach (var e in out_.EdgeIds(v0, v0))
                        {
                            foreach (var id in g_.InputEdgeIds(e))
                            {
                                degenerate_ids.Add(id);
                            }
                            used_[e] = true;
                        }
                    }
                    // Because the edges were created in layer order, and all sorts used are
                    // stable, the edges are still in layer order.  Therefore we can simply
                    // merge together all the edges in the same relative position.
                    var j = 0;
                    foreach (var e in out_edges)
                    {
                        foreach (var id in g_.InputEdgeIds(e))
                        {
                            merged_input_ids[j].Add(id);
                        }
                        used_[e] = true;
                        ++j;
                    }
                    foreach (var e in in_edges)
                    {
                        foreach (var id in g_.InputEdgeIds(e))
                        {
                            merged_input_ids[j].Add(id);
                        }
                        used_[e] = true;
                        ++j;
                    }
                    Assert.True(merged_input_ids.Count == j);
                }
                if (degenerate_ids.Any())
                {
                    degenerate_ids.Sort();
                    AssignDegenerateEdges(degenerate_ids, merged_input_ids);
                }
                // Output the merged edges.
                int v0_ = vertices[0], v1_ = vertices[1], vb = vertices.Last();
                foreach (var e in out_.EdgeIds(v0_, v1_))
                {
                    new_edges_.Add(new(v0_, vb));
                    new_edge_layers_.Add(edge_layers_[e]);
                }
                foreach (var e in out_.EdgeIds(v1_, v0_))
                {
                    new_edges_.Add(new(vb, v0_));
                    new_edge_layers_.Add(edge_layers_[e]);
                }
                foreach (var ids in merged_input_ids)
                {
                    new_input_edge_ids_.Add(input_edge_id_set_lexicon_.Add(ids));
                }
            }
            // Given a list of the input edge ids associated with degenerate edges in the
            // interior of an edge chain, assigns each input edge id to one of the output
            // edges.
            private void AssignDegenerateEdges(List<int> degenerate_ids, List<List<int>> merged_input_ids)
            {
                // Each degenerate edge is assigned to an output edge in the appropriate
                // layer.  If there is more than one candidate, we use heuristics so that if
                // the input consists of a chain of edges provided in consecutive order
                // (some of which became degenerate), then all those input edges are
                // assigned to the same output edge.  For example, suppose that one output
                // edge is labeled with input edges 3,4,7,8, while another output edge is
                // labeled with input edges 50,51,54,57.  Then if we encounter degenerate
                // edges labeled with input edges 5 and 6, we would prefer to assign them to
                // the first edge (yielding the continuous range 3,4,5,6,7,8).
                //
                // The heuristic below is only smart enough to handle the case where the
                // candidate output edges have non-overlapping ranges of input edges.
                // (Otherwise there is probably not a good heuristic for assigning the
                // degenerate edges in any case.)

                // Duplicate edge ids don't affect the heuristic below, so we don't bother
                // removing them.  (They will be removed by IdSetLexicon.Add.)
                foreach (var ids in merged_input_ids) ids.Sort();

                // Sort the output edges by their minimum input edge id.  This is sufficient
                // for the purpose of determining which layer they belong to.  With
                // EdgeType.UNDIRECTED, some edges might not have any input edge ids (i.e.,
                // if they consist entirely of siblings of input edges).  We simply remove
                // such edges from the lists of candidates.
                var order = new List<int>(merged_input_ids.Count);
                for (var i = 0; i < merged_input_ids.Count; ++i)
                {
                    if (merged_input_ids[i].Any()) order.Add(i);
                }
                order.Sort((int i, int j) => merged_input_ids[i][0].CompareTo(merged_input_ids[j][0]));

                // Now determine where each degenerate edge should be assigned.
                var comp = new MergedIdsComp(merged_input_ids);
                foreach (var degenerate_id in degenerate_ids)
                {
                    var layer = InputEdgeLayer(degenerate_id);

                    // Find the first output edge whose range of input edge ids starts after
                    // "degenerate_id".  If the previous edge belongs to the correct layer,
                    // then we assign the degenerate edge to it.
                    var index = order.GetUpperBound(degenerate_id, comp);
                    if (index < order.Count)
                    {
                        if (merged_input_ids[(int)order[index - 1]][0] >= layer_begins_[layer]) --index;
                    }
                    Assert.True(layer == InputEdgeLayer(merged_input_ids[(int)order[index]][0]));
                    merged_input_ids[(int)order[index]].Add(degenerate_id);
                }
            }

            private class MergedIdsComp : IComparer<int>
            {
                private readonly List<List<int>> merged_input_ids_;
                public MergedIdsComp(List<List<int>> merged_input_ids) { merged_input_ids_ = merged_input_ids; }

                public int Compare(int x, int y)
                {
                    return x.CompareTo(merged_input_ids_[y][0]);
                }
            }

            private readonly S2Builder builder_;
            private readonly Graph g_;
            private readonly Graph.VertexInMap in_;
            private readonly Graph.VertexOutMap out_;
            private readonly List<int> edge_layers_;
            private readonly List<List<int>> site_vertices_;
            private readonly List<List<OutputEdge>> layer_edges_;
            private readonly List<List<int>> layer_input_edge_ids_;
            private readonly IdSetLexicon input_edge_id_set_lexicon_;

            // Convenience member copied from builder_.
            private readonly List<int> layer_begins_;

            // is_interior_[v] indicates that VertexId "v" is eligible to be an interior
            // vertex of a simplified edge chain.  You can think of it as vertex whose
            // indegree and outdegree are both 1 (although the actual definition is a
            // bit more complicated because of duplicate edges and layers).
            private readonly List<bool> is_interior_;

            // used_[e] indicates that EdgeId "e" has already been processed.
            private readonly List<bool> used_;

            // Temporary vectors, declared here to avoid repeated allocation.
            private readonly List<int> tmp_vertices_ = new List<int>();

            // The output edges after simplification.
            private readonly List<OutputEdge> new_edges_ = new();
            private readonly List<int> new_input_edge_ids_ = new List<int>();
            private readonly List<int> new_edge_layers_ = new List<int>();
        }

        // This class is only needed by S2Builder.Layer implementations.  A layer is
        // responsible for assembling an S2Builder.Graph of snapped edges into the
        // desired output format (e.g., an S2Polygon).  The GraphOptions class allows
        // each Layer type to specify requirements on its input graph: for example, if
        // DegenerateEdges.DISCARD is specified, then S2Builder will ensure that all
        // degenerate edges are removed before passing the graph to the S2Layer.Build
        // method.
        public class GraphOptions : IEquatable<GraphOptions>
        {
            #region Fields, Constants

            // Specifies whether the S2Builder input edges should be treated as
            // undirected.  If true, then all input edges are duplicated into pairs
            // consisting of an edge and a sibling (reverse) edge.  The layer
            // implementation is responsible for ensuring that exactly one edge from
            // each pair is used in the output, i.e. *only half* of the graph edges will
            // be used.  (Note that some values of the sibling_pairs() option
            // automatically take care of this issue by removing half of the edges and
            // changing edge_type() to DIRECTED.)
            //
            // DEFAULT: EdgeType.DIRECTED
            public EdgeType EdgeType_ { get; set; }

            public DegenerateEdges DegenerateEdges_ { get; set; }

            public DuplicateEdges DuplicateEdges_ { get; set; }

            public SiblingPairs SiblingPairs_ { get; set; }

            // This is a specialized option that is only needed by clients want to work
            // with the graphs for multiple layers at the same time (e.g., in order to
            // check whether the same edge is present in two different graphs).  [Note
            // that if you need to do this, usually it is easier just to build a single
            // graph with suitable edge labels.]
            //
            // When there are a large number of layers, then by default S2Builder builds
            // a minimal subgraph for each layer containing only the vertices needed by
            // the edges in that layer.  This ensures that layer types that iterate over
            // the vertices run in time proportional to the size of that layer rather
            // than the size of all layers combined.  (For example, if there are a
            // million layers with one edge each, then each layer would be passed a
            // graph with 2 vertices rather than 2 million vertices.)
            //
            // If this option is set to false, this optimization is disabled.  Instead
            // the graph passed to this layer will contain the full set of vertices.
            // (This is not recommended when the number of layers could be large.)
            //
            // DEFAULT: true
            public bool AllowVertexFiltering { get; set; }

            #endregion

            #region Constructors
            
            // All S2Builder.Layer subtypes should specify GraphOptions explicitly
            // using this constructor, rather than relying on default values.
            public GraphOptions(EdgeType edge_type, DegenerateEdges degenerate_edges,
                           DuplicateEdges duplicate_edges, SiblingPairs sibling_pairs)
            {
                EdgeType_ = edge_type; DegenerateEdges_ = degenerate_edges;
                DuplicateEdges_ = duplicate_edges; SiblingPairs_ = sibling_pairs;
                AllowVertexFiltering = true;
            }

            // The default options specify that all edges should be kept, since this
            // produces the least surprising output and makes it easier to diagnose the
            // problem when an option is left unspecified.
            public GraphOptions()
            {
                EdgeType_ = EdgeType.DIRECTED;
                DegenerateEdges_ = DegenerateEdges.KEEP;
                DuplicateEdges_ = DuplicateEdges.KEEP;
                SiblingPairs_ = SiblingPairs.KEEP;
                AllowVertexFiltering = true;
            } 
            
            #endregion

            #region IEquatable

            public static bool operator ==(GraphOptions x, GraphOptions y) => Equals(x, y);

            public static bool operator !=(GraphOptions x, GraphOptions y) => !Equals(x, y);
            public override bool Equals(object other) => other is GraphOptions go && Equals(go);
            public bool Equals(GraphOptions other)
            {
                return EdgeType_ == other.EdgeType_ &&
                    DegenerateEdges_ == other.DegenerateEdges_ &&
                    DuplicateEdges_ == other.DuplicateEdges_ &&
                    SiblingPairs_ == other.SiblingPairs_ &&
                    AllowVertexFiltering == other.AllowVertexFiltering;
            }
            public override int GetHashCode()
            {
                return HashCode.Combine(EdgeType_, DegenerateEdges_, DuplicateEdges_, SiblingPairs_, AllowVertexFiltering);
            }

            #endregion

            // Controls how degenerate edges (i.e., an edge from a vertex to itself) are
            // handled.  Such edges may be present in the input, or they may be created
            // when both endpoints of an edge are snapped to the same output vertex.
            // The options available are:
            //
            // DISCARD: Discards all degenerate edges.  This is useful for layers that
            //          do not support degeneracies, such as S2PolygonLayer.
            //
            // DISCARD_EXCESS: Discards all degenerate edges that are connected to
            //                 non-degenerate edges.  (Any remaining duplicate edges can
            //                 be merged using DuplicateEdges.MERGE.)  This is useful
            //                 for simplifying polygons while ensuring that loops that
            //                 collapse to a single point do not disappear.
            //
            // KEEP: Keeps all degenerate edges.  Be aware that this may create many
            //       redundant edges when simplifying geometry (e.g., a polyline of the
            //       form AABBBBBCCCCCCDDDD).  DegenerateEdges.KEEP is mainly useful
            //       for algorithms that require an output edge for every input edge.
            //
            // DEFAULT: DegenerateEdges.KEEP
            public enum DegenerateEdges { DISCARD, DISCARD_EXCESS, KEEP }

            // Controls how duplicate edges (i.e., edges that are present multiple
            // times) are handled.  Such edges may be present in the input, or they can
            // be created when vertices are snapped together.  When several edges are
            // merged, the result is a single edge labelled with all of the original
            // input edge ids.
            //
            // DEFAULT: DuplicateEdges.KEEP
            public enum DuplicateEdges { MERGE, KEEP }

            // Controls how sibling edge pairs (i.e., pairs consisting of an edge and
            // its reverse edge) are handled.  Layer types that define an interior
            // (e.g., polygons) normally discard such edge pairs since they do not
            // affect the result (i.e., they define a "loop" with no interior).  The
            // various options include:
            //
            // DISCARD: Discards all sibling edge pairs.
            //
            // DISCARD_EXCESS: Like DISCARD, except that a single sibling pair is kept
            //                 if the result would otherwise be empty.  This is useful
            //                 for polygons with degeneracies (S2LaxPolygonShape), and
            //                 for simplifying polylines while ensuring that they are
            //                 not split into multiple disconnected pieces.
            //
            // KEEP: Keeps sibling pairs.  This can be used to create polylines that
            //       double back on themselves, or degenerate loops (with a layer type
            //       such as S2LaxPolygonShape).
            //
            // REQUIRE: Requires that all edges have a sibling (and returns an error
            //          otherwise).  This is useful with layer types that create a
            //          collection of adjacent polygons (a polygon mesh).
            //
            // CREATE: Ensures that all edges have a sibling edge by creating them if
            //         necessary.  This is useful with polygon meshes where the input
            //         polygons do not cover the entire sphere.  Such edges always
            //         have an empty set of labels.
            //
            // If edge_type() is EdgeType.UNDIRECTED, a sibling edge pair is considered
            // to consist of four edges (two duplicate edges and their siblings), since
            // only two of these four edges will be used in the final output.
            //
            // Furthermore, since the options REQUIRE and CREATE guarantee that all
            // edges will have siblings, S2Builder implements these options for
            // undirected edges by discarding half of the edges in each direction and
            // changing the edge_type() to EdgeType.DIRECTED.  For example, two
            // undirected input edges between vertices A and B would first be converted
            // into two directed edges in each direction, and then one edge of each pair
            // would be discarded leaving only one edge in each direction.
            //
            // Degenerate edges are considered not to have siblings.  If such edges are
            // present, they are passed through unchanged by SiblingPairs.DISCARD.  For
            // SiblingPairs.REQUIRE or SiblingPairs.CREATE with undirected edges, the
            // number of copies of each degenerate edge is reduced by a factor of two.
            //
            // Any of the options that discard edges (DISCARD, DISCARD_EXCESS, and
            // REQUIRE/CREATE in the case of undirected edges) have the side effect that
            // when duplicate edges are present, all of the corresponding edge labels
            // are merged together and assigned to the remaining edges.  (This avoids
            // the problem of having to decide which edges are discarded.)  Note that
            // this merging takes place even when all copies of an edge are kept, and
            // that even labels attached to duplicate degenerate edges are merged.  For
            // example, consider the graph {AB1, AB2, BA3, CD4, CD5} (where XYn denotes
            // an edge from X to Y with label "n").  With SiblingPairs.DISCARD, we need
            // to discard one of the copies of AB.  But which one?  Rather than choosing
            // arbitrarily, instead we merge the labels of all duplicate edges (even
            // ones where no sibling pairs were discarded), yielding {AB12, CD45, CD45}
            // (assuming that duplicate edges are being kept).
            //
            // DEFAULT: SiblingPairs.KEEP
            public enum SiblingPairs { DISCARD, DISCARD_EXCESS, KEEP, REQUIRE, CREATE } 
        }

        // For output layers that represent polygons, there is an ambiguity inherent
        // in spherical geometry that does not exist in planar geometry.  Namely, if
        // a polygon has no edges, does it represent the empty polygon (containing
        // no points) or the full polygon (containing all points)?  This ambiguity
        // also occurs for polygons that consist only of degeneracies, e.g. a
        // degenerate loop with only two edges could be either a degenerate shell in
        // the empty polygon or a degenerate hole in the full polygon.
        //
        // To resolve this ambiguity, an IsFullPolygonPredicate may be specified for
        // each output layer (see AddIsFullPolygonPredicate below).  If the output
        // after snapping consists only of degenerate edges and/or sibling pairs
        // (including the case where there are no edges at all), then the layer
        // implementation calls the given predicate to determine whether the polygon
        // is empty or full except for those degeneracies.  The predicate is given
        // an S2Builder.Graph containing the output edges, but note that in general
        // the predicate must also have knowledge of the input geometry in order to
        // determine the correct result.
        //
        // This predicate is only needed by layers that are assembled into polygons.
        // It is not used by other layer types.
        public delegate bool IsFullPolygonPredicate(Graph g, out S2Error error);

        // Default constructor; requires Init() to be called.
        public S2Builder() { }

        // Convenience constructor that calls Init().  Note that to use the default
        // options, C++ syntax requires an extra layer of parentheses:
        //
        //   S2Builder builder{S2Builder.Options()};
        public S2Builder(Options options)
        {
            Init(options);
        }

        // Initializes an S2Builder with the given options.
        public void Init(Options options)
        {
            Options_ = options;
            var snap_function = options.SnapFunction;
            var snap_radius = snap_function.SnapRadius;
            Assert.True(snap_radius <= SnapFunction.kMaxSnapRadius);

            // Convert the snap radius to an S1ChordAngle.  This is the "true snap
            // radius" used when evaluating exact predicates (s2predicates.h).
            site_snap_radius_ca_ = new S1ChordAngle(snap_radius);

            // When split_crossing_edges() is true, we need to use a larger snap radius
            // for edges than for vertices to ensure that both edges are snapped to the
            // edge intersection location.  This is because the computed intersection
            // point is not exact; it may be up to kIntersectionError away from its true
            // position.  The computed intersection point might then be snapped to some
            // other vertex up to snap_radius away.  So to ensure that both edges are
            // snapped to a common vertex, we need to increase the snap radius for edges
            // to at least the sum of these two values (calculated conservatively).
            var edge_snap_radius = snap_radius;
            if (!options.SplitCrossingEdges)
            {
                edge_snap_radius_ca_ = site_snap_radius_ca_;
            }
            else
            {
                edge_snap_radius += S2EdgeCrossings.kIntersectionErrorS1Angle;
                edge_snap_radius_ca_ = RoundUp(edge_snap_radius);
            }
            snapping_requested_ = (edge_snap_radius > S1Angle.Zero);

            // Compute the maximum distance that a vertex can be separated from an
            // edge while still affecting how that edge is snapped.
            max_edge_deviation_ = snap_function.MaxEdgeDeviation();
            edge_site_query_radius_ca_ = new S1ChordAngle(
                max_edge_deviation_ + snap_function.MinEdgeVertexSeparation());

            // Compute the maximum edge length such that even if both endpoints move by
            // the maximum distance allowed (i.e., snap_radius), the center of the edge
            // will still move by less than max_edge_deviation().  This saves us a lot
            // of work since then we don't need to check the actual deviation.
            min_edge_length_to_split_ca_ = S1ChordAngle.FromRadians(
                2 * Math.Acos(snap_radius.Sin / max_edge_deviation_.Sin));

            // If the condition below is violated, then AddExtraSites() needs to be
            // modified to check that snapped edges pass on the same side of each "site
            // to avoid" as the input edge.  Currently it doesn't need to do this
            // because the condition below guarantees that if the snapped edge passes on
            // the wrong side of the site then it is also too close, which will cause a
            // separation site to be added.
            //
            // Currently max_edge_deviation() is at most 1.1 * snap_radius(), whereas
            // min_edge_vertex_separation() is at least 0.219 * snap_radius() (based on
            // S2CellIdSnapFunction, which is currently the worst case).
            Assert.True(snap_function.MaxEdgeDeviation() <=
                      snap_function.SnapRadius+
                      snap_function.MinEdgeVertexSeparation());

            // To implement idempotency, we check whether the input geometry could
            // possibly be the output of a previous S2Builder invocation.  This involves
            // testing whether any site/site or edge/site pairs are too close together.
            // This is done using exact predicates, which require converting the minimum
            // separation values to an S1ChordAngle.
            min_site_separation_ = snap_function.MinVertexSeparation();
            min_site_separation_ca_ = new S1ChordAngle(min_site_separation_);
            min_edge_site_separation_ca_ =
                new S1ChordAngle(snap_function.MinEdgeVertexSeparation());

            // This is an upper bound on the distance computed by S2ClosestPointQuery
            // where the true distance might be less than min_edge_site_separation_ca_.
            min_edge_site_separation_ca_limit_ =
                AddPointToEdgeError(min_edge_site_separation_ca_);

            // Compute the maximum possible distance between two sites whose Voronoi
            // regions touch.  (The maximum radius of each Voronoi region is
            // edge_snap_radius_.)  Then increase this bound to account for errors.
            max_adjacent_site_separation_ca_ =
                AddPointToPointError(RoundUp(2 * edge_snap_radius));

            // Finally, we also precompute sin^2(edge_snap_radius), which is simply the
            // squared distance between a vertex and an edge measured perpendicular to
            // the plane containing the edge, and increase this value by the maximum
            // error in the calculation to compare this distance against the bound.
            var d = edge_snap_radius.Sin;
            edge_snap_radius_sin2_ = d * d;
            edge_snap_radius_sin2_ += ((9.5 * d + 2.5 + 2 * Math.Sqrt(3)) * d + 9 * S2Constants.DoubleEpsilon) * S2Constants.DoubleEpsilon;

            // Initialize the current label set.
            label_set_id_ = IdSetLexicon.kEmptySetId;
            label_set_modified_ = false;

            // If snapping was requested, we try to determine whether the input geometry
            // already meets the output requirements.  This is necessary for
            // idempotency, and can also save work.  If we discover any reason that the
            // input geometry needs to be modified, snapping_needed_ is set to true.
            snapping_needed_ = false;
        }

        // Starts a new output layer.  This method must be called before adding any
        // edges to the S2Builder.  You may call this method multiple times to build
        // multiple geometric objects that are snapped to the same set of sites.
        //
        // For example, if you have a set of contour lines, then you could put each
        // contour line in a separate layer.  This keeps the contour lines separate
        // from each other, while also ensuring that no crossing edges are created
        // when they are snapped and/or simplified.  (This is not true if the
        // contour lines are snapped or simplified independently.)
        //
        // Similarly, if you have a set of polygons that share common boundaries
        // (e.g., countries), you can snap and/or simplify them at the same time by
        // putting them in different layers, while ensuring that their boundaries
        // remain consistent (i.e., no crossing edges or T-vertices are introduced).
        //
        // Ownership of the layer is transferred to the S2Builder.  Example usage:
        //
        // S2Polyline line1, line2;
        // builder.StartLayer(new s2builderutil.S2PolylineLayer(line1)));
        // ... Add edges using builder.AddEdge(), etc ...
        // builder.StartLayer(new s2builderutil.S2PolylineLayer(line2)));
        // ... Add edges using builder.AddEdge(), etc ...
        // S2Error error;
        // Assert.True(builder.Build(out error)) << error;  // Builds "line1" & "line2"
        public void StartLayer(S2Builder.Layer layer)
        {
            layer_options_.Add(layer.GraphOptions_());
            layer_begins_.Add(input_edges_.Count);
            layer_is_full_polygon_predicates_.Add(IsFullPolygon(false));
            layers_.Add(layer);
        }

        // Adds a degenerate edge (representing a point) to the current layer.
        public void AddPoint(S2Point v)
        {
            AddEdge(v, v);
        }

        // Adds the given edge to the current layer.
        public void AddEdge(S2Point v0, S2Point v1)
        {
            Assert.True(layers_.Any()); // Call StartLayer before adding any edges

            if (v0 == v1 && (layer_options_.Last().DegenerateEdges_==
                             GraphOptions.DegenerateEdges.DISCARD))
            {
                return;
            }
            var j0 = AddVertex(v0);
            var j1 = AddVertex(v1);
            input_edges_.Add(new(j0, j1));

            // If there are any labels, then attach them to this input edge.
            if (label_set_modified_)
            {
                if (!label_set_ids_.Any())
                {
                    // Populate the missing entries with empty label sets.
                    label_set_ids_.Fill(label_set_id_, input_edges_.Count - 1);
                }
                label_set_id_ = label_set_lexicon_.Add(label_set_);
                label_set_ids_.Add(label_set_id_);
                label_set_modified_ = false;
            }
            else if (label_set_ids_.Any())
            {
                label_set_ids_.Add(label_set_id_);
            }
        }

        // Adds the edges in the given polyline.  (Note that if the polyline
        // consists of 0 or 1 vertices, this method does nothing.)
        public void AddPolyline(S2Polyline polyline)
        {
            var n = polyline.NumVertices;
            for (var i = 1; i < n; ++i)
            {
                AddEdge(polyline.Vertex(i - 1), polyline.Vertex(i));
            }
        }

        // Adds the edges in the given loop.  If the sign() of the loop is negative
        // (i.e. this loop represents a hole within a polygon), the edge directions
        // are automatically reversed to ensure that the polygon interior is always
        // to the left of every edge.
        public void AddLoop(S2Loop loop)
        {
            // Ignore loops that do not have a boundary.
            if (loop.IsEmptyOrFull) return;

            // For loops that represent holes, we add the edge from vertex n-1 to vertex
            // n-2 first.  This is because these edges will be assembled into a
            // clockwise loop, which will eventually be normalized in S2Polygon by
            // calling S2Loop.Invert().  S2Loop.Invert() reverses the order of the
            // vertices, so to end up with the original vertex order (0, 1, ..., n-1) we
            // need to build a clockwise loop with vertex order (n-1, n-2, ..., 0).
            // This is done by adding the edge (n-1, n-2) first, and then ensuring that
            // Build() assembles loops starting from edges in the order they were added.
            var n = loop.NumVertices;
            for (var i = 0; i < n; ++i)
            {
                AddEdge(loop.OrientedVertex(i), loop.OrientedVertex(i + 1));
            }
        }

        // Adds the loops in the given polygon.  Loops representing holes have their
        // edge directions automatically reversed as described for AddLoop().  Note
        // that this method does not distinguish between the empty and full polygons,
        // i.e. adding a full polygon has the same effect as adding an empty one.
        public void AddPolygon(S2Polygon polygon)
        {
            for (var i = 0; i < polygon.NumLoops(); ++i)
            {
                AddLoop(polygon.Loop(i));
            }
        }

        // Adds the edges of the given shape to the current layer.
        public void AddShape(S2Shape shape)
        {
            for (int e = 0, n = shape.NumEdges; e < n; ++e)
            {
                var edge = shape.GetEdge(e);
                AddEdge(edge.V0, edge.V1);
            }
        }

        // For layers that are assembled into polygons, this method specifies a
        // predicate that is called when the output consists entirely of degenerate
        // edges and/or sibling pairs.  The predicate is given an S2Builder.Graph
        // containing the output edges (if any) and is responsible for deciding
        // whether this graph represents the empty polygon (possibly with degenerate
        // shells) or the full polygon (possibly with degenerate holes).  Note that
        // this cannot be determined from the output edges alone; it also requires
        // knowledge of the input geometry.  (Also see IsFullPolygonPredicate above.)
        //
        // This method should be called at most once per layer; additional calls
        // simply overwrite the previous value for the current layer.
        //
        // The default predicate simply returns false (i.e., degenerate polygons are
        // assumed to be empty).  Arguably it would better to return an error in
        // this case, but the fact is that relatively few clients need to be able to
        // construct full polygons, and it is unreasonable to expect all such
        // clients to supply an appropriate predicate.
        //
        // The reason for having a predicate rather than a boolean value is that the
        // predicate is responsible for determining whether the output polygon is
        // empty or full.  In general the input geometry is not degenerate, but
        // rather collapses into a degenerate configuration due to snapping and/or
        // simplification.
        //
        // TODO(ericv): Provide standard predicates to handle common cases,
        // e.g. valid input geometry that becomes degenerate due to snapping.
        public void AddIsFullPolygonPredicate(IsFullPolygonPredicate predicate)
        {
            layer_is_full_polygon_predicates_[^1] = predicate;
        }

        // A predicate that returns an error indicating that no polygon predicate
        // has been specified.
        public static bool IsFullPolygonUnspecified(out S2Error error)
        {
            error = new(S2ErrorCode.BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED, "A degenerate polygon was found, but no predicate was specified to determine whether the polygon is empty or full.  Call S2Builder.AddIsFullPolygonPredicate() to fix this problem.");
            // Assumes the polygon is empty.
            return false;
        }

        // Returns a predicate that returns a constant value (true or false);
        public static IsFullPolygonPredicate IsFullPolygon(bool is_full)
        {
            return (Graph g, out S2Error error) => {
                error = S2Error.OK;
                return is_full;
            };
        }

        // Forces a vertex to be located at the given position.  This can be used to
        // prevent certain input vertices from moving.  However if you are trying to
        // preserve part of the input boundary, be aware that this option does not
        // prevent edges from being split by new vertices.
        //
        // Forced vertices are never snapped; if this is desired then you need to
        // call options().snap_function().SnapPoint() explicitly.  Forced vertices
        // are also never simplified away (if simplify_edge_chains() is used).
        //
        // Caveat: Since this method can place vertices arbitrarily close together,
        // S2Builder makes no minimum separation guaranteees with forced vertices.
        public void ForceVertex(S2Point vertex)
        {
            sites_.Add(vertex);
        }

        // Every edge can have a set of non-negative integer labels attached to it.
        // When used with an appropriate layer type, you can then retrieve the
        // labels associated with each output edge.  This can be useful when merging
        // or combining data from several sources.  (Note that in many cases it is
        // easier to use separate output layers rather than labels.)
        //
        // Labels are 32-bit non-negative integers.  To support other label types,
        // you can use ValueLexicon to store the set of unique labels seen so far:
        //
        //   ValueLexicon<MyLabel> my_label_lexicon;
        //   builder.set_label(my_label_lexicon.Add(label));
        //
        // The current set of labels is represented as a stack.  This makes it easy
        // to add and remove labels hierarchically (e.g., polygon 5, loop 2).  Use
        // set_label() and clear_labels() if you need at most one label per edge.
        //

        // Clear the stack of labels.
        public void ClearLabels()
        {
            label_set_.Clear();
            label_set_modified_ = true;
        }

        // Add a label to the stack.
        // REQUIRES: label >= 0.
        public void PushLabel(int label)
        {
            Assert.True(label >= 0);
            label_set_.Add(label);
            label_set_modified_ = true;
        }

        // Remove a label from the stack.
        public void PopLabel()
        {
            label_set_.RemoveAt(label_set_.Count - 1);
            label_set_modified_ = true;
        }

        // Convenience function that clears the stack and adds a single label.
        // REQUIRES: label >= 0.
        public void SetLabel(int label)
        {
            Assert.True(label >= 0);
            label_set_.Capacity = 1;
            label_set_[0] = label;
            label_set_modified_ = true;
        }

        // Performs the requested edge splitting, snapping, simplification, etc, and
        // then assembles the resulting edges into the requested output layers.
        //
        // Returns true if all edges were assembled; otherwise sets "error"
        // appropriately.  Depending on the error, some or all output layers may
        // have been created.  Automatically resets the S2Builder state so that it
        // can be reused.
        //
        // REQUIRES: error != null.
        public bool Build(out S2Error error)
        {
            // Mark the end of the last layer.
            layer_begins_.Add(input_edges_.Count);

            // See the algorithm overview at the top of this file.
            if (snapping_requested_ && !Options_.Idempotent)
            {
                snapping_needed_ = true;
            }
            ChooseSites();
            BuildLayers();
            Reset();
            error = error_;
            return error.IsOk;
        }

        // Clears all input data and resets the builder state.  Any options
        // specified are preserved.
        public void Reset()
        {
            input_vertices_.Clear();
            input_edges_.Clear();
            layers_.Clear();
            layer_options_.Clear();
            layer_begins_.Clear();
            layer_is_full_polygon_predicates_.Clear();
            label_set_ids_.Clear();
            label_set_lexicon_.Clear();
            label_set_.Clear();
            label_set_modified_ = false;
            sites_.Clear();
            edge_sites_.Clear();
            snapping_needed_ = false;
        }

        //////////////////////  Input Types  /////////////////////////
        // All types associated with the S2Builder inputs are prefixed with "Input".

        // Input vertices are stored in a vector, with some removal of duplicates.
        // Edges are represented as (VertexId, VertexId) pairs.  All edges are stored
        // in a single vector; each layer corresponds to a contiguous range.

        private int AddVertex(S2Point v)
        {
            // Remove duplicate vertices that follow the pattern AB, BC, CD.  If we want
            // to do anything more sophisticated, either use a ValueLexicon, or sort the
            // vertices once they have all been added, remove duplicates, and update the
            // edges.
            if (!input_vertices_.Any() || v != input_vertices_.Last())
            {
                input_vertices_.Add(v);
            }
            return input_vertices_.Count - 1;
        }
        private void ChooseSites()
        {
            if (!input_vertices_.Any()) return;

            MutableS2ShapeIndex input_edge_index = new();
            input_edge_index.Add(new VertexIdEdgeVectorShape(input_edges_, input_vertices_));
            if (Options_.SplitCrossingEdges)
            {
                AddEdgeCrossings(input_edge_index);
            }
            if (snapping_requested_)
            {
                var site_index = new S2PointIndex<int>();
                AddForcedSites(site_index);
                ChooseInitialSites(site_index);
                CollectSiteEdges(site_index);
            }
            if (snapping_needed_)
            {
                AddExtraSites(input_edge_index);
            }
            else
            {
                CopyInputEdges();
            }
        }
        private void CopyInputEdges()
        {
            // Sort the input vertices, discard duplicates, and update the input edges
            // to refer to the pruned vertex list.  (We sort in the same order used by
            // ChooseInitialSites() to avoid inconsistencies in tests.)
            var sorted = SortInputVertices();
            var vmap = new int[input_vertices_.Count];
            sites_.Clear();
            sites_.Capacity = input_vertices_.Count;
            for (var in_ = 0; in_ < sorted.Count; ) {
                var site = input_vertices_[sorted[in_].Item2];
                vmap[sorted[in_].Item2] = sites_.Count;
                while (++in_ < sorted.Count && input_vertices_[sorted[in_].Item2] == site) {
                    vmap[sorted[in_].Item2] = sites_.Count;
                }
                sites_.Add(site);
            }
            input_vertices_ = sites_;
            for (var i = 0; i < input_edges_.Count; i++)
            {
                var e = input_edges_[i];
                input_edges_[i] = new(vmap[e.Item1], vmap[e.Item2]);
            }
        }
        private List<InputVertexKey> SortInputVertices()
        {
            // Sort all the input vertices in the order that we wish to consider them as
            // candidate Voronoi sites.  Any sort order will produce correct output, so
            // we have complete flexibility in choosing the sort key.  We could even
            // leave them unsorted, although this would have the disadvantage that
            // changing the order of the input edges could cause S2Builder to snap to a
            // different set of Voronoi sites.
            //
            // We have chosen to sort them primarily by S2CellId since this improves the
            // performance of many S2Builder phases (due to better spatial locality).
            // It also allows the possibility of replacing the current S2PointIndex
            // approach with a more efficient recursive divide-and-conquer algorithm.
            //
            // However, sorting by leaf S2CellId alone has two small disadvantages in
            // the case where the candidate sites are densely spaced relative to the
            // snap radius (e.g., when using the IdentitySnapFunction, or when snapping
            // to E6/E7 near the poles, or snapping to S2CellId/E6/E7 using a snap
            // radius larger than the minimum value required):
            //
            //  - First, it tends to bias the Voronoi site locations towards points that
            //    are earlier on the S2CellId Hilbert curve.  For example, suppose that
            //    there are two parallel rows of input vertices on opposite sides of the
            //    edge between two large S2Cells, and the rows are separated by less
            //    than the snap radius.  Then only vertices from the cell with the
            //    smaller S2CellId are selected, because they are considered first and
            //    prevent us from selecting the sites from the other cell (because they
            //    are closer than "snap_radius" to an existing site).
            //
            //  - Second, it tends to choose more Voronoi sites than necessary, because
            //    at each step we choose the first site along the Hilbert curve that is
            //    at least "snap_radius" away from all previously selected sites.  This
            //    tends to yield sites whose "coverage discs" overlap quite a bit,
            //    whereas it would be better to cover all the input vertices with a
            //    smaller set of coverage discs that don't overlap as much.  (This is
            //    the "geometric set cover problem", which is NP-hard.)
            //
            // It is not worth going to much trouble to fix these problems, because they
            // really aren't that important (and don't affect the guarantees made by the
            // algorithm), but here are a couple of heuristics that might help:
            //
            // 1. Sort the input vertices by S2CellId at a coarse level (down to cells
            // that are O(snap_radius) in size), and then sort by a fingerprint of the
            // S2Point coordinates (i.e., quasi-randomly).  This would retain most of
            // the advantages of S2CellId sorting, but makes it more likely that we will
            // select sites that are further apart.
            //
            // 2. Rather than choosing the first uncovered input vertex and snapping it
            // to obtain the next Voronoi site, instead look ahead through following
            // candidates in S2CellId order and choose the furthest candidate whose
            // snapped location covers all previous uncovered input vertices.
            //
            // TODO(ericv): Experiment with these approaches.

            var keys = new List<InputVertexKey>(input_vertices_.Count);
            for (var i = 0; i < input_vertices_.Count; ++i)
            {
                keys.Add(new InputVertexKey(new S2CellId(input_vertices_[i]), i));
            }
            keys.Sort((InputVertexKey a, InputVertexKey b) => {
                if (a.Item1 < b.Item1) return -1;
                if (b.Item1 < a.Item1) return 1;
                return input_vertices_[a.Item2].CompareTo(input_vertices_[b.Item2]);
            });
            return keys;
        }

        // Check all edge pairs for crossings, and add the corresponding intersection
        // points to input_vertices_.  (The intersection points will be snapped and
        // merged with the other vertices during site selection.)
        private void AddEdgeCrossings(MutableS2ShapeIndex input_edge_index)
        {
            // We need to build a list of intersections and add them afterwards so that
            // we don't reallocate vertices_ during the VisitCrossings() call.
            var new_vertices = new List<S2Point>();
            input_edge_index.VisitCrossingEdgePairs(CrossingType.INTERIOR,
                (ShapeEdge a, ShapeEdge b, bool c) =>
                {
                    new_vertices.Add(S2EdgeCrossings.GetIntersection(a.V0, a.V1, b.V0, b.V1, null));
                    return true;  // Continue visiting.
                });
            if (new_vertices.Any())
            {
                snapping_needed_ = true;
                foreach (var vertex in new_vertices) AddVertex(vertex);
            }
        }
        private void AddForcedSites(S2PointIndex<int> site_index)
        {
            // Sort the forced sites and remove duplicates.
            var tmp = new SortedSet<S2Point>(sites_);
            sites_.Clear();
            sites_.AddRange(tmp);
            // Add the forced sites to the index.
            for (var id = 0; id < sites_.Count; ++id)
            {
                site_index.Add(sites_[id], id);
            }
            num_forced_sites_ = sites_.Count;
        }
        private bool IsForced(int v) => v < num_forced_sites_;
        private void ChooseInitialSites(S2PointIndex<int> site_index)
        {
            // Find all points whose distance is <= min_site_separation_ca_.
            var options = new S2ClosestPointQuery<int>.Options();
            options.ConservativeMaxDistance = (min_site_separation_ca_);
            var site_query = new S2ClosestPointQuery<int>(site_index, options);
            var results = new List<S2ClosestPointQueryBase<S1ChordAngle, int>.Result>();

            // Apply the snap_function() to each input vertex, then check whether any
            // existing site is closer than min_vertex_separation().  If not, then add a
            // new site.
            //
            // NOTE(ericv): There are actually two reasonable algorithms, which we call
            // "snap first" (the one above) and "snap last".  The latter checks for each
            // input vertex whether any existing site is closer than snap_radius(), and
            // only then applies the snap_function() and adds a new site.  "Snap last"
            // can yield slightly fewer sites in some cases, but it is also more
            // expensive and can produce surprising results.  For example, if you snap
            // the polyline "0:0, 0:0.7" using IntLatLngSnapFunction(0), the result is
            // "0:0, 0:0" rather than the expected "0:0, 0:1", because the snap radius
            // is approximately Math.Sqrt(2) degrees and therefore it is legal to snap both
            // input points to "0:0".  "Snap first" produces "0:0, 0:1" as expected.
            foreach (var key in SortInputVertices())
            {
                var vertex = input_vertices_[key.Item2];
                var site = SnapSite(vertex);
                // If any vertex moves when snapped, the output cannot be idempotent.
                snapping_needed_ = snapping_needed_ || site != vertex;

                // FindClosestPoints() measures distances conservatively, so we need to
                // recheck the distances using exact predicates.
                //
                // NOTE(ericv): When the snap radius is large compared to the average
                // vertex spacing, we could possibly avoid the call the FindClosestPoints
                // by checking whether sites_.Last() is close enough.
                var target = new S2ClosestPointQuery<int>.PointTarget(site);
                site_query.FindClosestPoints(target, results);
                var add_site = true;
                foreach (var result in results)
                {
                    if (S2Pred.CompareDistance(site, result.Point,
                                                min_site_separation_ca_) <= 0)
                    {
                        add_site = false;
                        // This pair of sites is too close.  If the sites are distinct, then
                        // the output cannot be idempotent.
                        snapping_needed_ = snapping_needed_ || site != result.Point;
                    }
                }
                if (add_site)
                {
                    site_index.Add(site, sites_.Count);
                    sites_.Add(site);
                    site_query.ReInit();
                }
            }
        }
        private S2Point SnapSite(S2Point point)
        {
            if (!snapping_requested_) return point;
            var site = Options_.SnapFunction.SnapPoint(point);
            var dist_moved = new S1ChordAngle(site, point);
            if (dist_moved > site_snap_radius_ca_)
            {
                error_ = new(S2ErrorCode.BUILDER_SNAP_RADIUS_TOO_SMALL, "Snap function moved vertex ({point.X:g}, {point.Y:g}, {point.Z:g}) by {dist_moved.ToAngle().Radians:g}, which is more than the specified snap radius of {site_snap_radius_ca_.ToAngle().Radians:g}");
            }
            return site;
        }

        // For each edge, find all sites within min_edge_site_query_radius_ca_ and
        // store them in edge_sites_.  Also, to implement idempotency this method also
        // checks whether the input vertices and edges may already satisfy the output
        // criteria.  If any problems are found then snapping_needed_ is set to true.
        private void CollectSiteEdges(S2PointIndex<int> site_index)
        {
            // Find all points whose distance is <= edge_site_query_radius_ca_.
            var options = new S2ClosestPointQuery<int>.Options();
            options.ConservativeMaxDistance = (edge_site_query_radius_ca_);
            var site_query= new S2ClosestPointQuery<int>(site_index, options);
            var results = new List<S2ClosestPointQueryBase<S1ChordAngle, int>.Result>();
            edge_sites_.Resize(input_edges_.Count, () => new List<int>());
            for (var e = 0; e < input_edges_.Count; ++e)
            {
                var edge = input_edges_[e];
                var v0 = input_vertices_[edge.Item1];
                var v1 = input_vertices_[edge.Item2];
#if s2builder_verbose
                System.Diagnostics.Debug.WriteLine($"S2Polyline: {v0.ToDebugString()}, {v1.ToDebugString()}");
#endif
                var target = new S2ClosestPointQuery<int>.EdgeTarget(v0, v1);
                site_query.FindClosestPoints(target, results);
                var sites = edge_sites_[e];
                sites.Capacity = Math.Max(results.Count, sites.Count);
                foreach (var result in results)
                {
                    sites.Add(result.Data);
                    if (!snapping_needed_ &&
                        result.Distance.IsLessThan(min_edge_site_separation_ca_limit_) &&
                        result.Point!= v0 && result.Point!= v1 &&
                        S2Pred.CompareEdgeDistance(result.Point, v0, v1, min_edge_site_separation_ca_) < 0)
                    {
                        snapping_needed_ = true;
                    }
                }
                SortSitesByDistance(v0, sites);
            }
        }
        private void SortSitesByDistance(S2Point x, List<int> sites) =>
            // Sort sites in increasing order of distance to X.
            sites.Sort((int i, int j) => S2Pred.CompareDistances(x, sites_[i], sites_[j]));

        // There are two situatons where we need to add extra Voronoi sites in order
        // to ensure that the snapped edges meet the output requirements:
        //
        //  (1) If a snapped edge deviates from its input edge by more than
        //      max_edge_deviation(), we add a new site on the input edge near the
        //      middle of the snapped edge.  This causes the snapped edge to split
        //      into two pieces, so that it follows the input edge more closely.
        //
        //  (2) If a snapped edge is closer than min_edge_vertex_separation() to any
        //      nearby site (the "site to avoid"), then we add a new site (the
        //      "separation site") on the input edge near the site to avoid.  This
        //      causes the snapped edge to follow the input edge more closely and is
        //      guaranteed to increase the separation to the required distance.
        //
        // We check these conditions by snapping all the input edges to a chain of
        // Voronoi sites and then testing each edge in the chain.  If a site needs to
        // be added, we mark all nearby edges for re-snapping.
        private void AddExtraSites(MutableS2ShapeIndex input_edge_index)
        {
            // When options_.split_crossing_edges() is true, this function may be called
            // even when site_snap_radius_ca_ == 0 (because edge_snap_radius_ca_ > 0).
            // However neither of the conditions above apply in that case.
            if (site_snap_radius_ca_ == S1ChordAngle.Zero) return;

            var chain = new List<int>();  // Temporary
            var snap_queue = new List<int>();
            for (var max_e = 0; max_e < input_edges_.Count; ++max_e)
            {
                snap_queue.Add(max_e);
                while (snap_queue.Any())
                {
                    var e = snap_queue.Last();
                    snap_queue.RemoveAt(snap_queue.Count - 1);
                    SnapEdge(e, chain);
                    // We could save the snapped chain here in a snapped_chains_ vector, to
                    // avoid resnapping it in AddSnappedEdges() below, however currently
                    // SnapEdge only accounts for less than 5% of the runtime.
                    MaybeAddExtraSites(e, max_e, chain, input_edge_index, snap_queue);
                }
            }
        }
        private void MaybeAddExtraSites(int edge_id, int max_edge_id, List<int> chain, MutableS2ShapeIndex input_edge_index, List<int> snap_queue)
        {
            // The snapped chain is always a *subsequence* of the nearby sites
            // (edge_sites_), so we walk through the two arrays in parallel looking for
            // sites that weren't snapped.  We also keep track of the current snapped
            // edge, since it is the only edge that can be too close.
            var i = 0;
            foreach (var id in edge_sites_[edge_id])
            {
                if (id == chain[i])
                {
                    if (++i == chain.Count) break;
                    // Check whether this snapped edge deviates too far from its original
                    // position.  If so, we split the edge by adding an extra site.
                    var v0 = sites_[chain[i - 1]];
                    var v1 = sites_[chain[i]];
                    if (new S1ChordAngle(v0, v1) < min_edge_length_to_split_ca_) continue;

                    var edge = input_edges_[edge_id];
                    var a0 = input_vertices_[edge.Item1];
                    var a1 = input_vertices_[edge.Item2];
                    if (!S2EdgeDistances.IsEdgeBNearEdgeA(a0, a1, v0, v1, max_edge_deviation_))
                    {
                        // Add a new site on the input edge, positioned so that it splits the
                        // snapped edge into two approximately equal pieces.  Then we find all
                        // the edges near the new site (including this one) and add them to
                        // the snap queue.
                        //
                        // Note that with large snap radii, it is possible that the snapped
                        // edge wraps around the sphere the "wrong way".  To handle this we
                        // find the preferred split location by projecting both endpoints onto
                        // the input edge and taking their midpoint.
                        var mid = (S2EdgeDistances.Project(v0, a0, a1) +
                                       S2EdgeDistances.Project(v1, a0, a1)).Normalized;
                        var new_site = GetSeparationSite(mid, v0, v1, edge_id);
                        AddExtraSite(new_site, max_edge_id, input_edge_index, snap_queue);
                        return;
                    }
                }
                else if (i > 0 && id >= num_forced_sites_)
                {
                    // Check whether this "site to avoid" is closer to the snapped edge than
                    // min_edge_vertex_separation().  Note that this is the only edge of the
                    // chain that can be too close because its vertices must span the point
                    // where "site_to_avoid" projects onto the input edge XY (this claim
                    // relies on the fact that all sites are separated by at least the snap
                    // radius).  We don't try to avoid sites added using ForceVertex()
                    // because we don't guarantee any minimum separation from such sites.
                    var site_to_avoid = sites_[id];
                    var v0 = sites_[chain[i - 1]];
                    var v1 = sites_[chain[i]];
                    if (S2Pred.CompareEdgeDistance(
                            site_to_avoid, v0, v1, min_edge_site_separation_ca_) < 0)
                    {
                        // A snapped edge can only approach a site too closely when there are
                        // no sites near the input edge near that point.  We fix that by
                        // adding a new site along the input edge (a "separation site"), then
                        // we find all the edges near the new site (including this one) and
                        // add them to the snap queue.
                        var new_site = GetSeparationSite(site_to_avoid, v0, v1, edge_id);
                        Assert.True(site_to_avoid != new_site);
                        AddExtraSite(new_site, max_edge_id, input_edge_index, snap_queue);
                        return;
                    }
                }
            }
        }

        // Adds a new site, then updates "edge_sites"_ for all edges near the new site
        // and adds them to "snap_queue" for resnapping (unless their edge id exceeds
        // "max_edge_id", since those edges have not been snapped the first time yet).
        private void AddExtraSite(S2Point new_site, int max_edge_id, MutableS2ShapeIndex input_edge_index, List<int> snap_queue)
        {
            var new_site_id = sites_.Count;
            sites_.Add(new_site);
            // Find all edges whose distance is <= edge_site_query_radius_ca_.
            var options = new S2ClosestEdgeQuery.Options();
            options.ConservativeMaxDistance = (edge_site_query_radius_ca_);
            options.IncludeInteriors = (false);
            var query = new S2ClosestEdgeQuery(input_edge_index, options);
            var target = new S2ClosestEdgeQuery.PointTarget(new_site);
            foreach (var result in query.FindClosestEdges(target))
            {
                var e = result.EdgeId;
                var site_ids = edge_sites_[e];
                site_ids.Add(new_site_id);
                SortSitesByDistance(input_vertices_[input_edges_[e].Item1], site_ids);
                if (e <= max_edge_id) snap_queue.Add(e);
            }
        }
        private S2Point GetSeparationSite(S2Point site_to_avoid, S2Point v0, S2Point v1, int input_edge_id)
        {
            // Define the "coverage disc" of a site S to be the disc centered at S with
            // radius "snap_radius".  Similarly, define the "coverage interval" of S for
            // an edge XY to be the intersection of XY with the coverage disc of S.  The
            // SnapFunction implementations guarantee that the only way that a snapped
            // edge can be closer than min_edge_vertex_separation() to a non-snapped
            // site (i.e., site_to_avoid) if is there is a gap in the coverage of XY
            // near this site.  We can fix this problem simply by adding a new site to
            // fill this gap, located as closely as possible to the site to avoid.
            //
            // To calculate the coverage gap, we look at the two snapped sites on
            // either side of site_to_avoid, and find the endpoints of their coverage
            // intervals.  The we place a new site in the gap, located as closely as
            // possible to the site to avoid.  Note that the new site may move when it
            // is snapped by the snap_function, but it is guaranteed not to move by
            // more than snap_radius and therefore its coverage interval will still
            // intersect the gap.
            var edge = input_edges_[input_edge_id];
            var x = input_vertices_[edge.Item1];
            var y = input_vertices_[edge.Item2];
            var xy_dir = y - x;
            var n = S2PointUtil.RobustCrossProd(x, y);
            var new_site = S2EdgeDistances.Project(site_to_avoid, x, y, n);
            var gap_min = GetCoverageEndpoint(v0, x, y, n);
            var gap_max = GetCoverageEndpoint(v1, y, x, -n);
            if ((new_site - gap_min).DotProd(xy_dir) < 0)
            {
                new_site = gap_min;
            }
            else if ((gap_max - new_site).DotProd(xy_dir) < 0)
            {
                new_site = gap_max;
            }
            new_site = SnapSite(new_site);
            Assert.True(v0 != new_site);
            Assert.True(v1 != new_site);
            return new_site;
        }

        // Given a site P and an edge XY with normal N, intersect XY with the disc of
        // radius snap_radius() around P, and return the intersection point that is
        // further along the edge XY toward Y.
#pragma warning disable IDE0060 // Quitar el parámetro no utilizado
        private S2Point GetCoverageEndpoint(S2Point p, S2Point x, S2Point y, S2Point n)
#pragma warning restore IDE0060 // Quitar el parámetro no utilizado
        {
            // Consider the plane perpendicular to P that cuts off a spherical cap of
            // radius snap_radius().  This plane intersects the plane through the edge
            // XY (perpendicular to N) along a line, and that line intersects the unit
            // sphere at two points Q and R, and we want to return the point R that is
            // further along the edge XY toward Y.
            //
            // Let M be the midpoint of QR.  This is the point along QR that is closest
            // to P.  We can now express R as the sum of two perpendicular vectors OM
            // and MR in the plane XY.  Vector MR is in the direction N x P, while
            // vector OM is in the direction (N x P) x N, where N = X x Y.
            //
            // The length of OM can be found using the Pythagorean theorem on triangle
            // OPM, and the length of MR can be found using the Pythagorean theorem on
            // triangle OMR.
            //
            // In the calculations below, we save some work by scaling all the vectors
            // by n.CrossProd(p).Norm2, and normalizing at the end.
            var n2 = n.Norm2;
            var nDp = n.DotProd(p);
            var nXp = n.CrossProd(p);
            var nXpXn = n2 * p - nDp * n;
            var om = Math.Sqrt(1 - edge_snap_radius_sin2_) * nXpXn;
            var mr2 = edge_snap_radius_sin2_ * n2 - nDp * nDp;

            // MR is constructed so that it points toward Y (rather than X).
            var mr = Math.Sqrt(Math.Max(0.0, mr2)) * nXp;
            return (om + mr).Normalized;
        }
        private void SnapEdge(int e, List<int> chain)
        {
            chain.Clear();
            var edge = input_edges_[e];
            if (!snapping_needed_)
            {
                chain.Add(edge.Item1);
                chain.Add(edge.Item2);
                return;
            }

            var x = input_vertices_[edge.Item1];
            var y = input_vertices_[edge.Item2];

            // Optimization: if there is only one nearby site, return.
            // Optimization: if there are exactly two nearby sites, and one is close
            // enough to each vertex, then return.

            // Now iterate through the sites.  We keep track of the sequence of sites
            // that are visited.
            var candidates = edge_sites_[e];
            foreach (var site_id in candidates)
            {
                var c = sites_[site_id];
                // Skip any sites that are too far away.  (There will be some of these,
                // because we also keep track of "sites to avoid".)  Note that some sites
                // may be close enough to the line containing the edge, but not to the
                // edge itself, so we can just use the dot product with the edge normal.
                if (S2Pred.CompareEdgeDistance(c, x, y, edge_snap_radius_ca_) > 0)
                {
                    continue;
                }
                // Check whether the new site C excludes the previous site B.  If so,
                // repeat with the previous site, and so on.
                var add_site_c = true;
                for (; chain.Any(); chain.RemoveAt(chain.Count-1))
                {
                    var b = sites_[chain.Last()];

                    // First, check whether B and C are so far apart that their clipped
                    // Voronoi regions can't intersect.
                    var bc = new S1ChordAngle(b, c);
                    if (bc >= max_adjacent_site_separation_ca_) break;

                    // Otherwise, we want to check whether site C prevents the Voronoi
                    // region of B from intersecting XY, or vice versa.  This can be
                    // determined by computing the "coverage interval" (the segment of XY
                    // intersected by the coverage disc of radius snap_radius) for each
                    // site.  If the coverage interval of one site contains the coverage
                    // interval of the other, then the contained site can be excluded.
                    var result = S2Pred.GetVoronoiSiteExclusion(
                        b, c, x, y, edge_snap_radius_ca_);
                    if (result == S2Pred.Excluded.FIRST) continue;  // Site B excluded by C
                    if (result == S2Pred.Excluded.SECOND)
                    {
                        add_site_c = false;  // Site C is excluded by B.
                        break;
                    }
                    Assert.True(S2Pred.Excluded.NEITHER == result);

                    // Otherwise check whether the previous site A is close enough to B and
                    // C that it might further clip the Voronoi region of B.
                    if (chain.Count < 2) break;
                    var a = sites_[chain[^2]];
                    var ac = new S1ChordAngle(a, c);
                    if (ac >= max_adjacent_site_separation_ca_) break;

                    // If triangles ABC and XYB have the same orientation, the circumcenter
                    // Z of ABC is guaranteed to be on the same side of XY as B.
                    var xyb = S2Pred.Sign(x, y, b);
                    if (S2Pred.Sign(a, b, c) == xyb)
                    {
                        break;  // The circumcenter is on the same side as B but further away.
                    }
                    // Other possible optimizations:
                    //  - if AB > max_adjacent_site_separation_ca_ then keep B.
                    //  - if d(B, XY) < 0.5 * Math.Min(AB, BC) then keep B.

                    // If the circumcenter of ABC is on the same side of XY as B, then B is
                    // excluded by A and C combined.  Otherwise B is needed and we can exit.
                    if (S2Pred.EdgeCircumcenterSign(x, y, a, b, c) != xyb) break;
                }
                if (add_site_c)
                {
                    chain.Add(site_id);
                }
            }
            Assert.True(chain.Any());
#if s2builder_verbose
            var sb = new StringBuilder($"({edge.Item1},{edge.Item2}): ");
            foreach (var id in chain) sb.Append(id + " ");
            System.Diagnostics.Debug.WriteLine(sb.ToString());
#endif
        }
        private const int kMinLayersForVertexFiltering = 10;

        private void BuildLayers()
        {
            // Each output edge has an "input edge id set id" (an int) representing
            // the set of input edge ids that were snapped to this edge.  The actual
            // ints can be retrieved using "input_edge_id_set_lexicon".
            var layer_edges = new List<List<OutputEdge>>();
            var layer_input_edge_ids = new List<List<int>>();
            var input_edge_id_set_lexicon = new IdSetLexicon();
            BuildLayerEdges(layer_edges, layer_input_edge_ids, input_edge_id_set_lexicon);

            // At this point we have no further need for the input geometry or nearby
            // site data, so we clear those fields to save space.
            input_vertices_.Clear();
            input_edges_.Clear();
            edge_sites_.Clear();

            // If there are a large number of layers, then we build a minimal subset of
            // vertices for each layer.  This ensures that layer types that iterate over
            // vertices will run in time proportional to the size of that layer rather
            // than the size of all layers combined.
            var layer_vertices = new List<List<S2Point>>();
            if (layers_.Count >= kMinLayersForVertexFiltering)
            {
                // Disable vertex filtering if it is disallowed by any layer.  (This could
                // be optimized, but in current applications either all layers allow
                // filtering or none of them do.)
                var allow_vertex_filtering = false;
                foreach (var options in layer_options_)
                {
                    allow_vertex_filtering &= options.AllowVertexFiltering;
                }
                if (allow_vertex_filtering)
                {
                    var filter_tmp = new List<int>();  // Temporary used by FilterVertices.
                    layer_vertices.Capacity = layers_.Count;
                    for (var i = 0; i < layers_.Count; ++i)
                    {
                        layer_vertices[i] = Graph.FilterVertices(sites_, layer_edges[i], filter_tmp);
                    }
                    sites_.Clear();  // Release memory
                }
            }
            for (var i = 0; i < layers_.Count; ++i)
            {
                var vertices = (!layer_vertices.Any() ? sites_ : layer_vertices[i]);
                var graph = new Graph(layer_options_[i], vertices, layer_edges[i],
                layer_input_edge_ids[i], input_edge_id_set_lexicon,
                label_set_ids_, label_set_lexicon_,
                layer_is_full_polygon_predicates_[i]);
                layers_[i].Build(graph, out error_);
                // Don't free the layer data until all layers have been built, in order to
                // support building multiple layers at once (e.g. ClosedSetNormalizer).
            }
        }

        // Snaps and possibly simplifies the edges for each layer, populating the
        // given output arguments.  The resulting edges can be used to construct an
        // S2Builder.Graph directly (no further processing is necessary).
        //
        // This method is not "const" because Graph.ProcessEdges can modify
        // layer_options_ in some cases (changing undirected edges to directed ones).
        private void BuildLayerEdges(List<List<OutputEdge>> layer_edges, List<List<int>> layer_input_edge_ids, IdSetLexicon input_edge_id_set_lexicon)
        {
            // OutputEdge chains are simplified only when a non-zero snap radius is specified.
            // If so, we build a map from each site to the set of input vertices that
            // snapped to that site.
            List<List<int>> site_vertices = new();
            var simplify = snapping_needed_ && Options_.SimplifyEdgeChains;
            if (simplify) site_vertices.Capacity = sites_.Count;

            while (layer_edges.Count < layers_.Count) layer_edges.Add(new());
            while (layer_input_edge_ids.Count < layers_.Count) layer_input_edge_ids.Add(new List<int>());
            for (var i = 0; i < layers_.Count; ++i)
            {
                AddSnappedEdges(layer_begins_[i], layer_begins_[i + 1], layer_options_[i],
                                layer_edges[i], layer_input_edge_ids[i],
                                input_edge_id_set_lexicon, site_vertices);
            }
            if (simplify)
            {
                SimplifyEdgeChains(site_vertices, layer_edges, layer_input_edge_ids,
                                   input_edge_id_set_lexicon);
            }
            // We simplify edge chains before processing the per-layer GraphOptions
            // because simplification can create duplicate edges and/or sibling edge
            // pairs which may need to be removed.
            for (var i = 0; i < layers_.Count; ++i)
            {
                // The errors generated by ProcessEdges are really warnings, so we simply
                // record them and continue.
                Graph.ProcessEdges(layer_options_[i], layer_edges[i],
                                    layer_input_edge_ids[i],
                                    input_edge_id_set_lexicon, out error_);
            }
        }

        // Snaps all the input edges for a given layer, populating the given output
        // arguments.  If (*site_vertices) is non-empty then it is updated so that
        // (*site_vertices)[site] contains a list of all input vertices that were
        // snapped to that site.
        private void AddSnappedEdges(int begin, int end, GraphOptions options, List<OutputEdge> edges, List<int> input_edge_ids, IdSetLexicon input_edge_id_set_lexicon, List<List<int>> site_vertices)
        {
            var discard_degenerate_edges = options.DegenerateEdges_== GraphOptions.DegenerateEdges.DISCARD;
            var chain = new List<int>();
            for (var e = begin; e < end; ++e)
            {
                var id = IdSetLexicon.AddSingleton(e);
                SnapEdge(e, chain);
                MaybeAddInputVertex(input_edges_[e].Item1, chain[0], site_vertices);
                if (chain.Count == 1)
                {
                    if (discard_degenerate_edges) continue;
                    AddSnappedEdge(chain[0], chain[0], id, options.EdgeType_,
                                   edges, input_edge_ids);
                }
                else
                {
                    MaybeAddInputVertex(input_edges_[e].Item2, chain.Last(), site_vertices);
                    for (var i = 1; i < chain.Count; ++i)
                    {
                        AddSnappedEdge(chain[i - 1], chain[i], id, options.EdgeType_,
                                       edges, input_edge_ids);
                    }
                }
            }
#if s2builder_verbose
            DumpEdges(edges, sites_);
#endif
        }
        private static void DumpEdges(List<OutputEdge> edges, List<S2Point> vertices)
        {
            foreach (var e in edges)
            {
                var v = new S2Point[2] { vertices[e.Item1], vertices[e.Item2] };
#if s2debug
                System.Diagnostics.Debug.WriteLine($"S2Polyline: {S2TextFormat.ToDebugString(v)}({e.Item1},{e.Item2})");
#endif
            }
        }

        // If "site_vertices" is non-empty, ensures that (*site_vertices)[id] contains
        // "v".  Duplicate entries are allowed.
        private static void MaybeAddInputVertex(int v, int id, List<List<int>> site_vertices)
        {
            if (!site_vertices.Any()) return;

            // Optimization: check if we just added this vertex.  This is worthwhile
            // because the input edges usually form a continuous chain, i.e. the
            // destination of one edge is the same as the source of the next edge.
            var vertices = site_vertices[id];
            if (!vertices.Any() || vertices.Last() != v)
            {
                vertices.Add(v);
            }
        }

        // Adds the given edge to "edges" and "input_edge_ids".  If undirected edges
        // are being used, also adds an edge in the opposite direction.
        private static void AddSnappedEdge(int src, int dst, int id, EdgeType edge_type, List<OutputEdge> edges, List<int> input_edge_ids)
        {
            edges.Add(new(src, dst));
            input_edge_ids.Add(id);
            if (edge_type == EdgeType.UNDIRECTED)
            {
                edges.Add(new(dst, src));
                // Automatically created edges do not have input edge ids or labels.  This
                // can be used to distinguish the original direction of the undirected edge.
                input_edge_ids.Add(IdSetLexicon.kEmptySetId);
            }
        }
        // Simplifies edge chains, updating its input/output arguments as necessary.
        private void SimplifyEdgeChains(List<List<int>> site_vertices, List<List<OutputEdge>> layer_edges, List<List<int>> layer_input_edge_ids, IdSetLexicon input_edge_id_set_lexicon)
        {
            if (!layers_.Any()) return;

            // Merge the edges from all layers (in order to build a single graph).
            List<OutputEdge> merged_edges = new();
            List<int> merged_input_edge_ids = new();
            List<int> merged_edge_layers = new();
            MergeLayerEdges(layer_edges, layer_input_edge_ids,
                            merged_edges, merged_input_edge_ids, merged_edge_layers);

            // The following fields will be reconstructed by EdgeChainSimplifier.
            foreach (var edges in layer_edges) edges.Clear();
            foreach (var input_edge_ids in layer_input_edge_ids) input_edge_ids.Clear();

            // The graph options are irrelevant for edge chain simplification, but we
            // try to set them appropriately anyway.
            var graph_options = new GraphOptions(EdgeType.DIRECTED,
                                                  GraphOptions.DegenerateEdges.KEEP,
                                                  GraphOptions.DuplicateEdges.KEEP,
                                                  GraphOptions.SiblingPairs.KEEP);
            var graph = new Graph(graph_options, sites_, merged_edges, merged_input_edge_ids,
              input_edge_id_set_lexicon, null, null,
              null);
            var simplifier = new EdgeChainSimplifier(
                this, graph, merged_edge_layers, site_vertices,
                layer_edges, layer_input_edge_ids, input_edge_id_set_lexicon);
            simplifier.Run();
        }
        // Merges the edges from all layers and sorts them in lexicographic order so
        // that we can construct a single graph.  The sort is stable, which means that
        // any duplicate edges within each layer will still be sorted by int.
        private static void MergeLayerEdges(List<List<OutputEdge>> layer_edges, List<List<int>> layer_input_edge_ids, List<OutputEdge> edges, List<int> input_edge_ids, List<int> edge_layers)
        {
            List<LayerEdgeId> order = new();
            for (var i = 0; i < layer_edges.Count; ++i)
            {
                for (var e = 0; e < layer_edges[i].Count; ++e)
                {
                    order.Add(new(i, e));
                }
            }
            order.Sort((LayerEdgeId ai, LayerEdgeId bi) => {
                return StableLessThan(layer_edges[ai.Item1][ai.Item2],
                                      layer_edges[bi.Item1][bi.Item2], ai, bi);
            });
            edges.Capacity = order.Count;
            input_edge_ids.Capacity = order.Count;
            edge_layers.Capacity = order.Count;
            foreach (var id in order)
            {
                edges.Add(layer_edges[id.Item1][id.Item2]);
                input_edge_ids.Add(layer_input_edge_ids[id.Item1][id.Item2]);
                edge_layers.Add(id.Item1);
            }
        }
        // A comparison function that allows stable sorting with sort (which is
        // fast but not stable).  It breaks ties between equal edges by comparing
        // their LayerEdgeIds.
        private static int StableLessThan(OutputEdge a, OutputEdge b, LayerEdgeId ai, LayerEdgeId bi)
        {
            // The compiler doesn't optimize this as well as it should:
            //   return (a, ai) < (b, bi);
            if (a.Item1.CompareTo(b.Item1) != 0) return a.Item1.CompareTo(b.Item1);
            if (a.Item2.CompareTo(b.Item2) != 0) return a.Item2.CompareTo(b.Item2);
            if (b.Item2.CompareTo(a.Item2) < 0) return 1;
            return ai.CompareTo(bi);  // Stable sort.
        }

        // Helper functions for computing error bounds:
        private static S1ChordAngle RoundUp(S1Angle a)
        {
            var ca = new S1ChordAngle(a);
            return ca.PlusError(ca.S1AngleConstructorMaxError);
        }
        private static S1ChordAngle AddPointToPointError(S1ChordAngle ca) => ca.PlusError(ca.GetS2PointConstructorMaxError());
        private static S1ChordAngle AddPointToEdgeError(S1ChordAngle ca) => ca.PlusError(S2EdgeDistances.GetUpdateMinDistanceMaxError(ca));

        //////////// Parameters /////////////

        // S2Builder options.
        public Options Options_ { get; set; }

        // The maximum distance (inclusive) that a vertex can move when snapped,
        // equal to S1ChordAngle(options_.snap_function().snap_radius()).
        private S1ChordAngle site_snap_radius_ca_;

        // The maximum distance (inclusive) that an edge can move when snapping to a
        // snap site.  It can be slightly larger than the site snap radius when
        // edges are being split at crossings.
        private S1ChordAngle edge_snap_radius_ca_;

        private S1Angle max_edge_deviation_;
        private S1ChordAngle edge_site_query_radius_ca_;
        private S1ChordAngle min_edge_length_to_split_ca_;

        private S1Angle min_site_separation_;
        private S1ChordAngle min_site_separation_ca_;
        private S1ChordAngle min_edge_site_separation_ca_;
        private S1ChordAngle min_edge_site_separation_ca_limit_;

        private S1ChordAngle max_adjacent_site_separation_ca_;

        // The squared sine of the edge snap radius.  This is equivalent to the snap
        // radius (squared) for distances measured through the interior of the
        // sphere to the plane containing an edge.  This value is used only when
        // interpolating new points along edges (see GetSeparationSite).
        private double edge_snap_radius_sin2_;

        // A copy of the argument to Build().
        private S2Error error_;

        // True if snapping was requested.  This is true if either snap_radius() is
        // positive, or split_crossing_edges() is true (which implicitly requests
        // snapping to ensure that both crossing edges are snapped to the
        // intersection point).
        private bool snapping_requested_;

        // Initially false, and set to true when it is discovered that at least one
        // input vertex or edge does not meet the output guarantees (e.g., that
        // vertices are separated by at least snap_function.min_vertex_separation).
        private bool snapping_needed_;

        //////////// Input Data /////////////

        // A flag indicating whether label_set_ has been modified since the last
        // time label_set_id_ was computed.
        private bool label_set_modified_;

        private List<S2Point> input_vertices_ = new();
        private readonly List<InputEdge> input_edges_ = new();

        private readonly List<Layer> layers_ = new();
        private readonly List<GraphOptions> layer_options_ = new();
        private readonly List<int> layer_begins_ = new();
        private readonly List<IsFullPolygonPredicate> layer_is_full_polygon_predicates_ = new();

        // Each input edge has "label set id" (an int) representing the set of
        // labels attached to that edge.  This vector is populated only if at least
        // one label is used.
        private readonly List<int> label_set_ids_ = new();
        private readonly IdSetLexicon label_set_lexicon_ = new();

        // The current set of labels (represented as a stack).
        private readonly List<int> label_set_ = new();

        // The LabelSetId corresponding to the current label set, computed on demand
        // (by adding it to label_set_lexicon()).
        private int label_set_id_;

        ////////////// Data for Snapping and Simplifying //////////////

        // The number of sites specified using ForceVertex().  These sites are
        // always at the beginning of the sites_ vector.
        private int num_forced_sites_;

        // The set of snapped vertex locations ("sites").
        private readonly List<S2Point> sites_ = new();

        // A map from each input edge to the set of sites "nearby" that edge,
        // defined as the set of sites that are candidates for snapping and/or
        // avoidance.  Note that compact_array will inline up to two sites, which
        // usually takes care of the vast majority of edges.  Sites are kept sorted
        // by increasing distance from the origin of the input edge.
        //
        // Once snapping is finished, this field is discarded unless edge chain
        // simplification was requested, in which case instead the sites are
        // filtered by removing the ones that each edge was snapped to, leaving only
        // the "sites to avoid" (needed for simplification).
        private readonly List<List<int>> edge_sites_ = new();

        // Indicates whether the input edges are undirected.  Typically this is
        // specified for each output layer (e.g., s2builderutil.S2PolygonLayer).
        //
        // Directed edges are preferred, since otherwise the output is ambiguous.
        // For example, output polygons may be the *inverse* of the intended result
        // (e.g., a polygon intended to represent the world's oceans may instead
        // represent the world's land masses).  Directed edges are also somewhat
        // more efficient.
        //
        // However even with undirected edges, most S2Builder layer types try to
        // preserve the input edge direction whenever possible.  Generally, edges
        // are reversed only when it would yield a simpler output.  For example,
        // S2PolygonLayer assumes that polygons created from undirected edges should
        // cover at most half of the sphere.  Similarly, S2PolylineVectorLayer
        // assembles edges into as few polylines as possible, even if this means
        // reversing some of the "undirected" input edges.
        //
        // For shapes with interiors, directed edges should be oriented so that the
        // interior is to the left of all edges.  This means that for a polygon with
        // holes, the outer loops ("shells") should be directed counter-clockwise
        // while the inner loops ("holes") should be directed clockwise.  Note that
        // S2Builder.AddPolygon() follows this convention automatically.
        public enum EdgeType { DIRECTED, UNDIRECTED };
    }
}
