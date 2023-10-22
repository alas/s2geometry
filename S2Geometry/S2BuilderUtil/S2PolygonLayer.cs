// Note that there are two supported output types for polygons: S2Polygon and
// S2LaxPolygonShape.  Use S2Polygon if you need the full range of operations
// that S2Polygon implements.  Use S2LaxPolygonShape if you want to represent
// polygons with zero-area degenerate regions, or if you need a type that has
// low memory overhead and fast initialization.  However, be aware that to
// convert from a S2LaxPolygonShape to an S2Polygon you will need to use
// S2Builder again.
//
// Similarly, there are two supported output formats for polygon meshes:
// S2LaxPolygonShapeVector and S2PolygonMesh.  Use S2PolygonMesh if you need
// to be able to determine which polygons are adjacent to each edge or vertex;
// otherwise use S2LaxPolygonShapeVector, which uses less memory and is faster
// to construct.

// A layer type that assembles edges (directed or undirected) into an
// S2Polygon.  Returns an error if the edges cannot be assembled into loops.
//
// If the input edges are directed, they must be oriented such that the
// polygon interior is to the left of all edges.  Directed edges are always
// preferred (see S2Builder.EdgeType).
//
// Before the edges are assembled into loops, "sibling pairs" consisting of an
// edge and its reverse edge are automatically removed.  Such edge pairs
// represent zero-area degenerate regions, which S2Polygon does not allow.
// (If you need to build polygons with degeneracies, use LaxPolygonLayer
// instead.)
//
// S2PolygonLayer is implemented such that if the input to S2Builder is a
// polygon and is not modified, then the output has the same cyclic ordering
// of loop vertices and the same loop ordering as the input polygon.
//
// If the polygon has no edges, then the graph's IsFullPolygonPredicate is
// called to determine whether the output polygon should be empty (containing
// no points) or full (containing all points).  This predicate can be
// specified as part of the S2Builder input geometry.

namespace S2Geometry.S2BuilderUtil;

using static S2Builder;
using static S2Builder.GraphOptions;
using LoopType = S2Builder.Graph.LoopType;
using LoopMap = Dictionary<S2Loop, (int, bool)>; // gtl.btree_map

public class S2PolygonLayer : Layer
{
    public class Options
    {
        // Constructor that uses the default options (listed below).
        public Options()
        {
            EdgeType = EdgeType.DIRECTED; Validate = false;
        }

        // Constructor that specifies the edge type.
        public Options(EdgeType edge_type)
        {
            EdgeType = edge_type; Validate = false;
        }

        // Indicates whether the input edges provided to S2Builder are directed or
        // undirected.  Directed edges should be used whenever possible (see
        // S2Builder.EdgeType for details).
        //
        // If the input edges are directed, they should be oriented so that the
        // polygon interior is to the left of all edges.  This means that for a
        // polygon with holes, the outer loops ("shells") should be directed
        // counter-clockwise while the inner loops ("holes") should be directed
        // clockwise.  Note that S2Builder.AddPolygon() does this automatically.
        //
        // DEFAULT: S2Builder.EdgeType.DIRECTED
        public EdgeType EdgeType { get; set; }

        // If true, calls FindValidationError() on the output polygon.  If any
        // error is found, it will be returned by S2Builder.Build().
        //
        // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
        // order to turn off the default error checking in debug builds.
        //
        // DEFAULT: false
        public bool Validate { get; set; }
    }

    // Specifies that a polygon should be constructed using the given options.
    public S2PolygonLayer(S2Polygon polygon, Options? options = null)
    {
        Init(polygon, null, null, options ?? new Options());
    }

    // Specifies that a polygon should be constructed using the given options,
    // and that any labels attached to the input edges should be returned in
    // "label_set_ids" and "label_set_lexicion".
    //
    // The labels associated with the edge "polygon.loop(i).vertex({j, j+1})"
    // can be retrieved as follows:
    //
    //   for (Int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
    public S2PolygonLayer(S2Polygon polygon, LabelSetIds label_set_ids,
                   IdSetLexicon label_set_lexicon,
                   Options? options = null)
    {
        Init(polygon, label_set_ids, label_set_lexicon, options ?? new Options());
    }

    // Layer interface:
    public override GraphOptions GraphOptions_()
    {
        // Prevent degenerate edges and sibling edge pairs.  There should not be any
        // duplicate edges if the input is valid, but if there are then we keep them
        // since this tends to produce more comprehensible errors.
        return new GraphOptions(options_.EdgeType, DegenerateEdges.DISCARD,
                            DuplicateEdges.KEEP, SiblingPairs.DISCARD);
    }

    public override void Build(Graph g, out S2Error error)
    {
        label_set_ids_?.Clear();

        // It's tricky to compute the edge labels for S2Polygons because the
        // S2Polygon.Init methods can reorder and/or invert the loops.  We handle
        // this by remembering the original vector index of each loop and whether or
        // not the loop contained S2.Origin.  By comparing this with the final
        // S2Polygon loops we can fix up the edge labels appropriately.
        var loop_map = new LoopMap();
        if (g.NumEdges == 0)
        {
            // The polygon is either full or empty.
            if (g.IsFullPolygon(out error))
            {
                polygon_ = new S2Polygon(S2Loop.KFull);
            }
            else
            {
                polygon_ = new S2Polygon([]);
            }
        }
        else if (g.Options.EdgeType_ == EdgeType.DIRECTED)
        {
            var edge_loops = new EdgeLoops();
            if (!g.GetDirectedLoops(LoopType.SIMPLE, edge_loops, out error))
            {
                return;
            }
            var loops = new List<S2Loop>();
            AppendS2Loops(g, edge_loops, loops);
            AppendEdgeLabels(g, edge_loops);
            edge_loops.Clear();  // Release memory
            InitLoopMap(loops, loop_map);
            polygon_ = new S2Polygon(loops);
        }
        else
        {
            var components = new List<UndirectedComponent>();
            if (!g.GetUndirectedComponents(LoopType.SIMPLE, components, out error))
            {
                return;
            }
            // It doesn't really matter which complement of each component we use,
            // since below we normalize all the loops so that they enclose at most
            // half of the sphere (to ensure that the loops can always be nested).
            //
            // The only reason to prefer one over the other is that when there are
            // multiple loops that touch, only one of the two complements matches the
            // structure of the input loops.  GetUndirectedComponents() tries to
            // ensure that this is always complement 0 of each component.
            var loops = new List<S2Loop>();
            foreach (var component in components)
            {
                AppendS2Loops(g, component[0], loops);
                AppendEdgeLabels(g, component[0]);
            }
            components.Clear();  // Release memory
            InitLoopMap(loops, loop_map);
            foreach (var loop in loops) loop.Normalize();
            polygon_ = new S2Polygon(loops);
        }
        ReorderEdgeLabels(loop_map);
        if (options_.Validate)
        {
            polygon_.FindValidationError(out error);
        }
    }

    private void Init(S2Polygon polygon, LabelSetIds? label_set_ids, IdSetLexicon? label_set_lexicon, Options options)
    {
        MyDebug.Assert((label_set_ids is null) == (label_set_lexicon is null));
        polygon_ = polygon;
        label_set_ids_ = label_set_ids;
        label_set_lexicon_ = label_set_lexicon;
        options_ = options;

        if (options_.Validate)
        {
            polygon_.S2DebugOverride = S2Debug.DISABLE;
        }
    }

    private static void AppendS2Loops(Graph g, EdgeLoops edge_loops, List<S2Loop> loops)
    {
        var vertices = new List<S2Point>();
        foreach (var edge_loop in edge_loops)
        {
            vertices.Capacity = edge_loop.Count;
            foreach (var edge_id in edge_loop)
            {
                vertices.Add(g.Vertex(g.GetEdge(edge_id).ShapeId));
            }
            loops.Add(new S2Loop(vertices));
            vertices.Clear();
        }
    }
    private void AppendEdgeLabels(Graph g, EdgeLoops edge_loops)
    {
        if (label_set_ids_ is null) return;

        var labels = new InputEdgeLoop();  // Temporary storage for labels.
        var fetcher = new Graph.LabelFetcher(g, options_.EdgeType);
        foreach (var edge_loop in edge_loops)
        {
            var loop_label_set_ids = new InputEdgeLoop(edge_loop.Count);
            foreach (var edge_id in edge_loop)
            {
                fetcher.Fetch(edge_id, labels);
                loop_label_set_ids.Add(label_set_lexicon_.Add(labels));
            }
            label_set_ids_.Add(loop_label_set_ids);
        }
    }
    private void InitLoopMap(List<S2Loop> loops, LoopMap loop_map)
    {
        if (label_set_ids_ is null) return;
        for (int i = 0; i < loops.Count; i++)
        {
            S2Loop loop = loops[i];
            loop_map[loop] = (i, loop.ContainsOrigin);
        }
    }
    private void ReorderEdgeLabels(LoopMap loop_map)
    {
        if (label_set_ids_ is null) return;

        var new_ids = new LabelSetIds(label_set_ids_.Count);
        for (int i = 0; i < polygon_.NumLoops(); ++i)
        {
            var loop = polygon_.Loop(i);
            var old = loop_map[loop];
            (label_set_ids_[old.Item1], new_ids[i])=(new_ids[i], label_set_ids_[old.Item1]);
            if (loop.ContainsOrigin != old.Item2)
            {
                // S2Loop.Invert() reverses the order of the vertices, which leaves
                // the last edge unchanged.  For example, the loop ABCD (with edges
                // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
                new_ids[i].Reverse(0, new_ids[i].Count - 1);
            }
        }
        label_set_ids_ = new_ids;
    }

    private S2Polygon polygon_;
    private LabelSetIds? label_set_ids_;
    private IdSetLexicon? label_set_lexicon_;
    private Options options_;
}

// Like S2PolygonLayer, but adds the polygon to a MutableS2ShapeIndex (if the
// polygon is non-empty).
public class IndexedS2PolygonLayer : Layer
{
    public IndexedS2PolygonLayer(MutableS2ShapeIndex index, S2PolygonLayer.Options? options = null)
    {
        index_ = index; polygon_ = new S2Polygon();
        layer_ = new S2PolygonLayer(polygon_, options ?? new S2PolygonLayer.Options());
    }

    public override GraphOptions GraphOptions_()
    {
        return layer_.GraphOptions_();
    }

    public override void Build(Graph g, out S2Error error)
    {
        layer_.Build(g, out error);
        if (error.IsOk() && !polygon_.IsEmpty())
        {
            index_.Add(new S2Polygon.OwningShape(polygon_));
        }
    }

    private readonly MutableS2ShapeIndex index_;
    private readonly S2Polygon polygon_;
    private readonly S2PolygonLayer layer_;
}
