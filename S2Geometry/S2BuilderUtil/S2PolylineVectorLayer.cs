using System;
using System.Collections.Generic;
namespace S2Geometry.S2BuilderUtil
{
    using static S2Builder;
    using static S2Builder.GraphOptions;

    using EdgeId = Int32;
    using Label = Int32;
    using LabelSetIds = List<List<int>>;

    // A layer type that assembles edges (directed or undirected) into multiple
    // S2Polylines.  Returns an error if S2Builder found any problem with the
    // input edges; this layer type does not generate any errors of its own.
    //
    // Duplicate edges are handled correctly (e.g., if a polyline backtracks on
    // itself, or loops around and retraces some of its previous edges.)  The
    // implementation attempts to preserve the order of the input edges whenever
    // possible, so that if the input is a polyline and it is not modified by
    // S2Builder, then the output will be the same polyline even if the polyline
    // forms a loop.  However, note that this is not guaranteed when undirected
    // edges are used: for example, if the input consists of a single undirected
    // edge, then either directed edge may be returned.
    public class S2PolylineVectorLayer : S2Builder.Layer
    {
        public class Options
        {
            // Constructor that uses the default options (listed below).
            public Options()
            {
                EdgeType_ = EdgeType.DIRECTED;
                PolylineType_ = Graph.PolylineType.PATH;
                DuplicateEdges_ = DuplicateEdges.KEEP;
                sibling_pairs_ = GraphOptions.SiblingPairs.KEEP;
                Validate = false;
            }

            // Constructor that specifies the edge type.
            public Options(EdgeType edge_type)
            {
                EdgeType_ = edge_type;
                PolylineType_ = Graph.PolylineType.PATH;
                DuplicateEdges_ = DuplicateEdges.KEEP;
                sibling_pairs_ = GraphOptions.SiblingPairs.KEEP;
                Validate = false;
            }

            // Indicates whether the input edges provided to S2Builder are directed or
            // undirected.
            //
            // Directed edges should be used whenever possible to avoid ambiguity.
            // The implementation attempts to preserve the structure of directed input
            // edges whenever possible, so that if the input is a vector of disjoint
            // polylines and none of them need to be modified then the output will be
            // the same polylines in the same order.  With undirected edges, there are
            // no such guarantees.
            //
            // DEFAULT: S2Builder.EdgeType.DIRECTED
            public EdgeType EdgeType_ { get; set; }

            // Indicates whether polylines should be "paths" (which don't allow
            // duplicate vertices, except possibly the first and last vertex) or
            // "walks" (which allow duplicate vertices and edges).
            //
            // If your input consists of polylines, and you want to split them into
            // separate pieces whenever they self-intersect or cross each other, then
            // use PolylineType.PATH (and probably use split_crossing_edges()).  If
            // you don't mind if your polylines backtrack or contain loops, then use
            // PolylineType.WALK.
            //
            // DEFAULT: PolylineType.PATH
            public Graph.PolylineType PolylineType_ { get; set; }

            // Indicates whether duplicate edges in the input should be kept (KEEP) or
            // merged together (MERGE).  Note you can use edge labels to determine
            // which input edges were merged into a given output edge.
            //
            // DEFAULT: DuplicateEdges.KEEP
            public DuplicateEdges DuplicateEdges_ { get; set; }

            // Indicates whether sibling edge pairs (i.e., pairs consisting of an edge
            // and its reverse edge) should be kept (KEEP) or discarded (DISCARD).
            // For example, if a polyline backtracks on itself, the DISCARD option
            // would cause this section of the polyline to be removed.  Note that this
            // option may cause a single polyline to split into several pieces (e.g.,
            // if a polyline has a "lollipop" shape).
            //
            // REQUIRES: sibling_pairs == { DISCARD, KEEP }
            //           (the CREATE and REQUIRE options are not allowed)
            //
            // DEFAULT: SiblingPairs.KEEP
            public SiblingPairs SiblingPairs
            {
                get => sibling_pairs_;
                set
                {
                    Assert.True(value == GraphOptions.SiblingPairs.KEEP ||
                           value == GraphOptions.SiblingPairs.DISCARD);
                    sibling_pairs_ = value;
                }
            }
            private SiblingPairs sibling_pairs_;

            // If true, calls FindValidationError() on each output polyline.  If any
            // error is found, it will be returned by S2Builder.Build().
            //
            // DEFAULT: false
            public bool Validate { get; set; }
        }

        // Specifies that a vector of polylines should be constructed using the
        // given options.
        public S2PolylineVectorLayer(List<S2Polyline> polylines, Options options = null)
        {
            Init(polylines, null, null, options ?? new Options());
        }

        // Specifies that a vector of polylines should be constructed using the
        // given options, and that any labels attached to the input edges should be
        // returned in "label_set_ids" and "label_set_lexicion".
        //
        // The labels associated with the edge "polyline[i].vertex({j, j+1})" can be
        // retrieved as follows:
        //
        //   for (Int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
        public S2PolylineVectorLayer(List<S2Polyline> polylines,
                              LabelSetIds label_set_ids,
                              IdSetLexicon label_set_lexicon,
                              Options options = null)
        {
            Init(polylines, label_set_ids, label_set_lexicon, options ?? new Options());
        }

        // Layer interface:
        public override GraphOptions GraphOptions_()
        {
            return new GraphOptions(options_.EdgeType_, DegenerateEdges.DISCARD,
                                options_.DuplicateEdges_, options_.SiblingPairs);
        }
        public override void Build(Graph g, out S2Error error)
        {
            error = S2Error.OK;
            var edge_polylines = g.GetPolylines(
                options_.PolylineType_);
            polylines_.Capacity = edge_polylines.Count;
            if (label_set_ids_ != null) label_set_ids_.Capacity = edge_polylines.Count;
            var vertices = new List<S2Point>();  // Temporary storage for vertices.
            var labels = new List<Label>();  // Temporary storage for labels.
            foreach (var edge_polyline in edge_polylines)
            {
                vertices.Add(g.Vertex(g.GetEdge(edge_polyline[0]).Item1));
                foreach (EdgeId e in edge_polyline)
                {
                    vertices.Add(g.Vertex(g.GetEdge(e).Item2));
                }
                var polyline = new S2Polyline(vertices.ToArray());
                vertices.Clear();
                if (options_.Validate)
                {
                    polyline.FindValidationError(out error);
                }
                polylines_.Add(polyline);
                if (label_set_ids_ != null)
                {
                    var fetcher = new Graph.LabelFetcher(g, options_.EdgeType_);
                    var polyline_labels = new List<Int32>(edge_polyline.Count);
                    foreach (EdgeId e in edge_polyline)
                    {
                        fetcher.Fetch(e, labels);
                        polyline_labels.Add(label_set_lexicon_.Add(labels));
                    }
                    label_set_ids_.Add(polyline_labels);
                }
            }
        }

        private void Init(List<S2Polyline> polylines,
            LabelSetIds label_set_ids, IdSetLexicon label_set_lexicon,
            Options options)
        {
            Assert.True((label_set_ids == null) == (label_set_lexicon == null));
            polylines_ = polylines;
            label_set_ids_ = label_set_ids;
            label_set_lexicon_ = label_set_lexicon;
            options_ = options;
        }

        private List<S2Polyline> polylines_;
        private LabelSetIds label_set_ids_;
        private IdSetLexicon label_set_lexicon_;
        private Options options_;
    }

    // Like S2PolylineVectorLayer, but adds the polylines to a MutableS2ShapeIndex.
    public class IndexedS2PolylineVectorLayer : Layer
    {

        public IndexedS2PolylineVectorLayer(MutableS2ShapeIndex index, S2PolylineVectorLayer.Options options = null)
        { index_ = index; layer_ = new S2PolylineVectorLayer(polylines_, options ?? new S2PolylineVectorLayer.Options()); }

        public override GraphOptions GraphOptions_()
        {
            return layer_.GraphOptions_();
        }

        public override void Build(Graph g, out S2Error error)
        {
            layer_.Build(g, out error);
            if (error.IsOk)
            {
                foreach (var polyline in polylines_)
                {
                    index_.Add(new S2Polyline.OwningShape(polyline));
                }
            }
        }

        private readonly MutableS2ShapeIndex index_;
        private readonly List<S2Polyline> polylines_ = new List<S2Polyline>();
        private readonly S2PolylineVectorLayer layer_;
    }
}
