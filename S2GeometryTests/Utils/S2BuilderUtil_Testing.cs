using System;
using System.Collections.Generic;
using System.Linq;
using Edge = S2Geometry.KeyKey<System.Int32, System.Int32>;                // Defines an output edge.

namespace S2Geometry
{
    // A class that copies an Graph and owns the underlying data
    // (unlike Graph, which is just a view).
    public class GraphClone
    {
        public GraphClone() { }  // Must call Init().
        public GraphClone(S2Builder.Graph g) { Init(g); }
        public void Init(S2Builder.Graph g)
        {
            options_ = g.Options;
            vertices_ = g.Vertices;
            edges_ = g.Edges;
            input_edge_id_set_ids_ = g.InputEdgeIdSetIds;
            input_edge_id_set_lexicon_ = g.InputEdgeIdSetLexicon;
            label_set_ids_ = g.LabelSetIds;
            label_set_lexicon_ = g.LabelSetLexicon;
            is_full_polygon_predicate_ = g.IsFullPolygonPredicate();
            g_ = new S2Builder.Graph(
                options_, vertices_, edges_, input_edge_id_set_ids_,
                input_edge_id_set_lexicon_, label_set_ids_, label_set_lexicon_,
                is_full_polygon_predicate_);
        }
        public S2Builder.Graph Graph() => g_;

        private S2Builder.GraphOptions options_;
        private List<S2Point> vertices_;
        private List<Edge> edges_;
        private List<Int32> input_edge_id_set_ids_;
        private IdSetLexicon input_edge_id_set_lexicon_;
        private List<Int32> label_set_ids_;
        private IdSetLexicon label_set_lexicon_;
        private S2Builder.IsFullPolygonPredicate is_full_polygon_predicate_;
        private S2Builder.Graph g_;
    }

    // A layer type that copies an Graph into a GraphClone object
    // (which owns the underlying data, unlike Graph itself).
    public class GraphCloningLayer : S2Builder.Layer
    {
        private readonly GraphClone gc_;

        public GraphCloningLayer(S2Builder.GraphOptions graph_options, GraphClone gc)
        { graph_options_ = graph_options; gc_ = gc; }

        public override S2Builder.GraphOptions GraphOptions_() => graph_options_;
        private readonly S2Builder.GraphOptions graph_options_;

        public override void Build(S2Builder.Graph g, out S2Error error)
        { error = S2Error.OK; gc_.Init(g); }
    }

    // A layer type that copies an Graph and appends it to a vector,
    // and appends the corresponding GraphClone object (which owns the Graph data)
    // to a separate vector.
    public class GraphAppendingLayer : S2Builder.Layer
    {
        private readonly List<S2Builder.Graph> graphs_;
        private readonly List<GraphClone> clones_;

        public GraphAppendingLayer(S2Builder.GraphOptions graph_options, List<S2Builder.Graph> graphs, List<GraphClone> clones)
        { graph_options_ = graph_options; graphs_ = graphs; clones_ = clones; }

        public override S2Builder.GraphOptions GraphOptions_() => graph_options_;
        private readonly S2Builder.GraphOptions graph_options_;

        public override void Build(S2Builder.Graph g, out S2Error error)
        {
            error = S2Error.OK;
            clones_.Add(new GraphClone(g));
            graphs_.Add(clones_.Last().Graph());
        }
    }
}
