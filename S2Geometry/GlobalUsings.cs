global using System.Diagnostics;

global using Edge = S2Geometry.S2ShapeUtil.ShapeEdgeId;                         // Defines an output edge.
global using OutputEdge = S2Geometry.S2ShapeUtil.ShapeEdgeId;                   // Defines an output edge.   

global using EdgeLoop = System.Collections.Generic.List<int>;          // A loop consisting of a sequence of edges.
// Identifies an edge in the graph.  Edges are numbered sequentially
// starting from zero.
global using EdgeId = int;
global using EdgePolyline = System.Collections.Generic.List<int>;
global using EdgeIdSet = System.Collections.Generic.List<int>;
global using EdgeLoops = System.Collections.Generic.List<System.Collections.Generic.List<int>>;
global using InputEdge = S2Geometry.Int32Int32;                                 // Defines an input edge.
global using InputEdgeId = int;
global using InputEdgeIdSetId = int;
global using InputEdgeLoop = System.Collections.Generic.List<int>;
global using DirectedComponent = System.Collections.Generic.List<System.Collections.Generic.List<int>>;
global using UndirectedComponent = S2Geometry.Array2<System.Collections.Generic.List<System.Collections.Generic.List<int>>>;
global using LayerEdgeId = S2Geometry.Int32Int32;                               // Identifies an output edge in a particular layer.   

global using InputVertexKey = S2Geometry.S2CellIdInt32;                         // Sort key for prioritizing input vertices.
global using VertexId = int;

global using Label = int;
global using LabelSet = System.Collections.Generic.List<int>;
global using LabelSetIds = System.Collections.Generic.List<System.Collections.Generic.List<int>>;

// InputEdgeCrossings represents all pairs of intersecting input edges and
// also certain GraphEdgeClipper state modifications (kSetInside, etc).
// It is sorted lexicographically except for entries representing state
// modifications, which are sorted by the first InputEdgeId only.
global using InputEdgeCrossings = System.Collections.Generic.List<(int, S2Geometry.CrossingInputEdge)>;
global using CrossingGraphEdgeVector = S2Geometry.Array2<S2Geometry.CrossingGraphEdge>;
global using SourceEdgeCrossing = System.ValueTuple<S2Geometry.SourceId, bool>;
// SourceEdgeCrossing represents an input edge that crosses some other
// edge; it crosses the edge from left to right iff the second parameter
// is "true".
global using SourceEdgeCrossings = System.Collections.Generic.List<(int, (S2Geometry.SourceId, bool))>; //SourceEdgeCrossing
global using SourceIdLexicon = System.Collections.Generic.Dictionary<S2Geometry.SourceId, int>; // ValueLexicon<SourceId>

global using S2ShapeIndexIdCell = S2Geometry.KeyData<S2Geometry.S2CellId, S2Geometry.S2ShapeIndexCell>;
