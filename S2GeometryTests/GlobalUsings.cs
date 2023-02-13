global using Xunit;
global using Xunit.Abstractions;

global using S2Geometry.S2BuilderUtil;

global using static S2Geometry.S2TextFormat;
global using static S2Geometry.S2Builder;

global using S2GeometryTests.Utils;

global using Edge = S2Geometry.S2ShapeUtil.ShapeEdgeId; // Defines an output edge.

global using EdgeId = System.Int32;
global using InputEdgeId = System.Int32;
global using VertexId = System.Int32;
global using EdgeVector = System.Collections.Generic.List<S2Geometry.S2ShapeUtil.ShapeEdgeId>; //Edge
global using InputEdgeIdSetId = System.Int32;

// A set of (edge string, Int32[]) pairs representing the
// InputEdgeIds attached to the edges of a graph.  Edges are in
// S2TextFormat.ToDebugString() format, such as "1:3, 4:5".
global using EdgeInputEdgeIds = System.Collections.Generic.List<(string, System.Int32[])>;

// Since we don't expect to have any crossing edges, the key for each edge is
// simply the sum of its endpoints.  This key has the advantage of being
// unchanged when the endpoints of an edge are swapped.
global using EdgeLabelMap = System.Collections.Generic.Dictionary<S2Geometry.Vector3<double>, System.Collections.Generic.List<System.Int32>>;

global using DirectedComponent = System.Collections.Generic.List<System.Collections.Generic.List<System.Int32>>;
global using UndirectedComponent = S2Geometry.Array2<System.Collections.Generic.List<System.Collections.Generic.List<System.Int32>>>;

global using Label = System.Int32;
global using LabelSetId = System.Int32;
global using LabelSet = System.Collections.Generic.List<System.Int32>;
global using LabelSetIds = System.Collections.Generic.List<System.Collections.Generic.List<int>>;

global using S2Point = S2Geometry.Vector3<double>;
//global using S2Point_Coder = S2Geometry.S2Coder_Testing.S2BasicCoder<S2Geometry.Vector3<double>>;

global using S2EdgeCrosser = S2Geometry.S2EdgeCrosserBase;
global using S2CopyingEdgeCrosser = S2Geometry.S2EdgeCrosserBase;

global using S2MinDistance = S2Geometry.S1ChordAngle;
global using S2MinDistanceTarget = S2Geometry.S2DistanceTarget<S2Geometry.S1ChordAngle>;
