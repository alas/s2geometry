namespace S2Geometry;

using DegenerateBoundaries = LaxPolygonLayer.Options.DegenerateBoundaries;
using EdgeLabelMap = Dictionary<S2Shape.Edge, LabelSet>;

public class S2BuilderUtil_LaxPolygonLayerTests(ITestOutputHelper logger)
{
    private readonly ITestOutputHelper _logger = logger;

    private static readonly DegenerateBoundaries[] kAllDegenerateBoundaries = [
        DegenerateBoundaries.DISCARD, DegenerateBoundaries.DISCARD_HOLES,
        DegenerateBoundaries.DISCARD_SHELLS, DegenerateBoundaries.KEEP];

    [Fact]
    internal void Test_LaxPolygonLayer_Empty()
    {
        foreach (var db in kAllDegenerateBoundaries)
        {
            TestLaxPolygonUnchanged("", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_Full()
    {
        foreach (var db in kAllDegenerateBoundaries)
        {
            TestLaxPolygonUnchanged("full", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_OneNormalShell()
    {
        foreach (var db in kAllDegenerateBoundaries)
        {
            TestLaxPolygonUnchanged("0:0, 0:1, 1:1", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_IsFullPolygonPredicateNotCalled()
    {
        // Test that the IsFullPolygonPredicate is not called when at least one
        // non-degenerate loop is present.
        foreach (var degenerate_boundaries in kAllDegenerateBoundaries)
        {
            S2Builder builder = new(new Options());
            S2LaxPolygonShape output = new();
            LaxPolygonLayer.Options options = new()
            {
                EdgeType = EdgeType.DIRECTED,
                DegenerateBoundaries_ = degenerate_boundaries
            };
            builder.StartLayer(new LaxPolygonLayer(output, options));
            var polygon = MakeLaxPolygonOrDie("0:0, 0:1, 1:1");
            builder.AddShape(polygon);
            // If the predicate is called, it will return an error.
            builder.AddIsFullPolygonPredicate((Graph g, out S2Error error) => IsFullPolygonUnspecified(out error));
            Assert.True(builder.Build(out var error));
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_TwoNormalShellsOneNormalHole()
    {
        // The second two loops are nested.  Note that S2LaxPolygon and S2Polygon
        // require opposite vertex orderings for holes.
        foreach (var db in kAllDegenerateBoundaries)
        {
            TestLaxPolygonUnchanged("0:1, 1:1, 0:0; " +
                                    "3:3, 3:6, 6:6, 6:3; " +
                                    "4:4, 5:4, 5:5, 4:5", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_AllDegenerateShells()
    {
        foreach (var db in new[]{DegenerateBoundaries.KEEP,
      DegenerateBoundaries.DISCARD_HOLES})
        {
            TestLaxPolygonUnchanged("1:1; 2:2, 3:3", db);
        }
        foreach (var db in new[]{DegenerateBoundaries.DISCARD,
      DegenerateBoundaries.DISCARD_SHELLS})
        {
            TestLaxPolygon("1:1; 2:2, 3:3", "", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_AllDegenerateHoles()
    {
        foreach (var db in new[]{DegenerateBoundaries.KEEP,
      DegenerateBoundaries.DISCARD_SHELLS})
        {
            TestLaxPolygonUnchanged("full; 1:1; 2:2, 3:3", db);
        }
        foreach (var db in new[]{DegenerateBoundaries.DISCARD,
      DegenerateBoundaries.DISCARD_HOLES})
        {
            TestLaxPolygon("full; 1:1; 2:2, 3:3", "full", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_SomeDegenerateShells()
    {
        string kNormal = "0:0, 0:9, 9:0; 1:1, 7:1, 1:7";
        string kInput = kNormal + "; 3:2; 2:2, 2:3";
        TestLaxPolygonUnchanged(kInput, DegenerateBoundaries.KEEP);
        TestLaxPolygonUnchanged(kInput, DegenerateBoundaries.DISCARD_HOLES);
        TestLaxPolygon(kInput, kNormal, DegenerateBoundaries.DISCARD);
        TestLaxPolygon(kInput, kNormal, DegenerateBoundaries.DISCARD_SHELLS);
    }

    [Fact]
    internal void Test_LaxPolygonLayer_SomeDegenerateHoles()
    {
        foreach (var db in new[]{DegenerateBoundaries.KEEP,
      DegenerateBoundaries.DISCARD_SHELLS})
        {
            TestLaxPolygonUnchanged("0:0, 0:9, 9:0; 1:1; 2:2, 3:3", db);
        }
        foreach (var db in new[]{DegenerateBoundaries.DISCARD,
      DegenerateBoundaries.DISCARD_HOLES})
        {
            TestLaxPolygon("0:0, 0:9, 9:0; 1:1; 2:2, 3:3", "0:0, 0:9, 9:0", db);
        }
    }

    [Fact]
    internal void Test_LaxPolygonLayer_NormalAndDegenerateShellsAndHoles()
    {
        // We start with two normal shells and one normal hole.
        string kNormal = "0:0, 0:9, 9:9, 9:0; " +
                               "0:10, 0:19, 9:19, 9:10; 1:11, 8:11, 8:18, 1:18";
        // These are the same loops augmented with degenerate interior filaments
        // (holes).  Note that one filament connects the second shell and hole
        // above, transforming them into a single loop.
        string kNormalWithDegenHoles =
            "0:0, 0:9, 1:8, 1:7, 1:8, 0:9, 9:9, 9:0; " +
            "0:10, 0:19, 9:19, 9:10, 0:10, 1:11, 8:11, 8:18, 1:18, 1:11";
        // Then we add other degenerate shells and holes, including a sibling pair
        // that connects the two shells above.
        string kDegenShells = "0:9, 0:10; 2:12; 3:13, 3:14; 20:20; 10:0, 10:1";
        string kDegenHoles = "2:5; 3:6, 3:7; 8:8";
        string kInput = kNormalWithDegenHoles + "; " +
                              kDegenShells + "; " + kDegenHoles;
        TestLaxPolygon(kInput, kNormal, DegenerateBoundaries.DISCARD);
        TestLaxPolygon(kInput, kNormal + "; " + kDegenShells,
                       DegenerateBoundaries.DISCARD_HOLES);
        TestLaxPolygon(kInput, kNormalWithDegenHoles + "; " + kDegenHoles,
                       DegenerateBoundaries.DISCARD_SHELLS);
        TestLaxPolygon(kInput, kInput, DegenerateBoundaries.KEEP);
    }

    [Fact]
    internal void Test_LaxPolygonLayer_PartialLoop()
    {
        S2Builder builder = new(new Options());
        S2LaxPolygonShape output = new();
        builder.StartLayer(new LaxPolygonLayer(output));
        builder.AddPolyline(MakePolylineOrDie("0:1, 2:3, 4:5"));
        Assert.False(builder.Build(out var error));
        Assert.Equal(S2ErrorCode.BUILDER_EDGES_DO_NOT_FORM_LOOPS, error.Code);
        Assert.True(output.IsEmpty());
    }

#if false
    // TODO(ericv): Implement validation of S2LaxPolygonShape.
    [Fact]
    internal void Test_LaxPolygonLayer_InvalidPolygon() {
        S2Builder builder = new(new S2Builder.Options());
        S2LaxPolygonShape output = new();
        LaxPolygonLayer.Options options = new();
        options.Validate = true;
        builder.StartLayer(new LaxPolygonLayer(output, options));
        builder.AddPolyline(S2TextFormat.MakePolylineOrDie("0:0, 0:10, 10:0, 10:10, 0:0"));
        Assert.False(builder.Build(out var error));
        Assert.Equal(S2ErrorCode.LOOP_SELF_INTERSECTION, error.Code);
    }
#endif

    [Fact]
    internal void Test_LaxPolygonLayer_DuplicateInputEdges()
    {
        // Check that LaxPolygonLayer removes duplicate edges in such a way that
        // degeneracies are not lost.
        S2Builder builder = new(new Options());
        S2LaxPolygonShape output = new();
        LaxPolygonLayer.Options options = new()
        {
            DegenerateBoundaries_ = DegenerateBoundaries.KEEP
        };
        builder.StartLayer(new LaxPolygonLayer(output, options));
        builder.AddShape(MakeLaxPolygonOrDie("0:0, 0:5, 5:5, 5:0"));
        builder.AddPoint(MakePointOrDie("0:0"));
        builder.AddPoint(MakePointOrDie("1:1"));
        builder.AddPoint(MakePointOrDie("1:1"));
        builder.AddShape(MakeLaxPolygonOrDie("2:2, 2:3"));
        builder.AddShape(MakeLaxPolygonOrDie("2:2, 2:3"));
        Assert.True(builder.Build(out _));
        Assert.Equal("0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 2:3",
        output.ToDebugString("; "));
    }

    [Fact]
    internal void Test_LaxPolygonLayer_EdgeLabels()
    {
        // TODO(ericv): Implement EdgeType.UNDIRECTED.
        foreach (var edge_type in new[] { EdgeType.DIRECTED })
        {
            foreach (var db in kAllDegenerateBoundaries)
            {
                // Test a polygon with normal and degenerate shells and holes.  Note
                // that this S2LaxPolygonShape has duplicate edges and is therefore not
                // valid in most contexts.
                TestEdgeLabels("1:1, 1:2; 0:0, 0:9, 9:9, 9:0; 1:2, 1:1; " +
                               "3:3, 8:3, 8:8, 3:8; 4:4; 4:5, 5:5; 4:4", edge_type, db);
            }
        }
    }

    [Fact]
    internal void Test_IndexedLaxPolygonLayer_AddsShape()
    {
        S2Builder builder = new(new Options());
        MutableS2ShapeIndex index = [];
        builder.StartLayer(new IndexedLaxPolygonLayer(index));
        string polygon_str = "0:0, 0:10, 10:0";
        builder.AddPolygon(MakePolygonOrDie(polygon_str));
        Assert.True(builder.Build(out _));
        Assert.Equal(1, index.NumShapeIds());
        var polygon = (S2LaxPolygonShape)index.Shape(0);
        Assert.Equal(polygon_str, polygon.ToDebugString());
    }

    [Fact]
    internal void Test_IndexedLaxPolygonLayer_IgnoresEmptyShape()
    {
        S2Builder builder = new(new Options());
        MutableS2ShapeIndex index = [];
        builder.StartLayer(new IndexedLaxPolygonLayer(index));
        Assert.True(builder.Build(out _));
        Assert.Equal(0, index.NumShapeIds());
    }

    private void TestLaxPolygon(string input_str, string expected_str,
        EdgeType edge_type, DegenerateBoundaries degenerate_boundaries)
    {
        _logger.WriteLine(degenerate_boundaries.ToString());
        S2Builder builder = new(new Options());
        S2LaxPolygonShape output = new();
        LaxPolygonLayer.Options options = new()
        {
            EdgeType = edge_type,
            DegenerateBoundaries_ = degenerate_boundaries
        };
        builder.StartLayer(new LaxPolygonLayer(output, options));

        var polygon = MakeLaxPolygonOrDie(input_str);
        builder.AddShape(polygon);

        // In order to construct polygons that are full except possibly for a
        // collection of degenerate holes, we must supply S2Builder with a predicate
        // that distinguishes empty polygons from full ones (modulo degeneracies).
        bool has_full_loop = false;
        for (int i = 0; i < polygon.NumLoops; ++i)
        {
            if (polygon.NumLoopVertices(i) == 0) has_full_loop = true;
        }
        builder.AddIsFullPolygonPredicate(IsFullPolygon(has_full_loop));
        Assert.True(builder.Build(out _));
        string actual_str = output.ToDebugString("; ");
        Assert.Equal(expected_str, actual_str);
    }

    private void TestLaxPolygon(string input_str,
                        string expected_str,
                        DegenerateBoundaries degenerate_boundaries)
    {
        TestLaxPolygon(input_str, expected_str, EdgeType.DIRECTED,
                       degenerate_boundaries);
#if false
// TODO(ericv): Implement.
TestLaxPolygon(input_str, expected_str, EdgeType.UNDIRECTED,
             degenerate_boundaries);
#endif
    }

    private void TestLaxPolygonUnchanged(string input_str,
                                 DegenerateBoundaries degenerate_boundaries)
    {
        TestLaxPolygon(input_str, input_str, degenerate_boundaries);
    }


    private static S2Shape.Edge GetKey(S2Shape.Edge edge, EdgeType edge_type)
    {
        // For undirected edges, sort the vertices in lexicographic order.
        if (edge_type == EdgeType.UNDIRECTED && edge.V0 > edge.V1)
        {
            edge = new S2Shape.Edge(edge.V1, edge.V0);
        }
        return edge;
    }

    private static void AddShapeWithLabels(S2Shape shape, EdgeType edge_type,
        S2Builder builder, EdgeLabelMap edge_label_map)
    {
        const int kLabelBegin = 1234;  // Arbitrary.
        for (int e = 0; e < shape.NumEdges(); ++e)
        {
            Int32 label = kLabelBegin + e;
            builder.SetLabel(label);
            // For undirected edges, reverse the direction of every other input edge.
            S2Shape.Edge edge = shape.GetEdge(e);
            if (edge_type == EdgeType.UNDIRECTED && ((e & 1) != 0))
            {
                edge = new S2Shape.Edge(edge.V1, edge.V0);
            }
            builder.AddEdge(edge.V0, edge.V1);
            edge_label_map[GetKey(edge, edge_type)].Add(label);
        }
    }

    // Converts "input_str" to an S2LaxPolygonShape, assigns labels to its edges,
    // then uses LaxPolygonLayer with the given arguments to build a new
    // S2LaxPolygonShape and verifies that all edges have the expected labels.
    // (This function does not test whether the output edges are correct.)
    private static void TestEdgeLabels(string input_str, EdgeType edge_type,
        DegenerateBoundaries degenerate_boundaries)
    {
        S2Builder builder = new(new Options());
        S2LaxPolygonShape output = new();
        LabelSetIds label_set_ids = [];
        IdSetLexicon label_set_lexicon = new();
        LaxPolygonLayer.Options options = new()
        {
            EdgeType = edge_type,
            DegenerateBoundaries_ = degenerate_boundaries
        };
        builder.StartLayer(new LaxPolygonLayer(
            output, label_set_ids, label_set_lexicon, options));

        EdgeLabelMap edge_label_map = [];
        AddShapeWithLabels(MakeLaxPolygonOrDie(input_str), edge_type,
                           builder, edge_label_map);
        Assert.True(builder.Build(out _));
        for (int i = 0; i < output.NumChains(); ++i)
        {
            for (int j = 0; j < output.GetChain(i).Length; ++j)
            {
                S2Shape.Edge edge = output.ChainEdge(i, j);
                var expected_labels = edge_label_map[GetKey(edge, edge_type)];
                Assert.Equal(expected_labels.Count,
                          label_set_lexicon.IdSet_(label_set_ids[i][j]).Count);
                Assert.True(expected_labels.SequenceEqual(label_set_lexicon.IdSet_(label_set_ids[i][j])));
            }
        }
    }
}
