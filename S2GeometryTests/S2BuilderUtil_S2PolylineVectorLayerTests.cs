namespace S2Geometry
{
    public class S2BuilderUtil_S2PolylineVectorLayerTests
    {
        [Fact]
        public void Test_S2PolylineVectorLayer_NoEdges()
        {
            TestS2PolylineVectorUnchanged(new List<string>());
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_TwoPolylines()
        {
            TestS2PolylineVectorUnchanged(new List<string> { "0:0, 1:1, 2:2", "4:4, 3:3" });
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_JoiningPolylines()
        {
            // Check that polylines are joined together when possible, even if they were
            // not adjacent in the input.  For undirected edges, the polyline direction
            // should be chosen such that the first edge of the polyline was added to
            // S2Builder before the last edge of the polyline.
            TestS2PolylineVector(new List<string> { "1:1, 2:2", "3:3, 2:2", "0:0, 1:1" },
                                 new List<string> { "3:3, 2:2", "0:0, 1:1, 2:2" }, EdgeType.DIRECTED);
            TestS2PolylineVector(new List<string> { "1:1, 2:2", "3:3, 2:2", "0:0, 1:1" },
                                 new List<string> { "3:3, 2:2, 1:1, 0:0" }, EdgeType.UNDIRECTED);
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_SegmentNetwork()
        {
            // Test a complex network of polylines that meet at shared vertices.
            TestS2PolylineVectorUnchanged(new List<string>{
      "0:0, 1:1, 2:2",
      "2:2, 2:3, 2:4",
      "2:4, 3:4, 4:4",
      "2:2, 3:2, 4:2",
      "4:2, 4:3, 4:4",
      "1:0, 2:2",
      "0:1, 2:2",
      "5:4, 4:4",
      "4:5, 4:4",
      "2:4, 2:5, 1:5, 1:4, 2:4",
      "4:2, 6:1, 5:0",  // Two nested loops
      "4:2, 7:0, 6:-1",
      "11:1, 11:0, 10:0, 10:1, 11:1"  // Isolated loop
    });
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_MultipleIntersectingWalks()
        {
            // This checks idempotency for directed edges in the case of several
            // polylines that share edges (and that even share loops).  The test
            // happens to pass for undirected edges as well.
            S2PolylineVectorLayer.Options layer_options = new();
            layer_options.PolylineType_ = (Graph.PolylineType.WALK);
            var input = new List<string>
  {
    "5:5, 5:6, 6:5, 5:5, 5:4, 5:3",
    "4:4, 5:5, 6:5, 5:6, 5:5, 5:6, 6:5, 5:5, 4:5",
    "3:5, 5:5, 5:6, 6:5, 5:5, 5:6, 6:6, 7:7",
  };
            TestS2PolylineVector(input, input, layer_options);
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_EarlyWalkTermination()
        {
            // This checks idempotency for cases where earlier polylines in the input
            // happen to terminate in the middle of later polylines.  This requires
            // building non-maximal polylines.
            S2PolylineVectorLayer.Options layer_options = new();
            layer_options.PolylineType_ = (Graph.PolylineType.WALK);
            var input = new List<string>
  {
    "0:1, 1:1",
    "1:0, 1:1, 1:2",
    "0:2, 1:2, 2:2",
    "2:1, 2:2, 2:3"
  };
            TestS2PolylineVector(input, input, layer_options);
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_InputEdgeStartsMultipleLoops()
        {
            // A single input edge is split into several segments by removing portions
            // of it, and then each of those segments becomes one edge of a loop.
            S2PolylineVectorLayer.Options layer_options = new();
            layer_options.PolylineType_ = (Graph.PolylineType.WALK);
            layer_options.SiblingPairs = (GraphOptions.SiblingPairs.DISCARD);
            Options builder_options = new();
            builder_options.SnapFunction = new IntLatLngSnapFunction(7);
            var input = new List<string>
                {
                    "0:10, 0:0",
                    "0:6, 1:6, 1:7, 0:7, 0:8",
                    "0:8, 1:8, 1:9, 0:9, 0:10",
                    "0:2, 1:2, 1:3, 0:3, 0:4",
                    "0:0, 1:0, 1:1, 0:1, 0:2",
                    "0:4, 1:4, 1:5, 0:5, 0:6",
                };
            var expected = new List<string>
                {
                    "0:1, 0:0, 1:0, 1:1, 0:1",
                    "0:3, 0:2, 1:2, 1:3, 0:3",
                    "0:5, 0:4, 1:4, 1:5, 0:5",
                    "0:7, 0:6, 1:6, 1:7, 0:7",
                    "0:9, 0:8, 1:8, 1:9, 0:9",
                };
            TestS2PolylineVector(input, expected, layer_options, builder_options);
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_ValidateFalse()
        {
            // Verifies that calling set_validate(false) does not turn off s2 debugging.
            S2PolylineVectorLayer.Options layer_options=new();
            layer_options.Validate = false;
            Assert.Equal(layer_options.s2debug_override_, S2Debug.ALLOW);
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_ValidateTrue()
        {
            // Verifies that the validate() option works.
            S2PolylineVectorLayer.Options layer_options=new();
            layer_options.Validate = true;
            Assert.Equal(layer_options.s2debug_override_, S2Debug.DISABLE);
            S2Builder builder=new(new S2Builder.Options());
            List<S2Polyline> output=new();
            builder.StartLayer(
                new S2PolylineVectorLayer(output, layer_options));
            builder.AddEdge(new S2Point(1, 0, 0), new S2Point(-1, 0, 0));
            S2Error error;
            Assert.False(builder.Build(out error));
            Assert.Equal(error.Code, S2ErrorCode.ANTIPODAL_VERTICES);
        }

        [Fact]
        public void Test_S2PolylineVectorLayer_SimpleEdgeLabels()
        {
            S2Builder builder = new(new Options());
            List<S2Polyline> output = new();
            LabelSetIds label_set_ids = new();
            IdSetLexicon label_set_lexicon = new();
            S2PolylineVectorLayer.Options layer_options = new();
            layer_options.EdgeType_ = (EdgeType.UNDIRECTED);
            layer_options.DuplicateEdges_ = (GraphOptions.DuplicateEdges.MERGE);
            builder.StartLayer(new S2PolylineVectorLayer(
                output, label_set_ids, label_set_lexicon, layer_options));
            builder.SetLabel(1);
            builder.AddPolyline(MakePolylineOrDie("0:0, 0:1, 0:2"));
            builder.SetLabel(2);
            builder.AddPolyline(MakePolylineOrDie("0:3, 0:2, 0:1"));
            builder.ClearLabels();
            builder.AddPolyline(MakePolylineOrDie("0:4, 0:5"));
            Assert.True(builder.Build(out _));
            var expected = new List<LabelSetIds>
                {
                    new LabelSetIds{
                        new LabelSet {1},
                        new LabelSet {1, 2},
                        new LabelSet {2},
                    },
                    new LabelSetIds{ new LabelSet { } },
                };
            Assert.Equal(expected.Count, label_set_ids.Count);
            for (int i = 0; i < expected.Count; ++i)
            {
                Assert.Equal(expected[i].Count, label_set_ids[i].Count);
                for (int j = 0; j < expected[i].Count; ++j)
                {
                    Assert.Equal(expected[i][j].Count,
                              label_set_lexicon.IdSet_(label_set_ids[i][j]).Count);
                    int k = 0;
                    foreach (var label in label_set_lexicon.IdSet_(label_set_ids[i][j]))
                    {
                        Assert.Equal(expected[i][j][k++], label);
                    }
                }
            }
        }

        [Fact]
        public void Test_IndexedS2PolylineVectorLayer_AddsShapes()
        {
            S2Builder builder = new(new Options());
            MutableS2ShapeIndex index = new();
            builder.StartLayer(new IndexedS2PolylineVectorLayer(index));
            string polyline0_str = "0:0, 1:1";
            string polyline1_str = "2:2, 3:3";
            builder.AddPolyline(MakePolylineOrDie(polyline0_str));
            builder.AddPolyline(MakePolylineOrDie(polyline1_str));
            Assert.True(builder.Build(out _));
            Assert.Equal(2, index.NumShapeIds());
            var polyline0 = ((S2Polyline.Shape)index.Shape(0)).Polyline;
            var polyline1 = ((S2Polyline.Shape)index.Shape(1)).Polyline;
            Assert.Equal(polyline0_str, polyline0.ToDebugString());
            Assert.Equal(polyline1_str, polyline1.ToDebugString());
        }

        // Convenience function that tests both directed and undirected edges.
        private static void TestS2PolylineVector(
            List<string> input_strs,
            List<string> expected_strs,
            S2PolylineVectorLayer.Options? layer_options = null,
            Options? builder_options = null)
        {
            layer_options ??= new S2PolylineVectorLayer.Options();
            builder_options ??= new Options();
            TestS2PolylineVector(input_strs, expected_strs, EdgeType.DIRECTED,
                                 layer_options, builder_options);
            TestS2PolylineVector(input_strs, expected_strs, EdgeType.UNDIRECTED,
                                 layer_options, builder_options);
        }

        private static void TestS2PolylineVectorUnchanged(List<string> input_strs)
        {
            TestS2PolylineVector(input_strs, input_strs);
        }

        private static void TestS2PolylineVector(
            List<string> input_strs,
            List<string> expected_strs,
            EdgeType edge_type,
            S2PolylineVectorLayer.Options? layer_options = null, // by value
            Options? builder_options = null)
        {
            layer_options ??= new S2PolylineVectorLayer.Options();
            builder_options ??= new Options();
            layer_options.EdgeType_ = (edge_type);
            S2Builder builder = new(builder_options);
            List<S2Polyline> output = new();
            builder.StartLayer(new S2PolylineVectorLayer(output, layer_options));
            foreach (var input_str in input_strs)
            {
                builder.AddPolyline(MakePolylineOrDie(input_str));
            }
            Assert.True(builder.Build(out _));
            List<string> output_strs = new();
            foreach (var polyline in output)
            {
                output_strs.Add(polyline.ToDebugString());
            }
            Assert.Equal(string.Join("; ", expected_strs),
                         string.Join("; ", output_strs));
        }
    }
}
