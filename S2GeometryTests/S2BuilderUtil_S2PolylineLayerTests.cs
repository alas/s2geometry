namespace S2Geometry
{
    public class S2BuilderUtil_S2PolylineLayerTests
    {
        [Fact]
        public void Test_S2PolylineLayer_NoEdges()
        {
            TestS2Polyline(Array.Empty<string>(), "");
        }

        [Fact]
        public void Test_S2PolylineLayer_OneEdge()
        {
            // Even with undirected edges, S2PolylineLayer prefers to reconstruct edges
            // in their original direction.
            TestS2PolylineUnchanged("3:4, 1:1");
            TestS2PolylineUnchanged("1:1, 3:4");
        }

        [Fact]
        public void Test_S2PolylineLayer_StraightLineWithBacktracking()
        {
            TestS2PolylineUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0, 4:0");
        }

        [Fact]
        public void Test_S2PolylineLayer_EarlyWalkTerminationWithEndLoop1()
        {
            // Test that the "early walk termination" code (which is needed by
            // S2PolylineVectorLayer in order to implement idempotency) does not create
            // two polylines when it is possible to assemble the edges into one.
            //
            // This example tests a code path where the early walk termination code
            // should not be triggered at all (but was at one point due to a bug).
            Options options = new();
            options.SnapFunction = (new IntLatLngSnapFunction(2));
            TestS2Polyline(new[] { "0:0, 0:2, 0:1" }, "0:0, 0:1, 0:2, 0:1", options);
        }

        [Fact]
        public void Test_S2PolylineLayer_EarlyWalkTerminationWithEndLoop2()
        {
            // This tests a different code path where the walk is terminated early
            // (yield a polyline with one edge), and then the walk is "maximimzed" by
            // appending a two-edge loop to the end.
            TestS2Polyline(new[] { "0:0, 0:1", "0:2, 0:1", "0:1, 0:2" },
                 "0:0, 0:1, 0:2, 0:1");
        }

        [Fact]
        public void Test_S2PolylineLayer_SimpleLoop()
        {
            TestS2PolylineUnchanged("0:0, 0:5, 5:5, 5:0, 0:0");
        }

        [Fact]
        public void Test_S2PolylineLayer_ManyLoops()
        {
            // This polyline consists of many overlapping loops that keep returning to
            // the same starting vertex (2:2).  This tests whether the implementation is
            // able to assemble the polyline in the original order.
            TestS2PolylineUnchanged(
                "0:0, 2:2, 2:4, 2:2, 2:4, 4:4, 4:2, 2:2, 4:4, 4:2, 2:2, 2:0, 2:2, " +
                "2:0, 4:0, 2:2, 4:2, 2:2, 0:2, 0:4, 2:2, 0:4, 0:2, 2:2, 0:4, 2:2, " +
                "0:2, 2:2, 0:0, 0:2, 2:2, 0:0");
        }

        [Fact]
        public void Test_S2PolylineLayer_UnorderedLoops()
        {
            // This test consists of 5 squares that touch diagonally, similar to the 5
            // white squares of a 3x3 chessboard.  The edges of these squares need to be
            // reordered to assemble them into a single unbroken polyline.
            TestS2Polyline(new[]{
                "3:3, 3:2, 2:2, 2:3, 3:3",
      "1:0, 0:0, 0:1, 1:1, 1:0",
      "3:1, 3:0, 2:0, 2:1, 3:1",
      "1:3, 1:2, 0:2, 0:1, 1:3",
      "1:1, 1:2, 2:2, 2:1, 1:1",  // Central square
      },
    "3:3, 3:2, 2:2, 2:1, 3:1, 3:0, 2:0, 2:1, 1:1, 1:0, 0:0, " +
    "0:1, 1:1, 1:2, 0:2, 0:1, 1:3, 1:2, 2:2, 2:3, 3:3");
        }

        [Fact]
        public void Test_S2PolylineLayer_SplitEdges()
        {
            // Test reconstruction of a polyline where two edges have been split into
            // many pieces by crossing edges.  This example is particularly challenging
            // because (1) the edges form a loop, and (2) the first and last edges are
            // identical (but reversed).  This is designed to test the heuristics that
            // attempt to find the first edge of the input polyline.
            Options options = new();
            options.SplitCrossingEdges = (true);
            options.SnapFunction = (new IntLatLngSnapFunction(7));
            TestS2Polyline(
      new[] { "0:10, 0:0, 1:0, -1:2, 1:4, -1:6, 1:8, -1:10, -5:0, 0:0, 0:10" },
      "0:10, 0:9, 0:7, 0:5, 0:3, 0:1, 0:0, 1:0, 0:1, -1:2, 0:3, 1:4, 0:5, " +
      "-1:6, 0:7, 1:8, 0:9, -1:10, -5:0, 0:0, 0:1, 0:3, 0:5, 0:7, 0:9, 0:10",
      options);
        }

        [Fact]
        public void Test_S2PolylineLayer_SimpleEdgeLabels()
        {
            S2Builder builder = new(new Options());
            S2Polyline output = new();
            LabelSet label_set_ids = new();
            IdSetLexicon label_set_lexicon = new();
            builder.StartLayer(new S2PolylineLayer(
                output, label_set_ids, label_set_lexicon,
                new S2PolylineLayer.Options(EdgeType.UNDIRECTED)));
            builder.SetLabel(5);
            builder.AddPolyline(MakePolylineOrDie("0:0, 0:1, 0:2"));
            builder.PushLabel(7);
            builder.AddPolyline(MakePolylineOrDie("0:3, 0:2"));
            builder.ClearLabels();
            builder.AddPolyline(MakePolylineOrDie("0:3, 0:4, 0:5"));
            builder.SetLabel(11);
            builder.AddPolyline(MakePolylineOrDie("0:6, 0:5"));
            Assert.True(builder.Build(out _));
            var expected = new List<LabelSet> {
                new LabelSet { 5 }, new LabelSet { 5 }, new LabelSet { 5, 7 },
                new LabelSet { }, new LabelSet { }, new LabelSet { 11 } };
            Assert.Equal(expected.Count, label_set_ids.Count);
            for (int i = 0; i < expected.Count; ++i)
            {
                Assert.Equal(expected[i].Count,
                          label_set_lexicon.IdSet_(label_set_ids[i]).Count);
                int j = 0;
                foreach (Int32 label in label_set_lexicon.IdSet_(label_set_ids[i]))
                {
                    Assert.Equal(expected[i][j++], label);
                }
            }
        }

        [Fact]
        public void Test_S2PolylineLayer_InvalidPolyline()
        {
            S2Builder builder = new(new Options());
            S2Polyline output = new();
            S2PolylineLayer.Options options = new();
            options.Validate = (true);
            builder.StartLayer(new S2PolylineLayer(output, options));
            var vertices = new[]
            {
                new S2Point(1, 0, 0),
                new S2Point(-1, 0, 0)
            };
            S2Polyline input = new(vertices, S2Debug.DISABLE);
            builder.AddPolyline(input);
            Assert.False(builder.Build(out var error));
            Assert.Equal(S2ErrorCode.ANTIPODAL_VERTICES, error.Code);
        }

        [Fact]
        public void Test_IndexedS2PolylineLayer_AddsShape()
        {
            S2Builder builder = new(new Options());
            MutableS2ShapeIndex index = new();
            builder.StartLayer(new IndexedS2PolylineLayer(index));
            string polyline_str = "0:0, 0:10";
            builder.AddPolyline(MakePolylineOrDie(polyline_str));
            Assert.True(builder.Build(out _));
            Assert.Equal(1, index.NumShapeIds());
            S2Polyline polyline = ((S2Polyline.Shape)
                index.Shape(0)).Polyline;
            Assert.Equal(polyline_str, polyline.ToDebugString());
        }

        [Fact]
        public void Test_IndexedS2PolylineLayer_AddsEmptyShape()
        {
            S2Builder builder = new(new Options());
            MutableS2ShapeIndex index = new();
            builder.StartLayer(new IndexedS2PolylineLayer(index));
            S2Polyline line = new();
            builder.AddPolyline(line);
            Assert.True(builder.Build(out _));
            Assert.Equal(0, index.NumShapeIds());
        }

        private static void TestS2Polyline(
            string[] input_strs,
            string expected_str, EdgeType edge_type,
            Options options = null)
        {
            S2Builder builder = new(options ?? new Options());
            S2Polyline output = new();
            builder.StartLayer(new S2PolylineLayer(
                output, new S2PolylineLayer.Options(edge_type)));
            foreach (var input_str in input_strs)
            {
                builder.AddPolyline(MakePolylineOrDie(input_str));
            }
            Assert.True(builder.Build(out _));
            Assert.Equal(expected_str, output.ToDebugString());
        }

        // Convenience function that tests both directed and undirected edges.
        private static void TestS2Polyline(string[] input_strs,
            string expected_str, Options options = null)
        {
            options ??= new Options();
            TestS2Polyline(input_strs, expected_str, EdgeType.DIRECTED, options);
            TestS2Polyline(input_strs, expected_str, EdgeType.UNDIRECTED, options);
        }

        private static void TestS2PolylineUnchanged(string input_str)
        {
            TestS2Polyline(new[] { input_str }, input_str);
        }
    }
}
