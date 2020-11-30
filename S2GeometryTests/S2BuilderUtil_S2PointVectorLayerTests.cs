using System;
using System.Collections.Generic;
using Xunit;
using S2Geometry.S2BuilderUtil;
using static S2Geometry.S2TextFormat;

namespace S2Geometry
{
    using LabelSet = List<int>;

    public class S2BuilderUtil_S2PointVectorLayerTests
    {
        [Fact]
        public void Test_S2PointVectorLayer_MergeDuplicates() {
            S2Builder builder=new(new S2Builder.Options());
            List<S2Point> output=new();
            IdSetLexicon label_set_lexicon = new();
            LabelSet label_set_ids = new();
            builder.StartLayer(new S2PointVectorLayer(
                output, label_set_ids, label_set_lexicon,
                new S2PointVectorLayer.Options(
                    S2Builder.GraphOptions.DuplicateEdges.MERGE)));

            builder.SetLabel(1);
            AddPoint(MakePointOrDie("0:1"), builder);
            AddPoint(MakePointOrDie("0:2"), builder);
            builder.SetLabel(2);
            AddPoint(MakePointOrDie("0:1"), builder);
            AddPoint(MakePointOrDie("0:4"), builder);
            AddPoint(MakePointOrDie("0:5"), builder);
            builder.ClearLabels();
            AddPoint(MakePointOrDie("0:5"), builder);
            AddPoint(MakePointOrDie("0:6"), builder);
            Assert.True(builder.Build(out _));

            Int32[][] expected_labels = { new[]{ 1, 2 }, new[] { 1 }, new[] { 2 }, new[] { 2 }, Array.Empty<Int32>() };
            string expected_points = "0:1, 0:2, 0:4, 0:5, 0:6";

            VerifyS2PointVectorLayerResults(label_set_ids, label_set_lexicon, output,
                                            expected_points, expected_labels);
        }

        [Fact]
        public void Test_S2PointVectorLayer_KeepDuplicates() {
            S2Builder builder=new(new S2Builder.Options());
            List<S2Point> output=new();
            IdSetLexicon label_set_lexicon=new();
            LabelSet label_set_ids=new();
            builder.StartLayer(new S2PointVectorLayer(
                output, label_set_ids, label_set_lexicon,
                new S2PointVectorLayer.Options(
                    S2Builder.GraphOptions.DuplicateEdges.KEEP)));

            builder.SetLabel(1);
            AddPoint(MakePointOrDie("0:1"), builder);
            AddPoint(MakePointOrDie("0:2"), builder);
            builder.SetLabel(2);
            AddPoint(MakePointOrDie("0:1"), builder);
            AddPoint(MakePointOrDie("0:4"), builder);
            AddPoint(MakePointOrDie("0:5"), builder);
            builder.ClearLabels();
            AddPoint(MakePointOrDie("0:5"), builder);
            AddPoint(MakePointOrDie("0:6"), builder);
            Assert.True(builder.Build(out _));

            Int32[][] expected_labels = { new[]{ 1 }, new[] { 2 }, new[] { 1 }, new[] { 2 }, new[] { 2 }, Array.Empty<int>(), Array.Empty<int>() };
            string expected_points = "0:1, 0:1, 0:2, 0:4, 0:5, 0:5, 0:6";

            VerifyS2PointVectorLayerResults(label_set_ids, label_set_lexicon, output,
                                            expected_points, expected_labels);
        }

        [Fact]
        public void Test_S2PointVectorLayer_Error() {
            S2Builder builder=new(new S2Builder.Options());
            List<S2Point> output=new();
            builder.StartLayer(new S2PointVectorLayer(
                output, new S2PointVectorLayer.Options(
                             S2Builder.GraphOptions.DuplicateEdges.KEEP)));

            AddPoint(MakePointOrDie("0:1"), builder);
            builder.AddEdge(MakePointOrDie("0:3"), MakePointOrDie("0:4"));
            AddPoint(MakePointOrDie("0:5"), builder);
            Assert.False(builder.Build(out var error));
            Assert.Equal(S2ErrorCode.INVALID_ARGUMENT, error.Code);
            Assert.Equal("Found non-degenerate edges", error.Text);

            Assert.Equal(2, output.Count);
            Assert.Equal(MakePointOrDie("0:1"), output[0]);
            Assert.Equal(MakePointOrDie("0:5"), output[1]);
        }

        [Fact]
        public void Test_IndexedS2PointVectorLayer_AddsShapes() {
            S2Builder builder=new(new S2Builder.Options());
            MutableS2ShapeIndex index=new();
            builder.StartLayer(new IndexedS2PointVectorLayer(index));
            string point0_str = "0:0";
            string point1_str = "2:2";
            builder.AddPoint(S2TextFormat.MakePointOrDie(point0_str));
            builder.AddPoint(S2TextFormat.MakePointOrDie(point1_str));
            Assert.True(builder.Build(out _));
            Assert.Equal(1, index.NumShapeIds());
            var shape = (S2PointVectorShape) index.Shape(0);
            Assert.Equal(2, shape.NumPoints);
            Assert.Equal(point0_str, shape.Point(0).ToDebugString());
            Assert.Equal(point1_str, shape.Point(1).ToDebugString());
        }

        [Fact]
        public void Test_IndexedS2PointVectorLayer_AddsEmptyShape() {
            S2Builder builder=new(new S2Builder.Options());
            MutableS2ShapeIndex index=new();
            builder.StartLayer(new IndexedS2PointVectorLayer(index));
            Assert.True(builder.Build(out _));
            Assert.Equal(0, index.NumShapeIds());
        }

        private static void VerifyS2PointVectorLayerResults(
            LabelSet label_set_ids, IdSetLexicon label_set_lexicon, List<S2Point> output,
            string str_expected_points, Int32[][] expected_labels)
        {
            var expected_points = S2TextFormat.ParsePointsOrDie(str_expected_points);

            Assert.Equal(expected_labels.Length, label_set_ids.Count);
            for (int i = 0; i < output.Count; ++i)
            {
                Assert.Equal(expected_points[i], output[i]);
                Assert.Equal(expected_labels[i].Length,
                          label_set_lexicon.IdSet_(label_set_ids[i]).Count);
                int k = 0;
                foreach (var label in label_set_lexicon.IdSet_(label_set_ids[i]))
                {
                    Assert.Equal(expected_labels[i][k++], label);
                }
            }
        }

        private static void AddPoint(S2Point p, S2Builder builder) { builder.AddEdge(p, p); }
    } 
}
