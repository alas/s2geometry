using System;
using Xunit;
using static S2Geometry.S2TextFormat;

namespace S2Geometry
{
    // Note that the "real" testing of these methods is in s2loop_measures_test
    // and s2polyline_measures_test.  This file only checks the handling of
    // multiple shapes and shapes of different dimensions.
    public class S2ShapeIndexMeasuresTests
    {
        [Fact]
        public void Test_GetDimension_Empty() {
            Assert.Equal(-1, S2ShapeIndexMeasures.GetDimension(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetDimension_Points() {
            Assert.Equal(0, S2ShapeIndexMeasures.GetDimension(MakeIndexOrDie("0:0 # #")));

            // Create an index with an empty point set.
            MutableS2ShapeIndex index = new();
            index.Add(new S2PointVectorShape());
            Assert.Equal(0, S2ShapeIndexMeasures.GetDimension(index));
        }

        [Fact]
        public void Test_GetDimension_PointsAndLines() {
            Assert.Equal(1, S2ShapeIndexMeasures.GetDimension(MakeIndexOrDie("0:0 # 1:1, 1:2 #")));

            // Note that a polyline with one vertex has no edges, so it is effectively
            // empty for the purpose of testing GetDimension().
            Assert.Equal(1, S2ShapeIndexMeasures.GetDimension(MakeIndexOrDie("0:0 # 1:1 #")));
        }

        [Fact]
        public void Test_GetDimension_PointsLinesAndPolygons() {
            Assert.Equal(2, S2ShapeIndexMeasures.GetDimension(MakeIndexOrDie(
                "0:0 # 1:1, 2:2 # 3:3, 3:4, 4:3")));

            Assert.Equal(2, S2ShapeIndexMeasures.GetDimension(MakeIndexOrDie("# # empty")));
        }

        [Fact]
        public void Test_GetNumPoints_Empty() {
            Assert.Equal(0, S2ShapeIndexMeasures.GetNumPoints(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetNumPoints_TwoPoints() {
            Assert.Equal(2, S2ShapeIndexMeasures.GetNumPoints(MakeIndexOrDie("0:0 | 1:0 # #")));
        }

        [Fact]
        public void Test_GetNumPoints_LineAndPolygon() {
            Assert.Equal(0, S2ShapeIndexMeasures.GetNumPoints(MakeIndexOrDie(
                "# 1:1, 1:2 # 0:3, 0:5, 2:5")));
        }

        [Fact]
        public void Test_GetLength_Empty() {
            Assert.Equal(S1Angle.Zero, S2ShapeIndexMeasures.GetLength(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetLength_TwoLines() {
            Assert.Equal(S1Angle.FromDegrees(2), S2ShapeIndexMeasures.GetLength(MakeIndexOrDie(
                "4:4 # 0:0, 1:0 | 1:0, 2:0 # 5:5, 5:6, 6:5")));
        }

        [Fact]
        public void Test_GetPerimeter_Empty() {
            Assert.Equal(S1Angle.Zero, S2ShapeIndexMeasures.GetPerimeter(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetPerimeter_DegeneratePolygon() {
            Assert2.Near(4.0, S2ShapeIndexMeasures.GetPerimeter(MakeIndexOrDie(
                "4:4 # 0:0, 1:0 | 2:0, 3:0 # 0:1, 0:2, 0:3")).GetDegrees());
        }

        [Fact]
        public void Test_GetArea_Empty() {
            Assert.Equal(0.0, S2ShapeIndexMeasures.GetArea(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetArea_TwoFullPolygons() {
            Assert.Equal(8 * Math.PI, S2ShapeIndexMeasures.GetArea(MakeIndexOrDie("# # full | full")));
        }

        [Fact]
        public void Test_GetApproxArea_Empty() {
            Assert.Equal(0.0, S2ShapeIndexMeasures.GetApproxArea(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetApproxArea_TwoFullPolygons() {
            Assert.Equal(8 * Math.PI, S2ShapeIndexMeasures.GetApproxArea(MakeIndexOrDie("# # full | full")));
        }

        [Fact]
        public void Test_GetCentroid_Empty() {
            Assert.Equal(S2Point.Empty, S2ShapeIndexMeasures.GetCentroid(MakeIndexOrDie("# #")));
        }

        [Fact]
        public void Test_GetCentroid_Points() {
            Assert.Equal(new S2Point(1, 1, 0),
                      S2ShapeIndexMeasures.GetCentroid(MakeIndexOrDie("0:0 | 0:90 # #")));
        }

        [Fact]
        public void Test_GetCentroid_Polyline() {
            // Checks that points are ignored when computing the centroid.
            Assert.True(S2.ApproxEquals(
                new S2Point(1, 1, 0),
                S2ShapeIndexMeasures.GetCentroid(MakeIndexOrDie("5:5 | 6:6 # 0:0, 0:90 #"))));
        }

        [Fact]
        public void Test_GetCentroid_Polygon() {
            // Checks that points and polylines are ignored when computing the centroid.
            Assert.True(S2.ApproxEquals(
                new S2Point(S2.M_PI_4, S2.M_PI_4, S2.M_PI_4),
                S2ShapeIndexMeasures.GetCentroid(MakeIndexOrDie("5:5 # 6:6, 7:7 # 0:0, 0:90, 90:0"))));
        }
    }
}
