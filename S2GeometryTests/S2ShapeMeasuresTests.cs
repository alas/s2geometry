using System;
using System.Collections.Generic;
using Xunit;
using static S2Geometry.S2TextFormat;

namespace S2Geometry
{
    // Note that the "real" testing of these methods is in s2loop_measures_test
    // and s2polyline_measures_test.  This file only checks the handling of shapes
    // of different dimensions and shapes with multiple edge chains.
    public class S2ShapeMeasuresTests
    {
        [Fact]
        public void Test_GetLength_WrongDimension()
        {
            Assert.Equal(S1Angle.Zero, S2ShapeMeasures.GetLength(MakeIndexOrDie("0:0 # #").Shape(0)));
            Assert.Equal(S1Angle.Zero,
                      S2ShapeMeasures.GetLength(MakeLaxPolygonOrDie("0:0, 0:1, 1:0")));
        }

        [Fact]
        public void Test_GetLength_NoPolylines()
        {
            Assert.Equal(S1Angle.Zero, S2ShapeIndexMeasures.GetLength(MakeIndexOrDie("")));
        }

        [Fact]
        public void Test_GetLength_ThreePolylinesInOneShape()
        {
            // S2EdgeVectorShape is the only standard S2Shape that can have more than
            // one edge chain of dimension 1.
            var p = ParsePointsOrDie("0:0, 1:0, 2:0, 3:0");
            S2EdgeVectorShape shape = new S2EdgeVectorShape(new List<(S2Point, S2Point)> { (p[0], p[1]), (p[0], p[2]), (p[0], p[3]) });
            Assert.Equal(S1Angle.FromDegrees(6), S2ShapeMeasures.GetLength(shape));
        }

        [Fact]
        public void Test_GetPerimeter_WrongDimension()
        {
            Assert.Equal(S1Angle.Zero,
                      S2ShapeMeasures.GetPerimeter(MakeIndexOrDie("0:0 # #").Shape(0)));
            Assert.Equal(S1Angle.Zero,
                      S2ShapeIndexMeasures.GetPerimeter(MakeIndexOrDie("0:0, 0:1, 1:0")));
        }

        [Fact]
        public void Test_GetPerimeter_EmptyPolygon()
        {
            Assert.Equal(S1Angle.Zero, S2ShapeMeasures.GetPerimeter(MakeLaxPolygonOrDie("empty")));
        }

        [Fact]
        public void Test_GetPerimeter_FullPolygon()
        {
            Assert.Equal(S1Angle.Zero, S2ShapeMeasures.GetPerimeter(MakeLaxPolygonOrDie("full")));
        }

        [Fact]
        public void Test_GetPerimeter_TwoLoopPolygon()
        {
            // To ensure that all edges are 1 degree long, we use degenerate loops.
            Assert.Equal(S1Angle.FromDegrees(6),
                      S2ShapeMeasures.GetPerimeter(MakeLaxPolygonOrDie("0:0, 1:0; 0:1, 0:2, 0:3")));
        }

        [Fact]
        public void Test_GetArea_WrongDimension()
        {
            Assert.Equal(0.0, S2ShapeMeasures.GetArea(MakeIndexOrDie("0:0 # #").Shape(0)));
            Assert.Equal(0.0, S2ShapeIndexMeasures.GetArea(MakeIndexOrDie("0:0, 0:1, 1:0")));
        }

        [Fact]
        public void Test_GetArea_EmptyPolygon()
        {
            Assert.Equal(0.0, S2ShapeMeasures.GetArea(MakeLaxPolygonOrDie("empty")));
        }

        [Fact]
        public void Test_GetArea_FullPolygon()
        {
            Assert.Equal(4 * Math.PI, S2ShapeMeasures.GetArea(MakeLaxPolygonOrDie("full")));
            Assert.Equal(4 * Math.PI,
                      S2ShapeMeasures.GetArea(new S2Polygon.OwningShape(MakePolygonOrDie("full"))));
        }

        [Fact]
        public void Test_GetArea_TwoTinyShells()
        {
            // Two small squares with sides about 10 um (micrometers) long.
            double side = S1Angle.FromDegrees(1e-10).Radians;
            Assert.Equal(2 * side * side, S2ShapeMeasures.GetArea(MakeLaxPolygonOrDie(
                "0:0, 0:1e-10, 1e-10:1e-10, 1e-10:0; " +
                "0:0, 0:-1e-10, -1e-10:-1e-10, -1e-10:0")));
        }

        [Fact]
        public void Test_GetArea_TinyShellAndHole()
        {
            // A square about 20 um on each side with a hole 10 um on each side.
            double side = S1Angle.FromDegrees(1e-10).Radians;
            Assert.Equal(3 * side * side, S2ShapeMeasures.GetArea(MakeLaxPolygonOrDie(
                "0:0, 0:2e-10, 2e-10:2e-10, 2e-10:0; " +
                "0.5e-10:0.5e-10, 1.5e-10:0.5e-10, 1.5e-10:1.5e-10, 0.5e-10:1.5e-10")));
        }

        [Fact]
        public void Test_GetApproxArea_LargeShellAndHolePolygon()
        {
            // Make sure that GetApproxArea works well for large polygons.
            Assert2.Near(S2ShapeMeasures.GetApproxArea(MakeLaxPolygonOrDie(
                "0:0, 0:90, 90:0; 0:22.5, 90:0, 0:67.5")),
                S2Constants.M_PI_4, 1e-12);
        }

        [Fact]
        public void Test_GetCentroid_Points()
        {
            // GetCentroid() returns the centroid multiplied by the number of points.
            Assert.Equal(new S2Point(1, 1, 0),
                S2ShapeMeasures.GetCentroid(MakeIndexOrDie("0:0 | 0:90 # #").Shape(0)));
        }

        [Fact]
        public void Test_GetCentroid_Polyline()
        {
            // GetCentroid returns the centroid multiplied by the length of the polyline.
            Assert.True(S2PointUtil.ApproxEquals(
                new S2Point(1, 1, 0),
                S2ShapeIndexMeasures.GetCentroid(MakeIndexOrDie("0:0, 0:90"))));
        }

        [Fact]
        public void Test_GetCentroid_Polygon()
        {
            // GetCentroid returns the centroid multiplied by the area of the polygon.
            Assert.True(S2PointUtil.ApproxEquals(
                new S2Point(S2Constants.M_PI_4, S2Constants.M_PI_4, S2Constants.M_PI_4),
                S2ShapeMeasures.GetCentroid(MakeLaxPolygonOrDie("0:0, 0:90, 90:0"))));
        }
    }
}
