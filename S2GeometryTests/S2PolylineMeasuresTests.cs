using System;
using System.Collections.Generic;
using Xunit;

namespace S2Geometry
{
    public class S2PolylineMeasuresTests
    {
        [Fact]
        public void Test_GetLengthAndCentroid_GreatCircles()
        {
            // Construct random great circles and divide them randomly into segments.
            // Then make sure that the length and centroid are correct.  Note that
            // because of the way the centroid is computed, it does not matter how
            // we split the great circle into segments.

            for (int iter = 0; iter < 100; ++iter)
            {
                // Choose a coordinate frame for the great circle.
                S2Testing.GetRandomFrame(out var x, out var y, out _);

                var lineList = new List<S2Point>();
                double theta = 0;
                while (theta < S2Constants.M_2_PI)
                {
                    lineList.Add(Math.Cos(theta) * x + Math.Sin(theta) * y);
                    theta += S2Testing.Random.RandDouble();
                }
                // Close the circle.
                lineList.Add(lineList[0]);
                var line = lineList.ToArray();
                S1Angle length = S2PolylineMeasures.GetLength(line);
                Assert.True(Math.Abs(length.Radians - S2Constants.M_2_PI) <= 2e-14);
                S2Point centroid = S2PolylineMeasures.GetCentroid(line);
                Assert.True(centroid.Norm <= 2e-14);
            }
        }
    }
}
