using System;
using Xunit;

namespace S2Geometry
{
    public class S2TestingTests
    {
        [Fact]
        public void Test_TriangleFractal() => TestFractal(7, 7, 1.0);

        [Fact]
        public void Test_TriangleMultiFractal() => TestFractal(2, 6, 1.0);

        [Fact]
        public void Test_SpaceFillingFractal() => TestFractal(4, 4, 1.999);

        [Fact]
        public void Test_KochCurveFractal() => TestFractal(7, 7, Math.Log(4) / Math.Log(3));

        [Fact]
        public void Test_KochCurveMultiFractal() => TestFractal(4, 8, Math.Log(4) / Math.Log(3));

        [Fact]
        public void Test_CesaroFractal() => TestFractal(7, 7, 1.8);

        [Fact]
        public void Test_CesaroMultiFractal() => TestFractal(3, 6, 1.8);

        private static void TestFractal(int min_level, int max_level, double dimension)
        {
            // This function constructs a fractal and then computes various metrics
            // (number of vertices, total length, minimum and maximum radius) and
            // verifies that they are within expected tolerances.  Essentially this
            // directly verifies that the shape constructed *is* a fractal, i.e. the
            // total length of the curve increases exponentially with the level, while
            // the area bounded by the fractal is more or lessant.

            // The radius needs to be fairly small to avoid spherical distortions.
            const double nominal_radius = 0.001;  // radians, or about 6km
            const double kDistortionError = 1e-5;

            var fractal = new S2Testing.Fractal();
            fractal.MinLevel = (min_level);
            fractal.MaxLevel = (max_level);
            fractal.FractalDimension = (dimension);
            var frame = S2Testing.GetRandomFrame();
            var loop = fractal.MakeLoop(frame, S1Angle.FromRadians(nominal_radius));
            Assert.True(loop.IsValid());

            // If min_level and max_level are not equal, then the number of vertices and
            // the total length of the curve are subject to random variation.  Here we
            // compute an approximation of the standard deviation relative to the mean,
            // noting that most of the variance is due to the random choices about
            // whether to stop subdividing at "min_level" or not.  (The random choices
            // at higher levels contribute progressively less and less to the variance.)
            // The "relative_error" below corresponds to *one* standard deviation of
            // error; it can be increased to a higher multiple if necessary.
            //
            // Details: Let n=3*(4**min_level) and k=(max_level-min_level+1).  Each of
            // the "n" edges at min_level stops subdividing at that level with
            // probability (1/k).  This gives a binomial distribution with mean u=(n/k)
            // and standard deviation s=Math.Sqrt((n/k)(1-1/k)).  The relative error (s/u)
            // can be simplified to Math.Sqrt((k-1)/n).
            var num_levels = max_level - min_level + 1;
            var min_vertices = NumVerticesAtLevel(min_level);
            var relative_error = Math.Sqrt((num_levels - 1.0) / min_vertices);

            // "expansion_factor" is the total fractal length at level "n+1" divided by
            // the total fractal length at level "n".
            var expansion_factor = Math.Pow(4, 1 - 1 / dimension);
            var expected_num_vertices = 0;
            var expected_length_sum = 0.0;

            // "triangle_perim" is the perimeter of the original equilateral triangle
            // before any subdivision occurs.
            var triangle_perim = 3 * Math.Sqrt(3) * Math.Tan(nominal_radius);
            var min_length_sum = triangle_perim * Math.Pow(expansion_factor, min_level);
            for (var level = min_level; level <= max_level; ++level)
            {
                expected_num_vertices += NumVerticesAtLevel(level);
                expected_length_sum += Math.Pow(expansion_factor, level);
            }
            expected_num_vertices /= num_levels;
            expected_length_sum *= triangle_perim / num_levels;

            Assert.True(loop.NumVertices >= min_vertices);
            Assert.True(loop.NumVertices <= NumVerticesAtLevel(max_level));
            Assert2.Near(expected_num_vertices, loop.NumVertices,
                        relative_error * (expected_num_vertices - min_vertices));

            var center = frame.Col(2);
            var min_radius = S2.M_2_PI;
            double max_radius = 0;
            double length_sum = 0;
            for (int i = 0; i < loop.NumVertices; ++i)
            {
                // Measure the radius of the fractal in the tangent plane at "center".
                double r = Math.Tan(center.Angle(loop.Vertex(i)));
                min_radius = Math.Min(min_radius, r);
                max_radius = Math.Max(max_radius, r);
                length_sum += loop.Vertex(i).Angle(loop.Vertex(i + 1));
            }
            // kVertexError is an approximate bound on the error when computing vertex
            // positions of the fractal (due to S2.FromFrame, trig calculations, etc).
            double kVertexError = 1e-14;

            // Although min_radius_factor() is only a lower bound in general, it happens
            // to be exact (to within numerical errors) unless the dimension is in the
            // range (1.0, 1.09).
            if (dimension == 1.0 || dimension >= 1.09)
            {
                // Expect the min radius to match very closely.
                Assert2.Near(min_radius, fractal.MinRadiusRactor() * nominal_radius,
                            kVertexError);
            }
            else
            {
                // Expect the min radius to satisfy the lower bound.
                Assert.True(min_radius >=
                          fractal.MinRadiusRactor() * nominal_radius - kVertexError);
            }
            // max_radius_factor() is exact (modulo errors) for all dimensions.
            Assert2.Near(max_radius, fractal.MaxRadiusFactor() * nominal_radius,
                        kVertexError);

            Assert2.Near(expected_length_sum, length_sum,
                        relative_error * (expected_length_sum - min_length_sum) +
                        kDistortionError * length_sum);
        }

        private static int NumVerticesAtLevel(int level)
        {
            Assert.True(level >= 0 && level <= 14);  // Sanity / overflow check
            return 3 * (1 << (2 * level));      // 3*(4**level)
        }
    }
}