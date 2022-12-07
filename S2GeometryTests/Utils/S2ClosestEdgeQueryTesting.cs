namespace S2Geometry;

// An abstract class that adds edges to a MutableS2ShapeIndex for benchmarking.
internal interface IShapeIndexFactory
{
    // Requests that approximately "num_edges" edges located within the given
    // S2Cap bound should be added to "index".
    void AddEdges(S2Cap index_cap, int num_edges, MutableS2ShapeIndex index);
}

// Generates a regular loop that approximately fills the given S2Cap.
//
// Regular loops are nearly the worst case for distance calculations, since
// many edges are nearly equidistant from any query point that is not
// immediately adjacent to the loop.
internal class RegularLoopShapeIndexFactory : IShapeIndexFactory
{
    public void AddEdges(S2Cap index_cap, int num_edges, MutableS2ShapeIndex index)
    {
        index.Add(new S2Loop.Shape(S2Loop.MakeRegularLoop(
            index_cap.Center, index_cap.RadiusAngle(), num_edges)));
    }
}

// Generates a fractal loop that approximately fills the given S2Cap.
internal class FractalLoopShapeIndexFactory : IShapeIndexFactory
{
    public void AddEdges(S2Cap index_cap, int num_edges, MutableS2ShapeIndex index)
    {
        var fractal = new S2Testing.Fractal();
        fractal.SetLevelForApproxMaxEdges(num_edges);
        index.Add(new S2Loop.Shape(
            fractal.MakeLoop(S2Testing.GetRandomFrameAt(index_cap.Center),
                             index_cap.RadiusAngle())));
    }
}

// Generates a cloud of points that approximately fills the given S2Cap.
internal class PointCloudShapeIndexFactory : IShapeIndexFactory
{
    public void AddEdges(S2Cap index_cap, int num_edges, MutableS2ShapeIndex index)
    {
        var points = new List<S2Point>();
        for (int i = 0; i < num_edges; ++i)
        {
            points.Add(S2Testing.SamplePoint(index_cap));
        }
        index.Add(new S2PointVectorShape(points.ToArray()));
    }
}
