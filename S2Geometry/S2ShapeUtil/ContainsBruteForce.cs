namespace S2Geometry.S2ShapeUtil
{
	public static partial class S2ShapeX
	{
		// Returns true if the given shape contains the given point.  Most clients
		// should not use this method, since its running time is linear in the number
		// of shape edges.  Instead clients should create an S2ShapeIndex and use
		// S2ContainsPointQuery, since this strategy is much more efficient when many
		// points need to be tested.
		//
		// Polygon boundaries are treated as being semi-open (see S2ContainsPointQuery
		// and S2VertexModel for other options).
		//
		// CAVEAT: Typically this method is only used internally.  Its running time is
		//         linear in the number of shape edges.
		public static bool ContainsBruteForce(this S2Shape shape, S2Point point)
		{
			if (shape.Dimension() < 2) return false;

			var ref_point = shape.GetReferencePoint();
			if (ref_point.Point == point) return ref_point.Contained;

			var crosser = new S2CopyingEdgeCrosser(ref_point.Point, point);
			bool inside = ref_point.Contained;
			for (int e = 0; e < shape.NumEdges; ++e) {
				var edge = shape.GetEdge(e);
				inside ^= crosser.EdgeOrVertexCrossing(edge.V0, edge.V1);
			}
			return inside;
		}
	}
}
