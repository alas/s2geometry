using System.Collections.Generic;
using System.Linq;

namespace S2Geometry
{
	// S2ConvexHullQuery builds the convex hull of any collection of points,
	// polylines, loops, and polygons.  It returns a single convex loop.
	//
	// The convex hull is defined as the smallest convex region on the sphere that
	// contains all of your input geometry.  Recall that a region is "convex" if
	// for every pair of points inside the region, the straight edge between them
	// is also inside the region.  In our case, a "straight" edge is a geodesic,
	// i.e. the shortest path on the sphere between two points.
	//
	// Containment of input geometry is defined as follows:
	//
	//  - Each input loop and polygon is contained by the convex hull exactly
	//    (i.e., according to S2Polygon.Contains(S2Polygon)).
	//
	//  - Each input point is either contained by the convex hull or is a vertex
	//    of the convex hull. (Recall that S2Loops do not necessarily contain their
	//    vertices.)
	//
	//  - For each input polyline, the convex hull contains all of its vertices
	//    according to the rule for points above.  (The definition of convexity
	//    then ensures that the convex hull also contains the polyline edges.)
	//
	// To use this class, call the Add*() methods to add your input geometry, and
	// then call GetConvexHull().  Note that GetConvexHull() does *not* reset the
	// state; you can continue adding geometry if desired and compute the convex
	// hull again.  If you want to start from scratch, simply declare a new
	// S2ConvexHullQuery object (they are cheap to create).
	//
	// This class is not thread-safe.  There are no "const" methods.
	// 
	// This implement Andrew's monotone chain algorithm, which is a variant of the
	// Graham scan (see https://en.wikipedia.org/wiki/Graham_scan).  The time
	// complexity is O(n log n), and the space required is O(n).  In fact only the
	// call to "sort" takes O(n log n) time; the rest of the algorithm is linear.
	//
	// Demonstration of the algorithm and code:
	// en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
	public class S2ConvexHullQuery
	{
		public S2ConvexHullQuery()
		{
			bound_ = S2LatLngRect.Empty; points_ = new List<S2Point>();
		}

		// Add a point to the input geometry.
		public void AddPoint(S2Point point)
		{
			bound_.AddPoint(point);
			points_.Add(point);
		}

		// Add a polyline to the input geometry.
		public void AddPolyline(S2Polyline polyline)
		{
			bound_ = bound_.Union(polyline.GetRectBound());
			for (int i = 0; i < polyline.NumVertices; ++i)
			{
				points_.Add(polyline.Vertex(i));
			}
		}

		// Add a loop to the input geometry.
		public void AddLoop(S2Loop loop)
		{
			bound_ = bound_.Union(loop.GetRectBound());
			if (loop.IsEmptyOrFull)
			{
				// The empty and full loops consist of a single fake "vertex" that should
				// not be added to our point collection.
				return;
			}
			for (int i = 0; i < loop.NumVertices; ++i)
			{
				points_.Add(loop.Vertex(i));
			}
		}

		// Add a polygon to the input geometry.
		public void AddPolygon(S2Polygon polygon)
		{
			for (int i = 0; i < polygon.NumLoops(); ++i)
			{
				var loop = polygon.Loop(i);
				// Only loops at depth 0 can contribute to the convex hull.
				if (loop.Depth == 0)
				{
					AddLoop(loop);
				}
			}
		}

		// Compute a bounding cap for the input geometry provided.
		//
		// Note that this method does not clear the geometry; you can continue
		// adding to it and call this method again if desired.
		public S2Cap GetCapBound()
		{
			// We keep track of a rectangular bound rather than a spherical cap because
			// it is easy to compute a tight bound for a union of rectangles, whereas it
			// is quite difficult to compute a tight bound around a union of caps.
			// Also, polygons and polylines implement GetCapBound() in terms of
			// GetRectBound() for this same reason, so it is much better to keep track
			// of a rectangular bound as we go along and convert it at the end.
			//
			// TODO(ericv): We could compute an optimal bound by implementing Welzl's
			// algorithm.  However we would still need to have special handling of loops
			// and polygons, since if a loop spans more than 180 degrees in any
			// direction (i.e., if it contains two antipodal points), then it is not
			// enough just to bound its vertices.  In this case the only convex bounding
			// cap is S2Cap.Full(), and the only convex bounding loop is the full loop.
			return bound_.GetCapBound();
		}

		// Compute the convex hull of the input geometry provided.
		//
		// If there is no geometry, this method returns an empty loop containing no
		// points (see S2Loop.IsEmpty).
		//
		// If the geometry spans more than half of the sphere, this method returns a
		// full loop containing the entire sphere (see S2Loop.IsFull).
		//
		// If the geometry contains 1 or 2 points, or a single edge, this method
		// returns a very small loop consisting of three vertices (which are a
		// superset of the input vertices).
		//
		// Note that this method does not clear the geometry; you can continue
		// adding to it and call this method again if desired.
		public S2Loop GetConvexHull()
		{
			S2Cap cap = GetCapBound();
			if (cap.Height>= 1)
			{
				// The bounding cap is not convex.  The current bounding cap
				// implementation is not optimal, but nevertheless it is likely that the
				// input geometry itself is not contained by any convex polygon.  In any
				// case, we need a convex bounding cap to proceed with the algorithm below
				// (in order to construct a point "origin" that is definitely outside the
				// convex hull).
				return S2Loop.kFull;
			}
			// This code implements Andrew's monotone chain algorithm, which is a simple
			// variant of the Graham scan.  Rather than sorting by x-coordinate, instead
			// we sort the points in CCW order around an origin O such that all points
			// are guaranteed to be on one side of some geodesic through O.  This
			// ensures that as we scan through the points, each new point can only
			// belong at the end of the chain (i.e., the chain is monotone in terms of
			// the angle around O from the starting point).
			S2Point origin = cap.Center.Ortho;
			points_.Sort(new OrderedCcwAround(origin));

			// Remove duplicates.  We need to do this before checking whether there are
			// fewer than 3 points.
			var tmp = points_.Distinct().ToList();
			points_.Clear();
			points_.AddRange(tmp);

			// Special cases for fewer than 3 points.
			if (points_.Count < 3)
			{
				if (!points_.Any())
				{
					return S2Loop.kEmpty;
				}
				else if (points_.Count == 1)
				{
					return GetSinglePointLoop(points_[0]);
				}
				else
				{
					return GetSingleEdgeLoop(points_[0], points_[1]);
				}
			}

			// Verify that all points lie within a 180 degree span around the origin.
			Assert.True(S2Pred.Sign(origin, points_.First(), points_.Last()) >= 0);

			// Generate the lower and upper halves of the convex hull.  Each half
			// consists of the maximal subset of vertices such that the edge chain makes
			// only left (CCW) turns.
			var lower = new List<S2Point>();
			var upper = new List<S2Point>();
			GetMonotoneChain(lower);
			points_.Reverse();
			GetMonotoneChain(upper);

			// Remove the duplicate vertices and combine the chains.
			Assert.True(lower.First() == upper.Last());
			Assert.True(lower.Last() == upper.First());
			lower.RemoveAt(lower.Count - 1);
			upper.RemoveAt(lower.Count - 1);
			lower.AddRange(upper);
			return new S2Loop(lower);
		}

		// Iterate through the given points, selecting the maximal subset of points
		// such that the edge chain makes only left (CCW) turns.
		private void GetMonotoneChain(List<S2Point> output)
		{
			Assert.True(!output.Any());
			foreach (S2Point p in points_)
			{
				// Remove any points that would cause the chain to make a clockwise turn.
				while (output.Count >= 2 && S2Pred.Sign(output[^2], output.Last(), p) <= 0)
				{
					output.RemoveAt(output.Count - 1);
				}
				output.Add(p);
			}
		}
		private S2Loop GetSinglePointLoop(S2Point p)
		{
			// Construct a 3-vertex polygon consisting of "p" and two nearby vertices.
			// Note that Contains(p) may be false for the resulting loop (see comments
			// in header file).
			var d0 = S2PointUtil.Ortho(p);
			var d1 = p.CrossProd(d0);
			var vertices = new S2Point[3];
			vertices[0] = p;

			const double kOffset = 1e-15;

			vertices[1] = (p + kOffset * d0).Normalized;
			vertices[2] = (p + kOffset * d1).Normalized;
			return new S2Loop(vertices);
		}
		private S2Loop GetSingleEdgeLoop(S2Point a, S2Point b)
		{
			// Construct a loop consisting of the two vertices and their midpoint.
			var vertices = new S2Point[3];
			vertices[0] = a;
			vertices[1] = b;
			vertices[2] = (a + b).Normalized;
			var loop = new S2Loop(vertices);
			// The resulting loop may be clockwise, so invert it if necessary.
			loop.Normalize();
			return loop;
		}

		// A comparator for sorting points in CCW around a central point "center".
		private class OrderedCcwAround : IComparer<S2Point>
		{
			private readonly S2Point center_;
			public OrderedCcwAround(S2Point center) { center_ = center; }
			public int Compare(S2Point a, S2Point b)
			{
				// If X and Y are equal, this will return false (as desired).
				return S2Pred.Sign(center_, a, b); // > 0;
			}
		}

		private S2LatLngRect bound_;
		private readonly List<S2Point> points_;
	}
}
