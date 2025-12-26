namespace S2Geometry;

public static partial class S2ShapeUtil
{
    // The purpose of this function is to construct polygons consisting of
    // multiple loops.  It takes as input a collection of loops whose boundaries
    // do not cross, and groups them into polygons whose interiors do not
    // intersect (where the boundary of each polygon may consist of multiple
    // loops).
    //
    // some of those islands have lakes, then the input to this function would
    // islands, and their lakes.  Each loop would actually be present twice, once
    // in each direction (see below).  The output would consist of one polygon
    // representing each lake, one polygon representing each island not including
    // islands or their lakes, and one polygon representing the rest of the world
    //
    // This method is intended for internal use; external clients should use
    // S2Builder, which has more convenient interface.
    //
    // The input consists of a set of connected components, where each component
    // consists of one or more loops.  The components must satisfy the following
    // properties:
    //
    //  - The loops in each component must form a subdivision of the sphere (i.e.,
    //    they must cover the entire sphere without overlap), except that a
    //    component may consist of a single loop if and only if that loop is
    //    degenerate (i.e., its interior is empty).
    //
    //  - The boundaries of different components must be disjoint (i.e. no
    //    crossing edges or shared vertices).
    //
    //  - No component should be empty, and no loop should have zero edges.
    //
    // The output consists of a set of polygons, where each polygon is defined by
    // the collection of loops that form its boundary.  This function does not
    // actually construct any S2Shapes; it simply identifies the loops that belong
    // to each polygon.
    public static void BuildPolygonBoundaries(List<List<S2Shape>> components, List<List<S2Shape>> polygons)
    {
        polygons.Clear();
        if (components.Count==0) return;

        // Since the loop boundaries do not cross, a loop nesting hierarchy can be
        // defined by choosing any point on the sphere as the "point at infinity".
        // Loop A then contains loop B if (1) A contains the boundary of B and (2)
        // loop A does not contain the point at infinity.
        //
        // We choose S2.Origin for this purpose.  The loop nesting hierarchy then
        // determines the face structure.  Here are the details:
        //
        // 1. Build an S2ShapeIndex of all loops that do not contain S2.Origin.
        //    This leaves at most one unindexed loop per connected component
        //    (the "outer loop").
        //
        // 2. For each component, choose a representative vertex and determine
        //    which indexed loops contain it.  The "depth" of this component is
        //    defined as the number of such loops.
        //
        // 3. Assign the outer loop of each component to the containing loop whose
        //    depth is one less.  This generates a set of multi-loop polygons.
        //
        // 4. The outer loops of all components at depth 0 become a single face.

        var index = new MutableS2ShapeIndex();
        // A map from shape.id() to the corresponding component number.
        var component_ids = new List<int>();
        var outer_loops = new List<S2Shape>();
        for (int i = 0; i < components.Count; ++i)
        {
            var component = components[i];
            foreach (var loop in component)
            {
                if (component.Count > 1 &&
                    !loop.ContainsBruteForce(S2.Origin))
                {
                    // Ownership is transferred back at the end of this function.
                    index.Add(loop);
                    component_ids.Add(i);
                }
                else
                {
                    outer_loops.Add(loop);
                }
            }
            // Check that there is exactly one outer loop in each component.
            MyDebug.Assert(i + 1 == outer_loops.Count, "Component is not a subdivision");
        }
        // Find the loops containing each component.
        var ancestors = new List<List<S2Shape>>(components.Count);
        var contains_query = index.MakeS2ContainsPointQuery();
        for (int i = 0; i < outer_loops.Count; ++i)
        {
            var loop = outer_loops[i];
            MyDebug.Assert(loop.NumEdges() > 0);
            ancestors[i] = contains_query.GetContainingShapes(loop.GetEdge(0).V0);
        }
        // Assign each outer loop to the component whose depth is one less.
        // Components at depth 0 become a single face.
        Dictionary<ValueTuple<S2Shape?>, List<S2Shape>> children = []; // btree_map
        for (int i = 0; i < outer_loops.Count; ++i)
        {
            S2Shape? ancestor = null;
            int depth = ancestors[i].Count;
            if (depth > 0)
            {
                foreach (var candidate in ancestors[i])
                {
                    if (ancestors[component_ids[candidate.Id]].Count == depth - 1)
                    {
                        MyDebug.Assert(ancestor is null);
                        ancestor = candidate;
                    }
                }
                MyDebug.Assert(ancestor is not null);
            }
            ValueTuple<S2Shape?> notNullAncestor = new(ancestor!);
            if (!children.ContainsKey(notNullAncestor)) children.Add(notNullAncestor, []);
            children[notNullAncestor].Add(outer_loops[i]);
        }
        // There is one face per loop that is not an outer loop, plus one for the
        // outer loops of components at depth 0.
        polygons.Resize(index.NumShapeIds() + 1, () => []);
        for (int i = 0; i < index.NumShapeIds(); ++i)
        {
            var polygon = polygons[i];
            var loop = index.Shape(i)!;
            ValueTuple<S2Shape?> loopKey = new(loop);
            var itr = children.TryGetValue(loopKey, out List<S2Shape>? value) ? value : null;
            if (itr is not null)
            {
                polygon = itr;
            }
            polygon.Add(loop);
        }
        polygons[^1] = children[new(null)];

        // Explicitly release the shapes from the index so they are not deleted.
        index.ReleaseAll();
    }
}
