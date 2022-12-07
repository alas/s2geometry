namespace S2Geometry;

internal class S2ShapeUtil_Testing
{
    /// <summary>
    /// Verifies that all methods of the two S2Shapes return identical results,
    /// except for id() and type_tag().
    /// </summary>
    internal static void ExpectEqual(S2Shape a, S2Shape b)
    {
        Assert.True(a.NumEdges() == b.NumEdges());
        for (int i = 0; i < a.NumEdges(); ++i)
        {
            Assert.Equal(a.GetEdge(i), b.GetEdge(i));
            Assert.True(a.GetChainPosition(i) == b.GetChainPosition(i));
        }
        Assert.True(a.Dimension() == b.Dimension());
        Assert.True(a.GetReferencePoint() == b.GetReferencePoint());
        Assert.True(a.NumChains() == b.NumChains());
        for (int i = 0; i < a.NumChains(); ++i)
        {
            Assert.True(a.GetChain(i) == b.GetChain(i));
            int chain_length = a.GetChain(i).Length;
            for (int j = 0; j < chain_length; ++j)
            {
                Assert.True(a.ChainEdge(i, j) == b.ChainEdge(i, j));
            }
        }
    }

    /// <summary>
    /// Verifies that two S2ShapeIndexes have identical contents (including all the
    /// S2Shapes in both indexes).
    /// </summary>
    internal static void ExpectEqual(S2ShapeIndex a, S2ShapeIndex b)
    {
        // Check that both indexes have identical shapes.
        Assert.True(a.NumShapeIds() == b.NumShapeIds());
        for (int shape_id = 0; shape_id < a.NumShapeIds(); ++shape_id)
        {
            var a_shape = a.Shape(shape_id);
            var b_shape = b.Shape(shape_id);
            if (a_shape == null || b_shape == null)
            {
                Assert.True(a_shape == b_shape);
            }
            else
            {
                Assert.True(a_shape.Id == b_shape.Id);
                Assert.True(a_shape == b_shape);
            }
        }

        // Check that both indexes have identical cell contents.
        var a_it = a.GetNewEnumerator();
        var b_it = b.GetNewEnumerator();
        var aHasNext = a_it.MoveNext();
        var bHasNext = b_it.MoveNext();
        while (aHasNext && bHasNext)
        {
            Assert.True(a_it.Current.Item1 == b_it.Current.Item1);
            var a_cell = a_it.Current.Item2;
            var b_cell = b_it.Current.Item2;
            Assert.True(a_cell.NumClipped() == b_cell.NumClipped());
            for (var i = 0; i < a_cell.NumClipped(); ++i)
            {
                var a_clipped = a_cell.Clipped(i);
                var b_clipped = b_cell.Clipped(i);
                Assert.True(a_clipped.ShapeId == b_clipped.ShapeId);
                Assert.True(a_clipped.ContainsCenter == b_clipped.ContainsCenter);
                Assert.True(a_clipped.NumEdges == b_clipped.NumEdges);
                for (int j = 0; j < a_clipped.NumEdges; ++j)
                {
                    Assert.True(a_clipped.Edge(j) == b_clipped.Edge(j));
                }
            }
            aHasNext = a_it.MoveNext();
            bHasNext = b_it.MoveNext();
        }
        Assert.True(!bHasNext);

        // Spot-check the other iterator methods.  (We know that both indexes have
        // the same contents, so any differences are due to implementation bugs.)
        a_it.Reset();
        b_it.Reset();
        a_it.MoveNext();
        b_it.MoveNext();
        Assert.True(a_it.Current.Item1 == b_it.Current.Item1);
        aHasNext = a_it.MoveNext();
        bHasNext = b_it.MoveNext();
        if (aHasNext)
        {
            Assert.True(a_it.Current.Item1 == b_it.Current.Item1);
            Assert.True(aHasNext == bHasNext);
            // Assert.True(a_it.MovePrevious());
            // Assert.True(b_it.MovePrevious());
            Assert.True(a_it.Current.Item1 == b_it.Current.Item1);
        }
        // Assert.False(a_it.MovePrevious());
        // Assert.False(b_it.MovePrevious());
        // a_it.Finish();
        // b_it.Finish();
        // Assert.True(a_it.id() == b_it.id());
        // a_it.Seek(a_it.id().Next);
        // b_it.Seek(b_it.id().Next);
        // Assert.True(a_it.id() == b_it.id());
    }
}
