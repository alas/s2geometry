namespace S2Geometry;

using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Text;

public static class S2TextFormat
{
    private static List<string> SplitString(string s, char separator) => (
        from t in s.Split(separator,
            StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
        select t).ToList();

    #region ParseLatLngs

    // Parses a string of one or more latitude-longitude coordinates in degrees,
    // and return the corresponding vector of S2LatLng points.
    // Examples of the input format:
    //     ""                            // no points
    //     "-20:150"                     // one point
    //     "-20:150, -20:151, -19:150"   // three points
    public static List<S2LatLng> ParseLatLngsOrDie(string str)
    {
        var latlngs = new List<S2LatLng>();
        MyDebug.Assert(ParseLatLngs(str, latlngs));
        return latlngs;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool ParseLatLngs(string str, List<S2LatLng> latlngs)
    {
        if (!ParseUtils.DictionaryParse(str, out var ps)) return false;

        foreach (var (lat, lng) in ps)
        {
            if (!double.TryParse(lat, NumberStyles.Any, CultureInfo.InvariantCulture, out var latDeg)) return false;

            if (!double.TryParse(lng, NumberStyles.Any, CultureInfo.InvariantCulture, out var lngDeg)) return false;

            latlngs.Add(S2LatLng.FromDegrees(latDeg, lngDeg));
        }
        return true;
    }

    #endregion

    // Parses a string in the same format as ParseLatLngs, and return the
    // corresponding vector of S2Point values.
    public static List<S2Point> ParsePointsOrDie(string str)
    {
        var vertices = new List<S2Point>();
        MyDebug.Assert(ParsePoints(str, vertices));
        return vertices;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool ParsePoints(string str, List<S2Point> vertices)
    {
        vertices.Clear();
        List<S2LatLng> latlngs = [];
        if (!ParseLatLngs(str, latlngs)) return false;

        vertices.AddRange(latlngs.Select(t => t.ToPoint()).ToArray());
        return true;
    }

    public static S2Point MakePointOrDie(string str)
    {
        MyDebug.Assert(MakePoint(str, out var point));
        return point;
    }

    // As above, but do not Assert.True-fail on invalid input. Returns true if conversion
    // is successful.
    public static bool MakePoint(string str, out S2Point point)
    {
        point = S2Point.Empty;
        var vertices = new List<S2Point>();
        if (!ParsePoints(str, vertices) || vertices.Count != 1) return false;

        point = vertices[0];
        return true;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeLatLng(string str, out S2LatLng latlng)
    {
        latlng = S2LatLng.Invalid;
        var latlngs = new List<S2LatLng>();
        if (!ParseLatLngs(str, latlngs) || latlngs.Count != 1) return false;

        latlng = latlngs[0];
        return true;
    }

    // Given a string in the same format as ParseLatLngs, returns a single S2LatLng.
    public static S2LatLng MakeLatLngOrDie(string str)
    {
        MyDebug.Assert(MakeLatLng(str, out var latlng));
        return latlng;
    }

    // Given a string in the same format as ParseLatLngs, returns the minimal
    // bounding S2LatLngRect that contains the coordinates.
    public static S2LatLngRect? MakeLatLngRectOrDie(string str)
    {
        MyDebug.Assert(MakeLatLngRect(str, out var rect));
        return rect;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeLatLngRect(string str, [NotNullWhen(true)] out S2LatLngRect? rect)
    {
        rect = null;
        var latlngs = new List<S2LatLng>();
        if (!ParseLatLngs(str, latlngs) || latlngs.Count==0) return false;

        rect = S2LatLngRect.FromPoint(latlngs[0]);
        foreach (var ll in latlngs.Skip(1))
        {
            rect = rect.Value.AddPoint(ll);
        }
        return true;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeCellId(string str, out S2CellId cellId)
    {
        cellId = S2CellId_FromDebugString(str);
        return cellId != S2CellId.None;
    }

    // Parses an S2CellId in the format "f/dd..d" where "f" is a digit in the
    // range [0-5] representing the S2CellId face, and "dd..d" is a string of
    // digits in the range [0-3] representing each child's position with respect
    // to its parent.  (Note that the latter string may be empty.)
    //
    // For example "4/" represents S2CellId.FromFace(4), and "3/02" represents
    // S2CellId.FromFace(3).child(0).child(2).
    //
    // This function is a wrapper for S2CellId.FromDebugString().
    public static S2CellId MakeCellIdOrDie(string str)
    {
        MyDebug.Assert(MakeCellId(str, out var cell_id));
        return cell_id;
    }

    /// <summary>
    /// Converts a string in the format returned by ToString() to an S2CellId.
    /// Returns S2CellId.None if the string could not be parsed.
    /// 
    /// The method name includes "Debug" in order to avoid possible confusion
    /// with FromToken() above.
    /// </summary>
    private static S2CellId S2CellId_FromDebugString(string str)
    {
        // This function is reasonably efficient, but is only intended for use in
        // tests.
        var level = str.Length - 2;
        if (level < 0 || level > S2.kMaxCellLevel) return S2CellId.None;

        var face = (int)char.GetNumericValue(str[0]);
        if (face < 0 || face > 5 || str[1] != '/') return S2CellId.None;

        var id = S2CellId.FromFace(face);
        for (var i = 2; i < str.Length; i++)
        {
            var child_pos = (int)char.GetNumericValue(str[i]);
            if (child_pos < 0 || child_pos > 3) return S2CellId.None;
            id = id.Child(child_pos);
        }
        return id;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeCellUnion(string str, [NotNullWhen(true)] out S2CellUnion? cell_union)
    {
        cell_union = null;
        var cellIds = new List<S2CellId>();
        foreach (var cell_str in SplitString(str, ','))
        {
            if (!MakeCellId(cell_str, out var cellId)) return false;
            cellIds.Add(cellId);
        }
        cell_union = new S2CellUnion(cellIds);
        return true;
    }

    // Parses a comma-separated list of S2CellIds in the format above, and returns
    // the corresponding S2CellUnion.  (Note that S2CellUnions are automatically
    // normalized by sorting, removing duplicates, and replacing groups of 4 child
    // cells by their parent cell.)
    public static S2CellUnion MakeCellUnionOrDie(string str)
    {
        MyDebug.Assert(MakeCellUnion(str, out var cell_union));
        return cell_union;
    }

    // Given a string of latitude-longitude coordinates in degrees,
    // returns a newly allocated loop.  Example of the input format:
    //     "-20:150, 10:-120, 0.123:-170.652"
    // The strings "empty" or "full" create an empty or full loop respectively.
    public static S2Loop MakeLoopOrDie(string str, S2Debug override_ = S2Debug.ALLOW)
    {
        MyDebug.Assert(MakeLoop(str, out var loop, override_));
        return loop!;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeLoop(string str, [NotNullWhen(true)] out S2Loop? loop, S2Debug override_ = S2Debug.ALLOW)
    {
        if (str == "empty")
        {
            loop = S2Loop.KEmpty;
            return true;
        }

        if (str == "full")
        {
            loop = S2Loop.KFull;
            return true;
        }

        loop = null;
        List<S2Point> vertices = [];
        if (!ParsePoints(str, vertices)) return false;
        loop = new S2Loop(vertices, override_);
        return true;
    }

    // Similar to MakeLoop(), but returns an S2Polyline rather than an S2Loop.
    public static S2Polyline MakePolylineOrDie(string str, S2Debug override_ = S2Debug.ALLOW)
    {
        MyDebug.Assert(MakePolyline(str, out var polyline, override_));
        return polyline!.Value;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakePolyline(string str, [NotNullWhen(true)] out S2Polyline? polyline, S2Debug override_ = S2Debug.ALLOW)
    {
        polyline = null;
        var vertices = new List<S2Point>();
        if (!ParsePoints(str, vertices)) return false;
        polyline = new S2Polyline([.. vertices], override_);
        return true;
    }

    // Like MakePolyline, but returns an S2LaxPolylineShape instead.
    public static S2LaxPolylineShape MakeLaxPolylineOrDie(string str)
    {
        MyDebug.Assert(MakeLaxPolyline(str, out var lax_polyline));
        return lax_polyline;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeLaxPolyline(string str, [NotNullWhen(true)] out S2LaxPolylineShape? lax_polyline)
    {
        lax_polyline = null;
        List<S2Point> vertices = [];
        if (!ParsePoints(str, vertices)) return false;
        lax_polyline = new S2LaxPolylineShape(vertices);
        return true;
    }

    private static bool InternalMakePolygon(string str, bool normalize_loops, [NotNullWhen(true)]out S2Polygon? polygon)
    {
        polygon = null;
        if (str == "empty") str = "";
        var loop_strs = SplitString(str, ';');
        List<S2Loop> loops = [];
        foreach (var loop_str in loop_strs)
        {
            if (!MakeLoop(loop_str, out var loop)) return false;
            // Don't normalize loops that were explicitly specified as "full".
            if (normalize_loops && !loop.IsFull()) loop.Normalize();
            loops.Add(loop);
        }
        polygon = new S2Polygon(loops);
        return true;
    }

    // Given a sequence of loops separated by semicolons, returns a newly
    // allocated polygon.  Loops are automatically normalized by inverting them
    // if necessary so that they enclose at most half of the unit sphere.
    // (Historically this was once a requirement of polygon loops.  It also
    // hides the problem that if the user thinks of the coordinates as X:Y
    // rather than LAT:LNG, it yields a loop with the opposite orientation.)
    //
    // Examples of the input format:
    //     "10:20, 90:0, 20:30"                                  // one loop
    //     "10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2:20.3"   // two loops
    //     ""       // the empty polygon (consisting of no loops)
    //     "empty"  // the empty polygon (consisting of no loops)
    //     "full"   // the full polygon (consisting of one full loop).
    public static S2Polygon MakePolygonOrDie(string str)
    {
        MyDebug.Assert(MakePolygon(str, out var polygon));
        return polygon;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakePolygon(string str, [NotNullWhen(true)]out S2Polygon? polygon)
    {
        return InternalMakePolygon(str, true, out polygon);
    }

    // Like MakePolygon(), except that it does not normalize loops (i.e., it
    // gives you exactly what you asked for).
    public static S2Polygon MakeVerbatimPolygonOrDie(string str)
    {
        MyDebug.Assert(MakeVerbatimPolygon(str, out var polygon));
        return polygon;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeVerbatimPolygon(string str, [NotNullWhen(true)] out S2Polygon? polygon)
    {
        return InternalMakePolygon(str, false, out polygon);
    }

    // Parses a string in the same format as MakePolygon, except that loops must
    // be oriented so that the interior of the loop is always on the left, and
    // polygons with degeneracies are supported.  As with MakePolygon, "full" and
    // denotes the full polygon and "" or "empty" denote the empty polygon.
    public static S2LaxPolygonShape MakeLaxPolygonOrDie(string str)
    {
        MyDebug.Assert(MakeLaxPolygon(str, out var lax_polygon));
        return lax_polygon;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeLaxPolygon(string str, [NotNullWhen(true)] out S2LaxPolygonShape? lax_polygon)
    {
        lax_polygon = null;
        var loop_strs = SplitString(str, ';');
        List<List<S2Point>> loops = [];
        foreach (var loop_str in loop_strs)
        {
            if (loop_str == "full")
            {
                loops.Add(new(0));
            }
            else if (loop_str != "empty")
            {
                List<S2Point> points = [];
                if (!ParsePoints(loop_str, points)) return false;
                loops.Add(points);
            }
        }
        lax_polygon = new S2LaxPolygonShape(loops);
        return true;
    }

    // Returns a MutableS2ShapeIndex containing the points, polylines, and loops
    // (in the form of one polygon for each group of loops) described by the
    // following format:
    //
    //   point1|point2|... # line1|line2|... # polygon1|polygon2|...
    //
    // Examples:
    //   1:2 | 2:3 # #                     // Two points (one S2PointVectorShape)
    //   # 0:0, 1:1, 2:2 | 3:3, 4:4 #      // Two polylines
    //   # # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2  // Two nested loops (one polygon)
    //   5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0    // One of each point, line, and polygon
    //   # # empty                         // One empty polygon
    //   # # empty | full                  // One empty polygon, one full polygon
    //
    // All the points, if any, are stored as a single S2PointVectorShape in the
    // index.  Polylines are stored as individual S2LaxPolylineShapes.  Polygons
    // are separated by '|', with distinct loops for a polygon separated by ';'.
    // Each group of loops is stored as an individual S2LaxPolygonShape.
    //
    // Loops should be directed so that the region's interior is on the left.
    // Loops can be degenerate (they do not need to meet S2Loop requirements).
    //
    // CAVEAT: Because whitespace is ignored, empty polygons must be specified
    //         as the string "empty" rather than as the empty string ("").
    public static MutableS2ShapeIndex MakeIndexOrDie(string str)
    {
        MyDebug.Assert(MakeIndex(str, out var index));
        return index;
    }

    // As above, but does not Assert.True-fail on invalid input. Returns true if
    // conversion is successful.
    public static bool MakeIndex(string str, [NotNullWhen(true)] out MutableS2ShapeIndex? index)
    {
        index = null;
        var result = new MutableS2ShapeIndex();
        var strs = str.Split('#');
        MyDebug.Assert(3 == strs.Length);

        var points = new List<S2Point>();
        foreach (var point_str in SplitString(strs[0], '|'))
        {
            if (!MakePoint(point_str, out var point)) return false;
            points.Add(point);
        }
        if (points.Count!=0)
        {
            result.Add(new S2PointVectorShape([.. points]));
        }
        foreach (var line_str in SplitString(strs[1], '|'))
        {
            if (!MakeLaxPolyline(line_str, out var lax_polyline)) return false;
            result.Add(lax_polyline);
        }
        foreach (var polygon_str in SplitString(strs[2], '|'))
        {
            if (!MakeLaxPolygon(polygon_str, out var lax_polygon)) return false;
            result.Add(lax_polygon);
        }
        index = result;
        return true;
    }

    private static string AppendVertex(S2LatLng ll)
    {
        var norm = ll.Normalized();
        return $"{norm.LatDegrees():g15}:{norm.LngDegrees():g15}";
    }

    private static string AppendVertex(S2Point p)
    {
        return AppendVertex(new S2LatLng(p));
    }

    private static string AppendVertices(IEnumerable<S2Point> v, int n)
    {
        return string.Join(", ", v.Take(n).Select(t => AppendVertex(t)));
    }

    // Convert an S2Point, S2LatLng, S2LatLngRect, S2CellId, S2CellUnion, loop,
    // polyline, or polygon to the string format above.
    public static string ToDebugString(this S2Point point)
    {
        return AppendVertex(point);
    }

    public static string ToDebugString(this S2LatLng latlng)
    {
        return AppendVertex(latlng);
    }

    public static string ToDebugString(this S2LatLngRect rect)
    {
        return $"{AppendVertex(rect.Lo())}, {AppendVertex(rect.Hi())}";
    }

    public static string ToDebugString(this S2CellUnion cell_union)
    {
        return string.Join(", ", cell_union.Select(t => t.ToString()));
    }

    public static string ToDebugString(this S2Loop loop)
    {
        if (loop.IsEmpty())
        {
            return "empty";
        }
        else if (loop.IsFull())
        {
            return "full";
        }

        if (loop.NumVertices > 0)
        {
            return AppendVertices(loop.CloneVertices(), loop.NumVertices);
        }
        return "";
    }

    public static string ToDebugString(this S2Polyline polyline)
    {
        if (polyline.NumVertices() > 0)
        {
            return AppendVertices(polyline.VerticesClone(), polyline.NumVertices());
        }
        return "";
    }

    public static string ToDebugString(this S2Polygon polygon, string loop_separator = ";\n")
    {
        if (polygon.IsEmpty())
        {
            return "empty";
        }
        else if (polygon.IsFull())
        {
            return "full";
        }
        var loops = polygon.Loops().Select(t => AppendVertices(t.Vertices, t.NumVertices));
        return string.Join(loop_separator, loops);
    }

    public static string ToDebugString(this S2Builder.Graph graph)
    {
        StringBuilder sb = new();
        foreach (var edge in graph.Edges)
        {
            S2Point[] vertices = [ graph.Vertex(edge.ShapeId),
                        graph.Vertex(edge.EdgeId) ];
            sb.Append(vertices.ToDebugString());
            sb.Append("; ");
        }
        return sb.ToString();
    }

    public static string ToDebugString(this S2Point[] points) => AppendVertices(points.ToArray(), points.Length);
    
    public static string ToDebugString(this S2PointLoopSpan points) => AppendVertices(points.ToArray(), points.Count);

    public static string ToDebugString(this List<S2LatLng> latlngs)
    {
        var lls = latlngs.Select(t => AppendVertex(t));
        return string.Join(", ", lls);
    }

    public static string ToDebugString(this S2LaxPolylineShape polyline)
    {
        if (polyline.NumVertices() > 0)
        {
            return AppendVertices(polyline!.Vertices()!, polyline.NumVertices());
        }
        return "";
    }

    public static string ToDebugString(this S2LaxPolygonShape polygon, string loop_separator = ";\n")
    {
        var sb = new List<string>();
        for (var i = 0; i < polygon.NumLoops; ++i)
        {
            var n = polygon.NumLoopVertices(i);
            if (n == 0)
            {
                sb.Add("full");
            }
            else
            {
                sb.Add(AppendVertices(polygon.LoopVertices(i, 0), n));
            }
        }
        return string.Join(loop_separator, sb);
    }

    // Convert the contents of an S2ShapeIndex to the format above.  The index may
    // contain S2Shapes of any type.  Shapes are reordered if necessary so that
    // all point geometry (shapes of dimension 0) are first, followed by all
    // polyline geometry, followed by all polygon geometry.
    public static string ToDebugString(this S2ShapeIndex index)
    {
        StringBuilder sb = new();
        for (var dim = 0; dim < 3; ++dim)
        {
            if (dim > 0) sb.Append('#');
            var count = 0;
            foreach (var shape in index)
            {
                if (shape is null || shape.Dimension() != dim) continue;

                sb.Append((count > 0) ? " | " : (dim > 0) ? " " : "");
                for (var i = 0; i < shape.NumChains(); ++i, ++count)
                {
                    if (i > 0) sb.Append((dim == 2) ? "; " : " | ");
                    var chain = shape.GetChain(i);
                    if (chain.Length == 0)
                    {
                        MyDebug.Assert(dim == 2);
                        sb.Append("full");
                    }
                    else
                    {
                        sb.Append(AppendVertex(shape.GetEdge(chain.Start).V0));
                    }
                    var limit = chain.Start + chain.Length;
                    if (dim != 1) --limit;
                    for (var e = chain.Start; e < limit; ++e)
                    {
                        sb.Append(", ");
                        sb.Append(AppendVertex(shape.GetEdge(e).V1));
                    }
                }
            }
            // Example output: "# #", "0:0 # #", "# # 0:0, 0:1, 1:0"
            if (dim == 1 || (dim == 0 && count > 0)) sb.Append(' ');
        }
        return sb.ToString();
    }
}
