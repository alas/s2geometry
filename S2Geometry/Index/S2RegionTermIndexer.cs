// S2RegionTermIndexer is a helper class for adding spatial data to an
// information retrieval system.  Such systems work by converting documents
// into a collection of "index terms" (e.g., representing words or phrases),
// and then building an "inverted index" that maps each term to a list of
// documents (and document positions) where that term occurs.
//
// This class deals with the problem of converting spatial data into index
// terms, which can then be indexed along with the other document information.
//
// Spatial data is represented using the S2Region type.  Useful S2Region
// subtypes include:
//
//   S2Cap
//    - a disc-shaped region
//
//   S2LatLngRect
//    - a rectangle in latitude-longitude coordinates
//
//   S2Polyline
//    - a polyline
//
//   S2Polygon
//    - a polygon, possibly with multiple holes and/or shells
//
//   S2CellUnion
//    - a region approximated as a collection of S2CellIds
//
//   S2ShapeIndexRegion
//    - an arbitrary collection of points, polylines, and polygons
//
//   S2ShapeIndexBufferedRegion
//    - like the above, but expanded by a given radius
//
//   S2RegionUnion, S2RegionIntersection
//    - the union or intersection of arbitrary other regions
//
// So for example, if you want to query documents that are within 500 meters
// of a polyline, you could use an S2ShapeIndexBufferedRegion containing the
// polyline with a radius of 500 meters.
//
// Example usage:
//
//   // This class is intended to be used with an external key-value store,
//   // but for this example will we use an unordered_map.  The key is an
//   // index term, and the value is a set of document ids.
//   unordered_map<string, int[]> index;
//
//   // Create an indexer that uses up to 10 cells to approximate each region.
//   S2RegionTermIndexer.Options options;
//   options.set_max_cells(10);
//   S2RegionTermIndexer indexer(options);
//
//   // For this example, we index a disc-shaped region with a 10km radius.
//   S2LatLng center = S2LatLng.FromDegrees(44.1, -56.235);
//   S1Angle radius = S2Earth.ToAngle(util.units.Kilometers(10.0));
//   S2Cap cap(center.ToPoint(), radius);
//
//   // Add the terms for this disc-shaped region to the index.
//   foreach (var term in indexer.GetIndexTerms(cap)) {
//     index[term].Add(kSomeDocumentId);
//   }
//
//   // And now at query time: build a latitude-longitude rectangle.
//   S2LatLngRect rect(S2LatLng.FromDegrees(-12.1, 10.2),
//                     S2LatLng.FromDegrees(-9.2,  120.5));
//
//   // Convert the query region to a set of terms, and compute the union
//   // of the document ids associated with those terms.
//   set<int> doc_ids;
//   foreach (var term in indexer.GetQueryTerms(rect)) {
//     doc_ids.insert(index[term].begin(), index[term].end());
//   }
//
//   // "doc_ids" now contains all documents that intersect the query region,
//   // along with some documents that nearly intersect it.  The results can
//   // be further pruned if desired by retrieving the original regions that
//   // were indexed (i.e., the document contents) and checking for exact
//   // intersection with the query region.
//
// Indexing Strategy
// -----------------
//
// Given a query region, we want to find all of the document regions that
// intersect it.  The first step is to represent all the regions as S2Cell
// coverings (see S2RegionCoverer).  We then split the problem into two parts,
// namely finding the document regions that are "smaller" than the query
// region and those that are "larger" than the query region.
//
// We do this by defining two terms for each S2CellId: a "covering term" and
// an "ancestor term".  (In the implementation below, covering terms are
// distinguished by prefixing a '$' to them.)  For each document region, we
// insert a covering term for every cell in the region's covering, and we
// insert an ancestor term for these cells *and* all of their ancestors.
//
// Then given a query region, we can look up all the document regions that
// intersect its covering by querying the union of the following terms:
//
// 1. An "ancestor term" for each cell in the query region.  These terms
//    ensure that we find all document regions that are "smaller" than the
//    query region, i.e. where the query region contains a cell that is either
//    a cell of a document region or one of its ancestors.
//
// 2. A "covering term" for every ancestor of the cells in the query region.
//    These terms ensure that we find all the document regions that are
//    "larger" than the query region, i.e. where document region contains a
//    cell that is a (proper) ancestor of a cell in the query region.
//
// Together, these terms find all of the document regions that intersect the
// query region.  Furthermore, the number of terms to be indexed and queried
// are both fairly small, and can be bounded in terms of max_cells() and the
// number of cell levels used.
//
// Optimizations
// -------------
//
// + Cells at the maximum level being indexed (max_level()) have the special
//   property that they will never be an ancestor of a cell in the query
//   region.  Therefore we can safely skip generating "covering terms" for
//   these cells (see query step 2 above).
//
// + If the index will contain only points (rather than general regions), then
//   we can skip all the covering terms mentioned above because there will
//   never be any document regions larger than the query region.  This can
//   significantly reduce the size of queries.
//
// + If it is more important to optimize index size rather than query speed,
//   the number of index terms can be reduced by creating ancestor terms only
//   for the *proper* ancestors of the cells in a document region, and
//   compensating for this by including covering terms for all cells in the
//   query region (in addition to their ancestors).
//
//   Effectively, when the query region and a document region contain exactly
//   the same cell, we have a choice about whether to treat this match as a
//   "covering term" or an "ancestor term".  One choice minimizes query size
//   while the other minimizes index size.

namespace S2Geometry;

public class S2RegionTermIndexer
{
    // The following parameters control the tradeoffs between index size, query
    // size, and accuracy (see s2region_coverer.h for details).
    //
    // IMPORTANT: You must use the same values for min_level(), max_level(), and
    // level_mod() for both indexing and queries, otherwise queries will return
    // incorrect results.  However, max_cells() can be changed as often as
    // desired -- you can even change this parameter for every region.
    public class Options : S2RegionCoverer.Options
    {
        public Options()
        {
            // Override the S2RegionCoverer defaults.
            MaxCells = 8;
            MinLevel = 4;
            MaxLevel = 16;
            LevelMod = 1;
        }

        ///////////////// Options Inherited From S2RegionCoverer ////////////////

        // max_cells() controls the maximum number of cells when approximating
        // each region.  This parameter value may be changed as often as desired
        // (using Options(), see below), e.g. to approximate some regions
        // more accurately than others.
        //
        // Increasing this value during indexing will make indexes more accurate
        // but larger.  Increasing this value for queries will make queries more
        // accurate but slower.  (See s2region_coverer.h for details on how this
        // parameter affects accuracy.)  For example, if you don't mind large
        // indexes but want fast serving, it might be reasonable to set
        // max_cells() == 100 during indexing and max_cells() == 8 for queries.
        //
        // DEFAULT: 8  (coarse approximations)
        // using S2RegionCoverer.Options.max_cells;
        // using S2RegionCoverer.Options.set_max_cells;

        // min_level() and max_level() control the minimum and maximum size of the
        // S2Cells used to approximate regions.  Setting these parameters
        // appropriately can reduce the size of the index and speed up queries by
        // reducing the number of terms needed.  For example, if you know that
        // your query regions will rarely be less than 100 meters in width, then
        // you could set max_level() as follows:
        //
        //   options.set_max_level(S2.kAvgEdge.GetClosestLevel(
        //       S2Earth.MetersToRadians(100)));
        //
        // This restricts the index to S2Cells that are approximately 100 meters
        // across or larger.  Similar, if you know that query regions will rarely
        // be larger than 1000km across, then you could set min_level() similarly.
        //
        // If min_level() is set too high, then large regions may generate too
        // many query terms.  If max_level() is set too low, then small query
        // regions will not be able to discriminate which regions they intersect
        // very precisely and may return many more candidates than necessary.
        //
        // If you have no idea about the scale of the regions being queried,
        // it is perfectly fine to set min_level() == 0 and max_level() == 30
        // (== S2.kMaxLevel).  The only drawback is that may result in a larger
        // index and slower queries.
        //
        // The default parameter values are suitable for query regions ranging
        // from about 100 meters to 3000 km across.
        //
        // DEFAULT: 4  (average cell width == 600km)
        // using S2RegionCoverer.Options.min_level;
        // using S2RegionCoverer.Options.set_min_level;

        // DEFAULT: 16 (average cell width == 150m)
        // using S2RegionCoverer.Options.max_level;
        // using S2RegionCoverer.Options.set_max_level;

        // Setting level_mod() to a value greater than 1 increases the effective
        // branching factor of the S2Cell hierarchy by skipping some levels.  For
        // example, if level_mod() == 2 then every second level is skipped (which
        // increases the effective branching factor to 16).  You might want to
        // consider doing this if your query regions are typically very small
        // (e.g., single points) and you don't mind increasing the index size
        // (since skipping levels will reduce the accuracy of cell coverings for a
        // given max_cells() limit).
        //
        // DEFAULT: 1  (don't skip any cell levels)
        // using S2RegionCoverer.Options.level_mod;
        // using S2RegionCoverer.Options.set_level_mod;

        // If your index will only contain points (rather than regions), be sure
        // to set this flag.  This will generate smaller and faster queries that
        // are specialized for the points-only case.
        //
        // With the default quality settings, this flag reduces the number of
        // query terms by about a factor of two.  (The improvement gets smaller
        // as max_cells() is increased, but there is really no reason not to use
        // this flag if your index consists entirely of points.)
        //
        // DEFAULT: false
        public bool IndexContainsPointsOnly { get; set; } = false;

        // If true, the index will be optimized for space rather than for query
        // time.  With the default quality settings, this flag reduces the number
        // of index terms and increases the number of query terms by the same
        // factor (approximately 1.3).  The factor increases up to a limiting
        // ratio of 2.0 as max_cells() is increased.
        //
        // CAVEAT: This option has no effect if the index contains only points.
        //
        // DEFAULT: false
        public bool OptimizeForSpace { get; set; } = false;

        // A non-alphanumeric character that is used internally to distinguish
        // between two different types of terms (by adding this character).
        //
        // REQUIRES: "ch" is non-alphanumeric.
        // DEFAULT: '$'
        public string Marker() => marker_;
        public char MarkerCharacter
        {
            get => marker_[0];
            set
            {
                System.Diagnostics.Debug.Assert(!char.IsLetterOrDigit(value));
                marker_ = new string(value, 1);
            }
        }

        private string marker_ = new('$', 1);
    }

    // Constructs an S2RegionTermIndexer with the given options.
    public S2RegionTermIndexer(Options options) => Options_ = options;

    // Returns the current options.  Options can be modifed between calls.
    public Options Options_ { get; set; }

    // Converts the given region into a set of terms for indexing.  Terms
    // consist of lowercase letters, numbers, '$', and an optional prefix.
    //
    // "prefix" is a unique prefix used to distinguish S2 terms from other terms
    // in the repository.  The prefix may also be used to index documents with
    // multiple types of location information (e.g. store footprint, entrances,
    // parking lots, etc).  The prefix should be kept short since it is
    // prepended to every term.
    public List<string> GetIndexTerms(IS2Region region, string prefix)
    {
        // Note that options may have changed since the last call.
        coverer_.Options_ = Options_;
        var covering = coverer_.GetCovering(region);
        return GetIndexTermsForCanonicalCovering(covering, prefix);
    }

    // Converts a given query region into a set of terms.  If you compute the
    // union of all the documents associated with these terms, the result will
    // include all documents whose index region intersects the query region.
    //
    // "prefix" should match the corresponding value used when indexing.
    public List<string> GetQueryTerms(IS2Region region, string prefix)
    {
        // Note that options may have changed since the last call.
        coverer_.Options_ = Options_;
        var covering = coverer_.GetCovering(region);
        return GetQueryTermsForCanonicalCovering(covering, prefix);
    }

    // Convenience methods that accept an S2Point rather than S2Region.  (These
    // methods are also faster.)
    //
    // Note that you can index an S2LatLng by converting it to an S2Point first:
    //     var terms = GetIndexTerms(S2Point(latlng), ...);
    public List<string> GetIndexTerms(S2Point point, string prefix)
    {
        // See the top of this file for an overview of the indexing strategy.
        //
        // The last cell generated by this loop is effectively the covering for
        // the given point.  You might expect that this cell would be indexed as a
        // covering term, but as an optimization we always index these cells as
        // ancestor terms only.  This is possible because query regions will never
        // contain a descendant of such cells.  Note that this is true even when
        // max_level() != true_max_level() (see S2RegionCoverer.Options).

        var id = new S2CellId(point);
        var terms = new List<string>();
        for (var level = Options_.MinLevel; level <= Options_.MaxLevel;
             level += Options_.LevelMod)
        {
            terms.Add(GetTerm(TermType.ANCESTOR, id.Parent(level), prefix));
        }
        return terms;
    }
    public List<string> GetQueryTerms(S2Point point, string prefix)
    {
        // See the top of this file for an overview of the indexing strategy.

        var id = new S2CellId(point);
        var terms = new List<string>();
        // Recall that all true_max_level() cells are indexed only as ancestor terms.
        var level = Options_.TrueMaxLevel;
        terms.Add(GetTerm(TermType.ANCESTOR, id.Parent(level), prefix));
        if (Options_.IndexContainsPointsOnly) return terms;

        // Add covering terms for all the ancestor cells.
        for (; level >= Options_.MinLevel; level -= Options_.LevelMod)
        {
            terms.Add(GetTerm(TermType.COVERING, id.Parent(level), prefix));
        }
        return terms;
    }

    // Low-level methods that accept an S2CellUnion covering of the region to be
    // indexed or queried.
    //
    // REQUIRES: "covering" satisfies the S2RegionCoverer.Options for this
    //           class (i.e., max_cells, min_level, max_level, and level_mod).
    //
    // If you have a covering that was computed using different options, then
    // you can either call the regular S2Region methods (since S2CellUnion is a
    // type of S2Region), or "canonicalize" the covering first by calling
    // S2RegionCoverer.CanonicalizeCovering() with the same options.
    public List<string> GetIndexTermsForCanonicalCovering(S2CellUnion covering, string prefix)
    {
        // See the top of this file for an overview of the indexing strategy.
        //
        // Cells in the covering are normally indexed as covering terms.  If we are
        // optimizing for query time rather than index space, they are also indexed
        // as ancestor terms (since this lets us reduce the number of terms in the
        // query).  Finally, as an optimization we always index true_max_level()
        // cells as ancestor cells only, since these cells have the special property
        // that query regions will never contain a descendant of these cells.

        System.Diagnostics.Debug.Assert(!Options_.IndexContainsPointsOnly);
#if s2debug
        coverer_.Options_ = Options_;
        System.Diagnostics.Debug.Assert(coverer_.IsCanonical(covering));
#endif
        var terms = new List<string>();
        var prev_id = S2CellId.None;
        var true_max_level = Options_.TrueMaxLevel;
        foreach (var id in covering)
        {
            // IsCanonical() already checks the following conditions, but we repeat
            // them here for documentation purposes.
            var level = id.Level();
            System.Diagnostics.Debug.Assert(level >= Options_.MinLevel);
            System.Diagnostics.Debug.Assert(level <= Options_.MaxLevel);
            System.Diagnostics.Debug.Assert(0 == (level - Options_.MinLevel) % Options_.LevelMod);

            if (level < true_max_level)
            {
                // Add a covering term for this cell.
                terms.Add(GetTerm(TermType.COVERING, id, prefix));
            }
            if (level == true_max_level || !Options_.OptimizeForSpace)
            {
                // Add an ancestor term for this cell at therained level.
                terms.Add(GetTerm(TermType.ANCESTOR, id.Parent(level), prefix));
            }
            // Finally, add ancestor terms for all the ancestors of this cell.
            while ((level -= Options_.LevelMod) >= Options_.MinLevel)
            {
                var ancestor_id = id.Parent(level);
                if (prev_id != S2CellId.None && prev_id.Level() > level &&
                    prev_id.Parent(level) == ancestor_id)
                {
                    break;  // We have already processed this cell and its ancestors.
                }
                terms.Add(GetTerm(TermType.ANCESTOR, ancestor_id, prefix));
            }
            prev_id = id;
        }
        return terms;
    }
    public List<string> GetQueryTermsForCanonicalCovering(S2CellUnion covering, string prefix)
    {
        // See the top of this file for an overview of the indexing strategy.

#if s2debug
        coverer_.Options_ = Options_;
        System.Diagnostics.Debug.Assert(coverer_.IsCanonical(covering));
#endif
        var terms = new List<string>();
        var prev_id = S2CellId.None;
        var true_max_level = Options_.TrueMaxLevel;
        foreach (var id in covering)
        {
            // IsCanonical() already checks the following conditions, but we repeat
            // them here for documentation purposes.
            var level = id.Level();
            System.Diagnostics.Debug.Assert(level >= Options_.MinLevel);
            System.Diagnostics.Debug.Assert(level <= Options_.MaxLevel);
            System.Diagnostics.Debug.Assert(0 == (level - Options_.MinLevel) % Options_.LevelMod);

            // Cells in the covering are always queried as ancestor terms.
            terms.Add(GetTerm(TermType.ANCESTOR, id, prefix));

            // If the index only contains points, there are no covering terms.
            if (Options_.IndexContainsPointsOnly) continue;

            // If we are optimizing for index space rather than query time, cells are
            // also queried as covering terms (except for true_max_level() cells,
            // which are indexed and queried as ancestor cells only).
            if (Options_.OptimizeForSpace && level < true_max_level)
            {
                terms.Add(GetTerm(TermType.COVERING, id, prefix));
            }
            // Finally, add covering terms for all the ancestors of this cell.
            while ((level -= Options_.LevelMod) >= Options_.MinLevel)
            {
                var ancestor_id = id.Parent(level);
                if (prev_id != S2CellId.None && prev_id.Level() > level &&
                    prev_id.Parent(level) == ancestor_id)
                {
                    break;  // We have already processed this cell and its ancestors.
                }
                terms.Add(GetTerm(TermType.COVERING, ancestor_id, prefix));
            }
            prev_id = id;
        }
        return terms;
    }

    private enum TermType { ANCESTOR, COVERING };

    // There are generally more ancestor terms than covering terms, so we add
    // the extra "marker" character to the covering terms to distinguish them.
    private string GetTerm(TermType term_type, S2CellId id, string prefix) => term_type == TermType.ANCESTOR
            ? prefix + id.ToToken()
            : prefix + Options_.Marker() + id.ToToken();

    private readonly S2RegionCoverer coverer_ = new();
}
