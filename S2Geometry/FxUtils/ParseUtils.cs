namespace S2Geometry;

public static class ParseUtils
{
    public static bool DictionaryParse(string encoded_str, out List<(string, string)> items)
    {
        items = [];
        if (string.IsNullOrEmpty(encoded_str)) return true;

        var entries = (
            from entry in encoded_str.Split(',')
            let fields = entry.Split(':')
            select (fields.Length != 2) ? (null, null) : (fields[0].Trim(), fields[1].Trim()))
            .ToList();

        if (entries.Any(t => t == (null, null))) return false;

        items = entries;
        return true;
    }
}
