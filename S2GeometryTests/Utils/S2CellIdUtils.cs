namespace S2Geometry;

internal class S2CellIdUtils
{
    public static S2CellId FromDebugString(string str)
    {
        // This function is reasonably efficient, but is only intended for use in
        // tests.
        int level = (int)str.Length - 2;
        if (level < 0 || level > S2.kMaxCellLevel) return S2CellId.None;
        int face = str[0] - '0';
        if (face < 0 || face > 5 || str[1] != '/') return S2CellId.None;
        S2CellId id = S2CellId.FromFace(face);
        for (int i = 2; i < str.Length; ++i)
        {
            int child_pos = str[i] - '0';
            if (child_pos < 0 || child_pos > 3) return S2CellId.None;
            id = id.Child(child_pos);
        }
        return id;
    }
}
