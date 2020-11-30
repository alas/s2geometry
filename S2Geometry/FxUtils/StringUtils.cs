using System.Text;

namespace S2Geometry
{
    public static class StringUtils
    {
        // HexEncodeStr returns the data in str in hex encoded form.
        public static string ToHexa(this string str)
        {
            var sb = new StringBuilder(str.Length * 2);

            foreach (var b in System.Text.Encoding.ASCII.GetBytes(str))
            {
                sb.AppendFormat($"{b:X}");
            }
            return sb.ToString();
        }
        // HexEncodeStr returns the data in str in hex encoded form.
        public static string ToHexa(this byte[] ba) => ToHexa(ba, ba.Length);
        public static string ToHexa(this byte[] ba, int count)
        {
            var sb = new StringBuilder(ba.Length * 2);

            for (var i = 0; i < count; i++)
            {
                sb.AppendFormat($"{ba[i]:X2}");
            }
            return sb.ToString();
        }

        public static string NoCR(this string s) => s.Replace("\r", "");
    }
}
