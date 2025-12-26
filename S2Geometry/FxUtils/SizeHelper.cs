namespace S2Geometry;

using System.Reflection.Emit;
using System.Runtime.InteropServices;

public static class SizeHelper
{
    public static int SizeOf<T>()
    {
        return SizeOf(typeof(T));
    }

    public static int SizeOf<T>(T _)
    {
        return SizeOf<T>();
    }

    public static int SizeOf(this Type type)
    {
        var dynamicMethod = new DynamicMethod("SizeOf", typeof(int), Type.EmptyTypes);
        var generator = dynamicMethod.GetILGenerator();

        generator.Emit(OpCodes.Sizeof, type);
        generator.Emit(OpCodes.Ret);

        var function = (Func<int>)dynamicMethod.CreateDelegate(typeof(Func<int>));
        return function();
    }

    public static LayoutKind AlignOf(this Type type)
    {
        // The common language runtime uses the Auto layout value by default.
        // To reduce layout-related problems associated with the Auto value,
        // C#, Visual Basic, and C++ compilers specify Sequential layout for value types.

        var sla = type.StructLayoutAttribute;
        return sla is not null
            ? sla.Value
            : type.IsValueType
                ? LayoutKind.Sequential
                : LayoutKind.Auto;
    }
}
