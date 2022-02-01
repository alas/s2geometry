namespace S2Geometry;

public class S2ErrorTests
{
    [Fact]
    public void Test_S2Error_Basic()
    {
        S2Error error = new(S2ErrorCode.DUPLICATE_VERTICES,
            $"Vertex {23} is the same as vertex {47}");
        // Prepend additional context to the message.
        error = new(error.Code, $"Loop {5}: {error.Text}");
        Assert.Equal(S2ErrorCode.DUPLICATE_VERTICES, error.Code);
        Assert.Equal("Loop 5: Vertex 23 is the same as vertex 47", error.Text);
    }

    [Fact]
    public void Test_S2Error_Constructor()
    {
        S2Error error = new(S2ErrorCode.RESOURCE_EXHAUSTED,
            $"Memory limit exceeded ({100} vs {50})");
        Assert.Equal(S2ErrorCode.RESOURCE_EXHAUSTED, error.Code);
        Assert.Equal("Memory limit exceeded (100 vs 50)", error.Text);
    }
}
