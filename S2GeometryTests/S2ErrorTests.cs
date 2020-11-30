using Xunit;

namespace S2Geometry
{
    public class S2ErrorTests
    {
        [Fact]
        public void Test_S2Error_Basic()
        {
            S2Error error = new(S2ErrorCode.DUPLICATE_VERTICES, $"Vertex {23} is the same as vertex {47}");
            // Prepend additional context to the message.
            error = new(error.Code, $"Loop {5}: {error.Text}");
            Assert.Equal(S2ErrorCode.DUPLICATE_VERTICES, error.Code);
            Assert.Equal("Loop 5: Vertex 23 is the same as vertex 47", error.Text);
        }
    }
}
