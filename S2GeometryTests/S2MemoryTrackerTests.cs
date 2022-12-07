namespace S2Geometry;

public class S2MemoryTrackerTests
{
    [Fact]
    internal void Test_S2MemoryTracker_PeriodicCallback()
    {
        S2MemoryTracker tracker=new();
        int callback_count = 0;
        S2MemoryTracker.Client client=new(tracker);

        // Test that a callback interval of 0 bytes invokes the callback every time.
        tracker.SetPeriodicCallback(0, () => ++callback_count);
        client.Tally(0);
        Assert.Equal(callback_count, 1);
        client.Tally(-10);
        Assert.Equal(callback_count, 2);

        // Test that the callback interval is based on total allocated bytes rather
        // than current usage.
        tracker.SetPeriodicCallback(100, () => ++callback_count);
        client.Tally(99);
        Assert.Equal(callback_count, 2);
        client.Tally(1);
        Assert.Equal(callback_count, 3);
        client.Tally(-50);
        client.Tally(50);
        client.Tally(-50);
        client.Tally(49);
        Assert.Equal(callback_count, 3);
        client.Tally(1);
        Assert.Equal(callback_count, 4);
    }
}
