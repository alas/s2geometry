// S2MemoryTracker is a helper class for tracking and limiting the memory
// usage of S2 operations.  It provides the following functionality:
//
//  - Tracks the current and maximum memory usage of certain S2 classes
//    (including S2Builder, S2BooleanOperation, and S2BufferOperation).
//
//  - Supports cancelling the current operation if a given memory limit would
//    otherwise be exceeded.
//
//  - Invokes an optional callback after every N bytes of memory allocation,
//    and periodically within certain calculations that might take a long
//    time.  This gives the client an opportunity to cancel the current
//    operation for any reason, e.g. because the memory usage of the entire
//    thread or process is too high, a deadline was exceeded, an external
//    cancellation request was received, etc.
//
// To use it, clients simply create an S2MemoryTracker object and pass it to
// the desired S2 operations.  For example:
//
//   S2MemoryTracker tracker;
//   tracker.set_limit_bytes(500 << 20);     // 500 MB limit
//   S2Builder::Options options;
//   options.set_memory_tracker(&tracker);
//   S2Builder builder{options};
//   ...
//   S2Error error;
//   if (!builder.Build(&error)) {
//     if (error.code() == S2Error::RESOURCE_EXHAUSTED) {
//       S2_LOG(ERROR) << error;  // Memory limit exceeded
//     }
//   }
//
// Here is an example showing how to invoke a callback after every 10 MB of
// memory allocation:
//
//   tracker.set_periodic_callback(10 << 20 /*10 MB*/, [&]() {
//       if (MyCancellationCheck()) {
//         tracker.SetError(S2Error::CANCELLED, "Operation cancelled");
//       }
//     });
//
// Note that the callback is invoked based on cumulative allocation rather
// than current usage, e.g. a loop that repeatedly allocates and frees 1 MB
// would invoke the callback every 10 iterations.  Also note that the callback
// has control over what type of error is generated.
//
// This class is not thread-safe and therefore all objects associated with a
// single S2MemoryTracker should be accessed using a single thread.
//
// Implementation Notes
// --------------------
//
// In order to write a new class that tracks memory using S2MemoryTracker,
// users must analyze their data structures and make appropriate method calls.
// The major drawback to this approach is that it is fragile, since users
// might change their code but not their memory tracking.  The only way to
// avoid this problem is through rigorous testing.  See s2builder_test.cc
// and s2memory_tracker_testing.h for useful techniques.
//
// Note that malloc hooks are not a good solution for memory tracking within
// the S2 library.  The reasons for this include: (1) malloc hooks are
// program-wide and affect all threads, (2) the S2 library is used on many
// platforms (and by open source projects) and cannot depend on the features
// of specific memory allocators, and (3) certain S2 code paths can allocate a
// lot of memory at once, so it is better to predict and avoid such
// allocations rather than detecting them after the fact (as would happen with
// malloc hooks).

namespace S2Geometry;

using System.Runtime.InteropServices;

public class S2MemoryTracker
{
    // Indicates that memory usage is unlimited.
    public const long kNoLimit = long.MaxValue;

    // The current tracked memory usage.
    //
    // CAVEAT: When an operation is cancelled (e.g. due to a memory limit being
    // exceeded) the value returned may be wildly inaccurate.  This is because
    // this method reports attempted rather than actual memory allocation, and
    // S2 operations often attempt to allocate memory even on their failure /
    // early exit code paths.
    public long UsageBytes { get; private set; } = 0;

    // The maximum tracked memory usage.
    //
    // CAVEAT: When an operation is cancelled the return value may be wildly
    // inaccurate (see usage() for details).
    public long MaxUsageBytes  { get; private set; } = 0;

    // Specifies a memory limit in bytes.  Whenever the tracked memory usage
    // would exceed this value, an error of type S2Error::RESOURCE_EXHAUSTED is
    // generated and the current operation will be cancelled.  If the value is
    // kNoLimit then memory usage is tracked but not limited.
    //
    // DEFAULT: kNoLimit
    public long LimitBytes { get; set; } = kNoLimit;

    private long AllocBytes = 0;

    // Returns true if no memory tracking errors have occurred.  If this method
    // returns false then the current S2 operation will be cancelled.
    public bool IsOk() { return Error.IsOk(); }

    // Returns the tracker's current error status.  Whenever an error exists
    // the current S2 operation will be cancelled.
    //
    // Sets the error status of the memory tracker.  Typically this method is
    // called from the periodic callback (see below).  Setting the error code to
    // anything other than S2Error::OK requests cancellation of the current S2
    // operation.
    //
    // CAVEAT: Do not use this method to clear an existing error unless you know
    // what you're doing.  Clients are not required to track memory accurately
    // once an operation has been cancelled, and therefore the only safe way to
    // reset the error status is to delete the S2MemoryTracker::Client object
    // where the error occurred (which will free all its memory and restore the
    // S2MemoryTracker to an accurate state).
    public S2Error Error { get; set; }

    // A function that is called periodically to check whether the current
    // S2 operation should be cancelled.
    //public Action PeriodicCallback;

    // Sets a function that is called after every "callbackAllocDeltaBytes"
    // of tracked memory allocation to check whether the current operation
    // should be cancelled.  The callback may also be called periodically during
    // calculations that take a long time.  Once an error has occurred, further
    // callbacks are suppressed.
    public void SetPeriodicCallback(long callbackAllocDeltaBytes, Action periodicCallback)
    {
        CallbackAllocDeltaBytes = callbackAllocDeltaBytes;
        PeriodicCallback = periodicCallback;
        CallbackAllocLimitBytes = AllocBytes + CallbackAllocDeltaBytes;
    }
    public long CallbackAllocDeltaBytes { get; private set; } = 0;
    public Action? PeriodicCallback { get; private set; }
    private long CallbackAllocLimitBytes = kNoLimit;

    // Resets usage() and max_usage() to zero and clears any error.  Leaves all
    // other parameters unchanged.
    public void Reset()
    {
        Error = S2Error.OK;
        UsageBytes = MaxUsageBytes = AllocBytes = 0;
        CallbackAllocLimitBytes = CallbackAllocDeltaBytes;
    }

    //////////////////////////////////////////////////////////////////////
    //
    // Everything below this point is only needed by classes that use
    // S2MemoryTracker to track their memory usage.

    // S2MemoryTracker::Client is used to track the memory used by a given S2
    // operation.  It points to an S2MemoryTracker tracker object and updates
    // the tracker's state as memory is allocated and freed.  Several client
    // objects can point to the same memory tracker; this makes it easier for
    // one S2 operation to use another S2 operation in its implementation.  (For
    // example, S2BooleanOperation is implemented using S2Builder.)
    //
    // The client object keeps track of its own memory usage in addition to
    // updating the shared S2MemoryTracker state.  This allows the memory used
    // by a particular client to be automatically subtracted from the total when
    // that client is destroyed.
    public class Client : IDisposable
    {
        // Specifies the S2MemoryTracker that will be used to track the memory
        // usage of this client.  Several S2 operations can use the same memory
        // tracker by creating different Client objects, e.g. S2BooleanOperation
        // has a client to track its memory usage, and it also uses S2Builder
        // which creates its own client.  This allows the total memory usage of
        // both classes to be tracked and controlled.

        // Default constructor.  Note that this class can be used without calling
        // Init(), but in that case memory usage is not tracked.
        public Client() { }

        // Convenience constructor that calls Init().
        public Client(S2MemoryTracker? tracker)
        {
            Init(tracker);
        }

        // Initializes this client to use the given memory tracker.  This function
        // may be called more than once (which is equivalent to destroying this
        // client and transferring the current memory usage to a new client).
        public void Init(S2MemoryTracker? tracker)
        {
            long usage_bytes = ClientUsageBytes;
            Tally(-usage_bytes);
            Tracker = tracker;
            Tally(usage_bytes);
        }

        // Returns the memory tracker associated with this client object.
        public S2MemoryTracker? Tracker { get; private set; } = null;

        // Returns true if this client has been initialized.
        public bool IsActive() { return Tracker is not null; }

        // When a Client object is destroyed, any remaining memory is subtracted
        // from the associated S2MemoryTracker (under the assumption that the
        // associated S2 operation has been destroyed as well).
        public void Dispose()
        {
            Tally(-ClientUsageBytes);
            GC.SuppressFinalize(this);
        }

        // Returns the current tracked memory usage.
        // XXX(ericv): Return 0 when not active.
        public long UsageBytes() => Tracker?.UsageBytes ?? 0;

        // Returns the current tracked memory usage of this client only.
        // Returns zero if tracker() has not been initialized.
        public long ClientUsageBytes { get; private set; } = 0;

        // Returns the tracker's current error status.
        public S2Error Error() => Tracker?.Error ?? S2Error.OK;

        // Returns true if no memory tracking errors have occurred.  If this method
        // returns false then the current S2 operation will be cancelled.
        public bool Ok() => Tracker?.IsOk() ?? true;

        // Records a "delta" bytes of memory use (positive or negative), returning
        // false if the current operation should be cancelled.
        public bool Tally(long delta_bytes)
        {
            if (!IsActive()) return true;
            ClientUsageBytes += delta_bytes;
            return Tracker!.Tally(delta_bytes);
        }

        // Specifies that "delta" bytes of memory will be allocated and then later
        // freed.  Returns false if the current operation should be cancelled.
        public bool TallyTemp(long delta_bytes)
        {
            Tally(delta_bytes);
            return Tally(-delta_bytes);
        }

        // Adds the memory used by the given vector to the current tally.  Returns
        // false if the current operation should be cancelled.
        public /*inline*/ bool Tally<T>(List<T> v)
        {
            return Tally(v.Capacity * SizeHelper.SizeOf(typeof(T)));
        }

        // Subtracts the memory used by the given vector from the current tally.
        // Returns false if the current operation should be cancelled.
        public /*inline*/ bool Untally<T>(List<T> v)
        {
            return Tally(-v.Capacity * SizeHelper.SizeOf(typeof(T)));
        }

        // Ensures that the given vector has space for "n" additional elements and
        // tracks any additional memory that was allocated for this purpose.  This
        // method should be called before actually adding any elements to the
        // vector, otherwise memory will not be tracked correctly.  Returns false
        // if the current operation should be cancelled.
        public /*inline*/ bool AddSpace<T>(List<T> v, int n)
        {
            return AddSpaceInternal<T>(v, n, false);
        }

        // Like AddSpace() except that if more memory is needed the vector is
        // sized to hold exactly "n" additional elements.  Returns false if the
        // current operation should be cancelled.
        public /*inline*/ bool AddSpaceExact<T>(List<T> v, int n)
        {
            return AddSpaceInternal(v, n, true);
        }

        public /*inline*/ bool AddSpaceInternal<T>(List<T> v, int n, bool exact)
        {
            var new_size = v.Count + n;
            var old_capacity = v.Capacity;
            if (new_size <= old_capacity) return true;
            var new_capacity = exact ? new_size : Math.Max(new_size, 2 * old_capacity);
            // Note that reserve() allocates new storage before freeing the old storage.
            if (!Tally(new_capacity * SizeHelper.SizeOf(typeof(T)))) return false;
            v.Capacity = new_capacity;
            MyDebug.Assert(v.Capacity == new_capacity);
            return Tally(-old_capacity * SizeHelper.SizeOf(typeof(T)));
        }

        // Deallocates storage for the given vector and updates the memory
        // tracking accordingly.  Returns false if the current operation should be
        // cancelled.
        public /*inline*/ bool Clear<T>(List<T> v)
        {
            var old_capacity = v.Capacity;
            v.Clear();
            return Tally(-old_capacity * SizeHelper.SizeOf(typeof(T)));
        }

        // Returns the number of allocated bytes used by gtl::compact_array<T>.
        // (This class is similar to std::vector but stores a small number of
        // elements inline.)
        public static int GetCompactArrayAllocBytes<T>(List<T> array)
        {
            // Unfortunately this information isn't part of the public API.
            var kMaxInlinedBytes = 11;
            var kInlined = kMaxInlinedBytes / Marshal.SizeOf(typeof(T));
            int n = array.Capacity;
            return (n <= kInlined) ? 0 : n * Marshal.SizeOf(typeof(T));
        }

        // Returns the estimated minimum number of allocated bytes for each
        // additional entry in an absl::btree_* container (e.g. absl::btree_map)
        // including overhead.  This estimate is accurate only if the container
        // nodes are nearly full, as when elements are added in sorted order.
        // Otherwise nodes are expected to be 75% full on average.
        //
        // This function can be used to approximately track memory usage while
        // such a container is being built.  Once construction is complete, the
        // exact memory usage can be determined using "container.bytes_used()".
        public static int GetBtreeMinBytesPerEntry<T>()
        {
            return (int)(1.12 * SizeHelper.SizeOf(typeof(T)));
        }
    }

    private bool Tally(long delta_bytes)
    {
        UsageBytes += delta_bytes;
        AllocBytes += Math.Max(0L, delta_bytes);
        MaxUsageBytes = Math.Max(MaxUsageBytes, UsageBytes);
        if (UsageBytes > LimitBytes && IsOk()) SetLimitExceededError();
        if (PeriodicCallback is not null && AllocBytes >= CallbackAllocLimitBytes)
        {
            CallbackAllocLimitBytes = AllocBytes + CallbackAllocDeltaBytes;
            if (IsOk()) PeriodicCallback();
        }
        return IsOk();
    }

    private void SetLimitExceededError()
    {
        Error = new(S2ErrorCode.RESOURCE_EXHAUSTED,
            $"Memory limit exceeded (tracked usage {UsageBytes} bytes, limit {LimitBytes} bytes)");
    }
}
