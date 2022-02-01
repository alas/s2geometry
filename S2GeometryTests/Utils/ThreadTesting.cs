// A test harness for verifying the thread safety properties of a class.  It
// creates one writer thread and several reader threads, and repeatedly
// alternates between executing some operation in the writer thread and then
// executing some other operation in parallel in the reader threads.
//
// It is intended for testing thread-compatible classes, i.e. those where
// const methods are thread safe and non-const methods are not thread safe.

namespace S2GeometryTests.Utils;

public abstract class ReaderWriterTest
{
    // The main loop of each reader thread.
    public void ReaderLoop()
    {
        lock_.WaitOne();
        for (int last_write = 0; ; last_write = num_writes_)
        {
            while (num_writes_ == last_write)
            {
                lock (lock_)
                {
                    Monitor.Wait(write_ready_);
                }
            }
            if (num_writes_ < 0) break;

            // Release the lock first so that all reader threads can run in parallel.
            lock_.ReleaseMutex();
            ReadOp();
            lock_.WaitOne();
            if (--num_readers_left_ == 0)
            {
                Monitor.Pulse(all_readers_done_);
            }
        }
        lock_.ReleaseMutex();
    }

    // Create the given number of reader threads and execute the given number of
    // (write, read) iterations.
    public void Run(int num_readers, int iters)
    {
        using ReaderThreadPool pool = new(ReaderLoop, num_readers);
        lock_.WaitOne();
        for (int iter = 0; iter < iters; ++iter)
        {
            // Loop invariant: lock_ is held and num_readers_left_ == 0.
            Assert.Equal(0, num_readers_left_);
            WriteOp();

            // Now set the readers loose.
            num_readers_left_ = num_readers;
            ++num_writes_;
            Monitor.PulseAll(write_ready_);
            while (num_readers_left_ > 0)
            {
                lock (lock_)
                {
                    Monitor.Wait(all_readers_done_);
                }
            }
        }
        // Signal the readers to exit.
        num_writes_ = -1;
        Monitor.PulseAll(write_ready_);
        lock_.ReleaseMutex();
        // ReaderThreadPool destructor waits for all threads to complete.
    }

    // The writer thread calls the following function once between reads.
    public abstract void WriteOp();

    // Each reader thread calls the following function once between writes.
    public abstract void ReadOp();

    // The following fields are guarded by lock_.
    private readonly Mutex lock_ = new();
    private int num_writes_ = 0;
    private int num_readers_left_ = 0;

    // TODO(ericv): Consider removing these condition variables now that the

    // Signalled when a new write is ready to be processed.
    private readonly object write_ready_ = new();
    // Signalled when all readers have processed the latest write.
    private readonly object all_readers_done_ = new();
}

internal class ReaderThreadPool : IDisposable
{
    public ReaderThreadPool(Action callback, int num_threads)
    {
        threads_ = new Thread[num_threads];
        num_threads_ = num_threads;

        for (int i = 0; i < num_threads_; ++i)
        {
            threads_[i] = new Thread(new ThreadStart(callback));
        }
    }

    public void Dispose()
    {
        for (int i = 0; i < num_threads_; ++i) threads_[i].Join();
    }

    private readonly Thread[] threads_;
    private readonly int num_threads_;
};
