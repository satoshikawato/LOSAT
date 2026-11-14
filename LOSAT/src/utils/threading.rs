//! Search-scoped execution; NCBI owns search semantics, Rayon only schedules.
use anyhow::{ensure, Context, Result};

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
// This borrow cannot outlive the one search, or the run (`with_run_pool`), that owns the pool.
pub struct SearchPool<'a> {
    requested: &'a usize,
    #[cfg(feature = "parallel")]
    pool: Option<&'a rayon::ThreadPool>,
    // TBLASTN subject plans (was LOSAT_X_TBNSSIDE): state the owner of the search attaches for
    // the stages below it, which all receive the pool.
    // No NCBI counterpart: read-only shared state (for example a subject translation) that the
    // owner of the search attaches for the stages below. It only lets stages share a value instead
    // of computing it again; it does not change any value NCBI computes.
    x_ext: Option<&'a (dyn std::any::Any + Send + Sync)>,
}

impl SearchPool<'static> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:253-261
    // ```c
    // GetDbIndexRunSearchFn()(
    //         chunk_queries, lut_options, word_options );
    // ...
    // if (IsMultiThreaded()) {
    //      x_LaunchMultiThreadedSearch(*chunk_data);
    // } else {
    // ```
    // NCBI runs a search without the thread launch when it is not multithreaded. This pool is what a
    // search that is itself one task of a parallel loop gets: it runs on one thread, in the NCBI
    // order. Scheduling only.
    /// The pool a one-thread search gets, for a search that is itself one task of a parallel loop.
    pub fn x_serial() -> Self {
        SearchPool {
            requested: &1,
            #[cfg(feature = "parallel")]
            pool: None,
            x_ext: None,
        }
    }
}

impl<'a> SearchPool<'a> {
    // No NCBI counterpart: it returns the same pool (same thread count, same rayon pool) with a
    // reference to shared state attached. It does not change any value NCBI computes.
    /// TBLASTN subject plans (was LOSAT_X_TBNSSIDE): this pool with `ext` attached.
    pub fn x_with_ext<'b>(&'b self, ext: &'b (dyn std::any::Any + Send + Sync)) -> SearchPool<'b>
    where
        'a: 'b,
    {
        SearchPool {
            requested: self.requested,
            #[cfg(feature = "parallel")]
            pool: self.pool,
            x_ext: Some(ext),
        }
    }

    // No NCBI counterpart: reads back the shared state attached above; it does not change any value
    // NCBI computes.
    /// TBLASTN subject plans (was LOSAT_X_TBNSSIDE): the attached state, if it is a `T`.
    pub fn x_ext<T: std::any::Any>(&self) -> Option<&'a T> {
        self.x_ext.and_then(|ext| ext.downcast_ref::<T>())
    }
}

impl SearchPool<'_> {
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:147
    // TBlastThreads the_threads(GetNumberOfThreads());
    pub fn threads(&self) -> usize {
        *self.requested
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
    // (*thread)->Run(); (*thread)->Join(&result);
    pub fn enabled(&self) -> bool {
        #[cfg(feature = "parallel")]
        return self.pool.is_some();
        #[cfg(not(feature = "parallel"))]
        false
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
    // (*thread)->Run(); (*thread)->Join(&result);
    // Nested DP/linking/redo work uses this same pool, never the global pool.
    pub fn install<F, R>(&self, work: F) -> R
    where
        F: FnOnce() -> R + Send,
        R: Send,
    {
        #[cfg(feature = "parallel")]
        if let Some(pool) = self.pool {
            return pool.install(work);
        }
        work()
    }
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187,3205-3236
// arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
// NCBI warns when restricting threads. LOSAT's explicit-N contract instead
// rejects unsupported requests, without silently reducing the requested pool.
pub fn validate_threads(requested: usize) -> Result<()> {
    ensure!(requested > 0, "num_threads must be greater than zero");
    if let Some(raw) = std::env::var_os("LOSAT_WASI_THREAD_CAP") {
        let cap = raw
            .to_str()
            .context("LOSAT_WASI_THREAD_CAP must be a positive integer")?
            .parse::<usize>()
            .context("LOSAT_WASI_THREAD_CAP must be a positive integer")?;
        ensure!(cap > 0, "LOSAT_WASI_THREAD_CAP must be a positive integer");
        ensure!(
            requested <= cap,
            "requested {requested} threads exceeds LOSAT_WASI_THREAD_CAP={cap}"
        );
    }
    if requested == 1 {
        return Ok(());
    }
    ensure!(
        cfg!(all(
            feature = "parallel",
            any(
                not(target_arch = "wasm32"),
                all(losat_wasi_threads, feature = "wasm-threads")
            )
        )),
        "unsupported num_threads={requested}: this build does not support parallel search"
    );
    #[cfg(feature = "parallel")]
    ensure!(
        requested <= rayon::max_num_threads(),
        "requested {requested} threads exceeds Rayon maximum {}, which is not supported by LOSAT",
        rayon::max_num_threads()
    );
    Ok(())
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads());
// (*thread)->Run(); (*thread)->Join(&result);
// The caller occupies slot zero, so N includes the caller and N-1 children.
// Rayon-core 1.13.0 registry.rs:910-938 runs start_handler inside the worker
// loop; WorkerThread::drop at 687-693 clears its TLS registration on return.
// Unlike use_current_thread(), running that loop does not leak caller state.
pub fn with_search_pool<F, R>(requested: usize, program: &str, work: F) -> Result<R>
where
    F: FnOnce(&SearchPool<'_>) -> Result<R> + Send,
    R: Send,
{
    validate_threads(requested)?;
    #[cfg(feature = "parallel")]
    if requested > 1 {
        // The run's pool was reported once, when the run built it.
        if let Some(pool) = run_pool(requested) {
            return work(&SearchPool {
                requested: &requested,
                pool: Some(pool),
                x_ext: None,
            });
        }
        ensure!(
            rayon::current_thread_index().is_none(),
            "cannot start a search from an already registered Rayon worker"
        );
        let pool_slot = std::cell::RefCell::new(None::<rayon::ThreadPool>);
        let result_slot = std::cell::RefCell::new(None);
        let mut work = Some(work);
        let mut body = || {
            let pool = pool_slot
                .borrow_mut()
                .take()
                .expect("constructed search pool");
            // Catch before Rayon's start_handler catch boundary, which otherwise
            // aborts. Drop the pool even after an error/panic to release workers.
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                ensure!(
                    pool.current_num_threads() == requested,
                    "search pool size mismatch"
                );
                report_pool(program, requested, pool.current_num_threads());
                let search = SearchPool {
                    requested: &requested,
                    pool: Some(&pool),
                    x_ext: None,
                };
                search.install(|| work.take().expect("one search per pool")(&search))
            }));
            *result_slot.borrow_mut() = Some(result);
            drop(pool);
        };
        // SAFETY: only captured slot zero invokes this callback, synchronously
        // on this caller below. body lives until its loop and all scoped children
        // have returned; failed construction never starts slot zero.
        let callback = unsafe { CallerCallback::new(&mut body) };
        std::thread::scope(|scope| -> Result<()> {
            let mut caller = None;
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(requested)
                .start_handler(move |index| {
                    #[cfg(feature = "blastx-worker-probe")]
                    probe_lifecycle("started", index);
                    if index == 0 {
                        // SAFETY: slot zero is run only by the owning caller.
                        unsafe { callback.call() };
                    }
                })
                .exit_handler(|_index| {
                    #[cfg(feature = "blastx-worker-probe")]
                    probe_lifecycle("exited", _index);
                })
                .spawn_handler(|thread| {
                    #[cfg(feature = "blastx-worker-probe")]
                    probe_lifecycle("spawn_attempt", thread.index());
                    if thread.index() == 0 {
                        caller = Some(thread);
                    } else {
                        std::thread::Builder::new().spawn_scoped(scope, || thread.run())?;
                    }
                    Ok(())
                })
                .build()
                .with_context(|| {
                    format!("failed to build {program} pool with {requested} threads; a thread count that the system cannot start is not supported by LOSAT")
                })?;
            *pool_slot.borrow_mut() = Some(pool);
            caller.expect("caller owns pool slot zero").run();
            Ok(())
        })?;
        #[cfg(feature = "blastx-worker-probe")]
        for index in 0..requested {
            probe_lifecycle("joined", index);
        }
        // The caller loop cleared its TLS and scope joined all children before
        // returning or resuming the original user panic.
        return match result_slot.into_inner().expect("caller completed search") {
            Ok(result) => result,
            Err(panic) => std::panic::resume_unwind(panic),
        };
    }
    report_pool(program, requested, 0);
    work(&SearchPool {
        requested: &requested,
        #[cfg(feature = "parallel")]
        pool: None,
        x_ext: None,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149,242-262
// ```c
// TBlastThreads the_threads(GetNumberOfThreads());
// ...
// _TRACE("Launching BLAST with " << GetNumberOfThreads() << " threads");
// ...
//             GetDbIndexRunSearchFn()(
//                     chunk_queries, lut_options, word_options );
//
//             if (IsMultiThreaded()) {
//                  x_LaunchMultiThreadedSearch(*chunk_data);
// ```
// NCBI launches its search threads for each query chunk (lines 253-257) and joins them before
// the next chunk. Here one pool is built for the whole loop of query batches instead of one per
// batch. Each batch is still searched by the same function with the same inputs, and the results
// are used in the input order, so only the placement of work on threads changes.
/// `with_search_pool` for work that is not `Send` (it borrows the caller's output stream, for
/// instance), so that one pool can span a whole loop of query batches instead of being rebuilt,
/// threads and all, for every batch.
pub fn x_with_search_pool_local<F, R>(requested: usize, program: &str, work: F) -> Result<R>
where
    F: FnOnce(&SearchPool<'_>) -> Result<R>,
{
    struct Local<T>(T);
    // SAFETY: `with_search_pool` calls `work` exactly once, synchronously, on
    // the thread that called it: without a pool directly, with a pool as slot
    // zero, which is the caller's own thread running the Rayon worker loop
    // (`caller.run()` above), from where `ThreadPool::install` runs its
    // closure in place, or inside a run (`with_run_pool`) on the run's pool
    // from this same thread. The result is written to and read from a slot on that
    // same thread's stack. Neither value is ever touched by another thread.
    unsafe impl<T> Send for Local<T> {}
    let work = Local(work);
    let result = with_search_pool(requested, program, move |pool| {
        // Use the wrapper as a whole so that the closure captures it, not its field.
        let work: Local<F> = work;
        Ok(Local((work.0)(pool)))
    })?;
    result.0
}

// No NCBI counterpart: the pool of a run that searches more than once. The Web adapter's run
// searches once for each hit list size of its formats (docs/web/abi_v2.md §7), and BLASTP's
// block search can search again when it does not apply; each search builds its own pool in the
// CLI. In the browser a thread worker is held for each thread that may start, and the threads of
// one pool may not have returned when the next pool spawns (they need the lock that wasi-libc's
// spawn holds), so every pool after the first needed thread workers of its own.
/// Runs `work`, which may start several searches of `requested` threads, with one pool for all of
/// them: inside `work`, `with_search_pool(requested, ..)` on this thread searches on this pool
/// instead of building one. Each search is the same function on the same inputs, so only the
/// threads that run its work change, as for the batches of `x_with_search_pool_local`. When the
/// pool cannot be built (a thread count that `validate_threads` rejects, or threads that do not
/// start), `work` runs without it, and its first search fails as it does without a run.
pub fn with_run_pool<F, R>(requested: usize, program: &str, work: F) -> R
where
    F: FnOnce() -> R,
{
    let mut work = Some(work);
    #[cfg(feature = "parallel")]
    if requested > 1 {
        let pooled = x_with_search_pool_local(requested, program, |search| {
            let work = work.take().expect("one run per pool");
            let Some(pool) = search.pool else {
                return Ok(work());
            };
            struct Restore(Option<(usize, *const rayon::ThreadPool)>);
            impl Drop for Restore {
                fn drop(&mut self) {
                    RUN_POOL.with(|slot| slot.set(self.0));
                }
            }
            let _restore =
                Restore(RUN_POOL.with(|slot| {
                    slot.replace(Some((requested, pool as *const rayon::ThreadPool)))
                }));
            Ok(work())
        });
        // `work` is left only when the pool was not built.
        if let Ok(done) = pooled {
            return done;
        }
    }
    let _ = (requested, program);
    work.take().expect("work not run")()
}

#[cfg(feature = "parallel")]
thread_local! {
    /// The pool of the `with_run_pool` running on this thread, with its size.
    static RUN_POOL: std::cell::Cell<Option<(usize, *const rayon::ThreadPool)>> =
        const { std::cell::Cell::new(None) };
}

/// The pool of the `with_run_pool` running on this thread, if it has `requested` threads.
#[cfg(feature = "parallel")]
fn run_pool<'a>(requested: usize) -> Option<&'a rayon::ThreadPool> {
    let (size, pool) = RUN_POOL.with(std::cell::Cell::get)?;
    // SAFETY: `with_run_pool` sets the slot on this thread only while its `work` runs inside
    // the pool, which `with_search_pool` owns until `work` has returned, and restores it before
    // then. A search that finds the pool here runs and returns within that `work`, and its
    // `SearchPool` does not outlive the search.
    (size == requested).then(|| unsafe { &*pool })
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:172-180
// (*thread)->Run(); (*thread)->Join(&result);
// Erase only the scoped caller callback's lifetime for Rayon's 'static startup
// hook. No borrowed data is dereferenced by a child, and scope joins all children.
#[cfg(feature = "parallel")]
struct CallerCallback {
    data: *mut (),
    invoke: unsafe fn(*mut ()),
}

// SAFETY: the hook checks slot zero; it is captured instead of spawned and runs
// on the borrowing caller. Other threads only retain the pointer without use.
#[cfg(feature = "parallel")]
unsafe impl Send for CallerCallback {}
// SAFETY: the callback is invoked exactly once, by that same caller.
#[cfg(feature = "parallel")]
unsafe impl Sync for CallerCallback {}

#[cfg(feature = "parallel")]
impl CallerCallback {
    // NCBI reference: prelim_stage.cpp:172-180; Run precedes Join.
    // SAFETY: body must remain live until slot zero and all children return.
    unsafe fn new<F: FnMut()>(body: &mut F) -> Self {
        unsafe fn invoke<F: FnMut()>(data: *mut ()) {
            // SAFETY: new stored this exact F, kept live by the enclosing scope.
            unsafe { (&mut *data.cast::<F>())() };
        }
        Self {
            data: std::ptr::from_mut(body).cast(),
            invoke: invoke::<F>,
        }
    }

    // NCBI reference: prelim_stage.cpp:172-180; Run precedes Join.
    // SAFETY: invoke once, on the caller which owns the still-live body.
    unsafe fn call(&self) {
        unsafe { (self.invoke)(self.data) };
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1659-1668
// /* Use a local diagnostics structure ... shared between multiple threads */
// LOSAT_WASI_THREADS_DEBUG is read once per process, as LOSAT's other switches are: the stage
// reports run once per subject chunk, and a read of the environment there cost about 340
// instructions each (3 % of a megablast search of 20,000 short subjects).
pub fn diagnostics_enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var("LOSAT_WASI_THREADS_DEBUG")
            .is_ok_and(|v| v == "1" || v.eq_ignore_ascii_case("true"))
    })
}

// Diagnostic lines go out as one write, so a line another thread or the host
// writes at the same time (stderr is unbuffered) never lands inside it.
pub(crate) fn write_diagnostic_line(line: String) {
    use std::io::Write;
    let _ = std::io::stderr().write_all(line.as_bytes());
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:147,177-188
// TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Join(&result);
fn report_pool(program: &str, requested: usize, pool: usize) {
    if diagnostics_enabled() {
        write_diagnostic_line(format!("[losat-thread-pool] program={program} requested_threads={requested} pool_threads={pool} caller_participates={} measured_activity=null\n", pool > 0));
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1659-1668
// /* Use a local diagnostics structure ... shared between multiple threads */
// Record stage selection separately from the already verified pool size.
pub fn report_stage(program: &str, stage: &str, work_items: usize, parallel: bool) {
    if diagnostics_enabled() {
        write_diagnostic_line(format!("[losat-thread-stage] program={program} stage={stage} work_items={work_items} parallel_selected={parallel} measured_activity=null\n"));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1659-1664
// ```c++
//     /* Use a local diagnostics structure, because the one passed in an input
//       argument can be shared between multiple threads, so we don't want to pass
//       it to the engine and have a lot of mutex contention. */
//     BlastDiagnostics* local_diagnostics = Blast_DiagnosticsInit();
//
//     if ((status =
// ```
// Diagnostic-only host clock/worker identity; never read by search or reporting.
#[cfg(feature = "blastx-worker-probe")]
fn probe_lifecycle(event: &str, index: usize) {
    probe_write(format!(
        "[blastx-pool-probe] event={event} index={index} thread={:?}\n",
        std::thread::current().id()
    ));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1660-1664
// ```c++
//       argument can be shared between multiple threads, so we don't want to pass
//       it to the engine and have a lot of mutex contention. */
//     BlastDiagnostics* local_diagnostics = Blast_DiagnosticsInit();
//
//     if ((status =
// ```
// Diagnostic-only host recording uses a separate append file. Native reporting
// owns stderr's lock throughout search; workers must never acquire that lock.
#[cfg(feature = "blastx-worker-probe")]
pub(crate) fn probe_write(message: String) {
    use std::io::Write;
    let path = std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
        .expect("diagnostic probe requires LOSAT_BLASTX_WORKER_LOG");
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .expect("diagnostic worker log must be writable");
    file.write_all(message.as_bytes())
        .expect("diagnostic worker log write");
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187
    // arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
    #[test]
    fn invalid_request_never_enters_search() {
        assert!(with_search_pool(0, "test", |_| -> Result<()> {
            panic!("invalid request entered search")
        })
        .is_err());
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Run(); (*thread)->Join(&result);
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    #[test]
    fn repeated_pools_include_caller_and_preserve_errors() {
        let caller = std::thread::current().id();
        for n in [1, 2, 4, 2, 1] {
            let error = with_search_pool(n, "test", |search| -> Result<()> {
                assert_eq!(search.threads(), n);
                assert_eq!(std::thread::current().id(), caller);
                if let Some(pool) = search.pool {
                    let workers = pool.broadcast(|_| std::thread::current().id());
                    assert_eq!(workers.len(), n);
                    assert_eq!(workers.iter().filter(|id| **id == caller).count(), 1);
                }
                anyhow::bail!("original search error")
            })
            .unwrap_err();
            assert_eq!(format!("{error:#}"), "original search error");
            assert_eq!(rayon::current_thread_index(), None);
        }
    }

    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:172-180
    // (*thread)->Run(); (*thread)->Join(&result);
    #[test]
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    fn caller_registration_is_cleared_after_panic() {
        let failed = std::panic::catch_unwind(|| {
            let _ = with_search_pool(2, "test", |_| -> Result<()> { panic!("search panic") });
        });
        assert!(failed.is_err());
        assert_eq!(rayon::current_thread_index(), None);
        with_search_pool(4, "test", |_| Ok(())).unwrap();
        assert_eq!(rayon::current_thread_index(), None);
    }

    // The searches of a run share its pool: the same worker threads, no pool of their own; a
    // search of another size cannot start inside it, and after the run a search builds its pool.
    #[test]
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    fn searches_of_a_run_share_its_pool() {
        let workers = |search: &SearchPool<'_>| {
            let mut ids = search
                .pool
                .expect("a pool")
                .broadcast(|_| format!("{:?}", std::thread::current().id()));
            ids.sort();
            ids
        };
        let (first, second) = with_run_pool(3, "test", || {
            let first = with_search_pool(3, "test", |search| Ok(workers(search))).unwrap();
            let second = with_search_pool(3, "test", |search| Ok(workers(search))).unwrap();
            assert!(with_search_pool(2, "test", |_| Ok(())).is_err());
            assert!(!with_search_pool(1, "test", |search| Ok(search.enabled())).unwrap());
            (first, second)
        });
        assert_eq!(first.len(), 3);
        assert_eq!(first, second);
        // A search after the run builds its pool: new threads besides the caller.
        let caller = format!("{:?}", std::thread::current().id());
        let after = with_search_pool(3, "test", |search| Ok(workers(search))).unwrap();
        assert!(after
            .iter()
            .filter(|id| **id != caller)
            .all(|id| !first.contains(id)));
        // A run of one thread has no pool, and its searches none either.
        assert!(
            !with_run_pool(1, "test", || with_search_pool(1, "test", |search| Ok(
                search.enabled()
            )))
            .unwrap()
        );
        // A run whose pool cannot be built runs without it: its search fails as on its own.
        let too_many = rayon::max_num_threads() + 1;
        let mut ran = false;
        let error = with_run_pool(too_many, "test", || {
            ran = true;
            with_search_pool(too_many, "test", |_| Ok(()))
        })
        .unwrap_err();
        assert!(ran);
        assert_eq!(
            format!("{error:#}"),
            format!(
                "{:#}",
                with_search_pool(too_many, "test", |_| Ok(())).unwrap_err()
            )
        );
    }
}
