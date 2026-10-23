//! EXPERIMENT (LOSAT_X_BNSTREAM): the subjects of a subject set searched by the threads of
//! the search pool, each subject claimed once in OID order, and the results handed over by
//! the owner thread in OID order while the threads go on with the next subjects. A thread
//! claims a subject only while it is fewer than `window` subjects past the last one handed
//! over, so that at most `window` subjects' results wait for the owner (no block ends at
//! which the threads wait for the slowest subject of the block).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427
//! ```c
//!     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//!        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//!            continue;
//!        }
//! ```
//! The hand-over keeps the order of NCBI's one-thread loop, so the hit lists are the same.
//!
//! With `assist` (the `-db` path, LOSAT_X_BNDBSTREAM, `database/x_db_stream.rs`), a thread that
//! waits (a helper for room in the window, the owner for the next result) first runs a pending
//! job of the pool (`rayon::yield_now`), so that the chunks of a long subject that another
//! thread's search hands out (a nested parallel stage) are searched by the waiting threads too, as
//! the idle workers of a `par_iter` would; it spins or sleeps only when no job is pending. The jobs
//! it can run are those nested stages' (the loops are broadcast jobs, which no other thread runs),
//! and they end on their own.
//!
//! No item is searched before every thread has started its loop (`wait_for_every_loop`): a
//! thread also runs its own broadcast job inside a wait for a nested job, so a loop started after
//! the thread took a job of an item's nested stage would run above that unfinished job.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;

/// No NCBI counterpart (scheduling only): items `0..n`, the next one to claim, the count the
/// owner has handed over, and a ring of `window` slots for the results waiting for it (item
/// `i` in slot `i % window`).
struct Window<R> {
    n: usize,
    window: usize,
    assist: bool,
    next: AtomicUsize,
    handed: AtomicUsize,
    slots: Vec<Mutex<Option<R>>>,
    ready: Vec<AtomicBool>,
    failed: AtomicBool,
}

/// No NCBI counterpart: what a claim found.
enum Claim {
    Item(usize),
    Full,
    Done,
}

// No NCBI counterpart: marks the window as failed if a search unwinds, so that the owner stops
// waiting for its result.
struct Searching<'a, R> {
    window: &'a Window<R>,
    done: bool,
}

impl<R> Drop for Searching<'_, R> {
    fn drop(&mut self) {
        if !self.done {
            self.window.failed.store(true, Ordering::Release);
        }
    }
}

// No NCBI counterpart: when the owner leaves its loop (at the end, or unwinding), the helpers
// that wait for room in the window stop waiting.
struct Closing<'a, R>(&'a Window<R>);

impl<R> Drop for Closing<'_, R> {
    fn drop(&mut self) {
        self.0.failed.store(true, Ordering::Release);
    }
}

impl<R> Window<R> {
    fn new(n: usize, window: usize, assist: bool) -> Self {
        let window = window.clamp(1, n.max(1));
        Self {
            n,
            window,
            assist,
            next: AtomicUsize::new(0),
            handed: AtomicUsize::new(0),
            slots: (0..window).map(|_| Mutex::new(None)).collect(),
            ready: (0..window).map(|_| AtomicBool::new(false)).collect(),
            failed: AtomicBool::new(false),
        }
    }

    /// The next item, if it is fewer than `window` items past the `handed` items handed over
    /// (its slot is then empty: the owner emptied it before handing over item `i - window`).
    fn claim(&self, handed: usize) -> Claim {
        let mut next = self.next.load(Ordering::Relaxed);
        loop {
            if next >= self.n {
                return Claim::Done;
            }
            if next >= handed + self.window {
                return Claim::Full;
            }
            match self.next.compare_exchange_weak(
                next,
                next + 1,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => return Claim::Item(next),
                Err(now) => next = now,
            }
        }
    }

    fn search<S>(&self, i: usize, state: &mut S, search: &impl Fn(&mut S, usize) -> R) {
        let mut guard = Searching {
            window: self,
            done: false,
        };
        let result = search(state, i);
        let slot = i % self.window;
        *self.slots[slot].lock().unwrap_or_else(|e| e.into_inner()) = Some(result);
        self.ready[slot].store(true, Ordering::Release);
        guard.done = true;
    }

    /// Helper: searches the items it claims until none is left or the owner stops.
    fn helper<S>(&self, init: &impl Fn() -> S, search: &impl Fn(&mut S, usize) -> R) {
        let mut state: Option<S> = None;
        // spin a little, then sleep 20 us, doubling up to 1 ms (a helper waits only while the
        // owner has not yet handed over the subject at the start of the window)
        let mut idle = 0u32;
        let mut nap = 20u64;
        loop {
            if self.failed.load(Ordering::Acquire) {
                return;
            }
            match self.claim(self.handed.load(Ordering::Acquire)) {
                Claim::Item(i) => {
                    idle = 0;
                    nap = 20;
                    let state = state.get_or_insert_with(init);
                    self.search(i, state, search);
                }
                Claim::Done => return,
                Claim::Full => {
                    if self.assist && run_pending_job() {
                        idle = 0;
                        nap = 20;
                        continue;
                    }
                    idle = idle.saturating_add(1);
                    if idle < 128 {
                        std::hint::spin_loop();
                    } else {
                        pause(nap);
                        nap = (nap * 2).min(1000);
                    }
                }
            }
        }
    }

    /// Owner: hands the results over in item order; while the next one is not there, it
    /// searches an item of its own if the window has room.
    fn owner<S>(
        &self,
        init: &impl Fn() -> S,
        search: &impl Fn(&mut S, usize) -> R,
        hand_over: &mut dyn FnMut(usize, R),
    ) {
        let _closing = Closing(self);
        let mut state: Option<S> = None;
        let mut idle = 0u32;
        let mut handed = 0usize;
        while handed < self.n {
            let slot = handed % self.window;
            if self.ready[slot].load(Ordering::Acquire) {
                let result = self.slots[slot]
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .take()
                    .expect("search_in_order: a ready slot holds its result");
                self.ready[slot].store(false, Ordering::Relaxed);
                self.handed.store(handed + 1, Ordering::Release);
                hand_over(handed, result);
                handed += 1;
                idle = 0;
                continue;
            }
            assert!(
                !self.failed.load(Ordering::Acquire),
                "search_in_order: a search of the pass failed"
            );
            match self.claim(handed) {
                Claim::Item(i) => {
                    idle = 0;
                    let state = state.get_or_insert_with(init);
                    self.search(i, state, search);
                }
                Claim::Full | Claim::Done => {
                    if self.assist && run_pending_job() {
                        idle = 0;
                        continue;
                    }
                    idle = idle.saturating_add(1);
                    if idle < 128 {
                        std::hint::spin_loop();
                    } else {
                        std::thread::yield_now();
                    }
                }
            }
        }
    }
}

/// No NCBI counterpart (scheduling only): with `assist`, a waiting thread runs one pending job
/// of the pool it runs in, if there is one (`true`).
fn run_pending_job() -> bool {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    return matches!(rayon::yield_now(), Some(rayon::Yield::Executed));
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    false
}

/// No NCBI counterpart (scheduling only): each thread of the pass counts itself in `started` and
/// waits until all `threads` have, so that no item is searched (and no nested job of an item can
/// be taken) before every thread has started its loop. rayon-core 1.13 runs a thread's own
/// broadcast job wherever the thread looks for work, also inside a join wait
/// (`WorkerThread::take_local_job` pops the broadcast deque after the local one, and
/// `join_context` calls it): a thread that stole a job of an item's nested parallel stage before
/// it started its loop ran the loop nested above that job, the loop waited for room in the window,
/// the window for that item, and the item for the job below it (the assist test hung for hours on a
/// loaded host). Once every thread has started, none has a loop left to start.
fn wait_for_every_loop(started: &AtomicUsize, threads: usize) {
    started.fetch_add(1, Ordering::AcqRel);
    let mut idle = 0u32;
    while started.load(Ordering::Acquire) < threads {
        idle = idle.saturating_add(1);
        if idle < 128 {
            std::hint::spin_loop();
        } else {
            std::thread::yield_now();
        }
    }
}

/// No NCBI counterpart: a helper's wait for room in the window. A native helper sleeps
/// `micros`; a wasm helper yields instead, because `std::thread::sleep` panics under the
/// browser's WASI shim (@bjorn3/browser_wasi_shim 0.4.2 returns from `poll_oneoff` without
/// writing the event count that std checks: "thread::sleep(): unexpected result of
/// poll_oneoff").
fn pause(micros: u64) {
    #[cfg(not(target_arch = "wasm32"))]
    std::thread::sleep(std::time::Duration::from_micros(micros));
    #[cfg(target_arch = "wasm32")]
    {
        let _ = micros;
        std::thread::yield_now();
    }
}

/// No NCBI counterpart (scheduling only): runs `search` for items `0..n` on every thread of
/// the pool this thread runs in (each thread with its own state from `init`, made when it
/// first claims an item), and `hand_over` for every item in item order on one of them (the
/// owner). At most `window` items past the last one handed over are claimed; with `assist`, a
/// thread that waits runs the pool's pending jobs first (module documentation). Call it from a
/// thread of the search pool (`SearchPool::install`); outside a pool, or in a pool of one
/// thread, the items are searched and handed over one at a time on this thread.
///
/// The threads run their loops as one broadcast job each, not as jobs that a thread may steal:
/// a thread that waits inside a search (a nested parallel stage) runs only stolen jobs that end
/// on their own, never a second loop that would wait for the search it interrupted. No item is
/// searched before every thread has started its loop (`wait_for_every_loop`), so that no thread
/// starts its own loop inside a wait for an item's nested job either. A pool runs one pass at a
/// time: a pass started inside an item, or beside another pass (from the items of a parallel
/// stage), could start its loops above another pass's items and wait for them there.
pub(super) fn search_in_order<S, R: Send>(
    n: usize,
    window: usize,
    assist: bool,
    init: &(impl Fn() -> S + Sync),
    search: &(impl Fn(&mut S, usize) -> R + Sync),
    hand_over: &mut (dyn FnMut(usize, R) + Send),
) {
    if n == 0 {
        return;
    }
    let shared = Window::new(n, window, assist);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if rayon::current_thread_index().is_some() && rayon::current_num_threads() > 1 {
        let shared = &shared;
        // The first thread to start its job is the owner.
        let owner = Mutex::new(Some(hand_over));
        let started = AtomicUsize::new(0);
        rayon::broadcast(|ctx| {
            let hand_over = owner.lock().unwrap_or_else(|e| e.into_inner()).take();
            wait_for_every_loop(&started, ctx.num_threads());
            match hand_over {
                Some(hand_over) => shared.owner(init, search, hand_over),
                None => shared.helper(init, search),
            }
        });
        return;
    }
    shared.owner(init, search, hand_over);
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod tests {
    use super::search_in_order;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn pool(threads: usize) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a test pool")
    }

    // Uneven work per item (a few slow items), so that results arrive out of order.
    fn work(i: usize) -> u64 {
        let rounds = if i % 37 == 5 {
            20_000
        } else {
            50 + (i * 7919) % 300
        };
        (0..rounds as u64).fold(i as u64, |acc, x| acc.wrapping_mul(31).wrapping_add(x))
    }

    #[test]
    fn hands_over_every_item_once_in_item_order() {
        for (threads, n, window) in [
            (1, 10, 4),
            (4, 0, 8),
            (4, 1, 8),
            (4, 500, 1),
            (4, 500, 16),
            (8, 1000, 64),
            (3, 7, 100),
        ] {
            for assist in [false, true] {
                let handed = pool(threads).install(|| {
                    let mut handed: Vec<(usize, u64)> = Vec::new();
                    search_in_order(
                        n,
                        window,
                        assist,
                        &|| 0usize,
                        &|calls: &mut usize, i| {
                            *calls += 1;
                            work(i)
                        },
                        &mut |i, value| handed.push((i, value)),
                    );
                    handed
                });
                let want: Vec<(usize, u64)> = (0..n).map(|i| (i, work(i))).collect();
                assert_eq!(
                    handed, want,
                    "threads {threads} n {n} window {window} assist {assist}"
                );
            }
        }
    }

    // A search with a nested parallel stage (as the speculative traceback has): a thread that
    // waits for its nested work must not pick up a second loop of this pass.
    #[test]
    fn a_nested_parallel_stage_does_not_block_the_pass() {
        use rayon::prelude::*;
        let threads = 4;
        let nested = |i: usize| -> u64 {
            (0..64usize)
                .into_par_iter()
                .map(|x| work(i + x))
                .reduce(|| 0, u64::wrapping_add)
        };
        let want: Vec<(usize, u64)> = (0..400).map(|i| (i, nested(i))).collect();
        for assist in [false, true] {
            let handed = pool(threads).install(|| {
                let mut handed: Vec<(usize, u64)> = Vec::new();
                search_in_order(
                    400,
                    8,
                    assist,
                    &|| (),
                    &|_: &mut (), i| nested(i),
                    &mut |i, value| handed.push((i, value)),
                );
                handed
            });
            assert_eq!(handed, want, "assist {assist}");
        }
    }

    // A thread of the pool that is busy when the pass starts holds back every item until it has
    // started its loop (`wait_for_every_loop`): before then it could take a nested job of an item
    // and start its loop inside the wait for that job.
    #[test]
    fn no_item_is_searched_before_every_thread_has_started_its_loop() {
        use std::sync::atomic::AtomicBool;
        use std::sync::{Arc, Mutex};
        use std::time::{Duration, Instant};
        let pool = pool(4);
        let busy = Arc::new(AtomicBool::new(false));
        let release = Arc::new(AtomicBool::new(false));
        {
            let (busy, release) = (Arc::clone(&busy), Arc::clone(&release));
            pool.spawn(move || {
                busy.store(true, Ordering::SeqCst);
                while !release.load(Ordering::SeqCst) {
                    std::hint::spin_loop();
                }
            });
        }
        while !busy.load(Ordering::SeqCst) {
            std::thread::yield_now();
        }
        let released_at: Arc<Mutex<Option<Instant>>> = Arc::new(Mutex::new(None));
        let releaser = {
            let (release, released_at) = (Arc::clone(&release), Arc::clone(&released_at));
            std::thread::spawn(move || {
                std::thread::sleep(Duration::from_millis(200));
                *released_at.lock().unwrap() = Some(Instant::now());
                release.store(true, Ordering::SeqCst);
            })
        };
        let first: Mutex<Option<Instant>> = Mutex::new(None);
        let handed = pool.install(|| {
            let mut handed = Vec::new();
            search_in_order(
                32,
                4,
                true,
                &|| (),
                &|_: &mut (), i| {
                    first.lock().unwrap().get_or_insert_with(Instant::now);
                    work(i)
                },
                &mut |i, _| handed.push(i),
            );
            handed
        });
        releaser.join().unwrap();
        assert_eq!(handed, (0..32).collect::<Vec<_>>());
        let first = first
            .into_inner()
            .unwrap()
            .expect("the items were searched");
        let released = released_at
            .lock()
            .unwrap()
            .expect("the busy thread was released");
        assert!(
            first >= released,
            "an item was searched before the busy thread started its loop"
        );
    }

    // With `assist`, the threads that wait for the window run the chunks that a slow item's
    // search hands out: a pass whose first item is long and parallel (the long subject of a `-db`
    // search) ends with the other threads' help, in item order.
    //
    // The first chunk (run inline by the thread that searches item 0, while the other chunks wait
    // in its deque) holds until another thread has run a chunk, so that the test does not depend
    // on the other threads being scheduled before item 0 ends on a loaded host (CI ran it at ~50 %
    // failures). Without the assist the other threads never take a chunk, and the deadline ends
    // the wait so that the assertion fails instead of the test hanging.
    #[test]
    fn waiting_threads_run_a_long_items_nested_jobs_with_assist() {
        use rayon::prelude::*;
        use std::collections::HashSet;
        use std::sync::Mutex;
        use std::time::{Duration, Instant};
        let threads = 4;
        let runners: Mutex<HashSet<usize>> = Mutex::new(HashSet::new());
        let search = |_: &mut (), i: usize| -> u64 {
            if i == 0 {
                (0..256usize)
                    .into_par_iter()
                    .map(|x| {
                        runners
                            .lock()
                            .unwrap()
                            .insert(rayon::current_thread_index().unwrap());
                        if x == 0 {
                            let deadline = Instant::now() + Duration::from_secs(30);
                            while runners.lock().unwrap().len() < 2 && Instant::now() < deadline {
                                std::thread::yield_now();
                            }
                        }
                        (0..20_000u64).fold(x as u64, |a, y| a.wrapping_mul(31).wrapping_add(y))
                    })
                    .reduce(|| 0, u64::wrapping_add)
            } else {
                work(i)
            }
        };
        let handed = pool(threads).install(|| {
            let mut handed: Vec<usize> = Vec::new();
            search_in_order(64, 8, true, &|| (), &search, &mut |i, _| handed.push(i));
            handed
        });
        assert_eq!(handed, (0..64).collect::<Vec<_>>());
        assert!(
            runners.lock().unwrap().len() > 1,
            "the long item ran on one thread"
        );
    }

    #[test]
    fn claims_stay_within_the_window() {
        let threads = 6;
        let window = 12;
        let handed = AtomicUsize::new(0);
        let worst = AtomicUsize::new(0);
        pool(threads).install(|| {
            search_in_order(
                2000,
                window,
                false,
                &|| (),
                &|_: &mut (), i| {
                    let ahead = i - handed.load(Ordering::SeqCst).min(i);
                    worst.fetch_max(ahead, Ordering::SeqCst);
                    work(i)
                },
                &mut |i, _| handed.store(i + 1, Ordering::SeqCst),
            );
        });
        // The owner counts an item as handed over just before `hand_over` runs.
        assert!(worst.load(Ordering::SeqCst) <= window);
    }

    #[test]
    fn makes_one_state_per_thread_at_most() {
        let threads = 4;
        let states = AtomicUsize::new(0);
        pool(threads).install(|| {
            search_in_order(
                300,
                32,
                false,
                &|| {
                    states.fetch_add(1, Ordering::SeqCst);
                },
                &|_: &mut (), i| work(i),
                &mut |_, _| {},
            );
        });
        let made = states.load(Ordering::SeqCst);
        assert!((1..=threads).contains(&made), "{made} states");
    }

    #[test]
    fn a_failed_search_stops_the_pass() {
        for failing in [0usize, 3, 299] {
            let outcome = std::panic::catch_unwind(|| {
                pool(4).install(|| {
                    search_in_order(
                        300,
                        16,
                        false,
                        &|| (),
                        &|_: &mut (), i| {
                            assert!(i != failing, "item {i} fails");
                            work(i)
                        },
                        &mut |_, _| {},
                    );
                })
            });
            assert!(outcome.is_err(), "item {failing}");
        }
    }

    #[test]
    fn runs_on_one_thread_outside_a_pool() {
        let mut handed = Vec::new();
        search_in_order(
            20,
            4,
            false,
            &|| (),
            &|_: &mut (), i| i * 2,
            &mut |i, value| handed.push((i, value)),
        );
        assert_eq!(handed, (0..20).map(|i| (i, i * 2)).collect::<Vec<_>>());
    }

    // B17 probe (run by hand with `--ignored`): many passes whose first item hands out nested jobs
    // at once, on more threads than CPUs. A pool thread that steals one of those jobs before it
    // has started its loop runs the loop nested inside the stolen job's join wait (rayon-core 1.13
    // `take_local_job` pops the thread's own broadcast deque there), above unfinished work of item
    // 0. A pass that makes no progress for LOSAT_B17_STALL_S seconds ends the process (exit 3).
    #[test]
    #[ignore]
    fn probe_b17_many_passes_with_a_nested_first_item() {
        use rayon::prelude::*;
        use std::time::{Duration, Instant};
        let var = |name: &str, default: usize| {
            std::env::var(name)
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(default)
        };
        let threads = var("LOSAT_B17_THREADS", 32);
        let passes = var("LOSAT_B17_PASSES", 100_000);
        let stall = var("LOSAT_B17_STALL_S", 20) as u64;
        let assist = var("LOSAT_B17_ASSIST", 1) == 1;
        let chunks = var("LOSAT_B17_CHUNKS", 256);
        let done = AtomicUsize::new(0);
        let start = Instant::now();
        let light =
            |x: usize| (0..200u64).fold(x as u64, |a, y| a.wrapping_mul(31).wrapping_add(y));
        std::thread::scope(|s| {
            s.spawn(|| {
                let mut last = 0;
                let mut since = Instant::now();
                loop {
                    std::thread::sleep(Duration::from_millis(200));
                    let now = done.load(Ordering::Relaxed);
                    if now >= passes {
                        return;
                    }
                    if now != last {
                        last = now;
                        since = Instant::now();
                    } else if since.elapsed() > Duration::from_secs(stall) {
                        eprintln!(
                            "B17 HANG at pass {now} of {passes} after {:.1} s (threads {threads} assist {assist})",
                            start.elapsed().as_secs_f64()
                        );
                        std::process::exit(3);
                    }
                }
            });
            pool(threads).install(|| {
                for p in 0..passes {
                    search_in_order(
                        16,
                        2,
                        assist,
                        &|| (),
                        &|_: &mut (), i: usize| -> u64 {
                            if i == 0 {
                                (0..chunks)
                                    .into_par_iter()
                                    .map(light)
                                    .reduce(|| 0, u64::wrapping_add)
                            } else {
                                light(i)
                            }
                        },
                        &mut |_, _| {},
                    );
                    done.store(p + 1, Ordering::Relaxed);
                }
            });
        });
        eprintln!(
            "B17 {passes} passes in {:.3} s (threads {threads} assist {assist})",
            start.elapsed().as_secs_f64()
        );
    }
}
