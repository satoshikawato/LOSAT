//! EXPERIMENT (LOSAT_X_BNSTREAM): the subjects of a subject set searched by the threads of
//! the search pool, each subject claimed once in OID order, and the results handed over by
//! the owner thread in OID order while the threads go on with the next subjects. A thread
//! claims a subject only while it is fewer than `window` subjects past the last one handed
//! over, so that at most `window` subjects' results wait for the owner (no block ends at
//! which the threads wait for the slowest subject of the block).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427
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

/// No NCBI counterpart (scheduling only): what narrows a window while many bytes of results wait
/// for the owner. Past `floor` items after the last one handed over, an item is claimed only
/// while `waiting_bytes` (kept by the caller) is at most `max_bytes`. The item at the start of the
/// window is always within `floor`, so the owner always gets to hand it over.
pub(super) struct WaitingBound<'a> {
    pub(super) floor: usize,
    pub(super) waiting_bytes: &'a AtomicUsize,
    pub(super) max_bytes: usize,
}

/// No NCBI counterpart (scheduling only): items `0..n`, the next one to claim, the count the
/// owner has handed over, and a ring of `window` slots for the results waiting for it (item
/// `i` in slot `i % window`).
struct Window<'a, R> {
    n: usize,
    window: usize,
    assist: bool,
    bound: Option<WaitingBound<'a>>,
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
struct Searching<'a, 'b, R> {
    window: &'a Window<'b, R>,
    done: bool,
}

impl<R> Drop for Searching<'_, '_, R> {
    fn drop(&mut self) {
        if !self.done {
            self.window.failed.store(true, Ordering::Release);
        }
    }
}

// No NCBI counterpart: when the owner leaves its loop (at the end, or unwinding), the helpers
// that wait for room in the window stop waiting.
struct Closing<'a, 'b, R>(&'a Window<'b, R>);

impl<R> Drop for Closing<'_, '_, R> {
    fn drop(&mut self) {
        self.0.failed.store(true, Ordering::Release);
    }
}

impl<'a, R> Window<'a, R> {
    fn new(n: usize, window: usize, assist: bool, bound: Option<WaitingBound<'a>>) -> Self {
        let window = window.clamp(1, n.max(1));
        Self {
            n,
            window,
            assist,
            bound,
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
            if let Some(bound) = &self.bound {
                if next >= handed + bound.floor
                    && bound.waiting_bytes.load(Ordering::Relaxed) > bound.max_bytes
                {
                    return Claim::Full;
                }
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
/// starts its own loop inside a wait for an item's nested job either. Callers must run one pass
/// at a time per pool (the code does not check it): a pass started inside an item, or beside
/// another pass (from the items of a parallel stage), could start its loops above the other
/// pass's items and wait there, at its start or for its window (a pass inside an item of a pass
/// without `assist` never fills its start, because the outer loops run no pending job).
pub(super) fn search_in_order<S, R: Send>(
    n: usize,
    window: usize,
    assist: bool,
    init: &(impl Fn() -> S + Sync),
    search: &(impl Fn(&mut S, usize) -> R + Sync),
    hand_over: &mut (dyn FnMut(usize, R) + Send),
) {
    search_in_order_bounded(n, window, assist, None, init, search, hand_over);
}

/// No NCBI counterpart (scheduling only): `search_in_order` with a window that `bound` narrows
/// (`WaitingBound`).
pub(super) fn search_in_order_bounded<S, R: Send>(
    n: usize,
    window: usize,
    assist: bool,
    bound: Option<WaitingBound<'_>>,
    init: &(impl Fn() -> S + Sync),
    search: &(impl Fn(&mut S, usize) -> R + Sync),
    hand_over: &mut (dyn FnMut(usize, R) + Send),
) {
    if n == 0 {
        return;
    }
    let shared = Window::new(n, window, assist, bound);
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
#[path = "x_bnstream_tests.rs"]
mod tests;
