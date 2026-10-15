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

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;

/// No NCBI counterpart (scheduling only): items `0..n`, the next one to claim, the count the
/// owner has handed over, and a ring of `window` slots for the results waiting for it (item
/// `i` in slot `i % window`).
struct Window<R> {
    n: usize,
    window: usize,
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
    fn new(n: usize, window: usize) -> Self {
        let window = window.clamp(1, n.max(1));
        Self {
            n,
            window,
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
                    idle = idle.saturating_add(1);
                    if idle < 128 {
                        std::hint::spin_loop();
                    } else {
                        std::thread::sleep(std::time::Duration::from_micros(nap));
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
                    .expect("LOSAT_X_BNSTREAM: a ready slot holds its result");
                self.ready[slot].store(false, Ordering::Relaxed);
                self.handed.store(handed + 1, Ordering::Release);
                hand_over(handed, result);
                handed += 1;
                idle = 0;
                continue;
            }
            assert!(
                !self.failed.load(Ordering::Acquire),
                "LOSAT_X_BNSTREAM: a subject search failed"
            );
            match self.claim(handed) {
                Claim::Item(i) => {
                    idle = 0;
                    let state = state.get_or_insert_with(init);
                    self.search(i, state, search);
                }
                Claim::Full | Claim::Done => {
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

/// No NCBI counterpart (scheduling only): runs `search` for items `0..n` on every thread of
/// the pool this thread runs in (each thread with its own state from `init`, made when it
/// first claims an item), and `hand_over` for every item in item order on one of them (the
/// owner). At most `window` items past the last one handed over are claimed. Call it from a
/// thread of the search pool (`SearchPool::install`); outside a pool, or in a pool of one
/// thread, the items are searched and handed over one at a time on this thread.
///
/// The threads run their loops as one broadcast job each, not as jobs that a thread may steal:
/// a thread that waits inside a search (a nested parallel stage) runs only stolen jobs that end
/// on their own, never a second loop that would wait for the search it interrupted.
pub(super) fn search_in_order<S, R: Send>(
    n: usize,
    window: usize,
    init: &(impl Fn() -> S + Sync),
    search: &(impl Fn(&mut S, usize) -> R + Sync),
    hand_over: &mut (dyn FnMut(usize, R) + Send),
) {
    if n == 0 {
        return;
    }
    let shared = Window::new(n, window);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if rayon::current_thread_index().is_some() && rayon::current_num_threads() > 1 {
        let shared = &shared;
        // The first thread to start its job is the owner.
        let owner = Mutex::new(Some(hand_over));
        rayon::broadcast(|_| {
            let hand_over = owner.lock().unwrap_or_else(|e| e.into_inner()).take();
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
            let handed = pool(threads).install(|| {
                let mut handed: Vec<(usize, u64)> = Vec::new();
                search_in_order(
                    n,
                    window,
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
            assert_eq!(handed, want, "threads {threads} n {n} window {window}");
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
        let handed = pool(threads).install(|| {
            let mut handed: Vec<(usize, u64)> = Vec::new();
            search_in_order(
                400,
                8,
                &|| (),
                &|_: &mut (), i| nested(i),
                &mut |i, value| handed.push((i, value)),
            );
            handed
        });
        let want: Vec<(usize, u64)> = (0..400).map(|i| (i, nested(i))).collect();
        assert_eq!(handed, want);
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
        search_in_order(20, 4, &|| (), &|_: &mut (), i| i * 2, &mut |i, value| {
            handed.push((i, value))
        });
        assert_eq!(handed, (0..20).map(|i| (i, i * 2)).collect::<Vec<_>>());
    }
}
