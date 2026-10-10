//! EXPERIMENT (LOSAT_X_AHEAD[=window]): an ordered loop whose expensive,
//! pure step is evaluated ahead of the loop by the other threads of the
//! search pool.
//!
//! The pattern in BLAST's gapped stages is
//!
//! ```text
//! for i in 0..n {                      // score order, must stay ordered
//!     if tree.contains(item[i]) { continue; }
//!     let value = align(item[i]);      // expensive, depends on item[i] only
//!     tree.add(value); ...             // changes what later items see
//! }
//! ```
//!
//! Only the owner of the loop (the thread that runs it) looks at the tree and
//! decides what happens to each index, exactly as in the serial loop. What the
//! other threads may do is evaluate `align(item[j])` for indexes the owner has
//! not reached yet. Because that value is a function of the index alone, it
//! does not matter who evaluates it, when, or on which scratch memory: the
//! owner either takes a value that is already there or evaluates it itself,
//! and a value nobody takes is dropped. The schedule therefore changes the
//! running time only.
//!
//! Compared with evaluating a fixed batch in a fork-join region and then
//! consuming it, the helpers here never stop between batches (no wake-up per
//! batch), an expensive index does not hold up the rest of its batch, and the
//! owner never waits idle: while a value it needs is still being produced it
//! evaluates later indexes itself.

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard, OnceLock};

/// Look-ahead distance in indexes; `None` when the switch is off.
pub(crate) fn window() -> Option<usize> {
    static WINDOW: OnceLock<Option<usize>> = OnceLock::new();
    *WINDOW.get_or_init(|| {
        let raw = std::env::var("LOSAT_X_AHEAD").ok()?;
        Some(
            raw.parse::<usize>()
                .ok()
                .filter(|n| (2..=65_536).contains(n))
                .unwrap_or(16),
        )
    })
}

/// `LOSAT_X_AHEADSHADOW`: the owner evaluates every value it takes once more
/// on its own scratch and compares.
pub(crate) fn shadow() -> bool {
    static SHADOW: OnceLock<bool> = OnceLock::new();
    *SHADOW.get_or_init(|| std::env::var_os("LOSAT_X_AHEADSHADOW").is_some())
}

enum State<R> {
    /// Not offered to the helpers, or withdrawn by the owner.
    Closed,
    /// Offered; nobody has started.
    Open,
    /// Being evaluated.
    Running,
    /// Evaluated, not taken yet.
    Ready(R),
}

struct Cell<R> {
    index: usize,
    state: State<R>,
}

#[repr(align(128))]
struct Padded<T>(T);

pub(crate) struct Ahead<R> {
    /// Ring indexed by `index % window`; the tag says which index a cell is for.
    cells: Box<[Mutex<Cell<R>>]>,
    /// Indexes below this have been offered or closed by the owner.
    announced: Padded<AtomicUsize>,
    /// The next index a helper looks at.
    cursor: Padded<AtomicUsize>,
    finished: AtomicBool,
    evaluated: AtomicU64,
    taken: AtomicU64,
}

/// Puts an index back to `Open` if its evaluation unwinds, so that the owner
/// evaluates it itself (and meets the same panic) instead of waiting forever.
struct Reopen<'a, R> {
    cell: &'a Mutex<Cell<R>>,
    index: usize,
    armed: bool,
}

impl<R> Drop for Reopen<'_, R> {
    fn drop(&mut self) {
        if self.armed {
            let mut cell = lock(self.cell);
            if cell.index == self.index && matches!(cell.state, State::Running) {
                cell.state = State::Open;
            }
        }
    }
}

fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    // No lock is held while user code runs, so poisoning carries no meaning.
    mutex
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

impl<R> Ahead<R> {
    pub(crate) fn new(window: usize) -> Self {
        let window = window.max(2);
        Self {
            cells: (0..window)
                .map(|_| {
                    Mutex::new(Cell {
                        index: usize::MAX,
                        state: State::Closed,
                    })
                })
                .collect(),
            announced: Padded(AtomicUsize::new(0)),
            cursor: Padded(AtomicUsize::new(0)),
            finished: AtomicBool::new(false),
            evaluated: AtomicU64::new(0),
            taken: AtomicU64::new(0),
        }
    }

    fn cell(&self, index: usize) -> &Mutex<Cell<R>> {
        &self.cells[index % self.cells.len()]
    }

    /// Owner, at index `at`: offer the indexes of `[at, at + window)` below
    /// `len` that have not been offered yet. `wanted(j)` is the owner's guess,
    /// with what it knows now, whether index `j` will need its value.
    pub(crate) fn announce(&self, at: usize, len: usize, mut wanted: impl FnMut(usize) -> bool) {
        let end = at.saturating_add(self.cells.len()).min(len);
        let mut next = self.announced.0.load(Ordering::Relaxed);
        while next < end {
            let state = if wanted(next) {
                State::Open
            } else {
                State::Closed
            };
            // The cell's previous index is below `at`: the owner is done with it.
            *lock(self.cell(next)) = Cell { index: next, state };
            next += 1;
            self.announced.0.store(next, Ordering::Release);
        }
    }

    fn claim(&self) -> Option<usize> {
        let mut next = self.cursor.0.load(Ordering::Relaxed);
        loop {
            if next >= self.announced.0.load(Ordering::Acquire) {
                return None;
            }
            match self.cursor.0.compare_exchange_weak(
                next,
                next + 1,
                Ordering::AcqRel,
                Ordering::Relaxed,
            ) {
                Ok(_) => return Some(next),
                Err(now) => next = now,
            }
        }
    }

    /// Evaluates `index` if it is still on offer.
    fn evaluate<S>(&self, index: usize, scratch: &mut S, compute: &impl Fn(usize, &mut S) -> R) {
        let cell = self.cell(index);
        {
            let mut cell = lock(cell);
            if cell.index != index || !matches!(cell.state, State::Open) {
                return;
            }
            cell.state = State::Running;
        }
        let mut reopen = Reopen {
            cell,
            index,
            armed: true,
        };
        let value = compute(index, scratch);
        reopen.armed = false;
        self.evaluated.fetch_add(1, Ordering::Relaxed);
        let mut cell = lock(cell);
        // The owner may have moved on and reused the cell for a later index.
        if cell.index == index && matches!(cell.state, State::Running) {
            cell.state = State::Ready(value);
        }
    }

    /// Helper thread: evaluate offered indexes until the owner is done.
    pub(crate) fn work<S>(&self, scratch: &mut S, compute: &impl Fn(usize, &mut S) -> R) {
        let mut idle = 0u32;
        while !self.finished.load(Ordering::Acquire) {
            match self.claim() {
                Some(index) => {
                    idle = 0;
                    self.evaluate(index, scratch, compute);
                }
                None => {
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

    /// Owner: index `index` does not need its value.
    pub(crate) fn skip(&self, index: usize) {
        let mut cell = lock(self.cell(index));
        if cell.index == index && matches!(cell.state, State::Open) {
            cell.state = State::Closed;
        }
    }

    /// Owner: the value of `index` if a helper has produced it or is producing
    /// it. `None` means nobody has started; the owner then evaluates it itself.
    pub(crate) fn take<S>(
        &self,
        index: usize,
        scratch: &mut S,
        compute: &impl Fn(usize, &mut S) -> R,
    ) -> Option<R> {
        let cell = self.cell(index);
        let mut idle = 0u32;
        loop {
            {
                let mut cell = lock(cell);
                if cell.index != index {
                    return None;
                }
                match std::mem::replace(&mut cell.state, State::Closed) {
                    State::Ready(value) => {
                        self.taken.fetch_add(1, Ordering::Relaxed);
                        return Some(value);
                    }
                    State::Running => cell.state = State::Running,
                    State::Open | State::Closed => return None,
                }
            }
            // A helper is producing this value: do useful work meanwhile.
            match self.claim() {
                Some(other) => {
                    idle = 0;
                    self.evaluate(other, scratch, compute);
                }
                None => {
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

    fn finish(&self) {
        self.finished.store(true, Ordering::Release);
    }

    /// (values evaluated ahead of the owner, values the owner took)
    pub(crate) fn counts(&self) -> (u64, u64) {
        (
            self.evaluated.load(Ordering::Relaxed),
            self.taken.load(Ordering::Relaxed),
        )
    }
}

struct Finish<'a, R>(&'a Ahead<R>);

impl<R> Drop for Finish<'_, R> {
    fn drop(&mut self) {
        self.0.finish();
    }
}

/// Runs `body` (the ordered loop) on this thread while one helper per scratch
/// runs on the other threads of the pool this thread belongs to. Without a
/// pool, or with `ahead == None`, this is just `body()`.
pub(crate) fn with_helpers<R, S, T>(
    ahead: Option<&Ahead<R>>,
    scratches: impl IntoIterator<Item = S>,
    compute: &(impl Fn(usize, &mut S) -> R + Sync),
    body: impl FnOnce() -> T,
) -> T
where
    R: Send,
    S: Send,
{
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if let Some(ahead) = ahead {
        // Outside a pool `in_place_scope` would start Rayon's global pool.
        if rayon::current_thread_index().is_some() {
            return rayon::in_place_scope(|scope| {
                // Set when the loop is over, also when it unwinds: the scope
                // waits for the helpers.
                let _finish = Finish(ahead);
                for mut scratch in scratches {
                    scope.spawn(move |_| ahead.work(&mut scratch, compute));
                }
                body()
            });
        }
    }
    let _ = (ahead, scratches, compute);
    body()
}

#[cfg(test)]
mod tests {
    use super::*;

    // The owner's result must not depend on what the helpers do.
    fn value_of(index: usize) -> u64 {
        let mut x = index as u64 ^ 0x9E37_79B9_7F4A_7C15;
        for _ in 0..(index % 7) * 3000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
        }
        x
    }

    fn ordered_sum(n: usize, window: usize, helpers: usize) -> (u64, Vec<usize>) {
        let compute = |index: usize, scratch: &mut u64| -> u64 {
            *scratch += 1;
            value_of(index)
        };
        let ahead = Ahead::new(window);
        let mut kept: Vec<usize> = Vec::new();
        let mut sum = 0u64;
        let mut owner_scratch = 0u64;
        let scratches: Vec<u64> = vec![0; helpers];
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(helpers + 1)
            .build()
            .unwrap();
        pool.install(|| {
            with_helpers(Some(&ahead), scratches, &compute, || {
                for index in 0..n {
                    // What the owner knows changes with every value it keeps.
                    let live = |j: usize, sum: u64| (j as u64).wrapping_add(sum) % 3 != 0;
                    ahead.announce(index, n, |j| live(j, sum));
                    if !live(index, sum) {
                        ahead.skip(index);
                        continue;
                    }
                    let value = ahead
                        .take(index, &mut owner_scratch, &compute)
                        .unwrap_or_else(|| compute(index, &mut owner_scratch));
                    sum = sum.wrapping_mul(31).wrapping_add(value);
                    kept.push(index);
                }
            })
        });
        let (evaluated, taken) = ahead.counts();
        assert!(taken <= evaluated);
        eprintln!("window {window} helpers {helpers}: evaluated ahead {evaluated}, taken {taken}, kept {}", kept.len());
        (sum, kept)
    }

    #[test]
    fn x_ahead_owner_result_is_independent_of_helpers() {
        let serial = {
            let mut kept = Vec::new();
            let mut sum = 0u64;
            for index in 0..5000usize {
                if (index as u64).wrapping_add(sum) % 3 == 0 {
                    continue;
                }
                sum = sum.wrapping_mul(31).wrapping_add(value_of(index));
                kept.push(index);
            }
            (sum, kept)
        };
        for &(window, helpers) in &[(2, 1), (8, 1), (64, 1), (64, 3), (5, 2), (256, 7)] {
            for _ in 0..3 {
                assert_eq!(ordered_sum(5000, window, helpers), serial);
            }
        }
    }

    #[test]
    fn x_ahead_without_pool_is_the_plain_loop() {
        let compute = |index: usize, _: &mut ()| index * 2;
        let ahead = Ahead::new(8);
        let total = with_helpers(Some(&ahead), vec![(), ()], &compute, || {
            let mut total = 0;
            for index in 0..100 {
                ahead.announce(index, 100, |_| true);
                total += ahead
                    .take(index, &mut (), &compute)
                    .unwrap_or_else(|| compute(index, &mut ()));
            }
            total
        });
        assert_eq!(total, 9900);
        assert_eq!(ahead.counts(), (0, 0));
    }
}
