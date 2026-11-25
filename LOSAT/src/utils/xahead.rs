//! No NCBI counterpart: an ordered loop whose expensive,
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
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919,4045-4046,4087-4088
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375,404-405,509-512,600-601
//! These are the two ordered loops of NCBI (`BLAST_GetGappedScore`, the preliminary
//! gapped extension, and `Blast_TracebackFromHSPList`). Each iteration tests the tree,
//! aligns, and adds to the tree. This module does not change that order. The owner thread
//! still tests and adds in the NCBI order; the other threads only evaluate the
//! alignment call for indexes the owner has not reached, which depends on the index
//! alone. A value nobody takes is dropped. This changes which thread computes a value and
//! when, and nothing else.

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard};

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919,4045-4046
/// The number of upcoming indexes of the loop above that may be evaluated ahead (the built-in
/// value of the former `LOSAT_X_AHEAD`, which graduated on 2026-10-22). Placement of work only.
pub(crate) const WINDOW: usize = 16;

// No NCBI counterpart: the states of one slot of the look-ahead ring; the slot holds a
// value that NCBI computes inside its loop, computed earlier.
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

// No NCBI counterpart: cache-line padding against false sharing; it does not change any
// value NCBI computes.
#[repr(align(128))]
struct Padded<T>(T);

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919,4045-4046
// One slot per upcoming index of this loop (and of the loop of `Blast_TracebackFromHSPList`).
// A slot holds the value of the alignment call for that index. The loop itself and its tree
// operations stay on the owner thread.
pub(crate) struct Ahead<R> {
    /// Ring indexed by `index % window`; the tag says which index a cell is for.
    cells: Box<[Mutex<Cell<R>>]>,
    /// Indexes below this have been offered or closed by the owner.
    announced: Padded<AtomicUsize>,
    /// The loop's length as the owner last gave it (`announce`), `usize::MAX` before.
    total: Padded<AtomicUsize>,
    /// Whether a helper leaves once every index of the loop is claimed (`x_new_leaving`), instead
    /// of waiting for the owner to finish.
    leave: bool,
    /// Whether a wait with nothing to do sleeps after spinning and yielding (`x_new_napping`).
    naps: bool,
    /// The next index a helper looks at.
    cursor: Padded<AtomicUsize>,
    finished: AtomicBool,
    evaluated: AtomicU64,
    taken: AtomicU64,
}

/// Puts an index back to `Open` if its evaluation unwinds, so that the owner
/// evaluates it itself (and meets the same panic) instead of waiting forever.
// No NCBI counterpart: puts a slot back on offer if its evaluation panics, so that the
// owner evaluates it itself; it does not change any value NCBI computes.
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

// No NCBI counterpart: a mutex lock that ignores poisoning; it does not change any value
// NCBI computes.
fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    // No lock is held while user code runs, so poisoning carries no meaning.
    mutex
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

impl<R> Ahead<R> {
    // No NCBI counterpart: makes an empty ring of `window` slots; it does not change any value
    // NCBI computes.
    pub(crate) fn new(window: usize) -> Self {
        Self::x_with(window, false, false)
    }

    // No NCBI counterpart: scheduling only. A helper that leaves frees its thread for other work
    // of the pool (another loop's tasks) while the owner takes the last values.
    /// `new(window)` whose helpers leave once every index of the loop is claimed.
    #[cfg_attr(
        not(feature = "parallel"),
        expect(dead_code, reason = "used only with the parallel feature")
    )]
    pub(crate) fn x_new_leaving(window: usize) -> Self {
        Self::x_with(window, true, false)
    }

    // No NCBI counterpart: scheduling only. For a loop whose values take milliseconds (BLASTX's
    // composition redo, `blastx/runtime/kappa_batch.rs`): a helper or the owner that waits sleeps
    // 50 us at a time once it has spun and yielded (`x_idle`), so that a long wait does not keep
    // a core in `sched_yield` (6 to 13 s of system time in 13 to 18 s of wall time, q200k
    // `-use_sw_tback` at BATCH_SIZE 40000 and 32 threads), while a value is ready again within
    // 50 us, little beside one value's time.
    /// `x_new_leaving(window)` whose waits nap.
    #[cfg_attr(
        not(feature = "parallel"),
        expect(dead_code, reason = "used only with the parallel feature")
    )]
    pub(crate) fn x_new_napping(window: usize) -> Self {
        Self::x_with(window, true, true)
    }

    // No NCBI counterpart: the constructor of the three (scheduling only).
    fn x_with(window: usize, leave: bool, naps: bool) -> Self {
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
            total: Padded(AtomicUsize::new(usize::MAX)),
            leave,
            naps,
            cursor: Padded(AtomicUsize::new(0)),
            finished: AtomicBool::new(false),
            evaluated: AtomicU64::new(0),
            taken: AtomicU64::new(0),
        }
    }

    fn cell(&self, index: usize) -> &Mutex<Cell<R>> {
        &self.cells[index % self.cells.len()]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3878
    /// The owner offers the indexes after the one it is at, in the order of the NCBI loop.
    /// `wanted` is the owner's guess whether an index will reach the alignment call; a wrong
    /// guess costs time only, because the owner decides for each index itself.
    /// Owner, at index `at`: offer the indexes of `[at, at + window)` below
    /// `len` that have not been offered yet. `wanted(j)` is the owner's guess,
    /// with what it knows now, whether index `j` will need its value.
    pub(crate) fn announce(&self, at: usize, len: usize, mut wanted: impl FnMut(usize) -> bool) {
        self.x_announce_while(at, len, |next| Some(wanted(next)));
    }

    // No NCBI counterpart (scheduling only): `announce` of a loop of `len` indexes that offers none
    // at or past `until` yet (the flow across BLASTX chunks offers the tasks of the loop's chunk and
    // the next one only, `blastx/x_subject_flow.rs` `run_chunks`). Helpers that leave do so once
    // every index below `len` is claimed, as with `announce`.
    #[cfg_attr(
        not(any(feature = "parallel", test)),
        expect(dead_code, reason = "used only with the parallel feature or tests")
    )]
    pub(crate) fn x_announce_until(
        &self,
        at: usize,
        until: usize,
        len: usize,
        mut wanted: impl FnMut(usize) -> bool,
    ) {
        self.x_announce_while(at, len, |next| (next < until).then(|| wanted(next)));
    }

    // No NCBI counterpart (scheduling only): `announce` whose guess may also say "not yet"
    // (`guess(j)` = `None`): the offers stop at that index, which with the later ones stays
    // unoffered until a later call, made when the owner knows more, offers or closes it (the
    // composition redo's matches redone ahead, `blastx/runtime/kappa_batch.rs`, which stops where
    // only its matches still on offer would end the redo early). The loop's length stays `len`, so
    // helpers that leave wait while the offers are stopped. `announce` and `x_announce_until` are
    // this loop with a guess that says "not yet" never or from `until` on.
    pub(crate) fn x_announce_while(
        &self,
        at: usize,
        len: usize,
        mut guess: impl FnMut(usize) -> Option<bool>,
    ) {
        self.total.0.store(len, Ordering::Release);
        let end = at.saturating_add(self.cells.len()).min(len);
        let mut next = self.announced.0.load(Ordering::Relaxed);
        while next < end {
            let state = match guess(next) {
                Some(true) => State::Open,
                Some(false) => State::Closed,
                None => return,
            };
            // The cell's previous index is below `at`: the owner is done with it.
            *lock(self.cell(next)) = Cell { index: next, state };
            next += 1;
            self.announced.0.store(next, Ordering::Release);
        }
    }

    // No NCBI counterpart: a helper takes the next offered index; which thread evaluates an
    // index does not change its value.
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

    /// No NCBI counterpart: runs the alignment call of the NCBI loop for an index ahead of the
    /// owner. The call depends on the index alone, so the value is the one the owner would get.
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

    /// No NCBI counterpart: the helper loop (scheduling of independent work); it does not
    /// change any value NCBI computes.
    /// Helper thread: evaluate offered indexes until the owner is done.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        expect(
            dead_code,
            reason = "used only where the helper pool runs (parallel, and wasm-threads on wasm32)"
        )
    )]
    pub(crate) fn work<S>(&self, scratch: &mut S, compute: &impl Fn(usize, &mut S) -> R) {
        let mut idle = 0u32;
        while !self.finished.load(Ordering::Acquire) {
            match self.claim() {
                Some(index) => {
                    idle = 0;
                    self.evaluate(index, scratch, compute);
                }
                None => {
                    // Every index below `total` is claimed: no later one will be offered.
                    if self.leave
                        && self.cursor.0.load(Ordering::Acquire)
                            >= self.total.0.load(Ordering::Acquire)
                    {
                        return;
                    }
                    self.x_idle(&mut idle);
                }
            }
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3918-3919
    /// The tree test said "contained", so NCBI does not call the aligner for this index. The
    /// owner withdraws the index from the helpers.
    /// Owner: index `index` does not need its value.
    pub(crate) fn skip(&self, index: usize) {
        let mut cell = lock(self.cell(index));
        if cell.index == index && matches!(cell.state, State::Open) {
            cell.state = State::Closed;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3918-3921,4045-4046
    /// Called by the owner at the point where NCBI calls the aligner, after the tree test has
    /// said "not contained". The owner uses a helper's value for this index, or computes it
    /// itself.
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
                None => self.x_idle(&mut idle),
            }
        }
    }

    // No NCBI counterpart (scheduling only): one step of a wait with nothing to claim, of a helper
    // (`work`) or of the owner (`take`): spin, then yield; with `naps` (`x_new_napping`), after 64
    // yields, sleep 50 us at a time. A wasm thread keeps yielding: `std::thread::sleep` panics
    // under the browser's WASI shim (`blastn/blast_engine/x_bnstream.rs`, `pause`).
    fn x_idle(&self, idle: &mut u32) {
        *idle = idle.saturating_add(1);
        if *idle < 128 {
            std::hint::spin_loop();
        } else if !self.naps || *idle < 128 + 64 || cfg!(target_arch = "wasm32") {
            std::thread::yield_now();
        } else {
            std::thread::sleep(std::time::Duration::from_micros(50));
        }
    }

    // No NCBI counterpart: tells the helpers that the loop has ended; it does not change any
    // value NCBI computes.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn finish(&self) {
        self.finished.store(true, Ordering::Release);
    }

    /// No NCBI counterpart: counters for LOSAT_X_STATS; they do not change any value NCBI
    /// computes.
    /// (values evaluated ahead of the owner, values the owner took)
    pub(crate) fn counts(&self) -> (u64, u64) {
        (
            self.evaluated.load(Ordering::Relaxed),
            self.taken.load(Ordering::Relaxed),
        )
    }
}

// No NCBI counterpart: ends the helpers when the loop exits, also on a panic.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Finish<'a, R>(&'a Ahead<R>);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<R> Drop for Finish<'_, R> {
    fn drop(&mut self) {
        self.0.finish();
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850,3918-3919,4087-4088
/// `body` is this loop, run in the NCBI order on the calling thread. The other threads of
/// the pool only evaluate values ahead; they never touch the tree or the result list.
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
#[path = "xahead_tests.rs"]
mod tests;
