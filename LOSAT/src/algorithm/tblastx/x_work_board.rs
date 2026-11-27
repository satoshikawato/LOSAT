//! EXPERIMENT (LOSAT_X_TBXBOARD): one board per TBLASTX search on which the subjects' units and
//! the linking calls' groups are handed out, groups first, then units by batch and subject order.
//!
//! No NCBI counterpart: NCBI searches its subjects one after the other on each thread
//! (blast_engine.c:1409-1413) and links each subject's HSP list on that thread
//! (blast_engine.c:870-873); LOSAT splits both into independent items (the units of PAIRPAR, the
//! groups of LINKSLICES) whose results do not depend on the thread that runs them or on when it
//! runs them (x_pair_par.rs `run_units`, x_link_slices.rs items 2-4). This file only chooses which
//! thread runs which item, and when.
//!
//! # Why
//!
//! With rayon's work stealing every subject of a batch (TBXNEXT) runs its units at once, so the
//! subjects of a batch end their seed stage together and their linking chains all start at the
//! end of the batch; a chain's link groups, pushed on its thread's deque, are stolen only by a
//! thread that has run out of its own units. On the genome 3 x 3 input at 16 threads the search
//! ends with one thread linking a self-comparison's largest groups after every unit has run. Here:
//!
//! 1. Units are taken in batch order, then subject order: a batch's first subject ends its seed
//!    stage first, and its chain links while the later subjects' units run.
//! 2. A thread that helps takes groups before units: a chain's groups run on the threads that take
//!    items as soon as they are free, while the remaining units wait, so the chain ends early.
//! 3. A subject's thread takes its own units only. When it took groups and earlier subjects' units
//!    first, the other threads took most of its units meanwhile and it then blocked on them: with
//!    the genome 3 x 3 input's subjects cut into 18,838 pieces of 100 bp the search took longer at
//!    16 threads than without the board.
//!
//! # How (no thread waits in the pool)
//!
//! An owner (a subject's thread with its units, a chain's thread with its groups) puts an item
//! set on the board and takes its set's items only, until none is left to take, then blocks until
//! the items others took have run. A helper (a thread of the batch's subject loop without a subject)
//! takes items until there is nothing left it may take, and sleeps while a subject of its batch has
//! not yet put its units on the board; a `rayon::spawn` job of the board takes one group.
//! Items are pure work: a unit's scan and a group's kernel make no pool call and wait for nothing,
//! so an item that was taken always ends, and a blocked owner always wakes. Taking an item never
//! runs a pool job, so no subject or chain runs nested inside another (the problem of
//! x_batch_pipe.rs item 6).
//!
//! Lifetimes: an item set's closure borrows its owner's frame; the board holds it as a raw pointer
//! that is read only for an index below the set's length, after a successful take, and the owner
//! returns only when every taken index has run (`XSet::wait`).

use std::cell::Cell;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock};

// The item sets and their owner's hold (`utils/xitems.rs`), shared with BLASTX.
pub(crate) use crate::utils::xitems::{in_item, XTaker};
use crate::utils::xitems::{Run, XOwned, XSet};

// No NCBI counterpart: the variable is read once per process; it chooses which thread runs which
// unit and group, not a value NCBI computes.
/// `LOSAT_X_TBXBOARD` (strict set, `utils/xswitch.rs`): native builds with the parallel feature.
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        cfg!(all(feature = "parallel", not(target_arch = "wasm32")))
            && crate::utils::xswitch::on("LOSAT_X_TBXBOARD")
    })
}

/// Where a thread of a subject loop sleeps while no item is left to take and a subject has not yet
/// put its units on the board. A ring that comes just before the sleeper counts itself is missed,
/// so the sleep is bounded: such a thread wakes 1 ms later at most, and so does a thread whose loop
/// stopped (`XBoard::help`).
#[derive(Default)]
struct XBell {
    lock: Mutex<()>,
    rung: Condvar,
    sleepers: AtomicUsize,
}

impl XBell {
    fn ring(&self) {
        if self.sleepers.load(Ordering::SeqCst) > 0 {
            let _lock = self.lock.lock().unwrap_or_else(|e| e.into_inner());
            self.rung.notify_all();
        }
    }

    fn sleep(&self) {
        let lock = self.lock.lock().unwrap_or_else(|e| e.into_inner());
        self.sleepers.fetch_add(1, Ordering::SeqCst);
        let waited = self
            .rung
            .wait_timeout(lock, std::time::Duration::from_millis(1));
        self.sleepers.fetch_sub(1, Ordering::SeqCst);
        drop(waited);
    }
}

/// The units of one batch's subjects, one slot per subject index.
pub(crate) struct XBatch {
    bell: Arc<XBell>,
    slots: Box<[Mutex<Option<Arc<XSet>>>]>,
    /// The subject has put its units on the board, or it will put none.
    closed: Box<[AtomicBool]>,
    /// No subject below this one has units left to take (a hint).
    lowest: AtomicUsize,
    /// Every subject below this one is drained (`drained`): a slot that is closed and has no unit
    /// left to take stays so, so the scan never goes back.
    drained_below: AtomicUsize,
}

impl XBatch {
    /// Subject `s` has put its units on the board (or will put none) and none is left to take.
    fn slot_drained(&self, s: usize) -> bool {
        // `closed` before the slot: `run_units` stores its set before it closes the slot.
        self.closed[s].load(Ordering::Acquire)
            && self.slots[s]
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .as_ref()
                .is_none_or(|set| set.taken())
    }

    /// Every subject has put its units on the board (or will put none) and none is left to take.
    /// Each call reads the slots from the first one not yet known to be drained and stops at the
    /// first that is not, so the calls of a batch read each slot about once.
    pub(crate) fn drained(&self) -> bool {
        let n = self.slots.len();
        let mut s = self.drained_below.load(Ordering::Acquire);
        while s < n {
            if !self.slot_drained(s) {
                return false;
            }
            s += 1;
            self.drained_below.fetch_max(s, Ordering::AcqRel);
        }
        true
    }

    /// Takes and runs one unit of the lowest subject that has one left.
    fn run_one(&self, taker: &mut XTaker) -> bool {
        let n = self.slots.len();
        let mut s = self.lowest.load(Ordering::Acquire);
        let mut all_below_taken = true;
        while s < n {
            // `closed` before the slot (`slot_drained`).
            let closed = self.closed[s].load(Ordering::Acquire);
            let set = self.slots[s]
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .clone();
            match set {
                Some(set) if !set.taken() => {
                    if set.run_one(taker) {
                        return true;
                    }
                }
                Some(_) if all_below_taken => {
                    self.lowest.fetch_max(s + 1, Ordering::AcqRel);
                }
                Some(_) => {}
                // Closed without units, or after its units ran (`run_units` clears its slot).
                None if closed => {
                    if all_below_taken {
                        self.lowest.fetch_max(s + 1, Ordering::AcqRel);
                    }
                }
                None => {
                    all_below_taken = false;
                }
            }
            s += 1;
        }
        false
    }

    /// The subject `s` will put no units on the board (it returned before its units).
    pub(crate) fn close_without_units(&self, s: usize) {
        self.closed[s].store(true, Ordering::Release);
        self.bell.ring();
    }
}

/// The group sets of one search. The board's helper job holds it, so it outlives the board.
#[derive(Default)]
struct XGroups {
    sets: Mutex<Vec<Arc<XSet>>>,
    /// A helper job waits in a thread's queue (`x_spawn_group_helper`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    queued: AtomicBool,
}

impl XGroups {
    /// Takes and runs one item of any set; false when there was none.
    fn run_one(&self, taker: &mut XTaker) -> bool {
        loop {
            let set = {
                let mut sets = self.sets.lock().unwrap_or_else(|e| e.into_inner());
                sets.retain(|set| !set.taken());
                sets.first().cloned()
            };
            match set {
                Some(set) => {
                    if set.run_one(taker) {
                        return true;
                    }
                }
                None => return false,
            }
        }
    }

    /// An item of some set is left to take.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn any_left(&self) -> bool {
        self.sets
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .iter()
            .any(|set| !set.taken())
    }
}

/// The board of one search.
#[derive(Default)]
pub(crate) struct XBoard {
    groups: Arc<XGroups>,
    batches: Mutex<Vec<Arc<XBatch>>>,
    bell: Arc<XBell>,
}

impl XBoard {
    /// A batch of `subjects` subjects, taken after the batches before it.
    pub(crate) fn batch(&self, subjects: usize) -> Arc<XBatch> {
        let batch = Arc::new(XBatch {
            bell: Arc::clone(&self.bell),
            slots: (0..subjects).map(|_| Mutex::new(None)).collect(),
            closed: (0..subjects).map(|_| AtomicBool::new(false)).collect(),
            lowest: AtomicUsize::new(0),
            drained_below: AtomicUsize::new(0),
        });
        self.batches
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(Arc::clone(&batch));
        batch
    }

    /// The batch has ended its subject loop.
    pub(crate) fn end_batch(&self, batch: &Arc<XBatch>) {
        self.batches
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .retain(|b| !Arc::ptr_eq(b, batch));
    }

    /// Takes and runs one item: a group of any linking call first, then a unit by batch and
    /// subject order. False when there was none.
    fn run_one(&self, taker: &mut XTaker) -> bool {
        self.groups.run_one(taker) || {
            let batches: Vec<Arc<XBatch>> = self
                .batches
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .clone();
            batches.iter().any(|batch| batch.run_one(taker))
        }
    }

    /// Takes and runs one item (groups first, then units by batch and subject order); false when
    /// there was none. For a thread that waits outside the pool.
    pub(crate) fn take_one(&self, taker: &mut XTaker) -> bool {
        self.run_one(taker)
    }

    /// A subject's units: `run(i, taker)` for every `i < n`, on this thread (which takes this
    /// set's units only, Why 3.), on the board's takers and, for a batch's only subject, on up to
    /// `threads - 1` helper jobs of the set. Returns when every unit has run.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        expect(dead_code, reason = "used only with `parallel` threads")
    )]
    pub(crate) fn run_units(&self, batch: &XBatch, s: usize, n: usize, run: &Run<'_>) {
        let set = XSet::new(run, n);
        let owned = XOwned(&set);
        *batch.slots[s].lock().unwrap_or_else(|e| e.into_inner()) = Some(Arc::clone(&set));
        batch.closed[s].store(true, Ordering::Release);
        self.bell.ring();
        // The set's own helpers, for the batch's only subject: it has no subject-loop threads to
        // take its units. In a subject loop (TBXNEXT's broadcast) every pool thread is in the loop
        // and takes units from the board (its subjects', then `help`), so helpers there would only
        // wait in the threads' queues until the loop ends. A helper takes this set's units only,
        // so a thread that runs it inside a pool wait holds it for this subject's units at most.
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if batch.slots.len() == 1 {
            for _ in 0..rayon::current_num_threads().min(n).saturating_sub(1) {
                let set = Arc::clone(&set);
                rayon::spawn(move || {
                    let mut taker = XTaker::default();
                    while set.run_one(&mut taker) {}
                });
            }
        }
        let mut taker = XTaker::default();
        while set.run_one(&mut taker) {}
        drop(taker);
        set.wait();
        drop(owned);
        // The slot stays closed and drained; the set is not kept until the batch ends.
        *batch.slots[s].lock().unwrap_or_else(|e| e.into_inner()) = None;
    }

    /// A linking call's groups: `run(i, taker)` for every `i < n` (the largest group first is the
    /// caller's order), on this thread, on the board's takers before any unit, and on the board's
    /// helper jobs. Returns when every group has run.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        expect(dead_code, reason = "used only with `parallel` threads")
    )]
    pub(crate) fn run_groups(&self, n: usize, run: &Run<'_>) {
        self.run_set(n, run, true);
    }

    /// `run_groups`, with or without the helper job (a panic in a `rayon::spawn` job aborts, so the
    /// tests of a panicking item run without it).
    fn run_set(
        &self,
        n: usize,
        run: &Run<'_>,
        #[cfg_attr(
            not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
            expect(
                unused_variables,
                reason = "the helper job is spawned only with `parallel` threads"
            )
        )]
        helper: bool,
    ) {
        let set = XSet::new(run, n);
        let owned = XOwned(&set);
        self.groups
            .sets
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(Arc::clone(&set));
        self.bell.ring();
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if helper && n > 1 {
            x_spawn_group_helper(&self.groups);
        }
        let mut taker = XTaker::default();
        while set.run_one(&mut taker) {}
        set.wait();
        drop(owned);
        // Its set leaves the list here: with every thread an owner nothing else may take from the
        // list, and a batch of N subjects would keep N sets until a helper job runs.
        self.groups
            .sets
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .retain(|other| !Arc::ptr_eq(other, &set));
    }

    /// A thread of `batch`'s subject loop without a subject: takes items until the batch's units
    /// are all taken (it waits, without the pool, while a subject has not yet put its units), or
    /// the loop stopped (`stopped`: a subject's search panicked, so the subjects not yet taken will
    /// never close).
    #[cfg_attr(
        not(all(feature = "parallel", not(target_arch = "wasm32"))),
        expect(
            dead_code,
            reason = "used only by the subject loop's helper, with `parallel` off wasm32"
        )
    )]
    pub(crate) fn help(&self, batch: &XBatch, stopped: &AtomicBool) {
        let mut taker = XTaker::default();
        loop {
            if self.run_one(&mut taker) {
                continue;
            }
            if batch.drained() || stopped.load(Ordering::Relaxed) {
                return;
            }
            self.bell.sleep();
        }
    }
}

/// The board's helper job, one in the threads' queues at a time: a thread that runs it spawns the
/// next while groups are left, so idle threads join one after the other, and sets whose items
/// were all taken by others leave at most this one job behind, however many subjects a batch has.
/// One item per job: a job that a thread runs inside a pool wait (nested) holds that thread for one
/// item at most.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_spawn_group_helper(groups: &Arc<XGroups>) {
    if groups.queued.swap(true, Ordering::AcqRel) {
        return;
    }
    let groups = Arc::clone(groups);
    rayon::spawn(move || {
        groups.queued.store(false, Ordering::Release);
        if groups.any_left() {
            x_spawn_group_helper(&groups);
        }
        let mut taker = XTaker::default();
        groups.run_one(&mut taker);
    });
}

thread_local! {
    // The board and batch of the subject this thread searches (`enter`).
    static CURRENT: Cell<Option<(*const XBoard, *const XBatch)>> = const { Cell::new(None) };
}

/// The board and batch of this thread's subject until the guard drops.
pub(crate) struct XEntered(Option<(*const XBoard, *const XBatch)>);

impl Drop for XEntered {
    fn drop(&mut self) {
        CURRENT.with(|current| current.set(self.0));
    }
}

/// This thread searches a subject of `batch` on `board` until the guard drops.
pub(crate) fn enter(board: &XBoard, batch: &XBatch) -> XEntered {
    XEntered(CURRENT.with(|current| current.replace(Some((board, batch)))))
}

/// The board and batch of this thread's subject (`enter`), if any.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn current<'a>() -> Option<(&'a XBoard, &'a XBatch)> {
    // SAFETY: the pointers are set by `enter` from references that outlive its guard, and the
    // guard restores the value before them on drop; the subject's work that reads them runs
    // inside the guard's life on this thread.
    CURRENT
        .with(Cell::get)
        .map(|(board, batch)| unsafe { (&*board, &*batch) })
}

#[cfg(all(test, feature = "parallel", not(target_arch = "wasm32")))]
#[path = "x_work_board_tests.rs"]
mod tests;
