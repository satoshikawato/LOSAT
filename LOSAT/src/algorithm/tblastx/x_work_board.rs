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
//! 2. Groups are taken before units: a chain's groups run on the threads that take items as soon as
//!    they are free, while the remaining units wait, so the chain ends early.
//!
//! # How (no thread waits in the pool)
//!
//! An owner (a subject's thread with its units, a chain's thread with its groups) puts an item
//! set on the board and takes items until its set has none left to take, then blocks until the
//! items others took have run. A helper (a thread of the batch's subject loop without a subject,
//! or a `rayon::spawn` job of a group set) takes items until there is nothing left it may take.
//! Items are pure work: a unit's scan and a group's kernel make no pool call and wait for nothing,
//! so an item that was taken always ends, and a blocked owner always wakes. Taking an item never
//! runs a pool job, so no subject or chain runs nested inside another (the problem of
//! x_batch_pipe.rs item 6).
//!
//! Lifetimes: an item set's closure borrows its owner's frame; the board holds it as a raw pointer
//! that is read only for an index below the set's length, after a successful take, and the owner
//! returns only when every taken index has run (`XSet::wait`).

use std::any::Any;
use std::cell::Cell;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock};

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

/// What a taker lends the items it runs: the unit path keeps its scratch here, made by the first
/// unit the taker runs and dropped when the taker stops (as rayon's `map_init` made one per split).
#[derive(Default)]
pub(crate) struct XTaker {
    pub(crate) scratch: Option<Box<dyn Any>>,
}

type Run<'a> = dyn Fn(usize, &mut XTaker) + Sync + 'a;

/// A set of independent items: a subject's units, or a linking call's groups.
pub(crate) struct XSet {
    run: *const Run<'static>,
    n: usize,
    next: AtomicUsize,
    done: Mutex<usize>,
    all_done: Condvar,
}

// SAFETY: `run` is a `Sync` closure; it is called only for a taken index below `n`, and its owner
// does not return before every such call has returned (`XSet::wait`).
unsafe impl Send for XSet {}
unsafe impl Sync for XSet {}

impl XSet {
    fn new(run: &Run<'_>, n: usize) -> Arc<Self> {
        // SAFETY: the lifetime is erased for the takers; see the type's `Send` / `Sync`.
        let run: *const Run<'static> = unsafe { std::mem::transmute(run) };
        Arc::new(Self {
            run,
            n,
            next: AtomicUsize::new(0),
            done: Mutex::new(0),
            all_done: Condvar::new(),
        })
    }

    /// Takes and runs one item; false when none was left to take.
    fn run_one(&self, taker: &mut XTaker) -> bool {
        if self.next.load(Ordering::Acquire) >= self.n {
            return false;
        }
        let i = self.next.fetch_add(1, Ordering::AcqRel);
        if i >= self.n {
            return false;
        }
        // The item counts as run when it returns or panics, so its owner never waits for it forever.
        struct Ran<'a>(&'a XSet, bool);
        impl Drop for Ran<'_> {
            fn drop(&mut self) {
                IN_ITEM.with(|flag| flag.set(self.1));
                let mut done = self.0.done.lock().unwrap_or_else(|e| e.into_inner());
                *done += 1;
                if *done == self.0.n {
                    self.0.all_done.notify_all();
                }
            }
        }
        let _ran = Ran(self, IN_ITEM.with(|flag| flag.replace(true)));
        // SAFETY: see the type's `Send` / `Sync`.
        unsafe { (*self.run)(i, taker) };
        true
    }

    /// Every item has been taken.
    fn taken(&self) -> bool {
        self.next.load(Ordering::Acquire) >= self.n
    }

    /// Blocks until every item has run.
    fn wait(&self) {
        self.wait_for(self.n);
    }

    /// Blocks until `n` items have run.
    fn wait_for(&self, n: usize) {
        let mut done = self.done.lock().unwrap_or_else(|e| e.into_inner());
        while *done < n {
            done = self.all_done.wait(done).unwrap_or_else(|e| e.into_inner());
        }
    }
}

/// The owner's hold on its set: if the owner unwinds (a panic in an item it ran), no item is taken
/// any more and the guard waits for the items others took, which borrow the owner's frame, before
/// the panic goes on.
struct XOwned<'a>(&'a XSet);

impl Drop for XOwned<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            let taken = self
                .0
                .next
                .fetch_max(self.0.n, Ordering::AcqRel)
                .min(self.0.n);
            self.0.wait_for(taken);
        }
    }
}

/// The units of one batch's subjects, one slot per subject index.
pub(crate) struct XBatch {
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
    }
}

/// The board of one search.
#[derive(Default)]
pub(crate) struct XBoard {
    groups: Mutex<Vec<Arc<XSet>>>,
    batches: Mutex<Vec<Arc<XBatch>>>,
}

impl XBoard {
    /// A batch of `subjects` subjects, taken after the batches before it.
    pub(crate) fn batch(&self, subjects: usize) -> Arc<XBatch> {
        let batch = Arc::new(XBatch {
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
        loop {
            let set = {
                let mut groups = self.groups.lock().unwrap_or_else(|e| e.into_inner());
                groups.retain(|set| !set.taken());
                groups.first().cloned()
            };
            match set {
                Some(set) => {
                    if set.run_one(taker) {
                        return true;
                    }
                }
                None => break,
            }
        }
        let batches: Vec<Arc<XBatch>> = self
            .batches
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .clone();
        batches.iter().any(|batch| batch.run_one(taker))
    }

    /// Takes and runs one item (groups first, then units by batch and subject order); false when
    /// there was none. For a thread that waits outside the pool.
    pub(crate) fn take_one(&self, taker: &mut XTaker) -> bool {
        self.run_one(taker)
    }

    /// A subject's units: `run(i, taker)` for every `i < n`, on this thread, on the board's
    /// takers and on up to `threads - 1` helper jobs of the set, units of earlier subjects and
    /// groups first. Returns when every unit has run.
    pub(crate) fn run_units(&self, batch: &XBatch, s: usize, n: usize, run: &Run<'_>) {
        let set = XSet::new(run, n);
        let owned = XOwned(&set);
        *batch.slots[s].lock().unwrap_or_else(|e| e.into_inner()) = Some(Arc::clone(&set));
        batch.closed[s].store(true, Ordering::Release);
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
        while !set.taken() {
            if !self.run_one(&mut taker) {
                // Nothing before this subject's units is left and they are taken.
                set.run_one(&mut taker);
            }
        }
        drop(taker);
        set.wait();
        drop(owned);
        // The slot stays closed and drained; the set is not kept until the batch ends.
        *batch.slots[s].lock().unwrap_or_else(|e| e.into_inner()) = None;
    }

    /// A linking call's groups: `run(i, taker)` for every `i < n` (the largest group first is the
    /// caller's order), on this thread, on the board's takers before any unit, and on up to
    /// `helpers` pool jobs. Returns when every group has run.
    pub(crate) fn run_groups(&self, n: usize, run: &Run<'_>, helpers: usize) {
        let set = XSet::new(run, n);
        let owned = XOwned(&set);
        self.groups
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(Arc::clone(&set));
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        x_spawn_group_helper(Arc::clone(&set), helpers.min(n.saturating_sub(1)));
        let _ = helpers;
        let mut taker = XTaker::default();
        while set.run_one(&mut taker) {}
        set.wait();
        drop(owned);
    }

    /// A thread of `batch`'s subject loop without a subject: takes items until the batch's units
    /// are all taken (it waits, without the pool, while a subject has not yet put its units), or
    /// the loop stopped (`stopped`: a subject's search panicked, so the subjects not yet taken will
    /// never close).
    pub(crate) fn help(&self, batch: &XBatch, stopped: &AtomicBool) {
        let mut taker = XTaker::default();
        let mut idle = 0u32;
        loop {
            if self.run_one(&mut taker) {
                idle = 0;
                continue;
            }
            if batch.drained() || stopped.load(Ordering::Relaxed) {
                return;
            }
            idle += 1;
            if idle < 64 {
                std::hint::spin_loop();
            } else {
                std::thread::yield_now();
            }
        }
    }
}

/// One helper job of a group set, which spawns the next before it takes its item while items are
/// left, `left` helpers in all: idle threads join one after the other, and a set whose items were
/// all taken by others leaves at most one job in a thread's queue. One item per helper: a helper
/// that a thread runs inside a pool wait (nested) holds that thread for one item at most.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_spawn_group_helper(set: Arc<XSet>, left: usize) {
    if left == 0 || set.taken() {
        return;
    }
    rayon::spawn(move || {
        if !set.taken() {
            x_spawn_group_helper(Arc::clone(&set), left - 1);
        }
        let mut taker = XTaker::default();
        set.run_one(&mut taker);
    });
}

thread_local! {
    // This thread runs a board item: the item's sorts stay on the thread (no pool, no board).
    static IN_ITEM: Cell<bool> = const { Cell::new(false) };
}

/// Whether this thread is running a board item (an item makes no pool or board call).
pub(crate) fn in_item() -> bool {
    IN_ITEM.with(Cell::get)
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
pub(crate) fn current<'a>() -> Option<(&'a XBoard, &'a XBatch)> {
    // SAFETY: the pointers are set by `enter` from references that outlive its guard, and the
    // guard restores the value before them on drop; the subject's work that reads them runs
    // inside the guard's life on this thread.
    CURRENT
        .with(Cell::get)
        .map(|(board, batch)| unsafe { (&*board, &*batch) })
}

#[cfg(all(test, feature = "parallel", not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicU64;

    fn work(s: usize, i: usize) -> u64 {
        let rounds = 2_000 + ((s * 7_919 + i * 104_729) % 20_000);
        (0..rounds as u64).fold((s * 1_000 + i) as u64, |acc, x| {
            acc.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(x)
        })
    }

    // Subjects that put their units on the board from the threads of a broadcast (the TBXNEXT
    // loop), the other threads helping, and each subject linking a few groups on the board after its
    // units: every unit and group runs once, each subject gets its own results in unit order, and
    // the run ends at 1, 2, 4 and 8 threads with fewer, as many and more subjects than threads.
    #[test]
    fn every_item_runs_once_and_every_owner_gets_its_results() {
        for threads in [1usize, 2, 4, 8] {
            for subjects in [1usize, 2, 3, threads, threads + 3] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("test pool");
                let board = XBoard::default();
                let batch = board.batch(subjects);
                let next = AtomicUsize::new(0);
                let ran = AtomicU64::new(0);
                let checked = AtomicUsize::new(0);
                pool.install(|| {
                    rayon::broadcast(|_| {
                        loop {
                            let s = next.fetch_add(1, Ordering::Relaxed);
                            if s >= subjects {
                                break;
                            }
                            let _entered = enter(&board, &batch);
                            if s % 5 == 4 {
                                // A subject too short to search: no units.
                                batch.close_without_units(s);
                                checked.fetch_add(1, Ordering::Relaxed);
                                continue;
                            }
                            let units = 3 + (s * 11) % 40;
                            let out: Vec<OnceLock<u64>> =
                                (0..units).map(|_| OnceLock::new()).collect();
                            board.run_units(&batch, s, units, &|i, taker: &mut XTaker| {
                                taker.scratch.get_or_insert_with(|| Box::new(0u64));
                                ran.fetch_add(1, Ordering::Relaxed);
                                out[i].set(work(s, i)).expect("a unit runs once");
                            });
                            for (i, cell) in out.iter().enumerate() {
                                assert_eq!(cell.get(), Some(&work(s, i)), "subject {s} unit {i}");
                            }
                            let groups: Vec<OnceLock<u64>> =
                                (0..4).map(|_| OnceLock::new()).collect();
                            board.run_groups(
                                4,
                                &|g, _| {
                                    groups[g].set(work(s + 100, g)).expect("a group runs once");
                                },
                                threads - 1,
                            );
                            for (g, cell) in groups.iter().enumerate() {
                                assert_eq!(cell.get(), Some(&work(s + 100, g)));
                            }
                            checked.fetch_add(1, Ordering::Relaxed);
                        }
                        board.help(&batch, &AtomicBool::new(false));
                    });
                });
                board.end_batch(&batch);
                let expected: u64 = (0..subjects)
                    .filter(|s| s % 5 != 4)
                    .map(|s| (3 + (s * 11) % 40) as u64)
                    .sum();
                assert_eq!(
                    ran.load(Ordering::Relaxed),
                    expected,
                    "{threads} threads, {subjects} subjects"
                );
                assert_eq!(checked.load(Ordering::Relaxed), subjects);
            }
        }
    }

    // `drained` reads each slot once per batch (its cursor never goes back) and is true only when
    // every slot is closed and taken; a slot closed without units does not hold the unit order's
    // hint (`lowest`) back.
    #[test]
    fn drained_and_lowest_pass_closed_slots_once() {
        let board = XBoard::default();
        let batch = board.batch(4);
        batch.close_without_units(0);
        batch.close_without_units(2);
        assert!(!batch.drained());
        assert_eq!(batch.drained_below.load(Ordering::Relaxed), 1);
        let ran = AtomicUsize::new(0);
        board.run_units(&batch, 1, 3, &|_, _: &mut XTaker| {
            ran.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(ran.load(Ordering::Relaxed), 3);
        assert!(
            batch.slots[1].lock().unwrap().is_none(),
            "the slot is cleared"
        );
        assert!(!batch.drained());
        assert_eq!(batch.drained_below.load(Ordering::Relaxed), 3);
        assert!(!batch.run_one(&mut XTaker::default()));
        assert_eq!(batch.lowest.load(Ordering::Relaxed), 3);
        batch.close_without_units(3);
        assert!(batch.drained());
        assert_eq!(batch.drained_below.load(Ordering::Relaxed), 4);
    }

    // An item that panics on another thread counts as run, so its owner returns, and the panic
    // reaches the caller (a hang here would stop the test at the watchdog).
    #[test]
    fn a_panicking_item_on_another_thread_does_not_hang_its_owner() {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(2)
                .build()
                .expect("test pool");
            let board = XBoard::default();
            let other_took = AtomicBool::new(false);
            let owner_done = AtomicBool::new(false);
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                pool.install(|| {
                    rayon::broadcast(|ctx| {
                        if ctx.index() == 0 {
                            board.run_groups(
                                4,
                                &|_, _| {
                                    if rayon::current_thread_index() != Some(0) {
                                        other_took.store(true, Ordering::SeqCst);
                                        panic!("an item panics on the other thread");
                                    }
                                    let start = std::time::Instant::now();
                                    while !other_took.load(Ordering::SeqCst)
                                        && start.elapsed() < std::time::Duration::from_secs(5)
                                    {
                                        std::hint::spin_loop();
                                    }
                                },
                                0,
                            );
                            owner_done.store(true, Ordering::SeqCst);
                        } else {
                            let mut taker = XTaker::default();
                            while !owner_done.load(Ordering::SeqCst) {
                                if !board.take_one(&mut taker) {
                                    std::thread::yield_now();
                                }
                            }
                        }
                    });
                });
            }));
            let _ = tx.send((result.is_err(), other_took.load(Ordering::SeqCst)));
        });
        let (panicked, other_took) = rx
            .recv_timeout(std::time::Duration::from_secs(60))
            .expect("the owner returned");
        if other_took {
            assert!(panicked, "the item's panic reaches the caller");
        }
    }
}
