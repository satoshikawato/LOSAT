//! No NCBI counterpart: a set of independent items that its owner runs with the idle threads of the
//! pool, where the owner never runs a pool job while it waits.
//!
//! The items are independent pieces of one NCBI loop (the units of a TBLASTX subject and the groups
//! of a linking call, `tblastx/x_work_board.rs`; the contexts of a BLASTX query batch that SEG filters
//! and whose statistical blocks are computed, and the chunks of a split query that are set up,
//! `blastx/query_setup.rs`, `blastx/parameters.rs`, `blastx/split.rs`), whose results do not depend
//! on the thread that runs them or on when it runs them. This file only chooses which thread runs
//! which item, and when.
//!
//! # Why
//!
//! A thread that waits in a pool call (`join`, `par_iter`, the end of a `scope`) runs other jobs of
//! the pool meanwhile, from its own queue or stolen from another thread's, and returns only when the
//! job it took has ended. Inside a flow of long tasks (BLASTX searches every query batch as a task,
//! `blastx/native.rs` `x_flow`, and its composition redo hands out helpers that run until the loop is
//! claimed) such a job is another batch's search or helper: the 120 contexts of one 80,000-letter
//! batch, 5 ms of SEG on one thread, took seconds to filter at 16 threads (`q200k`, `BATCH_SIZE`
//! 40000, `-use_sw_tback`: 19.6 to 52.8 s of SEG for 5 batches), and the batch started its search
//! that much later.
//!
//! # How (no thread waits in the pool)
//!
//! An owner (a subject's thread with its units, a chain's thread with its groups, a batch's thread
//! with its contexts) puts an item set on the board and takes its set's items only, until none is
//! left to take, then blocks until the items others took have run. A helper of the board (a thread
//! of the batch's subject loop without a subject) takes items until there is nothing left it may
//! take, and a `rayon::spawn` job of the board takes one group; a helper job of an `x_each` set
//! queues the set's next helper job and takes one item. Items are pure work: they make no pool call
//! and wait for nothing, so an item that was taken always ends, and a blocked owner always wakes.
//! Taking an item never runs a pool job, so no long task runs nested inside an owner's wait (the
//! problem of `tblastx/x_batch_pipe.rs` item 6).
//!
//! Lifetimes: an item set's closure borrows its owner's frame; the set holds it as a raw pointer
//! that is read only for an index below the set's length, after a successful take, and the owner
//! returns only when every taken index has run (`XSet::wait`). A helper job that starts after its
//! owner returned finds no index left to take and reads nothing of the frame.

use std::any::Any;
use std::cell::Cell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};

/// What a taker lends the items it runs: the unit path keeps its scratch here, made by the first
/// unit the taker runs and dropped when the taker stops (as rayon's `map_init` made one per split).
#[derive(Default)]
pub(crate) struct XTaker {
    pub(crate) scratch: Option<Box<dyn Any>>,
}

pub(crate) type Run<'a> = dyn Fn(usize, &mut XTaker) + Sync + 'a;

/// A set of independent items: a subject's units, a linking call's groups, a batch's contexts.
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
    pub(crate) fn new(run: &Run<'_>, n: usize) -> Arc<Self> {
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

    /// Takes the next index; None when none was left to take.
    fn take(&self) -> Option<usize> {
        if self.next.load(Ordering::Acquire) >= self.n {
            return None;
        }
        let i = self.next.fetch_add(1, Ordering::AcqRel);
        (i < self.n).then_some(i)
    }

    /// Takes and runs one item; false when none was left to take.
    pub(crate) fn run_one(&self, taker: &mut XTaker) -> bool {
        let Some(i) = self.take() else {
            return false;
        };
        let _ran = XRan::start(self);
        // SAFETY: see the type's `Send` / `Sync`.
        unsafe { (*self.run)(i, taker) };
        true
    }

    /// `run_one` for a helper job, which must not unwind: a panic of the item goes to `on_panic`
    /// before the item counts as run, so an owner that wakes because every item has run finds it.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn run_one_catching(
        &self,
        taker: &mut XTaker,
        on_panic: impl FnOnce(Box<dyn Any + Send>),
    ) -> bool {
        let Some(i) = self.take() else {
            return false;
        };
        let _ran = XRan::start(self);
        if let Err(panic) = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            // SAFETY: see the type's `Send` / `Sync`.
            unsafe { (*self.run)(i, taker) }
        })) {
            on_panic(panic);
        }
        true
    }

    /// Every item has been taken.
    pub(crate) fn taken(&self) -> bool {
        self.next.load(Ordering::Acquire) >= self.n
    }

    /// Blocks until every item has run.
    pub(crate) fn wait(&self) {
        let mut done = self.done.lock().unwrap_or_else(|e| e.into_inner());
        while *done < self.n {
            done = self.all_done.wait(done).unwrap_or_else(|e| e.into_inner());
        }
    }

    /// Blocks until `n` items have run, for an owner that unwinds: an item that ends notifies only
    /// when it is the set's last (a syscall per unit otherwise), so this wait reads the count again.
    fn wait_unwinding(&self, n: usize) {
        let mut done = self.done.lock().unwrap_or_else(|e| e.into_inner());
        while *done < n {
            done = self
                .all_done
                .wait_timeout(done, std::time::Duration::from_millis(1))
                .unwrap_or_else(|e| e.into_inner())
                .0;
        }
    }
}

/// An item being run: it counts as run when it returns or panics (the drop), so its owner never
/// waits for it forever, and `in_item` holds on its thread meanwhile.
struct XRan<'a>(&'a XSet, bool);

impl<'a> XRan<'a> {
    fn start(set: &'a XSet) -> Self {
        Self(set, IN_ITEM.with(|flag| flag.replace(true)))
    }
}

impl Drop for XRan<'_> {
    fn drop(&mut self) {
        IN_ITEM.with(|flag| flag.set(self.1));
        let mut done = self.0.done.lock().unwrap_or_else(|e| e.into_inner());
        *done += 1;
        if *done == self.0.n {
            self.0.all_done.notify_all();
        }
    }
}

/// The owner's hold on its set: if the owner unwinds (a panic in an item it ran), no item is taken
/// any more and the guard waits for the items others took, which borrow the owner's frame, before
/// the panic goes on.
pub(crate) struct XOwned<'a>(pub(crate) &'a XSet);

impl Drop for XOwned<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            let taken = self
                .0
                .next
                .fetch_max(self.0.n, Ordering::AcqRel)
                .min(self.0.n);
            self.0.wait_unwinding(taken);
        }
    }
}

thread_local! {
    // This thread runs an item: the item's own parallel steps stay on the thread (no pool, no
    // board, no set of its own).
    static IN_ITEM: Cell<bool> = const { Cell::new(false) };
}

/// Whether this thread is running an item (an item makes no pool call).
pub(crate) fn in_item() -> bool {
    IN_ITEM.with(Cell::get)
}

/// A set of `x_each` and its helper jobs.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XHelped {
    set: Arc<XSet>,
    /// A helper job waits in a thread's queue.
    queued: std::sync::atomic::AtomicBool,
    /// The first panic of an item that a helper job ran, raised again by the owner (a panic in a
    /// `rayon::spawn` job would end the process).
    panic: Mutex<Option<Box<dyn Any + Send>>>,
}

/// The set's helper job, one in the threads' queues at a time: a thread that runs it spawns the
/// next while items are left, so idle threads join one after the other, and a set whose items its
/// owner ran alone leaves one job behind, which ends at once. One item per job: a job that a thread
/// runs inside a pool wait holds that thread for one item at most.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_spawn_helper(helped: &Arc<XHelped>) {
    if helped.set.taken() || helped.queued.swap(true, Ordering::AcqRel) {
        return;
    }
    let helped = Arc::clone(helped);
    rayon::spawn(move || {
        helped.queued.store(false, Ordering::Release);
        x_spawn_helper(&helped);
        let mut taker = XTaker::default();
        // The payload is stored before the item counts as run (`run_one_catching`), so the owner,
        // which reads it once every item has run, cannot miss it.
        helped.set.run_one_catching(&mut taker, |panic| {
            helped
                .panic
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .get_or_insert(panic);
        });
    });
}

/// `run(i)` for every `i < n`; returns when every item has run. In a pool of 2+ threads this thread
/// takes the items in index order and the pool's idle threads take the next ones (the set's helper
/// jobs); this thread never runs a pool job, and once no item is left to take it blocks until the
/// items other threads took have run. Outside a pool, in a pool of one thread, and inside an item
/// (`in_item`) the items run here one after the other. A panic of an item is raised here.
pub(crate) fn x_each(n: usize, run: &(dyn Fn(usize) + Sync)) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if n > 1
        && !in_item()
        && rayon::current_thread_index().is_some()
        && rayon::current_num_threads() > 1
    {
        let item = |i: usize, _: &mut XTaker| run(i);
        let set = XSet::new(&item, n);
        let owned = XOwned(&set);
        let helped = Arc::new(XHelped {
            set: Arc::clone(&set),
            queued: std::sync::atomic::AtomicBool::new(false),
            panic: Mutex::new(None),
        });
        x_spawn_helper(&helped);
        let mut taker = XTaker::default();
        while set.run_one(&mut taker) {}
        set.wait();
        drop(owned);
        let panic = helped
            .panic
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .take();
        if let Some(panic) = panic {
            std::panic::resume_unwind(panic);
        }
        return;
    }
    for i in 0..n {
        run(i);
    }
}

/// `f(i)` for every `i < n` (`x_each`), in index order.
#[cfg_attr(
    not(feature = "parallel"),
    expect(dead_code, reason = "its callers are the parallel steps of BLASTX")
)]
pub(crate) fn x_map<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let slots: Vec<Mutex<Option<T>>> = (0..n).map(|_| Mutex::new(None)).collect();
    x_each(n, &|i| {
        let value = f(i);
        *slots[i].lock().unwrap_or_else(|e| e.into_inner()) = Some(value);
    });
    slots
        .into_iter()
        .map(|slot| {
            slot.into_inner()
                .unwrap_or_else(|e| e.into_inner())
                .expect("every item ran")
        })
        .collect()
}

#[cfg(all(test, feature = "parallel", not(target_arch = "wasm32")))]
#[path = "xitems_tests.rs"]
mod tests;
