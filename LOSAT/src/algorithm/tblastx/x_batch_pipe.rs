//! EXPERIMENT (LOSAT_X_TBXPIPE): the next query batch of a TBLASTX
//! search starts when the batch before it has ended its seed stage, while that batch still links
//! its subjects.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_subject_next.rs` is).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-195
//! ```c
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!     ...
//!                 results = lcl_blast.Run();
//! ```
//! NCBI searches its query batches one after the other, and so does `run_in_pool`. With TBXNEXT
//! every subject of a batch starts at once, one pool thread each, and its units run on the whole
//! pool, so the pool is full while units run. After a batch's last unit only its subjects' chains
//! are left (merge, both linking calls, the re-evaluation): on the genome 3 x 3 input at 16
//! threads that is 3.6 s of an 11.2 s search, with 45 of its 179 thread-seconds idle; the first
//! batch's tail alone is one subject's linking (2.0 s) while 13 threads wait. The work that can
//! fill it is the next batch's units.
//!
//! # What changes
//!
//! `run_in_pool` keeps its loop: its checks, its errors, its writes and their order. Only where a
//! batch's result comes from changes (`XPipe::take`):
//!
//! 1. A plan lists the batches the loop will search unless a search fails (`run_impl.rs`
//!    `x_pipe_plan`: the loop's boundaries and the checks that need no search).
//! 2. The loop searches a batch itself when it reaches it, unless the batch before already
//!    started it. A batch starts the next one at its *seed-stage end*: when the last of its
//!    subjects has ended its seed stage (`XSeedEnd`, at dispatch point 5 of `process_subject`,
//!    where the seed release (`x_seed_release.rs`) has freed the batch's lookup and unit tables),
//!    if the admission rule lets it. That search is a job of the loop's scope
//!    (`rayon::in_place_scope`) and runs on the pool like any other. With `LOSAT_X_TBXBOARD`
//!    (tblastx/x_work_board.rs) a batch whose units are on the board starts the next one earlier,
//!    when every one of its units has been taken (`XSeedEnd::units_taken`, while the HSPs held are
//!    below the limit): the next batch builds its tables while the last units run, then waits at a
//!    gate before its subjects (`XGateView`) until the batch before has ended its seed stage and
//!    the HSPs held leave room, or every unit of the batch before has run and the records they
//!    returned (a bound on the HSPs that batch's merges keep) fit beside the HSPs held, or the
//!    batch before has returned, or the pipe stopped, taking the board's items meanwhile.
//! 3. The loop takes a started batch's result from its slot. Until the slot is filled the loop's
//!    thread runs the pool's jobs (`rayon::yield_now`), and sleeps 50 us when there is none, so
//!    no thread waits for a condition outside rayon's latches and no job waits at all (plan v4
//!    §5: a waiting helper that a nested wait runs deadlocks). With the board the loop's thread
//!    and a gate's thread take the board's items instead, then their own pool jobs and their
//!    part of a broadcast (`rayon::yield_local`; a later batch's subject loop is a
//!    `rayon::broadcast`, which ends only after every thread has run its part), so they never
//!    run another thread's jobs, such as a batch's search, nested.
//! 4. Admission: a subject holds its HSP list from its seed-stage end until it returns (its chain
//!    links it, re-evaluates it and converts it). The next batch starts early only while the
//!    subjects past their seed stage hold at most `HSP_LIMIT` HSPs (4 Mi: 352 MiB of 88-byte
//!    records); otherwise
//!    the loop starts it when it reaches it, as before. The genome 3 x 3 batches hold 1.2-1.4 M
//!    HSPs; one 5 Mb by 5 Mb pair holds about 23 M, so a query file of such genomes keeps today's
//!    peak. Without the board the tables of two batches are never held at once: the next batch
//!    builds its own after the batch before has freed its tables. With the board the next batch's
//!    tables are built while the batch before's last units run, so both batches' tables are held
//!    for that time (one batch ahead at most: a batch starts the next only once its own units are
//!    all taken).
//! 5. A search that returns an error stops the early starts. The loop then meets that batch's
//!    error where it met it before; batches started after it are discarded unread. So does the
//!    loop's end (its own errors, as a database's titles that cannot be read): the scope then
//!    waits only for the batches already started. A search that panics fills its slot with the
//!    panic, which the loop raises again when it takes the batch (release builds abort first).
//! 6. While a later batch may run, a batch's subjects link without the pool (their two linking
//!    calls get `allow_parallel` false, the path of one thread, which gives the same lists): a
//!    thread waiting for its pool tasks runs whatever the pool holds, so a chain that waited for
//!    its link groups ran the next batch's subjects inside the wait and ended after them, and the
//!    chains of three batches ended one after the other on one thread at the end of the search
//!    (genome 3 x 3, 16 threads: the last batch's tail 1.43 s, against 0.91 s without the pipe).
//!    A chain without a wait never runs other work. The other threads run the next batch's units
//!    meanwhile. The subjects of the plan's last batch, and of a batch whose subjects hold more
//!    HSPs than the limit (no batch starts early then), link on the pool as before. With the
//!    board, such a chain puts its link groups and sorts on the board instead of running them on
//!    its own thread (`x_link_slices.rs` `x_helped`), where no thread waits in the pool.
//!
//! # Why the output is the loop's
//!
//! A batch's result is a function of its queries, the subjects and the run's arguments: every
//! table, cutoff and counter its search reads is a local of `search_query_batch` (the
//! process-wide values are switches read once and diagnostic counters), and its result does not
//! depend on the thread count or on which thread searched which subject (the threads tests, 1 to
//! 16 threads). The loop takes the batches it would have searched, in its order, under the key it
//! computed itself (`take` compares it with the plan's; a difference, which no input is known to
//! give, searches the batch on the loop and stops the early starts). The `LOSAT_TIMING` block of
//! a batch is kept as text in its result and printed when the loop takes it, in the batch order.
//!
//! A thread that waits inside a batch (a subject's units, a link1 group) may run the next batch's
//! jobs nested; that batch's subject then ends later, which only delays the loop's take and the
//! admission count. No thread-local scratch of TBLASTX is borrowed across a rayon wait (the link
//! kernels make no rayon call; the subject buffers are taken and put back).
//!
//! Off (the loop as before): the switch off; a build without `parallel`, or wasm32 (the loop is
//! kept there); one thread, or a loop not on a thread of its pool; fewer than two planned
//! batches; any debugging output of the search (`LOSAT_DEBUG*`, the HSP trace, the diagnostics,
//! the stage dumps), whose lines would interleave.
//!
//! Diagnostics: with `LOSAT_X_STATS`, `[X_STATS] TBXPIPE planned= taken= early= loop= waited=
//! key_mismatch= max_held_hsps=` at the end of the search (early: started at the batch before's
//! seed-stage end, or with the board when its units were all taken; loop: searched by the loop
//! when it reached them; waited: taken after the loop waited for the slot).

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use std::sync::Mutex;
use std::sync::OnceLock;

use anyhow::Result;

// No NCBI counterpart: the variable is read once per process; it chooses when a batch's search
// starts, not a value NCBI computes.
/// LOSAT_X_TBXPIPE: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXPIPE"))
}

#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
/// The HSPs the batches past their seed-stage end may hold when the next batch starts early
/// (module header, item 4).
const HSP_LIMIT: u64 = 1 << 22;

#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
// No NCBI counterpart: a tuning and test override of `HSP_LIMIT` (`LOSAT_X_TBXPIPEHSPS=<n>`).
fn hsp_limit() -> u64 {
    static LIMIT: OnceLock<u64> = OnceLock::new();
    *LIMIT.get_or_init(|| {
        std::env::var("LOSAT_X_TBXPIPEHSPS")
            .ok()
            .and_then(|value| value.trim().parse::<u64>().ok())
            .unwrap_or(HSP_LIMIT)
    })
}

// No NCBI counterpart: the pipe's count of the HSPs held (module header, item 4).
/// The HSPs that subjects past their seed stage hold, the most held at once, and the admission
/// limit.
#[derive(Default)]
pub(super) struct XHeld {
    now: AtomicU64,
    most: AtomicU64,
    limit: u64,
}

impl XHeld {
    fn with_limit(limit: u64) -> Self {
        Self {
            limit,
            ..Self::default()
        }
    }

    /// Adds `hsps`; the HSPs held now.
    fn add(&self, hsps: u64) -> u64 {
        let now = self.now.fetch_add(hsps, Ordering::SeqCst) + hsps;
        self.most.fetch_max(now, Ordering::Relaxed);
        now
    }

    fn release(&self, hsps: u64) {
        self.now.fetch_sub(hsps, Ordering::SeqCst);
    }

    fn room(&self) -> bool {
        self.now.load(Ordering::SeqCst) <= self.limit
    }
}

// No NCBI counterpart: the stage of one batch search that `process_subject` reports to the pipe.
/// The seed-stage end of one batch search: the subjects whose seed stage has not ended, the HSPs
/// of those that have, and what to do when the last one ends.
pub(super) struct XSeedEnd<'a> {
    left: AtomicUsize,
    hsps: AtomicU64,
    on_end: &'a (dyn Fn(u64) + Sync),
    held: &'a XHeld,
    /// The batch is the plan's last: no batch can start while it links.
    last: bool,
    // EXPERIMENT (LOSAT_X_TBXBOARD): what to do once every unit of the batch has been taken
    // (start the next batch's setup), once; and the gate this batch's subjects wait at when the
    // batch was started that early.
    units_taken: Option<&'a (dyn Fn() + Sync)>,
    units_taken_once: AtomicBool,
    gate: Option<XGateView<'a>>,
    // EXPERIMENT (LOSAT_X_TBXBOARD): the subjects whose units have not all run, the records their
    // units returned (a bound on the HSPs the merges keep), and what to do when the last one ends.
    units_left: AtomicUsize,
    unit_records: AtomicU64,
    units_done: Option<&'a (dyn Fn(u64) + Sync)>,
}

// EXPERIMENT (LOSAT_X_TBXBOARD): a batch started when the batch before had all its units taken
// builds its tables, then its subjects wait here until the batch before has ended its seed stage
// and the HSPs held leave room (module header, item 4), or every unit of the batch before has run
// and its units' records fit beside the HSPs held, or the batch before has returned, or the pipe
// stopped. No NCBI counterpart (scheduling only).
#[derive(Default)]
pub(super) struct XGate {
    prev_seed_done: AtomicBool,
    prev_returned: AtomicBool,
    /// Every unit of the batch before has run; `prev_bound` is the records they returned, which
    /// bounds the HSPs the batch before will hold from its seed-stage end.
    prev_units_done: AtomicBool,
    prev_bound: AtomicU64,
}

#[derive(Clone, Copy)]
pub(super) struct XGateView<'a> {
    gate: &'a XGate,
    held: &'a XHeld,
    stop: &'a AtomicBool,
}

impl XGateView<'_> {
    fn open(&self) -> bool {
        self.stop.load(Ordering::SeqCst)
            || self.gate.prev_returned.load(Ordering::SeqCst)
            || (self.gate.prev_seed_done.load(Ordering::SeqCst) && self.held.room())
            || (self.gate.prev_units_done.load(Ordering::SeqCst)
                && self.held.now.load(Ordering::SeqCst)
                    + self.gate.prev_bound.load(Ordering::SeqCst)
                    <= self.held.limit)
    }
}

impl<'a> XSeedEnd<'a> {
    /// `subjects` seed stages counted into `held`; `on_end` gets the batch's HSPs when the last
    /// one ends; `last` when no batch follows in the plan.
    pub(super) fn new(
        subjects: usize,
        on_end: &'a (dyn Fn(u64) + Sync),
        held: &'a XHeld,
        last: bool,
    ) -> Self {
        Self {
            left: AtomicUsize::new(subjects),
            hsps: AtomicU64::new(0),
            on_end,
            held,
            last,
            units_taken: None,
            units_taken_once: AtomicBool::new(false),
            gate: None,
            units_left: AtomicUsize::new(subjects),
            unit_records: AtomicU64::new(0),
            units_done: None,
        }
    }

    // EXPERIMENT (LOSAT_X_TBXBOARD): see the fields.
    fn with_early(
        mut self,
        units_taken: &'a (dyn Fn() + Sync),
        units_done: &'a (dyn Fn(u64) + Sync),
        gate: Option<XGateView<'a>>,
    ) -> Self {
        self.units_taken = Some(units_taken);
        self.units_done = Some(units_done);
        self.gate = gate;
        self
    }

    /// One subject's units have all run and returned `records` records (its HSPs for a subject
    /// without units); the last one tells the pipe the batch's bound.
    fn subject_units_done(&self, records: u64) {
        let total = self.unit_records.fetch_add(records, Ordering::SeqCst) + records;
        if self.units_left.fetch_sub(1, Ordering::SeqCst) == 1 {
            if let Some(units_done) = self.units_done {
                units_done(total);
            }
        }
    }

    /// Every unit of this batch has been taken (from the board): the next batch may build its
    /// tables. Once per batch.
    pub(super) fn units_taken(&self) {
        if let Some(units_taken) = self.units_taken {
            if !self.units_taken_once.swap(true, Ordering::SeqCst) {
                units_taken();
            }
        }
    }

    /// Before this batch's subjects: waits at the gate when the batch started early, running the
    /// board's items meanwhile (no pool wait).
    pub(super) fn wait_gate(
        &self,
        board: Option<&crate::algorithm::tblastx::x_work_board::XBoard>,
    ) {
        let Some(gate) = self.gate else {
            return;
        };
        let mut taker = crate::algorithm::tblastx::x_work_board::XTaker::default();
        while !gate.open() {
            if board.is_some_and(|board| board.take_one(&mut taker)) {
                continue;
            }
            // No unit to run now: a unit's scratch is not held through the chains' peak.
            taker.scratch = None;
            // This thread's own pool jobs and broadcast parts (see `XPipe::wait`).
            #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
            if matches!(rayon::yield_local(), Some(rayon::Yield::Executed)) {
                continue;
            }
            std::thread::sleep(std::time::Duration::from_micros(50));
        }
    }

    /// One subject's seed stage has ended with `hsps` HSPs, which it holds until it returns.
    /// True when its chain may use the pool (module header, item 6).
    fn subject_ended(&self, hsps: u64) -> bool {
        let held = self.held.add(hsps);
        let total = self.hsps.fetch_add(hsps, Ordering::SeqCst) + hsps;
        if self.left.fetch_sub(1, Ordering::SeqCst) == 1 {
            (self.on_end)(total);
        }
        self.last || held > self.held.limit
    }
}

// No NCBI counterpart: one subject's part of `XSeedEnd`.
/// The seed stage of one subject: `end` with its HSP count, or the drop (a subject that returns
/// before the seed stage ends, as one too short to search) with none; the drop releases the HSPs
/// from the count.
pub(super) struct XSeedGuard<'a> {
    seed_end: Option<&'a XSeedEnd<'a>>,
    ended: bool,
    hsps: u64,
    units_reported: bool,
}

impl<'a> XSeedGuard<'a> {
    pub(super) fn new(seed_end: Option<&'a XSeedEnd<'a>>) -> Self {
        Self {
            seed_end,
            ended: false,
            hsps: 0,
            units_reported: false,
        }
    }

    // EXPERIMENT (LOSAT_X_TBXBOARD): this subject's units have all run and returned `records`
    // records (`XSeedEnd::subject_units_done`); once.
    pub(super) fn units_done(&mut self, records: usize) {
        if self.units_reported {
            return;
        }
        self.units_reported = true;
        if let Some(seed_end) = self.seed_end {
            seed_end.subject_units_done(records as u64);
        }
    }

    /// The subject's seed stage has ended with `hsps` HSPs. True when the subject's chain (its
    /// linking calls) may use the pool: always without the pipe; with it, in the plan's last batch
    /// and when the HSPs held exceed the admission limit (no batch starts early then), else false:
    /// a later batch runs meanwhile, and a chain waiting for its pool tasks would run that batch's
    /// subjects inside the wait and end after them.
    pub(super) fn end(&mut self, hsps: usize) -> bool {
        if self.ended {
            return true;
        }
        // A subject that ran no units reports its HSPs as its bound.
        self.units_done(hsps);
        self.ended = true;
        match self.seed_end {
            Some(seed_end) => {
                self.hsps = hsps as u64;
                seed_end.subject_ended(self.hsps)
            }
            None => true,
        }
    }
}

impl Drop for XSeedGuard<'_> {
    fn drop(&mut self) {
        // An unwinding subject starts nothing.
        if std::thread::panicking() {
            return;
        }
        self.end(0);
        if let Some(seed_end) = self.seed_end {
            seed_end.held.release(self.hsps);
        }
    }
}

#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
// No NCBI counterpart: see `XPipe::closer`.
/// Stops the early starts when dropped: when the loop returns, and when it unwinds.
pub(super) struct XPipeCloser<'a>(&'a AtomicBool);

#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
impl Drop for XPipeCloser<'_> {
    fn drop(&mut self) {
        self.0.store(true, Ordering::SeqCst);
    }
}

/// What the pipe counts (`[X_STATS] TBXPIPE`).
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
#[derive(Default)]
struct Counts {
    taken: AtomicU64,
    early: AtomicU64,
    loop_searched: AtomicU64,
    waited: AtomicU64,
    key_mismatch: AtomicU64,
}

/// A batch search: the batch's query range and its seed-stage end, if the pipe follows it.
pub(super) type XSearch<'s, R> =
    dyn Fn(std::ops::Range<usize>, Option<&XSeedEnd<'_>>) -> Result<R> + Sync + 's;

#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
// No NCBI counterpart: see the module header.
/// The planned batches of one search, their slots, and the admission count.
pub(super) struct XPipe<'s, R> {
    plan: Vec<std::ops::Range<usize>>,
    subjects: usize,
    search: &'s XSearch<'s, R>,
    started: Vec<AtomicBool>,
    /// A batch's result, or the panic of its search (`run_early`).
    slots: Vec<Mutex<Option<std::thread::Result<Result<R>>>>>,
    next: AtomicUsize,
    stop: AtomicBool,
    held: XHeld,
    counts: Counts,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): one gate per planned batch; `gated[i]` when batch i was
    /// started at the batch before's units-taken point.
    gates: Vec<XGate>,
    gated: Vec<AtomicBool>,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): the search's board; the loop waits for a slot by taking its
    /// items instead of the pool's jobs (it never runs a batch's search nested).
    board: Option<&'s crate::algorithm::tblastx::x_work_board::XBoard>,
}

#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
impl<'s, R: Send> XPipe<'s, R> {
    /// The pipe of `plan` (the batches' query ranges, in the loop's order) over `subjects`
    /// subjects; `None` when the plan has fewer than two batches.
    pub(super) fn new(
        plan: Vec<std::ops::Range<usize>>,
        subjects: usize,
        search: &'s XSearch<'s, R>,
        board: Option<&'s crate::algorithm::tblastx::x_work_board::XBoard>,
    ) -> Option<Self> {
        Self::with_limit(plan, subjects, search, hsp_limit()).map(|mut pipe| {
            pipe.board = board;
            pipe
        })
    }

    fn with_limit(
        plan: Vec<std::ops::Range<usize>>,
        subjects: usize,
        search: &'s XSearch<'s, R>,
        limit: u64,
    ) -> Option<Self> {
        if plan.len() < 2 {
            crate::utils::xreach::x_reach!("LOSAT_X_TBXPIPE", false);
            return None;
        }
        let n = plan.len();
        Some(Self {
            plan,
            subjects,
            search,
            started: (0..n).map(|_| AtomicBool::new(false)).collect(),
            slots: (0..n).map(|_| Mutex::new(None)).collect(),
            next: AtomicUsize::new(0),
            stop: AtomicBool::new(false),
            held: XHeld::with_limit(limit),
            counts: Counts::default(),
            gates: (0..n).map(|_| XGate::default()).collect(),
            gated: (0..n).map(|_| AtomicBool::new(false)).collect(),
            board: None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:193-195
    // ```c
    //                 results = lcl_blast.Run();
    // ```
    /// The result of the loop's next batch, `range`: from its slot if the batch before started
    /// it, else searched here. Call it from the loop, inside `scope`, once per batch in order.
    pub(super) fn take(
        &'s self,
        scope: &rayon::Scope<'s>,
        range: std::ops::Range<usize>,
    ) -> Result<R> {
        let i = self.next.fetch_add(1, Ordering::SeqCst);
        if self.plan.get(i) != Some(&range) {
            self.counts.key_mismatch.fetch_add(1, Ordering::Relaxed);
            self.stop.store(true, Ordering::SeqCst);
            return (self.search)(range, None);
        }
        self.counts.taken.fetch_add(1, Ordering::Relaxed);
        if self.started[i].swap(true, Ordering::SeqCst) {
            self.wait(i)
        } else {
            self.counts.loop_searched.fetch_add(1, Ordering::Relaxed);
            self.run(scope, i)
        }
    }

    /// No batch starts early any more (the loop has ended); dropped by the loop's scope.
    pub(super) fn closer(&self) -> XPipeCloser<'_> {
        XPipeCloser(&self.stop)
    }

    /// The slot of batch `i`, once its search has filled it; the thread runs the pool's jobs
    /// meanwhile. A panic of that search is raised again here.
    fn wait(&self, i: usize) -> Result<R> {
        let mut waited = false;
        // One taker for the whole wait: the units it runs keep one scratch (`XTaker`).
        let mut taker = crate::algorithm::tblastx::x_work_board::XTaker::default();
        loop {
            if let Some(result) = self.slots[i]
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .take()
            {
                if waited {
                    self.counts.waited.fetch_add(1, Ordering::Relaxed);
                }
                return result.unwrap_or_else(|panic| std::panic::resume_unwind(panic));
            }
            waited = true;
            if let Some(board) = self.board {
                // The board's items, then this thread's own pool jobs and its part of a broadcast
                // (a later batch's subject loop waits until every thread has run its part), never
                // another thread's jobs.
                if board.take_one(&mut taker) {
                    continue;
                }
                // No unit to run now: a unit's scratch is not held through the chains' peak.
                taker.scratch = None;
                if !matches!(rayon::yield_local(), Some(rayon::Yield::Executed)) {
                    std::thread::sleep(std::time::Duration::from_micros(50));
                }
                continue;
            }
            match rayon::yield_now() {
                Some(rayon::Yield::Executed) => {}
                _ => std::thread::sleep(std::time::Duration::from_micros(50)),
            }
        }
    }

    /// Searches batch `i` with its seed-stage end, which may start batch `i + 1`.
    fn run(&'s self, scope: &rayon::Scope<'s>, i: usize) -> Result<R> {
        let on_end = |_: u64| self.seed_ended(scope, i);
        let units_taken = || self.units_taken(scope, i);
        let units_done = |bound: u64| {
            if let Some(next) = self.gates.get(i + 1) {
                next.prev_bound.store(bound, Ordering::SeqCst);
                next.prev_units_done.store(true, Ordering::SeqCst);
            }
        };
        let last = i + 1 == self.plan.len();
        let gate = self.gated[i].load(Ordering::SeqCst).then(|| XGateView {
            gate: &self.gates[i],
            held: &self.held,
            stop: &self.stop,
        });
        let seed_end = XSeedEnd::new(self.subjects, &on_end, &self.held, last).with_early(
            &units_taken,
            &units_done,
            gate,
        );
        let result = (self.search)(self.plan[i].clone(), Some(&seed_end));
        if result.is_err() {
            self.stop.store(true, Ordering::SeqCst);
        }
        if let Some(next) = self.gates.get(i + 1) {
            next.prev_returned.store(true, Ordering::SeqCst);
        }
        result
    }

    // EXPERIMENT (LOSAT_X_TBXBOARD): every unit of batch `i` has been taken: batch i + 1 starts
    // now, builds its tables and waits at its gate (`XGateView::open`) before its subjects.
    fn units_taken(&'s self, scope: &rayon::Scope<'s>, i: usize) {
        let next = i + 1;
        // Batch i's own HSPs are not counted yet (they are at its seed-stage end): the start needs
        // room left below the limit (item 4).
        if next >= self.plan.len()
            || self.stop.load(Ordering::SeqCst)
            || self.started[next].load(Ordering::SeqCst)
            || self.held.now.load(Ordering::SeqCst) >= self.held.limit
        {
            return;
        }
        self.gated[next].store(true, Ordering::SeqCst);
        if self.started[next].swap(true, Ordering::SeqCst) {
            return;
        }
        crate::utils::xreach::x_reach!("LOSAT_X_TBXPIPE", true);
        self.counts.early.fetch_add(1, Ordering::Relaxed);
        self.spawn(scope, next);
    }

    /// Batch `i` has ended its seed stage: batch `i + 1` starts now if the HSPs held leave room
    /// and nothing started it yet.
    fn seed_ended(&'s self, scope: &rayon::Scope<'s>, i: usize) {
        let next = i + 1;
        if next >= self.plan.len() {
            return;
        }
        self.gates[next]
            .prev_seed_done
            .store(true, Ordering::SeqCst);
        if !self.held.room() {
            // The HSPs held leave no room: the loop starts the next batch itself.
            crate::utils::xreach::x_reach!("LOSAT_X_TBXPIPE", false);
            return;
        }
        if self.stop.load(Ordering::SeqCst) || self.started[next].swap(true, Ordering::SeqCst) {
            return;
        }
        crate::utils::xreach::x_reach!("LOSAT_X_TBXPIPE", true);
        self.counts.early.fetch_add(1, Ordering::Relaxed);
        self.spawn(scope, next);
    }

    fn spawn(&'s self, scope: &rayon::Scope<'s>, next: usize) {
        scope.spawn(move |scope| {
            // The loop waits for this slot: a panic fills it too (item 5).
            let result =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| self.run(scope, next)));
            if result.is_err() {
                self.stop.store(true, Ordering::SeqCst);
            }
            *self.slots[next]
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(result);
        });
    }

    /// The diagnostic lines of the module header, as text.
    pub(super) fn lines(&self) -> Vec<String> {
        let c = &self.counts;
        let mut lines = Vec::new();
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            lines.push(format!(
                "[X_STATS] TBXPIPE planned={} taken={} early={} loop={} waited={} key_mismatch={} max_held_hsps={}",
                self.plan.len(),
                c.taken.load(Ordering::Relaxed),
                c.early.load(Ordering::Relaxed),
                c.loop_searched.load(Ordering::Relaxed),
                c.waited.load(Ordering::Relaxed),
                c.key_mismatch.load(Ordering::Relaxed),
                self.held.most.load(Ordering::Relaxed),
            ));
        }
        lines
    }
}

#[cfg(all(test, feature = "parallel", not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    fn pool(threads: usize) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a test pool")
    }

    /// A fake batch search: `subjects` "subjects" on the pool, each ending its seed stage with
    /// `range.start + 1` HSPs, then a value made from the range; it records how many searches run
    /// at once.
    fn fake_search<'a>(
        running: &'a AtomicUsize,
        most: &'a AtomicUsize,
        subjects: usize,
        fail_at: Option<usize>,
    ) -> impl Fn(std::ops::Range<usize>, Option<&XSeedEnd<'_>>) -> Result<Vec<usize>> + Sync + 'a
    {
        move |range, seed_end| {
            let now = running.fetch_add(1, Ordering::SeqCst) + 1;
            most.fetch_max(now, Ordering::SeqCst);
            use rayon::prelude::*;
            (0..subjects).into_par_iter().for_each(|s| {
                let mut guard = XSeedGuard::new(seed_end);
                let mut x = 0u64;
                for k in 0..20_000u64 {
                    x = x.wrapping_mul(31).wrapping_add(k ^ s as u64);
                }
                std::hint::black_box(x);
                guard.end(range.start + 1);
                // the "chain" after the seed stage
                for k in 0..20_000u64 {
                    x = x.wrapping_mul(17).wrapping_add(k);
                }
                std::hint::black_box(x);
            });
            running.fetch_sub(1, Ordering::SeqCst);
            if fail_at == Some(range.start) {
                anyhow::bail!("batch {} failed", range.start);
            }
            Ok(range.collect())
        }
    }

    fn run_loop(
        threads: usize,
        batches: usize,
        limit: u64,
        fail_at: Option<usize>,
    ) -> (Vec<Result<Vec<usize>>>, usize, u64) {
        let plan: Vec<_> = (0..batches).map(|b| b * 2..b * 2 + 2).collect();
        let running = AtomicUsize::new(0);
        let most = AtomicUsize::new(0);
        let search = fake_search(&running, &most, 3, fail_at);
        let pipe = XPipe::with_limit(plan.clone(), 3, &search, limit).expect("two batches or more");
        let results = pool(threads).install(|| {
            rayon::in_place_scope(|scope| {
                let _closer = pipe.closer();
                let mut out = Vec::new();
                for range in plan.iter().cloned() {
                    let result = pipe.take(scope, range);
                    let failed = result.is_err();
                    out.push(result);
                    if failed {
                        break;
                    }
                }
                out
            })
        });
        (
            results,
            most.load(Ordering::SeqCst),
            pipe.counts.early.load(Ordering::SeqCst),
        )
    }

    #[test]
    fn the_loop_gets_every_batch_in_order_at_any_thread_count() {
        for threads in [1, 2, 4, 8, 16] {
            let (results, most, early) = run_loop(threads, 6, u64::MAX, None);
            let values: Vec<Vec<usize>> =
                results.into_iter().map(|r| r.expect("no error")).collect();
            let expected: Vec<Vec<usize>> = (0..6).map(|b| vec![b * 2, b * 2 + 1]).collect();
            assert_eq!(values, expected, "threads {threads}");
            assert!(most <= 6, "threads {threads}");
            assert_eq!(
                early, 5,
                "every batch after the first starts at the seed end before it"
            );
        }
    }

    #[test]
    fn a_limit_below_one_batch_keeps_one_search_at_a_time() {
        for threads in [2, 8] {
            let (results, most, early) = run_loop(threads, 4, 0, None);
            assert_eq!(results.len(), 4);
            assert!(results.iter().all(Result::is_ok));
            assert_eq!(early, 0);
            assert_eq!(most, 1, "threads {threads}");
        }
    }

    #[test]
    fn an_error_stops_the_early_starts_and_is_taken_in_order() {
        let (results, _, _) = run_loop(4, 6, u64::MAX, Some(4));
        assert_eq!(results.len(), 3);
        assert!(results[0].is_ok() && results[1].is_ok());
        assert_eq!(
            format!("{:#}", results[2].as_ref().unwrap_err()),
            "batch 4 failed"
        );
    }

    #[test]
    fn a_key_the_plan_does_not_hold_is_searched_on_the_loop() {
        let plan = vec![0..2, 2..4, 4..6];
        let running = AtomicUsize::new(0);
        let most = AtomicUsize::new(0);
        let search = fake_search(&running, &most, 2, None);
        let pipe = XPipe::with_limit(plan, 2, &search, u64::MAX).expect("three batches");
        let values = pool(4).install(|| {
            rayon::in_place_scope(|scope| {
                vec![
                    pipe.take(scope, 0..2).expect("planned"),
                    pipe.take(scope, 2..3).expect("not planned"),
                ]
            })
        });
        assert_eq!(values, vec![vec![0, 1], vec![2]]);
        assert_eq!(pipe.counts.key_mismatch.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn a_chain_links_on_the_pool_in_the_last_batch_or_over_the_limit() {
        let noop = |_: u64| {};
        let held = XHeld::with_limit(10);
        {
            let seed_end = XSeedEnd::new(3, &noop, &held, false);
            let mut a = XSeedGuard::new(Some(&seed_end));
            let mut b = XSeedGuard::new(Some(&seed_end));
            assert!(!a.end(3), "room left, a later batch may run");
            assert!(b.end(8), "11 HSPs held, over the limit");
            assert_eq!(held.now.load(Ordering::SeqCst), 11);
            drop(a);
            assert_eq!(held.now.load(Ordering::SeqCst), 8);
            // A subject that returns before its seed stage ends holds nothing.
            drop(XSeedGuard::new(Some(&seed_end)));
            drop(b);
        }
        assert_eq!(held.now.load(Ordering::SeqCst), 0);
        assert_eq!(held.most.load(Ordering::SeqCst), 11);
        let last = XSeedEnd::new(1, &noop, &held, true);
        assert!(XSeedGuard::new(Some(&last)).end(1), "the plan's last batch");
        assert!(XSeedGuard::new(None).end(1), "no pipe");
    }

    #[test]
    fn a_panic_in_an_early_batch_is_raised_on_the_loop() {
        for threads in [1, 4] {
            let plan = vec![0..2, 2..4, 4..6];
            let running = AtomicUsize::new(0);
            let most = AtomicUsize::new(0);
            let inner = fake_search(&running, &most, 3, None);
            let search = |range: std::ops::Range<usize>, seed_end: Option<&XSeedEnd<'_>>| {
                let start = range.start;
                let result = inner(range, seed_end);
                assert!(start != 2, "batch {start} panicked");
                result
            };
            let pipe =
                XPipe::with_limit(plan.clone(), 3, &search, u64::MAX).expect("three batches");
            let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                pool(threads).install(|| {
                    rayon::in_place_scope(|scope| {
                        let _closer = pipe.closer();
                        for range in plan.iter().cloned() {
                            let _ = pipe.take(scope, range);
                        }
                    })
                })
            }));
            let panic = caught.expect_err("the loop raises the batch's panic");
            let message = panic.downcast_ref::<String>().cloned().unwrap_or_default();
            assert!(
                message.contains("batch 2 panicked"),
                "threads {threads}: {message}"
            );
            assert_eq!(
                pipe.counts.taken.load(Ordering::SeqCst),
                2,
                "threads {threads}"
            );
            // Batch 1 may have started batch 2 before it panicked.
            assert!(
                pipe.counts.early.load(Ordering::SeqCst) >= 1,
                "threads {threads}"
            );
        }
    }

    #[test]
    fn after_the_loop_ends_no_batch_starts_early() {
        let plan = vec![0..2, 2..4, 4..6];
        let running = AtomicUsize::new(0);
        let most = AtomicUsize::new(0);
        let search = fake_search(&running, &most, 2, None);
        let pipe = XPipe::with_limit(plan, 2, &search, u64::MAX).expect("three batches");
        pool(2).install(|| {
            rayon::in_place_scope(|scope| {
                drop(pipe.closer());
                pipe.seed_ended(scope, 0);
            })
        });
        assert_eq!(pipe.counts.early.load(Ordering::SeqCst), 0);
        assert!(!pipe.started[1].load(Ordering::SeqCst));
    }

    #[test]
    fn a_plan_of_one_batch_makes_no_pipe() {
        let running = AtomicUsize::new(0);
        let most = AtomicUsize::new(0);
        let search = fake_search(&running, &most, 2, None);
        assert!(XPipe::with_limit(vec![0..1], 2, &search, u64::MAX).is_none());
    }
}
