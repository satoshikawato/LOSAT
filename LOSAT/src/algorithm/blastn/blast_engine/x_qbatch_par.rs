//! EXPERIMENT (LOSAT_X_QBATCHPAR, shadow LOSAT_X_QBATCHPARSHADOW): blastn's query batches
//! searched side by side, ahead of the batch loop of `run_in_pool`.
//!
//! NCBI's blastn reads its queries in batches (`CBatchSizeMixer`, blastn_app.cpp:261-300): the
//! size of each batch after the first comes from the successful initial extensions
//! (`good_init_extends`, "gie") of the batch before, so batch i+1 cannot be read before batch i
//! is searched. LOSAT's loop does the same and searches one batch at a time on the whole search
//! pool. A search of many queries against one long subject or a few has little parallel work
//! inside a batch (one subject chunk at a time), so most of the pool waits.
//!
//! With the switch on, the loop stays the only place that decides the batch boundaries, the
//! errors, the warnings, the appends and the order of the output. What changes is where a
//! batch's result comes from:
//! - a planner keeps a copy of the loop's mixer, advanced with *predicted* gie (the hit rate per
//!   letter of the last batch of the same size, times the planned batch's letters; 0 for a batch
//!   that will be split into query chunks, as the split batch reports), and plans the keys
//!   `(first record, end record, batch size)` of the next `2 k` batches;
//! - `k` slots search the planned batches ahead of the loop, each with `m` threads (`k x m <= n`):
//!   with `m = 1` on a plain thread with a one-thread pool (`SearchPool::x_serial`), with `m > 1`
//!   on a private pool of `m` threads whose first thread is the slot's (`with_search_pool`, built
//!   once per slot and search), so a batch sees only its own slot's pool;
//! - the loop computes each batch's key from the true mixer, as before, and takes a slot's result
//!   only under that same key (`Driver::take`); after the batch it reports the true gie
//!   (`Driver::commit`) and the planner plans again from the true mixer. A planned key that the
//!   new plan no longer holds is dropped (queued), cancelled (running) or discarded (finished).
//!
//! Why the output is the loop's: a batch's result is a function of its key (its records, its
//! batch size, which enters only an error message) and of the search's fixed inputs; no state
//! crosses batches (`-max_target_seqs`, culling and the hit lists are per query, and a query
//! never spans batches), thread-local scratch is reset per call, and the result does not depend
//! on the thread count. The first batch is never searched ahead (it is the only one that writes
//! to the outputs, `BatchOutputs`), and the loop takes a result only under the key it computed
//! itself, so it takes NCBI's batches in NCBI's order, with their errors at the same batch.
//! A slot's error is a returned value, handled by the loop's error paths when it takes the batch;
//! a slot's panic is kept and resumed only if the loop takes that batch (its message is printed
//! then, not when it happens). That needs unwinding: the release profile has `panic = "abort"`,
//! so there a panic in a batch a slot searches ends the process with its message. To bound that,
//! no slot runs while the loop writes output: `Driver::finish` cancels and joins every slot before
//! the reports and before each error path writes anything, and once a slot's batch returns an
//! error the batches planned after it are no longer started. The batches still exposed past an
//! error at batch j: before the error is known, the slots may have searched up to the plan's
//! window (2k - 1 batches) past j; once it is known, at most k - 1 still run, and they are
//! cancelled when the loop meets j. Slot panics are counted (`panicked=`), and in the shadow
//! mode any slot panic, taken or not, is a difference. `LOSAT_TIMING` blocks are text in
//! the batch's result, printed when the loop takes it.
//!
//! Cancellation: a slot's search gets a flag through its pool (`SearchPool::x_with_ext`,
//! `Cancel`), which `search_query_batch` reads at its start, between subjects, between subject
//! chunks of 5 Mb, between query chunks, and before its traceback; set, the search ends with
//! `Cancelled`, which is never taken. The loop's own searches have no flag. The flag is set for a
//! running batch that a new plan drops, when the speculation stops, and in `finish`, so the
//! wait for a cancelled batch is at most one subject chunk's search (or one query chunk's) and
//! its batch's set-up.
//!
//! Memory: each planned batch carries an estimate of its memory while it is searched
//! (`batch_estimate_bytes`), and a slot starts a batch only while the estimates of the batches
//! running plus the result sizes of the finished batches the loop has not taken yet stay within
//! the budget plus the loop's next batch, which is always admitted (the batch the loop needs).
//!
//! Waste: two mispredicted next keys within `GUARD_COMMITS` batches, or more batches discarded
//! (dropped finished results and cancelled searches) than `WASTE_SLOTS x k` plus a
//! `1 / WASTE_SHARE` share of the batches taken, stop the speculation; the loop then searches the
//! rest itself, on the whole pool.
//!
//! Active threads: the slots use `k x m <= n` threads. The loop's thread waits while a slot
//! searches its batch. The loop searches a batch on the whole pool only when no slot is running
//! and none can start (the first batch, and every batch once the search does not speculate).
//! `LOSAT_X_QBATCHPARSHADOW` searches every batch on the loop as well, on the whole pool, while
//! the slots run (a test mode; more than n threads are active then).
//!
//! Off (the loop as before): switch not set; a build without `parallel`, or wasm32; one thread;
//! a BLAST database (`-db`); `-verbose`, `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`,
//! `LOSAT_DEBUG_COORDS[_START]`, the blastn trace or `LOSAT_WASI_THREADS_DEBUG`; subjects of n or
//! more work units (5 Mb subject chunks; today's subject-parallel batch already uses the pool);
//! fewer than two batches after the first; a memory budget of fewer than two slots.
//!
//! Slots: per batch in flight `4^w x 4 + 8 x letters` bytes of lookup table (w: the batch's
//! lookup width, two contexts per letter), 13 bytes per query letter of query-side copies, and
//! 220 bytes per predicted initial extension for the hits (calibrated to the measured 50 MiB
//! per 1 Mb blastn batch and 116 MiB per 5 Mb dc-megablast batch). The budget for the batches in
//! flight beyond the loop's one is `LOSAT_X_QBATCHPARMEM` MiB, or by default 96 MiB (blastn,
//! blastn-short) or 28 MiB (megablast, dc-megablast) per thread beyond the first: 0.8 of NCBI
//! BLAST+ 2.17.0's measured growth of max RSS per thread on the audit input.
//! `k = min(1 + budget / estimate, n, R)` for R batch units of work left (the planned letters
//! over the largest planned batch), then lowered to `ceil(R / ceil(R / k))` (as many rounds,
//! fewer batches in flight); `m = n / k`. `LOSAT_X_QBATCHPARSLOTS=k` (tuning) sets k directly.
//!
//! Diagnostics: with the switch and `LOSAT_X_STATS` set, one line on standard error at the end of
//! the search, `[X_STATS] QBATCHPAR batches= precomputed= loop= mispredicted= discarded=
//! cancelled= panicked= waste_bound= k= m= idle= budget_mib= estimate_mib= peak_mib= off=` (precomputed:
//! batches taken from a slot; loop: batches the loop searched itself; discarded: finished or
//! cancelled batches dropped, `cancelled` of them cancelled; waste_bound: the discards that stop
//! the speculation, at the end; idle: threads no slot uses; peak_mib: the most memory admitted
//! at once; off: why the search did not, or stopped to, speculate, `none` if it never stopped).
//! Without `LOSAT_X_STATS` the switch writes nothing of its own. With the shadow,
//! `[X_SHADOW] QBATCHPAR batches= compared= diffs= panicked= mispredicted= wasted=`; the first difference
//! prints it and panics. `LOSAT_X_QBATCHPARMISS=N` (test only) corrupts every N-th prediction so
//! that the discard, replan and guard paths run.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-300
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!     ...
//!                 results = lcl_blast.Run();
//!                 if (!batch_size)
//!                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
//! ```

#![cfg_attr(
    not(all(feature = "parallel", not(target_arch = "wasm32"))),
    allow(dead_code, unused_imports, unused_variables)
)]

use std::collections::{HashMap, HashSet, VecDeque};
use std::ops::Range;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, OnceLock, PoisonError};

use crate::algorithm::blastn::query_split::{calculate_num_chunks, SplitSizes};
use crate::blastinput::query_batch::BatchSizeMixer;
use crate::blastinput::seq_range::{next_ranged_batch, QueryInput};
use crate::utils::threading::SearchPool;

/// Bytes in a MiB.
pub(super) const MIB: u64 = 1 << 20;
/// Query-side bytes per query letter of a batch in flight (the batch's copy of its records,
/// reverse complements, contexts with and without sentinels, masks).
const QUERY_BYTES_PER_LETTER: u64 = 13;
/// Hit bytes per predicted successful initial extension of a batch in flight.
const HIT_BYTES_PER_EXTENSION: u64 = 220;
/// Mispredicted next keys within `GUARD_COMMITS` batches that stop the speculation.
const GUARD_COMMITS: u64 = 8;
/// Discarded batches that stop the speculation: `WASTE_SLOTS x k` (a plan's window) plus one
/// `WASTE_SHARE`-th of the batches taken from the slots.
const WASTE_SLOTS: u64 = 2;
const WASTE_SHARE: u64 = 4;
/// The stack of a slot's thread: that of the main thread, where the loop's searches (and their
/// top level) ran before.
const SLOT_STACK_BYTES: usize = 8 << 20;

/// No NCBI counterpart: the switches of this path, read once per search; they choose a LOSAT
/// code path and change no value NCBI computes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Switches {
    /// `LOSAT_X_QBATCHPAR` (any value), or the shadow.
    pub on: bool,
    /// `LOSAT_X_QBATCHPARSHADOW` (any value; implies the switch).
    pub shadow: bool,
    /// `LOSAT_X_STATS` (any value): the summary line.
    pub stats: bool,
    /// `LOSAT_X_QBATCHPARMEM`: the budget of the batches in flight beyond the loop's, in MiB.
    pub mem_mib: Option<u64>,
    /// `LOSAT_X_QBATCHPARSLOTS` (tuning): the number of slots, at most the thread count.
    pub slots: Option<usize>,
    /// `LOSAT_X_QBATCHPARMISS` (test only): corrupt every N-th prediction.
    pub miss_every: Option<u64>,
}

impl Switches {
    /// The switches as the environment sets them. A value variable that is not a number (or is
    /// 0 for the slots and the corruption) counts as not set.
    pub(super) fn from_env() -> Self {
        let number = |name: &str| -> Option<u64> {
            std::env::var(name)
                .ok()
                .and_then(|value| value.trim().parse::<u64>().ok())
        };
        let shadow = std::env::var_os("LOSAT_X_QBATCHPARSHADOW").is_some();
        Switches {
            on: shadow || std::env::var_os("LOSAT_X_QBATCHPAR").is_some(),
            shadow,
            stats: std::env::var_os("LOSAT_X_STATS").is_some(),
            mem_mib: number("LOSAT_X_QBATCHPARMEM"),
            slots: number("LOSAT_X_QBATCHPARSLOTS")
                .filter(|&slots| slots > 0)
                .map(|slots| usize::try_from(slots).unwrap_or(usize::MAX)),
            miss_every: number("LOSAT_X_QBATCHPARMISS").filter(|&every| every > 0),
        }
    }
}

/// Why a search does not speculate, or stopped to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The switch is not set.
    NotSet,
    /// A build without the `parallel` feature, or wasm32.
    Build,
    /// One thread.
    Threads,
    /// A BLAST database (`-db`).
    Database,
    /// Verbose, debug or trace output.
    Debug,
    /// Subjects of as many work units as threads.
    Subjects,
    /// The search ended before its first batch was searched.
    Unstarted,
    /// Fewer than two batches after the first.
    Single,
    /// A budget of fewer than two slots.
    Budget,
    /// A slot's thread or pool could not be started.
    Spawn,
    /// Two mispredicted next keys within `GUARD_COMMITS` batches.
    Guard,
    /// More discarded batches than the waste bound.
    Waste,
    /// The loop's next key was not planned (not expected to happen).
    Unplanned,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Build => "build",
            Off::Threads => "threads",
            Off::Database => "db",
            Off::Debug => "debug",
            Off::Subjects => "subjects",
            Off::Unstarted => "unstarted",
            Off::Single => "single",
            Off::Budget => "budget",
            Off::Spawn => "spawn",
            Off::Guard => "guard",
            Off::Waste => "waste",
            Off::Unplanned => "unplanned",
        }
    }
}

/// The key of a query batch: its input records `input_start..input_end` and the batch size the
/// loop read them with (it enters only an error message). The loop computes it from the true
/// mixer; a result is used only under the key it was searched for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) struct BatchKey {
    pub input_start: usize,
    pub input_end: usize,
    pub batch_size: u32,
}

/// A batch for a slot: its key, its searched queries, and its estimated memory while it is
/// searched (`batch_estimate_bytes`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Job {
    pub key: BatchKey,
    pub searched: Range<usize>,
    pub bytes: u64,
}

/// The query input as the batch loop reads it.
#[derive(Clone, Copy)]
pub(super) struct PlanInput<'a> {
    pub query_input: &'a QueryInput,
    /// The length of each searched query.
    pub lengths: &'a [usize],
    /// The reading ends with a reader error after the last record (`QueryReading::Error`).
    pub reading_error: bool,
    /// `BATCH_SIZE` is set: every batch has that size.
    pub fixed_batch_size: bool,
    pub split: SplitSizes,
}

/// A planned batch.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Planned {
    pub key: BatchKey,
    pub searched: Range<usize>,
    /// Its searched letters, and its longest query.
    pub letters: u64,
    pub max_len: usize,
    /// The predicted successful initial extensions.
    pub gie: u64,
}

/// No NCBI counterpart: the speculation of NCBI's `CBatchSizeMixer` (scheduling only).
#[derive(Debug, Default)]
pub(super) struct Planner {
    /// The hit rate per letter of the last committed batch of each size.
    rates: HashMap<u32, f64>,
    /// That of the last committed batch.
    last_rate: Option<f64>,
    /// Predictions made (for `miss_every`).
    predictions: u64,
    miss_every: Option<u64>,
}

impl Planner {
    pub(super) fn new(miss_every: Option<u64>) -> Self {
        Planner {
            miss_every,
            ..Planner::default()
        }
    }

    /// A committed batch of `batch_size` with `letters` searched letters had `gie`.
    pub(super) fn record(&mut self, batch_size: u32, gie: u64, letters: u64) {
        let rate = gie as f64 / letters.max(1) as f64;
        self.rates.insert(batch_size, rate);
        self.last_rate = Some(rate);
    }

    /// The predicted gie of a batch of `batch_size` and `letters`; 0 for a batch that NCBI
    /// splits into query chunks (its own counter is not used; the chunks have theirs).
    fn predict(&self, batch_size: u32, letters: u64, split: bool) -> u64 {
        if split {
            return 0;
        }
        let rate = self
            .rates
            .get(&batch_size)
            .copied()
            .or(self.last_rate)
            .unwrap_or(0.0);
        (rate * letters as f64) as u64
    }

    /// Up to `window` batches from input record `start`, the first of `batch_size` (the true
    /// next size), the next ones of the sizes that `mixer` (the loop's mixer after the last
    /// batch; `None` with `BATCH_SIZE`) gives for the predicted gie. Stops after the batch that
    /// reaches the last record, and before a batch that ends the run with an error in the loop
    /// (only skipped records, only queries without letters, or the reader error after the last
    /// record). `speculate` counts the predictions (and corrupts every `miss_every`-th).
    pub(super) fn plan(
        &mut self,
        input: &PlanInput<'_>,
        start: usize,
        batch_size: u32,
        mixer: Option<&BatchSizeMixer>,
        window: usize,
        speculate: bool,
    ) -> Vec<Planned> {
        let records = input.query_input.input_lengths.len();
        let mut mixer = mixer.cloned();
        let mut planned = Vec::new();
        let (mut input_start, mut size) = (start, batch_size);
        while planned.len() < window && input_start < records {
            let (input_end, reached_size) = next_ranged_batch(
                &input.query_input.input_lengths,
                &input.query_input.skipped,
                input_start,
                size,
            );
            if input_end == records && !reached_size && input.reading_error {
                break;
            }
            let searched = input.query_input.searched_in(input_start..input_end);
            let lengths = &input.lengths[searched.clone()];
            if lengths.iter().all(|&length| length == 0) {
                break;
            }
            let letters: u64 = lengths.iter().map(|&length| length as u64).sum();
            let max_len = lengths.iter().copied().max().unwrap_or(0);
            let split = calculate_num_chunks(input.split, letters as usize).0 > 1;
            let gie = self.predict(size, letters, split);
            planned.push(Planned {
                key: BatchKey {
                    input_start,
                    input_end,
                    batch_size: size,
                },
                searched,
                letters,
                max_len,
                gie,
            });
            if input_end == records || planned.len() == window {
                break;
            }
            // The loop's `mixer.batch_size(Some(gie as u32 as i32))`.
            let mut next = match mixer.as_mut() {
                Some(mixer) if !input.fixed_batch_size => {
                    mixer.batch_size(Some(gie as u32 as i32)) as u32
                }
                _ => size,
            };
            if speculate {
                self.predictions += 1;
                if self
                    .miss_every
                    .is_some_and(|every| self.predictions % every == 0)
                {
                    next = if next > 1 { next / 2 } else { 2 };
                }
            }
            input_start = input_end;
            size = next;
        }
        planned
    }
}

/// The default budget per thread beyond the first, in MiB (0.8 x NCBI BLAST+ 2.17.0's max-RSS
/// growth per thread on the audit input: blastn (1888 - 127) / 15, dc-megablast (679 - 150) /
/// 15; megablast uses dc-megablast's).
pub(super) fn allowance_mib(task: &str) -> u64 {
    match task {
        "megablast" | "dc-megablast" => 28,
        _ => 96,
    }
}

/// The estimated memory of one batch in flight: lookup table of width `lut_width` over the two
/// contexts of `letters`, the query-side copies, and the hits of `gie` initial extensions.
pub(super) fn batch_estimate_bytes(lut_width: usize, letters: u64, gie: u64) -> u64 {
    let table = 4u64.saturating_pow(lut_width as u32).saturating_mul(4);
    table
        .saturating_add(letters.saturating_mul(2 * 4 + QUERY_BYTES_PER_LETTER))
        .saturating_add(gie.saturating_mul(HIT_BYTES_PER_EXTENSION))
}

/// The slots and the threads per slot for `threads` threads, a budget of `budget_bytes` for the
/// batches in flight beyond the loop's one, `estimate_bytes` per batch, `units` batch units of
/// work left, and the slot count of `LOSAT_X_QBATCHPARSLOTS`.
pub(super) fn choose_slots(
    threads: usize,
    budget_bytes: u64,
    estimate_bytes: u64,
    slots: Option<usize>,
    units: usize,
) -> Result<(usize, usize), Off> {
    if let Some(slots) = slots {
        let k = slots.clamp(1, threads.max(1));
        return Ok((k, threads / k));
    }
    let by_budget = 1 + budget_bytes / estimate_bytes.max(1);
    if by_budget < 2 {
        return Err(Off::Budget);
    }
    let k = (by_budget.min(threads as u64) as usize).min(units);
    if k < 2 {
        return Err(Off::Single);
    }
    let rounds = units.div_ceil(k);
    let k = units.div_ceil(rounds);
    Ok((k, threads / k))
}

/// The discarded batches beyond which the speculation stops, with `k` slots after `taken`
/// batches taken from them: a plan's window of waste, plus a quarter of the work taken.
pub(super) fn waste_bound(k: usize, taken: u64) -> u64 {
    WASTE_SLOTS * k as u64 + taken / WASTE_SHARE
}

/// What the executor needs to know of a slot's search result.
pub(super) trait SlotResult {
    /// The search ended with an error (not a cancellation): the loop stops at this batch if it
    /// takes it, so the batches planned after it are not started.
    fn failed(&self) -> bool;
    /// The search was cancelled (`Cancelled`); such a result is never taken.
    fn cancelled(&self) -> bool;
    /// The memory the result holds until the loop takes it.
    fn held_bytes(&self) -> u64;
}

/// No NCBI counterpart: the error of a slot's search that its cancellation ended (`Cancel`). It
/// is never taken by the loop.
#[derive(Debug)]
pub(super) struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a query batch searched ahead of the loop was cancelled")
    }
}

impl std::error::Error for Cancelled {}

/// No NCBI counterpart: the cancellation flag of a slot, attached to the pool its searches get
/// (`SearchPool::x_with_ext`); `search_query_batch` reads it (`cancelled`).
#[derive(Debug, Default)]
pub(super) struct Cancel(AtomicBool);

impl Cancel {
    pub(super) fn is_set(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }

    fn set(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    fn clear(&self) {
        self.0.store(false, Ordering::Relaxed);
    }
}

/// Whether the search on `pool` is a slot's whose batch was cancelled (never for the loop's own
/// searches, whose pool has no flag).
pub(super) fn cancelled(pool: &SearchPool<'_>) -> bool {
    pool.x_ext::<Cancel>().is_some_and(Cancel::is_set)
}

/// What a slot made of a batch: its value, or its panic with the messages the panic hook
/// kept instead of printing them.
pub(super) enum Outcome<T> {
    Value(T),
    Panic(Box<dyn std::any::Any + Send>, Vec<String>),
}

impl<T: SlotResult> Outcome<T> {
    /// The value; a panic is printed (its kept messages) and resumed on this thread.
    fn resume(self) -> T {
        match self {
            Outcome::Value(value) => value,
            Outcome::Panic(payload, messages) => {
                for message in &messages {
                    eprint!("{message}");
                }
                std::panic::resume_unwind(payload)
            }
        }
    }

    /// An error or a panic: the loop stops at this batch if it takes it.
    fn failed(&self) -> bool {
        match self {
            Outcome::Value(value) => value.failed(),
            Outcome::Panic(..) => true,
        }
    }
}

/// The batches of the slots.
struct State<T> {
    /// Planned batches that no slot has started, in batch order.
    queue: VecDeque<Job>,
    /// The batch each slot is searching, with its estimate.
    running: Vec<Option<(BatchKey, u64)>>,
    /// Finished batches of the current plan, with the memory their results hold.
    done: HashMap<BatchKey, (Outcome<T>, u64)>,
    /// The keys of the current plan.
    wanted: HashSet<BatchKey>,
    /// The loop's next key and its estimate: always admitted.
    next: Option<(BatchKey, u64)>,
    /// The memory the batches in flight may use beyond the loop's next batch.
    budget: u64,
    shutdown: bool,
    /// Finished or cancelled batches dropped unused, and how many of them were cancelled.
    discarded: u64,
    cancelled: u64,
    /// Slot batches that panicked (kept for the loop or dropped), and the first one's messages.
    panicked: u64,
    first_panic: Option<String>,
    /// The most slots searching at once, and the most memory admitted at once.
    peak_busy: usize,
    peak_bytes: u64,
    /// Slot threads that have not ended.
    alive: usize,
}

impl<T> State<T> {
    fn new() -> Self {
        State {
            queue: VecDeque::new(),
            running: Vec::new(),
            done: HashMap::new(),
            wanted: HashSet::new(),
            next: None,
            budget: 0,
            shutdown: false,
            discarded: 0,
            cancelled: 0,
            panicked: 0,
            first_panic: None,
            peak_busy: 0,
            peak_bytes: 0,
            alive: 0,
        }
    }

    /// The memory of the batches running and of the results not taken yet.
    fn in_use(&self) -> u64 {
        let running: u64 = self.running.iter().flatten().map(|(_, bytes)| bytes).sum();
        let held: u64 = self.done.values().map(|(_, bytes)| bytes).sum();
        running + held
    }

    /// Whether a slot may start `job` now: the loop's next batch always, any other while the
    /// memory in use with it stays within the budget plus the next batch's estimate.
    fn admits(&self, job: &Job) -> bool {
        match self.next {
            Some((key, _)) if key == job.key => true,
            next => self.in_use() + job.bytes <= self.budget + next.map_or(0, |(_, bytes)| bytes),
        }
    }

    fn note_peak(&mut self) {
        self.peak_bytes = self.peak_bytes.max(self.in_use());
        let busy = self.running.iter().flatten().count();
        self.peak_busy = self.peak_busy.max(busy);
    }
}

/// The counters of the summary lines.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Stats {
    pub batches: u64,
    pub precomputed: u64,
    pub loop_searched: u64,
    pub mispredicted: u64,
    pub k: usize,
    pub m: usize,
    pub idle: usize,
    pub budget_mib: u64,
    pub estimate_mib: u64,
    pub waste_bound: u64,
    pub off: Option<Off>,
    pub compared: u64,
    pub diffs: u64,
}

/// What the run's counters were at the end, with the executor's own.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Summary {
    pub stats: Stats,
    pub discarded: u64,
    pub cancelled: u64,
    pub panicked: u64,
    pub peak_busy: usize,
    pub peak_bytes: u64,
}

/// A slot's kept panic messages.
type PanicSink = Mutex<Vec<String>>;

/// What the loop and the slots share: declared before the thread scope of the loop, and the
/// summary lines at its end.
pub(super) struct Shared<T> {
    switches: Switches,
    state: Mutex<State<T>>,
    wake: Condvar,
    sinks: OnceLock<Vec<Arc<PanicSink>>>,
    cancels: OnceLock<Vec<Cancel>>,
    stats: Mutex<Stats>,
}

impl<T> Shared<T> {
    pub(super) fn new(switches: Switches) -> Self {
        Shared {
            switches,
            state: Mutex::new(State::new()),
            wake: Condvar::new(),
            sinks: OnceLock::new(),
            cancels: OnceLock::new(),
            stats: Mutex::new(Stats {
                off: Some(if switches.on {
                    Off::Unstarted
                } else {
                    Off::NotSet
                }),
                ..Stats::default()
            }),
        }
    }

    pub(super) fn switches(&self) -> Switches {
        self.switches
    }

    fn lock(&self) -> MutexGuard<'_, State<T>> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    fn wait<'g>(&self, guard: MutexGuard<'g, State<T>>) -> MutexGuard<'g, State<T>> {
        self.wake
            .wait(guard)
            .unwrap_or_else(PoisonError::into_inner)
    }

    /// Set the cancellation flag of every slot whose running batch `drop` selects.
    fn cancel_running(&self, state: &State<T>, drop: impl Fn(&BatchKey) -> bool) {
        if let Some(cancels) = self.cancels.get() {
            for (slot, running) in state.running.iter().enumerate() {
                if running.as_ref().is_some_and(|(key, _)| drop(key)) {
                    cancels[slot].set();
                }
            }
        }
    }

    /// The counters so far.
    pub(super) fn summary(&self) -> Summary {
        let stats = *self.stats.lock().unwrap_or_else(PoisonError::into_inner);
        let state = self.lock();
        Summary {
            stats,
            discarded: state.discarded,
            cancelled: state.cancelled,
            panicked: state.panicked,
            peak_busy: state.peak_busy,
            peak_bytes: state.peak_bytes,
        }
    }

    /// `[X_STATS] QBATCHPAR ...`
    pub(super) fn stats_line(&self) -> String {
        let summary = self.summary();
        let stats = summary.stats;
        format!(
            "[X_STATS] QBATCHPAR batches={} precomputed={} loop={} mispredicted={} discarded={} cancelled={} panicked={} waste_bound={} k={} m={} idle={} budget_mib={} estimate_mib={} peak_mib={} off={}",
            stats.batches,
            stats.precomputed,
            stats.loop_searched,
            stats.mispredicted,
            summary.discarded,
            summary.cancelled,
            summary.panicked,
            stats.waste_bound,
            stats.k,
            stats.m,
            stats.idle,
            stats.budget_mib,
            stats.estimate_mib,
            summary.peak_bytes.div_ceil(MIB),
            stats.off.map_or("none", Off::name)
        )
    }

    /// `[X_SHADOW] QBATCHPAR ...`
    pub(super) fn shadow_line(&self) -> String {
        let summary = self.summary();
        let stats = summary.stats;
        format!(
            "[X_SHADOW] QBATCHPAR batches={} compared={} diffs={} panicked={} mispredicted={} wasted={}",
            stats.batches,
            stats.compared,
            stats.diffs,
            summary.panicked,
            stats.mispredicted,
            summary.discarded
        )
    }
}

impl<T> Drop for Shared<T> {
    // No NCBI counterpart: the summary lines at the end of a search: the counters with
    // LOSAT_X_STATS (as the other counters of the tree), the shadow's line with the shadow.
    fn drop(&mut self) {
        if !self.switches.on || std::thread::panicking() {
            return;
        }
        if self.switches.stats {
            eprintln!("{}", self.stats_line());
        }
        if self.switches.shadow {
            eprintln!("{}", self.shadow_line());
        }
    }
}

/// Whether the search speculates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Phase {
    /// Before the first batch is committed.
    Pending,
    Active,
    Off(Off),
}

/// The loop's side: the planner, the plan, and the slots' threads (started in the thread scope
/// of the loop when the first batch is committed).
pub(super) struct Driver<'scope, 'env, T, F> {
    shared: &'env Shared<T>,
    search: &'env F,
    scope: &'scope std::thread::Scope<'scope, 'env>,
    input: PlanInput<'env>,
    threads: usize,
    allowance_mib: u64,
    /// The lookup width of a batch of these letters and longest query.
    lut_width: Box<dyn Fn(u64, usize) -> usize + 'env>,
    planner: Planner,
    phase: Phase,
    /// The keys of the current plan; the first is the next key the loop takes.
    plan: Vec<BatchKey>,
    window: usize,
    commits: u64,
    last_miss: Option<u64>,
    stats: Stats,
    slots: Vec<std::thread::ScopedJoinHandle<'scope, ()>>,
    finished: bool,
}

impl<'scope, 'env, T, F> Driver<'scope, 'env, T, F>
where
    T: SlotResult + Send + 'env,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync + 'env,
{
    /// The driver of one search; `off` is why the search cannot speculate at all (`None` when
    /// it may, once its first batch is committed).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        shared: &'env Shared<T>,
        search: &'env F,
        scope: &'scope std::thread::Scope<'scope, 'env>,
        input: PlanInput<'env>,
        threads: usize,
        allowance_mib: u64,
        lut_width: Box<dyn Fn(u64, usize) -> usize + 'env>,
        off: Option<Off>,
    ) -> Self {
        let switches = shared.switches;
        let off = if switches.on { off } else { Some(Off::NotSet) };
        let phase = match off {
            Some(off) => Phase::Off(off),
            None => Phase::Pending,
        };
        Driver {
            shared,
            search,
            scope,
            input,
            threads,
            allowance_mib,
            lut_width,
            planner: Planner::new(switches.miss_every),
            phase,
            plan: Vec::new(),
            window: 0,
            commits: 0,
            last_miss: None,
            stats: Stats {
                off: off.or(Some(Off::Unstarted)),
                ..Stats::default()
            },
            slots: Vec::new(),
            finished: false,
        }
    }

    /// Whether the loop also searches every batch it takes from a slot, and compares.
    pub(super) fn shadow(&self) -> bool {
        self.shared.switches.shadow && self.phase == Phase::Active
    }

    /// The result of batch `key` if a slot searched it (waiting for it if a slot is searching
    /// it or will); `None` when the loop searches it itself, on the whole pool: then no slot is
    /// running and none can start. The first batch is never searched ahead.
    pub(super) fn take(&mut self, key: BatchKey) -> Option<T> {
        if key.input_start != 0 && self.phase == Phase::Active {
            let mut state = self.shared.lock();
            let off = loop {
                if let Some((outcome, _)) = state.done.remove(&key) {
                    drop(state);
                    // Its memory is free for the next batches.
                    self.shared.wake.notify_all();
                    self.stats.precomputed += 1;
                    self.flush();
                    return Some(outcome.resume());
                }
                let planned = state
                    .running
                    .iter()
                    .flatten()
                    .any(|(running, _)| *running == key)
                    || state.queue.iter().any(|job| job.key == key);
                if !planned {
                    break Off::Unplanned;
                }
                // Every slot's thread ended outside a batch (after a panic it kept): nothing
                // will search the queued key.
                if state.alive == 0 {
                    break Off::Spawn;
                }
                state = self.shared.wait(state);
            };
            drop(state);
            self.stop(off);
        }
        self.wait_idle();
        self.stats.loop_searched += 1;
        self.flush();
        None
    }

    /// The shadow: the loop's own result of a batch taken from a slot equals the slot's
    /// (`compare` gives the first difference); the first difference prints the shadow line and
    /// panics.
    pub(super) fn shadow_check(
        &mut self,
        key: BatchKey,
        precomputed: &T,
        reference: &T,
        compare: impl Fn(&T, &T) -> Result<(), String>,
    ) {
        self.stats.compared += 1;
        if let Err(difference) = compare(precomputed, reference) {
            self.stats.diffs += 1;
            self.flush();
            eprintln!("{}", self.shared.shadow_line());
            panic!("[X_SHADOW] QBATCHPAR batch {key:?}: {difference}");
        }
        self.flush();
    }

    /// The shadow: a slot's batch that panicked, kept for the loop or dropped, is a difference;
    /// the first one prints the shadow line and panics with the batch's message (in the release
    /// profile, which aborts on a panic, the slot's panic has already ended the process).
    fn shadow_panics(&mut self) {
        if !self.shared.switches.shadow || self.stats.diffs > 0 {
            return;
        }
        let first = {
            let state = self.shared.lock();
            (state.panicked > 0).then(|| state.first_panic.clone().unwrap_or_default())
        };
        if let Some(message) = first {
            self.stats.diffs += 1;
            self.flush();
            eprintln!("{}", self.shared.shadow_line());
            panic!("[X_SHADOW] QBATCHPAR a batch searched by a slot panicked: {message}");
        }
    }

    /// Batch `key` is done in the loop with `gie` successful initial extensions; `mixer` is the
    /// loop's mixer after it (`None` with `BATCH_SIZE`) and `next_size` the next batch's size.
    pub(super) fn commit(
        &mut self,
        key: BatchKey,
        gie: u64,
        mixer: Option<&BatchSizeMixer>,
        next_size: u32,
    ) {
        self.stats.batches += 1;
        self.commits += 1;
        self.shadow_panics();
        let searched = self
            .input
            .query_input
            .searched_in(key.input_start..key.input_end);
        let letters = self.input.lengths[searched]
            .iter()
            .map(|&length| length as u64)
            .sum();
        self.planner.record(key.batch_size, gie, letters);
        if self.phase == Phase::Pending {
            self.activate(key, mixer, next_size);
        }
        if self.phase == Phase::Active {
            let next = self.planner.plan(
                &self.input,
                key.input_end,
                next_size,
                mixer,
                self.window,
                true,
            );
            let predicted = if self.plan.first() == Some(&key) {
                self.plan.get(1).copied()
            } else {
                None
            };
            let mut off = None;
            if let (Some(predicted), Some(actual)) = (predicted, next.first()) {
                if predicted != actual.key {
                    self.stats.mispredicted += 1;
                    let commit = self.commits;
                    if self
                        .last_miss
                        .is_some_and(|last| commit - last < GUARD_COMMITS)
                    {
                        off = Some(Off::Guard);
                    }
                    self.last_miss = Some(commit);
                }
            }
            match off {
                Some(off) => self.stop(off),
                None => self.replan(next),
            }
            // The waste bound: discarded (dropped or cancelled) batches.
            self.stats.waste_bound = waste_bound(self.stats.k, self.stats.precomputed);
            if self.phase == Phase::Active && self.shared.lock().discarded > self.stats.waste_bound
            {
                self.stop(Off::Waste);
            }
        }
        self.flush();
    }

    /// The end of the speculation before the loop writes output (the reports, or what an
    /// error writes): nothing more starts, running batches are cancelled, and every slot's
    /// thread is joined. Also run when the driver is dropped.
    pub(super) fn finish(&mut self) {
        if self.finished {
            return;
        }
        self.finished = true;
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            state.shutdown = true;
            state.queue.clear();
            state.wanted.clear();
            state.next = None;
            let unused: Vec<(Outcome<T>, u64)> = state.done.drain().map(|(_, done)| done).collect();
            state.discarded += unused.len() as u64;
            self.shared.cancel_running(&state, |_| true);
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
        let started = !self.slots.is_empty();
        for slot in self.slots.drain(..) {
            // A slot's thread catches the panics of its batches; a panic outside them (not
            // expected) has ended the thread, and the batches were dropped above.
            let _ = slot.join();
            UNJOINED.with(|unjoined| unjoined.set(unjoined.get().saturating_sub(1)));
        }
        if started && !std::thread::panicking() {
            check_no_global_pool();
            // The shadow: every panic of a slot's batch, taken or not, is a difference.
            self.shadow_panics();
        }
        self.flush();
    }

    fn flush(&self) {
        *self
            .shared
            .stats
            .lock()
            .unwrap_or_else(PoisonError::into_inner) = self.stats;
    }

    fn set_off(&mut self, off: Off) {
        self.phase = Phase::Off(off);
        self.stats.off = Some(off);
    }

    /// After the first batch: plan the whole input once, choose the slots, start them.
    fn activate(&mut self, key: BatchKey, mixer: Option<&BatchSizeMixer>, next_size: u32) {
        let full = self.planner.plan(
            &self.input,
            key.input_end,
            next_size,
            mixer,
            usize::MAX,
            false,
        );
        if full.len() < 2 {
            self.set_off(Off::Single);
            return;
        }
        let total: u64 = full.iter().map(|planned| planned.letters).sum();
        let largest = full
            .iter()
            .map(|planned| planned.letters)
            .max()
            .unwrap_or(1)
            .max(1);
        let units = usize::try_from(total.div_ceil(largest)).unwrap_or(usize::MAX);
        let estimate = full
            .iter()
            .map(|planned| self.estimate(planned))
            .max()
            .unwrap_or(0);
        let switches = self.shared.switches;
        let budget_mib = switches
            .mem_mib
            .unwrap_or(self.allowance_mib * (self.threads as u64).saturating_sub(1));
        self.stats.budget_mib = budget_mib;
        self.stats.estimate_mib = estimate.div_ceil(MIB);
        let (k, m) = match choose_slots(
            self.threads,
            budget_mib.saturating_mul(MIB),
            estimate,
            switches.slots,
            units,
        ) {
            Ok(slots) => slots,
            Err(off) => {
                self.set_off(off);
                return;
            }
        };
        self.stats.k = k;
        self.stats.m = m;
        self.stats.idle = self.threads.saturating_sub(k * m);
        self.shared.lock().budget = budget_mib.saturating_mul(MIB);
        match self.start_slots(k, m) {
            Ok(()) => {
                self.phase = Phase::Active;
                self.stats.off = None;
                self.window = 2 * k;
            }
            Err(off) => self.stop(off),
        }
    }

    /// The memory estimate of a planned batch.
    fn estimate(&self, planned: &Planned) -> u64 {
        batch_estimate_bytes(
            (self.lut_width)(planned.letters, planned.max_len),
            planned.letters,
            planned.gie,
        )
    }

    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    fn start_slots(&mut self, _k: usize, _m: usize) -> Result<(), Off> {
        Err(Off::Build)
    }

    /// `k` slot threads (named `losat-qbp-<slot>`), each with a pool of `m` threads when `m > 1`
    /// (`slot_main`).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    fn start_slots(&mut self, k: usize, m: usize) -> Result<(), Off> {
        install_panic_hook();
        let sinks: Vec<Arc<PanicSink>> = (0..k).map(|_| Arc::default()).collect();
        let cancels: Vec<Cancel> = (0..k).map(|_| Cancel::default()).collect();
        if self.shared.sinks.set(sinks).is_err() || self.shared.cancels.set(cancels).is_err() {
            return Err(Off::Spawn);
        }
        self.shared.lock().running = vec![None; k];
        for slot in 0..k {
            let (shared, search) = (self.shared, self.search);
            self.shared.lock().alive += 1;
            let spawned = std::thread::Builder::new()
                .name(format!("losat-qbp-{slot}"))
                .stack_size(SLOT_STACK_BYTES)
                .spawn_scoped(self.scope, move || slot_main(shared, search, slot, m));
            match spawned {
                Ok(handle) => {
                    self.slots.push(handle);
                    UNJOINED.with(|unjoined| unjoined.set(unjoined.get() + 1));
                }
                Err(_) => {
                    self.shared.lock().alive -= 1;
                    return Err(Off::Spawn);
                }
            }
        }
        Ok(())
    }

    /// Stop speculating: nothing more is started, running batches are cancelled, finished
    /// results are dropped.
    fn stop(&mut self, off: Off) {
        self.set_off(off);
        self.plan.clear();
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            state.queue.clear();
            state.wanted.clear();
            state.next = None;
            let unused: Vec<(Outcome<T>, u64)> = state.done.drain().map(|(_, done)| done).collect();
            state.discarded += unused.len() as u64;
            self.shared.cancel_running(&state, |_| true);
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
    }

    /// The plan `next` replaces the current one: its keys that no slot has started or finished
    /// are queued in order (up to a finished batch that failed: the loop stops there if it
    /// reaches it); finished batches outside it are dropped and running ones cancelled.
    fn replan(&mut self, next: Vec<Planned>) {
        let jobs: Vec<Job> = next
            .iter()
            .map(|planned| Job {
                key: planned.key,
                searched: planned.searched.clone(),
                bytes: self.estimate(planned),
            })
            .collect();
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            let mut jobs = jobs;
            if let Some(failed) = jobs.iter().position(|job| {
                state
                    .done
                    .get(&job.key)
                    .is_some_and(|(outcome, _)| outcome.failed())
            }) {
                jobs.truncate(failed + 1);
            }
            self.plan = jobs.iter().map(|job| job.key).collect();
            let wanted: HashSet<BatchKey> = self.plan.iter().copied().collect();
            let stale: Vec<BatchKey> = state
                .done
                .keys()
                .filter(|key| !wanted.contains(key))
                .copied()
                .collect();
            let unused: Vec<(Outcome<T>, u64)> = stale
                .iter()
                .filter_map(|key| state.done.remove(key))
                .collect();
            state.discarded += unused.len() as u64;
            self.shared
                .cancel_running(&state, |key| !wanted.contains(key));
            state.next = jobs.first().map(|job| (job.key, job.bytes));
            let queue: VecDeque<Job> = jobs
                .into_iter()
                .filter(|job| {
                    !state
                        .running
                        .iter()
                        .flatten()
                        .any(|(running, _)| *running == job.key)
                        && !state.done.contains_key(&job.key)
                })
                .collect();
            state.queue = queue;
            state.wanted = wanted;
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
    }

    /// Wait until no slot is searching (none starts: nothing is queued when this is called).
    fn wait_idle(&self) {
        let mut state = self.shared.lock();
        while state.running.iter().any(Option::is_some) {
            state = self.shared.wait(state);
        }
    }
}

impl<T, F> Drop for Driver<'_, '_, T, F> {
    // The end of the loop without `finish` (an early return or a panic): the same steps.
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        self.finished = true;
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            state.shutdown = true;
            state.queue.clear();
            state.wanted.clear();
            state.next = None;
            let unused: Vec<(Outcome<T>, u64)> = state.done.drain().map(|(_, done)| done).collect();
            state.discarded += unused.len() as u64;
            self.shared.cancel_running(&state, |_| true);
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
        for slot in self.slots.drain(..) {
            let _ = slot.join();
            UNJOINED.with(|unjoined| unjoined.set(unjoined.get().saturating_sub(1)));
        }
        *self
            .shared
            .stats
            .lock()
            .unwrap_or_else(PoisonError::into_inner) = self.stats;
    }
}

/// A slot's thread: search queued batches until the loop ends. With `threads > 1` the slot's
/// thread is the first thread of a pool of `threads` (`with_search_pool`, as the search's own
/// pool has the loop's thread), so a batch's search runs on the slot's thread with that pool
/// around it for its parallel stages (one heap arena keeps a batch's large buffers, as the loop
/// keeps them in its own), and the pool's threads are joined when the slot ends. A pool that
/// cannot be built leaves the slot searching on its own thread. Either pool carries the slot's
/// cancellation flag.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn slot_main<T, F>(shared: &Shared<T>, search: &F, slot: usize, threads: usize)
where
    T: SlotResult + Send,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync,
{
    let sink = Arc::clone(
        &shared
            .sinks
            .get()
            .expect("slot sinks are set before the slots")[slot],
    );
    let cancel = &shared
        .cancels
        .get()
        .expect("slot flags are set before the slots")[slot];
    set_panic_sink(Some(Arc::clone(&sink)));
    // When the slot's thread ends (also by a panic outside a batch, not expected), the loop must
    // not wait for it.
    struct Ended<'a, T>(&'a Shared<T>, usize);
    impl<T> Drop for Ended<'_, T> {
        fn drop(&mut self) {
            let mut state = self.0.lock();
            state.alive -= 1;
            if let Some(running) = state.running.get_mut(self.1) {
                *running = None;
            }
            drop(state);
            self.0.wake.notify_all();
        }
    }
    let _ended = Ended(shared, slot);
    if threads > 1 {
        let served = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
            rayon::broadcast(|_| set_panic_sink(Some(Arc::clone(&sink))));
            slot_loop(
                shared,
                search,
                slot,
                &pool.x_with_ext(cancel),
                &sink,
                cancel,
            );
            Ok(())
        });
        if served.is_ok() {
            return;
        }
    }
    let serial = SearchPool::x_serial();
    slot_loop(
        shared,
        search,
        slot,
        &serial.x_with_ext(cancel),
        &sink,
        cancel,
    );
}

/// The batches of one slot, each searched with `pool`, until the loop ends or a batch panics
/// (its thread-local state is then not trusted further; the panic is kept for the loop).
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn slot_loop<T, F>(
    shared: &Shared<T>,
    search: &F,
    slot: usize,
    pool: &SearchPool<'_>,
    sink: &PanicSink,
    cancel: &Cancel,
) where
    T: SlotResult + Send,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync,
{
    loop {
        let job = {
            let mut state = shared.lock();
            loop {
                if state.shutdown {
                    return;
                }
                if state.queue.front().is_some_and(|job| state.admits(job)) {
                    let job = state.queue.pop_front().expect("a queued batch");
                    state.running[slot] = Some((job.key, job.bytes));
                    cancel.clear();
                    state.note_peak();
                    break job;
                }
                state = shared.wait(state);
            }
        };
        sink.lock().unwrap_or_else(PoisonError::into_inner).clear();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| search(&job, pool)));
        let panicked = result.is_err();
        let outcome = match result {
            Ok(value) => Outcome::Value(value),
            Err(payload) => Outcome::Panic(
                payload,
                std::mem::take(&mut *sink.lock().unwrap_or_else(PoisonError::into_inner)),
            ),
        };
        // The result's size, before the lock (it walks every HSP).
        let held = match &outcome {
            Outcome::Value(value) => value.held_bytes(),
            Outcome::Panic(..) => 0,
        };
        let unused = {
            let mut state = shared.lock();
            // The flag as it stands when the search has returned (it is set and cleared only
            // under this lock).
            let flagged = cancel.is_set();
            settle(&mut state, slot, job, outcome, held, flagged)
        };
        shared.wake.notify_all();
        drop(unused);
        if panicked {
            return;
        }
    }
}

/// What becomes of a slot's finished batch, under the state lock; returns what to drop outside
/// it. `held` is the result's size and `flagged` whether the slot's cancellation flag was set
/// when the search returned. A flagged search counts as cancelled whatever it returned (a search
/// may stop part-way at any of its checks, and a later edit could return a partial value
/// there), so it is never taken: it is counted, and searched again if the plan wants its key.
/// A panic is counted whether it is kept or dropped.
fn settle<T: SlotResult>(
    state: &mut State<T>,
    slot: usize,
    job: Job,
    outcome: Outcome<T>,
    held: u64,
    flagged: bool,
) -> Option<Outcome<T>> {
    state.running[slot] = None;
    if let Outcome::Panic(_, messages) = &outcome {
        state.panicked += 1;
        if state.first_panic.is_none() {
            state.first_panic = Some(messages.concat());
        }
    }
    let wanted =
        !state.shutdown && state.wanted.contains(&job.key) && !state.done.contains_key(&job.key);
    let cancelled = flagged || matches!(&outcome, Outcome::Value(value) if value.cancelled());
    if cancelled {
        state.discarded += 1;
        state.cancelled += 1;
        if wanted {
            let at = state
                .queue
                .iter()
                .position(|queued| queued.key.input_start > job.key.input_start)
                .unwrap_or(state.queue.len());
            state.queue.insert(at, job);
        }
        return Some(outcome);
    }
    if !wanted {
        state.discarded += 1;
        return Some(outcome);
    }
    if outcome.failed() {
        // The loop stops at this batch if it reaches it: nothing after it starts.
        state
            .queue
            .retain(|queued| queued.key.input_start < job.key.input_start);
    }
    state.done.insert(job.key, (outcome, held));
    state.note_peak();
    None
}

thread_local! {
    // No NCBI counterpart: the slot threads of the search on this thread (the loop's) that have
    // not been joined; the writers of `run_in_pool` check that there are none.
    static UNJOINED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// The writers of `run_in_pool` (the reports, and what a failed batch writes) run only after
/// `Driver::finish` has joined every slot of the search, so that no speculative batch runs
/// while output is written. With debug assertions this checks it, on the loop's thread.
pub(super) fn debug_assert_slots_joined() {
    debug_assert_eq!(
        UNJOINED.with(std::cell::Cell::get),
        0,
        "a writer of run_in_pool runs while LOSAT_X_QBATCHPAR slots may still run (Driver::finish was not called before it)"
    );
}

thread_local! {
    // No NCBI counterpart: where the panic hook keeps the message of a panic on a slot's
    // threads (diagnostics only).
    static PANIC_SINK: std::cell::RefCell<Option<Arc<PanicSink>>> =
        const { std::cell::RefCell::new(None) };
}

fn set_panic_sink(sink: Option<Arc<PanicSink>>) {
    PANIC_SINK.with(|current| *current.borrow_mut() = sink);
}

/// The panic hook of the process, once: on a slot's threads it keeps the message (printed if
/// the loop takes the batch); elsewhere it calls the hook it replaced. Not with `panic =
/// "abort"` (the release profile): a panic cannot be kept there, and the process ends with the
/// default hook's message.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn install_panic_hook() {
    if cfg!(panic = "abort") {
        return;
    }
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            let sink = PANIC_SINK
                .try_with(|current| current.borrow().clone())
                .ok()
                .flatten();
            let Some(sink) = sink else {
                previous(info);
                return;
            };
            let payload = info
                .payload()
                .downcast_ref::<&str>()
                .map(|text| (*text).to_string())
                .or_else(|| info.payload().downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "Box<dyn Any>".to_string());
            let location = info
                .location()
                .map(|location| {
                    format!(
                        "{}:{}:{}",
                        location.file(),
                        location.line(),
                        location.column()
                    )
                })
                .unwrap_or_default();
            let thread = std::thread::current();
            let message = format!(
                "thread '{}' panicked at {location}:\n{payload}\n",
                thread.name().unwrap_or("<unnamed>")
            );
            sink.lock()
                .unwrap_or_else(PoisonError::into_inner)
                .push(message);
        }));
    });
}

/// Design section 7.4: no slot reaches rayon's global pool. A slot of one thread runs on a
/// plain thread, where a `par_iter` or `rayon::join` without a check of
/// `rayon::current_thread_index()` would start the global pool with every CPU; the slots of
/// `m > 1` threads have their own pool. With debug assertions (and not in this crate's unit
/// tests, whose other tests may use the global pool), after the slots are joined: the global
/// pool must not exist yet, unless this check made it (a one-thread pool, so a later search's
/// slots are checked as well as the first one's are only if nothing else started it).
/// Caveats, all in debug builds only: `build_global` can also fail for another reason (a thread
/// that cannot be spawned), and a host that embeds LOSAT and has built the global pool itself
/// trips the check; the check leaves a one-thread global pool in the process. Fine for the
/// command line, whose blastn path never builds the global pool.
fn check_no_global_pool() {
    #[cfg(all(
        debug_assertions,
        not(test),
        feature = "parallel",
        not(target_arch = "wasm32")
    ))]
    {
        static OURS: AtomicBool = AtomicBool::new(false);
        let built = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .thread_name(|_| "losat-qbp-global-check".to_string())
            .build_global()
            .is_ok();
        if built {
            OURS.store(true, Ordering::Relaxed);
        }
        debug_assert!(
            built || OURS.load(Ordering::Relaxed),
            "a LOSAT_X_QBATCHPAR slot started rayon's global pool (a parallel iterator or join without a check of rayon::current_thread_index())"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::seq_range::Placements;

    fn input_of(lengths: &[usize], skipped: &[bool]) -> QueryInput {
        let ordinals = (0..lengths.len())
            .filter(|&index| !skipped[index])
            .collect();
        QueryInput {
            placements: Placements::default(),
            ordinals,
            input_lengths: lengths.to_vec(),
            skipped: skipped.to_vec(),
            range_length: None,
        }
    }

    fn searched_lengths(input: &QueryInput) -> Vec<usize> {
        input
            .ordinals
            .iter()
            .map(|&ordinal| input.input_lengths[ordinal])
            .collect()
    }

    fn plan_input<'a>(
        query_input: &'a QueryInput,
        lengths: &'a [usize],
        fixed: bool,
        megablast: bool,
    ) -> PlanInput<'a> {
        PlanInput {
            query_input,
            lengths,
            reading_error: false,
            fixed_batch_size: fixed,
            split: SplitSizes::default_for(megablast),
        }
    }

    /// One batch of the loop: its key, letters and gie, the mixer after it, and the next size.
    type TrueBatch = (BatchKey, u64, u64, BatchSizeMixer, u32);

    /// The loop's batches for `gie` (a function of a batch's searched letters), as `run_in_pool`
    /// reads them.
    fn true_batches(
        input: &PlanInput<'_>,
        first: u32,
        max: i32,
        gie: impl Fn(u64) -> u64,
    ) -> Vec<TrueBatch> {
        let mut mixer = BatchSizeMixer::new(max);
        assert_eq!(mixer.batch_size(None) as u32, first);
        let mut size = first;
        let mut start = 0;
        let mut out = Vec::new();
        let records = input.query_input.input_lengths.len();
        while start < records {
            let (end, _) = next_ranged_batch(
                &input.query_input.input_lengths,
                &input.query_input.skipped,
                start,
                size,
            );
            let searched = input.query_input.searched_in(start..end);
            let letters: u64 = input.lengths[searched].iter().map(|&l| l as u64).sum();
            let extensions = gie(letters);
            let key = BatchKey {
                input_start: start,
                input_end: end,
                batch_size: size,
            };
            size = mixer.batch_size(Some(extensions as u32 as i32)) as u32;
            out.push((key, letters, extensions, mixer.clone(), size));
            start = end;
        }
        out
    }

    /// Replays the loop with the planner: how many next keys the plan mispredicted.
    fn mispredictions(input: &PlanInput<'_>, batches: &[TrueBatch]) -> usize {
        let mut planner = Planner::new(None);
        let mut plan: Vec<BatchKey> = Vec::new();
        let mut misses = 0;
        for (index, (key, letters, gie, mixer, next)) in batches.iter().enumerate() {
            planner.record(key.batch_size, *gie, *letters);
            let next_plan = planner.plan(input, key.input_end, *next, Some(mixer), 8, true);
            if index > 0 && plan.first() == Some(key) {
                if let (Some(predicted), Some(actual)) = (plan.get(1), next_plan.first()) {
                    if predicted != &actual.key {
                        misses += 1;
                    }
                }
            }
            plan = next_plan.iter().map(|planned| planned.key).collect();
            // The plan's first key is the loop's next key.
            if let Some((next_key, ..)) = batches.get(index + 1) {
                assert_eq!(plan.first(), Some(next_key), "batch {index}");
            }
        }
        misses
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:65-81 (`GetBatchSize`).
    // dc-megablast on the audit input: a batch with few extensions per letter (the first one of
    // 5,000, or a full one of 4,999,000 with fewer than ~2,328) makes `(Int4)` of the next size
    // overflow to INT_MIN, so the next batch is the smallest (100, one record), and a small
    // batch is followed by a full one. The scaled rate of the first batch predicts it.
    #[test]
    fn the_planner_follows_the_mixer_alternation_of_tiny_and_full_batches() {
        let lengths: Vec<usize> = (0..12_000)
            .map(|index| 1_500 + (index * 37) % 2_000)
            .collect();
        let query_input = input_of(&lengths, &vec![false; lengths.len()]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, false, true);
        let max = SplitSizes::default_for(true).mixer_max_batch_size();
        let batches = true_batches(&input, 5_000, max, |letters| letters / 4_000);
        let sizes: Vec<u32> = batches.iter().map(|(key, ..)| key.batch_size).collect();
        assert_eq!(&sizes[..5], &[5_000, 100, 4_999_000, 100, 4_999_000]);
        assert!(sizes.len() > 5, "{sizes:?}");
        assert_eq!(mispredictions(&input, &batches), 0);
    }

    // blastn (word 11) on the audit input: every batch after the first has the maximum size
    // (999,000) and its gie stays in the clamped regime.
    #[test]
    fn the_planner_predicts_blastn_batches_without_misses() {
        let lengths: Vec<usize> = (0..3_000)
            .map(|index| 2_000 + (index * 53) % 3_000)
            .collect();
        let query_input = input_of(&lengths, &vec![false; lengths.len()]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, false, false);
        let max = SplitSizes::default_for(false).mixer_max_batch_size();
        let batches = true_batches(&input, 5_000, max, |letters| letters * 72 / 1_000);
        assert!(batches.len() > 6);
        assert!(batches[1..batches.len() - 1]
            .iter()
            .all(|(key, ..)| key.batch_size == 999_000));
        assert_eq!(mispredictions(&input, &batches), 0);
    }

    // `BATCH_SIZE`: every batch has that size, whatever its gie; no prediction can be wrong.
    #[test]
    fn a_fixed_batch_size_plans_exact_keys() {
        let lengths = [300, 300, 300, 300, 300, 300, 300];
        let query_input = input_of(&lengths, &[false; 7]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(None);
        planner.record(500, 10, 600);
        let plan = planner.plan(&input, 2, 500, None, 10, true);
        let keys: Vec<(usize, usize)> = plan
            .iter()
            .map(|planned| (planned.key.input_start, planned.key.input_end))
            .collect();
        assert_eq!(keys, [(2, 4), (4, 6), (6, 7)]);
        assert!(plan.iter().all(|planned| planned.key.batch_size == 500));
        // The window bounds the plan.
        assert_eq!(planner.plan(&input, 2, 500, None, 2, true).len(), 2);
    }

    // The plan stops before a batch that ends the run in the loop. A batch reads records until
    // their letters reach its size (skipped records and records without letters add none), so
    // such a batch is the last one.
    #[test]
    fn the_plan_stops_before_a_batch_that_ends_the_run() {
        // Records 3 and 4 are skipped (`-query_loc` past their end): the batch from record 3
        // has no searched query (`Empty CBlastQueryVector`).
        let lengths = [100, 100, 100, 50, 50];
        let skipped = [false, false, false, true, true];
        let query_input = input_of(&lengths, &skipped);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(None);
        let plan = planner.plan(&input, 1, 1, None, 10, true);
        assert_eq!(plan.len(), 2);
        assert_eq!(plan[1].key.input_end, 3);
        // A skipped record inside the input joins the next batch.
        let skipped = [false, true, false, false, false];
        let query_input = input_of(&lengths, &skipped);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let plan = planner.plan(&input, 1, 1, None, 10, true);
        assert_eq!(plan[0].key.input_end, 3);
        assert_eq!(plan[0].searched, 1..2);
        // Queries without letters only (`Sequence contains no data`).
        let lengths = [100, 100, 0, 0];
        let query_input = input_of(&lengths, &[false; 4]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let plan = planner.plan(&input, 1, 1, None, 10, true);
        assert_eq!(plan.len(), 1);
        // A batch with a query with letters is searched.
        let plan = planner.plan(&input, 1, 150, None, 10, true);
        assert_eq!(plan[0].key.input_end, 4);
        // The reader error after the last record comes in the batch that does not reach its
        // size; a batch that reaches it is searched, and the error comes after it.
        let lengths = [100, 100, 100];
        let query_input = input_of(&lengths, &[false; 3]);
        let searched = searched_lengths(&query_input);
        let mut input = plan_input(&query_input, &searched, true, false);
        input.reading_error = true;
        let plan = planner.plan(&input, 0, 150, None, 10, true);
        assert_eq!(plan.len(), 1, "{plan:?}");
        let plan = planner.plan(&input, 0, 100, None, 10, true);
        assert_eq!(plan.len(), 3, "{plan:?}");
    }

    // A batch that NCBI splits into query chunks reports no extensions of its own.
    #[test]
    fn a_split_batch_is_predicted_without_extensions() {
        let lengths = [2_500_000, 100, 100];
        let query_input = input_of(&lengths, &[false; 3]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(None);
        planner.record(1, 1_000, 1_000);
        let plan = planner.plan(&input, 0, 1, None, 3, true);
        assert_eq!(plan[0].gie, 0);
        assert_eq!(plan[1].gie, 100);
    }

    // `LOSAT_X_QBATCHPARMISS=N`: every N-th prediction gives another size.
    #[test]
    fn the_test_switch_corrupts_every_nth_prediction() {
        let lengths = [100; 12];
        let query_input = input_of(&lengths, &[false; 12]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(Some(2));
        let plan = planner.plan(&input, 0, 200, None, 4, true);
        let sizes: Vec<u32> = plan.iter().map(|planned| planned.key.batch_size).collect();
        assert_eq!(sizes, [200, 200, 100, 100]);
        // Without counting (the whole-input plan of `activate`) nothing is corrupted.
        let plan = planner.plan(&input, 0, 200, None, 4, false);
        assert!(plan.iter().all(|planned| planned.key.batch_size == 200));
    }

    // The estimate is calibrated to the measured memory of one batch in flight: ~50 MiB for a
    // 1 Mb blastn batch (lookup width 11, ~72,000 extensions) and ~116 MiB for a 5 Mb
    // dc-megablast batch (width 11, ~1,000 extensions).
    #[test]
    fn the_estimate_matches_the_measured_batches() {
        let blastn = batch_estimate_bytes(11, 999_193, 71_763) as f64 / MIB as f64;
        let dc = batch_estimate_bytes(11, 5_001_442, 1_000) as f64 / MIB as f64;
        assert!((47.0..=54.0).contains(&blastn), "{blastn}");
        assert!((110.0..=122.0).contains(&dc), "{dc}");
    }

    // Decision M-30 (4): blastn at 16 threads has 15 x 96 MiB, so k = 16 by the budget and 15
    // after round balancing of ~86 units; dc-megablast has 15 x 28 = 420 MiB: 4 slots x 4
    // threads.
    #[test]
    fn slots_follow_the_budget_and_the_rounds() {
        let blastn = batch_estimate_bytes(11, 1_000_000, 72_000);
        assert_eq!(
            choose_slots(16, 15 * 96 * MIB, blastn, None, 86),
            Ok((15, 1))
        );
        let dc = batch_estimate_bytes(11, 5_000_000, 1_000);
        assert_eq!(choose_slots(16, 15 * 28 * MIB, dc, None, 18), Ok((4, 4)));
        // A large budget: 18 dc units take two rounds at 16 slots, so 9 slots.
        assert_eq!(choose_slots(16, 4096 * MIB, dc, None, 18), Ok((9, 1)));
        // The tuning variable wins, within the thread count.
        assert_eq!(choose_slots(16, 0, dc, Some(8), 18), Ok((8, 2)));
        assert_eq!(choose_slots(4, 0, dc, Some(8), 18), Ok((4, 1)));
        // Too small a budget, or one unit of work left.
        assert_eq!(choose_slots(16, 50 * MIB, dc, None, 18), Err(Off::Budget));
        assert_eq!(choose_slots(16, 4096 * MIB, dc, None, 1), Err(Off::Single));
        assert_eq!(allowance_mib("blastn"), 96);
        assert_eq!(allowance_mib("blastn-short"), 96);
        assert_eq!(allowance_mib("dc-megablast"), 28);
        assert_eq!(allowance_mib("megablast"), 28);
    }

    // The waste bound: a plan's window (two batches per slot) plus a quarter of the batches
    // taken.
    #[test]
    fn the_waste_bound_grows_with_the_slots_and_the_work_taken() {
        assert_eq!(waste_bound(2, 0), 4);
        assert_eq!(waste_bound(4, 7), 9);
        assert_eq!(waste_bound(15, 100), 55);
    }

    /// A result for `settle`: `Err("cancelled")` is a cancelled search, another `Err` a failed
    /// one.
    type Settled = Result<u32, &'static str>;

    impl SlotResult for Settled {
        fn failed(&self) -> bool {
            matches!(self, Err(error) if *error != "cancelled")
        }

        fn cancelled(&self) -> bool {
            matches!(self, Err(error) if *error == "cancelled")
        }

        fn held_bytes(&self) -> u64 {
            7
        }
    }

    fn key_at(input_start: usize) -> BatchKey {
        BatchKey {
            input_start,
            input_end: input_start + 1,
            batch_size: 100,
        }
    }

    fn job_at(input_start: usize) -> Job {
        Job {
            key: key_at(input_start),
            searched: input_start..input_start + 1,
            bytes: 10,
        }
    }

    /// A state with one slot running batch 1, batches 1-4 wanted and 2-4 queued.
    fn settling_state() -> State<Settled> {
        let mut state = State::new();
        state.running = vec![Some((key_at(1), 10))];
        state.wanted = (1..5).map(key_at).collect();
        state.queue = (2..5).map(job_at).collect();
        state
    }

    // Review R2: the slot decides "cancelled" from its flag after the search, not from what the
    // search returned: a search that stopped part-way and returned a value that looks complete
    // is never taken, and its key, still wanted, is searched again.
    #[test]
    fn a_flagged_search_is_never_taken_whatever_it_returned() {
        let mut state = settling_state();
        let unused = settle(&mut state, 0, job_at(1), Outcome::Value(Ok(1)), 7, true);
        assert!(matches!(unused, Some(Outcome::Value(Ok(1)))));
        assert!(state.done.is_empty());
        assert_eq!((state.cancelled, state.discarded), (1, 1));
        let queued: Vec<usize> = state.queue.iter().map(|job| job.key.input_start).collect();
        assert_eq!(queued, [1, 2, 3, 4]);
        assert_eq!(state.running, [None]);
        // Not flagged: the same value is kept for the loop.
        let mut state = settling_state();
        assert!(settle(&mut state, 0, job_at(1), Outcome::Value(Ok(1)), 7, false).is_none());
        assert_eq!(state.done.get(&key_at(1)).map(|(_, held)| *held), Some(7));
        assert_eq!(state.cancelled, 0);
        // A key no plan wants is dropped, not queued again.
        let mut state = settling_state();
        state.wanted.remove(&key_at(1));
        assert!(settle(&mut state, 0, job_at(1), Outcome::Value(Ok(1)), 7, true).is_some());
        assert_eq!(state.queue.len(), 3);
    }

    // Review R3: a slot's panic is counted, kept for the loop when its key is wanted (it then
    // ends the batches after it, as a failure does), dropped otherwise.
    #[test]
    fn slot_panics_are_counted_kept_or_not() {
        let mut state = settling_state();
        let panic = Outcome::Panic(Box::new("boom"), vec!["message\n".to_string()]);
        assert!(settle(&mut state, 0, job_at(1), panic, 0, false).is_none());
        assert_eq!(state.panicked, 1);
        assert_eq!(state.first_panic.as_deref(), Some("message\n"));
        assert!(state.done.contains_key(&key_at(1)));
        assert!(
            state.queue.is_empty(),
            "nothing after a failed batch starts"
        );
        let mut state = settling_state();
        state.wanted.clear();
        let panic = Outcome::Panic(Box::new("boom"), Vec::new());
        assert!(settle(&mut state, 0, job_at(1), panic, 0, false).is_some());
        assert_eq!((state.panicked, state.discarded), (1, 1));
        // A failed (not cancelled) result ends the queue after it too.
        let mut state = settling_state();
        assert!(settle(&mut state, 0, job_at(1), Outcome::Value(Err("e")), 7, false).is_none());
        assert!(state.queue.is_empty());
    }

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    mod executor {
        use super::*;
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::time::{Duration, Instant};

        /// A fake batch result: its key and the thread count of its pool; `Err("cancelled")`
        /// for a search that saw its cancellation.
        type Fake = Result<(BatchKey, usize), String>;

        impl SlotResult for Fake {
            fn failed(&self) -> bool {
                matches!(self, Err(error) if error != "cancelled")
            }

            fn cancelled(&self) -> bool {
                matches!(self, Err(error) if error == "cancelled")
            }

            fn held_bytes(&self) -> u64 {
                1
            }
        }

        /// What the fake searches do besides returning their key.
        #[derive(Clone, Copy, Default)]
        struct Fakes {
            /// The batch of this input start returns an error (the loop stops there).
            fail: Option<usize>,
            /// The batch of this input start panics.
            panic: Option<usize>,
            /// The failing batch waits until a slot has started the panicking one.
            fail_waits_for_panic: bool,
            /// A batch of another size than `BATCH_SIZE` (a corrupted prediction) waits for its
            /// cancellation (at most 10 s), or returns at once.
            stale: Stale,
            /// The loop waits this long after each batch it takes, so the slots run ahead.
            loop_delay: Option<Duration>,
            /// The lookup width the estimates assume (8 if 0).
            lut: usize,
        }

        /// The thread count a search reports when it saw its cancellation and returned anyway
        /// (`Stale::PartialOnCancel`): such a result must never be taken.
        const PARTIAL: usize = 9_999;

        #[derive(Clone, Copy, Default, PartialEq, Eq)]
        enum Stale {
            #[default]
            Normal,
            WaitsForCancel,
            /// Waits for its cancellation, then returns a value as if it were complete.
            PartialOnCancel,
            Quick,
        }

        #[derive(Debug)]
        struct Run {
            summary: Summary,
            results: Vec<Fake>,
            /// The keys the searches started, in order.
            started: Vec<BatchKey>,
            /// Searches running after `finish` (0: every slot was joined).
            active_after_finish: usize,
            line: String,
        }

        /// Runs the loop over `lengths` (fixed `batch_size`, or the mixer with a gie of one per
        /// ten letters) with `threads` threads, as `run_in_pool` does.
        fn run(
            lengths: &[usize],
            batch_size: Option<u32>,
            threads: usize,
            switches: Switches,
            fakes: Fakes,
        ) -> Run {
            let query_input = input_of(lengths, &vec![false; lengths.len()]);
            let searched = searched_lengths(&query_input);
            let input = PlanInput {
                query_input: &query_input,
                lengths: &searched,
                reading_error: false,
                fixed_batch_size: batch_size.is_some(),
                split: SplitSizes::default_for(false),
            };
            let active = AtomicUsize::new(0);
            let started = Mutex::new(Vec::new());
            let panic_started = AtomicBool::new(false);
            let search = |job: &Job, pool: &SearchPool<'_>| -> Fake {
                started.lock().unwrap().push(job.key);
                active.fetch_add(pool.threads(), Ordering::SeqCst);
                let stale = batch_size.is_some_and(|size| job.key.batch_size != size);
                match fakes.stale {
                    Stale::WaitsForCancel | Stale::PartialOnCancel if stale => {
                        let since = Instant::now();
                        while !cancelled(pool) && since.elapsed() < Duration::from_secs(10) {
                            std::thread::yield_now();
                        }
                    }
                    Stale::Quick if stale => {}
                    _ => std::thread::sleep(Duration::from_millis(2)),
                }
                active.fetch_sub(pool.threads(), Ordering::SeqCst);
                if cancelled(pool) {
                    if fakes.stale == Stale::PartialOnCancel {
                        return Ok((job.key, PARTIAL));
                    }
                    return Err("cancelled".to_string());
                }
                if fakes.panic == Some(job.key.input_start) {
                    panic_started.store(true, Ordering::SeqCst);
                    panic!("fake batch panic at {}", job.key.input_start);
                }
                if fakes.fail == Some(job.key.input_start) && fakes.fail_waits_for_panic {
                    while !panic_started.load(Ordering::SeqCst) {
                        std::thread::yield_now();
                    }
                }
                if pool.threads() > 1 {
                    assert_eq!(rayon::current_num_threads(), pool.threads());
                }
                if fakes.fail == Some(job.key.input_start) {
                    return Err(format!("error at {}", job.key.input_start));
                }
                Ok((job.key, pool.threads()))
            };
            let shared: Shared<Fake> = Shared::new(switches);
            let (results, active_after_finish) = std::thread::scope(|scope| {
                let mut driver = Driver::new(
                    &shared,
                    &search,
                    scope,
                    input,
                    threads,
                    1024,
                    Box::new(move |_, _| if fakes.lut == 0 { 8 } else { fakes.lut }),
                    None,
                );
                let mut mixer =
                    BatchSizeMixer::new(SplitSizes::default_for(false).mixer_max_batch_size());
                let mut size = batch_size.unwrap_or_else(|| mixer.batch_size(None) as u32);
                let mut start = 0;
                let mut results = Vec::new();
                while start < lengths.len() {
                    let (end, _) = next_ranged_batch(lengths, &query_input.skipped, start, size);
                    let key = BatchKey {
                        input_start: start,
                        input_end: end,
                        batch_size: size,
                    };
                    let result = match driver.take(key) {
                        Some(result) => result,
                        None => {
                            assert_eq!(rayon::current_thread_index(), None);
                            let job = Job {
                                key,
                                searched: start..end,
                                bytes: 0,
                            };
                            search(&job, &SearchPool::x_serial())
                        }
                    };
                    let stop = result.is_err();
                    results.push(result);
                    if stop {
                        break;
                    }
                    if let Some(delay) = fakes.loop_delay {
                        std::thread::sleep(delay);
                    }
                    let letters: u64 = lengths[start..end].iter().map(|&l| l as u64).sum();
                    let gie = letters / 10;
                    if batch_size.is_none() {
                        size = mixer.batch_size(Some(gie as u32 as i32)) as u32;
                    }
                    driver.commit(key, gie, batch_size.is_none().then_some(&mixer), size);
                    start = end;
                }
                driver.finish();
                (results, active.load(Ordering::SeqCst))
            });
            let started = started.into_inner().unwrap();
            Run {
                summary: shared.summary(),
                results,
                started,
                active_after_finish,
                line: shared.stats_line(),
            }
        }

        fn on() -> Switches {
            Switches {
                on: true,
                ..Switches::default()
            }
        }

        fn keys(results: &[Fake]) -> Vec<Result<BatchKey, String>> {
            results
                .iter()
                .map(|result| result.as_ref().map(|(key, _)| *key).map_err(Clone::clone))
                .collect()
        }

        fn off_run(lengths: &[usize], batch_size: Option<u32>) -> Run {
            run(
                lengths,
                batch_size,
                4,
                Switches::default(),
                Fakes::default(),
            )
        }

        #[test]
        fn slots_give_the_loops_batches_in_order() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(300));
            let run4 = run(&lengths, Some(300), 4, on(), Fakes::default());
            let stats = run4.summary.stats;
            assert_eq!(keys(&run4.results), keys(&off.results));
            assert_eq!(stats.batches as usize, run4.results.len());
            assert_eq!(stats.precomputed + stats.loop_searched, stats.batches);
            assert!(stats.precomputed >= stats.batches - 1, "{}", run4.line);
            assert_eq!((stats.k, stats.m, stats.idle), (4, 1, 0), "{}", run4.line);
            assert!(run4.summary.peak_busy <= 4);
            assert_eq!(run4.active_after_finish, 0);
            assert_eq!(run4.summary.panicked, 0);
            assert_eq!(off.summary.stats.precomputed, 0);
            assert_eq!(off.summary.stats.off, Some(Off::NotSet));
        }

        #[test]
        fn sub_pools_give_the_same_batches() {
            let lengths = vec![100; 30];
            let switches = Switches {
                slots: Some(2),
                ..on()
            };
            let run2 = run(&lengths, Some(200), 4, switches, Fakes::default());
            let stats = run2.summary.stats;
            assert_eq!((stats.k, stats.m), (2, 2), "{}", run2.line);
            assert!(stats.precomputed > 0, "{}", run2.line);
            assert!(run2
                .results
                .iter()
                .skip(1)
                .any(|result| matches!(result, Ok((_, 2)))));
        }

        #[test]
        fn the_mixer_batches_are_taken_from_the_slots() {
            let lengths: Vec<usize> = (0..400).map(|index| 3_000 + (index * 17) % 1_000).collect();
            let off = off_run(&lengths, None);
            let run4 = run(&lengths, None, 4, on(), Fakes::default());
            assert_eq!(keys(&run4.results), keys(&off.results));
            assert_eq!(run4.summary.stats.mispredicted, 0, "{}", run4.line);
            assert!(run4.summary.stats.precomputed > 0, "{}", run4.line);
        }

        // Corrupted predictions: the slots' batches under other keys are discarded, the loop's
        // batches do not change, and two misses within eight batches (or the waste bound) stop
        // the speculation.
        #[test]
        fn mispredicted_batches_are_discarded() {
            let lengths = vec![100; 60];
            let off = off_run(&lengths, Some(200));
            for every in [1, 2, 5] {
                let switches = Switches {
                    slots: Some(2),
                    miss_every: Some(every),
                    ..on()
                };
                let missed = run(&lengths, Some(200), 4, switches, Fakes::default());
                assert_eq!(keys(&missed.results), keys(&off.results), "every {every}");
                assert!(
                    missed.summary.stats.mispredicted > 0,
                    "every {every}: {}",
                    missed.line
                );
                // Every prediction corrupted: the second batch's and third's next keys are
                // both wrong, before the waste bound can be reached. With fewer corrupted, the
                // discards may reach the waste bound first.
                let stops: &[Off] = if every == 1 {
                    &[Off::Guard]
                } else {
                    &[Off::Guard, Off::Waste]
                };
                assert!(
                    missed
                        .summary
                        .stats
                        .off
                        .is_some_and(|off| stops.contains(&off)),
                    "every {every}: {}",
                    missed.line
                );
                assert_eq!(missed.summary.panicked, 0);
            }
        }

        // A running batch that a new plan drops is cancelled, and its result is never taken.
        #[test]
        fn a_running_batch_that_the_plan_drops_is_cancelled() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(2),
                miss_every: Some(1),
                ..on()
            };
            let fakes = Fakes {
                stale: Stale::WaitsForCancel,
                ..Fakes::default()
            };
            let cancelled_run = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&cancelled_run.results), keys(&off.results));
            assert!(
                cancelled_run.summary.cancelled > 0,
                "{}",
                cancelled_run.line
            );
            assert_eq!(cancelled_run.active_after_finish, 0);
        }

        // Discarded batches beyond the waste bound stop the speculation, also when no next key
        // is mispredicted (four slots plan eight batches with seven predictions; every seventh
        // corrupted makes only the last planned batch wrong, and the loop waits after each batch
        // so that the slots search it before the next plan drops it).
        #[test]
        fn waste_beyond_the_bound_stops_the_speculation() {
            let lengths = vec![100; 180];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(4),
                miss_every: Some(7),
                ..on()
            };
            let fakes = Fakes {
                stale: Stale::Quick,
                loop_delay: Some(Duration::from_millis(5)),
                ..Fakes::default()
            };
            let wasted = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&wasted.results), keys(&off.results));
            assert_eq!(wasted.summary.stats.mispredicted, 0, "{}", wasted.line);
            assert_eq!(
                wasted.summary.stats.off,
                Some(Off::Waste),
                "{}",
                wasted.line
            );
        }

        // The memory budget: with none, a slot starts only the loop's next batch, so one batch
        // runs at a time.
        #[test]
        fn the_budget_limits_the_batches_in_flight() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(4),
                mem_mib: Some(0),
                ..on()
            };
            let tight = run(&lengths, Some(200), 4, switches, Fakes::default());
            assert_eq!(keys(&tight.results), keys(&off.results));
            assert_eq!(tight.summary.peak_busy, 1, "{}", tight.line);
            assert!(tight.summary.stats.precomputed > 0, "{}", tight.line);
        }

        // An error at batch j ends the loop there; batches the slots searched past j are
        // dropped, and so is a panic of one of them, which is never resumed (batch j's search
        // waits until a slot has started the batch that panics).
        #[test]
        fn an_error_ends_the_loop_and_drops_batches_past_it() {
            let lengths = vec![100; 40];
            let fakes = Fakes {
                fail: Some(10),
                panic: Some(12),
                fail_waits_for_panic: true,
                ..Fakes::default()
            };
            let failed = run(&lengths, Some(100), 4, on(), fakes);
            assert_eq!(failed.results.len(), 11);
            assert_eq!(failed.results[10], Err("error at 10".to_string()));
            assert!(failed.summary.discarded >= 1, "{}", failed.line);
            assert_eq!(failed.active_after_finish, 0);
            // The panic of batch 12, past the error, is counted though never resumed.
            assert_eq!(failed.summary.panicked, 1, "{}", failed.line);
        }

        // Once a slot's batch returns an error, the batches planned after it are not started:
        // one slot searches the batches in order and stops at the failing one.
        #[test]
        fn no_batch_after_a_failed_one_starts() {
            let lengths = vec![100; 40];
            let switches = Switches {
                slots: Some(1),
                ..on()
            };
            let fakes = Fakes {
                fail: Some(10),
                ..Fakes::default()
            };
            let failed = run(&lengths, Some(100), 4, switches, fakes);
            assert_eq!(failed.results.len(), 11);
            let last = failed
                .started
                .iter()
                .map(|key| key.input_start)
                .max()
                .unwrap();
            assert_eq!(last, 10, "{:?}", failed.started);
        }

        // A panic of a batch the loop takes is resumed on the loop's thread.
        #[test]
        fn a_panic_of_a_taken_batch_is_resumed() {
            let lengths = vec![100; 20];
            let caught = std::panic::catch_unwind(|| {
                let fakes = Fakes {
                    panic: Some(5),
                    ..Fakes::default()
                };
                run(&lengths, Some(100), 4, on(), fakes)
            });
            let payload = caught.expect_err("the panic of batch 5 is resumed");
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .unwrap_or_default();
            assert_eq!(message, "fake batch panic at 5");
        }

        // Review R2 on whole runs: stale searches that see their cancellation and return a value
        // as if complete are counted as cancelled and never taken.
        #[test]
        fn a_value_returned_after_a_cancel_is_never_taken() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(2),
                miss_every: Some(1),
                ..on()
            };
            let fakes = Fakes {
                stale: Stale::PartialOnCancel,
                ..Fakes::default()
            };
            let partial = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&partial.results), keys(&off.results));
            assert!(partial
                .results
                .iter()
                .all(|result| !matches!(result, Ok((_, PARTIAL)))));
            assert!(partial.summary.cancelled > 0, "{}", partial.line);
            assert_eq!(partial.summary.panicked, 0);
        }

        // A budget between none and ample: the estimates (lookup width 10, 4 MiB each) and a
        // budget of 5 MiB admit the loop's next batch and one more, so two of the four slots
        // search at a time while the queued batches wait.
        #[test]
        fn a_budget_admits_some_batches_and_holds_the_rest() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(4),
                mem_mib: Some(5),
                ..on()
            };
            let fakes = Fakes {
                lut: 10,
                ..Fakes::default()
            };
            let bounded = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&bounded.results), keys(&off.results));
            assert_eq!(bounded.summary.peak_busy, 2, "{}", bounded.line);
            let one = batch_estimate_bytes(10, 200, 20);
            assert!(
                bounded.summary.peak_bytes <= 5 * MIB + 2 * one,
                "{}",
                bounded.line
            );
            assert!(bounded.summary.stats.precomputed > 0, "{}", bounded.line);
        }

        // Review R3: in the shadow mode a slot's panic is a difference even when its batch is
        // dropped (here past an error at batch 10): the loop stops with the shadow's panic.
        #[test]
        fn the_shadow_fails_on_a_dropped_slot_panic() {
            let lengths = vec![100; 40];
            let switches = Switches {
                shadow: true,
                ..on()
            };
            let fakes = Fakes {
                fail: Some(10),
                panic: Some(12),
                fail_waits_for_panic: true,
                ..Fakes::default()
            };
            let caught = std::panic::catch_unwind(|| run(&lengths, Some(100), 4, switches, fakes));
            let payload = caught.expect_err("the shadow fails on the slot's panic");
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .unwrap_or_default();
            assert!(
                message.starts_with("[X_SHADOW] QBATCHPAR a batch searched by a slot panicked"),
                "{message}"
            );
            assert!(message.contains("fake batch panic at 12"), "{message}");
        }

        #[test]
        fn one_batch_after_the_first_does_not_speculate() {
            let lengths = vec![100; 2];
            let single = run(&lengths, Some(100), 4, on(), Fakes::default());
            assert_eq!(single.summary.stats.precomputed, 0);
            assert_eq!(
                single.summary.stats.off,
                Some(Off::Single),
                "{}",
                single.line
            );
        }
    }
}
