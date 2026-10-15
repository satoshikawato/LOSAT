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
//!   new plan no longer holds is dropped (queued) or discarded when it finishes (running). After
//!   two mispredicted next keys within eight batches the search stops speculating and the loop
//!   searches the rest itself, on the whole pool.
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
//! so there a panic in any batch a slot searches ends the process with its message, also in a
//! batch the loop would not take (a batch past an error, or of a mispredicted key); panics are
//! bugs, and the loop's own batches would panic the same. `LOSAT_TIMING` blocks are text in the
//! batch's result, printed when the loop takes it.
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
//! Diagnostics: with the switch set one line on standard error at the end of the search,
//! `[X_STATS] QBATCHPAR batches= precomputed= loop= mispredicted= discarded= k= m= budget_mib=
//! estimate_mib= off=` (precomputed: batches taken from a slot; loop: batches the loop searched
//! itself; off: why the search did not, or stopped to, speculate, `none` if it never stopped). With
//! the shadow, `[X_SHADOW] QBATCHPAR batches= compared= diffs= mispredicted= wasted=`; the first
//! difference prints it and panics. `LOSAT_X_QBATCHPARMISS=N` (test only) corrupts every N-th
//! prediction so that the discard, replan and guard paths run.
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

/// No NCBI counterpart: the switches of this path, read once per search; they choose a LOSAT
/// code path and change no value NCBI computes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Switches {
    /// `LOSAT_X_QBATCHPAR` (any value), or the shadow.
    pub on: bool,
    /// `LOSAT_X_QBATCHPARSHADOW` (any value; implies the switch).
    pub shadow: bool,
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
    /// The loop's next key was not planned (not expected to happen).
    Unplanned,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
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

/// A batch for a slot: its key and its searched queries.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Job {
    pub key: BatchKey,
    pub searched: Range<usize>,
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

/// What a slot made of a batch: its value, or its panic with the messages the panic hook
/// kept instead of printing them.
pub(super) enum Outcome<T> {
    Value(T),
    Panic(Box<dyn std::any::Any + Send>, Vec<String>),
}

impl<T> Outcome<T> {
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
}

/// The batches of the slots.
struct State<T> {
    /// Planned batches that no slot has started, in batch order.
    queue: VecDeque<Job>,
    /// The batch each slot is searching.
    running: Vec<Option<BatchKey>>,
    /// Finished batches of the current plan.
    done: HashMap<BatchKey, Outcome<T>>,
    /// The keys of the current plan.
    wanted: HashSet<BatchKey>,
    shutdown: bool,
    /// Finished batches dropped unused.
    discarded: u64,
    /// The most slots searching at once.
    peak_busy: usize,
    /// Slot threads that have not ended.
    alive: usize,
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
    pub budget_mib: u64,
    pub estimate_mib: u64,
    pub off: Option<Off>,
    pub compared: u64,
    pub diffs: u64,
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
    stats: Mutex<Stats>,
}

impl<T> Shared<T> {
    pub(super) fn new(switches: Switches) -> Self {
        Shared {
            switches,
            state: Mutex::new(State {
                queue: VecDeque::new(),
                running: Vec::new(),
                done: HashMap::new(),
                wanted: HashSet::new(),
                shutdown: false,
                discarded: 0,
                peak_busy: 0,
                alive: 0,
            }),
            wake: Condvar::new(),
            sinks: OnceLock::new(),
            stats: Mutex::new(Stats {
                off: Some(Off::Unstarted),
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

    /// The counters so far, with the discarded batches and the peak of busy slots.
    pub(super) fn stats(&self) -> (Stats, u64, usize) {
        let stats = *self.stats.lock().unwrap_or_else(PoisonError::into_inner);
        let state = self.lock();
        (stats, state.discarded, state.peak_busy)
    }

    /// `[X_STATS] QBATCHPAR ...`
    pub(super) fn stats_line(&self) -> String {
        let (stats, discarded, _) = self.stats();
        format!(
            "[X_STATS] QBATCHPAR batches={} precomputed={} loop={} mispredicted={} discarded={} k={} m={} budget_mib={} estimate_mib={} off={}",
            stats.batches,
            stats.precomputed,
            stats.loop_searched,
            stats.mispredicted,
            discarded,
            stats.k,
            stats.m,
            stats.budget_mib,
            stats.estimate_mib,
            stats.off.map_or("none", Off::name)
        )
    }

    /// `[X_SHADOW] QBATCHPAR ...`
    pub(super) fn shadow_line(&self) -> String {
        let (stats, discarded, _) = self.stats();
        format!(
            "[X_SHADOW] QBATCHPAR batches={} compared={} diffs={} mispredicted={} wasted={}",
            stats.batches, stats.compared, stats.diffs, stats.mispredicted, discarded
        )
    }
}

impl<T> Drop for Shared<T> {
    // No NCBI counterpart: the summary lines at the end of a search with the switch set.
    fn drop(&mut self) {
        if !self.switches.on || std::thread::panicking() {
            return;
        }
        eprintln!("{}", self.stats_line());
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
}

impl<'scope, 'env, T, F> Driver<'scope, 'env, T, F>
where
    T: Send + 'env,
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
        let off = if switches.on { off } else { Some(Off::Build) };
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
                if let Some(outcome) = state.done.remove(&key) {
                    drop(state);
                    self.stats.precomputed += 1;
                    self.flush();
                    return Some(outcome.resume());
                }
                let planned = state.running.contains(&Some(key))
                    || state.queue.iter().any(|job| job.key == key);
                if !planned {
                    break Off::Unplanned;
                }
                // Every slot's thread ended outside a batch (not expected): nothing will search
                // the queued key.
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
            if let (Some(predicted), Some(actual)) = (predicted, next.first()) {
                if predicted != actual.key {
                    self.stats.mispredicted += 1;
                    let commit = self.commits;
                    if self
                        .last_miss
                        .is_some_and(|last| commit - last < GUARD_COMMITS)
                    {
                        self.stop(Off::Guard);
                        self.flush();
                        return;
                    }
                    self.last_miss = Some(commit);
                }
            }
            self.replan(next);
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
            .map(|planned| {
                batch_estimate_bytes(
                    (self.lut_width)(planned.letters, planned.max_len),
                    planned.letters,
                    planned.gie,
                )
            })
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
        match self.start_slots(k, m) {
            Ok(()) => {
                self.phase = Phase::Active;
                self.stats.off = None;
                self.window = 2 * k;
            }
            Err(off) => self.stop(off),
        }
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
        if self.shared.sinks.set(sinks).is_err() {
            return Err(Off::Spawn);
        }
        self.shared.lock().running = vec![None; k];
        for slot in 0..k {
            let (shared, search) = (self.shared, self.search);
            self.shared.lock().alive += 1;
            let spawned = std::thread::Builder::new()
                .name(format!("losat-qbp-{slot}"))
                .spawn_scoped(self.scope, move || slot_main(shared, search, slot, m));
            if spawned.is_err() {
                self.shared.lock().alive -= 1;
                return Err(Off::Spawn);
            }
        }
        Ok(())
    }

    /// Stop speculating: nothing more is started; finished results are dropped, running ones
    /// when they finish.
    fn stop(&mut self, off: Off) {
        self.set_off(off);
        self.plan.clear();
        let unused: Vec<Outcome<T>> = {
            let mut state = self.shared.lock();
            state.queue.clear();
            state.wanted.clear();
            let unused: Vec<Outcome<T>> = state.done.drain().map(|(_, outcome)| outcome).collect();
            state.discarded += unused.len() as u64;
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
    }

    /// The plan `next` replaces the current one: its keys that no slot has started or finished
    /// are queued in order; finished batches outside it are dropped.
    fn replan(&mut self, next: Vec<Planned>) {
        self.plan = next.iter().map(|planned| planned.key).collect();
        let wanted: HashSet<BatchKey> = self.plan.iter().copied().collect();
        let unused: Vec<Outcome<T>> = {
            let mut state = self.shared.lock();
            let stale: Vec<BatchKey> = state
                .done
                .keys()
                .filter(|key| !wanted.contains(key))
                .copied()
                .collect();
            let unused: Vec<Outcome<T>> = stale
                .iter()
                .filter_map(|key| state.done.remove(key))
                .collect();
            state.discarded += unused.len() as u64;
            let queue: VecDeque<Job> = next
                .into_iter()
                .filter(|planned| {
                    !state.running.contains(&Some(planned.key))
                        && !state.done.contains_key(&planned.key)
                })
                .map(|planned| Job {
                    key: planned.key,
                    searched: planned.searched,
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
    // The end of the loop (also by an error or a panic): the slots stop after their batch, and
    // the thread scope joins them.
    fn drop(&mut self) {
        let unused: Vec<Outcome<T>> = {
            let mut state = self.shared.lock();
            state.shutdown = true;
            state.queue.clear();
            state.wanted.clear();
            let unused: Vec<Outcome<T>> = state.done.drain().map(|(_, outcome)| outcome).collect();
            state.discarded += unused.len() as u64;
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
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
/// cannot be built leaves the slot searching on its own thread.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn slot_main<T, F>(shared: &Shared<T>, search: &F, slot: usize, threads: usize)
where
    T: Send,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync,
{
    let sink = Arc::clone(
        &shared
            .sinks
            .get()
            .expect("slot sinks are set before the slots")[slot],
    );
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
            slot_loop(shared, search, slot, pool, &sink);
            Ok(())
        });
        if served.is_ok() {
            return;
        }
    }
    slot_loop(shared, search, slot, &SearchPool::x_serial(), &sink);
}

/// The batches of one slot, each searched with `pool`, until the loop ends.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn slot_loop<T, F>(
    shared: &Shared<T>,
    search: &F,
    slot: usize,
    pool: &SearchPool<'_>,
    sink: &PanicSink,
) where
    T: Send,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync,
{
    loop {
        let job = {
            let mut state = shared.lock();
            loop {
                if state.shutdown {
                    return;
                }
                if let Some(job) = state.queue.pop_front() {
                    state.running[slot] = Some(job.key);
                    let busy = state.running.iter().flatten().count();
                    state.peak_busy = state.peak_busy.max(busy);
                    break job;
                }
                state = shared.wait(state);
            }
        };
        sink.lock().unwrap_or_else(PoisonError::into_inner).clear();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| search(&job, pool)));
        let outcome = match result {
            Ok(value) => Outcome::Value(value),
            Err(payload) => Outcome::Panic(
                payload,
                std::mem::take(&mut *sink.lock().unwrap_or_else(PoisonError::into_inner)),
            ),
        };
        let unused = {
            let mut state = shared.lock();
            state.running[slot] = None;
            if !state.shutdown && state.wanted.contains(&job.key) {
                state.done.insert(job.key, outcome);
                None
            } else {
                state.discarded += 1;
                Some(outcome)
            }
        };
        shared.wake.notify_all();
        drop(unused);
    }
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

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    mod executor {
        use super::*;
        use std::sync::atomic::{AtomicUsize, Ordering};

        /// A fake batch result: its key and the thread count of its pool.
        type Fake = Result<(BatchKey, usize), String>;

        #[derive(Debug)]
        struct Run {
            stats: Stats,
            discarded: u64,
            peak_busy: usize,
            results: Vec<Fake>,
            line: String,
        }

        /// Runs the loop over `lengths` (fixed `batch_size`, or the mixer with `gie` per
        /// letter) with `threads` threads; `fail` makes the batch of that input start an error
        /// (the loop stops there), `panic` makes it panic.
        #[allow(clippy::too_many_arguments)]
        fn run(
            lengths: &[usize],
            batch_size: Option<u32>,
            threads: usize,
            switches: Switches,
            fail: Option<usize>,
            panic: Option<usize>,
            active: &AtomicUsize,
            peak: &AtomicUsize,
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
            let panic_started = std::sync::atomic::AtomicBool::new(false);
            let search = |job: &Job, pool: &SearchPool<'_>| -> Fake {
                let now = active.fetch_add(pool.threads(), Ordering::SeqCst) + pool.threads();
                peak.fetch_max(now, Ordering::SeqCst);
                std::thread::sleep(std::time::Duration::from_millis(2));
                active.fetch_sub(pool.threads(), Ordering::SeqCst);
                if panic == Some(job.key.input_start) {
                    panic_started.store(true, Ordering::SeqCst);
                    panic!("fake batch panic at {}", job.key.input_start);
                }
                // The failing batch past which a batch panics waits until a slot has started
                // that batch (with a panic before the failure, it is never reached).
                if fail == Some(job.key.input_start)
                    && panic.is_some_and(|at| at > job.key.input_start)
                {
                    while !panic_started.load(Ordering::SeqCst) {
                        std::thread::yield_now();
                    }
                }
                if pool.threads() > 1 {
                    assert_eq!(rayon::current_num_threads(), pool.threads());
                }
                if fail == Some(job.key.input_start) {
                    return Err(format!("error at {}", job.key.input_start));
                }
                Ok((job.key, pool.threads()))
            };
            let shared: Shared<Fake> = Shared::new(switches);
            let results = std::thread::scope(|scope| {
                let mut driver = Driver::new(
                    &shared,
                    &search,
                    scope,
                    input,
                    threads,
                    1024,
                    Box::new(|_, _| 8),
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
                            let pool = SearchPool::x_serial();
                            search(
                                &Job {
                                    key,
                                    searched: start..end,
                                },
                                &pool,
                            )
                        }
                    };
                    let stop = result.is_err();
                    results.push(result);
                    if stop {
                        break;
                    }
                    let letters: u64 = lengths[start..end].iter().map(|&l| l as u64).sum();
                    let gie = letters / 10;
                    if batch_size.is_none() {
                        size = mixer.batch_size(Some(gie as u32 as i32)) as u32;
                    }
                    driver.commit(key, gie, batch_size.is_none().then_some(&mixer), size);
                    start = end;
                }
                results
            });
            let (stats, discarded, peak_busy) = shared.stats();
            let line = shared.stats_line();
            Run {
                stats,
                discarded,
                peak_busy,
                results,
                line,
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

        #[test]
        fn slots_give_the_loops_batches_in_order() {
            let lengths = vec![100; 40];
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let off = run(
                &lengths,
                Some(300),
                4,
                Switches::default(),
                None,
                None,
                &active,
                &peak,
            );
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let run4 = run(&lengths, Some(300), 4, on(), None, None, &active, &peak);
            assert_eq!(keys(&run4.results), keys(&off.results));
            assert_eq!(run4.stats.batches as usize, run4.results.len());
            assert_eq!(
                run4.stats.precomputed + run4.stats.loop_searched,
                run4.stats.batches
            );
            assert!(
                run4.stats.precomputed >= run4.stats.batches - 1,
                "{}",
                run4.line
            );
            assert_eq!((run4.stats.k, run4.stats.m), (4, 1), "{}", run4.line);
            assert!(run4.peak_busy <= 4 && peak.load(Ordering::SeqCst) <= 4);
            assert_eq!(off.stats.precomputed, 0);
        }

        #[test]
        fn sub_pools_give_the_same_batches() {
            let lengths = vec![100; 30];
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let switches = Switches {
                slots: Some(2),
                ..on()
            };
            let run2 = run(&lengths, Some(200), 4, switches, None, None, &active, &peak);
            assert_eq!((run2.stats.k, run2.stats.m), (2, 2), "{}", run2.line);
            assert!(run2.stats.precomputed > 0, "{}", run2.line);
            assert!(run2
                .results
                .iter()
                .skip(1)
                .any(|result| matches!(result, Ok((_, 2)))));
            assert!(peak.load(Ordering::SeqCst) <= 4);
        }

        #[test]
        fn the_mixer_batches_are_taken_from_the_slots() {
            let lengths: Vec<usize> = (0..400).map(|index| 3_000 + (index * 17) % 1_000).collect();
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let off = run(
                &lengths,
                None,
                4,
                Switches::default(),
                None,
                None,
                &active,
                &peak,
            );
            let run4 = run(&lengths, None, 4, on(), None, None, &active, &peak);
            assert_eq!(keys(&run4.results), keys(&off.results));
            assert_eq!(run4.stats.mispredicted, 0, "{}", run4.line);
            assert!(run4.stats.precomputed > 0, "{}", run4.line);
        }

        // Corrupted predictions: the slots' batches under other keys are discarded, the loop's
        // batches do not change, and the guard stops the speculation.
        #[test]
        fn mispredicted_batches_are_discarded() {
            let lengths = vec![100; 60];
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let off = run(
                &lengths,
                Some(200),
                4,
                Switches::default(),
                None,
                None,
                &active,
                &peak,
            );
            for every in [1, 2, 5] {
                let switches = Switches {
                    slots: Some(2),
                    miss_every: Some(every),
                    ..on()
                };
                let missed = run(&lengths, Some(200), 4, switches, None, None, &active, &peak);
                assert_eq!(keys(&missed.results), keys(&off.results), "every {every}");
                assert!(
                    missed.stats.mispredicted > 0,
                    "every {every}: {}",
                    missed.line
                );
                assert_eq!(
                    missed.stats.off,
                    Some(Off::Guard),
                    "every {every}: {}",
                    missed.line
                );
            }
        }

        // An error at batch j ends the loop there; batches the slots searched past j are
        // dropped, and so is a panic of one of them, which is never resumed (batch j's search
        // waits until a slot has started the batch that panics).
        #[test]
        fn an_error_ends_the_loop_and_drops_batches_past_it() {
            let lengths = vec![100; 40];
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let failed = run(
                &lengths,
                Some(100),
                4,
                on(),
                Some(10),
                Some(12),
                &active,
                &peak,
            );
            assert_eq!(failed.results.len(), 11);
            assert_eq!(failed.results[10], Err("error at 10".to_string()));
            assert!(failed.discarded >= 1, "{}", failed.line);
        }

        // A panic of a batch the loop takes is resumed on the loop's thread.
        #[test]
        fn a_panic_of_a_taken_batch_is_resumed() {
            let lengths = vec![100; 20];
            let caught = std::panic::catch_unwind(|| {
                let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
                run(&lengths, Some(100), 4, on(), None, Some(5), &active, &peak)
            });
            let payload = caught.expect_err("the panic of batch 5 is resumed");
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .unwrap_or_default();
            assert_eq!(message, "fake batch panic at 5");
        }

        #[test]
        fn one_batch_after_the_first_does_not_speculate() {
            let lengths = vec![100; 2];
            let (active, peak) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let single = run(&lengths, Some(100), 4, on(), None, None, &active, &peak);
            assert_eq!(single.stats.precomputed, 0);
            assert_eq!(single.stats.off, Some(Off::Single), "{}", single.line);
        }
    }
}
