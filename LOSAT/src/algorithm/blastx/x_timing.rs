//! Diagnostic (LOSAT_TIMING): where the wall time of a BLASTX run goes, so that the serial parts
//! of a run with many threads have numbers.
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The module reads a monotonic clock and
//! adds differences to counters; it computes no value NCBI computes and no value the search
//! reads. With `LOSAT_TIMING` not set it never reads the clock.
//!
//! # What is measured
//!
//! `native::run` prints these lines after the epilog (`write`):
//!
//! - `[TIMING] blastx_run: startup .. prolog .. read_stores .. loop .. epilog .. total ..
//!   (search .. render .. titles .. input ..)`: wall times on the thread that runs `run`.
//!   `startup` is everything before the prolog (options, the `-subject` file, the database's
//!   index); `read_stores` reads the database's sequences; `loop` is the query loop (the pool's
//!   start and end included). Inside `loop`: `search` is the time the loop waits for searches,
//!   `render` the time it writes reports (of which `titles` reads the reported sequences' titles
//!   with `-db`), and `input = loop - search - render` the reader's own time (and the queueing of
//!   batches with `LOSAT_X_BXBATCH`).
//! - `[TIMING] blastx_waves: waves=N batches=N wall ..s busy ..s threads=N idle ..s tail=N
//!   tail_wall ..s`: with `LOSAT_X_BXBATCH`, the waves of batches searched side by side.
//!   `wall` is their wall time, `busy` the searches' own times added over the threads,
//!   `idle = wall x threads - busy` the thread time the waves left unused (a wave ends when its
//!   last batch does). `tail` counts the batches searched one after another, each with the whole
//!   pool, because fewer than the threads were queued (the end of the input, or a long query,
//!   which is searched at once); `tail_wall` is their wall time.
//! - `[TIMING] blastx_flow: batches=N alone=N main_jobs=N main_job_time ..s blocked ..s`: with
//!   `LOSAT_X_BXFLOW` (instead of the waves), the batches read, those searched alone (with the
//!   pool inside) at the end of the input, the tasks the reading thread ran itself while it
//!   waited (and their time), and the time it waited for a result with no task to run. Both
//!   waits count as `search` of `blastx_run`.
//! - `[TIMING] blastx_subject_waves: waves=N tasks=N wall ..s busy ..s longest ..s latest_start ..s
//!   last_end ..s threads=N idle ..s`: the parallel subject waves of `search::search_core` (the
//!   `scan_par` stage). A task is the subjects one worker searches in a wave (`LOSAT_X_BXCHUNK`
//!   of them, or one). `wall` adds the waves' wall times, `busy` the tasks' own times, `longest`
//!   each wave's longest task, `latest_start` each wave's latest task start and `last_end` each
//!   wave's latest task end, both after the wave began; `idle = wall x threads - busy`. A wave
//!   ends with its last task: `longest` near `wall` with `busy` well below `wall x threads` is
//!   uneven work, `latest_start` the time workers take to begin, `wall - last_end` the join.
//! - `[TIMING] blastx_subject_flow: loops=N tasks=N ahead=N taken=N wall ..s`: with
//!   `LOSAT_X_BXSUBJFLOW` (instead of the waves), the subject loops run as a flow, their tasks,
//!   the tasks searched ahead of the loop (by a helper, or by the loop's thread while it waited
//!   for a helper's task) and those of them the loop took (it searched the others itself at their
//!   turn), and the loops' wall time, the hand-over to the sink included (it is also the
//!   `scan_par` of the batch line, whose `replay` is then 0).
//! - `[TIMING] blastx_batch_<mode>: batches=N query_setup .. setup .. chunks_par .. lookup ..
//!   scan_par .. scan_inline .. replay .. chunk_end .. batch_end .. sum ..`: the stages of
//!   `search::search_core`, per query batch, in the order they run. Mode `side` is a batch
//!   searched on one thread beside other batches (`runtime::x_search_internal_serial`): its times
//!   are added over the threads. Mode `alone` is a batch searched by itself, with the pool if
//!   there is one: its times are the wall times of the thread that runs it. `query_setup`
//!   translates and masks the queries; `setup` is the batch's statistics and query chunks;
//!   `chunks_par` the parallel chunk path of a long query (`LOSAT_X_BXPAR`); `lookup` builds a
//!   chunk's lookup table and diagonal tables; `scan_par` is the wall time of the parallel
//!   subject waves (scan, ungapped and preliminary gapped extension); `replay` is the loop that
//!   hands the subjects to the batch's sink in subject order (linking, reaping, the HSP
//!   stream), and `scan_inline` the part of it that computes subjects inline (one thread);
//!   `batch_end` is the end of the batch (`runtime::Runtime::batch_end`); `sum` is the whole
//!   batch, so `sum` minus the stages is what lies between them.
//! - `[TIMING] blastx_batch_end_<mode>: init .. read .. kappa_spec .. kappa_serial .. loop ..
//!   finish .. (matches=N reused=N serial=N)`: inside `batch_end`. `init` is up to the
//!   composition state, `read` takes the matches from the stream, `loop` is the in-order loop
//!   over the matches, of which `kappa_spec` is the parallel redo made ahead (`LOSAT_X_BXPAR`)
//!   and `kappa_serial` the redo of the matches done in the loop itself; `finish` is from the
//!   heaps to the batch's results (hit lists, culling, sorting, the copy of the batch).
//!   `matches` counts the matches, `reused` those whose redo was made ahead and `serial` those
//!   redone in the loop.
//! - `[TIMING] blastx_diag: mode=.. subjects=N cells=N batched_hits=N candidates=N`: with
//!   `LOSAT_X_BXDIAG`, the fused seed loop's cell test and counts (`seed/x_seed_diag.rs`).

use std::io::{self, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;
use std::time::Instant;

// No NCBI counterpart: whether LOSAT_TIMING is set, read once.
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

// No NCBI counterpart: the clock, read only when LOSAT_TIMING is set.
#[inline]
pub(crate) fn lap() -> Option<Instant> {
    if on() {
        Some(Instant::now())
    } else {
        None
    }
}

// No NCBI counterpart: nanoseconds since `start` (0 without a start).
#[inline]
fn since(start: Option<Instant>) -> u64 {
    start.map_or(0, |t| {
        u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
    })
}

// No NCBI counterpart: names of the run's wall times and counts.
/// The intervals and counts of `blastx_run` and `blastx_waves`.
#[derive(Clone, Copy)]
pub(crate) enum Run {
    Startup,
    Prolog,
    ReadStores,
    Loop,
    Epilog,
    Total,
    Search,
    Render,
    Titles,
    Waves,
    WaveBatches,
    WaveWall,
    WaveBusy,
    TailBatches,
    TailWall,
    FlowBatches,
    FlowJobs,
    FlowJobTime,
    FlowBlocked,
    FlowAlone,
    SubjectWaves,
    SubjectTasks,
    SubjectWall,
    SubjectBusy,
    SubjectLongest,
    SubjectLatestStart,
    SubjectLastEnd,
    FlowLoops,
    FlowTasks,
    FlowAhead,
    FlowTaken,
    FlowWall,
}
const RUN_LEN: usize = Run::FlowWall as usize + 1;

// No NCBI counterpart: names of the stages of `search_core`.
/// The stages of one query batch in `search::search_core`, in the order they run.
#[derive(Clone, Copy)]
pub(crate) enum Batch {
    Batches,
    QuerySetup,
    Setup,
    ChunksPar,
    Lookup,
    ScanPar,
    ScanInline,
    Replay,
    ChunkEnd,
    BatchEnd,
    Sum,
}
const BATCH_LEN: usize = Batch::Sum as usize + 1;

// No NCBI counterpart: names of the intervals of `batch_end`.
/// The intervals and counts of `runtime::Runtime::batch_end`.
#[derive(Clone, Copy)]
pub(crate) enum End {
    Init,
    Read,
    KappaSpec,
    KappaSerial,
    Loop,
    Finish,
    Matches,
    Reused,
    Serial,
}
const END_LEN: usize = End::Serial as usize + 1;

// No NCBI counterpart: the counters (nanoseconds or counts); the batch counters per mode
// (0 = side by side, 1 = alone).
static RUN: [AtomicU64; RUN_LEN] = [const { AtomicU64::new(0) }; RUN_LEN];
static BATCH: [[AtomicU64; BATCH_LEN]; 2] = [const { [const { AtomicU64::new(0) }; BATCH_LEN] }; 2];
static END: [[AtomicU64; END_LEN]; 2] = [const { [const { AtomicU64::new(0) }; END_LEN] }; 2];

// No NCBI counterpart: the mode of the batch this thread runs.
#[inline]
fn mode() -> usize {
    usize::from(!super::runtime::x_inner_serial())
}

// No NCBI counterpart: adds the time since `start` to a run interval; returns it.
#[inline]
pub(crate) fn run(stage: Run, start: Option<Instant>) -> u64 {
    let ns = since(start);
    if start.is_some() {
        RUN[stage as usize].fetch_add(ns, Ordering::Relaxed);
    }
    ns
}

// No NCBI counterpart: adds `n` to a run count or interval.
#[inline]
pub(crate) fn run_add(stage: Run, n: u64) {
    if on() {
        RUN[stage as usize].fetch_add(n, Ordering::Relaxed);
    }
}

// No NCBI counterpart: the clocks of one wave of the parallel subject loop.
/// One wave of `search::search_core`'s parallel subject loop (`blastx_subject_waves`): when it
/// began and, over its tasks, the longest one, the latest start and end, and their sum.
pub(crate) struct SubjectWave {
    start: Option<Instant>,
    tasks: AtomicU64,
    busy: AtomicU64,
    longest: AtomicU64,
    latest_start: AtomicU64,
    last_end: AtomicU64,
}

impl SubjectWave {
    // No NCBI counterpart: a wave that begins now (no clock without LOSAT_TIMING).
    pub(crate) fn begin() -> Self {
        Self {
            start: lap(),
            tasks: AtomicU64::new(0),
            busy: AtomicU64::new(0),
            longest: AtomicU64::new(0),
            latest_start: AtomicU64::new(0),
            last_end: AtomicU64::new(0),
        }
    }

    // No NCBI counterpart: runs one task of the wave, timing it when the wave has a clock.
    #[inline]
    pub(crate) fn task<T>(&self, f: impl FnOnce() -> T) -> T {
        let Some(wave_start) = self.start else {
            return f();
        };
        let t0 = Instant::now();
        let result = f();
        let ns = |d: std::time::Duration| u64::try_from(d.as_nanos()).unwrap_or(u64::MAX);
        let began = ns(t0.duration_since(wave_start));
        let busy = ns(t0.elapsed());
        self.tasks.fetch_add(1, Ordering::Relaxed);
        self.busy.fetch_add(busy, Ordering::Relaxed);
        self.longest.fetch_max(busy, Ordering::Relaxed);
        self.latest_start.fetch_max(began, Ordering::Relaxed);
        self.last_end
            .fetch_max(began.saturating_add(busy), Ordering::Relaxed);
        result
    }

    // No NCBI counterpart: adds the wave to the run's counters (a wave without tasks: nothing).
    pub(crate) fn end(self) {
        let wall = since(self.start);
        let tasks = self.tasks.into_inner();
        if self.start.is_none() || tasks == 0 {
            return;
        }
        for (stage, value) in [
            (Run::SubjectWaves, 1),
            (Run::SubjectTasks, tasks),
            (Run::SubjectWall, wall),
            (Run::SubjectBusy, self.busy.into_inner()),
            (Run::SubjectLongest, self.longest.into_inner()),
            (Run::SubjectLatestStart, self.latest_start.into_inner()),
            (Run::SubjectLastEnd, self.last_end.into_inner()),
        ] {
            RUN[stage as usize].fetch_add(value, Ordering::Relaxed);
        }
    }
}

// No NCBI counterpart: adds one subject flow (`x_subject_flow::run`) to the run's counters.
pub(crate) fn subject_flow(start: Option<Instant>, tasks: u64, (ahead, taken): (u64, u64)) {
    if start.is_none() {
        return;
    }
    for (stage, value) in [
        (Run::FlowLoops, 1),
        (Run::FlowTasks, tasks),
        (Run::FlowAhead, ahead),
        (Run::FlowTaken, taken),
        (Run::FlowWall, since(start)),
    ] {
        RUN[stage as usize].fetch_add(value, Ordering::Relaxed);
    }
}

// No NCBI counterpart: adds the time since `start` to a stage of this thread's batch mode.
#[inline]
pub(crate) fn batch(stage: Batch, start: Option<Instant>) {
    if start.is_some() {
        BATCH[mode()][stage as usize].fetch_add(since(start), Ordering::Relaxed);
    }
}

// No NCBI counterpart: counts one batch of this thread's mode.
#[inline]
pub(crate) fn batch_count() {
    if on() {
        BATCH[mode()][Batch::Batches as usize].fetch_add(1, Ordering::Relaxed);
    }
}

// No NCBI counterpart: adds the time since `start` to an interval of `batch_end`.
#[inline]
pub(crate) fn end(stage: End, start: Option<Instant>) {
    if start.is_some() {
        END[mode()][stage as usize].fetch_add(since(start), Ordering::Relaxed);
    }
}

// No NCBI counterpart: adds `n` to a count of `batch_end`.
#[inline]
pub(crate) fn end_add(stage: End, n: u64) {
    if on() {
        END[mode()][stage as usize].fetch_add(n, Ordering::Relaxed);
    }
}

// No NCBI counterpart: seconds of a counter.
fn s(counter: &AtomicU64) -> f64 {
    counter.load(Ordering::Relaxed) as f64 / 1e9
}

// No NCBI counterpart: a count.
fn n(counter: &AtomicU64) -> u64 {
    counter.load(Ordering::Relaxed)
}

/// Writes the `[TIMING]` lines of the run (nothing without `LOSAT_TIMING`).
pub(crate) fn write(out: &mut dyn Write, threads: usize) -> io::Result<()> {
    if !on() {
        return Ok(());
    }
    let r = |stage: Run| &RUN[stage as usize];
    let input = (s(r(Run::Loop)) - s(r(Run::Search)) - s(r(Run::Render))).max(0.0);
    writeln!(
        out,
        "[TIMING] blastx_run: startup {:.3}s prolog {:.3}s read_stores {:.3}s loop {:.3}s epilog {:.3}s total {:.3}s (search {:.3}s render {:.3}s titles {:.3}s input {:.3}s)",
        s(r(Run::Startup)),
        s(r(Run::Prolog)),
        s(r(Run::ReadStores)),
        s(r(Run::Loop)),
        s(r(Run::Epilog)),
        s(r(Run::Total)),
        s(r(Run::Search)),
        s(r(Run::Render)),
        s(r(Run::Titles)),
        input,
    )?;
    let wall = s(r(Run::WaveWall));
    let busy = s(r(Run::WaveBusy));
    writeln!(
        out,
        "[TIMING] blastx_waves: waves={} batches={} wall {:.3}s busy {:.3}s threads={} idle {:.3}s tail={} tail_wall {:.3}s",
        n(r(Run::Waves)),
        n(r(Run::WaveBatches)),
        wall,
        busy,
        threads,
        (wall * threads as f64 - busy).max(0.0),
        n(r(Run::TailBatches)),
        s(r(Run::TailWall)),
    )?;
    writeln!(
        out,
        "[TIMING] blastx_flow: batches={} alone={} main_jobs={} main_job_time {:.3}s blocked {:.3}s",
        n(r(Run::FlowBatches)),
        n(r(Run::FlowAlone)),
        n(r(Run::FlowJobs)),
        s(r(Run::FlowJobTime)),
        s(r(Run::FlowBlocked)),
    )?;
    let wall = s(r(Run::SubjectWall));
    let busy = s(r(Run::SubjectBusy));
    writeln!(
        out,
        "[TIMING] blastx_subject_waves: waves={} tasks={} wall {:.3}s busy {:.3}s longest {:.3}s latest_start {:.3}s last_end {:.3}s threads={} idle {:.3}s",
        n(r(Run::SubjectWaves)),
        n(r(Run::SubjectTasks)),
        wall,
        busy,
        s(r(Run::SubjectLongest)),
        s(r(Run::SubjectLatestStart)),
        s(r(Run::SubjectLastEnd)),
        threads,
        (wall * threads as f64 - busy).max(0.0),
    )?;
    writeln!(
        out,
        "[TIMING] blastx_subject_flow: loops={} tasks={} ahead={} taken={} wall {:.3}s",
        n(r(Run::FlowLoops)),
        n(r(Run::FlowTasks)),
        n(r(Run::FlowAhead)),
        n(r(Run::FlowTaken)),
        s(r(Run::FlowWall)),
    )?;
    for (m, name) in [(0, "side"), (1, "alone")] {
        let b = |stage: Batch| &BATCH[m][stage as usize];
        writeln!(
            out,
            "[TIMING] blastx_batch_{name}: batches={} query_setup {:.3}s setup {:.3}s chunks_par {:.3}s lookup {:.3}s scan_par {:.3}s scan_inline {:.3}s replay {:.3}s chunk_end {:.3}s batch_end {:.3}s sum {:.3}s",
            n(b(Batch::Batches)),
            s(b(Batch::QuerySetup)),
            s(b(Batch::Setup)),
            s(b(Batch::ChunksPar)),
            s(b(Batch::Lookup)),
            s(b(Batch::ScanPar)),
            s(b(Batch::ScanInline)),
            s(b(Batch::Replay)),
            s(b(Batch::ChunkEnd)),
            s(b(Batch::BatchEnd)),
            s(b(Batch::Sum)),
        )?;
        let e = |stage: End| &END[m][stage as usize];
        writeln!(
            out,
            "[TIMING] blastx_batch_end_{name}: init {:.3}s read {:.3}s kappa_spec {:.3}s kappa_serial {:.3}s loop {:.3}s finish {:.3}s (matches={} reused={} serial={})",
            s(e(End::Init)),
            s(e(End::Read)),
            s(e(End::KappaSpec)),
            s(e(End::KappaSerial)),
            s(e(End::Loop)),
            s(e(End::Finish)),
            n(e(End::Matches)),
            n(e(End::Reused)),
            n(e(End::Serial)),
        )?;
    }
    super::seed::x_seed_diag::x_write_timing(out)?;
    Ok(())
}
