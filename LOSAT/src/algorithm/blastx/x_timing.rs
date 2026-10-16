//! Diagnostic (LOSAT_TIMING): where the wall time of a BLASTX run goes, so that the serial parts
//! of a run with many threads have numbers.
//!
//! No NCBI counterpart: NCBI BLAST+ has no such timer. The module reads a monotonic clock and
//! adds differences to counters; it computes no value NCBI computes and no value the search
//! reads. With `LOSAT_TIMING` not set it never reads the clock.
//!
//! # What is measured
//!
//! `native::run` prints four lines after the epilog (`write`):
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
}
const RUN_LEN: usize = Run::TailWall as usize + 1;

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
    Ok(())
}
