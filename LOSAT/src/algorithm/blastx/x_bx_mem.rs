//! BLASTX's memory at many threads: the diagonal tables of the pool's workers kept one per thread
//! (`LOSAT_X_BXDIAGKEEP`), and at most K query batches or chunks searched at a time
//! (`LOSAT_X_BXUNITCAP[=K]`). Both are in the strict set (on unless `=0` or `LOSAT_X_DEFAULT=0`);
//! the third part, the lookup tables' neighbour lists shared by all threads and their vectors
//! kept between builds, is `LOSAT_X_BXLUTPOOL` (`tblastx/lookup/x_lut_direct.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-271
//! ```c
//!         for (Uint4 i = 0; i < query_splitter->GetNumberOfChunks(); i++) {
//! ...
//!                 if (IsMultiThreaded()) {
//!                      x_LaunchMultiThreadedSearch(*chunk_data);
//!                 } else {
//!                     retval =
//!                         CPrelimSearchRunner(*chunk_data, opts_memento.get())();
//! ```
//! NCBI builds one lookup table per query chunk and searches the chunks one after another; with
//! `-subject` it runs one thread. LOSAT searches query batches (`LOSAT_X_BXFLOW`) and the chunks of
//! a split query (`LOSAT_X_BXPAR`) side by side, one per thread, each with a lookup table of its
//! own (about 3-4 MiB for a batch of 10,002 bases), and a batch searched with the pool inside makes
//! one diagonal table per pool thread. Every thread that ever built a table keeps that much in its
//! glibc arena after the table is freed, so the maximum RSS grew by 4-9 MiB per thread and passed
//! NCBI's `-subject` maximum at 17-20 threads (v4-MEM REPORT §9.3, §10.3).
//!
//! - `LOSAT_X_BXDIAGKEEP`: where the subject flow (`x_subject_flow.rs`) searches a chunk's
//!   subjects, each helper takes the diagonal table its thread keeps (one of the chunk's length and
//!   window, else a new one) and gives it back when it ends, instead of a table per helper made by
//!   the searching thread for each chunk; the chunks' subject tasks of `LOSAT_X_BXUNITCAP` do the
//!   same.
//! - `LOSAT_X_BXUNITCAP[=K]` (built-in K = 8): with more than K threads, where K units have tasks
//!   for all threads, at most K query batches are searched at a time by the batch flow, each with the
//!   pool inside (the subject flow, `flow_unit_cap`), instead of one batch per thread; the chunks of a
//!   split query are searched in waves of at most K chunks, each chunk's subjects as tasks of
//!   `UNIT_CAP_SUBJECTS_PER_TASK` on the pool (`chunk_unit_cap`), instead of waves of one chunk per
//!   thread. A wave of a split batch's chunks inside the capped batch flow has K + 1 chunks less the
//!   capped batches searched and the chunks other waves hold at its start (at least one; it holds its
//!   chunks until it ends, `chunk_wave`), so the lookup tables alive at once stay near K (at most
//!   about 2K) whatever the number of threads. The capped flow keeps the
//!   threads busy: it counts the batches being searched, not the batches not yet written (up to
//!   `flow_bound` of those, so a slow batch at the head of the output does not stop the others), a
//!   search that ends starts the next batch itself, and the subject flow's helpers inside a capped
//!   batch leave once every task of the chunk is claimed (`in_capped_batch`) instead of waiting for
//!   the batch's thread to finish. These counts are the run's: the native run attaches them to its
//!   pool (`XCapRun`), so searches side by side in one process (lib tests, an embedding host) do
//!   not see each other's batches; a search on a pool without them counts process-wide.
//! - `LOSAT_X_BXCAPFILL[=f]` (built-in f = 4; acts only where `LOSAT_X_BXUNITCAP` caps): the capped
//!   units keep the threads busy with few subjects too. A capped unit's subjects run as tasks of
//!   `capped_subjects_per_task` subjects: f tasks for each of the `ceil(threads / K)` threads a unit
//!   is left, at least one subject and at most the size without the cap (64 for the batch flow's
//!   subject flow, 16 for the chunk waves; a set keeps that size from 64 x f x ceil(threads / K)
//!   subjects, 512 at 16 threads), and the cap holds wherever there are more threads than K
//!   (`unit_cap`; without the fill it also needs K units with a task for every thread). And the
//!   capped batch flow prepares each query batch of at most two search units' letters
//!   (translation, filtering, Karlin blocks: `reader_batch::prepare_lone_unit`, which
//!   `search::search_core` runs when a batch prepares itself) as a task of its own when the reader
//!   delivers it, as it prepares the units of a reader batch of several (`native.rs`
//!   `XUnit::split`), so that the lookup table's build and the subject search are what the cap
//!   holds to K at a time (a batch spawned before its preparation ends holds one of the K places
//!   while it waits). Without it, a unit of 64 subjects or fewer had no subject flow, so the cap
//!   left the batch flow alone and every thread searched its own batch with its own table; forced
//!   on, the cap ran the preparation (about half of the work at one thread) on K threads. Measured
//!   (ws node 0, untimed, 800 cDNAs x 50 yeast proteins `-subject`, 16 threads): one batch per
//!   thread 0.17 s and 128 MiB; the cap with tasks of 4 subjects but the preparation inside it
//!   x1.59 and 68 MiB; this switch x0.98 and 71 MiB (NCBI `-subject`, one thread: 6.0 s, 80 MiB);
//!   200 proteins x0.74. With `LOSAT_TIMING` the preparation ahead is not in the batch's
//!   `QuerySetup` clock.
//! - `LOSAT_X_BXBUILDPAR` (acts where `LOSAT_X_BXUNITCAP` caps): the lookup table of a capped unit
//!   is built on `ceil(threads / K)` tasks of the pool (the threads a unit is left), its passes by
//!   ranges of source cells (`tblastx/lookup/x_lut_direct.rs` `x_partitioned_passes`, which gives
//!   the argument that the table is the one-thread build's), instead of on the unit's thread
//!   alone: with K units at a time the builds were the part of the unit no other thread could
//!   take (2.7 ms of a 3.6 ms batch at 32 threads, 800 cDNAs x 8 yeast proteins). Build alone (ws,
//!   a 10,002-base batch's six frames): 7.8 ms on one thread, 1.8 ms in 4 parts on 4 threads.
//!
//! Exactness. A subject's computation depends on the subject and the query chunk alone, not on the
//! diagonal table it runs on, while no diagonal offset passes 2^30 (`x_subject_flow.rs`, module
//! comment); both switches act only where `x_subject_flow::offsets_stay_below_2_30` holds. The batch
//! flow writes each batch's reports in the reader's order, and a capped batch is searched by the
//! same function as the batches the flow searches with the pool inside; the chunk waves hand the
//! sink every chunk and subject in the order of the serial chunk loop, as before. A batch prepared
//! ahead (`LOSAT_X_BXCAPFILL`) is prepared by the statements its search runs when it prepares
//! itself, from the same records and options, and its search takes that preparation instead (a
//! preparation that fails is dropped, and the search prepares again and meets the error where it
//! met it before). Only which thread computes what, and when, changes.

use super::{compute_subject, SubjectComputation};
use crate::algorithm::blastx::{
    args::ResolvedOptions,
    parameters::ContextParameters,
    query_setup::PreparedQueryBatch,
    seed::{Diagonals, Lookup},
    subjects::Subjects,
};
use anyhow::Result;
use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: whether `LOSAT_X_BXDIAGKEEP` is on.
fn diag_keep_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXDIAGKEEP"))
}

/// No NCBI counterpart: `diag_keep`'s rule, without recording a side.
fn keeps_tables(options: &ResolvedOptions, subjects: Subjects<'_>) -> bool {
    diag_keep_on()
        && super::x_subject_flow::offsets_stay_below_2_30(
            options.window_size,
            subjects.longest_length(),
        )
}

/// No NCBI counterpart: whether the pool's workers keep one diagonal table per thread
/// (`LOSAT_X_BXDIAGKEEP`, module comment), and so may run on a table another chunk or batch used.
pub(crate) fn diag_keep(options: &ResolvedOptions, subjects: Subjects<'_>) -> bool {
    if !diag_keep_on() {
        return false;
    }
    let own = keeps_tables(options, subjects);
    crate::utils::xreach::x_reach!("LOSAT_X_BXDIAGKEEP", own);
    own
}

/// No NCBI counterpart: the built-in K of `LOSAT_X_BXUNITCAP` (units searched at a time).
const UNIT_CAP: usize = 8;

/// No NCBI counterpart: subjects per task of a chunk searched with the pool inside under
/// `LOSAT_X_BXUNITCAP`; 543 subjects make 34 tasks, enough for the threads a chunk is left.
const UNIT_CAP_SUBJECTS_PER_TASK: usize = 16;

/// No NCBI counterpart: K of `LOSAT_X_BXUNITCAP[=K]`; `None` when the switch is off.
fn cap_k() -> Option<usize> {
    static K: OnceLock<Option<usize>> = OnceLock::new();
    *K.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BXUNITCAP", 2..=4_096, UNIT_CAP))
}

/// No NCBI counterpart: the built-in f of `LOSAT_X_BXCAPFILL` (tasks per thread of the capped units).
const CAP_FILL: usize = 4;

/// No NCBI counterpart: f of `LOSAT_X_BXCAPFILL[=f]`; `None` when the switch is off.
fn cap_fill() -> Option<usize> {
    static F: OnceLock<Option<usize>> = OnceLock::new();
    *F.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BXCAPFILL", 2..=64, CAP_FILL))
}

/// No NCBI counterpart (scheduling only, `LOSAT_X_BXCAPFILL`): subjects per task of a unit searched
/// under the cap K with `threads` threads: `most` (the size without the cap) when the switch is
/// off; otherwise few enough that the unit has f tasks for each of the `ceil(threads / K)` threads
/// it is left, at least one subject and at most `most`.
fn capped_subjects_per_task(threads: usize, k: usize, subjects: usize, most: usize) -> usize {
    cap_fill().map_or(most, |f| {
        filled_subjects_per_task(f, threads, k, subjects, most)
    })
}

/// No NCBI counterpart: `capped_subjects_per_task` with f tasks per thread.
fn filled_subjects_per_task(
    f: usize,
    threads: usize,
    k: usize,
    subjects: usize,
    most: usize,
) -> usize {
    subjects
        .div_ceil(f.saturating_mul(threads.div_ceil(k.max(1))).max(1))
        .clamp(1, most.max(1))
}

/// No NCBI counterpart: `Some(K)` where `LOSAT_X_BXUNITCAP[=K]` caps the units searched at a time,
/// each unit with `tasks(K)` tasks for the pool: the switch is on, there are more than K threads, K
/// units have a task for every thread (else threads would wait that one unit per thread keeps busy;
/// with `LOSAT_X_BXCAPFILL` this condition is dropped), and the workers keep their diagonal tables
/// (else each of the K units would make a table per thread). `LOSAT_X_BXCAPFILL` acts exactly where
/// the cap does.
fn unit_cap(
    threads: usize,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
    tasks: impl FnOnce(usize) -> usize,
) -> Option<usize> {
    let k = cap_k()?;
    // With LOSAT_X_BXCAPFILL the cap holds wherever there are more threads than K: K units of few
    // subjects leave threads without a subject task, but those threads take the parts of the K
    // units' lookup-table builds (LOSAT_X_BXBUILDPAR) and the query batches prepared ahead, and a
    // unit per thread would make a lookup table per thread.
    let filled = cap_fill().is_some();
    let own = threads > k
        && (filled || k.saturating_mul(tasks(k)) >= threads)
        && keeps_tables(options, subjects);
    crate::utils::xreach::x_reach!("LOSAT_X_BXUNITCAP", own);
    if cap_fill().is_some() {
        crate::utils::xreach::x_reach!("LOSAT_X_BXCAPFILL", own);
    }
    own.then_some(k)
}

/// No NCBI counterpart: `LOSAT_X_BXUNITCAP`'s K for the batch flow (`native.rs`), where each capped
/// batch's subjects would run as the subject flow (`x_subject_flow::fits`) in tasks of its `k`
/// (`LOSAT_X_BXCAPFILL`: of `capped_subjects_per_task`).
pub(crate) fn flow_unit_cap(
    pool: &crate::utils::threading::SearchPool<'_>,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<usize> {
    let threads = pool.threads();
    let cap = unit_cap(threads, options, subjects, |cap| {
        super::x_subject_flow::subjects_per_task()
            .map(|most| capped_subjects_per_task(threads, cap, subjects.len(), most))
            .filter(|&k| super::x_subject_flow::fits(k, pool, false, options, subjects))
            .map_or(0, |k| subjects.len().div_ceil(k))
    });
    cap
}

/// No NCBI counterpart: subjects per task of the subject flow in a chunk of a batch searched with
/// the pool inside: `most`, or inside the capped batch flow (`in_capped_batch`) the size of
/// `flow_unit_cap`'s tasks.
pub(crate) fn flow_subjects_per_task(
    pool: &crate::utils::threading::SearchPool<'_>,
    subjects: usize,
    most: usize,
) -> usize {
    match cap_k() {
        Some(k) if in_capped_batch(pool) => {
            capped_subjects_per_task(pool.threads(), k, subjects, most)
        }
        _ => most,
    }
}

/// No NCBI counterpart: `LOSAT_X_BXUNITCAP`'s K for the chunk waves of a split query
/// (`search::search_core`) and the subjects per task of each chunk on the pool:
/// `UNIT_CAP_SUBJECTS_PER_TASK` (`LOSAT_X_BXCAPFILL`: at most that, `capped_subjects_per_task`).
pub(crate) fn chunk_unit_cap(
    pool: &crate::utils::threading::SearchPool<'_>,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<(usize, usize)> {
    let threads = pool.threads();
    let per_task =
        |k| capped_subjects_per_task(threads, k, subjects.len(), UNIT_CAP_SUBJECTS_PER_TASK);
    let k = unit_cap(threads, options, subjects, |k| {
        subjects.len().div_ceil(per_task(k))
    })?;
    Some((k, per_task(k)))
}

/// No NCBI counterpart: whether the capped batch flow prepares each query batch as a task of its
/// own outside the cap (`LOSAT_X_BXCAPFILL`, module comment).
pub(crate) fn prepares_batches_ahead() -> bool {
    cap_fill().is_some()
}

/// No NCBI counterpart: whether `LOSAT_X_BXBUILDPAR` is on.
fn build_par_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXBUILDPAR"))
}

/// No NCBI counterpart (scheduling only, `LOSAT_X_BXBUILDPAR`): the parts of the lookup table of
/// a unit of `pool`'s search: `ceil(threads / K)` (the threads a unit is left) for a unit the cap K
/// holds (`capped`: a chunk of the capped chunk waves, or any chunk inside the capped batch flow),
/// else 1 (the build on one thread).
pub(crate) fn build_parts(pool: &crate::utils::threading::SearchPool<'_>, capped: bool) -> usize {
    if !build_par_on() {
        return 1;
    }
    let k = cap_k().filter(|_| capped || in_capped_batch(pool));
    let parts = k.map_or(1, |k| pool.threads().div_ceil(k));
    crate::utils::xreach::x_reach!("LOSAT_X_BXBUILDPAR", parts > 1);
    parts
}

/// No NCBI counterpart (scheduling only, `LOSAT_X_BXBUILDPAR`): the parts of the first chunk's lookup
/// table in the flow across chunks (`LOSAT_X_BXCHUNKCROSS`, `x_subject_flow::run_chunks`), built
/// before the helpers start, while the other threads of the pool are idle: one per thread (at most
/// `X_MOST_PARTS` = 8, `x_lut_direct.rs`). The later chunks' tables are built on one thread.
pub(crate) fn cross_build_parts(pool: &crate::utils::threading::SearchPool<'_>) -> usize {
    if !build_par_on() {
        return 1;
    }
    let parts = pool.threads();
    crate::utils::xreach::x_reach!("LOSAT_X_BXBUILDPAR", parts > 1);
    parts
}

// No NCBI counterpart: the counts of one run's capped batch flow (scheduling only): the batches it
// is searching now and the chunks the waves of its split batches hold. The native run attaches one
// to its pool (`XCapRun::attach`); a search on a pool without one counts in `X_PROCESS_RUN`.
#[derive(Default)]
pub(crate) struct XCapRun {
    batches: AtomicUsize,
    wave_chunks: AtomicUsize,
}

// No NCBI counterpart: the counts of the searches whose pool carries no `XCapRun`.
static X_PROCESS_RUN: XCapRun = XCapRun {
    batches: AtomicUsize::new(0),
    wave_chunks: AtomicUsize::new(0),
};

impl XCapRun {
    /// No NCBI counterpart: `pool` with `self` attached, so that the run's batches count here and
    /// not with another search in the process (lib tests in parallel, an embedding host).
    pub(crate) fn attach<'b>(
        &'b self,
        pool: &'b crate::utils::threading::SearchPool<'_>,
    ) -> crate::utils::threading::SearchPool<'b> {
        pool.x_with_ext(self)
    }

    /// No NCBI counterpart: the counts of the run searched on `pool`.
    fn of<'a>(pool: &crate::utils::threading::SearchPool<'a>) -> &'a XCapRun {
        pool.x_ext::<XCapRun>().unwrap_or(&X_PROCESS_RUN)
    }
}

/// No NCBI counterpart: counts one batch of the capped batch flow while it is searched.
pub(crate) struct XCappedBatch<'a>(&'a XCapRun);

impl<'a> XCappedBatch<'a> {
    pub(crate) fn enter(pool: &crate::utils::threading::SearchPool<'a>) -> Self {
        let run = XCapRun::of(pool);
        run.batches.fetch_add(1, Ordering::Relaxed);
        Self(run)
    }
}

impl Drop for XCappedBatch<'_> {
    fn drop(&mut self) {
        self.0.batches.fetch_sub(1, Ordering::Relaxed);
    }
}

/// No NCBI counterpart: the chunks of one wave of a split query, held until the wave ends.
pub(crate) struct XChunkWave<'a>(usize, &'a XCapRun);

impl XChunkWave<'_> {
    pub(crate) fn size(&self) -> usize {
        self.0
    }
}

impl Drop for XChunkWave<'_> {
    fn drop(&mut self) {
        self.1.wave_chunks.fetch_sub(self.0, Ordering::Relaxed);
    }
}

/// No NCBI counterpart: the next wave of a split query under the cap K: K + 1 chunks less the capped
/// batches searched now (the batch itself among them) and the chunks other waves hold, at least
/// one, at most K. Waves that start at once take their chunks one after another, and a wave keeps
/// its chunks until it ends, so a batch that starts later does not take K as well. The counts are
/// those of the run searched on `pool`.
pub(crate) fn chunk_wave<'a>(
    pool: &crate::utils::threading::SearchPool<'a>,
    k: usize,
) -> XChunkWave<'a> {
    let run = XCapRun::of(pool);
    let batches = run.batches.load(Ordering::Relaxed);
    let mut held = run.wave_chunks.load(Ordering::Relaxed);
    loop {
        let size = (k + 1).saturating_sub(batches + held).clamp(1, k);
        match run.wave_chunks.compare_exchange_weak(
            held,
            held + size,
            Ordering::Relaxed,
            Ordering::Relaxed,
        ) {
            Ok(_) => return XChunkWave(size, run),
            Err(now) => held = now,
        }
    }
}

/// No NCBI counterpart: whether a batch of the capped batch flow of the run searched on `pool` is
/// being searched now (every batch of that run then is one).
pub(crate) fn in_capped_batch(pool: &crate::utils::threading::SearchPool<'_>) -> bool {
    XCapRun::of(pool).batches.load(Ordering::Relaxed) > 0
}

/// No NCBI counterpart: batches the capped batch flow may have read and not yet written, under the
/// cap K: four per batch searched at a time.
pub(crate) fn flow_bound(k: usize) -> usize {
    4 * k
}

thread_local! {
    // No NCBI counterpart: the diagonal table this thread keeps between searches (buffer reuse).
    static X_DIAG: RefCell<Option<Diagonals>> = const { RefCell::new(None) };
}

/// No NCBI counterpart (buffer reuse, `LOSAT_X_BXDIAGKEEP`): a worker's diagonal table, taken on
/// first use from the thread that uses it (the table that thread keeps where its length and window
/// are the chunk's, else a new one) and given back to the thread that drops it.
pub(crate) struct XDiagSlot {
    table: Option<Diagonals>,
    length: usize,
    window: i32,
}

impl XDiagSlot {
    // No NCBI counterpart: buffer reuse only.
    /// A slot for a chunk of `length` letters and `window`; nothing is taken yet.
    pub(crate) fn new(length: usize, window: i32) -> Self {
        Self {
            table: None,
            length,
            window,
        }
    }

    // No NCBI counterpart: buffer reuse only (a new table is `Diagonals::new`, the NCBI port).
    /// The slot's table, taken from this thread on first use.
    pub(crate) fn get(&mut self) -> Result<&mut Diagonals> {
        if self.table.is_none() {
            let kept = X_DIAG.with(|kept| kept.borrow_mut().take());
            self.table = Some(match kept {
                Some(table) if table.x_fits(self.length, self.window) => table,
                _ => Diagonals::new(self.length, self.window)?,
            });
        }
        Ok(self.table.as_mut().expect("set above"))
    }
}

// No NCBI counterpart: buffer reuse only. A table left while a panic unwinds (written for a
// subject whose `finish_subject` did not run) is dropped, not kept.
impl Drop for XDiagSlot {
    fn drop(&mut self) {
        if let Some(table) = self.table.take() {
            if !std::thread::panicking() {
                let _ = X_DIAG.try_with(|kept| *kept.borrow_mut() = Some(table));
            }
        }
    }
}

/// No NCBI counterpart (buffer reuse): the diagonal table a task of the subject flow runs on: the
/// loop thread's own, a helper's own (made for the chunk), or a helper's thread's (`XDiagSlot`).
pub(crate) enum XScratch<'a> {
    Borrowed(&'a mut Diagonals),
    Owned(Diagonals),
    Slot(XDiagSlot),
}

impl XScratch<'_> {
    // No NCBI counterpart: buffer reuse only.
    /// The task's table.
    pub(crate) fn get(&mut self) -> Result<&mut Diagonals> {
        match self {
            Self::Borrowed(table) => Ok(table),
            Self::Owned(table) => Ok(table),
            Self::Slot(slot) => slot.get(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
// ```c
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// NCBI scans a chunk's subjects one after another on one diagonal table. Here they are tasks of
// consecutive subjects on the pool, each task on its thread's table (module comment: the table's
// history does not change a subject's computation), and the results are kept in subject order.
/// `search::x_compute_chunk`'s subject loop as tasks of `k` consecutive subjects on the pool
/// (`LOSAT_X_BXUNITCAP`, `chunk_unit_cap`): the subjects with something for the sink, in subject
/// order, or the error of the first subject that failed.
#[allow(clippy::too_many_arguments)]
pub(super) fn chunk_subjects_on_pool(
    k: usize,
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &Lookup,
    xdrop: i32,
    db_length: i64,
    x_containment: bool,
) -> Result<Vec<(usize, SubjectComputation)>> {
    use rayon::prelude::*;
    // The cap acts only where the workers keep their tables (`unit_cap`).
    crate::utils::xreach::x_reach!("LOSAT_X_BXDIAGKEEP", true);
    let last = batch.contexts.last().expect("BLASTX contexts");
    let length = last.offset + last.length;
    let n = subjects.len();
    // The fused loop's choice for the chunk's table, worked out once for every task (a thread
    // takes tasks of the wave's chunks in turn).
    let ratio = crate::algorithm::blastx::seed::x_seed_diag::x_shared_table_ratio(lookup, options);
    let parts: Vec<Result<Vec<(usize, SubjectComputation)>>> = (0..n.div_ceil(k))
        .into_par_iter()
        .with_max_len(1)
        .map(|task| {
            crate::algorithm::blastx::seed::x_seed_diag::x_take_table_ratio(lookup, ratio);
            let (first, end) = (task * k, n.min((task + 1) * k));
            // EXPERIMENT (LOSAT_X_BXDBSTREAM): the task's subjects read as one block.
            let block = subjects.x_read_range(first, end)?;
            let subjects = subjects.x_with_block(block.as_ref());
            let mut slot = XDiagSlot::new(length, options.window_size);
            let table = slot.get()?;
            let mut kept = Vec::new();
            for oid in first..end {
                let mut computed = compute_subject(
                    subjects,
                    batch,
                    params,
                    options,
                    lookup,
                    table,
                    xdrop,
                    false,
                    db_length,
                    oid,
                    x_containment,
                )?;
                computed.seeds = Vec::new();
                computed.initial = Vec::new();
                if !computed.is_empty() {
                    kept.push((oid, computed));
                }
            }
            Ok(kept)
        })
        .collect();
    let mut results = Vec::new();
    for part in parts {
        results.extend(part?);
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    // No NCBI counterpart (scheduling only): the capped units' task sizes of LOSAT_X_BXCAPFILL.
    #[test]
    fn capped_units_have_f_tasks_for_each_thread_they_are_left() {
        use super::filled_subjects_per_task as size;
        // 50 subjects, K = 8: 16 threads leave each unit 2, so 8 tasks of 7; 64 threads 8, so 32
        // tasks of 2; 6,600 subjects keep the size without the cap.
        assert_eq!(size(4, 16, 8, 50, 64), 7);
        assert_eq!(size(4, 64, 8, 50, 64), 2);
        assert_eq!(size(4, 16, 8, 6_600, 64), 64);
        assert_eq!(size(4, 16, 8, 6_600, 16), 16);
        // K units then hold a task for every thread (the rule `unit_cap` keeps without the fill)
        // wherever K units
        // of one subject each would.
        for threads in [9, 16, 17, 32, 48, 64, 128] {
            for subjects in [2, 10, 50, 200, 543, 6_600] {
                let per_task = size(4, threads, 8, subjects, 64);
                let tasks = subjects.div_ceil(per_task);
                assert!(
                    8 * tasks >= threads || 8 * subjects < threads,
                    "{threads} threads, {subjects} subjects: {per_task} per task"
                );
            }
        }
        // One subject at least, also for no subjects.
        assert_eq!(size(4, 16, 8, 0, 64), 1);
        assert_eq!(size(64, 4_096, 2, 1, 16), 1);
    }

    // No NCBI counterpart (scheduling only). The counts are the run's own (attached to its pool),
    // so other tests searching at the same time do not move them.
    #[test]
    fn a_wave_holds_its_chunks_until_it_ends() {
        let run = super::XCapRun::default();
        let pool = crate::utils::threading::SearchPool::x_serial();
        let pool = run.attach(&pool);
        assert!(!super::in_capped_batch(&pool));
        let batch = super::XCappedBatch::enter(&pool);
        assert!(super::in_capped_batch(&pool));
        // K + 1 less the batch itself, then one while the first wave holds its chunks.
        let first = super::chunk_wave(&pool, 8);
        assert_eq!(first.size(), 8);
        let second = super::chunk_wave(&pool, 8);
        assert_eq!(second.size(), 1);
        drop(first);
        assert_eq!(super::chunk_wave(&pool, 8).size(), 7);
        drop((second, batch));
        assert!(!super::in_capped_batch(&pool));
    }
}
