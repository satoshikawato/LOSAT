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
//!
//! Exactness. A subject's computation depends on the subject and the query chunk alone, not on the
//! diagonal table it runs on, while no diagonal offset passes 2^30 (`x_subject_flow.rs`, module
//! comment); both switches act only where `x_subject_flow::offsets_stay_below_2_30` holds. The batch
//! flow writes each batch's reports in the reader's order, and a capped batch is searched by the
//! same function as the batches the flow searches with the pool inside; the chunk waves hand the
//! sink every chunk and subject in the order of the serial chunk loop, as before. Only which thread
//! computes what, and when, changes.

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

/// No NCBI counterpart: `Some(K)` where `LOSAT_X_BXUNITCAP[=K]` caps the units searched at a time,
/// each unit with `tasks` tasks for the pool: the switch is on, there are more than K threads, K
/// units have a task for every thread (else threads would wait that one unit per thread keeps busy),
/// and the workers keep their diagonal tables (else each of the K units would make a table per
/// thread).
fn unit_cap(
    threads: usize,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
    tasks: usize,
) -> Option<usize> {
    static K: OnceLock<Option<usize>> = OnceLock::new();
    let k = (*K
        .get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BXUNITCAP", 2..=4_096, UNIT_CAP)))?;
    let own = threads > k && k.saturating_mul(tasks) >= threads && keeps_tables(options, subjects);
    crate::utils::xreach::x_reach!("LOSAT_X_BXUNITCAP", own);
    own.then_some(k)
}

/// No NCBI counterpart: `LOSAT_X_BXUNITCAP`'s K for the batch flow (`native.rs`), where each capped
/// batch's subjects would run as the subject flow (`x_subject_flow::fits`) in tasks of its `k`.
pub(crate) fn flow_unit_cap(
    pool: &crate::utils::threading::SearchPool<'_>,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<usize> {
    let tasks = super::x_subject_flow::subjects_per_task()
        .filter(|&k| super::x_subject_flow::fits(k, pool, false, options, subjects))
        .map_or(0, |k| subjects.len().div_ceil(k));
    unit_cap(pool.threads(), options, subjects, tasks)
}

/// No NCBI counterpart: `LOSAT_X_BXUNITCAP`'s K for the chunk waves of a split query
/// (`search::search_core`), each chunk's subjects in tasks of `UNIT_CAP_SUBJECTS_PER_TASK`.
pub(crate) fn chunk_unit_cap(
    pool: &crate::utils::threading::SearchPool<'_>,
    options: &ResolvedOptions,
    subjects: Subjects<'_>,
) -> Option<usize> {
    let tasks = subjects.len().div_ceil(UNIT_CAP_SUBJECTS_PER_TASK);
    unit_cap(pool.threads(), options, subjects, tasks)
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
/// `search::x_compute_chunk`'s subject loop as tasks of `UNIT_CAP_SUBJECTS_PER_TASK` consecutive
/// subjects on the pool (`LOSAT_X_BXUNITCAP`): the subjects with something for the sink, in subject
/// order, or the error of the first subject that failed.
#[allow(clippy::too_many_arguments)]
pub(super) fn chunk_subjects_on_pool(
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
    let k = UNIT_CAP_SUBJECTS_PER_TASK;
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
